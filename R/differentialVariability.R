########################################################################################################################
## differentialVaribility.R
## created: 2017-28-07
## creator: Michael Scherer
## ---------------------------------------------------------------------------------------------------------------------
## The differential varibility analysis methods between sample groups.
########################################################################################################################

#' diffVar
#' 
#' This routine applies the diffVar method from the \code{missMethyl} package that determines sites exhibiting
#' differential varibility between two sample groups
#' 
#' @author Michael Scherer
#' @param meth.matrix Matrix containing the methylation information used to calculate differentially variable sites 
#'                    between the two groups
#' @param inds.g1 Indices in the phenotypic table corresponding to the first group.
#' @param inds.g2 Indices in the phenotypic table corresponding to the second group.
#' @param adjustment.table A \code{data.frame} containing variables to adjust for in the testing
#' @return P-values as the result of the diffVar method not adjusted for multiple hypothesis testing.
#' 
#' @export
#' @references Phipson, Belinda, Oshlack, Alicia (2014)
#'             DiffVar: a new method for detecting differential variability with application to methylation in cancer and aging
#'             Genome Biology 15(9):465.
diffVar <- function(meth.matrix,inds.g1,inds.g2,adjustment.table=NULL){
  logger.start("Started diffVar method")
  #rnb.require('missMethyl')
  if(is.logical(inds.g1)) inds.g1 <- which(inds.g1)
  if(is.logical(inds.g2)) inds.g2 <- which(inds.g2)
  n.g1 <- length(inds.g1)
  n.g2 <- length(inds.g2)
  n <- n.g1 + n.g2
  if (!is.null(adjustment.table)){
    if (!(is.data.frame(adjustment.table) && nrow(adjustment.table)==n && (!any(is.na(adjustment.table))))) {
      stop("invalid value for adjustment.table")
    }
    m <- ncol(adjustment.table)
    if (m == 0) {
      adjustment.table <- NULL
    } else {
      colnames(adjustment.table) <- paste0("x",1:m,"x")
    }
  }
  
  ind.vec <- c(inds.g1,inds.g2)
  if (length(ind.vec) < 2) stop("need at least two samples indices to compare")
  meth.matrix <- rnb.beta2mval(meth.matrix[,ind.vec,drop=FALSE])
  df <- data.frame(xg = factor(rep(c("group1","group2"), c(n.g1,n.g2)), levels=c("group1","group2")))
  if (!is.null(adjustment.table)){
    df <- cbind(df,adjustment.table)
  }
  formula.text <- paste0(c("~0",colnames(df)),collapse="+")
  design.m <- model.matrix(as.formula(formula.text),data=df)
  colnames(design.m) <- make.names(colnames(design.m),unique=TRUE)
  colnames(design.m)[1:2] <- c("group1","group2")
  var.fit <- varFit(data=meth.matrix,design=design.m)
  contrasts.m <- makeContrasts(group1vs2=group1-group2,levels=design.m)
  var.fit <- contrasts.fit(var.fit,contrasts.m)
  var.fit <- eBayes(var.fit)
  logger.completed()
  return(var.fit$p.value[,"group1vs2"])
}

#' apply.iEVORA
#' 
#' This routine applies the iEVORA method created by Teschendorff et.al. to the supplied methylation matrix in a similar way
#' as the diffVar method.
#' 
#' @autor Michael Scherer
#' @param meth.matrix Matrix containing the methylation information used to calculate differentially variable sites 
#'                    between the two groups
#' @param inds.g1 Indices in the phenotypic table corresponding to the first group.
#' @param inds.g2 Indices in the phenotypic table corresponding to the second group.
#' @return Q-values as the result of applying the iEVORA method and then correct for multiple testing.
#' @export
apply.iEVORA <- function(meth.matrix,inds.g1,inds.g2){
  if(is.logical(inds.g1)) inds.g1 <- which(inds.g1)
  if(is.logical(inds.g2)) inds.g2 <- which(inds.g2)
  n <- length(inds.g1)+length(inds.g2)
  pheno.binary <- rep(0,n)
  pheno.binary[inds.g2] <- 1
  result <- suppressWarnings(iEVORA(data.m=meth.matrix,pheno.v=pheno.binary))
  p.values <- rep(1,dim(meth.matrix)[1])
  p.values[result[,"index"]] <- result[,"q(BT)"]
  return(p.values)
}

#' iEVORA
#' 
#' This method performs the iEVORA algorithm introduced by Teschendorff et.al. to detect sites showing differential
#' variability between two groups. The routine has been modified from the R script associated with the publication.
#' 
#' @author Andrew E. Teschendorff, with minor modifications by Michael Scherer
#' @param data.m Methylation matrix with rows labeling features, and columns labeling samples. 
#'               Rownames should be feature/probe IDs and should be provided.
#' @param pheno.v Phenotype vector with entries either 0 or 1 correpsonding to the group assignment.
#' @param thDV q-value threshold for the differntial varibility test. Default is 0.001.
#' @param thDM p-value threshold for differntial methylation means. Default is 0.05.
#' 
#' @return A matrix of ranked differentially variable (DV) and differentially methylated CpGs (DVMCs), 
#'          ranked according to the t-statistic p-value, but selected using the Bartlett's DV test. 
#'          Columns labels are the t-statistic, its P-value, the mean of phenotype-1, the mean of phenotype-0, 
#'          the log-ratio of the variances of phenotype-1 to phenotype-0,as well as the Bartlett's test P-value 
#'          and q-value.
#' @noRd
#' @references Teschendorff, Andrew E, Gao, Yang, Jones, Allison, Ruebner, Matthias, Beckmann, Matthias W., Wachter, David L.
#'              ,Fasching, Peter A., Widschwendter, Martin
#'             DNA methylation outliers in normal breast tissue identify field defects that are enriched in cancer
#'             Nature Communications 7:10478.
iEVORA <- function(data.m,pheno.v,thDV=0.001,thDM=0.05){
  logger.start("iEVORA method")
  #rnb.require("qvalue")
  statDVC.m <- t(apply(data.m,1,doDV,pheno.v));
  qvDVC.v <- qvalue(statDVC.m[,2])$qval;
  dvc.idx <- which(qvDVC.v < thDV);
  nDVC <- length(dvc.idx);
  if( nDVC > 0 ){
    statDMC.m <- t(apply(data.m[dvc.idx,],1,doTT,pheno.v))
    tmp.s <- sort(statDMC.m[,2],decreasing=FALSE,index.return=TRUE);
    pvDMC.v <- tmp.s$x;
    ntop <- length(which(pvDMC.v < 0.05));
    if(ntop > 0){
      topDVMC.m <- cbind(dvc.idx[tmp.s$ix[1:ntop]],statDMC.m[tmp.s$ix[1:ntop],],statDVC.m[dvc.idx[tmp.s$ix[1:ntop]],c(3:4,1:2)],qvDVC.v[dvc.idx[tmp.s$ix[1:ntop]]]);
      colnames(topDVMC.m) <- c("index","t","P(TT)","Av1","Av0","log[V1/V0]","P(BT)","q(BT)");
      rownames(topDVMC.m) <- rownames(statDMC.m)[tmp.s$ix[1:ntop]];
    }
    else {
      logger.info("No DVMCs found.")
    }
  }
  else {
    logger.info("No DVCs detected. Consider lowering the differential varibility threshold.");
  }
  logger.completed()
  return(topDVMC.m);
}

#' doDV
#' 
#' Helper function for the iEVORA method also supplied with the publication used to determine differentially variable
#' sites.
#' 
#' @autor Andrew Teschendorff, minor modifications by Michael Scherer
#' @param tmp.v A row from the methylation matrix corresponding to a single site.
#' @param pheno.v Vector containing the splitting into the two groups for which differential variability should be detected.
#' 
#' @return Vector containing the log odds of variances, the p-value from the Bartlett's test and the mean methylation 
#'          levels for the two groups.
#' @noRd
doDV <- function(tmp.v,pheno.v){
  co.idx <- which(pheno.v==0);
  ca.idx <- which(pheno.v==1);
  bt.o <- bartlett.test(x=tmp.v,g=pheno.v);
  pv <- bt.o$p.value;
  logR <- log2(var(tmp.v[ca.idx])/var(tmp.v[co.idx]));
  avCA <- mean(tmp.v[ca.idx]);
  avCO <- mean(tmp.v[co.idx]);
  out.v <- c(logR,pv,avCA,avCO);
  names(out.v) <- c("log(V1/V0)","P(BT)","Av1","Av0");
  return(out.v);
}

#' doTT
#' Helper function for the iEVORA method also supplied with the publication used to calculate t-statistics of differential
#' methylation between two sample groups.
#' 
#' 
#' @author Andrew Teschendorff, with minor modifications by Michael Scherer
#' @param tmp.v A row from the methylation matrix corresponding to a single site.
#' @param pheno.v Vector containing the splitting into the two groups for which differential variability should be detected.
#' 
#' @return Vector containing the t-statistics and the t-statistics p-value.
#' @noRd
doTT <- function(tmp.v,pheno.v){
  tt.o <- t.test(tmp.v ~ pheno.v);
  out.v <- c(-tt.o$stat,tt.o$p.val);
  names(out.v) <- c("t","P");
  return(out.v);
}

computeDiffVar.bin <- function(meth.matrix,inds.g1,inds.g2,adjustment.table=NULL){
  #ToDo
}

rnb.execute.diffVar <- function(rnb.set,pheno.cols,columns.adj=rnb.getOption("covariate.adjustment.columns"),
                                adjust.celltype=rnb.getOption("differential.adjustment.celltype"),
                                disk.dump=rnb.getOption("disk.dump.big.matrices"),
                                disk.dump.dir=tempfile(pattern="diffVarTables_")){
  
  logger.start("Retrieving comparison info")
  cmp.info <- get.comparison.info(rnb.set, pheno.cols=pheno.cols, 
                                  columns.adj=columns.adj,
                                  adjust.celltype=adjust.celltype)
  logger.completed()
  if (is.null(cmp.info)) {
    return(NULL)
  }
  
  variability.method <- rnb.getOption("variability.method")
  diff.meth <- new("RnBDiffMeth",variability.method=variability.method,disk.dump=disk.dump,disk.path=disk.dump.dir)
  for(i in 1:length(cmp.info)){
    cmp.info.cur <- cmp.info[[i]]
    logger.start(c("Comparing ",cmp.info.cur$comparison))
    meth.matrix <- meth(rnb.set)
    diffVar <- computeDiffVar.bin(meth.matrix=meth.matrix,inds.g1=cmp.info.cur$group.inds$group1,
                                  inds.g2=cmp.info.cur$group.inds$group2, adjustment.table=cmp.info.cur$adjustment.table)
    diff.meth <- addDiffVarTable(diff.meth,diffVar,cmp.info.cur$comparison,cmp.info.cur$group.names)
  }
  
  diff.meth <- addComparisonInfo(diff.meth,cmp.info)
  logger.completed()
  return(diff.meth)
}
