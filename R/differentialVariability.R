########################################################################################################################
## differentialVaribility.R
## created: 2017-28-07
## creator: Michael Scherer
## ---------------------------------------------------------------------------------------------------------------------
## The differential variability analysis methods between sample groups.
########################################################################################################################

P.VAL.CUT <- 0.05 #p-value cutoff to consider differentially variable sites as significant
DENS.SCATTER.SUBSAMPLE.THRES <- 2e6 #threshold to induce subsampling
DENS.SCATTER.SPARSE.POINTS.PERC <- 0.01 #percentage of points to plot in the sparsely populated regions
DENS.SCATTER.SPARSE.POINTS.MAX <- 1e4 #maximum number of points to plot in the sparsely populated regions
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
  logger.start("diffVar method")
  rnb.require('missMethyl')
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
#' @author Michael Scherer
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
  p.values[result[,"index"]] <- result[,"P(BT)"]
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
  rnb.require("qvalue")
  rnb.require("stats")
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
#' @author Andrew Teschendorff, minor modifications by Michael Scherer
#' @param tmp.v A row from the methylation matrix corresponding to a single site.
#' @param pheno.v Vector containing the splitting into the two groups for which differential variability should be detected.
#' 
#' @return Vector containing the log odds of variances, the p-value from the Bartlett's test and the mean methylation 
#'          levels for the two groups.
#' @noRd
doDV <- function(tmp.v,pheno.v){
  rnb.require('stats')
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

#'computeDiffVar.site
#'
#'This function performs the differential variability analysis on the given methylatio matrix and the group assignment
#'specified by the column indices.
#'@author Michael Scherer
#'@param meth.matrix Methylation matrix on which differential variability analysis should be conducted.
#'@param inds.g1 Column indices of the first group.
#'@param inds.g2 Column indices of the second group.
#'@param method Differrential analysis method to conduct. Possibilities are: "diffVar" and "iEVORA".
#'@param adjustment.table A \code{data.frame} containing variables to adjust for in the testing
#'@return Differential Variability data frame for the analysis.
#'@rdname computeDiffVar.site
#' @aliases computeDiffVar.site
#' @aliases computeDiffVar.default.site
#' @aliases computeDiffVar.extended.site
#'@export
computeDiffVar.default.site <- function(meth.matrix,inds.g1,inds.g2,method=rnb.getOption('differential.variability.method'),
                                adjustment.table=NULL,eps=0.01){
  if (!(method %in% c("diffVar","iEVORA"))) {
    stop("Invalid method for differential varibility method")
  }
  meth.g1 <- meth.matrix[,inds.g1]
  meth.g2 <- meth.matrix[,inds.g2]
  if(length(inds.g1)<2) {
    logger.info("Group 1 has less than 2 members")
    meth.g1 <- as.matrix(meth.g1)
  }
  if(length(inds.g2)<2) {
    logger.info("Group 2 has less than 2 members")
    meth.g2 <- as.matrix(meth.g2)
  }
  
  p.vals <- rep(as.double(NA),nrow(meth.matrix))
  do.p.vals <- ncol(meth.g1) > 1 || ncol(meth.g2) > 1
  if (do.p.vals) {
    if (method == "diffVar"){
      logger.info("Conducting differential variability using diffVar")
      tryCatch(
        p.vals <- diffVar(meth.matrix,inds.g1,inds.g2,adjustment.table=adjustment.table),
        error = function(ee) {
          logger.warning(c("Could not compute p-values using diffVar:",ee$message))
        }
      )
    } else if (method == "iEVORA"){
      logger.info("Conducting differential variability using iEVORA")
      tryCatch(
        p.vals <- apply.iEVORA(meth.matrix,inds.g1,inds.g2),
        error = function(ee) {
          logger.warning(c("Could not compute p-values using iEVORA:",ee$message))
        }
      )
    }
  } else {
    logger.warning("Skipping p-value computation due to insufficient sample numbers")
  }
  p.vals.is.na <- is.na(p.vals)
  if (!all(p.vals.is.na)){
    if (any(p.vals.is.na)){
      logger.info(c(sum(p.vals.is.na),"p-values are NA. They are treated as 1 in FDR adjustment"))
      p.vals[is.na(p.vals.t.na.adj)] <- 1
    }
    p.vals.adj <- p.adjust(p.vals, method = "fdr")
  } else {
    p.vals.adj <- rep(NA,length(p.vals))
  }

  var.g1 <- apply(meth.matrix[,inds.g1],1,var)
  var.g2 <- apply(meth.matrix[,inds.g2],1,var)
  var.diff <- var.g1-var.g2
  var.log.ratio <- log2(var.g1+eps/var.g2+eps)
  neg.log10.p <- -log10(p.vals)
  neg.log10.fdr <- -log10(p.vals.adj)
  tt <- data.frame(var.g1=var.g1,var.g2=var.g2,var.diff=var.diff,var.log.ratio=var.log.ratio,diffVar.p.val=p.vals,diffVar.p.adj.fdr=p.vals.adj,log10P=neg.log10.p,
                   log10FDR=neg.log10.fdr)
  return(tt)
}
#' @rdname computeDiffVar.site
#' @export
computeDiffVar.extended.site <- function(meth.matrix,inds.g1,inds.g2,method=rnb.getOption('differential.variability.method'),
                                        adjustment.table=NULL,covg=NULL,covg.thres=rnb.getOption("filtering.coverage.threshold")){
  matrix.default <- computeDiffVar.default.site(
    meth.matrix,inds.g1=inds.g1,inds.g2=inds.g2,
    method=method,adjustment.table=adjustment.table
  )
  #coverage information
  if (!is.null(covg) & all(dim(covg)==dim(meth.matrix))){
    covg[is.na(covg)] <- 0 #set NA to 0 coverage
    
    tab.covg.g1 <- covg[,inds.g1]
    tab.covg.g2 <- covg[,inds.g2]
    if(length(inds.g1)<2) {
      logger.info("Group 1 has less than 2 members")
      tab.covg.g1 <- as.matrix(tab.covg.g1)
    }
    if(length(inds.g2)<2) {
      logger.info("Group 2 has less than 2 members")
      tab.covg.g2 <- as.matrix(tab.covg.g2)
    }
    mean.covg.g1 <- rowMeans(tab.covg.g1)
    mean.covg.g2 <- rowMeans(tab.covg.g2)
    min.covg.g1 <- rowMins(tab.covg.g1)
    min.covg.g2 <- rowMins(tab.covg.g2)
    max.covg.g1 <- rowMaxs(tab.covg.g1)
    max.covg.g2 <- rowMaxs(tab.covg.g2)
    covg.thresh.nsamples.g1 <- rowSums(tab.covg.g1>=covg.thres)
    covg.thresh.nsamples.g2 <- rowSums(tab.covg.g2>=covg.thres)
    
    matrix.ext	<- data.frame(mean.covg.g1=mean.covg.g1,mean.covg.g2=mean.covg.g2,
                              min.covg.g1=min.covg.g1,min.covg.g2=min.covg.g2,
                              max.covg.g1=max.covg.g1,max.covg.g2=max.covg.g2,
                              covg.thresh.nsamples.g1=covg.thresh.nsamples.g1,covg.thresh.nsamples.g2=covg.thresh.nsamples.g2)
    matrix.default <- cbind(matrix.default,matrix.ext)	
  }
  return(matrix.default)  
}

#' addReportPlot.diffVar.scatter
#' 
#' This function creates a scatterplot for the given comparison comparing the difference in variances with the
#' usted) p-value as the result of the differentially variability analysis.
#' 
#' @param report Object of class \code{\linkS4class{Report}} to which the plot should be added
#' @param var.table Variability table containing the relevant information to be put into the scatterplot
#' @param comparison.name Comparison for which the scatterplot should be creaed
#' @param group.name1 Name of the first group of the comparison
#' @param group.name2 Name of the second group of the comparison
#' @return list of scatterplots
#' @noRd
addReportPlot.diffVar.scatter <- function (report, var.table, comparison.name,
                              group.name1="Group1",group.name2="Group2"){
  ret <- list()
  sparse.points <- DENS.SCATTER.SPARSE.POINTS.PERC
  if (DENS.SCATTER.SPARSE.POINTS.MAX < sparse.points*nrow(var.table)){
    sparse.points <- DENS.SCATTER.SPARSE.POINTS.MAX
  }
  dens.subsample <- FALSE
  if (nrow(var.table) > dens.subsample){
    dens.subsample <- DENS.SCATTER.SUBSAMPLE.THRES
  }
  
  if("diffVar.p.adj.fdr" %in% colnames(var.table)){
    var.sites <- var.table[,"diffVar.p.adj.fdr"] < P.VAL.CUT
    plot <- create.densityScatter(var.table[,c("var.g1","var.g2")],is.special=var.sites,dens.subsample=dens.subsample) +
        xlab(paste("Variance",group.name1)) + ylab(paste("Variance",group.name2))
    comp.type <- "fdrAdjPval"
    fig.name <- paste("diffVar",comparison.name,comp.type,sep = "_")
    plot <- createReportGgPlot(plot,fig.name,report=report,create.pdf = FALSE)
    plot <- off(plot,handle.errors=TRUE)
    ret <- c(ret,list(plot))
  }
  return(ret)
}

#' addReportPlot.diffVar.volcano
#' 
#' This function creates volcano plots for the given comparison comparing the difference in variances with the
#' negative decadic logarithm of the (adjusted) p-value as the result of the differentially variability analysis.
#' 
#' @param report Object of class \code{\linkS4class{Report}} to which the plot should be added
#' @param var.table Variability table containing the relevant information to be put into the scatterplot
#' @param comparison.name Comparison for which the scatterplot should be creaed
#' @param group.name1 Name of the first group of the comparison
#' @param group.name2 Name of the second group of the comparison
#' @return list of volcano plots
#' @noRd
addReportPlot.diffVar.volcano <- function (report, var.table, comparison.name,
                                           group.name1="Group1",group.name2="Group2"){
  ret <- list()
  dont.plot.p.val <- all(is.na(var.table[,"diffVar.p.val"]))
  
  fig.name <- paste("diffVar_volcano",comparison.name,"pVal",sep="_")
  if (!dont.plot.p.val){
    pp <- ggplot(var.table) + aes(x=var.diff,y=log10P,color=log10(combined.rank.var)) + 
      scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
      geom_point()
  } else {
    pp <- rnb.message.plot("No p-value available")
  }
  report.plot <- createReportGgPlot(pp,fig.name, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  ret <- c(ret,list(report.plot))
  
  fig.name <- paste("diffVar_volcano",comparison.name,"pValAdj",sep="_")
  pp <- ggplot(var.table) + aes(x=var.diff,y=log10FDR,color=log10(combined.rank.var)) +
    scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient")))+
    geom_point()
  report.plot <- createReportGgPlot(pp,fig.name, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  ret <- c(ret,list(report.plot))
  
  return(ret)
}

#' addReportPlot.diffVar.meth
#' 
#' This function compares differentially methylated sites (DMCs) with differentially variable sites (DMVs) in a scatterplot
#' and colors significantly differential sites with different colors.
#' 
#' @param report Object of class \code{\linkS4class{Report}} to which the plot should be added.
#' @param var.table Table containing the differentially variable sites.
#' @param meth.table Table containing the differentially methylated sites.
#' @param comparison.name Name of the comparison to be conducted.

addReportPlot.diffVar.meth <- function(report, var.table, meth.table, comparison.name,
                           group.name1,group.name2){
  ret <- list()
  dont.plot.p.val <- all(is.na(var.table[,"diffVar.p.adj.fdr"]))
  dont.plot.p.val <- dont.plot.p.val && all(is.na(meth.table[,"diffmeth.p.adj.fdr"]))
  
  fig.name <- paste("diffVarMeth",comparison.name,sep="_")
  if (!dont.plot.p.val){
    toPlot <- data.frame(Meth.p=meth.table[,"diffmeth.p.adj.fdr"],Var.p=var.table[,"diffVar.p.adj.fdr"])
    is.diff.meth <- toPlot$Meth.p < P.VAL.CUT
    is.diff.var <- toPlot$Var.p < P.VAL.CUT
    color.diff <- rep('Not Significant',dim(toPlot)[1])
    color.diff[is.diff.meth] <- "DMR"
    color.diff[is.diff.var] <- "DMV"
    color.diff[is.diff.meth&is.diff.var] <- "Both"
    toPlot <- data.frame(toPlot,Color.diff=color.diff)
    toPlot <- toPlot[is.diff.meth|is.diff.var,]
    pp <- ggplot(toPlot) + aes(x=-log10(Meth.p),y=-log10(Var.p),color=Color.diff) +
          scale_color_manual(values = c("#AF648C","#A6CEE3","#FB9A99")) +
          geom_vline(xintercept = -log10(P.VAL.CUT),linetype='dotted') + geom_hline(yintercept = -log10(P.VAL.CUT),linetype='dotted') +
          geom_point()
  } else {
    pp <- rnb.message.plot("No p-values available")
  }
  report.plot <- createReportGgPlot(pp,fig.name, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  ret <- c(ret,list(report.plot))
  
  return(ret)
}


#' create.diffvar.column.description
#' 
#' This function returns a detailed description of the column names for the differential variaiability table.
#' @param includeCovg Flag indicating if coverage information should also be explaines
#' @param covgThres Coverage threshold
#' @return List of column descriptions for each column in the differential variability table.
#' @noRd
create.diffvar.column.description <- function(includeCovg, covgThres=-1L){
  res <- list(
    "id: site id",
    "Chromosome: chromosome of the site",
    "Start: start position of the site",
    "Strand: strand of the site",
    c("var.g1, var.g2: (g1 and g2 are replaced by the corrspondinhg group used in the differentiality analysis) ",
        "the variances found in the groups"),
    "var.diff: difference in variance values between the two groups g1 and g2 (=var.g1-var.g2)",
    "var.log.ratio: Log2 of the ratio between the variances of the two groups g1 and g2 (=log2(var.g1+eps/var.g2+eps), default eps=0.01)",
    "diffvar.p.val: p-value resulting from applying the selected differentially variability method (diffVar or iEVORA)",
    "diffVar.p.adj.fdr: FDR-adjusted p-value for differential variability",
    "log10P: negative decadic logarithm of the p-value",
    "log10FDR: negative decadic logarithmn of the FDR-adjusted p-value",
    "combined.rank.var: combined rank consisting of the difference in variances, the log ration between the group variances and the p-value of the statistical test"
  )
  if (includeCovg){
    res <- c(res, list(
      paste0("mean.covg.g1,mean.covg.g2: mean coverage of groups 1 and 2 respectively (In case of Infinium array methylation data, coverage is defined as combined beadcount.)"),
      paste0("min.covg.g1,min.covg.g2: minimum coverage of groups 1 and 2 respectively"),
      paste0("max.covg.g1,max.covg.g2: maximum coverage of groups 1 and 2 respectively"),
      paste0("covg.thresh.nsamples.g1,covg.thresh.nsamples.g2: number of samples in group 1 and 2 respectively exceeding the coverage threshold (", covgThres, ").")
    ))
  }
  return(res)
}
#' get.diffvar.tab.annot.cols
#' 
#' This function returns the column names of the differential variability table.
#' @param includeCovg Flag indicating if coverage information is present in the data set
#' @param covgThres Coverage threshold
#' @return List of column names in the differential variability table
#' @noRd
get.diffvar.tab.annot.cols <- function(includeCovg, covgThres=-1L){
  res <- c()
  res <- c("var.g1","var.g2","var.diff",
           "var.log.ratio","diffVar.p.val",
           "diffVar.p.adj.fdr","log10P",
           "log10FDR","combined.rank.var")
  if (includeCovg){
    res <- c(res,c("mean.covg.g1","mean.covg.g2",
                   "min.covg.g1","min.covg.g2","max.covg.g1","max.covg.g2",
                   "covg.thresh.nsamples.g1","covg.thresh.nsamples.g2"))
  }
  return(res)
}

#' create.diffvar.tab.annot.colnames.pretty
#' 
#' This function returns a nice representation for the column names in the differential variaiability table.
#' @param grp.name1 Name of the first group in the comparison
#' @param grp.name2 Name of the second group in the comparison
#' @param includeCovg Flag indicating if coverage information should also be explaines
#' @param covgThres Coverage threshold
#' @return List of column descriptions for each column in the differential variability table.
#' @noRd  
get.diffvar.tab.annot.colnames.pretty <- function(grp.name1, grp.name2, includeCovg, covgThres=-1L){
  res <- c(paste("var",grp.name1,sep="."),paste("var",grp.name2,sep="."),"var.diff",
           "var.log.ratio","diffvar.p.val","diffvar.p.adj.fdr","log10P","log10FDR","combined.rank.var")
  if (includeCovg){
    res <- c(res,c(paste("mean.covg",grp.name1,sep="."),paste("mean.covg",grp.name2,sep="."),
                   paste("min.covg",grp.name1,sep="."),paste("min.covg",grp.name2,sep="."),
                   paste("max.covg",grp.name1,sep="."),paste("max.covg",grp.name2,sep="."),
                   paste("nsamples.covg",paste("thres",covgThres,sep=""),grp.name1,sep="."),
                   paste("nsamples.covg",paste("thres",covgThres,sep=""),grp.name2,sep=".")))
  }
  return(res)
}
#' rnb.section.diffVar
#'
#' This functions add information about the differential variability analysis to the specified report.
#' 
#' @param rnb.set Object of class \code{\linkS4class{RnBSet}} on which the analysis was conducted.
#' @param diff.meth Object of class \code{\linkS4class{RnBDiffMeth}} as the result of applying \code{rnb.execute.diffVar} 
#'                  containing the differentially variable sites
#' @param report Report object (\code{\linkS4class{Report}}) to which information should be added.
#' @param gzTable Flag indicating if the tables should ne gzipped
#' @return Modified report with section about result of differential variability analysis
#' @author Michael Scherer
#' @noRd
rnb.section.diffVar <- function(rnb.set,diff.meth,report,gzTable=FALSE,differentiality.method=rnb.getOption("differential.variability.method")){
  if (length(get.comparisons(diff.meth))<1){
    stop("no valid comparisons")
  }
  diff.var.cutoff <- 0.01
  logger.start("Adding Differential Variability Information")
  txt <- c("Differentially variable sites were computed with <code>", differentiality.method, "</code>.",
           " For more information about the method, have a look at the <a href=")
  if(differentiality.method=='diffVar'){
    txt <- c(txt,"https://bioconductor.org/packages/release/bioc/html/missMethyl.html>missMethyl</a> Bioconductor package.")
  }else if (differentiality.method=='iEVORA'){
    txt <- c(txt,"https://www.nature.com/articles/ncomms10478>iEVORA</a> method published in Nature.")
  }
  txt <- c(txt," This sections contains plots and tables describing the results of this test and further analyses ",
            "of the sites that were selected as differentially variable.")
  report <- rnb.add.section(report, 'Differential Variability', txt)
  
  #scatterplots
  logger.start("Adding scatterplots")
  rnb.cleanMem()
  plots <- list()
  if(parallel.isEnabled()){
    plots <- foreach(i=1:length(get.comparisons(diff.meth)),.combine="c") %dopar% {
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.variability.table(diff.meth,comp,return.data.frame=TRUE)
      comp.name.allowed <- ifelse(is.valid.fname(comp.name),comp.name,paste("cmp",i,sep=""))
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      res <- addReportPlot.diffVar.scatter(report,var.table,comp.name.allowed,
                                                      group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      res
    }
  } else {
    for (i in 1:length(get.comparisons(diff.meth))){
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.variability.table(diff.meth,comp,return.data.frame=TRUE)
      comp.name.allowed <- ifelse(is.valid.fname(comp.name),comp.name,paste("cmp",i,sep=""))
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      res <- addReportPlot.diffVar.scatter(report,var.table,comp.name.allowed,
                                           group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      plots <- c(plots,res)
    }
  }
  
  comps <- get.comparisons(diff.meth)
  diff.var.type = paste("FDR adjusted p-value &lt;",P.VAL.CUT)
  names(diff.var.type) = "fdrAdjPval"
  setting.names <- list(
    'comparison' = comps,
    'differential variability measure' = diff.var.type)
  description <- c('Scatterplot for differential variable sites.',
                   ' The transparency corresponds to point density. If the number of points exceeds ',DENS.SCATTER.SUBSAMPLE.THRES,
                   ' then the number of points for density estimation is reduced to that number by random sampling.',
                   ' The',round(DENS.SCATTER.SPARSE.POINTS.PERC*100),
                   '% of the points in the sparsest populated plot regions are drawn explicitly (up to a maximum of ',DENS.SCATTER.SPARSE.POINTS.MAX,
                   " points).",
                   ' Additionally, the colored points represent differentially variable sites (according to the selected criterion).')
  report <- rnb.add.figure(report, description, plots, setting.names)
  logger.completed()
  
  #volcano plots
  logger.start("Adding volcano plots")
  rnb.cleanMem()
  plots <- list()
  if(parallel.isEnabled()){
    plots <- foreach(i=1:length(get.comparisons(diff.meth)),.combine="c") %dopar% {
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.variability.table(diff.meth,comp,return.data.frame=TRUE)
      comp.name.allowed <- ifelse(is.valid.fname(comp.name),comp.name,paste("cmp",i,sep=""))
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      res <- addReportPlot.diffVar.volcano(report, var.table, comp.name.allowed,
                                                      group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      res
    }
  } else {
    for (i in 1:length(get.comparisons(diff.meth))){
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.variability.table(diff.meth,comp,return.data.frame=TRUE)
      comp.name.allowed <- comp.name
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      res <- addReportPlot.diffVar.volcano(report, var.table, comp.name.allowed,
                                           group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      plots <- c(plots,res)
    }
  }
  
  signif.measure <- c("pVal"="p-value","pValAdj"="adjusted p-value")
  setting.names <- list(
    'comparison' = comps,
    'significance metric' = signif.measure)
  description <- 'Volcano plot for differential variable sites.'
  report <- rnb.add.figure(report, description, plots, setting.names)
  logger.completed()
  
  #Comparing differentially methylated sites and differentially variable
  logger.start("Comparing DMCs and DMVs")
  rnb.cleanMem()
  plots <- list()
  if(parallel.isEnabled()){
    plots <- foreach(i=1:length(get.comparisons(diff.meth)),.combine="c") %dopar% {
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.variability.table(diff.meth,comp,return.data.frame=TRUE)
      meth.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
      comp.name.allowed <- ifelse(is.valid.fname(comp.name),comp.name,paste("cmp",i,sep=""))
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      res <- addReportPlot.diffVar.meth(report, var.table, meth.table, comp.name.allowed,
                                           group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      res
    }
  } else {
    for (i in 1:length(get.comparisons(diff.meth))){
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.variability.table(diff.meth,comp,return.data.frame=TRUE)
      meth.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
      comp.name.allowed <- ifelse(is.valid.fname(comp.name),comp.name,paste("cmp",i,sep=""))
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      res <- addReportPlot.diffVar.meth(report, var.table, meth.table, comp.name.allowed,
                                        group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      plots <- c(plots,res)
    }
  }
  setting.names <- list(
    'comparison' = comps)
  description <- 'Scatterplot comparing differentially methylated and variable sites.'
  report <- rnb.add.figure(report, description, plots, setting.names)
  logger.completed()
  
  
  #Tables
  logger.start("Adding tables")
  includeCovg <- hasCovg(rnb.set)
  
  txt <- c("Table describing the sites that were found to be differentially variable between two groups.
                 Below, a brief explanation of the different columns is shown:")
  report <- rnb.add.section(report, "Differential Variability Tables", txt, level = 2)
  
  column.description <- create.diffvar.column.description(includeCovg, covgThres=get.covg.thres(diff.meth))
  rnb.add.list(report, column.description,type="u")
  
  sectionText <- "The tables for the individual comparisons can be found here:\n<ul>\n"
  annot.cols <- c("Chromosome","Start","Strand")
  sites.info <- annotation(rnb.set,type="sites",add.names=FALSE)[, annot.cols]
  #add cg identifier for infinium datasets
  if (!inherits(rnb.set,"RnBiseqSet") && !is.null(rownames(sites.info))){
   sites.info <- data.frame(cgid=rownames(sites.info),sites.info,stringsAsFactors=FALSE)
  }
  grp.names <- get.comparison.grouplabels(diff.meth)
  for (i in 1:length(comps)){
    cc <- comps[i]
     
    annot.vec <- get.diffvar.tab.annot.cols(includeCovg)
    colname.vec <- get.diffvar.tab.annot.colnames.pretty(grp.names[i,1], grp.names[i,2], includeCovg, covgThres=get.covg.thres(diff.meth))
    var.table <- get.variability.table(diff.meth,cc,return.data.frame=TRUE)[,annot.vec]		
    colnames(var.table) <- colname.vec
    var.table <- cbind(rownames(var.table),sites.info,var.table)
    colnames(var.table)[1] <- "id"
    rownames(var.table) <- NULL
    ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
    fname <- paste("diffVarTable_",ccn,".csv",sep="")
    fname <- rnb.write.table(var.table,fname,fpath=rnb.get.directory(report, "data", absolute = TRUE),format="csv",gz=gzTable,row.names = FALSE,quote=FALSE)
    txt <- paste(c("<a href=\"", rnb.get.directory(report, "data"), "/", fname,"\">",cc,"</a>"),collapse="")
    sectionText <- paste(sectionText,"<li>",txt,"</li>\n",sep="")
  }
  sectionText <- paste(sectionText,"</ul>",sep="")
  rnb.add.paragraph(report, sectionText)
  
  logger.completed()
  logger.completed()
  return(report)
}
#' cols.to.rank.site
#' 
#' Return a matrix contianing the negative absolute values of the information used to rank the sites. Those are currently:
#' the variance difference, the log ratio in variances and the p-value from the statistical test.
#' 
#' @param diff.var A differential variability table.
#' @return A matrix with the absolute values of the relevant columns
#' @author Michael Scherer
#' @rdname cols.to.rank
#' @aliases cols.to.rank.site
#' @aliases cols.to.rank.region
cols.to.rank.site <- function(diff.var){
  return(cbind(-abs(diff.var$var.diff),-abs(diff.var$var.log.ratio),-abs(diff.var$diffVar.p.val)))
}

#' computeDiffVar.bin
#' 
#' This function computes differentially variable sites between the two groups specifified by the column indices on the given
#' methylation matrix for single sites.
#' @author Michael Scherer
#' @param meth.matrix Methylation Matrix with sites in the rows and samples as columns on which differential variability analysis
#'                    should be conducted
#' @param inds.g1 Sample indices in \code{rnbSet} of group 1 members
#' @param inds.g2 Sample indices in \code{rnbSet} of group 2 members
#' @param adjustment.table A \code{data.frame} containing variables to adjust for in the testing
#' @return A \code{data.frame} containing the information about the differential variability testing.
#' @noRd
computeDiffVar.bin <- function(meth.matrix,inds.g1,inds.g2,adjustment.table=NULL,...){
  if (length(union(inds.g1,inds.g2)) != (length(inds.g1)+length(inds.g2))){
    logger.error("Overlapping sample sets in differential variability analysis")
  }
  logger.start("Computing Differential Variability Table")
  diff.var <- computeDiffVar.extended.site(meth.matrix,inds.g1=inds.g1,inds.g2=inds.g2,...)
  to.rank <- cols.to.rank.site(diff.var = diff.var)
  ranks <- combinedRanking.tab(to.rank)
  diff.var <- data.frame(diff.var,combined.rank.var=ranks)
  logger.completed()
  return(diff.var)
  #Perhaps add permutation tests here
}

#' rnb.execute.diffVar
#' 
#' This routine computes sites that are differentially variable between two sample groups specified as the column name 
#' in the phenotypic table.
#' @author Michael Scherer
#' @param rnb.set Object of type \code{\linkS4class{RnBSet}} on which differential varibility analysis should be conducted
#' @param pheno.cols Column names used to define the classes, whose methylation variability should be compared with each other
#' @param columns.adj Column names or indices in the table of phenotypic information to be used for confounder adjustment in the
#'        differential variability analysis.
#' @param diff.meth Object of type \code{\linkS4class{RnBDiffMeth}}, possibly containing information about differentially methylated
#'        sites. The differentially variable sites then are also added to this object instead of creating a new object. 
#' @param adjust.celltype Flag indicating whether the resulting table should also contain estimated celltype contributions.
#' 				See \code{\link{rnb.execute.ct.estimation}} for details.
#' @param disk.dump Flag indicating whether the resulting differential methylation object should be file backed, ie.e the matrices dumped to disk
#' @param disk.dump.dir disk location for file backing of the resulting differential methylation object. Only meaningful if \code{disk.dump=TRUE}.
#' @return Object of type \code{\linkS4class{RnBDiffMeth}} containing information about the differential variability analysis.
#' @export
rnb.execute.diffVar <- function(rnb.set,pheno.cols,columns.adj=rnb.getOption("covariate.adjustment.columns"),
                                diff.meth=NULL,
                                adjust.celltype=rnb.getOption("differential.adjustment.celltype"),
                                disk.dump=rnb.getOption("disk.dump.big.matrices"),
                                disk.dump.dir=tempfile(pattern="diffVarTables_")){
  
  logger.start("Differential Variability")
  logger.start("Retrieving comparison info")
  cmp.info <- get.comparison.info(rnb.set, pheno.cols=pheno.cols, 
                                  columns.adj=columns.adj,
                                  adjust.celltype=adjust.celltype)
  logger.completed()
  if (is.null(cmp.info)) {
    return(NULL)
  }
  
  variability.method <- rnb.getOption("differential.variability.method")
  if(is.null(diff.meth)){
    diff.meth <- new("RnBDiffMeth",variability.method=variability.method,disk.dump=disk.dump,disk.path.DMVs=disk.dump.dir)
  }else{
    diff.meth <- prepare.diff.var(diff.meth,variability.method=variability.method,disk.dump=disk.dump,disk.path.DMVs=disk.dump.dir)
  }
  for(i in 1:length(cmp.info)){
    cmp.info.cur <- cmp.info[[i]]
    logger.start(c("Comparing ",cmp.info.cur$comparison))
    if(!isImputed(rnb.set)) rnb.set <- rnb.execute.imputation(rnb.set)
    meth.matrix <- meth(rnb.set)
    diffVar <- computeDiffVar.bin(meth.matrix=meth.matrix,inds.g1=cmp.info.cur$group.inds$group1,
                                  inds.g2=cmp.info.cur$group.inds$group2, adjustment.table=cmp.info.cur$adjustment.table,
                                  covg=covg(rnb.set))
    diff.meth <- addDiffVarTable(diff.meth,diffVar,cmp.info.cur$comparison,cmp.info.cur$group.names)
    logger.completed()
  }
  diff.meth <- addComparisonInfo(diff.meth,cmp.info)
  logger.completed()
  return(diff.meth)
}
