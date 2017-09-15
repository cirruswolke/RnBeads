########################################################################################################################
## differentialVaribility.R
## created: 2017-07-28
## creator: Michael Scherer
## ---------------------------------------------------------------------------------------------------------------------
## The differential variability analysis methods between sample groups.
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
      logger.info("No DVCs found.")
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

#' Returns the colum names of the differential variability table.
#' @param includeCovg Flag indicating if dataset contains coverage information
#' @return Column names of the differential variability table
get.table.ids <- function(includeCovg=FALSE){
  ret <- c("var.g1","var.g2","var.diff","var.log.ratio","diffVar.p.val")
  if(includeCovg){
    ret <- c(ret,"mean.covg.g1","covg.thresh.nsamples.g1","mean.covg.g2","covg.thresh.nsamples.g2")
  }
  return(ret)
}

#' addReportPlot.diffVar.scatter.site
#' 
#' This function creates a scatterplot for the given comparison comparing the difference in variances with the
#' usted) p-value as the result of the differentially variability analysis.
#' 
#' @param report Object of class \code{\linkS4class{Report}} to which the plot should be added
#' @param var.table Variability table containing the relevant information to be put into the scatterplot
#' @param comparison.name Comparison for which the scatterplot should be creaed
#' @param rank.cutoffs.numbers Numbers to be used for plotting the best ranked sites
#' @param auto.cutoff Auto selected cutoffs for the given comparison
#' @param group.name1 Name of the first group of the comparison
#' @param group.name2 Name of the second group of the comparison
#' @return list of scatterplots
#' @noRd
addReportPlot.diffVar.scatter.site <- function (report, var.table, comparison.name,
                              rank.cutoffs.numbers, auto.cutoff, rerank=TRUE,
                              group.name1="Group1",group.name2="Group2"){
  ret <- list()
  sparse.points <- DENS.SCATTER.SPARSE.POINTS.PERC
  if (DENS.SCATTER.SPARSE.POINTS.MAX < sparse.points*nrow(var.table)){
    sparse.points <- DENS.SCATTER.SPARSE.POINTS.MAX
  }
  dens.subsample <- FALSE
  if (nrow(var.table) > DENS.SCATTER.SUBSAMPLE.THRES){
    dens.subsample <- DENS.SCATTER.SUBSAMPLE.THRES
  }
  
  if("diffVar.p.adj.fdr" %in% colnames(var.table)){
    var.sites <- var.table[,"diffVar.p.adj.fdr"] < P.VAL.CUT
    plot <- create.densityScatter(var.table[,c("var.g1","var.g2")],is.special=var.sites,dens.subsample=dens.subsample) +
        xlab(paste("Variance",group.name1)) + ylab(paste("Variance",group.name2))
    comp.type <- "fdrAdjPval"
    fig.name <- paste("diffVar_site",comparison.name,comp.type,sep = "_")
    plot <- createReportGgPlot(plot,fig.name,report=report,create.pdf = FALSE, high.png = 200)
    plot <- off(plot,handle.errors=TRUE)
    ret <- c(ret,list(plot))
  }
  ranks <- var.table[,"combinedRank.var"]
  if (rerank)	ranks <- rank(ranks,na.last="keep",ties.method="min")
  for(i in 1:length(rank.cutoffs.numbers)){
    number <- rank.cutoffs.numbers[i]
    cutoff.name <- paste0("rc",i)
    var.table$isDVC <- ranks < number
    plot <- create.densityScatter(var.table[,c("var.g1","var.g2")], is.special = var.table$isDVC,dens.subsample=dens.subsample) +
      xlab(paste("Variance",group.name1)) + ylab(paste("Variance",group.name2))
    fig.name <- paste("diffVar_site",comparison.name,cutoff.name,sep="_")
    plot <- createReportGgPlot(plot,fig.name,report=report,create.pdf = FALSE, high.png = 200)
    plot <- off(plot,handle.errors=TRUE)
    ret <- c(ret,list(plot))
  }
  if(is.integer(auto.cutoff)){
    var.table$isDVC <- ranks < auto.cutoff
    plot <- create.densityScatter(var.table[,c("var.g1","var.g2")], is.special = var.table$isDVC,dens.subsample=dens.subsample) +
      xlab(paste("Variance",group.name1)) + ylab(paste("Variance",group.name2))
    fig.name <- paste("diffVar_site",comparison.name,"rcAuto",sep="_")
    plot <- createReportGgPlot(plot,fig.name,report=report,create.pdf = FALSE, high.png = 200)
    plot <- off(plot,handle.errors=TRUE)
    ret <- c(ret,list(plot))
  }
  return(ret)
}
#' @noRd
addReportPlot.diffVar.scatter.region <- function (report, var.table, comparison.name, region.name, rerank=TRUE,
                                           ranking.cutoffs, auto.cutoff=NULL, group.name1="Group1",group.name2="Group2"){
  ret <- list()
  
  al.x <- paste("Mean Variance",group.name1,sep=".")
  al.y <- paste("Mean Variance",group.name2,sep=".")
  
  dens.subsample <- FALSE
  if(nrow(var.table) > DENS.SCATTER.SUBSAMPLE.THRES){
    dens.subsample <- DENS.SCATTER.SUBSAMPLE.THRES
  }
  
  if("comb.p.adj.var.fdr" %in% colnames(var.table)){
    var.sites <- var.table[,"comb.p.adj.var.fdr"] < P.VAL.CUT
    plot <- create.densityScatter(var.table[,c("mean.var.g1","mean.var.g2")],is.special=var.sites,dens.subsample = dens.subsample) +
      xlab(al.x) + ylab(al.y)
    comp.type <- "fdrAdjPval"
    fig.name <- paste("diffVar_region",comparison.name,region.name,comp.type,sep = "_")
    plot <- createReportGgPlot(plot,fig.name,report=report,create.pdf = FALSE)
    plot <- off(plot,handle.errors=TRUE)
    ret <- c(ret,list(plot))
  }
  
  rrs <- var.table[,"combinedRank.var"]
  if (rerank)	rrs <- rank(rrs,na.last="keep",ties.method="min")
  for (i in 1:length(ranking.cutoffs)){
    rc <- ranking.cutoffs[i]
    cur.cut.name <- paste("rc",i,sep="")
    var.table$isDVR <- rrs < rc
    
    pp <- create.densityScatter(var.table[,c("mean.var.g1","mean.var.g2")],is.special=var.table$isDVR,
                                dens.subsample = dens.subsample, add.text.cor=TRUE) +
      labs(x=al.x, y=al.y) + coord_fixed()
    
    figName <- paste("diffVar_region",comparison.name,region.name,cur.cut.name,sep="_")
    report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
    report.plot <- off(report.plot,handle.errors=TRUE)
    ret <- c(ret,list(report.plot))
  }
  
  if (is.integer(auto.cutoff)){
    var.table$isDVR <- var.table[,"combinedRank.var"] <= auto.cutoff
    pp <- create.densityScatter(var.table[,c("mean.var.g1","mean.var.g2")],is.special=var.table$isDVR,
                                dens.subsample = dens.subsample, add.text.cor=TRUE) +
      labs(x=al.x, y=al.y) + coord_fixed()
    figName <- paste("diffVar_region",comparison.name,region.name,"rcAuto",sep="_")
    report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
    report.plot <- off(report.plot,handle.errors=TRUE)
    ret <- c(ret,list(report.plot))
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
#' @aliases addReportPlot.diffVar.volcano
#' @aliases addReportPlot.diffVar.volcano.region
#' @return list of volcano plots
#' @noRd
addReportPlot.diffVar.volcano <- function (report, var.table, comparison.name,
                                           group.name1="Group1",group.name2="Group2"){
  ret <- list()
  dont.plot.p.val <- all(is.na(var.table[,"diffVar.p.val"]))
  
  fig.name <- paste("diffVar_volcano_site",comparison.name,"pVal",sep="_")
  if (!dont.plot.p.val){
    pp <- ggplot(var.table) + aes(x=var.diff,y=log10P,color=log10(combinedRank.var)) + 
      scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
      geom_point()
  } else {
    pp <- rnb.message.plot("No p-value available")
  }
  report.plot <- createReportGgPlot(pp,fig.name, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  ret <- c(ret,list(report.plot))
  
  fig.name <- paste("diffVar_volcano_site",comparison.name,"pValAdj",sep="_")
  pp <- ggplot(var.table) + aes(x=var.diff,y=log10FDR,color=log10(combinedRank.var)) +
    scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient")))+
    geom_point()
  report.plot <- createReportGgPlot(pp,fig.name, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  ret <- c(ret,list(report.plot))
  
  return(ret)
}
#' addReportPlot.diffVar.volcano.region
#' @noRd
addReportPlot.diffVar.volcano.region <- function(report, var.table,comparison.name,region.type,
                                                  group.name1,group.name2){
  cn.d <- "mean.var.diff"
  cn.q <- "mean.var.log.ratio"
  cn.p <- "comb.p.val.var"
  cn.pa <- "comb.p.adj.var.fdr"

  figPlots <- list()
  dont.plot.p.val <- all(is.na(var.table[,cn.p]))
  
  figName <- paste("diffVar_region_volcano",comparison.name,region.type,"diff","pVal",sep="_")
  if (!dont.plot.p.val){
    pp <- ggplot(var.table) + aes_string(cn.d, paste0("-log10(",cn.p,")"), color="log10(combinedRank.var)") +
      scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
      geom_point(aes(order=plyr::desc(rank(combinedRank.var,ties.method="first",na.last=TRUE))))#(alpha=0.3)
  } else {
    pp <- rnb.message.plot("No p-value available")
  }
  report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  figPlots <- c(figPlots,list(report.plot))
  
  figName <- paste("diffVar_region_volcano",comparison.name,region.type,"diff","pValAdj",sep="_")
  pp <- ggplot(var.table) + aes_string(cn.d, paste0("-log10(",cn.pa,")"), color="log10(combinedRank.var)") +
    scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
    geom_point(aes(order=plyr::desc(rank(combinedRank.var,ties.method="first",na.last=TRUE))))#(alpha=0.3)
  report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  figPlots <- c(figPlots,list(report.plot))
  
  figName <- paste("diffVar_region_volcano",comparison.name,region.type,"quot","pVal",sep="_")
  if (!dont.plot.p.val){
    pp <- ggplot(var.table) + aes_string(cn.q, paste0("-log10(",cn.p,")"), color="log10(combinedRank.var)") +
      scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
      geom_point(aes(order=plyr::desc(rank(combinedRank.var,ties.method="first",na.last=TRUE))))#(alpha=0.3)
  } else {
    pp <- rnb.message.plot("No p-value available")
  }
  report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  figPlots <- c(figPlots,list(report.plot))
  
  figName <- paste("diffVar_region_volcano",comparison.name,region.type,"quot","pValAdj",sep="_")
  pp <- ggplot(var.table) + aes_string(cn.q, paste0("-log10(",cn.pa,")"), color="log10(combinedRank.var)") +
    scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
    geom_point(aes(order=plyr::desc(rank(combinedRank.var,ties.method="first",na.last=TRUE))))#(alpha=0.3)
  report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  figPlots <- c(figPlots,list(report.plot))
  
  figName <- paste("diffVar_region_volcano",comparison.name,region.type,"diff","quotSig",sep="_")
  pp <- ggplot(var.table) + aes_string(cn.d, cn.q, color="log10(combinedRank.var)") +
    scale_color_gradientn(colours=rev(rnb.getOption("colors.gradient"))) +
    geom_point(aes(order=plyr::desc(rank(combinedRank.var,ties.method="first",na.last=TRUE))))#(alpha=0.3)
  report.plot <- createReportGgPlot(pp, figName, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  figPlots <- c(figPlots,list(report.plot))
  
  figName <- paste("diffVar_region_volcano",comparison.name,region.type,"quot","quotSig",sep="_")
  pp <- rnb.message.plot("Quotient--Quotient scatterplot not available")
  report.plot <- createReportGgPlot(pp,figName, report,create.pdf=FALSE,high.png=200)
  report.plot <- off(report.plot,handle.errors=TRUE)
  figPlots <- c(figPlots,list(report.plot))
  
  return(figPlots)
}
#' addReportPlot.diffVar.meth
#' 
#' This function compares differentially methylated sites (DMCs) with differentially variable sites (DVCs) in a scatterplot
#' and colors significantly differential sites with different colors.
#' 
#' @param report Object of class \code{\linkS4class{Report}} to which the plot should be added.
#' @param var.table Table containing the differentially variable sites.
#' @param meth.table Table containing the differentially methylated sites.
#' @param comparison.name Name of the comparison to be conducted.
#' @param group.name1 Name of the first group in the comparison
#' @param group.name2 Name of the second group in the comparison
#' @return The modified report
#' @noRd
addReportPlot.diffVar.meth <- function(report, var.table, comparison.name,
                           group.name1,group.name2){
  ret <- list()
  dont.plot.p.val <- all(is.na(var.table[,"diffVar.p.adj.fdr"]))
  dont.plot.p.val <- dont.plot.p.val && all(is.na(meth.table[,"diffmeth.p.adj.fdr"]))
  
  fig.name <- paste("diffVarMeth_site",comparison.name,sep="_")
  if (!dont.plot.p.val){
    toPlot <- data.frame(Meth.p=var.table[,"diffmeth.p.adj.fdr"],Var.p=var.table[,"diffVar.p.adj.fdr"])
    is.diff.meth <- toPlot$Meth.p < P.VAL.CUT
    is.diff.var <- toPlot$Var.p < P.VAL.CUT
    color.diff <- rep('Not Significant',dim(toPlot)[1])
    color.diff[is.diff.meth] <- "DMR"
    color.diff[is.diff.var] <- "DVC"
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
  rank.cutoffs.numbers <- c(100,1000,10000,100000)
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
  
  diffRegionRankCut <- c(100,500,1000)
  comps <- get.comparisons(diff.meth)

  logger.start("Selection of rank cutoffs")
  rank.cuts.auto <- lapply(1:length(comps),FUN=function(i){
    comp.name <- names(comps)[i]
    comp <- comps[[i]]
    var.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
    res <- auto.select.rank.cut(var.table$diffVar.p.adj.fdr,var.table$combinedRank.var,alpha=0.1)
    return(as.integer(res))
  })
  txt <- paste("The following rank cutoffs have been automatically selected for the analysis of differentially",
               "variable sites:")
  rnb.add.paragraph(report, txt)
  
  tt <- data.frame(matrix(unlist(rank.cuts.auto),ncol=1,nrow=length(comps),byrow=TRUE))
  colnames(tt) <- "Rank Cutoff"
  rownames(tt) <- comps
  rnb.add.table(report,tt)
  logger.completed()
  
  #scatterplots
  logger.start("Adding scatterplots")
  rnb.cleanMem()
  plots <- list()
  if(parallel.isEnabled()){
    plots <- foreach(i=1:length(get.comparisons(diff.meth)),.combine="c") %dopar% {
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
      comp.name.allowed <- ifelse(is.valid.fname(comp.name),comp.name,paste("cmp",i,sep=""))
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      selected.cutoff <- rank.cuts.auto[[i]]
      res <- addReportPlot.diffVar.scatter.site(report,var.table,comp.name.allowed,
                                                      rank.cutoffs.numbers=rank.cutoffs.numbers,
                                                      auto.cutoff=selected.cutoff, rerank = TRUE,
                                                      group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      res
    }
  } else {
    for (i in 1:length(get.comparisons(diff.meth))){
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
      comp.name.allowed <- ifelse(is.valid.fname(comp.name),comp.name,paste("cmp",i,sep=""))
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      selected.cutoff <- rank.cuts.auto[[i]]
      res <- addReportPlot.diffVar.scatter.site(report,var.table,comp.name.allowed,
                                                rank.cutoffs.numbers=rank.cutoffs.numbers,
                                                auto.cutoff=selected.cutoff,
                                           group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      plots <- c(plots,res)
    }
  }
  
  comps <- get.comparisons(diff.meth)
  diff.var.type <- c(paste("FDR adjusted p-value &lt;",P.VAL.CUT),
                     paste("combined rank among the",rank.cutoffs.numbers,"best ranking sites"),
                     "automatically selected rank cutoff")
  names(diff.var.type) <-  c("fdrAdjPval",paste0("rc",1:length(rank.cutoffs.numbers)),"rcAuto")
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
      var.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
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
      var.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
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
  logger.start("Comparing DMCs and DVCs")
  rnb.cleanMem()
  plots <- list()
  if(parallel.isEnabled()){
    plots <- foreach(i=1:length(get.comparisons(diff.meth)),.combine="c") %dopar% {
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
      comp.name.allowed <- ifelse(is.valid.fname(comp.name),comp.name,paste("cmp",i,sep=""))
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      res <- addReportPlot.diffVar.meth(report, var.table, comp.name.allowed,
                                           group.name1=grp.names[1],group.name2=grp.names[2])
      rnb.cleanMem()
      res
    }
  } else {
    for (i in 1:length(get.comparisons(diff.meth))){
      comp.name <- names(get.comparisons(diff.meth))[i]
      comp <- get.comparisons(diff.meth)[comp.name]
      var.table <- get.table(diff.meth,comp,region.type="sites",return.data.frame=TRUE)
      comp.name.allowed <- ifelse(is.valid.fname(comp.name),comp.name,paste("cmp",i,sep=""))
      grp.names <- get.comparison.grouplabels(diff.meth)[comp,]
      res <- addReportPlot.diffVar.meth(report, var.table, comp.name.allowed,
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
  
  logger.completed()
  return(report)
}
#' rnb.section.diffVar.region
#'
#' Adds information for differentially variable regions to the report
#' @author Michael Scherer
#' @aliases rnb.section.diffMeth.region
#' @param rnb.set Object of type \code{\linkS4class{RnBSet}} containing methylation information
#' @param diff.meth RnBDiffMeth object. See \code{\link{RnBDiffMeth-class}} for details.
#' @param report Report object to which the content is added
#' @param gzTable Flag indicating if tables should be gzipped
#' @return The modified report object
rnb.section.diffVar.region <- function(rnb.set,diff.meth,report,gzTable=FALSE){
  if (length(get.comparisons(diff.meth))<1){
    stop("no valid comparisons")
  }
  if (length(get.region.types(diff.meth))<1){
    stop("no valid region types")
  }

  diffRegionRankCut <- c(100,500,1000)
  logger.start("Adding Region Level Information")

  sectionText <- c("Differential variability on the region level was computed based on a variety of metrics. ", 
                   "Those metrics were similar to the ones used for differential variability, but the mean of variances, ",
                   "the log-ratio of the quotient of variances as well as the p-values from the differentiality test were ",
                   "employed. Ranking was performed in line with the ranking of differential methylation.")
  report <- rnb.add.section(report, "Region Level", sectionText)
  
  comps <- get.comparisons(diff.meth)
  reg.types <- get.region.types(diff.meth)
  
  logger.start("Selection of rank cutoffs")
  rank.cuts.auto <- lapply(1:length(comps),FUN=function(i){
    lapply(1:length(reg.types),FUN=function(j){
      var.table <- get.table(diff.meth,comps[i],region.type=reg.types[j],return.data.frame=TRUE)
      res <- auto.select.rank.cut(var.table$comb.p.adj.var.fdr,var.table$combinedRank.var,alpha=0.1)
      return(as.integer(res))
    })
  })
  txt <- paste("The following rank cutoffs have been automatically selected for the analysis of differentially",
               "variable regions:")
  rnb.add.paragraph(report, txt)
  
  tt <- data.frame(matrix(unlist(rank.cuts.auto),ncol=length(reg.types),nrow=length(comps),byrow=TRUE))
  colnames(tt) <- reg.types
  rownames(tt) <- comps
  rnb.add.table(report,tt)
  logger.completed()
  
  #scatterplots
  logger.start("Adding scatterplots")
  rnb.cleanMem()
  addedPlots <- list()
  grp.labels <- get.comparison.grouplabels(diff.meth)
  if(parallel.isEnabled()){
    iis <- 1:length(comps)
    jjs <- 1:length(reg.types)
    pps <- expand.grid(iis,jjs)
    
    addedPlots <- foreach(k=1:nrow(pps),.combine="c") %dopar% {
      i <- pps[k,1]
      j <- pps[k,2]
      cc <- names(comps)[i]
      ccc <- comps[cc]
      ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
      rr <- reg.types[j]
      rrn <- ifelse(is.valid.fname(rr),rr,paste("reg",j,sep=""))
      auto.rank.cut <- rank.cuts.auto[[i]][[j]]
      var.table <- get.table(diff.meth,ccc,rr,return.data.frame=TRUE)
      res <- addReportPlot.diffVar.scatter.region(report,var.table,comparison.name=ccn,region.name=rrn,
                                                  ranking.cutoffs = diffRegionRankCut, rerank=TRUE,
                                                  auto.cutoff = auto.rank.cut,group.name1 = grp.labels[ccc,1],group.name2 = grp.labels[ccc,2])
      rnb.cleanMem()
      res
    }
    
  } else {
    for (i in 1:length(comps)){
      cc <- names(comps)[i]
      ccc <- comps[cc]
      ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
      for (j in 1:length(reg.types)){
        rr <- reg.types[j]
        rrn <- ifelse(is.valid.fname(rr),rr,paste("reg",j,sep=""))
        auto.rank.cut <- rank.cuts.auto[[i]][[j]]
        var.table <- get.table(diff.meth,ccc,rr,return.data.frame=TRUE)
        addedPlots <- c(addedPlots,addReportPlot.diffVar.scatter.region(report,var.table,comparison.name=ccn,region.name=rrn,
                                                    ranking.cutoffs = diffRegionRankCut,
                                                    auto.cutoff = auto.rank.cut,group.name1 = grp.labels[ccc,1],group.name2 = grp.labels[ccc,2])
        )
       rnb.cleanMem()
      }
    }
  }
  
  names(reg.types) <- ifelse(is.valid.fname(reg.types),reg.types,paste("reg",1:length(reg.types),sep=""))
  diffVarType = c(paste("FDR adjusted p-value <",P.VAL.CUT),
                   paste("combined rank among the ",diffRegionRankCut," best ranking regions",sep=""),
                   "automatically selected rank cutoff")
  names(diffVarType) = c("fdrAdjPval",paste("rc",1:length(diffRegionRankCut),sep=""),"rcAuto")
  setting.names <- list(
    'comparison' = comps,
    'regions' = reg.types ,
    'differential variability measure' = diffVarType)
  description <- 'Scatterplot for differential variable regions. The transparency corresponds to point density. The 1% of the points in the sparsest populated plot regions are drawn explicitly.
  Additionally, the colored points represent differentially methylated regions (according to the selected criterion).'
  report <- rnb.add.figure(report, description, addedPlots, setting.names)
  logger.completed()
  
  #volcano plots
  logger.start("Adding volcano plots")
  rnb.cleanMem()
  addedPlots <- list()
  if(parallel.isEnabled()){
    #generate pairs of comparison, region combinations with indices
    iis <- 1:length(comps)
    jjs <- 1:length(reg.types)
    pps <- expand.grid(iis,jjs)

    addedPlots <- foreach(k=1:nrow(pps),.combine="c") %dopar% {
      i <- pps[k,1]
      j <- pps[k,2]
      cc <- names(comps)[i]
      ccc <- comps[cc]
      ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
      rr <- reg.types[j]
      rrn <- ifelse(is.valid.fname(rr),rr,paste("reg",j,sep=""))
      var.table <- get.table(diff.meth,ccc,region.type=rr,return.data.frame=TRUE)
      res <- addReportPlot.diffVar.volcano.region(report,var.table,comparison.name=ccn,region.type=rrn,
                                                        group.name1=grp.labels[ccc,1],group.name2=grp.labels[ccc,2])
      rnb.cleanMem()
      res
    }
  } else {
    for (i in 1:length(comps)){
      cc <- names(comps)[i]
      ccc <- comps[cc]
      ccn <- ifelse(is.valid.fname(cc),cc,paste("cmp",i,sep=""))
      for (j in 1:length(reg.types)){
        rr <- reg.types[j]
        rrn <- ifelse(is.valid.fname(rr),rr,paste("reg",j,sep=""))
      var.table <- get.table(diff.meth,ccc,region.type=rr,return.data.frame=TRUE)
      addedPlots <- c(addedPlots,addReportPlot.diffVar.volcano.region(report,var.table,comparison.name=ccn,region.type=rrn,
                                                        group.name1=grp.labels[ccc,1],group.name2=grp.labels[ccc,2])
      )
      rnb.cleanMem()
      }
    }
  }
  diff.measure <- c("diff"="Difference","quot"="Quotient")
  signif.measure <- c("pVal"="combined p-value","pValAdj"="adjusted combined p-value","quotSig"="Quotient (only meaningful if 'Difference' is selected above)")
  setting.names <- list(
    'comparison' = comps,
    'regions' = reg.types ,
    'difference metric' = diff.measure,
    'significance metric' = signif.measure)
  description <- 'Volcano plot for differential variability quantified by various metrics. Color scale according to
  combined ranking.'
  report <- rnb.add.figure(report, description, addedPlots, setting.names)
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
  return(cbind(-abs(diff.var$var.diff), -abs(diff.var$var.log.ratio), diff.var$diffVar.p.val))
}
#' @rdname cols.to.rank
cols.to.rank.region <- function(diff.var){
  return(cbind(-abs(diff.var$mean.var.diff), -abs(diff.var$mean.var.log.ratio), diff.var$comb.p.adj.var.fdr))
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
rnb.execute.diffVar <- function(rnb.set,pheno.cols=rnb.getOption("differential.comparison.columns"),
                                region.types=rnb.region.types.for.analysis(rnb.set),
                                columns.adj=rnb.getOption("covariate.adjustment.columns"),
                                diff.meth=NULL,
                                adjust.celltype=rnb.getOption("differential.adjustment.celltype"),
                                disk.dump=rnb.getOption("disk.dump.big.matrices"),
                                disk.dump.dir=tempfile(pattern="diffVarTables_")){
  
  logger.start("Differential Variability")
  logger.start("Retrieving comparison info")
  cmp.info <- get.comparison.info(rnb.set, pheno.cols=pheno.cols,region.types = region.types, 
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
    diffVar <- computeDiffVar.bin.site(meth.matrix=meth.matrix,inds.g1=cmp.info.cur$group.inds$group1,
                                  inds.g2=cmp.info.cur$group.inds$group2, adjustment.table=cmp.info.cur$adjustment.table,
                                  covg=covg(rnb.set))
    diff.meth <- addDiffVarTable(diff.meth,diffVar,comparison=cmp.info.cur$comparison,grp.labs=cmp.info.cur$group.names)
    rnb.cleanMem()
    if (length(cmp.info.cur$region.types)>0){
      diffVar.region <- computeDiffVar.bin.region(rnb.set,diffVar,
                              inds.g1=cmp.info.cur$group.inds$group1,inds.g2=cmp.info.cur$group.inds$group2,
                              region.types=cmp.info.cur$region.types
        )
      for (rt in cmp.info.cur$region.types){
        diff.meth <- addDiffVarTable(diff.meth,diffVar.region[[rt]],comparison=cmp.info.cur$comparison, 
                                     region.type=rt, grp.labs=cmp.info.cur$group.names
        )
      }
    }
    logger.completed()
  }
  if(!all(get.comparisons(diff.meth)%in%names(cmp.info)))
    diff.meth <- addComparisonInfo(diff.meth,cmp.info)
  logger.completed()
  return(diff.meth)
}