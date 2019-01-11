########################################################################################################################
## segmentation.R
## created: 2019-01-04
## creator: Michael Scherer
## ---------------------------------------------------------------------------------------------------------------------
## Functions for segmentation of the methylome through MethylSeekR
########################################################################################################################

#' rnb.execute.segmentation
#' 
#' This function computes methylation segmentation by MethylSeekR into PMDs, UMRs/LMRs, and HMRs. It is recommened to only
#' execute this function on WGBS data  (with coverage >=10 according to the developer's recommendation), but could also be
#' used with RRBS_HaeIII without guarantee and the results should be interpreted carefully.
#' 
#' @param rnb.set An object of type \code{\link{RnBiseqSet-class}} containing methylation and coverage information.
#' @param sample.name The sample for which segmentation is to be executed. Segemntation can only be exectued for each sample
#'                     individually.
#' @param meth.level Methylation cutoff to be used in UMR/LMR computation
#' @param fdr False discovery rate cutoff to be used in percent
#' @param min.cover The coverage threshold
#' @param n.cores The number of cores available for analysis
#' @param plot.path Location on disk on which diagnostic plots are to be stored. Defaults to the working directory.
#' @param temp.dir The temporary directory. Defaults to the R temporary directory.
#' @return The input RnBSet object with segementation added as an additional region type. Furthermore, three new annotations
#'          are set globally containing segmentation into PMDs, UMRs/LMRs, and HMRs for the sample that was specified.
#' @details For further descriptions on the methods, see \code{MethylSeekR}-documentation. The new annotations can be accessed
#'          via \code{rnb.get.annotation("[PMDs,UMRsLMRs,HMRs]_[sample.name]")}.
#' @author Michael Scherer, based on a script by Abdulrahman Salhab
#' @export
rnb.execute.segmentation <- function(rnb.set,
                                     sample.name,
                                     meth.level=0.5,
                                     fdr=5,
                                     min.cover=5,
                                     n.cores=1,
                                     plot.path=getwd(),
                                     temp.dir=tempdir()){
  if(!inherits(rnb.set,"RnBiseqSet")){
    logger.error("Invalid value for rnb.set, needs to be RnBiseqSet")
  }
  if(length(sample.name)>1){
    logger.error("Only single sample can be analysed")
  }
  rnb.require("MethylSeekR")
  asb <- assembly(rnb.set)
  if(asb %in% c("hg19","hg38")){
    rnb.require(paste0("BSgenome.Hsapiens.UCSC.",asb))
    myGenomeSeq <- Hsapiens
  }else if(asb %in% "mm10"){
    rnb.require(paste0("BSgenome.Mmusculus.UCSC.",asb))
    myGenomeSeq <- Mmusculus
  }else{
    logger.info("Invalid assembly; Segmentation only supported for hg19, hg38 and mm10")
  }

  sLengths <- seqlengths(rnb.get.annotation("CpG",assembly = assembly(rnb.set)))
  CpGislands.gr <- rnb.get.annotation("cpgislands",assembly=assembly(rnb.set))
  CpGislands.gr <-suppressWarnings(resize(CpGislands.gr, 5000, fix="center"))
  
  if(any(!(rnb.region.types.for.analysis(rnb.set) %in% rnb.region.types(assembly(rnb.set))))){
    rnb.load.annotation.from.db(rnb.region.types.for.analysis(rnb.set)[!(rnb.region.types.for.analysis(rnb.set) %in% rnb.region.types(assembly(rnb.set)))],assembly(rnb.set))
  }
  sel.sample <- sample.name %in% samples(rnb.set)
  if(!sel.sample){
    logger.error(paste("Specified sample",sample.name,"not in RnBSet"))
  }
  sel.sample <- which(sel.sample)
  meth.rnb <- rnb.set@meth.sites[,sel.sample]
  is.na.meth <- is.na(meth.rnb)
  covg.rnb <- rnb.set@covg.sites[,sel.sample]
  meth.matrix <- round(meth.rnb*covg.rnb)
  nas <- is.na(covg.rnb) | is.na.meth
  anno.rnb <- data.frame2GRanges(annotation(rnb.set),assembly=assembly(rnb.set))
  tmp.anno <- anno.rnb
  values(tmp.anno) <- data.frame(T=covg.rnb,M=meth.matrix)
  tmp.anno <- tmp.anno[!nas]
  tmp.file.meth <- file.path(temp.dir,"GRanges_for_MethylSeekR.RDS")
  saveRDS(tmp.anno,tmp.file.meth)
  data.gr <- readMethylome(tmp.file.meth, seqLengths=seqlengths(rnb.get.annotation("CpG",assembly(rnb.set))), format="GRanges")
  
  ### read SNPs
  cat("Reading the snp file ......\n")
  # The same thing as for the Methylome holds here
  anno.rnb <- anno.rnb[!is.na(values(anno.rnb)$SNPs)]
  values(anno.rnb) <- values(anno.rnb)$SNPs
  tmp.file.snps <- file.path(temp.dir,"GRangesSNP_for_MethylSeekR.RDS")
  saveRDS(anno.rnb,tmp.file.snps)
  snps.gr <- readSNPTable(tmp.file.snps, seqLengths=seqlengths(rnb.get.annotation("CpG",assembly(rnb.set))), format="GRanges")
  
  ### removing SNPs
  if(length(snps.gr)>0){
    cat("Filtering the SNPs ......\n")
    data.gr <- removeSNPs(data.gr, snps.gr)
  }
  
  ### find PMD and plot
  logger.start("Detecting PMDs")
  PMDsegments.gr <- segmentPMDs(m=data.gr,pdfFilename=file.path(plot.path,paste(sample.name,"alpha.model.fit.pdf",sep=".")), chr.sel="chr2",seqLengths=sLengths, num.cores=n.cores)
  logger.completed()
  
  logger.start("Plotting alpha distribution")
  plotAlphaDistributionOneChr(m=data.gr, chr.sel="chr2",pdfFilename=file.path(plot.path,paste(sample.name,"alpha_distribution.pdf",sep=".")),num.cores=n.cores)
  logger.completed()
  
  ### FDR calculation
  logger.start("FDR calculation")
  stats <- calculateFDRs(m=data.gr, CGIs=CpGislands.gr,PMDs=PMDsegments.gr, nCpG.cutoffs =seq(1, 17, by=3),pdfFilename=file.path(plot.path,paste(sample.name,"FDR.pdf",sep=".")),num.cores=n.cores)
  FDR.cutoff <- as.numeric(fdr)
  m.sel <- as.numeric(meth.level)
  n.sel <- as.integer(names(stats$FDRs[as.character(m.sel), ][stats$FDRs[as.character(m.sel), ]<FDR.cutoff])[1])
  logger.info(paste("Minimum number of CpGs in LMRs:",n.sel,"CpGs"))
  logger.completed()
  
  ### find UMR and LMR
  logger.start("Detecting UMRs and LMRs")
  UMRLMRsegments.gr <- segmentUMRsLMRs(m=data.gr, meth.cutoff=m.sel,nCpG.cutoff=n.sel, PMDs=PMDsegments.gr,pdfFilename=file.path(plot.path,paste(sample.name,"UMR.LMR.scatter.pdf",sep=".")),num.cores=n.cores, myGenomeSeq=myGenomeSeq,seqLengths=sLengths, minCover = min.cover)
  logger.completed()
  
  # create final segmentation
  tile.type <- rnb.region.types.for.analysis(rnb.set)[grep("tiling",rnb.region.types.for.analysis(rnb.set))]
  hmr.segments <- annotation(rnb.set,tile.type)
  hmr.segments <- GRanges(Rle(hmr.segments$Chromosome),IRanges(start = hmr.segments$Start,end=hmr.segments$End))
  op.pmds <- findOverlaps(hmr.segments,PMDsegments.gr[values(PMDsegments.gr)$type %in% "notPMD"])
  hmr.segments <- hmr.segments[queryHits(op.pmds)]
  op.umr.lmr <- findOverlaps(hmr.segments,UMRLMRsegments.gr)
  hmr.segments <- hmr.segments[-queryHits(op.umr.lmr)]
  values(hmr.segments)$HMR <- rep("HMR",length(hmr.segments))
  
  # set new annotations PMD, UMR/LMR, HMR
  pmd.frame <- data.frame(Chromosome=seqnames(PMDsegments.gr),Start=start(PMDsegments.gr),End=end(PMDsegments.gr),
                           PMD=values(PMDsegments.gr)$type)
  rnb.set.annotation(paste0("PMDs_",sample.name),regions=pmd.frame,description = "Partially Methylated Domains by MethylSeekR",assembly = asb)
  umr.lmr.frame <- data.frame(Chromosome=seqnames(UMRLMRsegments.gr),Start=start(UMRLMRsegments.gr),End=end(UMRLMRsegments.gr),
                          UMR_LMR=values(UMRLMRsegments.gr)$type)
  rnb.set.annotation(paste0("UMRsLMRs_",sample.name),regions=umr.lmr.frame,description = "Unmethylated and Lowly Methylated Regions by MethylSeekR",assembly = asb)
  hmr.frame <- data.frame(Chromosome=seqnames(hmr.segments),Start=start(hmr.segments),End=end(hmr.segments),
                              UMR_LMR=values(hmr.segments)$HMR)
  rnb.set.annotation(paste0("HMRs_",sample.name),regions=hmr.frame,description = "Highly Methylated Regions by MethylSeekR",assembly = asb)
  
  rnb.set <- summarize.regions(rnb.set,paste("PMDs",sample.name,sep="_"))
  rnb.set <- summarize.regions(rnb.set,paste("UMRsLMRs",sample.name,sep="_"))
  rnb.set <- summarize.regions(rnb.set,paste("HMRs",sample.name,sep="_"))
  
  unlink(tmp.file.meth)
  unlink(tmp.file.snps)
  return(rnb.set)
}

#' rnb.bed.from.segmentation
#' 
#' This function creates a BED file from the segmentation result of \code{rnb.execute.segmentation} and stores it on disk.
#' 
#' @param rnb.set An \code{\link{RnBSet-class}} object obtained by executing \code{rnb.execute.segmentation}.
#' @param sample.name The sample name for which segmentation was computed.
#' @param type The type of segmentation (\code{PMDs}, \code{UMRsLMRs} or \code{HMRs}).
#' @param store.path Path to which the BED file is to be stored.
#' @author Michael Scherer
#' @export
rnb.bed.from.segmentation <- function(rnb.set,
                                      sample.name,
                                      type="PMD",
                                      store.path=getwd()){
  if(!type %in% c("PMDs","UMRsLMRs","HMRs")){
    logger.error("Invalid value for type, needs to be PMDs, UMRsLMRs or HMRs.")
  }
  region.name <- paste(type,sample.name,sep = "_")
  if(!(region.name %in% summarized.regions(rnb.set))){
    logger.error("Segmentation not yet available, execute rnb.execute.segementation first")
  }
  bed.frame <- annotation(rnb.set,region.name)
  meth.seg <- meth(rnb.set)[,sample.name]
  bed.frame <- data.frame(bed.frame[,c("chromosome","start","end")],AvgMeth=meth.seg)
  write.table(bed.frame,file.path(store.path,paste0(sample.name,"_",type,".bed")),sep="\t",row.names=F)
}

#' rnb.plot.segmentation.final
#' 
#' This functions plots the final segmentation result.
#' 
#' @param data.gr The data object as \code{\link{GRanges}} object
#' @param UMR.LMR.segments.gr Segmentation into UMR/LMR as GRanges
#' @param PMD.segments Segmentation into PMDs/notPMDs as GRanges
#' @param n.regions Number of regions
#' @param meth.cutoff The methylation cutoff
#' @author Michael Scherer
#' @noRd
rnb.plot.segmentation.final <- function(data.gr,
                                        UMR.LMR.segments.gr,
                                        PMD.segments,
                                        n.regions=4,
                                        meth.cutoff){
  logger.start("Plotting final segmentation")
  plotFinalSegmentation(m=data.gr, segs=UMR.LMR.segments.gr,PMDs=PMD.segments,numRegions = n.regions,meth.cutoff=meth.cutoff)
  logger.completed()
}

#'rnb.plot.segmentation.distributions
#'
#'Plots the distributions of methylation and coverage.
#'
#'@param data.gr The data object as \code{\link{GRanges}} object
#'@author Michael Scherer
#'@noRd
rnb.plot.segmentation.distributions <- function(data.gr){
  df <- as.data.frame(data.gr)
  df$meth <- df$M/df$T
  meth.plot <- ggplot(df[df$T>=5,], aes(x=df[df$T>=5,"meth"])) + geom_density(colour="dodgerblue1",size=1) +  ylab("density") + xlab("beta value") + ggtitle("Methylation level density")
  covg.plot <- ggplot(df, aes(x=df$T)) + geom_histogram(binwidth=1,alpha=.5,position="identity",colour = "dodgerblue1", fill = "dodgerblue1") + geom_vline(xintercept=mean(df$T),colour="black", linetype = "longdash")  + ylab("Frequency") + xlab("read coverage per CpG") + geom_text(aes(x2,y2,label = texthere),data.frame(x2=mean(df$T), y2=max(table(df$T)), texthere=round(mean(df$T),2)))
  return((list(Methylation=meth.plot,Coverage=covg.plot)))
}