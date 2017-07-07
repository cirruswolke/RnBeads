########################################################################################################################
## lolaUtils.R
## created: 2017-07-06
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Methods for LOLA enrichment analysis
########################################################################################################################

#' getCellTypesFromLolaDb
#'
#' retrieve or guess cell types from a LOLA DB object
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \link{\code{loadLolaDbs}}
#' @return character vector with cell types
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' lolaDb <- loadLolaDbs(lolaDirs[["hg19"]])
#' getCellTypesFromLolaDb(lolaDb)
#' }
getCellTypesFromLolaDb <- function(lolaDb){
	tt <- data.frame(lolaDb$regionAnno)
	res <- tt$cellType

	# special treatment for "sheffield_dnase" from LOLACore
	isInCollection <- tt$collection=="sheffield_dnase"
	if (sum(isInCollection) > 0) res[isInCollection] <- tt$description[isInCollection]

	# special treatment for "roadmap_epigenomics" from LOLAExt
	isInCollection <- tt$collection=="roadmap_epigenomics"
	if (sum(isInCollection) > 0){
		res[isInCollection] <- gsub("^(.+)-(.+)$", "\\1", tt$filename[isInCollection])
	}

	return(res)
}

#' getTargetFromLolaDb
#'
#' retrieve or guess the target from a LOLA DB object. Here, target typically
#' refers to antibodies for ChIP-seq experiments, but could also refer to other annotations
#' (e.g. motifs in TF motif databases, annotation according to UCSC features etc.)
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \link{\code{loadLolaDbs}}
#' @return character vector with targets
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' lolaDb <- loadLolaDbs(lolaDirs[["hg19"]])
#' getTargetFromLolaDb(lolaDb)
#' }
getTargetFromLolaDb <- function(lolaDb){
	tt <- data.frame(lolaDb$regionAnno)
	res <- tt$antibody

	# special treatment for "ucsc_features" from LOLACore
	isInCollection <- tt$collection=="ucsc_features"
	if (sum(isInCollection) > 0) res[isInCollection] <- gsub("^UCSC ", "", tt$description[isInCollection])

	# special treatment for "sheffield_dnase" from LOLACore
	isInCollection <- tt$collection=="sheffield_dnase"
	if (sum(isInCollection) > 0) res[isInCollection] <- "DNase"

	# special treatment for "encode_segmentation" from LOLACore
	isInCollection <- tt$collection=="encode_segmentation"
	if (sum(isInCollection) > 0){
		res[isInCollection] <- gsub(" Segments$", "", tt$description[isInCollection])
	}

	# special treatment for "jaspar_motifs" from LOLAExt
	isInCollection <- tt$collection=="jaspar_motifs"
	if (sum(isInCollection) > 0) res[isInCollection] <- gsub("\\.bed$", "", tt$filename[isInCollection])

	# special treatment for "roadmap_epigenomics" from LOLAExt
	isInCollection <- tt$collection=="roadmap_epigenomics"
	if (sum(isInCollection) > 0){
		res[isInCollection] <- gsub("(\\.(hotspot\\.(fdr0\\.01|all)|macs2))?\\.(peaks\\.bed|narrowPeak)$", "", 
			gsub("^(.+)-(.+)$", "\\2", tt$filename[isInCollection])
		)
	}

	return(res)
}

#' getNamesFromLolaDb
#'
#' get human readable names from a LOLA DB object
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \link{\code{loadLolaDbs}}
#' @param addCollectionNames  attach the name of the collection to the name
#' @param addDbId  attach the index of the item in the LOLA DB object to the name
#' @return character vector with human readable names
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' lolaDb <- loadLolaDbs(lolaDirs[["hg19"]])
#' getNamesFromLolaDb(lolaDb)
#' }
getNamesFromLolaDb <- function(lolaDb, addCollectionNames=FALSE, addDbId=TRUE){
	tt <- data.frame(lolaDb$regionAnno)
	res <- paste0(tt$description)

	cellTypes <- getCellTypesFromLolaDb(lolaDb)
	targets <- getTargetFromLolaDb(lolaDb)

	# if there is a valid "target" (see getTargetFromLolaDb) inferred from the annotation, construct a name from it
	hasDesc <- !is.na(targets)
	if (sum(hasDesc) > 0) res[hasDesc] <- targets[hasDesc]

	# if both cellType and target annotation were inferred, construct a name from them
	hasDesc <- !is.na(cellTypes) & !is.na(targets)
	if (sum(hasDesc) > 0) res[hasDesc] <- paste0(targets[hasDesc], " - ", cellTypes[hasDesc])

	# add suffixes if requested
	if (addCollectionNames){
		res <- paste0(res, " [", tt$collection, "]")
	}
	if (addDbId){
		res <- paste0(res, " [db", 1:nrow(tt), "]")
	}
	return(res)
}

#' lolaPrepareDataFrameForPlot
#'
#' Helper function for preparing a data.frame for multiple LOLA plots
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \link{\code{loadLolaDbs}}
#' @param lolaRes  LOLA enrichment result as returned by the \code{runLOLA} function from the \code{LOLA} package
#' @param scoreCol column name in \code{lolaRes} to be plotted
#' @param orderCol column name in \code{lolaRes} which is used for sorting the results
#' @param includedCollections vector of collection names to be included in the plot. If empty (default), all collections are used
#' @param pvalCut.neglog p-value cutoff (natural logarithm) to be employed for filtering the results
#' @param maxTerms maximum number of items to be included in the plot
#' @param colorpanel colors to be used for coloring the bars according to "target" (see \link{\code{getTargetFromLolaDb}}). An empty
#'                 vector indicates that black will be used for all bars.
#' @param groupByCollection facet the plot by collection
#' @param change.pval.base change p-value to base 10 (instead of base e)
#' @param orderDecreasing flag indicating whether the value in \code{orderCol} should be considered as decreasing (as opposed
#'                 to increasing). \code{NULL} (default) for automatic determination.
#' @return ggplot object containing the plot
#'
#' @author Fabian Mueller
#' @noRd
lolaPrepareDataFrameForPlot <- function(lolaDb, lolaRes, scoreCol="pValueLog", orderCol=scoreCol, includedCollections=c(), pvalCut.neglog=-log(0.01), maxTerms=50, perUserSet=FALSE, groupByCollection=TRUE, change.pval.base=TRUE, orderDecreasing=NULL){
	#dedect by column name whether decreasing order needs to be used
	if (is.null(orderDecreasing)){
		oset.dec <- c("pValueLog", "logOddsRatio")
		oset.inc <- c("maxRnk", "meanRnk")
		if (!is.element(orderCol, c(oset.dec, oset.inc))){
			logger.error(c("Could not determine whether to use increasing or decreasing order for column:", orderCol))
		}
		orderDecreasing <- is.element(orderCol, oset.dec)
	}

	lolaRes <- as.data.frame(lolaRes)
	lolaRes$name <- getNamesFromLolaDb(lolaDb, addCollectionNames=!groupByCollection)[lolaRes$dbSet]
	lolaRes$name <- trimChar(lolaRes$name, len.out=50, trim.str="...", len.pref=30)
	lolaRes$target <- getTargetFromLolaDb(lolaDb)[lolaRes$dbSet]

	if (!is.element(scoreCol, colnames(lolaRes))){
		logger.error(c("Missing score column from LOLA result:", scoreCol))
	}
	if (!is.element(orderCol, colnames(lolaRes))){
		logger.error(c("Missing order column from LOLA result:", orderCol))
	}

	if (length(includedCollections) > 0){
		lolaRes <- lolaRes[lolaRes[["collection"]] %in% includedCollections,]
	}

	if (change.pval.base) {
		# change the base from the natural logarithm to log10
		lolaRes$pValueLog <- lolaRes$pValueLog / log(10)
		pvalCut.neglog <- pvalCut.neglog / log(10)
	}
	calledSignif <- lolaRes$pValueLog > pvalCut.neglog
	lolaRes$pvalLab <- format(exp(-lolaRes$pValueLog), digits=2)

	if (sum(calledSignif) < 1){
		lolaRes.signif <- lolaRes
		logger.warning("No significant association found")
		return(rnb.message.plot("No significant association found"))
	} else {
		lolaRes.signif <- lolaRes[calledSignif,]
	}
	lolaRes.signif <- lolaRes.signif[order(lolaRes.signif[,orderCol], decreasing=orderDecreasing),]

	# data frame for plotting sorted according to significance/odds ratio
	df2p <- lolaRes.signif
	dbSets2keep <- unique(df2p[["dbSet"]])
	if (length(dbSets2keep) > maxTerms){
		if (perUserSet){
			logger.info(c("Reduced the number of region sets in the plot to", maxTerms, "per userSet"))
			uSetF <- factor(lolaRes[["userSet"]])
			dbSets2keep <- c()
			for (ss in levels(uSetF)){
				dbSets.cur <- unique(df2p[uSetF==ss,"dbSet"])
				dbSets2keep <- union(dbSets2keep, dbSets.cur[1:min(maxTerms, length(dbSets.cur))])
			}
		} else {
			logger.info(c("Reduced the number of region sets in the plot to", maxTerms))
			dbSets2keep <- dbSets2keep[1:maxTerms]
		}
		df2p <- df2p[df2p[["dbSet"]] %in% dbSets2keep, ]
	}

	df2p[["term"]] <- factor(df2p[["name"]], levels=unique(df2p[["name"]])) #sort factor by occurence (data.frame is ordered)
	df2p[["target"]] <- factor(df2p[["target"]])

	#bound score at maximum (excluding +/-Inf)
	scoresFinite <- df2p[[scoreCol]][is.finite(df2p[[scoreCol]])]
	cutVal <- max(scoresFinite, na.rm=TRUE)
	df2p[[scoreCol]][df2p[[scoreCol]]>cutVal] <- cutVal
	cutVal <- min(scoresFinite, na.rm=TRUE)
	df2p[[scoreCol]][df2p[[scoreCol]]<cutVal] <- cutVal

	return(df2p)
}

#' lolaBarPlot
#'
#' plot a barplot of LOLA enrichment results
#'
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \link{\code{loadLolaDbs}}
#' @param lolaRes  LOLA enrichment result as returned by the \code{runLOLA} function from the \code{LOLA} package
#' @param scoreCol column name in \code{lolaRes} to be plotted
#' @param orderCol column name in \code{lolaRes} which is used for sorting the results
#' @param includedCollections vector of collection names to be included in the plot. If empty (default), all collections are used
#' @param pvalCut.neglog p-value cutoff (natural logarithm) to be employed for filtering the results
#' @param maxTerms maximum number of items to be included in the plot
#' @param colorpanel colors to be used for coloring the bars according to "target" (see \link{\code{getTargetFromLolaDb}}). An empty
#'                 vector indicates that black will be used for all bars.
#' @param groupByCollection facet the plot by collection
#' @param change.pval.base change p-value to base 10 (instead of base e)
#' @param orderDecreasing flag indicating whether the value in \code{orderCol} should be considered as decreasing (as opposed
#'                 to increasing). \code{NULL} (default) for automatic determination.
#' @return ggplot object containing the plot
#'
#' @author Fabian Mueller
#' @export
#' @examples
#' \donttest{
#' library(RnBeads.hg19)
#' data(small.example.object)
#' logger.start(fname=NA)
#' # compute differential methylation
#' dm <- rnb.execute.computeDiffMeth(rnb.set.example,pheno.cols=c("Sample_Group","Treatment"))
#' # download LOLA DB
#' lolaDest <- tempfile()
#' dir.create(lolaDest)
#' lolaDirs <- downloadLolaDbs(lolaDest, dbs="LOLACore")
#' # perform enrichment analysis
#' res <- performLolaEnrichment.diffMeth(rnb.set.example,dm,lolaDirs[["hg19"]])
#' # select the 500 most hypermethylated tiling regions in ESCs compared to iPSCs
#' # in the example dataset
#' lolaRes <- res$region[["hESC vs. hiPSC (based on Sample_Group)"]][["tiling"]]
#' lolaRes <- lolaRes[lolaRes$userSet=="rankCut_500_hyper",]
#' # plot
#' lolaBarPlot(res$lolaDb, lolaRes, scoreCol="logOddsRatio", orderCol="maxRnk", pvalCut.neglog=-log(0.05))
#' }
lolaBarPlot <- function(lolaDb, lolaRes, scoreCol="pValueLog", orderCol=scoreCol, includedCollections=c(), pvalCut.neglog=-log(0.01), maxTerms=50, colorpanel=sample(rainbow(maxTerms,v=0.5)), groupByCollection=TRUE, change.pval.base=TRUE, orderDecreasing=NULL){
	if (length(unique(lolaRes[["userSet"]])) > 1){
		logger.warning("Multiple userSets contained in LOLA result object")
	}
	#prepare data.frame for plotting
	df2p <- lolaPrepareDataFrameForPlot(lolaDb, lolaRes, scoreCol=scoreCol, orderCol=orderCol, includedCollections=includedCollections, pvalCut.neglog=pvalCut.neglog, maxTerms=maxTerms, perUserSet=FALSE, groupByCollection=groupByCollection, change.pval.base=change.pval.base, orderDecreasing=orderDecreasing)

	aesObj <- aes_string("term", scoreCol)
	if (length(colorpanel) > 0) aesObj <- aes_string("term", scoreCol, fill="target")

	pp <- ggplot(df2p) + aesObj + geom_col() + 
		  scale_x_discrete(name="") + 
		  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
	if (length(colorpanel) > 0){
		cpanel <- rep(colorpanel, length.out=nlevels(df2p[["target"]]))
		names(cpanel) <- levels(df2p[["target"]])
		pp <- pp + scale_fill_manual(na.value="#C0C0C0", values=cpanel, guide=FALSE)
	} 
	if (groupByCollection){
		pp <- pp + facet_grid(. ~ collection, scales = "free", space = "free")
	}
	return(pp)
}
