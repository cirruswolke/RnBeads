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
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{RnBeads:::loadLolaDbs}
#' @return character vector with cell types
#'
#' @author Fabian Mueller
#' @noRd
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
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{RnBeads:::loadLolaDbs}
#' @return character vector with targets
#'
#' @author Fabian Mueller
#' @noRd
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
#' @param lolaDb   LOLA DB object as returned by \code{LOLA::loadRegionDB} or \code{RnBeads:::loadLolaDbs}
#' @param addCollectionNames  attach the name of the collection to the name
#' @param addDbId  attach the index of the item in the LOLA DB object to the name
#' @return character vector with human readable names
#'
#' @author Fabian Mueller
#' @noRd
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

