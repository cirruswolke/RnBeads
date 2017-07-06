########################################################################################################################
## lolaUtils.R
## created: 2017-07-06
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Methods for LOLA enrichment analysis
########################################################################################################################

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

	# trim the following names from the start of the description
	# trimNamesPattern.start <- paste0("^(", paste(c("ChIP", "UCSC", "cistrome_cistrome", "cistrome_epigenome"), collapse="|"), ")")
	trimNamesPattern.start <- paste0("^(", paste(c("UCSC"), collapse="|"), ")")
	res <- gsub(trimNamesPattern.start, "", res)

	# if both cellType and antibody annotation are present, construct a name from them
	hasDesc <- !is.na(tt$cellType) & !is.na(tt$antibody)
	if (sum(hasDesc) > 0) res[hasDesc] <- paste0(tt$antibody[hasDesc], " - ", tt$cellType[hasDesc])

	# special treatment for "encode_segmentation" from LOLACore
	isInCollection <- tt$collection=="encode_segmentation"
	if (sum(isInCollection) > 0){
		res[isInCollection] <- paste0(
			gsub(" Segments$", "", tt$description[isInCollection]),
			" - ", tt$cellType[isInCollection]
		)
	}

	# special treatment for "jaspar_motifs" from LOLAExt
	isInCollection <- tt$collection=="jaspar_motifs"
	if (sum(isInCollection) > 0) res[isInCollection] <- gsub("\\.bed$", "", tt$filename[isInCollection])

	# special treatment for "roadmap_epigenomics" from LOLAExt
	isInCollection <- tt$collection=="roadmap_epigenomics"
	if (sum(isInCollection) > 0) res[isInCollection] <- gsub("\\.(peaks\\.bed|narrowPeak)$", "", tt$filename[isInCollection])

	# add suffixes if requested
	if (addCollectionNames){
		res <- paste0(res, " [", tt$collection, "]")
	}
	if (addDbId){
		res <- paste0(res, " [db", 1:nrow(tt), "]")
	}
	return(res)
}

