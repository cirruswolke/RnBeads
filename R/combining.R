########################################################################################################################
## combining.R
## created: 2017-04-13
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Functions related to combining datasets of the same or different platforms, differential methylation results, etc.
########################################################################################################################

## TODO: Move the "combine" method here; use the function rnb.combine.arrays in it

## F U N C T I O N S ###################################################################################################

#' Combine dataset slots
#' 
#' Combines the data matrices of two datasets into a single one.
#' 
#' @param m1    Matrix of the first dataset. Set this to \code{NULL} if the first dataset does not contain this slot.
#' @param m2    Matrix of the second dataset. Set this to \code{NULL} if the second dataset does not contain this slot.
#' @param ii    Two-column \code{integer} \code{matrix} of indices in \code{m1} and \code{m2} to be considered for the
#'              combined matrix.
#' @param nn    Two-element \code{integer} \code{vector} listing the number of samples in the first and second dataset,
#'              respectively. If \code{m1} is provided, it should contain \code{nn[1]} columns.
#' @param useff Flag indicating if disk dumping for the output matrix must be enabled.
#' @return When \code{m1} and \code{m2} are both \code{NULL}, the return value is also \code{NULL}. Otherwise, the
#'         constructed combined matrix as an \code{ff} object (if \code{useff} is \code{TRUE}), or as an instance of
#'         \code{matrix} (if \code{useff} is \code{FALSE}). The dimensions of the matrix are \code{nrow(ii)} x
#'         \code{sum(nn)}.
#' 
#' @author Yassen Assenov
#' @noRd
rnb.combine.matrices <- function(m1, m2, ii, nn, useff = rnb.getOption("disk.dump.big.matrices")) {
	if (is.null(m1)) {
		if (is.null(m2)) {
			return(NULL)
		}
		if ("ff" %in% class(m2)) {
			d.type <- vmode(m2)
		} else {
			d.type <- typeof(m2)
		}
	} else if ("ff" %in% class(m1)) {
		d.type <- vmode(m1)
	} else {
		d.type <- typeof(m1)
	}
	na.value <- do.call(paste0("as.", d.type), list(x = NA))
	if (useff) {
		mm <- ff(na.value, dim = c(nrow(ii), sum(nn)), dimnames = list(NULL, NULL), vmode = d.type)
	} else {
		mm <- matrix(na.value, nrow = nrow(ii), ncol = sum(nn))
	}
	if (!is.null(m1)) {
		for (i in 1:(nn[1])) {
			mm[, i] <- m1[ii[, 1], i]
		}
	}
	if (!is.null(m2)) {
		for (i in 1:(nn[2])) {
			mm[, nn[1] + i] <- m2[ii[, 2], i]
		}
	}
	mm
}

########################################################################################################################

#' Combine pheno tables
#' 
#' Combines two sample annotation tables.
#' 
#' @param dataset1 First input dataset as an object of type inheriting \code{\linkS4class{RnBeadSet}}.
#' @param dataset2 Second input dataset as an object of type inheriting \code{\linkS4class{RnBeadSet}}.
#' @return Combined sample annotation table as a \code{data.frame}.
#' 
#' @author Yassen Assenov
#' @noRd
rnb.combine.pheno <- function(dataset1, dataset2) {
	tbl1 <- dataset1@pheno
	tbl2 <- dataset2@pheno
	if (!identical(sapply(tbl1, class), sapply(tbl2, class))) {
		stop("Mismatch in sample annotation tables")
	}
	for (i in which(sapply(tbl1, class) == "factor")) {
		if (!identical(levels(tbl1[, i]), levels(tbl2[, i]))) {
			stop("Mismatch in sample annotation tables")
		}
	}
	rbind(tbl1, tbl2)
}

########################################################################################################################

#' Combine array-based datasets
#'
#' Concatenates two array-based datasets focusing on the common probes.
#'
#' @param dataset1 First input dataset as an object of type inheriting \code{\linkS4class{RnBeadSet}}.
#' @param dataset2 Second input dataset as an object of type inheriting \code{\linkS4class{RnBeadSet}}.
#' @return Combined dataset as an object of type inheriting \code{\linkS4class{RnBeadSet}}.
#'
#' @details \describe{
#'   \item{Sample annotation tables}{This method expects that the sample annotation tables of the two datasets have
#'      identical structures.}
#'   \item{Genome assembly}{This method expects that the two datasets target the same genome assembly.}
#'   \item{Platform}{The platform of the combined dataset is the most recent among the platforms of the input datasets.}
#'   \item{Intensity values}{The combined dataset is of type \code{\linkS4class{RnBeadRawSet}} only when both input
#'      datasets are of this type. Otherwise, any intensity value data is ignored.}
#'   \item{Probes}{Only the common probes are included in the resulting dataset.}
#'   \item{Regions}{Regions summarized in any of the input datasets are ignored. In the resulting dataset, regions are
#'      summarized as specified in the analysis option \code{"region.types"}.}
#'   \item{Quality control data}{QC data in the input datasets is ignored. The combined dataset includes no data on QC
#'      probe intensities.}
#'   \item{Infered covariates}{Inferred covariates in the input datasets are ignored. The combined dataset includes no
#'      data on inferred covariates.}
#'   \item{Disk dumping}{The combined dataset stores big tables on disk when the analysis option
#'      \code{"disk.dump.big.matrices"} is enabled.}
#' }
#' 
#' @author Yassen Assenov
#' @export
rnb.combine.arrays <- function(dataset1, dataset2) {
	if (!inherits(dataset1, "RnBeadSet")) {
		stop("Invalid value for dataset1")
	}
	if (!inherits(dataset2, "RnBeadSet")) {
		stop("Invalid value for dataset2")
	}
	if (dataset1@assembly != dataset2@assembly) {
		stop("Incompatible assemblies")
	}
	i <- c(dataset1@target, dataset2@target)
	common.platform <- c("probesEPIC" = "EPIC", "probes450" = "450k", "probes27" = "27k")
	if (!(i[1] %in% names(common.platform))) {
		stop("Unsupported platform for dataset1")
	}
	if (!(i[2] %in% names(common.platform))) {
		stop("Unsupported platform for dataset2")
	}
	common.platform <- common.platform[intersect(names(common.platform), unique(i))]
	common.platform <- unname(common.platform[1])
	is.raw <- (inherits(dataset1, "RnBeadRawSet") && inherits(dataset2, "RnBeadRawSet"))

	## Combine sample annotation tables
	tbl <- rnb.combine.pheno(dataset1, dataset2)
	nn <- c(nrow(dataset1@pheno), nrow(dataset2@pheno))

	## Identify common probes
	common.sites <- intersect(rownames(dataset1@sites), rownames(dataset2@sites))
	if (length(common.sites) == 0) {
		stop("No common sites identified")
	}
	ii <- cbind( # probe names must be sorted in both datasets!
		which(rownames(dataset1@sites) %in% common.sites),
		which(rownames(dataset2@sites) %in% common.sites))
	useff <- rnb.getOption("disk.dump.big.matrices")
	
	## Combine the data matrices
	result <- list(
		"pheno" = tbl,
		"probes" = common.sites,
		"platform" = common.platform,
		"region.types" = rnb.region.types.for.analysis(dataset1@assembly),
		"qc" = NULL,
		"useff" = useff,
		"ffcleanup" = rnb.getOption("enforce.destroy.disk.dumps"))
	slot.names <- c("p.values")
	if (is.raw) {
		slot.names <- c(slot.names, "M", "U", "M0", "U0", "bead.counts.M", "bead.counts.U")
	} else {
		slot.names <- c(slot.names, "betas", "bead.counts")
	}
	for (slot.name in slot.names) {
		if (slot.name == "p.values") {
			s.name <- "pval.sites"
		} else if (slot.name == "betas") {
			s.name <- "meth.sites"
		} else {
			s.name <- slot.name
		}
		result[[slot.name]] <- rnb.combine.matrices(slot(dataset1, s.name), slot(dataset2, s.name), ii, nn, useff)
	}
	
	## Construct the resulting object
	result <- do.call(ifelse(is.raw, "RnBeadRawSet", "RnBeadSet"), result)
	rm(tbl, i, common.platform, is.raw, nn, common.sites, ii, useff, slot.names, slot.name, s.name)

	## Set normalization and background subtraction methods
	if (dataset1@status$normalized == dataset2@status$normalized) {
		result@status$normalized <- dataset1@status$normalized
	} else {
		result@status$normalized <- "none"
		warning("Inconsistent normalization methods; setting to none")
	}
	if (dataset1@status$background == dataset2@status$background) {
		result@status$background <- dataset1@status$background
	} else {
		result@status$background <- "none"
		warning("Inconsistent background subtraction methods; setting to none")
	}
	result
}
