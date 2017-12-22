########################################################################################################################
## immune.R
## created: 2017-11-20
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Functions related to the estimation of immune cell infiltration.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' Leukocytes unmethylation for purity
#'
#' Implementation of the LUMP (Leukocytes UnMethylation for Purity) algorithm for purity estimation on methylation
#' datasets.
#'
#' @param dataset Methylation dataset to study, provided as an object of type inheriting \code{RnBSet}.
#' @return Purity esimates provided as a \code{vector} of values in the range \code{[0, 1]}. The attribute
#'         \code{"sites"} contains the number of sites used in estimating the immune cell proportions. In case the
#'         dataset does not contain measurements for any of the sites on which LUMP focuses, the return values is
#'         \code{NULL}.
#'
#' @details The LUMP algorithm is developed by Dvir Aran, Marina Sirota and Atul J. Buttea.
#  It is described in their publication "Systematic pan-cancer analysis of tumour purity" in 2015.
#'
#' @author Yassen Assenov
#' @export
rnb.execute.lump <- function(dataset) {
	if (!inherits(dataset, "RnBSet")) {
		stop("Invalid value for dataset")
	}

	## Load the sites used by the LUMP algorithm
	load(system.file("data/lump.RData", package = "RnBeads"))
	x <- get0(paste0("lump.", dataset@assembly), inherits = FALSE)
	if (is.null(x)) {
		stop("Unsupported assembly")
	}
	x <- x[[dataset@target]]
	if (is.null(x)) {
		stop("Unsupported platform")
	}

	## Identify the indices of the sites in the dataset
	i.sites <- rep(0L, nrow(x)) # indices in the methylation matrix of the selected sites
	for (i.chrom in unique(x[, 1])) {
		i.x <- which(x[, 1] == i.chrom)
		i.d <- which(dataset@sites[, 2] == i.chrom)
		for (i in i.x) {
			j <- which(dataset@sites[i.d, 3] == x[i, 2])
			if (length(j) == 1) { i.sites[i] <- i.d[j] }
		}
	}
	rm(i.chrom, i.x, i.d, i, j)
	i.sites <- i.sites[i.sites != 0]
	if (length(i.sites) == 0) {
		stop(NULL)
	}

	result <- apply(dataset@meth.sites[i.sites, , drop = FALSE], 2, mean, na.rm = TRUE)
	result <- 1 - pmin(result / 0.85, 1)
	attr(result, "sites") <- length(i.sites)
	result
}

########################################################################################################################

#' rnb.section.lump
#'
#' Creates a report section dedicated to immune cell content estimation.
#'
#' @param report         Report on covariate inference to contain the immune cell content section. This must be an
#'                       object of type \code{\linkS4class{Report}}.
#' @param immune.content Estimation of the immune cell content, as returned by \code{\link{rnb.execute.lump}}.
#' @return The modified report.
#'
#' @author Yassen Assenov
#' @noRd
rnb.section.lump <- function(report, immune.content) {
	s.title <- "Immune Cell Content Estimation"
	txt <- NULL
	refText <- paste("Aran, D. et al. (2015) Systematic pan-cancer analysis of tumour purity.",
		"<i>Nature Communications</i> <b>6</b>:8971.")

	if (is.null(immune.content)) {
		report <- rnb.add.reference(report, refText)
		txt <- c("Immune cell content estimation was not performed because the dataset does not contain measurements ",
			"for any of the sites used in the LUMP algorithm ", rnb.get.reference(report, refText), ".")
		return(rnb.add.section(report, s.title, txt))
	} else if (is.character(immune.content)) {
		txt <- gsub("^Unsupported ", "", immune.content)
		txt <- c("Immune cell content estimation was not performed because the ", txt, " is not supported.")
		return(rnb.add.section(report, s.title, txt))
	}

	report <- rnb.add.reference(report, refText)
	iss <- length(immune.content) == 1
	txt <- c("Immune cell content estimation was performed using the LUMP algorithm ",
		rnb.get.reference(report, refText), ". ", ifelse(iss, "The estimate is a value", "Estimates are values"),
		" between 0 and 1 and ", ifelse(iss, "is", "are") , " based on ",
		ifelse(iss, "", "up to "), attr(immune.content, "sites"), " sites", ifelse(iss, "", " per sample"), ".")
	if (!iss) {
		txt <- c(txt, " The figure below shows the distribution of immune cell content values.")
	}
	report <- rnb.add.section(report, s.title, txt)
	if (!iss) {
		pp <- ggplot(data.frame(x = immune.content), aes_string("x")) +
			geom_histogram(breaks = seq(0, 1, length.out = 51)) + labs(x = "LUMP Estimate", y = "Frequency") +
			scale_x_continuous(breaks = seq(0, 1, length.out = 11), limits = c(0, 1), expand = c(0, 0)) +
			scale_y_continuous(expand = c(0, 0)) + theme(plot.margin = grid::unit(0.1 + c(0, 0.1, 0, 0), "in"))
		rplot <- createReportPlot("histogram_lump", report, width = 6.2, height = 5.2)
		print(pp)
		rplot <- off(rplot)
		txt <- paste("Histogram of all", length(immune.content), "estimates obtained after runnin LUMP.")
		report <- rnb.add.figure(report, txt, rplot)
	}

	return(report)
}
