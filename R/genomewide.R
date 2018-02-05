########################################################################################################################
## genomewide.R
## created: 2018-01-11
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Functions related to the calculation of genome-wide methylation.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' Genome-wide methylation level
#'
#' Computes genome-wide methylation levels per sample.
#'
#' @param dataset Methylation dataset to study, provided as an object of type inheriting \code{RnBSet}.
#' @return \code{vector} of values in the range \code{[0, 1]}, storing the average beta values per sample.
#'
#' @author Yassen Assenov
#' @export
rnb.execute.genomewide <- function(dataset) {
	if (!inherits(dataset, "RnBSet")) {
		stop("Invalid value for dataset")
	}

	N <- ncol(dataset@meth.sites)
	result <- rep(NA_real_, N)
	for (i in 1:N) {
		result[i] <- mean(dataset@meth.sites[, i], na.rm = TRUE)
	}
	result
}

########################################################################################################################

#' rnb.section.genomewide
#'
#' Creates a report section dedicated to the computed genome-wide methylation levels.
#'
#' @param report      Report on covariate inference to contain the genome-wide methylation section. This must be an
#'                    object of type \code{\linkS4class{Report}}.
#' @param meth.levels Calculated genome-wide methylation levels, as returned by \code{\link{rnb.execute.genomewide}}.
#' @return The modified report.
#'
#' @author Yassen Assenov
#' @noRd
rnb.section.genomewide <- function(report, meth.levels) {
	s.title <- "Genome-wide Methylation Levels"

	txt <- "because no beta value measurements are available"
	if (length(meth.levels) == 1) {
		if (is.na(meth.levels)) {
			txt <- paste0("Genome-wide methylation level could not be computed for the sample ", txt, ".")
		} else {
			txt <- c("Genome-wide methylation level was computed for the sample. The resulting value is ", meth.levels)
		}
		return(rnb.add.section(report, s.title, txt))
	}

	i.valid <- which(!is.na(meth.levels))
	M <- length(i.valid)
	if (M == 0) {
		txt <- paste0("Genome-wide methylation level could not be computed for any sample ", txt, ".")
		return(rnb.add.section(report, s.title, txt))
	}
	if (M == length(meth.levels)) {
		txt <- paste("Genome-wide methylation levels were computed for all", M, "samples.")
	} else {
		if (length(meth.levels) - M == 1L) {
			txt <- paste("It was not computed for the remaining sample", txt, "for it.")
		} else {
			txt <- paste("They were not computed for the remaining", length(meth.levels) - M, "samples", txt,
				"for them.")
		}
		if (M == 1) {
			txt <- paste("Genome-wide methylation level was computed for one sample only.", txt)
		} else {
			txt <- paste("Genome-wide methylation levels were computed for", M, "samples.", txt)
		}
	}
	txt <- c(txt, " The figure below shows ", ifelse(M == 1, "this value", "the distribution of these values"), ".")
	report <- rnb.add.section(report, s.title, txt)

	pp <- ggplot(data.frame(x = meth.levels[i.valid]), aes_string("x")) +
		geom_histogram(breaks = seq(0, 1, length.out = 51)) + labs(x = "Genome-wide methylation", y = "Frequency") +
		scale_x_continuous(breaks = seq(0, 1, length.out = 11), limits = c(0, 1), expand = c(0, 0)) +
		scale_y_continuous(expand = c(0, 0)) + theme(plot.margin = grid::unit(0.1 + c(0, 0.1, 0, 0), "in"))
	rplot <- createReportPlot("histogram_genomewide", report, width = 6.2, height = 5.2)
	print(pp)
	rplot <- off(rplot)
	txt <- c("Histogram of ", M, " genome-wide methylation level", ifelse(M == 1, "", "s"), ".")
	return(rnb.add.figure(report, txt, rplot))
}
