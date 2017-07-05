library(argparse)
suppressPackageStartupMessages(library(RnBeads))

ap <- ArgumentParser()
ap$add_argument("-x", "--xml", action="store", help="Configuration xml file")
ap$add_argument("-s", "--rnbSet", action="store", dest="rnbSetLoc", help="Location of stored RnBSet")
ap$add_argument("-d", "--diffMethChunkId", action="append", dest="diffMethChunkIds", help="Ids for chunks for which differential methylation has already been computed. Will be used to construct the input file locations.")
ap$add_argument("-o", "--output", action="store", help="Output directory")
ap$add_argument("-c", "--cores", action="store", type="integer", default=1, help="Number of cores used for the analysis")
cmdArgs <- ap$parse_args()
module.name <- "differential_wrapup"

logger.start(fname=NA)
logger.status(c("...Started module:",module.name))

logger.start("Configuring Analysis")
	rnb.settings <- rnb.xml2options(cmdArgs$xml,return.full.structure=TRUE)

	report.dir <- rnb.settings$analysis.params[["dir.reports"]]
	analysis.options <- rnb.settings$options
	dm.chunk.ids <- cmdArgs$diffMethChunkIds

	if ("preanalysis.script" %in% names(rnb.settings)){
		source(rnb.settings$preanalysis.script)
	} 
	## Set options
	if (length(analysis.options) != 0) {
		do.call(rnb.options, analysis.options)
	}

	logger.machine.name()

	if (cmdArgs$cores > 1) {
		parallel.setup(cmdArgs$cores)
	}

	aname <- rnb.getOption("analysis.name")
	if (!(is.null(aname) || is.na(aname) || nchar(aname) == 0)) {
		logger.info(c("Analysis Title:", aname))
	}
	ncores <- parallel.getNumWorkers()
	if (ncores == -1) {
		ncores <- 1L
	}
	logger.info(c("Number of cores:", ncores))
	rm(aname, ncores)
logger.completed()

logger.start("Loading RnBSet")
	rnb.set <- load.rnb.set(cmdArgs$rnbSetLoc)
logger.completed()

logger.start(fname=c(file.path(report.dir,paste0("analysis_",module.name,".log")),NA))

################################################################################
# main script
################################################################################
if (length(dm.chunk.ids)<1){
	stop("Specify at least one id for a differential methylation chunk")
}

disk.dump <- rnb.getOption("disk.dump.big.matrices")

logger.start("Loading and Combining DiffMeth Files")
	chunk.id <- dm.chunk.ids[1]
	diffmeth.path <- file.path(cmdArgs$output,paste0("differential_chunk_",chunk.id,"_rnbDiffMeth"))
	logger.info(c("Loading:",diffmeth.path))
	diffmeth <- load.rnb.diffmeth(diffmeth.path)
	diffmeth.go.enrichment <- NULL
	enrichment.go.file <- file.path(diffmeth.path, "enrichment_go.RData")
	if (file.exists(enrichment.go.file)){
		logger.info(c("Loading:",enrichment.go.file))
		load(enrichment.go.file) #loads 'diffmeth.go.enrichment' object
	}
	diffmeth.lola.enrichment <- NULL
	enrichment.lola.file <- file.path(diffmeth.path, "enrichment_lola.RData")
	if (file.exists(enrichment.lola.file)){
		logger.info(c("Loading:",enrichment.lola.file))
		load(enrichment.lola.file) #loads 'diffmeth.lola.enrichment' object
	}

	diffmeth.res <- list(report=NULL,diffmeth=diffmeth,dm.go.enrich=diffmeth.go.enrichment,dm.lola.enrich=diffmeth.lola.enrichment)

	rm(diffmeth,diffmeth.go.enrichment,diffmeth.lola.enrichment) # to make sure a result.diffmeth object is not stored twice if a following file does not contain one
	RnBeads:::rnb.cleanMem()

	if (length(dm.chunk.ids)>1){
		for (chunk.id in dm.chunk.ids[2:length(dm.chunk.ids)]){
			diffmeth.path <- file.path(cmdArgs$output,paste0("differential_chunk_",chunk.id,"_rnbDiffMeth"))
			logger.info(c("Loading:",diffmeth.path))
			diffmeth <- load.rnb.diffmeth(diffmeth.path)
			enrichment.go.file <- file.path(diffmeth.path, "enrichment_go.RData")
			diffmeth.go.enrichment <- NULL
			if (file.exists(enrichment.go.file)){
				logger.info(c("Loading:",enrichment.go.file))
				load(enrichment.go.file) #loads 'diffmeth.go.enrichment' object
			}
			enrichment.lola.file <- file.path(diffmeth.path, "enrichment_lola.RData")
			diffmeth.lola.enrichment <- NULL
			if (file.exists(enrichment.lola.file)){
				logger.info(c("Loading:",enrichment.lola.file))
				load(enrichment.lola.file) #loads 'diffmeth.lola.enrichment' object
			}
			diffmeth.res.cur <- list(report=NULL,diffmeth=diffmeth,dm.go.enrich=diffmeth.go.enrichment,dm.lola.enrich=diffmeth.lola.enrichment)

			ll <- list(diffmeth.res,diffmeth.res.cur)
			diffmeth.res <- RnBeads:::combine.diffMeth.objs(ll)

			destroy(diffmeth)
			rm(diffmeth.res.cur,diffmeth,diffmeth.go.enrichment,diffmeth.lola.enrichment) # to make sure a result.diffmeth object is not stored twice if a following file does not contain one
			RnBeads:::rnb.cleanMem()
		}
	}
	if (!is.valid(diffmeth.res$diffmeth,verbose=TRUE)){
		stop("RnBDiffMeth object is invalid")
	}
logger.completed()

logger.start("Saving combined data")
	diffmeth.path <- file.path(cmdArgs$output,paste0("differential_rnbDiffMeth"))
	save.rnb.diffmeth(diffmeth.res$diffmeth, diffmeth.path)
	diffmeth.go.enrichment <- diffmeth.res$dm.go.enrich
	if (!is.null(diffmeth.go.enrichment)){
		save(diffmeth.go.enrichment, file=file.path(diffmeth.path, "enrichment_go.RData"))
	}
	diffmeth.lola.enrichment <- diffmeth.res$dm.lola.enrich
	if (!is.null(diffmeth.lola.enrichment)){
		save(diffmeth.lola.enrichment, file=file.path(diffmeth.path, "enrichment_lola.RData"))
	}
logger.completed()

logger.start("Differential Methylation")
	logger.start("Report Generation")
		aname <- rnb.getOption("analysis.name")
		if (!(is.null(aname) || is.na(aname) || nchar(aname) == 0)) {
			logger.info(c("Analysis Title:", aname))
		}
		ncores <- parallel.getNumWorkers()
		if (ncores == -1) {
			ncores <- 1L
		}
		logger.info(c("Number of cores:", ncores))
		rm(aname, ncores)

		diffmeth <- diffmeth.res$diffmeth
		dm.go.enrich <- diffmeth.res$dm.go.enrich
		dm.lola.enrich <- <- diffmeth.res$dm.lola.enrich

		init.configuration <- !file.exists(file.path(report.dir, "configuration"))
		report <- RnBeads:::init.pipeline.report("differential_methylation", report.dir, init.configuration)
		optionlist <- rnb.options("analyze.sites","region.types", "differential.permutations", "differential.comparison.columns",
			"differential.comparison.columns.all.pairwise","columns.pairing","differential.site.test.method","covariate.adjustment.columns",
			"differential.adjustment.sva","differential.adjustment.celltype","differential.enrichment.go","differential.enrichment.lola")
		report <- RnBeads:::rnb.add.optionlist(report, optionlist)
		
		if (is.null(diffmeth)){
			txt <- "Differential methylation analyis was skipped because no valid grouping information could be found."
			report <- rnb.add.section(report, "Differential Methylation Analysis", txt)
		} else {
			gz <- rnb.getOption("gz.large.files")
			includeSites <- rnb.getOption("analyze.sites") && rnb.getOption("differential.report.sites")
			report <- RnBeads:::rnb.section.diffMeth.introduction(diffmeth,report)
			if (includeSites){
				report <- RnBeads:::rnb.section.diffMeth.site(rnb.set,diffmeth,report,gzTable=gz)
			}
			if (length(get.region.types(diffmeth))>0){
				report <- RnBeads:::rnb.section.diffMeth.region(rnb.set,diffmeth,report,dm.go.enrich=dm.go.enrich,dm.lola.enrich=dm.lola.enrich,gzTable=gz)
			}
		}

		off(report)

	logger.completed()
logger.completed()

logger.status(c("...Completed module:",module.name))
quit(save='no')
