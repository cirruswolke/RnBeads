################################################################################
# RnBeadsDJ - RnBeads Data Juggler
# A Shiny app for configuring and running RnBeads analyses
#-------------------------------------------------------------------------------
# Main app file
#-------------------------------------------------------------------------------
# created: 2017-10-04
# creator: Fabian Mueller
#-------------------------------------------------------------------------------
# Run using:
# shiny::runApp(file.path('RnBeadsDJ'))
################################################################################
library(shiny)
library(shinyjs)
# library(shinyFiles)
library(RnBeads)

################################################################################
# Options
################################################################################
options(shiny.maxRequestSize=50*1024^2)

################################################################################
# globals
################################################################################
REFRESH.TIME <- 5000
RNB.MODULES <- c(
	"Data Import"="data_import",
	"Quality Control"="quality_control",
	"Preprocessing"="preprocessing",
	"Tracks and Tables"="tracks_and_tables",
	"Covariate Inference"="covariate_inference",
	"Exploratory Analysis"="exploratory_analysis",
	"Differential Methylation"="differential_methylation"
)
RNB.MODULES.LOG.MSG <- c(
	"data_import"="Loading Data",
	"quality_control"="Quality Control",
	"preprocessing"="Preprocessing",
	"tracks_and_tables"="Tracks and Tables",
	"covariate_inference"="Covariate Inference",
	"exploratory_analysis"="Exploratory Analysis",
	"differential_methylation"="Differential Methylation"
)
RNB.PLATFORMS <- c("Bisulfite Sequencing"="biseq", "Illumina EPIC"="illEpic", "Illumina 450k"="ill450k", "Illumina27k"="ill27k")
RNB.ASSEMBLIES <- rnb.get.assemblies()
RNB.BED.STYLES <- c("BisSNP"="BisSNP", "ENCODE"="Encode", "EPP"="EPP", "Bismark cytosine"="bismarkCytosine", "Bismark coverage"="bismarkCov")
RNB.FILTERING.SNP <- c("No filtering"="no", "3 SNPs"="3", "5 SNPs"="5", "Any SNPs"="any")
RNB.NORMALIZATION.METHODS=c("none", "bmiq", "illumina", "swan", "minfi.funnorm", "wm.dasen", "wm.nasen", "wm.betaqn", "wm.naten", "wm.nanet", "wm.nanes", "wm.danes", "wm.danet", "wm.danen", "wm.daten1", "wm.daten2", "wm.tost", "wm.fuks", "wm.swan")
RNB.NORMALIZATION.BG.METHODS <- c("none", "methylumi.noob", "methylumi.goob", "enmix.oob")
RNB.IMPUTATION.METHODS <- c("none", "mean.cpgs", "mean.samples", "random", "knn")
RNB.TRACKHUB.FORMATS <- c("bigBed", "bigWig")
RNB.SVA.NUM.METHODS <- c("leek", "be")
RNB.DIFFMETH.TEST.METHODS <- c("limma", "refFreeEWAS")
RNB.DIFFVAR.METHODS <- c("diffVar", "iEVORA")
RNB.COLSCHEMES.CATEGORY <- list(
	default=c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"),
	extended=c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#2166AC","#B2182B","#00441B","#40004B","#053061","#003D7C","#D50911")
)
RNB.COLSCHEMES.METH <- list(
	default=c("#AD0021","#909090","#39278C"),
	YlBl=c("#EDF8B1","#41B6C4","#081D58")
)

RNB.OPTION.DESC.TAB <- RnBeads:::rnb.options.description.table()
RNB.OPTION.DESC <- sapply(names(rnb.options()), FUN=function(x){
	if (is.element(x, rownames(RNB.OPTION.DESC.TAB))){
		return(RNB.OPTION.DESC.TAB[x, "desc"])
	} else {
		return("See the help pages: '?rnb.options'")
	}
})
################################################################################
# Choose local file or directory
# adapted from https://github.com/wleepang/shiny-directory-input
################################################################################
# Interactive Choosers for MacOS and Linux (For Windos this is already implemented in the utils package)
#' Choose a Folder Interactively (Mac OS X)
#'
#' Display a folder selection dialog under Mac OS X
#'
#' @param default which folder to show initially
#' @param caption the caption on the selection dialog
#'
#' @details
#' Uses an Apple Script to display a folder selection dialog.  With \code{default = NA},
#' the initial folder selection is determined by default behavior of the
#' "choose folder" Apple Script command.  Otherwise, paths are expanded with
#' \link{path.expand}.
#'
#' @return
#' A length one character vector, character NA if 'Cancel' was selected.
#'
if (Sys.info()['sysname'] == 'Darwin') {
	choose.dir <- function(default = NA, caption = NA, ...) {
		command = 'osascript'
		args = '-e "POSIX path of (choose folder{{prompt}}{{default}})"'

		if (!is.null(caption) && !is.na(caption) && nzchar(caption)) {
			prompt = sprintf(' with prompt \\"%s\\"', caption)
		} else {
			prompt = ''
		}
		args = sub('{{prompt}}', prompt, args, fixed = T)

		if (!is.null(default) && !is.na(default) && nzchar(default)) {
			default = sprintf(' default location \\"%s\\"', path.expand(default))
		} else {
			default = ''
		}
		args = sub('{{default}}', default, args, fixed = T)

		suppressWarnings({
			path = system2(command, args = args, stderr = TRUE)
			path <- path[length(path)]
		})
		if (!file.exists(path)) {
			# user canceled
		 	path = NA
		}

		return(path)
	}

	choose.files <- function(default = NA, caption = NA, ...) {
		command = 'osascript'
		args = '-e "POSIX path of (choose file{{prompt}}{{default}})"'

		if (!is.null(caption) && !is.na(caption) && nzchar(caption)) {
			prompt = sprintf(' with prompt \\"%s\\"', caption)
		} else {
			prompt = ''
		}
		args = sub('{{prompt}}', prompt, args, fixed = T)

		if (!is.null(default) && !is.na(default) && nzchar(default)) {
			default = sprintf(' default location \\"%s\\"', path.expand(default))
		} else {
			default = ''
		}
		args = sub('{{default}}', default, args, fixed = T)

		suppressWarnings({
			path = system2(command, args = args, stderr = TRUE)
			path <- path[length(path)]
		})
		if (!file.exists(path)) {
			# user canceled
			path = NA
		}

		return(path)
	}
} else if (Sys.info()['sysname'] == 'Linux') {
	library(tcltk)
	choose.dir <- tk_choose.dir
	choose.files <- tk_choose.files
}

#' Directory Selection Control
#'
#' Create a directory selection control to select a directory on the server
#'
#' @param inputId The \code{input} slot that will be used to access the value
#' @param label Display label for the control, or NULL for no label
#' @param value Initial value.  Paths are exapnded via \code{\link{path.expand}}.
#'
#' @details
#' This widget relies on \link{\code{choose.dir}} to present an interactive
#' dialog to users for selecting a directory on the local filesystem.  Therefore,
#' this widget is intended for shiny apps that are run locally - i.e. on the
#' same system that files/directories are to be accessed - and not from hosted
#' applications (e.g. from shinyapps.io).
#'
#' @return
#' A directory input control that can be added to a UI definition.
#'
#' @seealso
#' \link{updateDirectoryInput}, \link{readDirectoryInput}, \link[utils]{choose.dir}
directoryInput = function(inputId, label, value = NULL) {
	if (!is.null(value) && !is.na(value)) {
		value = path.expand(value)
	}

	tagList(
		singleton(
			tags$head(
				tags$script(src = 'js/directory_input_binding.js')
			)
		),

		div(
			class = 'form-group directory-input-container',
			shiny:::`%AND%`(label, tags$label(label)),
			div(
				span(
					class = 'col-xs-9 col-md-11',
					style = 'padding-left: 0; padding-right: 5px;',
					div(
						class = 'input-group shiny-input-container',
						style = 'width:100%;',
						div(class = 'input-group-addon', icon('folder-o')),
						tags$input(
							id = sprintf('%s__chosen_dir', inputId),
							value = value,
							type = 'text',
							class = 'form-control directory-input-chosen-dir',
							readonly = 'readonly'
						)
					)
				),
				span(
					class = 'shiny-input-container',
					tags$button(
						id = inputId,
						class = 'btn btn-default directory-input',
						'...'
					)
				)
			)
		)
	)
}

#' Local File Selection Control
#'
#' Create a local file selection control to select a local file on the server
#'
#' @param inputId The \code{input} slot that will be used to access the value
#' @param label Display label for the control, or NULL for no label
#' @param value Initial value.  Paths are exapnded via \code{\link{path.expand}}.
#'
#' @details
#' This widget relies on \link{\code{choose.files}} to present an interactive
#' dialog to users for selecting files on the local filesystem.  Therefore,
#' this widget is intended for shiny apps that are run locally - i.e. on the
#' same system that files/directories are to be accessed - and not from hosted
#' applications (e.g. from shinyapps.io).
#'
#' @return
#' A local file input control that can be added to a UI definition.
#'
#' @seealso
#' \link{updateLocalFileInput}, \link{readLocalFileInput}, \link[utils]{choose.files}
localFileInput = function(inputId, label, value = NULL) {
	if (!is.null(value) && !is.na(value)) {
		value = path.expand(value)
	}

	tagList(
		singleton(
			tags$head(
				tags$script(src = 'js/localfile_input_binding.js')
			)
		),

		div(
			class = 'form-group localfile-input-container',
			shiny:::`%AND%`(label, tags$label(label)),
			div(
				span(
					class = 'col-xs-9 col-md-11',
					style = 'padding-left: 0; padding-right: 5px;',
					div(
						class = 'input-group shiny-input-container',
						style = 'width:100%;',
						div(class = 'input-group-addon', icon('file-o')),
						tags$input(
							id = sprintf('%s__chosen_file', inputId),
							value = value,
							type = 'text',
							class = 'form-control localfile-input-chosen-file',
							readonly = 'readonly'
						)
					)
				),
				span(
					class = 'shiny-input-container',
					tags$button(
						id = inputId,
						class = 'btn btn-default localfile-input',
						'...'
					)
				)
			)
		)
	)
}

#' Change the value of a directoryInput on the client
#'
#' @param session The \code{session} object passed to function given to \code{shinyServer}.
#' @param inputId The id of the input object.
#' @param value A directory path to set
#' @param ... Additional arguments passed to \link{\code{choose.dir}}.  Only used
#'    if \code{value} is \code{NULL}.
#'
#' @details
#' Sends a message to the client, telling it to change the value of the input
#' object.  For \code{directoryInput} objects, this changes the value displayed
#' in the text-field and triggers a client-side change event.  A directory
#' selection dialog is not displayed.
#'
updateDirectoryInput = function(session, inputId, value = NULL, ...) {
  if (is.null(value)) {
    value = choose.dir(...)
  }
  session$sendInputMessage(inputId, list(chosen_dir = value))
}

#' Change the value of a localFileInput on the client
#'
#' @param session The \code{session} object passed to function given to \code{shinyServer}.
#' @param inputId The id of the input object.
#' @param value A file path to set
#' @param ... Additional arguments passed to \link{\code{choose.files}}.  Only used
#'    if \code{value} is \code{NULL}.
#'
#' @details
#' Sends a message to the client, telling it to change the value of the input
#' object.  For \code{localFileInput} objects, this changes the value displayed
#' in the text-field and triggers a client-side change event.  A file
#' selection dialog is not displayed.
#'
updateLocalFileInput = function(session, inputId, value = NULL, ...) {
  if (is.null(value)) {
    value = choose.files(...)
  }
  session$sendInputMessage(inputId, list(chosen_file = value))
}

#' Read the value of a directoryInput
#'
#' @param session The \code{session} object passed to function given to \code{shinyServer}.
#' @param inputId The id of the input object
#'
#' @details
#' Reads the value of the text field associated with a \code{directoryInput}
#' object that stores the user selected directory path.
#'
readDirectoryInput = function(session, inputId) {
  session$input[[sprintf('%s__chosen_dir', inputId)]]
}

#' Read the value of a localFileInput
#'
#' @param session The \code{session} object passed to function given to \code{shinyServer}.
#' @param inputId The id of the input object
#'
#' @details
#' Reads the value of the text field associated with a \code{localFileInput}
#' object that stores the user selected file path.
#'
readLocalFileInput = function(session, inputId) {
  session$input[[sprintf('%s__chosen_file', inputId)]]
}

observeDirectoryInput <- function(input, session, inputId){
	observeEvent(ignoreNULL=TRUE, eventExpr={input[[inputId]]},
		handlerExpr={
			if (input[[inputId]] > 0) {	      
				# launch the directory selection dialog with initial path read from the widget
				path = choose.dir(default = readDirectoryInput(session, inputId))
				# update the widget value
				updateDirectoryInput(session, inputId, value=path)
			}
		}
	)
}

observeLocalFileInput <- function(input, session, inputId){
	observeEvent(ignoreNULL=TRUE, eventExpr={input[[inputId]]},
		handlerExpr={
			if (input[[inputId]] > 0) {
				# launch the file selection dialog with initial path read from the widget 
				filt <- matrix(c("comma-separated", ".csv", "tab-separated", ".tsv", "Text", ".txt", "All files", "*"), 4, 2, byrow = TRUE)
				path = choose.files(default = readLocalFileInput(session, inputId), filters=filt)
				# update the widget value
				updateLocalFileInput(session, inputId, value=path)
			}
		}
	)
}

################################################################################
# Little Helpers
################################################################################
plotColPal <- function(colpal){
	require(plotrix)
	par(mar=c(0,0,0,0))
	plot.new()
	gradient.rect(0,0,1,1,col=colpal,nslices=length(colpal),gradient="x",border=NA)
}

getRnbStatusFromLog <- function(logFile){
	if (!file.exists(logFile)) return(NULL)
	ll <- readLines(logFile)
	res <- data.frame(
		module=RNB.MODULES,
		status=NA,
		scheduled=NA
	)
	rownames(res) <- RNB.MODULES
	for (mm in RNB.MODULES){
		if(sum(grepl(paste0("COMPLETED ", RNB.MODULES.LOG.MSG[mm]), ll))>0){
			res[mm, "status"] <- "completed"
		} else if(sum(grepl(paste0("STARTED ", RNB.MODULES.LOG.MSG[mm]), ll))>0){
			res[mm, "status"] <- "started"
		}
	}
	#check the options if the anlysis is scheduled
	res["data_import", "scheduled"]              <- rnb.getOption("import")
	res["quality_control", "scheduled"]          <- rnb.getOption("qc")
	res["preprocessing", "scheduled"]            <- rnb.getOption("preprocessing")
	res["tracks_and_tables", "scheduled"]        <- rnb.getOption("export.to.csv") || rnb.getOption("export.to.bed") || length(rnb.getOption("export.to.trackhub")) > 0
	res["covariate_inference", "scheduled"]      <- rnb.getOption("inference")
	res["exploratory_analysis", "scheduled"]     <- rnb.getOption("exploratory")
	res["differential_methylation", "scheduled"] <- rnb.getOption("import")
	return(res)
}

checkReportDir <- function(repDir){
	res <- list(valid=FALSE)
	res$reportHtml <- rep(NA, length(RNB.MODULES))
	names(res$reportHtml) <- RNB.MODULES
	res$moduleStatus <- list(NULL)
	res$createdByDj <- FALSE # check if this report directory was created by RnBeadsDJ

	if (!dir.exists(repDir)) return(res)
	contents <- list.files(repDir, include.dirs=TRUE, all.files=TRUE)
	res$valid <- all(c("configuration", "analysis.log") %in% contents)
	if (!res$valid) return(res)
	res$createdByDj <- is.element(".RnBeadsDJ", contents) 
	res$valid <- res$createdByDj || all(c("index.html", "analysis_options.xml") %in% contents)
	if (!res$valid) return(res)
	htmlModules <- paste0(RNB.MODULES, ".html")
	htmlModules.exist <- htmlModules %in% contents
	res$reportHtml[htmlModules.exist] <- htmlModules[htmlModules.exist]
	res$moduleStatus <- getRnbStatusFromLog(file.path(repDir, "analysis.log"))
	return(res)
}
# mark a report directory as created by RnBeadsDJ
markDirDJ <- function(repDir){
	repStatus <- checkReportDir(repDir)
	if (!repStatus$createdByDj){
		file.create(file.path(repDir, ".RnBeadsDJ"))
	}
	invisible(NULL)
}
################################################################################
# UI configuration
################################################################################
ui <- tagList(useShinyjs(), navbarPage(
	windowTitle="RnBeadsDJ",
	tags$p(tags$a(href="http://rnbeads.mpi-inf.mpg.de", tags$img(width=145, height=50, src="img/rnbeads_logo.png")), "DJ"),
	# tabPanel("Sandbox", icon=icon("dropbox"),
	# 	verbatimTextOutput("sandboxOut"),
	# 	localFileInput('sandboxIn2', label='select a local file'),
	# 	verbatimTextOutput("sandboxOut2")
	# ),
	tabPanel("About", icon=icon("book"),
		tags$h1("Welcome to the RnBeads Data Juggler"),
		tags$p("Here, you can configure and run your RnBeads analyses. You can ..."),
		tags$ul(
			tags$li("... run new analyses by specifying a non-existing data directory in the", "'Analysis'", "tab."),
			tags$li("... view the status of an RnBeads analysis by specifying an existing report directory in the", "'Analysis'", "tab."),
			tags$li("... run individual RnBeads modules via the", "'Modules'", "tab.")
		)
	),
	tabPanel("Analysis", icon=icon("bar-chart"),
		sidebarPanel(
			directoryInput('outDir', label = 'Select analysis directory', value = '~'),
			textInput("reportSubDir", "Choose the name of the report directory", "rnbeads_report"),
			sliderInput('numCores', "Select the number of cores to use", min=1, max=detectCores(), value=1, step=1),
			selectInput('ggplotTheme', "Select a theme for the plots", c("Black & White"="bw", "Grey"="grey"))
		),
		mainPanel(
			uiOutput("anaStatus")
		)
	),
	tabPanel("Input", icon=icon("sign-in"),
		sidebarPanel(
			localFileInput("sampleAnnotFile", "Select sample annotation file"),
			selectInput("rnbOptsI.import.table.separator", "Separator:", c("comma" = ",", "tab"="\t")),
			selectInput("platform", "Platform", RNB.PLATFORMS),
			directoryInput('dataDir', label='Select input data directory')
		),
		mainPanel(
			uiOutput("inputStatus")
		)
	),
	tabPanel("Analysis Options", icon=icon("sliders"),
		tabsetPanel(
			tabPanel("General",
				tags$table(class="table table-hover",
					tags$thead(tags$tr(
						tags$th("Name"),
						tags$th("Setting"),
						tags$th("Value")
					)),
					tags$tbody(
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["analysis.name"], tags$code("analysis.name"))
							),
							tags$td(
								textInput("rnbOptsI.analysis.name", NULL, "RnBeads Analysis")
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.analysis.name")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["assembly"], tags$code("assembly"))
							),
							tags$td(
								selectInput("rnbOptsI.assembly", NULL, RNB.ASSEMBLIES, selected="hg38")
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.assembly", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["region.types"], tags$code("region.types"))
							),
							tags$td(
								uiOutput("selRegionTypes")
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.region.types", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["identifiers.column"], tags$code("identifiers.column"))
							),
							tags$td(
								uiOutput('selColumn.id')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.identifiers.column", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["min.group.size"], tags$code("min.group.size"))
							),
							tags$td(
								sliderInput("rnbOptsI.min.group.size", NULL, min=1, max=20, value=2)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.min.group.size")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["max.group.count"], tags$code("max.group.count"))
							),
							tags$td(
								sliderInput("rnbOptsI.max.group.count", NULL, min=2, max=20, value=20)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.max.group.count")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["colors.category"], tags$code("colors.category"))
							),
							tags$td(
								selectInput("rnbOptsI.colors.category", NULL, names(RNB.COLSCHEMES.CATEGORY))
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.colors.category"),
								plotOutput("rnbOptsOP.colors.category", height="30px")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["colors.meth"], tags$code("colors.meth"))
							),
							tags$td(
								selectInput("rnbOptsI.colors.meth", NULL, names(RNB.COLSCHEMES.METH))
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.colors.meth"),
								plotOutput("rnbOptsOP.colors.meth", height="30px")
							)
						)
					)
				)
			),
			tabPanel("Import",
				tags$table(class="table table-hover",
					tags$thead(tags$tr(
						tags$th("Name"),
						tags$th("Setting"),
						tags$th("Value")
					)),
					tags$tbody(
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["import.default.data.type"], tags$code("import.default.data.type"))
							),
							tags$td(
								"Defined per 'Platform' in the 'Input' section"
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.import.default.data.type")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["import.table.separator"], tags$code("import.table.separator"))
							),
							tags$td(
								"Defined per 'Separator' in the 'Input' section"
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.import.table.separator")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["import.bed.style"], tags$code("import.bed.style"))
							),
							tags$td(
								selectInput("rnbOptsI.import.bed.style", NULL, names(RNB.BED.STYLES))
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.import.import.bed.style")
							)
						)
					)
				)
			),
			tabPanel("Preprocessing",
				tags$table(class="table table-hover",
					tags$thead(tags$tr(
						tags$th("Name"),
						tags$th("Setting"),
						tags$th("Value")
					)),
					tags$tbody(
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.coverage.threshold"], tags$code("filtering.coverage.threshold"))
							),
							tags$td(
								sliderInput("rnbOptsI.filtering.coverage.threshold", NULL, min=1, max=100, value=5)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.coverage.threshold")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.low.coverage.masking"], tags$code("filtering.low.coverage.masking"))
							),
							tags$td(
								checkboxInput("rnbOptsI.filtering.low.coverage.masking", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.low.coverage.masking")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.high.coverage.outliers"], tags$code("filtering.high.coverage.outliers"))
							),
							tags$td(
								checkboxInput("rnbOptsI.filtering.high.coverage.outliers", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.high.coverage.outliers")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.missing.value.quantile"], tags$code("filtering.missing.value.quantile"))
							),
							tags$td(
								sliderInput("rnbOptsI.filtering.missing.value.quantile", NULL, min=0, max=1, value=1)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.missing.value.quantile")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.greedycut"], tags$code("filtering.greedycut"))
							),
							tags$td(
								checkboxInput("rnbOptsI.filtering.greedycut", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.greedycut")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.sex.chromosomes.removal"], tags$code("filtering.sex.chromosomes.removal"))
							),
							tags$td(
								checkboxInput("rnbOptsI.filtering.sex.chromosomes.removal", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.sex.chromosomes.removal")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.snp"], tags$code("filtering.snp"))
							),
							tags$td(
								selectInput("rnbOptsI.filtering.snp", NULL, RNB.FILTERING.SNP, selected="3")
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.snp")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["filtering.cross.reactive"], tags$code("filtering.cross.reactive"))
							),
							tags$td(
								checkboxInput("rnbOptsI.filtering.cross.reactive", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.filtering.cross.reactive")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["normalization.method"], tags$code("normalization.method"))
							),
							tags$td(
								selectInput("rnbOptsI.normalization.method", NULL, RNB.NORMALIZATION.METHODS, selected="swan")
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.normalization.method")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["normalization.background.method"], tags$code("normalization.background.method"))
							),
							tags$td(
								selectInput("rnbOptsI.normalization.background.method", NULL, RNB.NORMALIZATION.BG.METHODS, selected="methylumi.noob")
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.normalization.background.method")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["imputation.method"], tags$code("imputation.method"))
							),
							tags$td(
								selectInput("rnbOptsI.imputation.method", NULL, RNB.IMPUTATION.METHODS)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.imputation.method")
							)
						)
					)
				)
			),
			tabPanel("Export",
				tags$table(class="table table-hover",
					tags$thead(tags$tr(
						tags$th("Name"),
						tags$th("Setting"),
						tags$th("Value")
					)),
					tags$tbody(
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["export.to.trackhub"], tags$code("export.to.trackhub"))
							),
							tags$td(
								selectInput("rnbOptsI.export.to.trackhub", NULL, RNB.TRACKHUB.FORMATS, multiple=TRUE, selected=RNB.TRACKHUB.FORMATS)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.export.to.trackhub", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["export.types"], tags$code("export.types"))
							),
							tags$td(
								uiOutput('selRegions.export')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.export.types", placeholder=TRUE)
							)
						)
					)
				)
			),
			tabPanel("Inference",
				tags$table(class="table table-hover",
					tags$thead(tags$tr(
						tags$th("Name"),
						tags$th("Setting"),
						tags$th("Value")
					)),
					tags$tbody(
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["inference"], tags$code("inference"))
							),
							tags$td(
								checkboxInput("rnbOptsI.inference", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.inference")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["inference.age.prediction"], tags$code("inference.age.prediction"))
							),
							tags$td(
								checkboxInput("rnbOptsI.inference.age.prediction", "Enable", value=TRUE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.inference.age.prediction")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["inference.targets.sva"], tags$code("inference.targets.sva"))
							),
							tags$td(
								uiOutput('selColumn.sva')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.inference.targets.sva", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["inference.sva.num.method"], tags$code("inference.sva.num.method"))
							),
							tags$td(
								selectInput("rnbOptsI.inference.sva.num.method", NULL, RNB.SVA.NUM.METHODS)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.inference.sva.num.method", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["inference.reference.methylome.column"], tags$code("inference.reference.methylome.column"))
							),
							tags$td(
								uiOutput('selColumn.cellTypeRef')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.inference.reference.methylome.column", placeholder=TRUE)
							)
						)
					)
				)
			),
			tabPanel("Exploratory",
				tags$table(class="table table-hover",
					tags$thead(tags$tr(
						tags$th("Name"),
						tags$th("Setting"),
						tags$th("Value")
					)),
					tags$tbody(
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["exploratory.columns"], tags$code("exploratory.columns"))
							),
							tags$td(
								uiOutput('selColumn.ex')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.exploratory.columns", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["exploratory.intersample"], tags$code("exploratory.intersample"))
							),
							tags$td(
								checkboxInput("rnbOptsI.exploratory.intersample", "Enable", value=TRUE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.exploratory.intersample")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["exploratory.beta.distribution"], tags$code("exploratory.beta.distribution"))
							),
							tags$td(
								checkboxInput("rnbOptsI.exploratory.beta.distribution", "Enable", value=TRUE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.exploratory.beta.distribution")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["exploratory.correlation.qc"], tags$code("exploratory.correlation.qc"))
							),
							tags$td(
								checkboxInput("rnbOptsI.exploratory.correlation.qc", "Enable", value=TRUE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.exploratory.correlation.qc")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["exploratory.region.profiles"], tags$code("exploratory.region.profiles"))
							),
							tags$td(
								uiOutput('selRegionProfiles.ex')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.exploratory.region.profiles", placeholder=TRUE)
							)
						)
					)
				)
			),
			tabPanel("Differential",
				tags$table(class="table table-hover",
					tags$thead(tags$tr(
						tags$th("Name"),
						tags$th("Setting"),
						tags$th("Value")
					)),
					tags$tbody(
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["differential.comparison.columns"], tags$code("differential.comparison.columns"))
							),
							tags$td(
								uiOutput('selColumn.diff')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.comparison.columns", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["covariate.adjustment.columns"], tags$code("covariate.adjustment.columns"))
							),
							tags$td(
								uiOutput('selAdjColumns')
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.covariate.adjustment.columns", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["differential.site.test.method"], tags$code("differential.site.test.method"))
							),
							tags$td(
								selectInput("rnbOptsI.differential.site.test.method", NULL, RNB.DIFFMETH.TEST.METHODS)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.site.test.method", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["differential.report.sites"], tags$code("differential.report.sites"))
							),
							tags$td(
								checkboxInput("rnbOptsI.differential.report.sites", "Enable", value=TRUE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.report.sites")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["differential.variability"], tags$code("differential.variability"))
							),
							tags$td(
								checkboxInput("rnbOptsI.differential.variability", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.variability")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["differential.variability.method"], tags$code("differential.variability.method"))
							),
							tags$td(
								selectInput("rnbOptsI.differential.variability.method", NULL, RNB.DIFFVAR.METHODS)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.variability.method")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["differential.enrichment.go"], tags$code("differential.enrichment.go"))
							),
							tags$td(
								checkboxInput("rnbOptsI.differential.enrichment.go", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.enrichment.go")
							)
						),
						tags$tr(
							tags$td(
								tags$div(title=RNB.OPTION.DESC["differential.enrichment.lola"], tags$code("differential.enrichment.lola"))
							),
							tags$td(
								checkboxInput("rnbOptsI.differential.enrichment.lola", "Enable", value=FALSE)
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.enrichment.lola")
							)
						)
					)
				)
			),
			tabPanel("R Object",
				verbatimTextOutput("rnbOpts")
			)
		)
	),
	tabPanel("Run", icon=icon("play"),
		uiOutput("runHint"),
		actionButton("runRnb", "Run Analysis", class="btn-primary"),
		textOutput("runStatusMsg")
	),
	tabPanel("Modules", icon=icon("cubes"),
		tabsetPanel(
			tabPanel("Data Import",
				sidebarPanel(
					tags$h1("Import from ..."),
					wellPanel(
						tags$h4("... new Dataset"),
						tags$p("Import using the inputs and options specified in other tabs."),
						uiOutput("modImportNew.about"),
						actionButton("modImportNew", "Import New Data", class="btn-primary")
					),
					wellPanel(
						tags$h4("... existing Analysis Directory"),
						uiOutput("modImportAnaDir.about"),
						actionButton("modImportAnaDir", "Import from Analysis Directory", class="btn-primary")
					),
					wellPanel(
						tags$h4("... RnBeads Objects"),
						directoryInput('modImportRnBSetDir', label = 'Select directory containing an RnBSet object', value = '~'),
						# localFileInput('modImportOptionsFile', label = 'Select Options XML file (optional)'),
						actionButton("modImportRObjects", "Import", class="btn-primary")
					),
					wellPanel(
						tags$h4("Reset"),
						actionButton("modImportReset", "Unload dataset", class="btn-primary")
					)
				),
				mainPanel(
					tags$h1("Loading Status"),
					uiOutput("modImport.status"),
					uiOutput("rnbSetInfo")
				)
			),
			tabPanel("Quality Control",
				sidebarPanel(
					checkboxInput("modQC.overwrite", "Overwrite existing report", value=FALSE),
					actionButton("modQC.run", "Run Quality Control", class="btn-primary")
				),
				mainPanel(
					uiOutput("modQC.status")
				)
			),
			tabPanel("...",
				"TODO"
			)
		)
	)
))
################################################################################
# Server configuration
################################################################################
server <- function(input, output, session) {
	############################################################################
	# SANDBOX
	# file and directory handling
	# shinyDirChoose(input, "outDirSel", roots=c(disk="/", home="~", wd="."), session=session)

	# dynamic outputs
	observeLocalFileInput(input, session, 'sandboxIn2')
	output$sandboxOut2 <- renderPrint({readLocalFileInput(session, 'sandboxIn2')})
	############################################################################

	# output$outDirSel <- renderPrint({input$outDirSel})
	# output$outDirSel <- renderPrint({readDirectoryInput(session, 'outDirSel')})


	############################################################################
	# Analysis status
	############################################################################
	refreshTimer <- reactiveTimer(REFRESH.TIME)
	observeDirectoryInput(input, session, 'outDir')
	reportDir <- reactive({
		file.path(readDirectoryInput(session, 'outDir'), input$reportSubDir)
	})
	anaStatus <- reactive({
		refreshTimer()
		res <- list(status="invalid", statusTab=NULL, rnbSet.paths=c())
		if (dir.exists(reportDir())){
			reportStatus <- checkReportDir(reportDir())
			if (reportStatus$valid){
				#Valid existing report directory
				statusTab <- data.frame(
					Module=names(RNB.MODULES),
					Report=NA,
					Status=NA
				)
				rownames(statusTab) <- RNB.MODULES
				for (mm in RNB.MODULES){
					if (!is.na(reportStatus$reportHtml[mm])){
						statusTab[mm, "Report"] <- as.character(tags$a(href=paste0("file://", file.path(reportDir(), reportStatus$reportHtml[mm])), tags$code("html")))
						if (reportStatus$moduleStatus[mm, "status"]=="completed"){
							statusTab[mm, "Status"] <- as.character(tags$span(style="color:green", icon("check")))
						} else if (reportStatus$moduleStatus[mm, "status"]=="started"){
							statusTab[mm, "Status"] <- as.character(tags$span(style="color:orange", icon("play")))
						}
						
					} else {
						statusTab[mm, "Status"] <- as.character(tags$span(style="color:red", icon("times")))
					}
				}
				res$rnbSet.paths <- file.path(reportDir(), grep("rnbSet_", list.dirs(reportDir(), full.names=FALSE, recursive=FALSE), value=TRUE))
				if (!modImportStatus$dataset.loaded) shinyjs::enable("modImportAnaDir")
				res$status <- "reportDir"
				res$statusTab <- statusTab
			} else {
				# Existing directory, but no valid RnBeads report
				shinyjs::disable("modImportAnaDir")
			}	
		} else {
			# New RnBeads analysis
			shinyjs::disable("modImportAnaDir")
			res$status <- "new"
		}
		res
	})
	output$anaStatus <- renderUI({
		curStatus <- anaStatus()
		if (curStatus$status=="new"){
			tagList(
				tags$h1("Configure new analysis"),
				tags$p(
					"Configure the input and analysis options by clicking on the tabs above",
					"and selecting corresponding values. Finally go to 'Run' to start your analysis."
				),
				tags$p(
					"A new report will be created in", tags$code(reportDir())
				)
			)
		} else if (curStatus$status=="reportDir"){
			tagList(
				tags$h1("Existing RnBeads analysis"),
				renderTable(curStatus$statusTab, sanitize.text.function=function(x){x})
			)
		} else if (curStatus$status=="invalid"){
			tagList(
				tags$h1("Invalid analysis directory"),
				tags$p("Please select a non-existing analysis directory for new analyses or an existing RnBeads report directory.")
			)
		}
	})

	############################################################################
	# Input
	############################################################################
	observeLocalFileInput(input, session, 'sampleAnnotFile')
	observeDirectoryInput(input, session, 'dataDir')

	tabSep <- reactive({input$rnbOptsI.import.table.separator})
	sampleAnnotFile <- reactive({
		readLocalFileInput(session, 'sampleAnnotFile')
	})
	inputDataDir <- reactive({
		readDirectoryInput(session, 'dataDir')
	})

	sannot.fromFile <- reactive({
		sampleAnnotFn <- sampleAnnotFile()
		if (is.null(sampleAnnotFn) || sampleAnnotFn=="") return(NULL)
		tryCatch(
			{
				RnBeads::read.sample.annotation(sampleAnnotFn, sep=tabSep())
			},
			error = function(err) {
				data.frame(Error="Unable to load table", Why=err$message)
			}
		)
	})
	sannot <- reactive({
		# first look if an RnBSet object has been loaded. If so, take the sample annotation from there
		# otherwise require the user to upload an annotation table via file input
		if (!is.null(rnbData$rnbSet)){
			pheno(rnbData$rnbSet)
		} else {
			sannot.fromFile()
		}
	})
	sannot.nsamples <- reactive({
		nrow(sannot())
	})
	sannot.cols <- reactive({
		cnames <- colnames(sannot())
		if (length(cnames)==2 && cnames==c("Error", "Why")) cnames <- NULL #error --> no column names
		cnames
	})
	sannot.cols.grps <- reactive({
		if(length(sannot.cols())>0) {
			names(rnb.sample.groups(sannot()))
		} else {
			NULL
		}
	})
	output$sampleAnnotContent <- renderTable({sannot.fromFile()})

	isBiseq <- reactive({
		res <- input$platform == "biseq"
		if (res) {
			shinyjs::enable("rnbOptsI.assembly")
			shinyjs::enable("rnbOptsI.import.bed.style")
			shinyjs::enable("rnbOptsI.filtering.coverage.threshold")
			shinyjs::enable("rnbOptsI.filtering.low.coverage.masking")
			shinyjs::enable("rnbOptsI.filtering.high.coverage.outliers")
			updateCheckboxInput(session, "rnbOptsI.filtering.greedycut", value=FALSE)
			shinyjs::disable("rnbOptsI.filtering.greedycut")
			shinyjs::disable("rnbOptsI.filtering.cross.reactive")
			shinyjs::disable("rnbOptsI.normalization.method")
			shinyjs::disable("rnbOptsI.normalization.background.method")
			shinyjs::disable("rnbOptsI.exploratory.correlation.qc")
		} else {
			shinyjs::disable("rnbOptsI.assembly")
			shinyjs::disable("rnbOptsI.import.bed.style")
			shinyjs::disable("rnbOptsI.filtering.coverage.threshold")
			shinyjs::disable("rnbOptsI.filtering.low.coverage.masking")
			shinyjs::disable("rnbOptsI.filtering.high.coverage.outliers")
			updateCheckboxInput(session, "rnbOptsI.filtering.greedycut", value=TRUE)
			shinyjs::enable("rnbOptsI.filtering.greedycut")
			shinyjs::enable("rnbOptsI.filtering.cross.reactive")
			shinyjs::enable("rnbOptsI.normalization.method")
			shinyjs::enable("rnbOptsI.normalization.background.method")
			shinyjs::enable("rnbOptsI.exploratory.correlation.qc")
		}
		res
	})

	output$selColumn.id <- renderUI({
		selCols <- c("[None]", sannot.cols())
		selectInput('rnbOptsI.identifiers.column', NULL, selCols)
	})
	output$selColumn.sva<- renderUI({
		selCols <- c(sannot.cols.grps())
		selectInput('rnbOptsI.inference.targets.sva', NULL, selCols, multiple=TRUE)
	})
	output$selColumn.cellTypeRef <- renderUI({
		selCols <- c("[None]", sannot.cols())
		selectInput('rnbOptsI.inference.reference.methylome.column', NULL, selCols)
	})
	output$selColumn.ex <- renderUI({
		selCols <- c("[automatic]", sannot.cols.grps())
		selectInput('rnbOptsI.exploratory.columns', NULL, selCols, multiple=TRUE, selected="[automatic]")
	})
	output$selColumn.diff <- renderUI({
		selCols <- c("[automatic]", sannot.cols.grps())
		selectInput('rnbOptsI.differential.comparison.columns', NULL, selCols, multiple=TRUE, selected="[automatic]")
	})
	output$selAdjColumns <- renderUI({
		selCols <- sannot.cols()
		selectInput('rnbOptsI.covariate.adjustment.columns', NULL, selCols, multiple=TRUE, selected=c())
	})

	output$inputStatus <- renderUI({
		res <- list()
		showAnnotTab <- TRUE
		sampleAnnotFn <- sampleAnnotFile()
		if (is.null(sampleAnnotFn) || sampleAnnotFn=="") showAnnotTab <- FALSE
		if (showAnnotTab){
			res <- c(res,
				list(tags$h1("Preview of the sample annotation table")),
				list(tableOutput("sampleAnnotContent"))
			)
		} else {
			res <- c(res,
				list(tags$p("No sample annotation loaded"))
			)
		}
		tagList(res)
	})

	############################################################################
	# Analysis Options
	############################################################################
	# autoOpts <- c("analysis.name", "import.table.separator")
	# for (oo in autoOpts){
	# 	output[[paste0("rnbOptsO.", oo)]] <- renderText({
	# 		args <- list()
	# 		args[[oo]] <- input[[paste0("rnbOptsI.", oo)]]
	# 		do.call("rnb.options", args)
	# 		rnb.getOption(oo)
	# 	})
	# }

	output$rnbOpts <- renderPrint({
		depDummy <- assemblySel()
		# auto update on any option change to rnbOptsI
		for (oo in grep("^rnbOptsI.", names(input), value=TRUE)){
			input[[oo]]
		}
		rnb.options()
	})
	output$rnbOptsO.analysis.name <- renderText({
		rnb.options(analysis.name=input$rnbOptsI.analysis.name)
		rnb.getOption("analysis.name")
	})
	#default setting for assembly if platform is not bisulfite and thus no assembly input has been given
	assemblySel <- reactive({
		# print("HI")
		res <- "hg19"
		if (isBiseq() && !is.null(input$rnbOptsI.assembly)){
			res <- input$rnbOptsI.assembly
		}
		rnb.options(assembly=res)
		res
	})
	output$rnbOptsO.assembly <- renderText({
		depDummy <- assemblySel()
		rnb.getOption("assembly")
	})
	regTypes.all <- reactive({
		depDummy <- assemblySel() # dummy to update on dependency: assembly
		rnb.region.types(rnb.getOption("assembly"))
	})
	output$selRegionTypes <- renderUI({
		selectInput('rnbOptsI.region.types', NULL, regTypes.all(), multiple=TRUE, selected=regTypes.all())
	})
	regTypes <- reactive({
		input$rnbOptsI.region.types
	})
	output$rnbOptsO.region.types <- renderText({
		rnb.options(region.types=regTypes())
		rnb.getOption("region.types")
	})
	output$selRegionProfiles.ex <- renderUI({
		selRegs <- regTypes()
		defaultRegs <- intersect(c("genes", "promoters", "cpgislands"), selRegs)
		selectInput('rnbOptsI.exploratory.region.profiles', NULL, selRegs, selected=defaultRegs, multiple=TRUE)
	})
	output$selRegions.export <- renderUI({
		selRegs <- c("sites", regTypes())
		selectInput('rnbOptsI.export.types', NULL, selRegs, selected="sites", multiple=TRUE)
	})
	output$rnbOptsO.identifiers.column <- renderText({
		cname <- input$rnbOptsI.identifiers.column
		if (is.null(cname) || cname=="[None]") cname <- NULL
		rnb.options(identifiers.column=cname)
		rnb.getOption("identifiers.column")
	})
	output$rnbOptsO.min.group.size <- renderText({
		rnb.options(min.group.size=input$rnbOptsI.min.group.size)
		rnb.getOption("min.group.size")
	})
	observe({
		val.min.group.size <- input$rnbOptsI.min.group.size
		val.max.group.count <- input$rnbOptsI.max.group.count
		maxval <- 20
		if (!is.null(sannot.nsamples()) && sannot.nsamples() > 0) maxval <- sannot.nsamples()
		updateSliderInput(session, "rnbOptsI.min.group.size", value=val.min.group.size, min=1, max=maxval, step=1)
		updateSliderInput(session, "rnbOptsI.max.group.count", value=val.max.group.count, min=2, max=maxval, step=1)
	})
	output$rnbOptsO.max.group.count <- renderText({
		rnb.options(max.group.count=input$rnbOptsI.max.group.count)
		rnb.getOption("max.group.count")
	})
	output$rnbOptsO.colors.category <- renderText({
		cols <- RNB.COLSCHEMES.CATEGORY[[input$rnbOptsI.colors.category]]
		rnb.options(colors.category=cols)
		rnb.getOption("colors.category")
	})
	output$rnbOptsOP.colors.category <- renderPlot({
		cols <- RNB.COLSCHEMES.CATEGORY[[input$rnbOptsI.colors.category]]
		# rnb.options(colors.category=cols)
		# rnb.getOption("colors.category")
		plotColPal(cols)
	})
	output$rnbOptsO.colors.meth <- renderText({
		cols <- RNB.COLSCHEMES.METH[[input$rnbOptsI.colors.meth]]
		rnb.options(colors.meth=cols)
		rnb.getOption("colors.meth")
	})
	output$rnbOptsOP.colors.meth <- renderPlot({
		cols <- RNB.COLSCHEMES.METH[[input$rnbOptsI.colors.meth]]
		# rnb.options(colors.meth=cols)
		# rnb.getOption("colors.meth")
		plotColPal(cols)
	})
	output$rnbOptsO.import.default.data.type <- renderText({
		dt <- "bed.dir"
		if (!isBiseq()) dt <- "idat.dir"
		rnb.options(import.default.data.type=dt)
		rnb.getOption("import.default.data.type")
	})
	output$rnbOptsO.import.table.separator <- renderText({
		rnb.options(import.table.separator=tabSep())
		rnb.getOption("import.table.separator")
	})
	output$rnbOptsO.import.bed.style <- renderText({
		res <- rnb.getOption("import.bed.style")
		if (isBiseq()){
			res <- input$rnbOptsI.import.bed.style
			rnb.options(import.bed.style=res)
		}
		res
	})
	output$rnbOptsO.filtering.coverage.threshold <- renderText({
		res <- rnb.getOption("filtering.coverage.threshold")
		if (isBiseq()){
			res <- input$rnbOptsI.filtering.coverage.threshold
			rnb.options(filtering.coverage.threshold=res)
		}
		res
	})
	output$rnbOptsO.filtering.low.coverage.masking <- renderText({
		res <- rnb.getOption("filtering.low.coverage.masking")
		if (isBiseq()){
			res <- input$rnbOptsI.filtering.low.coverage.masking
			rnb.options(filtering.low.coverage.masking=res)
		}
		res
	})
	output$rnbOptsO.filtering.high.coverage.outliers <- renderText({
		res <- rnb.getOption("filtering.high.coverage.outliers")
		if (isBiseq()){
			res <- input$rnbOptsI.filtering.high.coverage.outliers
			rnb.options(filtering.high.coverage.outliers=res)
		}
		res
	})
	output$rnbOptsO.filtering.missing.value.quantile <- renderText({
		rnb.options(filtering.missing.value.quantile=input$rnbOptsI.filtering.missing.value.quantile)
		rnb.getOption("filtering.missing.value.quantile")
	})
	output$rnbOptsO.filtering.greedycut <- renderText({
		rnb.options(filtering.greedycut=input$rnbOptsI.filtering.greedycut)
		rnb.getOption("filtering.greedycut")
	})
	output$rnbOptsO.filtering.sex.chromosomes.removal <- renderText({
		rnb.options(filtering.sex.chromosomes.removal=input$rnbOptsI.filtering.sex.chromosomes.removal)
		rnb.getOption("filtering.sex.chromosomes.removal")
	})
	output$rnbOptsO.filtering.snp <- renderText({
		rnb.options(filtering.snp=input$rnbOptsI.filtering.snp)
		rnb.getOption("filtering.snp")
	})
	output$rnbOptsO.filtering.cross.reactive <- renderText({
		rnb.options(filtering.cross.reactive=input$rnbOptsI.filtering.cross.reactive)
		rnb.getOption("filtering.cross.reactive")
	})
	output$rnbOptsO.normalization.method <- renderText({
		res <- rnb.getOption("normalization.method")
		if (!isBiseq()){
			res <- input$rnbOptsI.normalization.method
			rnb.options(normalization.method=res)
		}
		res
	})
	output$rnbOptsO.normalization.background.method <- renderText({
		res <- rnb.getOption("normalization.background.method")
		if (!isBiseq()){
			res <- input$rnbOptsI.normalization.background.method
			rnb.options(normalization.background.method=res)
		}
		res
	})
	output$rnbOptsO.imputation.method <- renderText({
		rnb.options(imputation.method=input$rnbOptsI.imputation.method)
		rnb.getOption("imputation.method")
	})
	output$rnbOptsO.export.to.trackhub <- renderText({
		rnb.options(export.to.trackhub=input$rnbOptsI.export.to.trackhub)
		rnb.getOption("export.to.trackhub")
	})
	output$rnbOptsO.export.types <- renderText({
		rts <- input$rnbOptsI.export.types
		if (length(rts)<1) rts <- NULL
		rnb.options(export.types=rts)
		rnb.getOption("export.types")
	})
	doInference <- reactive({
		res <- input$rnbOptsI.inference
		inferenceOptNames <- grep("^rnbOptsI.inference.", names(input), value=TRUE)
		if (res){
			for (oo in inferenceOptNames){
				shinyjs::enable(oo)
			}
		} else {
			for (oo in inferenceOptNames){
				shinyjs::disable(oo)
			}
		}
		res
	})
	output$rnbOptsO.inference <- renderText({
		rnb.options(inference=doInference())
		rnb.getOption("inference")
	})
	output$rnbOptsO.inference.age.prediction <- renderText({
		rnb.options(inference.age.prediction=input$rnbOptsI.inference.age.prediction)
		rnb.getOption("inference.age.prediction")
	})
	output$rnbOptsO.inference.targets.sva <- renderText({
		cnames <- input$rnbOptsI.inference.targets.sva
		if (length(cnames)<1) cnames <- character(0)
		rnb.options(inference.targets.sva=cnames)
		rnb.getOption("inference.targets.sva")
	})
	output$rnbOptsO.inference.sva.num.method <- renderText({
		rnb.options(inference.sva.num.method=input$rnbOptsI.inference.sva.num.method)
		rnb.getOption("inference.sva.num.method")
	})
	output$rnbOptsO.inference.reference.methylome.column <- renderText({
		cname <- input$rnbOptsI.inference.reference.methylome.column
		if (is.null(cname) || cname=="[None]") cname <- NULL
		rnb.options(inference.reference.methylome.column=cname)
		rnb.getOption("inference.reference.methylome.column")
	})
	output$rnbOptsO.exploratory.columns <- renderText({
		cnames <- input$rnbOptsI.exploratory.columns
		if (length(cnames)<1) cnames <- NULL
		if (length(cnames)==1 && cnames=="[automatic]") cnames <- NULL
		cnames <- setdiff(cnames, "[automatic]")
		rnb.options(exploratory.columns=cnames)
		rnb.getOption("exploratory.columns")
	})
	output$rnbOptsO.exploratory.intersample <- renderText({
		rnb.options(exploratory.intersample=input$rnbOptsI.exploratory.intersample)
		rnb.getOption("exploratory.intersample")
	})
	output$rnbOptsO.exploratory.beta.distribution <- renderText({
		rnb.options(exploratory.beta.distribution=input$rnbOptsI.exploratory.beta.distribution)
		rnb.getOption("exploratory.beta.distribution")
	})
	output$rnbOptsO.exploratory.correlation.qc <- renderText({
		res <- rnb.getOption("exploratory.correlation.qc")
		if (!isBiseq()){
			res <- input$rnbOptsI.exploratory.correlation.qc
			rnb.options(exploratory.correlation.qc=res)
		}
		res
	})
	output$rnbOptsO.exploratory.region.profiles <- renderText({
		rts <- input$rnbOptsI.exploratory.region.profiles
		if (length(rts)<1) rts <- character(0)
		rnb.options(exploratory.region.profiles=rts)
		rnb.getOption("exploratory.region.profiles")
	})
	output$rnbOptsO.differential.comparison.columns <- renderText({
		cnames <- input$rnbOptsI.differential.comparison.columns
		if (length(cnames)<1) cnames <- NULL
		if (length(cnames)==1 && cnames=="[automatic]") cnames <- NULL
		cnames <- setdiff(cnames, "[automatic]")
		rnb.options(differential.comparison.columns=cnames)
		rnb.getOption("differential.comparison.columns")
	})
	output$rnbOptsO.covariate.adjustment.columns <- renderText({
		cnames <- input$rnbOptsI.covariate.adjustment.columns
		if (length(cnames)<1) cnames <- NULL
		rnb.options(covariate.adjustment.columns=cnames)
		rnb.getOption("covariate.adjustment.columns")
	})
	output$rnbOptsO.differential.site.test.method <- renderText({
		rnb.options(differential.site.test.method=input$rnbOptsI.differential.site.test.method)
		rnb.getOption("differential.site.test.method")
	})
	output$rnbOptsO.differential.report.sites <- renderText({
		rnb.options(differential.report.sites=input$rnbOptsI.differential.report.sites)
		rnb.getOption("differential.report.sites")
	})
	doDiffVar <- reactive({
		res <- input$rnbOptsI.differential.variability
		if (res){
			shinyjs::enable("rnbOptsI.differential.variability.method")
		} else {
			shinyjs::disable("rnbOptsI.differential.variability.method")
		}
		res
	})
	output$rnbOptsO.differential.variability <- renderText({
		rnb.options(differential.variability=doDiffVar())
		rnb.getOption("differential.variability")
	})
	output$rnbOptsO.differential.variability.method <- renderText({
		rnb.options(differential.variability.method=input$rnbOptsI.differential.variability.method)
		rnb.getOption("differential.variability.method")
	})
	output$rnbOptsO.differential.enrichment.go <- renderText({
		rnb.options(differential.enrichment.go=input$rnbOptsI.differential.enrichment.go)
		rnb.getOption("differential.enrichment.go")
	})
	output$rnbOptsO.differential.enrichment.lola <- renderText({
		rnb.options(differential.enrichment.lola=input$rnbOptsI.differential.enrichment.lola)
		rnb.getOption("differential.enrichment.lola")
	})
	############################################################################
	# Run
	############################################################################
	enableRun <- reactive({
		res <- TRUE
		if (dir.exists(reportDir())) res <- FALSE #directory exists
		if (!file.exists(sampleAnnotFile())) res <- FALSE #sample annotation file does not exist
		if (!dir.exists(inputDataDir())) res <- FALSE # data directory does not exist
		return(res)
	})
	output$runHint <- renderUI({
		if (enableRun()){
			shinyjs::enable("runRnb")
			tags$p(tags$span(style="color:green", icon("play"), "Ready to run"))
		} else {
			shinyjs::disable("runRnb")
			tags$p(tags$span(style="color:red", icon("warning"), "Unable to run the analysis. Please make sure that you selected a non-existing report directory ('Analysis' tab) and that you specified the sample annotation file and data directory correctly ('Input' tab)"))
		}
	})
	observeEvent(input$numCores, {
		if (input$numCores > 1){
			parallel.setup(input$numCores)
		} else {
			parallel.teardown()
		}
		logger.close()
	})
	observeEvent(input$ggplotTheme, {
		eval(parse(text=paste0("theme_set(theme_", input$ggplotTheme, "())")))
	})
	observeEvent(input$runRnb, {
		if(input$runRnb == 0) return()
        shinyjs::disable("runRnb")

    	logger.close() # close the server's current logger s.t. a new log file is created in the report directory
		withProgress({
			rnb.run.analysis(
				dir.reports=reportDir(),
				sample.sheet=sampleAnnotFile(),
				data.dir=inputDataDir(),
				data.type=rnb.getOption("import.default.data.type"),
				save.rdata=TRUE
			)
		}, message="Runnning RnBeads analysis")
	})
	# output$runStatusMsg <- renderText({
	# 	refreshTimer()
	# 	print(input$runRnb)
	# 	if(input$runRnb == 0) return("Analysis not started")
	# 	res <- "Analysis started"
	# 	# logFn <- file.path(reportDir(), "analysis.log")
	# 	# status <- getRnbStatusFromLog(logFn)
	# 	# status <- status[status$scheduled,]
	# 	# if (all(status$status=="completed")) return("Analysis completed")
	# 	# startedModules <- which(status$status=="started")
	# 	# if (length(startedModules) > 0){
	# 	# 	mm <- status$module[startedModules[length(startedModules)]] # pick last started module
	# 	# 	res <- paste0("Running ", mm)
	# 	# }
	# 	return(res)
	# })

	############################################################################
	# Modules
	############################################################################
	# IMPORT
	output$modImportNew.about <- renderUI({
		if (enableRun()){
			if (!modImportStatus$dataset.loaded){
				shinyjs::enable("modImportNew")
				shinyjs::enable("modImportNew.save")
			}
			tagList(
				tags$p(tags$span(icon("play"), "Ready to run")),
				checkboxInput("modImportNew.save", "Save RnBSet object", value=TRUE)
			)
		} else {
			shinyjs::disable("modImportNew")
			shinyjs::disable("modImportNew.save")
			tags$p(tags$span(icon("warning"), "Unable to run. Please make sure that you selected a non-existing report directory ('Analysis' tab) and that you specified the sample annotation file and data directory correctly ('Input' tab)"))
		}
	})
	output$modImportAnaDir.about <- renderUI({
		curStatus <- anaStatus()
		if (curStatus$status=="reportDir"){
			if (length(curStatus$rnbSet.paths) > 0){
				pp <- curStatus$rnbSet.paths
				names(pp) <- basename(pp)
				tagList(
					tags$p("An existing RnBeads report has been specified. Ready to import."),
					selectInput("modImportAnaDir.rnbSet", "Select RnBSet to load", pp)
				)
			} else {
				tagList(tags$p(icon("warning"), "The analysis directory does not contain RnBSet objects."))
			}
		} else {
			tagList(tags$p(icon("warning"), "Please specify an existing RnBeads report in the 'Analysis' tab to continue."))
		}
	})
	output$rnbSetInfo <- renderUI({
		if (is.null(rnbData$rnbSet)) return(NULL)
		tagList(
			tags$h2("RnBSet Object:"),
			renderPrint({methods::show(rnbData$rnbSet)}),
			tags$h2("Sample Annotation:"),
			renderTable({sannot()})
		)
	})
	modImportStatus <- reactiveValues(
		dataset.loaded=FALSE,
		dataset.loaded.nsamples=-1
	)
	rnbData <- reactiveValues(
		rnbSet=NULL
	)
	# watch if dataset has been loaded
	reactDataset <- observeEvent(modImportStatus$dataset.loaded, {
		if (modImportStatus$dataset.loaded){
			shinyjs::disable("modImportNew")
			shinyjs::disable("modImportNew.save")
			shinyjs::disable("modImportAnaDir")
			shinyjs::disable("modImportRObjects")
			shinyjs::enable("modImportReset")
		} else {
			if (enableRun()) {
				shinyjs::enable("modImportNew")
				shinyjs::enable("modImportNew.save")
			}
			if (anaStatus()$status=="reportDir" && length(anaStatus()$rnbSet.paths) > 0) shinyjs::enable("modImportAnaDir")
			shinyjs::enable("modImportRObjects")
			shinyjs::disable("modImportReset")
		}
	})
	output$modImport.status <- renderUI({
		# curStatus <- modImportStatus()
		res <- list()
		if (modImportStatus$dataset.loaded){
			res <- c(res, list(tags$p(tags$span(style="color:green", icon("check"), "Loaded dataset containing", modImportStatus$dataset.loaded.nsamples, "samples"))))
		} else {
			res <- c(res, list(tags$p(tags$span(style="color:red", icon("times"), "No dataset loaded"))))
		}
		res <- c(res, list(tags$p("The current report directory is ", tags$code(reportDir()), ". This can be configured in the 'Analysis' tab.")))
		tagList(res)
	})
	observeEvent(input$modImportNew, {
		modImportStatus$dataset.loaded <- TRUE
		withProgress({
			tryCatch({
					if(!dir.exists(reportDir())) rnb.initialize.reports(reportDir())
					logger.start(fname=file.path(reportDir(), "analysis.log"))
					res <- rnb.run.import(c(inputDataDir(), sampleAnnotFile()), data.type=rnb.getOption("import.default.data.type"), dir.reports=reportDir())
					logger.close()
					rnbData$rnbSet <- res$rnb.set
					markDirDJ(reportDir())
					if (input$modImportNew.save) save.rnb.set(rnbData$rnbSet, file.path(reportDir(), "rnbSet_import"), archive=FALSE)
					modImportStatus$dataset.loaded.nsamples <- length(samples(rnbData$rnbSet))
				},
				error = function(err) {
					rnbData$rnbSet <- NULL
					modImportStatus$dataset.loaded <- FALSE
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Failed to import data: ", err$message)))
				}
			)
		}, message="Importing dataset")
	})
	observeEvent(input$modImportAnaDir, {
		modImportStatus$dataset.loaded <- TRUE
		withProgress({
			tryCatch({
					rnbData$rnbSet <- load.rnb.set(input$modImportAnaDir.rnbSet)
					modImportStatus$dataset.loaded.nsamples <- length(samples(rnbData$rnbSet))
				},
				error = function(err) {
					rnbData$rnbSet <- NULL
					modImportStatus$dataset.loaded <- FALSE
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Failed to load RnBSet object: ", err$message)))
				}
			)
		}, message="Loading dataset")
		#TODO: import options
	})
	observeDirectoryInput(input, session, 'modImportRnBSetDir')
	# observeLocalFileInput(input, session, 'modImportOptionsFile')
	observeEvent(input$modImportRObjects, {
		modImportStatus$dataset.loaded <- TRUE
		withProgress({
			tryCatch({
					rnbData$rnbSet <- load.rnb.set(readDirectoryInput(session, 'modImportRnBSetDir'))
					modImportStatus$dataset.loaded.nsamples <- length(samples(rnbData$rnbSet))
				},
				error = function(err) {
					rnbData$rnbSet <- NULL
					modImportStatus$dataset.loaded <- FALSE
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Failed to load RnBSet object: ", err$message)))
				}
			)
		}, message="Loading dataset")
		#TODO: import options
	})
	observeEvent(input$modImportReset, {
		modImportStatus$dataset.loaded <- FALSE
		modImportStatus$dataset.loaded.nsamples <- -1
		rnbData$rnbSet <- NULL
	})

	# QUALITY CONTROL
	reportExists.quality_control <- reactive({
		curStatus <- anaStatus()
		if (is.null(curStatus$statusTab)) return(FALSE)
		!is.na(curStatus$statusTab["quality_control", "Report"])
	})
	output$modQC.status <- renderUI({
		reportStatus <- checkReportDir(reportDir())
		rdy <- FALSE
		res <- list()
		if (modImportStatus$dataset.loaded){
			rdy <- TRUE
			res <- c(res, list(tags$p(tags$span(style="color:green", icon("check"), "Dataset loaded."))))
		} else {
			res <- c(res, list(tags$p(tags$span(style="color:red", icon("times"), "No dataset loaded. Please load a dataset using the 'Data Import' tab."))))
		}
		res <- c(res, list(tags$p("The current report directory is ", tags$code(reportDir()), ". This can be configured in the 'Analysis' tab.")))
		if (reportExists.quality_control()){
			if (rdy){
				shinyjs::enable("modQC.overwrite")
			}
			if (!input$modQC.overwrite){
				rdy <- FALSE
			}
			res <- c(res, list(tags$div(title="If the link does not open use right-click, copy the link and paste the address in a new browser window", tags$p(tags$span(icon("search"), tags$a(href=paste0("file://", file.path(reportDir(), reportStatus$reportHtml["quality_control"])), "View report"))))))
		} else {
			shinyjs::disable("modQC.overwrite")
		}
		if (rdy){
			shinyjs::enable("modQC.run")
		} else {
			shinyjs::disable("modQC.run")
		}
		tagList(res)
	})
	observeEvent(input$modQC.run, {
		withProgress({
			tryCatch({
					if(!dir.exists(reportDir())) rnb.initialize.reports(reportDir())
					if (reportExists.quality_control() && input$modQC.overwrite){
						unlink(file.path(reportDir(), paste0("quality_control","*")), recursive=TRUE)
						showNotification(tags$span(icon("trash"), paste0("Deleted previous report")))
					}
					updateCheckboxInput(session, "modQC.overwrite", value=FALSE)
					logger.start(fname=file.path(reportDir(), "analysis.log"))
					res <- rnb.run.qc(rnbData$rnbSet, dir.reports=reportDir())
					logger.close()
					markDirDJ(reportDir())
					showNotification(tags$span(style="color:green", icon("check"), paste0("Analysis (Quality Control) completed")))
				},
				error = function(err) {
					showNotification(tags$span(style="color:red", icon("warning"), paste0("Analysis (Quality Control) failed:", err$message)))
				}
			)
		}, message="Performing analysis: Quality Control")
	})
}

################################################################################
# Main
################################################################################
shinyApp(ui = ui, server = server)

################################################################################
# Sandbox
################################################################################
# useful link for themes/layout
# https://shiny.rstudio.com/gallery/shiny-theme-selector.html

################################################################################
# TODOs:
# - help page parser for tooltips for options (is '?parse_Rd' useful?)
################################################################################

