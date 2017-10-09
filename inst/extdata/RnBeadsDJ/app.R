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
RNB.NORMALIZATION.METHODS=c("none", "bmiq", "illumina", "swan", "minfi.funnorm", "wm.dasen", "wm.nasen", "wm.betaqn", "wm.naten", "wm.nanet", "wm.nanes", "wm.danes", "wm.danet", "wm.danen", "wm.daten1", "wm.daten2", "wm.tost", "wm.fuks", "wm.swan")
RNB.NORMALIZATION.BG.METHODS <- c("none", "methylumi.noob", "methylumi.goob", "enmix.oob")
RNB.COLSCHEMES.CATEGORY <- list(
	default=c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"),
	extended=c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#2166AC","#B2182B","#00441B","#40004B","#053061","#003D7C","#D50911")
)
RNB.COLSCHEMES.METH <- list(
	default=c("#AD0021","#909090","#39278C"),
	YlBl=c("#EDF8B1","#41B6C4","#081D58")
)

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
	choose.dir <- function(default = NA, caption = NA) {
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

	choose.files <- function(default = NA, caption = NA) {
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
				path = choose.files(default = readLocalFileInput(session, inputId))
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

	if (!dir.exists(repDir)) return(res)
	contents <- list.files(repDir, include.dirs=TRUE)
	res$valid <- all(c("configuration", "index.html", "analysis_options.xml", "analysis.log") %in% contents)
	if (!res$valid) return(res)
	htmlModules <- paste0(RNB.MODULES, ".html")
	htmlModules.exist <- htmlModules %in% contents
	res$reportHtml[htmlModules.exist] <- htmlModules[htmlModules.exist]
	res$moduleStatus <- getRnbStatusFromLog(file.path(repDir, "analysis.log"))
	return(res)
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
			textInput("reportSubDir", "Choose the name of the report directory", "rnbeads_report")
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
			tags$h1("Preview of the sample annotation table"),
			tableOutput("sampleAnnotContent")
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
								tags$code("analysis.name")
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
								tags$code("assembly")
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
								tags$code("region.types")
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
								tags$code("identifiers.column")
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
								tags$code("min.group.size")
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
								tags$code("max.group.count")
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
								tags$code("colors.category")
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
								tags$code("colors.meth")
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
								tags$code("import.default.data.type")
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
								tags$code("import.table.separator")
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
								tags$code("import.bed.style")
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
								tags$code("filtering.coverage.threshold")
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
								tags$code("filtering.sex.chromosomes.removal")
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
								tags$code("normalization.method")
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
								tags$code("normalization.background.method")
							),
							tags$td(
								selectInput("rnbOptsI.normalization.background.method", NULL, RNB.NORMALIZATION.BG.METHODS, selected="methylumi.noob")
							),
							tags$td(
								verbatimTextOutput("rnbOptsO.normalization.background.method")
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
								tags$code("exploratory.columns")
							),
							tags$td(
								uiOutput('selColumn.ex')
							)
							,
							tags$td(
								verbatimTextOutput("rnbOptsO.exploratory.columns", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$code("exploratory.intersample")
							),
							tags$td(
								checkboxInput("rnbOptsI.exploratory.intersample", "Enable", value=TRUE)
							)
							,
							tags$td(
								verbatimTextOutput("rnbOptsO.exploratory.intersample")
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
								tags$code("differential.comparison.columns")
							),
							tags$td(
								uiOutput('selColumn.diff')
							)
							,
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.comparison.columns", placeholder=TRUE)
							)
						),
						tags$tr(
							tags$td(
								tags$code("differential.report.sites")
							),
							tags$td(
								checkboxInput("rnbOptsI.differential.report.sites", "Enable", value=TRUE)
							)
							,
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.report.sites")
							)
						),
						tags$tr(
							tags$td(
								tags$code("differential.enrichment.go")
							),
							tags$td(
								checkboxInput("rnbOptsI.differential.enrichment.go", "Enable", value=FALSE)
							)
							,
							tags$td(
								verbatimTextOutput("rnbOptsO.differential.enrichment.go")
							)
						),
						tags$tr(
							tags$td(
								tags$code("differential.enrichment.lola")
							),
							tags$td(
								checkboxInput("rnbOptsI.differential.enrichment.lola", "Enable", value=FALSE)
							)
							,
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
		tags$p(textOutput("runHint")),
		sliderInput('numCores', "Select the number of cores to use", min=1, max=detectCores(), value=1, step=1),
		selectInput('ggplotTheme', "Select a theme for the plots", c("Black & White"="bw", "Grey"="grey")),
		actionButton("runRnb", "Run Analysis", class="btn-primary"),
		textOutput("runStatusMsg")
	),
	tabPanel("Modules", icon=icon("cubes"),
		tabsetPanel(
			tabPanel("Data Import",
				sidebarPanel(
					wellPanel(
						tags$h4("New Dataset"),
						actionButton("modImportNew", "Import New Data", class="btn-primary")
					),
					wellPanel(
						tags$h4("Existing Analysis Directory"),
						uiOutput("modImportAnaDir.about"),
						actionButton("modImportAnaDir", "Import from Analysis Directory", class="btn-primary")
					),
					wellPanel(
						tags$h4("RnBeads Objects"),
						directoryInput('modImportRnBSetDir', label = 'Select directory containing an RnBSet object', value = '~'),
						localFileInput('modImportOptionsFile', label = 'Select Options XML file (optional)'),
						actionButton("modImportRObjects", "Import", class="btn-primary")
					),
					wellPanel(
						tags$h4("Reset"),
						actionButton("modImportReset", "Unload dataset", class="btn-primary")
					)
				),
				mainPanel(
					tags$h1("Loading Status"),
					uiOutput("modImport.status")
				)
			),
			tabPanel("Quality Control",
				"Blubb"
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
	refreshTimer <- reactiveTimer(5000)
	observeDirectoryInput(input, session, 'outDir')
	reportDir <- reactive({
		file.path(readDirectoryInput(session, 'outDir'), input$reportSubDir)
	})
	anaStatus <- reactive({
		res <- list(status="invalid", statusTab=NULL)
		if (dir.exists(reportDir())){
			reportStatus <- checkReportDir(reportDir())
			if (reportStatus$valid){
				#Valid existing report directory
				statusTab <- data.frame(
					Module=names(RNB.MODULES),
					Report=NA,
					Status=NA
				)
				for (i in 1:length(RNB.MODULES)){
					mm <- RNB.MODULES[i]
					if (!is.na(reportStatus$reportHtml[mm])){
						statusTab$Report[i] <- as.character(tags$a(href=paste0("file://", file.path(reportDir(), reportStatus$reportHtml[mm])), tags$code("html")))
						if (reportStatus$moduleStatus[mm, "status"]=="completed"){
							statusTab$Status[i] <- as.character(tags$span(style="color:green", icon("check")))
						} else if (reportStatus$moduleStatus[mm, "status"]=="started"){
							statusTab$Status[i] <- as.character(tags$span(style="color:orange", icon("play")))
						}
						
					} else {
						statusTab$Status[i] <- as.character(tags$span(style="color:red", icon("times")))
					}
				}
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

	sannot <- reactive({
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
	output$sampleAnnotContent <- renderTable({sannot()})

	isBiseq <- reactive({
		res <- input$platform == "biseq"
		if (res) {
			shinyjs::enable("rnbOptsI.assembly")
			shinyjs::enable("rnbOptsI.import.bed.style")
			shinyjs::enable("rnbOptsI.filtering.coverage.threshold")
			shinyjs::disable("rnbOptsI.normalization.method")
			shinyjs::disable("rnbOptsI.normalization.background.method")
		} else {
			shinyjs::disable("rnbOptsI.assembly")
			shinyjs::disable("rnbOptsI.import.bed.style")
			shinyjs::disable("rnbOptsI.filtering.coverage.threshold")
			shinyjs::enable("rnbOptsI.normalization.method")
			shinyjs::enable("rnbOptsI.normalization.background.method")
		}
		res
	})

	output$selColumn.id <- renderUI({
		selCols <- c("[None]", sannot.cols())
		selectInput('rnbOptsI.identifiers.column', NULL, selCols)
	})
	output$selColumn.ex <- renderUI({
		selCols <- c("[automatic]", sannot.cols.grps())
		selectInput('rnbOptsI.exploratory.columns', NULL, selCols, multiple=TRUE)
	})
	output$selColumn.diff <- renderUI({
		selCols <- c("[automatic]", sannot.cols.grps())
		selectInput('rnbOptsI.differential.comparison.columns', NULL, selCols, multiple=TRUE)
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
	regTypes <- reactive({
		depDummy <- assemblySel() # dummy to update on dependency: assembly
		rnb.region.types(rnb.getOption("assembly"))
	})
	output$selRegionTypes <- renderUI({
		selectInput('rnbOptsI.region.types', NULL, regTypes(), multiple=TRUE, selected=regTypes())
	})
	output$rnbOptsO.region.types <- renderText({
		rnb.options(region.types=regTypes())
		rnb.getOption("region.types")
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
	output$rnbOptsO.filtering.sex.chromosomes.removal <- renderText({
		rnb.options(filtering.sex.chromosomes.removal=input$rnbOptsI.filtering.sex.chromosomes.removal)
		rnb.getOption("filtering.sex.chromosomes.removal")
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
	output$rnbOptsO.exploratory.columns <- renderText({
		cnames <- input$rnbOptsI.exploratory.columns
		if (is.null(cnames)) cnames <- NULL
		if (length(cnames)==1 && cnames=="[automatic]") cnames <- NULL
		cnames <- setdiff(cnames, "[automatic]")
		rnb.options(exploratory.columns=cnames)
		rnb.getOption("exploratory.columns")
	})
	output$rnbOptsO.exploratory.intersample <- renderText({
		rnb.options(exploratory.intersample=input$rnbOptsI.exploratory.intersample)
		rnb.getOption("exploratory.intersample")
	})
	output$rnbOptsO.differential.comparison.columns <- renderText({
		cnames <- input$rnbOptsI.differential.comparison.columns
		if (is.null(cnames)) cnames <- NULL
		if (length(cnames)==1 && cnames=="[automatic]") cnames <- NULL
		cnames <- setdiff(cnames, "[automatic]")
		rnb.options(differential.comparison.columns=cnames)
		rnb.getOption("differential.comparison.columns")
	})
	output$rnbOptsO.differential.report.sites <- renderText({
		rnb.options(differential.report.sites=input$rnbOptsI.differential.report.sites)
		rnb.getOption("differential.report.sites")
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
	output$runHint <- renderText({
		if (enableRun()){
			shinyjs::enable("numCores")
			shinyjs::enable("ggplotTheme")
			shinyjs::enable("runRnb")
			"Ready to run"
		} else {
			shinyjs::disable("numCores")
			shinyjs::disable("ggplotTheme")
			shinyjs::disable("runRnb")
			"Unable to run the analysis. Please make sure that you selected a non-existing report directory ('Analysis' tab) and that you specified the sample annotation file and data directory correctly ('Input' tab)"
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
	rnbSet <- NULL
	output$modImportAnaDir.about <- renderUI({
		curStatus <- anaStatus()
		if (curStatus$status=="reportDir"){
			tagList(tags$p("An existing RnBeads report has been specified. Ready to import."))
		} else {
			tagList(tags$p(icon("warning"), "Please specify an existing RnBeads report in the 'Analysis' tab to continue."))
		}
	})
	modImportStatus <- reactiveValues(
		dataset.loaded=FALSE,
		dataset.loaded.nsamples=-1
	)
	# watch if dataset has been loaded
	reactDataset <- observeEvent(modImportStatus$dataset.loaded, {
		if (modImportStatus$dataset.loaded){
			shinyjs::disable("modImportNew")
			shinyjs::disable("modImportAnaDir")
			shinyjs::disable("modImportRObjects")
			shinyjs::enable("modImportReset")
		} else {
			shinyjs::enable("modImportNew")
			if (anaStatus()$status=="reportDir") shinyjs::enable("modImportAnaDir")
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
		tagList(res)
	})
	observeEvent(input$modImportNew, {
		modImportStatus$dataset.loaded <- TRUE
		#TODO
	})
	observeEvent(input$modImportAnaDir, {
		modImportStatus$dataset.loaded <- TRUE
		#TODO
	})
	observeDirectoryInput(input, session, 'modImportRnBSetDir')
	observeLocalFileInput(input, session, 'modImportOptionsFile')
	observeEvent(input$modImportRObjects, {
		modImportStatus$dataset.loaded <- TRUE
		withProgress({
			tryCatch({
					rnbSet <- load.rnb.set(readDirectoryInput(session, 'modImportRnBSetDir'))
					modImportStatus$dataset.loaded.nsamples <- length(samples(rnbSet))
				},
				error = function(err) {
					rnbSet <- NULL
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
		rnbSet <- NULL
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
