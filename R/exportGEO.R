########################################################################################################################
## exportGEO.R
## created: 2017-09-25
## creator: Pavlo Lutsik
## ---------------------------------------------------------------------------------------------------------------------
## Uploading 450k/EPIC data to GEO.
########################################################################################################################



# prepareGEOSampleInfoTemplate
# 
# Prepare a GEO sample information template for further use in \code{prepareSOFTfileForGEO}
# @param rnb.set                   Object inheriting from class \code{RnBeadSet}
#                                  with \code{"GSE"}.
# @param sample.source.col         \code{integer} specifying a column in the \code{pheno} slot of \code{rnb.set} 
#                                  containing information which will be written into the field \code{Sample_source_name_ch1}
#                                  of each sample record
# @param sample.description        \code{integer} specifying a column in the \code{pheno} slot of \code{rnb.set} 
#                                  containing information which will be written into the field \code{Sample_desription}
#                                  of each sample record
prepareGEOSampleInfoTemplate<-function(rnb.set, pheno.cols, sample.source.col, sample.description.col)
{
	pheno.data<-pheno(rnb.set)[pheno.cols]
	sample.names<-samples(rnb.set)
	sample.sources<-pheno(rnb.set)[[sample.source.col]]
	sample.descr<-pheno(rnb.set)[[sample.description.col]]
	
	platform <- c("HumanMethylation27", "HumanMethylation450", "HumanMethylationEPIC")[match(rnb.set@target, c('probes27','probes450','probesEPIC'))]
	chipVersion<-"unknown"
	organism <- c("H.sapiens", "M.muscullus", "R.norvegicus")[which(!is.na(pmatch( c('hg','mm','rn'), assembly(rnb.set))))]
	platformID<-c(
			"HumanMethylation27"="GPL8490",
			"HumanMethylation450"="GPL13534",
			"HumanMethylationEPIC"="GPL23976")[platform]
	
	templateTitle <- c(
			"Sample_title", 
			"Sample_channel_count",
			"Sample_source_name_ch1", 
			"Sample_organism_ch1", 
		rep("Sample_characteristics_ch1", ncol(pheno.data)),
			"Sample_molecule_ch1", 
			"Sample_extract_protocol_ch1",
			"Sample_label_ch1", 
			"Sample_label_protocol_ch1", 
			"Sample_hyb_protocol",
			"Sample_scan_protocol", 
			"Sample_description", 
			"Sample_data_processing",
			"Sample_platform_id", 
			"Sample_supplementary_file")
	
	templateContent <- c(
			"", 
			"2", 
			"", 
			organism, 
			colnames(pheno.data),
			"genomic DNA",
			"standard as recommended by Illumina", 
			"Cy3,Cy5", 
			"standard as recommended by Illumina",
			"standard as recommended by Illumina", 
			"standard as recommended by Illumina",
			"", 
			"", 
			platformID, 
			"none")
	
	norm_method<-rnb.set@status$normalized
	bgcorr_method<-rnb.set@status$background
	
	preprocessMethod <- paste("normalization: ", norm_method, "; ", "background: ", bgcorr_method, sep="")
	
	templateContent[templateTitle == "Sample_data_processing"] <- preprocessMethod
	template <- templateTitle
	pheno.cols<-grep("Sample_characteristics_ch1", template)
	for (i in seq(sample.names)) {
		templateContentCopy<-templateContent
		
		templateContentCopy[templateTitle == "Sample_title"] <- sample.names[i]
		templateContentCopy[templateTitle == "Sample_source_name_ch1"] <- sample.sources[i]
		templateContentCopy[templateTitle == "Sample_description"] <- sample.descr[i]
		
		templateContentCopy[pheno.cols] <- paste(templateContentCopy[pheno.cols], as.character(pheno.data[i,]), sep=": ")
		template <- rbind(template, templateContentCopy)
	}
	template <- cbind(c("sampleID", sample.names), template)
	
	colnames(template)<-template[1,]
	template<-template[-1,]
	
	return(template)
}

#'
#' prepareSOFTfileForGEO
#'
#' Starting from an \code{RnBeadSet} object generates a batch submission file for Gene Expression Omnibus series in SOFT format 
#'
#' @param rnb.set                   Object inheriting from class \code{RnBeadSet}
#'                                  with \code{"GSE"}.
#' @param filename                  Absolute path or a name of a SOFT file to be generated
#' @param sample.source.col         \code{integer} singleton specifying a column in the \code{pheno} slot of \code{rnb.set} 
#'                                  containing information which will be written into the field \code{Sample_source_name_ch1}
#'                                  of each sample record
#' @param sample.description        \code{integer} singleton specifying a column in the \code{pheno} slot of \code{rnb.set} 
#'                                  containing information which will be written into the field \code{Sample_desription}
#'                                  of each sample record
#' @param export.cols               \code{integer} vector specifying columns in the \code{pheno} slot of \code{rnb.set} 
#'                                  containing information which will be written into the fields \code{Sample_characteristics_ch1}
#'                                  of each sample record
#' @param series.info               A \code{list} with elements to be written to the series record. Elements should be character
#'                                  singletons named \code{SERIES} (contains a valid GSE identifier for updating an existing series)
#'                                  \code{Series_title}, \code{Series_summary}, \code{Series_type}, \code{Series_overall_design},
#'                                  \code{Series_contributor}, \code{Series_sample_id}
#' 
#' @return \code{TRUE} on success.
#'
#' @details The code was largely adapted from a similar function in package \code{lumi} which is due to Pan Du.
#' 
#' @author Pavlo Lutsik
#' @export
#'
prepareSOFTfileForGEO<-function(rnb.set, 
		filename, 
		sample.source.col, 
		sample.description.col,
		export.cols=seq(ncol(pheno(rnb.set))),
		rnb.set.raw=NULL, 
		series.info=NULL){
	
	if(!inherits(rnb.set, "RnBeadSet")){
		rnb.error("Wrong value for rnb.set supplied: not an RnBeadSet object")
	}
	
	if(!is.null(rnb.set.raw) && !inherits(rnb.set.raw, "RnBeadSet")){
		rnb.error("Wrong value for rnb.set.raw supplied: not an RnBeadSet object")
	}
	
	if(!(is.character(filename) && length(filename)==1)){
		rnb.error("Wrong value for rnb.set.raw supplied: not an RnBeadSet object")
	}
	
	expr.norm <- signif(meth(rnb.set),5)
	
	rnb.set.intens<-NULL
	if (is.null(rnb.set.raw)) {
		meth.raw <- NULL
		rnb.set.intens<-rnb.set
	}else{
		meth.raw <- signif(meth(rnb.set.raw),5)
		rnb.set.intens<-rnb.set.raw
	}
	
	if(!is.null(dpval(rnb.set.intens))){
		detect <- signif(dpval(rnb.set.intens), 5)
	}else{
		detect <- NULL
	}
	
	if(inherits(rnb.set.intens, "RnBeadRawSet") && !is.null(M(rnb.set.intens)) && !is.null(U(rnb.set.intens))){
		methyData <- signif(M(rnb.set.intens), 5)
		unmethyData <- signif(U(rnb.set.intens), 5)
	}else{
		methyData<-NULL
		unmethyData<-NULL
	}
	
	sampleInfo <- prepareGEOSampleInfoTemplate(rnb.set, 1:ncol(pheno(rnb.set)), 2, 2)
	
	sampleInfoTitle <- colnames(sampleInfo)
	if (any(sapply(sampleInfo[, -1, drop = F], nchar) == 0))
		stop("No blank fields are allowed in the sampleInfo table!\nYou can check some example submissions, like GSM296418, at the GEO website.\n")
	
	nuID <- rownames(annotation(rnb.set))
	probeId <- nuID
	
#	if (length(which(is.nuID(sample(nuID, 100)))) < 20) {
#		nuID <- NULL
#	}
	
	sampleID <- sampleInfo[, "sampleID"]
	sampleTitle <- sampleInfo[, "Sample_title"]
	outputFile <- filename
	for (i in seq(sampleID)) {
		cat("Processing sample", i, "\n")
		if (i == 1) {
			cat("^SAMPLE =", sampleTitle[i], "\n", sep = "",
					file = outputFile, append = FALSE)
		}
		else {
			cat("^SAMPLE =", sampleTitle[i], "\n", sep = "",
					file = outputFile, append = TRUE)
		}
		sampleInfo.i <- paste("!", sampleInfoTitle[-1], " = ",
				sampleInfo[i, -1], "\n", sep = "", collapse = "")
		sampleInfo.i <- gsub("'", "\\'", sampleInfo.i)
		cat(sampleInfo.i, file = outputFile, sep = "", append = TRUE)
		tableHead <- "ID_REF"
		cat("#ID_REF = Illumina ID\n", file = outputFile, append = TRUE)
#		if (!is.null(nuID)) {
#			cat("#nuID = nucleotide universal IDentifier (nuID), convertible to and from probe sequence. See Bioconductor lumi package for more details.\n",
#					file = outputFile, append = TRUE)
#			tableHead <- c(tableHead, "nuID")
#		}
		cat("#VALUE = Beta-value\n", file = outputFile, append = TRUE)
		if (!is.null(meth.raw))
			cat("#RAW_VALUE = raw Beta-value\n", file = outputFile,
					append = TRUE)
		tableHead <- c(tableHead, "VALUE")
		if (!is.null(meth.raw))
			tableHead <- c(tableHead, "RAW_VALUE")
		if (!is.null(methyData)) {
			cat("#METHYLATED = the intensities measured by methylated probes\n",
					file = outputFile, append = TRUE)
			tableHead <- c(tableHead, "METHYLATED")
		}
		if (!is.null(unmethyData)) {
			cat("#UNMETHYLATED = the intensities measured by unmethylated probes\n",
					file = outputFile, append = TRUE)
			tableHead <- c(tableHead, "UNMETHYLATED")
		}
		if (!is.null(detect)) {
			cat("#Detection_Pval = the detection p-value of the probe\n",
					file = outputFile, append = TRUE)
			tableHead <- c(tableHead, "Detection_Pval")
		}
		sampleTable.i <- probeId
#		if (!is.null(nuID))
#			sampleTable.i <- cbind(sampleTable.i, nuID)
		sampleTable.i <- cbind(sampleTable.i, expr.norm[, sampleID[i],
						drop = FALSE])
		if (!is.null(meth.raw))
			sampleTable.i <- cbind(sampleTable.i, meth.raw[, sampleID[i],
							drop = FALSE])
		if (!is.null(methyData))
			sampleTable.i <- cbind(sampleTable.i, methyData[,
							sampleID[i], drop = FALSE])
		if (!is.null(unmethyData))
			sampleTable.i <- cbind(sampleTable.i, unmethyData[,
							sampleID[i], drop = FALSE])
		if (!is.null(detect))
			sampleTable.i <- cbind(sampleTable.i, detect[, sampleID[i],
							drop = FALSE])
		sampleTable.i <- rbind(tableHead, sampleTable.i)
		cat("!sample_table_begin\n", file = outputFile, append = TRUE)
		write.table(sampleTable.i, sep = "\t", quote = FALSE,
				file = outputFile, append = TRUE, col.names = FALSE,
				row.names = FALSE)
		cat("!sample_table_end\n", file = outputFile, append = TRUE)
	}
	
	if(!is.null(series.info)){
		label_series<-c(
				'Series_title',
				'Series_summary',
				'Series_type',
				'Series_overall_design')
		
		cat(paste('^SERIES', series.info[["SERIES"]], sep=' = '), file = outputFile, append = TRUE)
		
		labelvalue_series<-paste('!',sapply(label_series,function(x){paste(x,series.info[[x]], sep=' = ')}), sep='')
		cat(labelvalue_series, file = outputFile, append = TRUE)
		
		##write.soft.append(paste.soft('!Series_summary',seriestable$Series_summary))
		
		seriescontributor<-strsplit(series.info$Series_contributor,';')[[1]]
		cat(paste('!Series_contributor',seriescontributor, sep=' = '), file = outputFile, append = TRUE)
		
		cat(paste('!Series_sample_id', paste(sampleID, sep=";") ,sep=' = '), file = outputFile, append = TRUE)
	
	}
	return TRUE
}

