#!/bin/env Rscript

#Add variables from command line tsv files of DSnumbers and Institute
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

option_list = list(
	make_option(c("--inputfile"), type="character", default=NULL,
		help="input file name. Optional tab-delimited file. Entries will be matched on DS column, and Age and Institution columns will be propagated to trackhub", metavar="character"),
	make_option(c("--out"), type="character", default=NULL,
		help="output file name", metavar="character"),
	make_option(c("--workingDir"), type="character", default=NULL, 
		help="full path working directory name", metavar="character"),
	make_option(c("--project"), type="character", default="",
		help="Enable custom directory search and group behavior: [byFC, CEGS_byLocus, humanENCODEdnase, mouseENCODEdnase, humanENCODEchipseq, mouseENCODEchipseq]", metavar="character"),
	make_option(c("--descend"), action="store_true", type="logical",
		help="Assume that sample folders are organized under two levels of subdirectories, e.g. FCxxx/dna")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(opt$project != "" && !opt$project %in% c("byFC", "CEGS_byLocus", "humanENCODEdnase", "mouseENCODEdnase", "humanENCODEchipseq", "mouseENCODEchipseq")) {
    message("[samplesforTrackhub] invalid option for --project ", opt$project)
    quit(save="no", status=1)
}

if(opt$project %in% c("byFC", "CEGS_byLocus") & ! "descend" %in% names(opt)) {
	stop("ERROR --project byFC or --project CEGS_byLocus require --descend for now")
}

if(!is.null(opt$out)) {
	message("[samplesforTrackhub] ", 'Output file: ', opt$out)
} else {
	message("[samplesforTrackhub] ", 'Output file: ', 'stdout')
}

if(!is.null(opt$workingDir)) {
	setwd(opt$workingDir)
	pwd <- opt$workingDir
} else {
	pwd <- getwd()
}

message("[samplesforTrackhub] ", 'Working Directory: ', pwd, "; Project: ", opt$project)

##############################################################
#for debug
#opt=list(file="/vol/isg/encode/dnase/SampleIDs.tsv", out="/vol/isg/encode/dnase/SamplesForTrackhub.tsv")
#opt=list(file="/vol/isg/encode/mouseencode_chipseq_2018/SampleIDs.tsv",out="/vol/isg/encode/mouseencode_chipseq_2018/SamplesForTrackhub.tsv")

if(!is.null(opt$inputfile)) {
	inputSampleIDs <- as.data.frame(fread(opt$inputfile))
	if(is.null(inputSampleIDs$DS)) {
		stop("DS column is required in provided inputfile", call.=FALSE)
	} else {
		inputSampleIDs <- inputSampleIDs[order(inputSampleIDs$DS),]
		message("[samplesforTrackhub] ", 'Number of unique DS numbers: ', length(unique(inputSampleIDs$DS)))
	}
	message("[samplesforTrackhub] ", 'Dimensions of Input file: ', nrow(inputSampleIDs), ' rows x ', ncol(inputSampleIDs), ' cols')
} else {
	inputSampleIDs <- NULL
}


# Get paths of directories to search for sample directories relative to pwd
if("descend" %in% names(opt)) {
	projectdirs <- NULL
	
	# Loop over all the flowcells to find sample directories.
	flowcell_dates = list()
	for(flowcell in list.dirs(path=pwd, full.names=F, recursive = F)) {
		if(flowcell=='src') next;
		if(flowcell=='trackhub') next;
		
		# To temporarily skip flowcell directories which are incomplete.
		# if(flowcell=='FCH2NNMBBXY') next;
		
		projectdirs <- append(projectdirs, paste0(flowcell, "/", list.dirs(path=paste0(pwd, "/", flowcell, "/"), full.names=F, recursive = F)))
		
		if(opt$project=="byFC") {
			# Get flowcell date
			infofile <- paste('/vol/mauranolab/flowcells/data/', flowcell, '/info.txt', sep="")
			if(file.exists(infofile)) {
				infofile_lines <- readLines(infofile)
				flowcell_dates[[flowcell]] <- gsub(unlist(strsplit(infofile_lines[grep("#Load date", infofile_lines)], '\t'))[2], pattern='-', replacement='')
			} else {
				message("[samplesforTrackhub] ", "WARNING No info.txt file found for ", infofile)
				flowcell_dates[[flowcell]] <- NA
			}
		}
	}
} else {
	projectdirs <- "."
}


# Get paths of sample directories relative to pwd
mappeddirs <- NULL
for(curdir in projectdirs) {
	thisProjectMappedDirs <- list.dirs(path=paste0(pwd, "/", curdir), full.names=F, recursive = F)
	if(curdir!=".") {
		thisProjectMappedDirs <- paste0(curdir, "/", thisProjectMappedDirs)
	}
	mappeddirs <- append(mappeddirs, thisProjectMappedDirs)
}

# Prune unwanted directories
mappeddirs <- mappeddirs[grep('/new', mappeddirs, invert=TRUE)]
mappeddirs <- mappeddirs[grep('/bak', mappeddirs, invert=TRUE)]
mappeddirs <- mappeddirs[grep('/trash', mappeddirs, invert=TRUE)]

#NB These show up in old pipeline
mappeddirs <- mappeddirs[grep('hotspots', mappeddirs, invert=TRUE)]
mappeddirs <- mappeddirs[grep('fastqc', mappeddirs, invert=TRUE)]

message("[samplesforTrackhub] ", 'Mapped directories to process: ', length(mappeddirs))


#Initialize color palette
set.seed(12345)
colorPalette <- c(brewer.pal(n=9, 'Set1'), brewer.pal(n=8, 'Dark2'), brewer.pal(n=12, 'Paired')[c(FALSE, TRUE)])
colorPalette <- gsub('#FFFF33', '#F0027F', colorPalette)
colorPalette <- gsub('#FFFFFF', '#BEAED4', colorPalette)
colorPalette <- sapply(colorPalette, FUN=function(x) { paste(col2rgb(x), collapse=',') })
#Store as 0-indexed so that %% works, and add 1 to access color
nextColorFromPalette <- 0
colorAssignments <- NULL


# The "pipelineParameters" parser function. Will be used in the below for loop.
pipelineParametersParser <- function(pipelineParameters, fieldName) {
	regex_out <- regexpr('^Running (?<processingCommand>[^,]+),(?<analysisCommand>[^ ]+) analysis [^\\-]+\\-[^ ]+ \\(.+\\) against genome (?<mappedgenome>[^ ]+)( \\(aka (?<annotationgenome>[^ ]+)\\))?', pipelineParameters, perl=T)

	regex_out.start <- attr(regex_out,"capture.start")[1, fieldName]
	regex_out.lngth <- attr(regex_out,"capture.length")[1, fieldName]
	return(substr(pipelineParameters, regex_out.start, regex_out.start+regex_out.lngth-1))
}

# Initialize "data" with just column names.  We'll be adding rows to this later on in the code.
outputCols <- c("Name", "SampleID", "Assay", "Group", "filebase", "Mapped_Genome", "Annotation_Genome", "Color", "analyzed_reads", "Genomic_coverage", "SPOT", "Num_hotspots", "Exclude", "Age", "Institution", "Replicate", "Bait_set")
if(opt$project == "CEGS_byLocus") {
    outputCols <- c(outputCols, "Study", "Project", "Assembly", "Type")
}

data <- data.frame(matrix(ncol=length(outputCols), nrow=1))
colnames(data) <- outputCols
i <- 0 # This will be our "data" output variable index.
for(curdir in mappeddirs) {
	message("[samplesforTrackhub] ", curdir)
	SampleID <- basename(curdir)
	SampleIDsplit <- unlist(strsplit(SampleID, "-"))
	
	#NB makeTracks is obsolete naming convention
	analysisFiles <- list.files(path=paste0(pwd, '/', curdir), pattern="^(makeTracks|analysis).*")
	
	if(length(analysisFiles)==0) {
		message("[samplesforTrackhub] ", "WARNING No analysis files found in ", curdir)
		next # Nothing follows here except for the long 'for(analysisFile...' loop.
	}
	
	for(analysisFile in analysisFiles) {
		analysisFileContents <- readLines(paste0(pwd, '/', curdir, '/', analysisFile), n=2000)
		
		if(tail(analysisFileContents, 2)[1] != "Done!"){
			message("[samplesforTrackhub] ", "WARNING ", analysisFile, " appears not to have completed successfully")
		} else {
			i <- i+1
			# We need to add a new row to "data". The values will be set within this for loop.
			data[i,] <- NA
			
			pipelineParameters <- analysisFileContents[grep('^Running [^,]+,[^,]+ analysis', analysisFileContents, perl=T)]
			if(length(pipelineParameters)>0) {
				sampleType <- pipelineParametersParser(pipelineParameters, "analysisCommand")
				mappedgenome <- pipelineParametersParser(pipelineParameters, "mappedgenome")
				annotationgenome <- pipelineParametersParser(pipelineParameters, "annotationgenome")
				
				data$Mapped_Genome[i] <- mappedgenome
				
				# Is the annotation geneome in the "Running ... analysis ..." line ?
 				if(nchar(annotationgenome) > 0) {
					# It is.
 					data$Annotation_Genome[i] <- annotationgenome
				} else {
					# It is not. Extract it from the "mappedgenome" field.
 					annotationgenome <- gsub("_.+$", "", mappedgenome)
 					data$Annotation_Genome[i] <- gsub("all$", "", annotationgenome)
				}
			} else {
				#assume dnase for old pipeline
				sampleType <- "dnase"
				data$Mapped_Genome[i] <- NA
				data$Annotation_Genome[i] <- NA
			}
			
			# Adding a new Assay type also requires changes to be made to MakeTrackhub.py
			# Look for the initialization of "assay_type" in MakeTrackhub.py for comments on this.
			if(sampleType=="dnase") {
				data$Assay[i] <- "DNase-seq"
			} else if(sampleType=="atac") {
				data$Assay[i] <- "ATAC-seq"
			} else if(sampleType=="none") {
				data$Assay[i] <- "None"
			} else if(sampleType=="dna" || sampleType=="callsnps") {
				data$Assay[i] <- "DNA"
			} else if(sampleType=="capture" || sampleType=="callsnpsCapture") {
				data$Assay[i] <- "DNA Capture"
			} else if(sampleType=="chipseq") {
				data$Assay[i] <- SampleIDsplit[2]
			} else {
				message("[samplesforTrackhub] ", "WARNING don't recognize sampleType: ", sampleType)
				data$Assay[i] <- NA
			}
			
			data$Name[i] <- SampleIDsplit[1]
			
			data[i, "SampleID"] <- SampleIDsplit[length(SampleIDsplit)]
			
			data$Replicate[i] <- NA
			
			#TODO parameterize setting a different key for color lookup
			if(is.null(colorAssignments) || !data$Name[i] %in% colorAssignments$group) {
				colorAssignments <- rbind(colorAssignments, data.frame(group=data$Name[i], rgb=colorPalette[nextColorFromPalette+1], stringsAsFactors=F))
				nextColorFromPalette <- (nextColorFromPalette + 1) %% length(colorPalette)
			}
			data$Color[i] <- colorAssignments[colorAssignments$group == data$Name[i], "rgb"]
			
			data$analyzed_reads[i] <- strsplit(analysisFileContents[grep('^Num_analyzed_(tags|reads)\t', analysisFileContents)], '\t')[[1]][2] #Tags is for old pipeline
			
			#Parse SampleAnnotation as a list of values with the key as name
			if(any(grep('^SampleAnnotation\t', analysisFileContents))) {
				SampleAnnotation <- strsplit(analysisFileContents[grep('^SampleAnnotation\t', analysisFileContents)], '\t')[[1]][2]
				#Split up comma-delimited key/value pairs
				SampleAnnotation <- strsplit(SampleAnnotation, ",")[[1]]
				#Split up key/value
				SampleAnnotation <- sapply(SampleAnnotation, USE.NAMES=F, FUN=function(x) {
					ret=list()
					for(cur in strsplit(x, "=")) { 
						ret[cur[1]] = cur[2] 
					}
					return(ret)
				})
			} else {
				SampleAnnotation <- list()
			}
			
			if( !is.null(SampleAnnotation[["Bait_set"]]) ) {
				data$Bait_set[i] <- SampleAnnotation[["Bait_set"]]
			} else {
				data$Bait_set[i] <- NA
			}
			
			if(any(grepl('^Genomic_coverage\t', analysisFileContents))) {
				data$Genomic_coverage[i] <- strsplit(analysisFileContents[grep('^Genomic_coverage\t', analysisFileContents)], '\t')[[1]][2]
			} else {
				data$Genomic_coverage[i] <- NA
			}
			
			if(any(grepl('^Num_hotspots2\t', analysisFileContents))) {
				data$Num_hotspots[i] <- strsplit(analysisFileContents[grep('^Num_hotspots2\t', analysisFileContents)], '\t')[[1]][2]
			} else {
				data$Num_hotspots[i] <- NA
			}
			
			if(any(grepl('^SPOT\t', analysisFileContents))) {
				data$SPOT[i] <- strsplit(analysisFileContents[grep('^SPOT\t', analysisFileContents)], '\t')[[1]][2]
			} else {
				data$SPOT[i] <- NA
			}
			
			data$Exclude[i] <- NA
			
			if(!is.null(inputSampleIDs)) {
				#Matching the DS number from the analysis file to inputSampleIDs can give multiple results, so just pick the first
				inputSampleIDrow <- which(inputSampleIDs[,"DS"] == data[i, "SampleID"])[1]
			}
			
			#TODO allow other columns to be taken from inputSampleIDs
			if("Age" %in% colnames(inputSampleIDs)) {
				data$Age[i] <- inputSampleIDs[inputSampleIDrow, "Age"]
			} else {
				data$Age[i] <- NA
			}
			
			if("Institution" %in% colnames(inputSampleIDs)) {
				data$Institution[i] <- inputSampleIDs[inputSampleIDrow, "Institution"]
				if(opt$project=="humanENCODEdnase") {
					if(data$Institution[i] == 'UMass') {data$Institution[i] <- 'UW'}
				}
			} else {
				data$Institution[i] <- NA
			}
			
			if(opt$project %in% c("byFC", "CEGS_byLocus")) {
				data$Institution[i] <- "NYU"
				
				if(opt$project=="byFC") {
					curFC <- unlist(strsplit(curdir, "/"))[1]
					data$Group[i] <- curFC
					if(curFC %in% names(flowcell_dates)) {
						#Group values will be in the form of YYYMMDD_<flowcell>
						data$Group[i] <- paste0(flowcell_dates[[curFC]] , "_" , data$Group[i])
					}
				} else if(opt$project=="CEGS_byLocus") {
					#Group values will be in the form of Study ID
					SampleNameSplit <- unlist(strsplit(SampleIDsplit[1], "_"))
					CEGSsampleType <- SampleNameSplit[length(SampleNameSplit)]
					if(CEGSsampleType %in% c("Yeast", "DNA", "BAC", "RepoBAC", "Ecoli", "Amplicon")) {
						data$Study[i] <- SampleNameSplit[1]
						data$Project[i] <- SampleNameSplit[2]
						data$Assembly[i] <- SampleNameSplit[3]
						data$Type[i] <- SampleNameSplit[4]
						data$Group[i] <- data$Study[i]
					}
					#TODO Find capture samples based on bait?
				} else {
					stop("ERROR Impossible!")
				}
			} else if(opt$project %in% c("humanENCODEdnase", "mouseENCODEdnase", "humanENCODEchipseq", "mouseENCODEchipseq")) {
				if(grepl('^f[A-Z]|^mfLiver_', data$Name[i])) {
					data$Group[i] <- "Fetal tissues"
				} else if(grepl('^(adipo|aggregated_lymphoid_nodule|adrenal|ammon|aorta|artery|astrocyte|bladder|body|bone|bowel|brain|breast|bronchial_epithelial_cell|cardia|cardiocyte|cerebellum|colon|coronary_artery|cortex|cortical_plate|dendritic_cell|derm|endothelial_cell_of_|epithelial_cell_of_choroid_plexus|erythroblast|esopha|eye|fetal_umbilical_cord|fibroblast|gast|gastro|globus|glom|gonad|gyrus|heart|hepatocyte|intest|keratinocyte|kidney|limb|liver|lung|mammary_epithelial_cell|medial_popliteal_nerve|medull|medulla|mesenchymal_stem_cell|mid_neurogenesis_radial_glial_cells|muscle|myotube|neural_cell|neural_progenitor_cell|neural_stem_progenitor_cell|neuroepithelial_stem_cell|neuron|nucleus|oesteoblast|olfact|fat_pad|osteo|ovary|pancrea|placenta|pons|prostate|psoas|putamen|radial_glial_cell|renal|retina|retinal_pigment_epithelial_cell|right_atrium_auricular_region|right_lobe_of_liver|skin|spinal|spleen|stomach|test[ei]s|thymus|thyroid|tibial_artery|tibial_nerve|tongue|trophoblast_cell|urothelia|uteru|vagina|ventriculus|amniotic_stem_cell|bipolar_spindle_neuron|caput_mediale_musculus_gastrocnemius|inferior_parietal_cortex|islet_precursor_cell|midbrain|middle_frontal_gyrus|pentadactyl_limb|ascending_aorta|bipolar_neuron|epithelial_cell_of_esophagus|epithelial_cell_of_prostate|foreskin_keratinocyte|lower_leg_skin|Peyers_patch|right_cardiac_atrium|sigmoid_colon|skeletal_muscle|small_intestine|smooth_muscle_cell|suprapubic_skin|thoracic_aorta|transverse_colon|upper_lobe_of_left_lung|urinary_bladder|brown_adipose_tissue|forebrain|hindbrain|myocyte|Muller_cell|telencephalon)', data$Name[i], ignore.case=T)) {
					data$Group[i] <- "Tissues"
				}
				if(grepl('^CD|^[him]?A?T[HhNnRr][0-9]*$|^GM[012][0-9][0-9][0-9][0-9]|m?B_?cell|neutrophil|natural_killer|regulatory_T_cell|^MEL$|^MEL_GATA1_ER$|macrophage|CH12LX|G1E|mononuclear|dendritic|leukemia_stem_cell', data$Name[i])) { data$Group[i] <- 'Hematopoietic cells' }
				if(grepl('ES|^H[0-9]|^iPS|E14TG2a4|^trophoblastic_cell$|^mesendoderm$|^endodermal_cell$|^ectodermal_cell$|^mesodermal_cell$|^WW6$|^ZHBTc4$', data$Name[i])) { data$Group[i] <- 'Pluripotent' }
				if(opt$project=="humanENCODEdnase") {
					if(is.na(data$Group[i]) || data$Institution[i] == "Duke") {
						data$Group[i] <- data$Institution[i]
					}
				} else if(opt$project %in% c("mouseENCODEdnase", "mouseENCODEchipseq", "humanENCODEchipseq")) {
					if(is.na(data$Group[i])) { data$Group[i] <- "Cell lines" }
					if(data$Group[i]=="GM12878") { data$Group[i] <- "Cell lines" }
					if(opt$project %in% c("mouseENCODEchipseq", "humanENCODEchipseq")) {
						if(grepl("[Tt]issues$", data$Group[i])) {
							if(grepl('^H[234][ABFK]', data$Assay[i])) {
								data$Group[i] <- paste0(data$Group[i], "-Histone marks")
							} else {
								data$Group[i] <- paste0(data$Group[i], "-TFs")
							}
						}
						if(grepl('eGFP|(3x)?FLAG', data$Assay[i], ignore.case=T)) { data$Group[i] <- 'Epitope-tagged TFs' }
						if(grepl('control', data$Assay[i], ignore.case=T)) { data$Group[i] <- 'Control' }
						if(data$Assay[i] == "CTCF") { data$Group[i] <- 'CTCF' }
					}
				}
			} else {
				data$Group[i] <- data$Institution[i]
			}
			
			data$filebase[i] <- paste0(curdir, "/", paste0(unlist(strsplit(basename(analysisFile), "\\."))[2:3], collapse="."))
		}
	}
}


if(opt$project %in% c("humanENCODEdnase", "mouseENCODEdnase", "humanENCODEchipseq", "mouseENCODEchipseq")) {
	#Fix sample age. 
	#Only keep first entry e.g. male (week 7) male (week8)
	data$Age <- gsub(').*', '', data$Age)
	data$Age[grep('^[0-9]+ days?$', data$Age)] <- paste0(round(as.numeric(gsub(' days?', '', data$Age[grep('^[0-9]+ days?$', data$Age)]))/7), ' weeks')
	
	#Seems to be obselete
	#Add replicate numbers based on highest number of non redundant reads
#	Replicates <- unique(subset(data, select=c(Name, Assay)))
#	Replicates$Name <- gsub('_L$|_R$', '', Replicates$Name)
#	for (i in 1:nrow(Replicates)){
#		#message("[samplesforTrackhub] ", Replicates$Name[i], Replicates$Group[i])
#		matchingRows <- gsub('_L$|_R$', '', data$Name)==Replicates$Name[i] & data$Assay==Replicates$Assay[i]
#		data[matchingRows,]$Replicate <- rank(-as.numeric(data[matchingRows,]$analyzed_reads))
#		#Fix color for Left and right tissues
#		data[matchingRows,]$Color <- names(tail(table(data[gsub('_L$|_R$', '', data$Name)==Replicates$Name[i],]$Color), 1))
#	}
	
	
	dataRep <- data
	dataRep$splitBy <- paste(gsub('_L$|_R$', '', dataRep$Name), dataRep$Assay, sep="-")
	dataRep <- split(dataRep, dataRep$splitBy)
	dataRepList <- lapply(dataRep, function(repDF) {
		if( nrow(repDF[grep('UW', repDF$Institution),]) >0 ) {
			repDF[grep('UW', repDF$Institution),]$Replicate <- rank(-as.numeric(repDF[grep('UW', repDF$Institution),]$analyzed_reads))
			repDF[grep('UW', repDF$Institution, invert=T),]$Replicate <- rank(-as.numeric(repDF[grep('UW', repDF$Institution , invert=T),]$analyzed_reads)) + nrow(repDF[grep('UW', repDF$Institution),])
		} else {
			repDF$Replicate <- rank(-as.numeric(repDF$analyzed_reads))
		}
		return(repDF)
	})
	data <- do.call(rbind, dataRepList)
	data <- data[,!colnames(data)=="splitBy"]
	
	data[data$Replicate>=3,]$Replicate <- 'Other'
	#data$Group[data$Group=='UMass'] <- 'UW'
	data$Replicate <- paste0('rep', data$Replicate)
}

if(opt$project=="CEGS_byLocus") {
	# Delete line items not associated with one of the enumerated "CEGSsampleType" sample types.
	data <- data[!is.na(data$Group),]
}

#Output file:
write.table(data, file=opt$out, sep='\t', col.names=T, row.names=F, quote=F)

# message("[samplesforTrackhub] ", warnings())

message("[samplesforTrackhub] ", "Done!!!")

