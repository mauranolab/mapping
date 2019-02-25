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
		help="Enable custom directory search and group behavior: [CEGS, humanENCODEdnase, mouseENCODEdnase, humanENCODEchipseq, mouseENCODEchipseq]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if(!is.null(opt$out)) {
	message('Output file: ', opt$out)
} else {
	message('Output file: ', 'stdout')
}

if(!is.null(opt$workingDir)) {
	setwd(opt$workingDir)
	pwd <- opt$workingDir
} else {
	pwd <- getwd()
}
message('Working Directory: ', pwd)

project=opt$project
message('Project: ', project)


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
		message('Number of unique DS numbers: ', length(unique(inputSampleIDs$DS)))
	}
	message('Dimensions of Input file: ', nrow(inputSampleIDs), ' rows x ', ncol(inputSampleIDs), ' cols')
} else {
	# This is the case for project==CEGS
	inputSampleIDs <- NULL
}


# Get paths of directories to search for sample directories relative to pwd
if(project=="CEGS") {
	projectdirs <- NULL
	
	# Loop over all the flowcells to find sample directories.
	flowcell_dates = list()
	for(flowcell in list.dirs(path=pwd, full.names=F, recursive = F)) {
		# To temporarily skip flowcell directories which are incomplete.
		# if(flowcell=='FCH2NNMBBXY') next;
		
		projectdirs <- append(projectdirs, paste0(flowcell, "/", list.dirs(path=paste0(pwd, "/", flowcell, "/"), full.names=F, recursive = F)))
		
		# Get flowcell date
		infofile <- paste('/vol/mauranolab/flowcells/data/', flowcell, '/info.txt', sep="")
		if(file.exists(infofile)) {
			infofile_lines <- readLines(infofile)
			flowcell_dates[[flowcell]] <- gsub(unlist(strsplit(infofile_lines[grep("#Load date", infofile_lines)], '\t'))[2], pattern='-', replacement='')
		} else { 
			message("WARNING No info.txt file found for ", infofile)
			flowcell_dates[[flowcell]] <- NA
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
mappeddirs <- mappeddirs[grep('Project_CEGS/new', mappeddirs, invert=TRUE)]
mappeddirs <- mappeddirs[grep('bak', mappeddirs, invert=TRUE)]
mappeddirs <- mappeddirs[grep('trash', mappeddirs, invert=TRUE)]

#NB These show up in old pipeline
mappeddirs <- mappeddirs[grep('hotspots', mappeddirs, invert=TRUE)]
mappeddirs <- mappeddirs[grep('fastqc', mappeddirs, invert=TRUE)]

message('Mapped directories to process: ', length(mappeddirs))


#Initialize color palette
set.seed(12345)
colorPalette <- c(brewer.pal(n=9, 'Set1'), brewer.pal(n=8, 'Dark2'), brewer.pal(n=12, 'Paired')[c(FALSE, TRUE)])
colorPalette <- gsub('#FFFF33', '#F0027F', colorPalette)
colorPalette <- gsub('#FFFFFF', '#BEAED4', colorPalette)
colorPalette <- sapply(colorPalette, FUN=function(x) { paste(col2rgb(x), collapse=',') })
#Store as 0-indexed so that %% works, and add 1 to access color
nextColorFromPalette <- 0
colorAssignments <- NULL


# Initialize "data" with just column names.  We'll be adding rows to this later on in the code.
data <- data.frame(matrix(ncol=14, nrow=1))
colnames(data) <- c("Name", "DS", "Replicate", "Color", "Assay", "analyzed_reads", "Genomic_coverage", "SPOT", "Num_hotspots", "Exclude", "Group", "Age", "Institution", "filebase")
i <- 0 # This will be our "data" output variable index.
for(curdir in mappeddirs){
	message("Working on ", curdir)
	SampleID <- basename(curdir)
	SampleIDsplit <- unlist(strsplit(SampleID, "-"))
	
	#NB makeTracks is obsolete naming convention
	analysisFiles <- list.files(path=paste0(pwd, '/', curdir), pattern="^(makeTracks|analysis).*")
	
	if(length(analysisFiles)==0) {
		message("WARNING No analysis files found in ", curdir)
		next # Nothing follows here except for the long 'for(analysisFile...' loop.
	}
	
	for(analysisFile in analysisFiles) {
		analysisFileContents <- readLines(paste0(pwd, '/', curdir, '/', analysisFile), n=2000)
		
		if(tail(analysisFileContents, 2)[1] != "Done!"){
			message("WARNING ", analysisFile, " appears not to have completed successfully")
		} else {
			i <- i+1
			# We need to add a new row to "data".  The values will be set within this for loop.
			data[i,] <- NA
			
			pipelineParameters <- analysisFileContents[grep('^Running [^,]+,[^,]+ analysis', analysisFileContents, perl=T)]
			if(length(pipelineParameters)>0) {
				pipelineParameters <- unlist(strsplit(pipelineParameters, " "))[2]
				analysisCommand <- unlist(strsplit(pipelineParameters, ","))[2]
			} else {
				#assume dnase for old pipeline
				analysisCommand <- "dnase"
			}
			
			if(analysisCommand=="dnase") {
				data$Assay[i] <- "DNase-seq"
			} else if(analysisCommand=="callsnps") {
				data$Assay[i] <- "DNA"
			} else if(analysisCommand=="chipseq") {
				data$Assay[i] <- SampleIDsplit[2]
			} else {
				message("WARNING don't recognize analysisCommand: ", analysisCommand)
				data$Assay[i] <- NA
			}
			
			data$Name[i] <- SampleIDsplit[1]
			
			data$DS[i] <- SampleIDsplit[length(SampleIDsplit)]
			
			data$Replicate[i] <- NA
			
			#TODO parameterize setting a different key for color lookup
			if(is.null(colorAssignments) || !data$Name[i] %in% colorAssignments$group) {
				colorAssignments <- rbind(colorAssignments, data.frame(group=data$Name[i], rgb=colorPalette[nextColorFromPalette+1], stringsAsFactors=F))
				nextColorFromPalette <- (nextColorFromPalette + 1) %% length(colorPalette)
			}
			data$Color[i] <- colorAssignments[colorAssignments$group == data$Name[i], "rgb"]
			
			data$analyzed_reads[i] <- strsplit(analysisFileContents[grep('^Num_analyzed_(tags|reads)\t', analysisFileContents)], '\t')[[1]][2] #Tags is for old pipeline
			
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
				inputSampleIDrow <- which(inputSampleIDs[,"DS"] == data$DS[i])[1]
			}
			
			if("Age" %in% colnames(inputSampleIDs)) {
				data$Age[i] <- inputSampleIDs[inputSampleIDrow, "Age"]
			} else {
				data$Age[i] <- NA
			}
			
			if("Institution" %in% colnames(inputSampleIDs)) {
				data$Institution[i] <- inputSampleIDs[inputSampleIDrow, "Institution"]
				if(project=="humanENCODEdnase") {
					if(data$Institution[i] == 'UMass') {data$Institution[i] <- 'UW'}
				}
			} else {
				data$Institution[i] <- NA
			}
			
			if(project=="CEGS") {
				data$Institution[i] <- "NYU" # Since we will stop using input file for project==CEGS
				
				curFC <- unlist(strsplit(curdir, "/"))[1] 
				data$Group[i] <- curFC
				if(curFC %in% names(flowcell_dates)) {
					data$Group[i] <- paste0(flowcell_dates[[curFC]] , "_" , data$Group[i]) #Group values will now be in the form of YYYMMDD_<flowcell>
				}
			} else if(project %in% c("humanENCODEdnase", "mouseENCODEdnase", "humanENCODEchipseq", "mouseENCODEchipseq")) {
				if(grepl('^f[A-Z]|^mfLiver_', data$Name[i])) {
					data$Group[i] <- "Fetal tissues"
				} else if(grepl('^(adipo|aggregated_lymphoid_nodule|adrenal|ammon|aorta|artery|astrocyte|bladder|body|bone|bowel|brain|breast|bronchial_epithelial_cell|cardia|cardiocyte|cerebellum|colon|coronary_artery|cortex|cortical_plate|dendritic_cell|derm|endothelial_cell_of_|epithelial_cell_of_choroid_plexus|erythroblast|esopha|eye|fetal_umbilical_cord|fibroblast|gast|gastro|globus|glom|gonad|gyrus|heart|hepatocyte|intest|keratinocyte|kidney|limb|liver|lung|mammary_epithelial_cell|medial_popliteal_nerve|medull|medulla|mesenchymal_stem_cell|mid_neurogenesis_radial_glial_cells|muscle|myotube|neural_cell|neural_progenitor_cell|neural_stem_progenitor_cell|neuroepithelial_stem_cell|neuron|nucleus|oesteoblast|olfact|fat_pad|osteo|ovary|pancrea|placenta|pons|prostate|psoas|putamen|radial_glial_cell|renal|retina|retinal_pigment_epithelial_cell|right_atrium_auricular_region|right_lobe_of_liver|skin|spinal|spleen|stomach|test[ei]s|thymus|thyroid|tibial_artery|tibial_nerve|tongue|trophoblast_cell|urothelia|uteru|vagina|ventriculus|amniotic_stem_cell|bipolar_spindle_neuron|caput_mediale_musculus_gastrocnemius|inferior_parietal_cortex|islet_precursor_cell|midbrain|middle_frontal_gyrus|pentadactyl_limb|ascending_aorta|bipolar_neuron|epithelial_cell_of_esophagus|epithelial_cell_of_prostate|foreskin_keratinocyte|lower_leg_skin|Peyers_patch|right_cardiac_atrium|sigmoid_colon|skeletal_muscle|small_intestine|smooth_muscle_cell|suprapubic_skin|thoracic_aorta|transverse_colon|upper_lobe_of_left_lung|urinary_bladder|brown_adipose_tissue|forebrain|hindbrain|myocyte|Muller_cell|telencephalon)', data$Name[i], ignore.case=T)) {
					data$Group[i] <- "Tissues"
				}
				if(grepl('^CD|^[him]?A?T[HhNnRr][0-9]*$|^GM[012][0-9][0-9][0-9][0-9]|m?B_?cell|neutrophil|natural_killer|regulatory_T_cell|^MEL$|^MEL_GATA1_ER$|macrophage|CH12LX|G1E|mononuclear|dendritic|leukemia_stem_cell', data$Name[i])) { data$Group[i] <- 'Hematopoietic cells' }
				if(grepl('ES|^H[0-9]|^iPS|E14TG2a4|^trophoblastic_cell$|^mesendoderm$|^endodermal_cell$|^ectodermal_cell$|^mesodermal_cell$|^WW6$|^ZHBTc4$', data$Name[i])) { data$Group[i] <- 'Pluripotent' }
				if(project=="humanENCODEdnase") {
					if(is.na(data$Group[i]) || data$Institution[i] == "Duke") {
						data$Group[i] <- data$Institution[i]
					}
				} else if(project %in% c("mouseENCODEdnase", "mouseENCODEchipseq", "humanENCODEchipseq")) {
					if(is.na(data$Group[i])) { data$Group[i] <- "Cell lines" }
					if(data$Group[i]=="GM12878") { data$Group[i] <- "Cell lines" }
					if(project %in% c("mouseENCODEchipseq", "humanENCODEchipseq")) {
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


if(project %in% c("humanENCODEdnase", "mouseENCODEdnase", "humanENCODEchipseq", "mouseENCODEchipseq")) {
	#Fix sample age. 
	#Only keep first entry e.g. male (week 7) male (week8)
	data$Age <- gsub(').*', '', data$Age)
	data$Age[grep('^[0-9]+ days?$', data$Age)] <- paste0(round(as.numeric(gsub(' days?', '', data$Age[grep('^[0-9]+ days?$', data$Age)]))/7), ' weeks')
	
	#Seems to be obselete
	#Add replicate numbers based on highest number of non redundant reads
#	Replicates <- unique(subset(data, select=c(Name, Assay)))
#	Replicates$Name <- gsub('_L$|_R$', '', Replicates$Name)
#	for (i in 1:nrow(Replicates)){
#		#message(Replicates$Name[i], Replicates$Group[i])
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


#Output file:
write.table(data, file=opt$out, sep='\t', col.names=T, row.names=F, quote=F)

# message(warnings())

message("Done!!!")
