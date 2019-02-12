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
		help="Enable custom directory search and group behavior: [CEGS, humanENCODEdnase, mouseENCODEchipseq]", metavar="character")
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
	message('Dimensions of Input file: ', nrow(unique(inputSampleIDs)), ' rows x ', nrow(unique(inputSampleIDs)))
} else {
	inputSampleIDs <- NULL
}


# Get paths of directories to search for sample directories relative to pwd
if(project=="CEGS") {
	projectdirs <- NULL
	
	# Loop over all the flowcells to find sample directories.
	flowcell_dates = list()
	for(flowcell in list.dirs(path=pwd, full.names=F, recursive = F)) {
		# To temporarily skip flowcell directories which are incomplete.
		# if(flowcell=='FCH75CLAFXY') next;
		
		projectdirs <- append(projectdirs, paste0(flowcell, "/", list.dirs(path=paste0(pwd, "/", flowcell, "/"), full.names=F, recursive = F)))
		
		# Get flowcell date
		infofile <- paste('/vol/mauranolab/flowcells/data/', flowcell, '/info.txt', sep="")
		if(exists(infofile)) {
			infofile_lines <- readLines(infofile)
			flowcell_dates[[flowcell]] <- gsub(unlist(strsplit(infofile_lines[grep("#Load date", infofile_lines)], '\t'))[2], pattern='-', replacement='')
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

#Thes should show up in old pipeline
mappeddirs <- mappeddirs[grep('hotspots', mappeddirs, invert=TRUE)]
mappeddirs <- mappeddirs[grep('fastqc', mappeddirs, invert=TRUE)]


message('Mapped directories to process: ', length(mappeddirs))


#pick colors
set.seed(12345)
col <- c(brewer.pal(n=9, 'Set1'), brewer.pal(n=8, 'Dark2'), brewer.pal(n=12, 'Paired')[c(FALSE, TRUE)])
col <- gsub('#FFFF33', '#F0027F', col)
col <- gsub('#FFFFFF', '#BEAED4', col)

#Color by groups
#TODO permit coloring by assay?
groupnames <- gsub("-.*", '', mappeddirs)
col <- rep(col, round(length(groupnames)/length(col), 2))[as.factor(unique(groupnames))]	 
col <- as.data.frame(col, unique(groupnames))


# Initialize "data" with just column names.  We'll be adding rows to this later on in the code.
data <- data.frame(matrix(ncol=14, nrow=1))
colnames(data) <- c('Name', 'DS', 'Replicate', 'Color', 'Assay', 'analyzed_reads', 'Genomic_coverage', 'SPOT', 'Num_hotspots', 'Exclude', 'Group', 'Age', 'Institution', "filebase")
i <- 0 # This will be our "data" output variable index.
for(curdir in mappeddirs){
	message("Working on ", curdir)
	SampleID <- basename(curdir)
	SampleIDsplit <- unlist(strsplit(SampleID, "-"))
	
	#NB makeTracks is obsolete naming convention
	analysisFiles <- list.files(path=paste0(pwd, '/', curdir), pattern="^(makeTracks|analysis).*")
	
	if(length(analysisFiles)==0) {
		message("WARNING No analysis files found in ", curdir)
	}
	
	for(analysisFile in analysisFiles) {
		analysisFileContents <- readLines(paste0(pwd, '/', curdir, '/', analysisFile), n=2000)
		
		# If Done! is not there, we move on to the next analysisFile in the for loop.
		if(tail(analysisFileContents, 2)[1] == 'Done!'){
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
				message("WARNING: can't parse SampleID properly")
				data$Assay[i] <- NA
			}
			
			colGroup <- col2rgb(as.character(col[rownames(col) %in% gsub("-.*", '', curdir),][1]))[,1]
			
			data$Name[i] <- SampleIDsplit[1]
			data$DS[i] <- SampleIDsplit[length(SampleIDsplit)]
			if(project=="CEGS") {
				data$Name[i] <- paste(data$Name[i], data$DS[i], sep="-")
			}
			data$Replicate[i] <- NA
			data$Color[i] <- paste(colGroup[1], colGroup[2], colGroup[3], sep=',')
			data$analyzed_reads[i] <- strsplit(analysisFileContents[grep('^Num_analyzed_(tags|reads)\t', analysisFileContents)], '\t')[[1]][2] #Tags is for old pipeline
			
			if(any(grepl('^Genomic_coverage', analysisFileContents))) {
				data$Genomic_coverage[i] <- strsplit(analysisFileContents[grep('^Genomic_coverage\t', analysisFileContents)], '\t')[[1]][2]
			} else {
				data$Genomic_coverage[i] <- NA
			}
			
			if(any(grepl('^Num_hotspots2', analysisFileContents))) {
				data$Num_hotspots[i] <- strsplit(analysisFileContents[grep('^Num_hotspots2\t', analysisFileContents)], '\t')[[1]][2]
			} else {
				data$Num_hotspots[i] <- NA
			}
			
			if(any(grepl('^SPOT2', analysisFileContents))) {
				data$SPOT[i] <- strsplit(analysisFileContents[grep('^SPOT2\t', analysisFileContents)], '\t')[[1]][2]
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
				curFC <- unlist(strsplit(curdir, "/"))[1] 
				data$Group[i] <- curFC
				if(curFC %in% names(flowcell_dates)) {
					data$Group[i] <- paste0(flowcell_dates[[curFC]] , "_" , data$Group[i]) #Group values will now be in the form of YYYMMDD_<flowcell>
				}
			} else if(project=="humanENCODEdnase") {
				if(data$Institution[i] != "Duke") {
					#Divide samples based on category
					if(length(grep('^CD|^[hi]?T[hHR][0-9]+$|^GM[012][0-9][0-9][0-9][0-9]|B_cell|neutrophil|natural_killer|regulatory_T_cell|macrophage|CH12LX|G1E|mononuclear|dendritic', data$Name[i]))>0){data$Group[i]<-'Hematopoietic cells'}
					if(length(grep('ES|^H[0-9]|^iPS|E14TG2a4',data$Name[i]))>0){data$Group[i]<-'Pluripotent'}
					if(length(grep('^f[A-Z]', data$Name[i]))>0){data$Group[i] <- 'Fetal-REMC'}
					if(length(grep('testis|spinal|Skin|bladder|urothelia|ventriculus|colon|limb|placenta|heart|cortex|kidney|bone|pancrea|cardia|eye|renal|gonad|muscle|osteo|medulla|brain|ovary|olfact|uteru|fibroblast|lung|tongue|bowel|putamen|esopha|gastro|ammon|derm|nucleus|gast|glom|gyrus|thyroid|adipo|neuron|prostate|intest|medull|[Ll]iver|aggregated_lymphoid_nodule|aorta|artery|psoas|stomach|testes|tibial_artery|vagina|omental_fat_depot|fetal_umbilical_cord|pons|medial_popliteal_nerve|globus|Spleen', data$Name[i][grep('^f[A-Z]', data$Name[i], invert=T)], ignore.case=T, invert=F))>0){data$Group[i] <- 'Tissues'}
				}
			} else if(project=="mouseENCODEchipseq"){
				if(length(grep('Broad|Duke|HAIB|HMS|Stanford|UMass|USC|UTA|UW|Yale|UCSD', data$Group[i], ignore.case=T))>0){data$Group[i]<-'Cell lines'}
				if(length(grep('^CD|^[hi]?T[hHR][0-9]+$|^GM[012][0-9][0-9][0-9][0-9]|B_cell|neutrophil|natural_killer|regulatory_T_cell|^MEL$|macrophage|CH12LX|G1E|mononuclear|dendritic', data$Name[i]))>0){data$Group[i]<-'Hematopoietic cells'}
				if(length(grep('ES|^H[0-9]|^iPS|E14TG2a4' ,data$Name[i]))>0){data$Group[i]<-'Pluripotent'}
				if(length(grep('HUES|H54|C2C12|myocyte', data$Name[i]))>0){data$Group[i]<-'Cell lines'}
				if(length(grep('mesenchymal_stem_cell|trophoblast_cell|testis|spinal|aorta|body|breast|coronary_artery|neural_cell|omental_fat_pad|right_atrium_auricular_region|right_lobe_of_liver|spleen|stomach|tibial_artery|tibial_nerve|vagina|Skin|bladder|urothelia|ventriculus|colon|limb|placenta|heart|cortex|kidney|bone|oesteoblast|pancrea|cardia|eye|renal|gonad|muscle|osteo|medulla|brain|ovary|olfact|uteru|fibroblast|lung|tongue|bowel|putamen|liver|esopha|gastro|ammon|derm|nucleus|gast|glom|gyrus|thyroid|adipo|neuron|dendritic_cell|hepatocyte|mid_neurogenesis_radial_glial_cells|neuroepithelial_stem_cell|radial_glial_cell|prostate|intest|medull|thymus|cerebellum|cortical_plate', data$Name[i],ignore.case=T,invert=F))>0 || length(grep('day',data$Age[i]))>0) {
					if(length(grep('^H3',data$Assay[i]))>0) {
						data$Group[i]<-'Tissues-Histone marks'
					} else {
						data$Group[i]<-'Tissues-TFs'
					}
				}
				if(length(grep('eGFP|FLAG|MCF10A_Er_Src', data$Name[i], ignore.case=T))>0){data$Group[i]<-'Epitope-tagged TFs'}
				if(length(grep('control', SampleID, ignore.case=T))>0){data$Group[i]<-'Control'}
				if(data$Assay[i] == "CTCF"){data$Group[i]<-'CTCF'}
			} else {
				#pre-populate group field with institution name
				data$Group[i] <- data$Institution[i]
			}
			
			data$filebase[i] <- paste0(curdir, "/", paste0(unlist(strsplit(basename(analysisFile), "\\."))[2:3], collapse="."))
		}
	}
}


if(project %in% c("humanENCODEdnase", "mouseENCODEchipseq")) {
	#Fix sample age. 
	#Only keep first entry e.g. male (week 7) male (week8)
	data$Age <- gsub(').*', '', data$Age)
	data$Age[grep('day', data$Age)] <- paste0(round(as.numeric(gsub(' day| days', '', data$Age[grep('day', data$Age)]))/7), ' weeks')
	
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
	dataReplist <- lapply(dataRep, function(repDF) {
		if( nrow(repDF[grep('UW', repDF$Institution),]) >0 ) {
			repDF[grep('UW', repDF$Institution),]$Replicate <- rank(-as.numeric(repDF[grep('UW', repDF$Institution),]$analyzed_reads))
			repDF[grep('UW', repDF$Institution, invert=T),]$Replicate <- rank(-as.numeric(repDF[grep('UW', repDF$Institution , invert=T),]$analyzed_reads)) + nrow(repDF[grep('UW', repDF$Institution),])
		} else {
			repDF$Replicate <- rank(-as.numeric(repDF$analyzed_reads))
		}
		return(repDF)
	})
	data <- do.call(rbind, dataReplist)
	
	data[data$Replicate>=3,]$Replicate <- 'Other'
	#data$Group[data$Group=='UMass'] <- 'UW'
	data$Replicate <- paste0('rep', data$Replicate)
}


#Output file:
write.table(data, file=opt$out, sep='\t', col.names=T, row.names=F, quote=F)

# message(warnings())

message("Done!!!")

