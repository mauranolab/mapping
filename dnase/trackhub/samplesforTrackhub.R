#Add variables from command line tsv files of DSnumbers and Institute
library(data.table)
library(RColorBrewer)
library(dplyr)
library(optparse)

option_list = list(
	make_option(c("-f", "--file"), type="character", default=NULL, 
		help="input file name. Tab-delimited file containing GroupID column with ID of samples to be included. Age and Institution columns will be propagated to trachub", metavar="character"),
	make_option(c("-o", "--out"), type="character", default=NULL, 
		help="output file name", metavar="character"),
	make_option(c("-w", "--workingDir"), type="character", default=NULL, 
		help="full path working directory name", metavar="character"),
	make_option(c("-p", "--project"), type="character", default=NULL, 
		help="Enable custom behavior: [CEGS, humanENCODEdnase, mouseENCODEchipseq]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
	print_help(opt_parser)
	stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

#for debug
#opt=list(file="/vol/isg/encode/dnase/SampleIDs.tsv", out="/vol/isg/encode/dnase/SamplesForTrackhub.tsv")
#opt=list(file="/vol/isg/encode/mouseencode_chipseq_2018/SampleIDs.tsv",out="/vol/isg/encode/mouseencode_chipseq_2018/SamplesForTrackhub.tsv")

message('Input file: ', opt$file)

if (!is.null(opt$out)) {
	message('Output file: ', opt$out)
} else {
	message('Output file: ', 'stdout')
}

if (!is.null(opt$workingDir)) {
	setwd(opt$workingDir)
	pwd <- opt$workingDir
} else {
	pwd <- getwd()
}
message('Working Directory: ', pwd)

project=opt$project
message('Project: ', project)
##############################################################

#inputSampleIDs = read.table(opt$file,sep='\t',header=T,stringsAsFactors=F,na.strings=c(""," ","NA"), fill = TRUE, check.names=FALSE,colClasses=c('character',NA))
#inputSampleIDs = read(opt$file,fill=TRUE,header=T, check.names=FALSE,colClasses=c('character',NA))

inputSampleIDs <- as.data.frame(fread(opt$file))

if (is.null(inputSampleIDs$GroupID)) {
	stop("GroupID column is required", call.=FALSE)
} else {
	inputSampleIDs <- inputSampleIDs[order(inputSampleIDs$GroupID),]
	message('Number of unique groups:', length(unique(inputSampleIDs$GroupID)))
}

#TODO for now we just use the DS number and Age so we don't have to manually rename all Names
#Age <- fread('../SampleIDs_SampleAge.tsv')
#Age <- as.data.frame(Age)
#Age <- Age[order(Age$GroupID),]
#
#Age$Name <- inputSampleIDs$Name
#inputSampleIDs <- Age

# Put a space between the 2 values
message('Dimensions of Input file: ', dim(unique(inputSampleIDs))[1], ' ', dim(unique(inputSampleIDs))[2])

# Get names of sample directories.
#TODO replace with --recursive command line option
if(project=="CEGS") {
	mappeddirs <- NULL  # We'll store sample directory path names in this.
	
	# Loop over all the flowcells to find sample directories.
	for(flowcell in list.dirs(path=pwd, full.names=F, recursive = F)) {
        # TEMP
        if(flowcell=='FCH75CLAFXY') next;

		# Get full path names of all the sample directories within <flowcell>/Project_CEGS
		mappeddirs <- append(mappeddirs, list.dirs(path=paste(pwd, "/", flowcell, "/Project_CEGS", sep=""), full.names = TRUE, recursive = FALSE))
	}

	# Check for unwanted directories
	mappeddirs <- mappeddirs[grep('Project_CEGS/fastqc', mappeddirs, invert=TRUE)]
	mappeddirs <- mappeddirs[grep('Project_CEGS/new', mappeddirs, invert=TRUE)]
} else {
	mappeddirs <- list.dirs(path=pwd, full.names=F, recursive = F)
	
	message('Number of unique groups that have a mapped directory: ', length(mappeddirs[grep(paste(unique(inputSampleIDs$GroupID),collapse="|"), mappeddirs)]))
	
	mappeddirs <- mappeddirs[grep(paste(unique(inputSampleIDs$GroupID), collapse="|"), mappeddirs)]
}

# Get the names into the form (for CEGS):  <flowcell>/Project_CEGS/<sample_directory>
# For non-CEGS, they are already like this:  <flowcell>
# So nothing happens in the non-CEGS case.
mappeddirs <- gsub(paste0(pwd, "/"), "", mappeddirs, fixed=T)

#get rid of anything in a bak or trash directory
mappeddirs <- mappeddirs[grep('bak', mappeddirs, invert=T)]
mappeddirs <- mappeddirs[grep('trash', mappeddirs, invert=T)]


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
colnames(data) <- c('Name', 'DS', 'Replicate', 'Color', 'Assay', 'analyzed_reads', 'SPOT', 'Num_hotspots', 'Exclude', 'Genomic_coverage', 'Group', 'Age', 'Uni', "filebase")

i <- 0            # This will be our "data" output variable index.
for(curdir in mappeddirs){
    #message(curdir)
	SampleID <- basename(curdir)
	SampleIDsplit <- unlist(strsplit(SampleID, "-"))
	
	#NB makeTracks is obsolete naming convention
	if(project=="CEGS") {
		mergedFiles <- list.files(path=paste0(pwd, '/', curdir), pattern="^(makeTracks|analysis).*")
	} else {
		mergedFiles <- list.files(path=paste0(pwd, '/', SampleID), pattern="^(makeTracks|analysis).*")
	}
	
	#message(mergedFiles)
	for(mergedFile in mergedFiles) {
		sampleFile <- readLines(paste0(pwd, '/', curdir, '/', mergedFile), n=2000)
		
		
		#TODO Use unlist(strsplit(sampleFile[grep('^Running ', sampleFile)], ","))[2] instead?
		if(project=="CEGS") {
			# Check for DNase analysis files.
			if(any(grepl('dnase analysis', sampleFile))) next;
		}
		if(project=="CEGS") {
			curAssay <- 'DNA';
		} else if(length(SampleIDsplit) == 2) {
			curAssay <- 'DNase'
		} else if(length(SampleIDsplit) == 3) {
			curAssay <- SampleIDsplit[2]
		} else {
			message("WARNING: can't parse SampleID properly")
			curAssay <- 'UNKNOWN'
		}
		
		
		# If Done! is not there, we move on to the next mergedFile in the for loop.
		if(tail(sampleFile, 2)[1] == 'Done!'){
			i <- i+1
			# We need to add a new row to "data".  The values will be set within this for loop.
			data[i,] <- NA
			
            # Need to test that this is working as intended in all cases.
			#colGroup <- col2rgb(as.character(col[rownames(col) %in% gsub("-.*", '', mappeddirs[j]),][1]))[,1]
			colGroup <- col2rgb(as.character(col[rownames(col) %in% gsub("-.*", '', curdir),][1]))[,1]
			
			data$Name[i] <- SampleIDsplit[1]
			data$DS[i] <- SampleIDsplit[length(SampleIDsplit)]
			data$Replicate[i] <- 1 #BUGBUG???
			data$Color[i] <- paste(colGroup[1], colGroup[2], colGroup[3], sep=',')
			data$Assay[i] <- curAssay
			data$analyzed_reads[i] <- strsplit(sampleFile[grep('Num_analyzed_reads\t', sampleFile)], '\t')[[1]][2]
			
			if(any(grepl('Num_hotspots2', sampleFile))) {
				data$Num_hotspots[i] <- strsplit(sampleFile[grep('Num_hotspots2\t', sampleFile)], '\t')[[1]][2]
			} else {
				data$Num_hotspots[i] <- NA
			}

			if(any(grepl('SPOT2', sampleFile))) {
				data$SPOT[i] <- strsplit(sampleFile[grep('SPOT2\t', sampleFile)], '\t')[[1]][2]
			} else {
				data$SPOT[i] <- NA
			}

			if(any(grepl('Genomic_coverage', sampleFile))) {
				data$Genomic_coverage[i] <- strsplit(sampleFile[grep('Genomic_coverage\t', sampleFile)], '\t')[[1]][2]
			} else {
				data$Genomic_coverage[i] <- NA
			}
			
			data$Exclude[i] <- NA
			
			#BUGBUG looking up by GroupID gives multiple results that triggers R warning when multiple fastq files exist per group (true for CEGS ?)
			if(project=="CEGS") {
				# Could put Age in the input file too, if needed.
				data$Age[i] <- "NoAge"
				data$Group[i] <- unlist(strsplit(curdir, "/"))[1] #use FC name
				data$Uni[i]  <- inputSampleIDs[inputSampleIDs[,"GroupID"]==data$Group[i], "Institution"]
			} else {
				data$Age[i] <- inputSampleIDs[grep(data$DS[i], inputSampleIDs$GroupID), "Age"]
				data$Uni[i] <- inputSampleIDs[grep(data$DS[i], inputSampleIDs$GroupID), "Institution"]
				
				#pre-populate group field with institution name
				if(project=="humanENCODEdnase") {
					if(data$Uni[i] == 'UMass') {data$Uni[i] <- 'UW'}
				}
				data$Group[i] <- data$Uni[i]
			}
			
			if(project=="CEGS") {
				data$filebase[i] <- paste0(curdir, "/", paste0(unlist(strsplit(basename(mergedFile), "\\."))[2:3], collapse="."))
			} else {
                # Needs to be tested
				data$filebase[i] <- paste0(SampleID, "/", paste0(unlist(strsplit(basename(mergedFile), "\\."))[2:3], collapse="."))
			}
			
			if(project=="humanENCODEdnase") {
				if(data$Group[i] != "Duke") {
					#Divide samples based on category
					if(length(grep('^CD|^[hi]?T[hHR][0-9]+$|^GM[012][0-9][0-9][0-9][0-9]$1|B_cell|neutrophil|natural_killer|regulatory_T_cell|^MEL$|macrophage|CH12LX|G1E|mononuclear|dendritic', data$Name[i]))>0){data$Group[i]<-'Hematopoietic cells'}
					if(length(grep('ES|^H[0-9]|^iPS|E14TG2a4',data$Name[i]))>0){data$Group[i]<-'Pluripotent'}
					if(length(grep('^f[A-Z]', data$Name[i]))>0){data$Group[i] <- 'Fetal-REMC'}
					if(length(grep('testis|spinal|Skin|bladder|urothelia|ventriculus|colon|limb|placenta|heart|cortex|kidney|bone|pancrea|cardia|eye|renal|gonad|muscle|osteo|medulla|brain|ovary|olfact|uteru|fibroblast|lung|tongue|bowel|putamen|esopha|gastro|ammon|derm|nucleus|gast|glom|gyrus|thyroid|adipo|neuron|prostate|intest|medull|[Ll]iver|aggregated_lymphoid_nodule|aorta|artery|psoas|stomach|testes|tibial_artery|vagina|omental_fat_depot|fetal_umbilical_cord|pons|medial_popliteal_nerve|globus|Spleen', data$Name[i][grep('^f[A-Z]', data$Name[i], invert=T)], ignore.case=T, invert=F))>0){data$Group[i] <- 'Tissues'}
				}
			} else if(project=="mouseENCODEchipseq"){
				if(length(grep('Broad|Duke|HAIB|HMS|Stanford|UMass|USC|UTA|UW|Yale|UCSD', data$Group[i], ignore.case=T))>0){data$Group[i]<-'Cell lines'}
				if(length(grep('^CD|^[hi]?T[hHR][0-9]+$|^GM[012][0-9][0-9][0-9][0-9]$1|B_cell|neutrophil|natural_killer|regulatory_T_cell|^MEL$|macrophage|CH12LX|G1E|mononuclear|dendritic', data$Name[i]))>0){data$Group[i]<-'Hematopoietic cells'}
				if(length(grep('ES|^H[0-9]|^iPS|E14TG2a4' ,data$Name[i]))>0){data$Group[i]<-'Pluripotent'}
				if(length(grep('HUES|H54|C2C12|myocyte', data$Name[i]))>0){data$Group[i]<-'Cell lines'}
				if(length(grep('mesenchymal_stem_cell|trophoblast_cell|testis|spinal|aorta|body|breast|coronary_artery|neural_cell|omental_fat_pad|right_atrium_auricular_region|right_lobe_of_liver|spleen|stomach|tibial_artery|tibial_nerve|vagina|Skin|bladder|urothelia|ventriculus|colon|limb|placenta|heart|cortex|kidney|bone|oesteoblast|pancrea|cardia|eye|renal|gonad|muscle|osteo|medulla|brain|ovary|olfact|uteru|fibroblast|lung|tongue|bowel|putamen|liver|esopha|gastro|ammon|derm|nucleus|gast|glom|gyrus|thyroid|adipo|neuron|dendritic_cell|hepatocyte|mid_neurogenesis_radial_glial_cells|neuroepithelial_stem_cell|radial_glial_cell|prostate|intest|medull|thymus|cerebellum|cortical_plate',
				data$Name[i],ignore.case=T,invert=F))>0 || length(grep('day',data$Age[i]))>0) {
					if(length(grep('^H3',data$Assay[i]))>0) {
						data$Group[i]<-'Tissues-Histone marks'
					} else {
						data$Group[i]<-'Tissues-TFs'
					}
				}
				if(length(grep('eGFP|FLAG|MCF10A_Er_Src', data$Name[i], ignore.case=T))>0){data$Group[i]<-'Epitope-tagged TFs'}
				if(length(grep('control', SampleID, ignore.case=T))>0){data$Group[i]<-'Control'}
				if(data$Assay[i] == "CTCF"){data$Group[i]<-'CTCF'}
			}
		}
	}
}


if(project %in% c("humanENCODEdnase", "mouseENCODEchipseq")) {
	#Fix sample age. 
	#Only keep first entry e.g. male (week 7) male (week8)
	data$Age <- gsub(').*', '', data$Age)
	data$Age[grep('day', data$Age)] <- paste0(round(as.numeric(gsub(' day| days', '', data$Age[grep('day', data$Age)]))/7), ' weeks')
	data$Age[grep('^8 ', data$Age)] <- '08 weeks' 

	##Rename the fetal tissues 
	#fetalRename <- data[grep('day|week', data$Age),] %>% 
	#	filter(!grepl('^f', Name)) %>%
	#	filter(!grepl('AG04449|AG04450|IMR_90', Name)) %>%
	#	select(DS)
	#
	#library(Hmisc)
	#data[data$DS%in%fetalRename$DS,]$Name <- gsub('^', 'f', capitalize(data[data$DS%in%fetalRename$DS,]$Name))
	#data[data$DS%in%fetalRename$DS,]$Group <- 'Fetal_roadmap'

	#Add replicate numbers based on highest number of non redundant reads
	data <- data[!is.na(data$Name),]
	Replicates <- unique(subset(data, select=c(Name, Group)))
	Replicates$Name <- gsub('_L$|_R$', '', Replicates$Name)
	for (i in 1:nrow(Replicates)){
		#message(Replicates$Name[i], Replicates$Group[i])
		data[gsub('_L$|_R$', '', data$Name)==Replicates$Name[i] & data$Group==Replicates$Group[i],]$Replicate <- rank(-as.numeric(data[gsub('_L$|_R$', '', data$Name)==Replicates$Name[i] & data$Group==Replicates$Group[i],]$analyzed_reads))
		#Fix color for Left and right tissues
		data[gsub('_L$|_R$', '', data$Name)==Replicates$Name[i] & data$Group==Replicates$Group[i],]$Color <- names(tail(table(data[gsub('_L$|_R$', '', data$Name)==Replicates$Name[i],]$Color), 1))
	}
	
	
	dataRep <- data
	dataRep$Name <- gsub('_L$|_R$', '', dataRep$Name)
	dataRep <- split(dataRep, dataRep$Name)
	#repDF <- dataRep[[3]]
	dataReplist <- lapply(dataRep, function(repDF) {
		head(repDF)
		#repDF$Replicate <- NA
		if ( nrow(repDF[grep('UW', repDF$Uni),]) >0 ) {
			 repDF[grep('UW', repDF$Uni),]$Replicate <- rank(-as.numeric(repDF[grep('UW', repDF$Uni),]$analyzed_reads))
			 repDF[grep('UW', repDF$Uni, invert=T),]$Replicate <- rank(-as.numeric(repDF[grep('UW', repDF$Uni , invert=T),]$analyzed_reads)) + nrow(repDF[grep('UW', repDF$Uni),])
		} else {
			repDF$Replicate <- rank(-as.numeric(repDF$analyzed_reads))
		}
		repDF
	})
	
	data <- do.call(rbind, dataReplist)
	
	data[data$Replicate>2,]$Replicate <- 'Other'
	#data$Group[data$Group=='UMass'] <- 'UW'
	data[is.na(data$Age),]$Age <- 'NoAge'
	data$Replicate <- paste0('rep', data$Replicate)
} else {
	data$Replicate <- 'REP?'
}


#Output file:
write.table(subset(data, select=-c(Uni)), file=opt$out, sep='\t', col.names=T, row.names=F, quote=F)

# message(warnings())

message("Done!!!")

