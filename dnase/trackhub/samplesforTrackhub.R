#Add variables from command line tsv files of DSnumbers and Institute
library("optparse")
library(data.table)
library(RColorBrewer)
library(dplyr)
option_list = list(
	make_option(c("-f", "--file"), type="character", default=NULL, 
		help="dataset file name", metavar="character"),
	make_option(c("-o", "--out"), type="character", default=NULL, 
		help="output file name", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$file)){
	print_help(opt_parser)
	stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#for debug
#opt=list(file="/vol/isg/encode/dnase201805/SampleIDs_20180502_MTM.tsv",out="/vol/isg/encode/dnase201805/SamplesForTrackhub.tsv")
#opt=list(file="/vol/isg/encode/mouseencode_chipseq_2018/SampleIDs.tsv",out="/vol/isg/encode/mouseencode_chipseq_2018/SamplesForTrackhub.tsv")

cat('Input file: ', opt$file,'\n')
cat('Output file: ', opt$out,'\n')

#df = read.table(opt$file,sep='\t',header=T,stringsAsFactors=F,na.strings=c(""," ","NA"), fill = TRUE, check.names=FALSE,colClasses=c('character',NA))
#df = read(opt$file,fill=TRUE,header=T, check.names=FALSE,colClasses=c('character',NA))
df <- fread(opt$file)
df <- as.data.frame(df)
df <- df[order(df$GroupID),]
#TODO for now we just use the DS number and Age so we don't have to manually rename all celltypes
#Age <- fread('../SampleIDs_SampleAge.tsv')
#Age <- as.data.frame(Age)
#Age <- Age[order(Age$GroupID),]
#
#Age$cellType <- df$cellType
#df <- Age
pwd <- getwd()
cat('Dimensions of Input file: ', dim(unique(df)), '\n')

if (is.null(df$GroupID)){
	stop("GroupID column is required", call.=FALSE)
} else {
	cat('Number of unique groups:', length(unique(df$GroupID)),'\n')
}


mappeddirs <- list.dirs(path='.', full.names=F, recursive = F)

cat('Number of unique groups that have a mapped directory: ', length(mappeddirs[grep(paste(unique(df$GroupID),collapse="|"), mappeddirs)]), '\n')
mappeddirs <- mappeddirs[grep(paste(unique(df$GroupID), collapse="|"), mappeddirs)]
#get rid of anything in a bak directory
mappeddirs <- mappeddirs[grep('bak', mappeddirs, invert=T)]

#pick colors
set.seed(12345)
col <- c(brewer.pal(n=9, 'Set1'), brewer.pal(n=8, 'Dark2'), brewer.pal(n=12, 'Paired')[c(FALSE, TRUE)])
col <- gsub('#FFFF33', '#F0027F', col)
col <- gsub('#FFFFFF', '#BEAED4', col)

groupnames <- gsub("-.*", '', mappeddirs)
col <- rep(col, round(length(groupnames)/length(col), 2))[as.factor(unique(groupnames))]	 
col <- as.data.frame(col, unique(groupnames))

data <- data.frame(matrix(ncol=12, nrow=length(mappeddirs)))
colnames(data) <- c('cellType', 'DSnumber', 'Replicate', 'Color', 'Assay', 'analyzed_reads', 'SPOT', 'Hotspots', 'Exclude', 'Variable', 'Age', 'Uni')

#Change cellType names for Fetal tissues


for(i in 1:length(mappeddirs)){
	SampleID <- mappeddirs[i]
	cat(SampleID, '\n')
	mergedFiles <- list.files(path=paste0(pwd, '/', SampleID), pattern="^(makeTracks|analysis).*")
	#cat(mergedFiles, '\n')
	if(length(mergedFiles) > 0)
		sampleFile <- readLines(paste0(pwd, '/', SampleID, '/', mergedFiles), n=2000)
		if(tail(sampleFile, 2)[1] == 'Done!'){
			SampleIDsplit <- unlist(strsplit(SampleID, "-"))
			
			colGroup <- col2rgb(as.character(col[rownames(col) %in% gsub("-.*", '', mappeddirs[i]),][1]))[,1]
			data$cellType[i] <- paste(SampleIDsplit[-length(SampleIDsplit)], collapse="-")
			data$DSnumber[i] <- SampleIDsplit[length(SampleIDsplit)]
			data$Replicate[i] <- 1 #BUGBUG???
			data$Color[i] <- paste(colGroup[1], colGroup[2], colGroup[3], sep=',')
			
			#TODO hardcoded. Use unlist(strsplit(sampleFile[grep('^Running ', sampleFile)], ","))[2]
#			data$Assay[i] <- 'DNase'
			data$Assay[i] <- SampleIDsplit[2]
			
			data$analyzed_reads[i] <- strsplit(sampleFile[grep('Num_analyzed_reads\t', sampleFile)], '\t')[[1]][2]
			data$Hotspots[i] <- strsplit(sampleFile[grep('Num_hotspots2\t', sampleFile)], '\t')[[1]][2]
			data$SPOT[i] <- strsplit(sampleFile[grep('SPOT2\t', sampleFile)], '\t')[[1]][2]
			data$Age[i] <- df[grep(data$DSnumber[i], df$GroupID), "Age"]
			data$Uni[i] <- df[grep(data$DSnumber[i], df$GroupID), "Institution"]
			#if there's data in the Variable column, add the information. 
			if(is.null(df$Variable)) {
				data$Variable[i] <- " "
			} else {
				#MTM -- what does this do? Is he filling in info from samples with same DS num?
				data$Variable[i] <- df[grep(data$DSnumber[i], df$GroupID),]$Variable[1]
			}
			if(data$Variable[i] != "Duke") {
				#Divide samples based on category
				if(length(grep('^f[A-Z]', SampleID))>0){data$Variable[i] <- 'Fetal_Roadmap'}
				if(length(grep('^CD|^Th|^TH|^hTH|^hTR|^iTH|^th|^GM1|^GM06990|^Jurkat', SampleID))>0){data$Variable[i] <- 'Hematopoietic_lineage'}
				if(length(grep('^H1|^H7|ES|^H[0-9]|^iPS', SampleID))>0){data$Variable[i] <- 'Pluripotent'}
				if(length(grep('testis|spinal|Skin|bladder|urothelia|ventriculus|colon|limb|placenta|heart|cortex|kidney|bone|pancrea|cardia|eye|renal|gonad|muscle|osteo|medulla|brain|ovary|olfact|uteru|fibroblast|lung|tongue|bowel|putamen|esopha|gastro|ammon|derm|nucleus|gast|glom|gyrus|thyroid|adipo|neuron|prostate|intest|medull|[Ll]iver|aggregated_lymphoid_nodule|aorta|artery|psoas|stomach|testes|tibial_artery|vagina|omental_fat_depot|fetal_umbilical_cord|pons|medial_popliteal_nerve|globus|Spleen', SampleID[grep('^f[A-Z]', SampleID, invert=T)], ignore.case=T, invert=F))>0){data$Variable[i] <- 'Tissues'}
			}
			#rm(sampleFile)
			#gc()
		}
}


#Fix sample age. 
#Only keep first entry e.g. male (week 7) male (week8)
data$Age <- gsub(').*', '', data$Age)
data$Age[grep('day', data$Age)] <- paste0(round(as.numeric(gsub(' day| days', '', data$Age[grep('day', data$Age)]))/7), ' weeks')
data$Age[grep('^8 ', data$Age)] <- '08 weeks' 

##Rename the fetal tissues 
#fetalRename <- data[grep('day|week', data$Age),] %>% 
#	 filter(!grepl('^f', cellType)) %>%
#	 filter(!grepl('AG04449|AG04450|IMR_90', cellType)) %>%
#	 select(DSnumber)
#	 
#library(Hmisc)
#data[data$DSnumber%in%fetalRename$DSnumber,]$cellType <- gsub('^', 'f', capitalize(data[data$DSnumber%in%fetalRename$DSnumber,]$cellType))
#data[data$DSnumber%in%fetalRename$DSnumber,]$Variable <- 'Fetal_roadmap'


#Add replicate numbers based on highest number of non redundant reads
data <- data[!is.na(data$cellType),]
Replicates <- unique(subset(data, select=c(cellType, Variable)))
Replicates$cellType <- gsub('_L$|_R$', '', Replicates$cellType)
for (i in 1:nrow(Replicates)){
	#cat(Replicates$cellType[i], Replicates$Variable[i])
	data[gsub('_L$|_R$', '', data$cellType)==Replicates$cellType[i] & data$Variable==Replicates$Variable[i],]$Replicate <- rank(-as.numeric(data[gsub('_L$|_R$', '', data$cellType)==Replicates$cellType[i] & data$Variable==Replicates$Variable[i],]$analyzed_reads))
	#Fix color for Left and right tissues
	data[gsub('_L$|_R$', '', data$cellType)==Replicates$cellType[i] & data$Variable==Replicates$Variable[i],]$Color <- names(tail(table(data[gsub('_L$|_R$', '', data$cellType)==Replicates$cellType[i],]$Color), 1))
}

data$Variable[data$Variable=='UMass'] <- 'UW'
data$Uni[data$Uni=='UMass'] <- 'UW'

#Delete to restore from 20180129
dataRep <- data
dataRep$cellType <- gsub('_L$|_R$', '', dataRep$cellType)
dataRep <- split(dataRep, dataRep$cellType)
#repDF <- dataRep[[3]]
dataReplist <- lapply(dataRep, function(repDF) {
	head(repDF)
	#repDF$Replicate <- NA
	if ( nrow(repDF[grep('UW', repDF$Uni),]) >0 ) {
		 repDF[grep('UW', repDF$Uni),]$Replicate <- rank(-as.numeric(repDF[grep('UW', repDF$Uni),]$analyzed_reads))
		 repDF[grep('UW', repDF$Uni, invert=T),]$Replicate <- rank(-as.numeric(repDF[grep('UW', repDF$Uni , invert=T),]$analyzed_reads)) + nrow(repDF[grep('UW', repDF$Uni),])
	} else { repDF$Replicate <- rank(-as.numeric(repDF$analyzed_reads))
	}
	repDF
}) 

data <- do.call(rbind, dataReplist)

data[data$Replicate>2,]$Replicate <- 'Other'
#data$Variable[data$Variable=='UMass'] <- 'UW'
data[is.na(data$Age),]$Age <- 'NoAge'
data$Replicate <- paste0('rep', data$Replicate)
#Output file
write.table(subset(data, select=-c(Uni)), file=opt$out, sep='\t', col.names=T, row.names=F, quote=F)


warnings()

cat("Done!!!")
