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
pwd<-getwd()
cat('Dimensions of Input file: ',dim(unique(df)),'\n')

if (is.null(df$GroupID)){
  stop("GroupID column is required", call.=FALSE)
} else {cat('Number of unique groups:', length(unique(df$GroupID)),'\n')}


bwfiles<-list.files(path='./',pattern='hg38')
#bwfiles<-bwfiles[grep('DS23661|DS24771',bwfiles,invert=T)]

cat('Number of unique group that has a bigWig file: ',length(bwfiles[grep(paste(unique(df$GroupID),collapse="|"),bwfiles)]),'\n')
bwfiles<-bwfiles[grep(paste(unique(df$GroupID),collapse="|"),bwfiles)]

bwfiles <- bwfiles[grep('bak',bwfiles,invert=T)]
#pick colours

set.seed(12345)
col<-c(brewer.pal(n=9,'Set1'),brewer.pal(n=8,'Dark2'),brewer.pal(n=12,'Paired')[c(FALSE, TRUE)])
col<-gsub('#FFFF33','#F0027F',col)
col<-gsub('#FFFFFF','#BEAED4',col)

Groups<-gsub("-.*",'',bwfiles)
col<-rep(col,round(length(Groups)/length(col),2))[as.factor(unique(Groups))]        
col<-as.data.frame(col,unique(Groups))

data<-data.frame(matrix(ncol=12,nrow=length(bwfiles)))
colnames(data)<-c('cellType','DSnumber','Replicate','Color','Assay','nonredundant_tags','SPOT','Hotspots','Exclude','Variable','Age', 'Uni')

#Change cellType names for Fetal tissues


for (i in 1:length(bwfiles)){
       SampleID<-bwfiles[i]
       cat(SampleID,'\n')
       mergedFiles<-list.files(path=paste0(pwd, '/',SampleID), pattern="^makeTracks.*")
       #cat(mergedFiles, '\n')
       if (length(mergedFiles)>0)
              sampleFile<-readLines(paste0(pwd, '/', SampleID,'/',mergedFiles), n=2000)
             # if(tail(sampleFile)[2]=='Done!'){
                     colGroup<-col2rgb(as.character(col[rownames(col)%in%gsub("-.*",'',bwfiles[i]),][1]))[,1]
                     data$cellType[i]<-gsub('-.*','',SampleID)
                     data$DSnumber[i]<-gsub('.hg38','',gsub('.*-','',SampleID))
                     data$Replicate[i]<-1
                     data$Color[i]<-paste(colGroup[1],colGroup[2],colGroup[3],sep=',')
                     data$Assay[i]<-'DNase'
                     data$nonredundant_tags[i]<-strsplit(sampleFile[grep('Num_uniquely_mapped_reads\t',sampleFile)],'\t')[[1]][2]
                     data$Hotspots[i]<-strsplit(sampleFile[grep('Num_hotspots2\t',sampleFile)],'\t')[[1]][2]
                     data$SPOT[i]<-strsplit(sampleFile[grep('SPOT2\t',sampleFile)],'\t')[[1]][2]
                     data$Age[i] <- df$Age[grep(gsub('.hg38','',gsub('.*-','',SampleID)),df$GroupID)]
                     data$Uni[i] <- df[grep(data$DSnumber[i],df$GroupID),]$Variable[1]#ADDED 20180129
                     #if there's data in the Variable column, add the information. 
                     if(is.null(df$Variable)){data$Variable[i]<-" "
                     } else {data$Variable[i]<-df[grep(data$DSnumber[i],df$GroupID),]$Variable[1]}
                     #Divide samples based on category
                     if(length(grep('^f[A-Z]',SampleID))>0){data$Variable[i]<-'Fetal_Roadmap'}
                     if(length(grep('^CD|^Th|^TH|^hTH|^iTH|^th|^GM1',SampleID))>0){data$Variable[i]<-'Hematopoietic_lineage'}
                     if(length(grep('^H1|^H7|ES|^H[0-9]|^iPS',SampleID))>0){data$Variable[i]<-'Pluripotent'}
                     if(length(grep('testis|spinal|Skin|bladder|urothelia|ventriculus|colon|limb|placenta|heart|cortex|kidney|bone|oesteoblast|pancrea|cardia|eye|renal|gonad|muscle|osteo|medulla|brain|ovary|olfact|uteru|fibroblast|lung|tongue|bowel|putamen|esopha|gastro|ammon|derm|nucleus|gast|glom|gyrus|thyroid|adipo|neuron|prostate|intest|medull',
                     SampleID[grep('^f[A-Z]',SampleID,invert=T)],ignore.case=T,invert=F))>0){data$Variable[i]<-'Tissues'}
                     #rm(sampleFile)
                     #gc()
                     
              #}      
}

#Fix sample age. 
#Only keep first entry e.g. male (week 7) male (week8)
data$Age <- gsub(').*','',data$Age)
data$Age[grep('day',data$Age)] <- paste0(round(as.numeric(gsub(' day| days','',data$Age[grep('day',data$Age)]))/7),' weeks')
data$Age[grep('^8 ',data$Age)] <- '08 weeks' 

##Rename the fetal tissues 
#fetalRename <- data[grep('day|week',data$Age),] %>% 
#       filter(!grepl('^f',cellType)) %>%
#       filter(!grepl('AG04449|AG04450|IMR_90',cellType)) %>%
#       select(DSnumber)
#       
#library(Hmisc)
#data[data$DSnumber%in%fetalRename$DSnumber,]$cellType <- gsub('^','f',capitalize(data[data$DSnumber%in%fetalRename$DSnumber,]$cellType))
#data[data$DSnumber%in%fetalRename$DSnumber,]$Variable <- 'Fetal_roadmap'


#Add replicate numbers based on highest number of non redundant tags
data<-data[!is.na(data$cellType),]
Replicates<-unique(subset(data,select=c(cellType,Variable)))
Replicates$cellType <- gsub('_L$|_R$','',Replicates$cellType)
for (i in 1:nrow(Replicates)){
       #cat(Replicates$cellType[i], Replicates$Variable[i])
       data[gsub('_L$|_R$','',data$cellType)==Replicates$cellType[i] & data$Variable==Replicates$Variable[i],]$Replicate <- rank(-as.numeric(data[gsub('_L$|_R$','',data$cellType)==Replicates$cellType[i] & data$Variable==Replicates$Variable[i],]$nonredundant_tags))
       #Fix color for Left and right tissues
       data[gsub('_L$|_R$','',data$cellType)==Replicates$cellType[i] & data$Variable==Replicates$Variable[i],]$Color <- names(tail(table(data[gsub('_L$|_R$','',data$cellType)==Replicates$cellType[i],]$Color),1))
}

data$Variable[data$Variable=='UMass']<-'UW'
data$Uni[data$Uni=='UMass']<-'UW'

#Delete to restore from 20180129
dataRep <- data
dataRep$cellType <- gsub('_L$|_R$','',dataRep$cellType)
dataRep <- split(dataRep, dataRep$cellType)
#repDF <- dataRep[[3]]
dataReplist <- lapply(dataRep, function(repDF) {
       head(repDF)
       #repDF$Replicate <- NA
       if ( nrow(repDF[grep('UW', repDF$Uni),]) >0 ) {
               repDF[grep('UW', repDF$Uni),]$Replicate <- rank(-as.numeric(repDF[grep('UW', repDF$Uni),]$nonredundant_tags))
               repDF[grep('UW', repDF$Uni, invert=T),]$Replicate <- rank(-as.numeric(repDF[grep('UW', repDF$Uni , invert=T),]$nonredundant_tags)) + nrow(repDF[grep('UW', repDF$Uni),])
       } else { repDF$Replicate <- rank(-as.numeric(repDF$nonredundant_tags))
       }
       repDF
})               

data <- do.call(rbind, dataReplist)

data[data$Replicate>2,]$Replicate<-'Other'
#data$Variable[data$Variable=='UMass']<-'UW'
data[is.na(data$Age),]$Age <- 'NoAge'
data$Replicate<-paste0('rep',data$Replicate)
#Output file
write.table(subset(data, select=-c(Uni)), file=opt$out,sep='\t',col.names=T,row.names=F,quote=F)
