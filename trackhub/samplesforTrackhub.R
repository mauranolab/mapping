#Add variables from command line tsv files of DSnumbers and Institute
library("optparse")
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

df = read.table(opt$file,sep='\t',header=T,stringsAsFactors=F,na.strings=c(""," ","NA",'\t'), fill = TRUE)
pwd<-getwd()
cat('Dimensions of Input file: ',dim(unique(df)),'\n')

if (is.null(df$GroupID)){
  stop("GroupID column is required", call.=FALSE)
} else {cat('Number of unique groups:', length(unique(df$GroupID)),'\n')}


bwfiles<-list.files(path=pwd,pattern='hg38.bw')

cat('Number of unique group that has a bigWig file: ',length(bwfiles[grep(paste(unique(df$GroupID),collapse="|"),bwfiles)]),'\n')
bwfiles<-bwfiles[grep(paste(unique(df$GroupID),collapse="|"),bwfiles)]
mergedFiles<-list.files(path=pwd, pattern="merge.*")



#pick colours
library(RColorBrewer)
set.seed(12345)
col<-c(brewer.pal(n=9,'Set1'),brewer.pal(n=8,'Dark2'),brewer.pal(n=12,'Paired')[c(FALSE, TRUE)])
col<-gsub('#FFFF33','#F0027F',col)
col<-gsub('#FFFFFF','#BEAED4',col)

Groups<-gsub("-.*",'',bwfiles)
col<-rep(col,round(length(Groups)/length(col),2))[as.factor(unique(Groups))]        
col<-as.data.frame(col,unique(Groups))

data<-data.frame(matrix(ncol=10,nrow=length(bwfiles)))
colnames(data)<-c('cellType','DSnumber','Replicate','Color','Assay','nonredundant_tags','SPOT','Hotspots','Exclude','Variable')
for (i in 1:length(bwfiles)){
       SampleID<-gsub('.hg38.bw','',bwfiles[i])
       if (length(grep(paste0('merge.',SampleID),mergedFiles))>0)
              sampleFile<-readLines(paste0('/vol/mauranolab/public_html/encode/dnase/mapped/',mergedFiles[grep(paste0('merge.',SampleID),mergedFiles)]))
              if(tail(sampleFile)[1]=='Done!'){
                     colGroup<-col2rgb(as.character(col[rownames(col)%in%gsub("-.*",'',bwfiles[i]),][1]))[,1]
                     data$cellType[i]<-gsub('-.*','',SampleID)
                     data$DSnumber[i]<-gsub('.*-','',SampleID)
                     data$Replicate[i]<-1
                     data$Color[i]<-paste(colGroup[1],colGroup[2],colGroup[3],sep=',')
                     data$Assay[i]<-'DNase'
                     data$nonredundant_tags[i]<-strsplit(sampleFile[grep('Num_nonredundant_tags\t',sampleFile)],'\t')[[1]][2]
                     data$Hotspots[i]<-strsplit(sampleFile[grep('Num_hotspots\t',sampleFile)],'\t')[[1]][2]
                     data$SPOT[i]<-strsplit(sampleFile[grep('SPOT\t',sampleFile)],'\t')[[1]][2]
                     
                     #if there's data in the Variable column, add the information. 
                     if(is.null(df$Variable)){data$Variable[i]<-" "
                     } else {data$Variable[i]<-df[grep(data$DSnumber[i],df$GroupID),]$Variable[1]}
                     #Divide samples based on category
                     if(length(grep('^f[A-Z]',SampleID))>0){data$Variable[i]<-'Fetal_roadmap'}
                     if(length(grep('^CD|^Th|^TH|^hTH|^iTH|^th',SampleID))>0){data$Variable[i]<-'Hematopoietic_lineage'}
                     if(length(grep('^H1|^H7|ES|^H[0-9]|^iPS',SampleID))>0){data$Variable[i]<-'Pluripotent'}
                     if(length(grep('testis|spinal|Skin|bladder|urothelia|ventriculus|colon|limb|placenta|heart|cortex|kidney|bone|oesteoblast|pancrea|cardia|eye|renal|gonad|muscle|osteo|medulla|brain|ovary|olfact|uteru|fibroblast|lung|tongue|bowel|putamen|esopha|gastro|ammon|derm|nucleus|gast|glom|gyrus|thyroid|adipo|neuron|prostate|intest|medull',
                     SampleID[grep('^f[A-Z]',SampleID,invert=T)],ignore.case=T,invert=F))>0){data$Variable[i]<-'Tissues'}
                     
                     
              }      
}


#Add replicate numbers based on highest number of non redundant tags
data<-data[!is.na(data$cellType),]
Replicates<-unique(data$cellType)
for (i in 1:length(Replicates)){
       Order<-order(data[grep(Replicates[i],data$cellType),]$nonredundant_tags,decreasing=T)
       data[grep(Replicates[i],data$cellType),]$Replicate<-rank(-as.numeric(data[grep(Replicates[i],data$cellType),]$nonredundant_tags))
       
}
data[data$Replicate>2,]$Replicate<-'Other'
data$Replicate<-paste0('rep',data$Replicate)
data$Replicate<-gsub('repOther','Other',data$Replicate)
data[data$Hotspots=='NA',]$Hotspots<-0
data$Variable[data$Variable=='UMass']<-'Roadmap'


#Output file
write.table(data, file=opt$out,sep='\t',col.names=T,row.names=F,quote=F)
