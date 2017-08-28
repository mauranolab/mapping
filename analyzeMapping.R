library(reshape)
library(ggplot2)
library(stringr)
library(RColorBrewer)

#list the flowcells based on aligned.
flowCells<-list.files('./',pattern='aligned')
flowCells<-flowCells[grep('bak',flowCells,invert=T)]

#Define functions for splitting the lines
splitLines <- function(Term,sep) {
data.frame(do.call('rbind', strsplit(as.character(Sample[grep(Term,Sample,fixed=T)]),sep,fixed=TRUE)),stringsAsFactors=F)}

#As above but +1 line
getLength<- function(Term) {
data.frame(do.call('rbind', strsplit(as.character(Sample[grep(Term,Sample,fixed=T)+1]),' ',fixed=TRUE)),stringsAsFactors=F)}

getUMIHist<- function(Sample) {
data.frame(do.call('rbind', strsplit(as.character(Sample[grep('Histogram of number of reads per barcode+UMI',Sample,fixed=T)+1:10]),' ',fixed=TRUE)),stringsAsFactors=F)$X1}

getBCHist<- function(Sample) {
data.frame(do.call('rbind', strsplit(as.character(Sample[grep('Histogram of number of reads per barcode',Sample,fixed=T)+1:10]),' ',fixed=TRUE)),stringsAsFactors=F)$X1}

getInsert<- function(Sample) {
data.frame(do.call('rbind', strsplit(as.character(Sample[grep('Histogram of number of insertion sites per barcode',Sample,fixed=T)+1:2]),' ',fixed=TRUE)),stringsAsFactors=F)$X1}

getChrom <- function(Sample) {
data.frame(do.call('rbind', strsplit(as.character(Sample[grep('Integration sites by chrom',Sample,fixed=T)+1:50]),' ',fixed=TRUE)),stringsAsFactors=F)[,c(1:2)]}

getReadsPerSite<-function(Sample) {
data.frame(do.call('rbind', strsplit(as.character(Sample[grep('Histogram of number of insertion sites per barcode',Sample,fixed=T)+4:15]),' ',fixed=TRUE)),stringsAsFactors=F)[,c(1:2)]}

getGeneLoc <- function(Sample) {
data.frame(do.call('rbind', strsplit(as.character(Sample[grep('doing GenicLocation',Sample,fixed=T)+1:50]),' ',fixed=TRUE)),stringsAsFactors=F)[,c(1:2)]}

SummarizeFlowcell<-list()
flowCells<-flowCells[grep('MaxDNA',flowCells,invert=T)]
#flowCells<-flowCells[8]
for (flow in 1:length(flowCells)){# change to 1:length(flowCells) when all flowcells are done
    cat(flowCells[flow],'\n') 
    SummarizeFlowcell[[flow]]<-flowCells[flow]
    FlowSamples<-list.dirs(flowCells[flow])
    FlowSamples<-FlowSamples[-1]
    Barcodes<-list.files(flowCells[flow],pattern='^BS00',include.dirs = FALSE)
    Barcodes<-Barcodes[grep('minReads|BS00118A_19_20_121A|iPCRvsDNA',Barcodes,invert=T)]#REMOVE WHEN FILES ARE FINISHED
    Barcodes<-Barcodes[grep('\\.o[0-9]', Barcodes)]
    Barcodes<-Barcodes[file.info(paste0(flowCells[flow],'/',Barcodes))$size>200]#Removes files with errors
    
    SumBarcodes <- list()
    sumBarcode<-as.data.frame(matrix(ncol=19,nrow=length(Barcodes))) 
    colnames(sumBarcode)<-c('Flowcell','BSnumber','Name','sampleType','Total reads','Total barcodes','BC+UMI','UMI length','Unique BC','BC length','Complexity','PF reads','QC mapped reads','Prop PF','Unique sites','Unique sites +-5bp','#BC 1 site','#BC 2+ sites','Prop 2+')
    #Create files for plotting
    allHistograms<-list()
    allChromosome<-list()
    allReadsperSite<-list()
    allDistToTSS<-list()
    allGeneLoc<-list()
    
    for (files in 1:length(Barcodes)){cat(Barcodes[files],'\n')
        #####
        #ANALYZE THE BARCODES
        #####
        #Read in Sample
		Sample<-readLines(paste0(flowCells[flow],'/',Barcodes[files]))
        #
        #Remove leading whitespace
        Sample<-sub("^\\s+","",Sample)
        
        #Sample name
       sumBarcode[files,]$BSnumber<-gsub('-.*','',Barcodes[files])
       sumBarcode[files,]$Name<-gsub('_RNA$|_DNA$|_iPCR$','',gsub('\\.o.*|.*-','',Barcodes[files]))
       if(length(grep('Merged',Barcodes[files]))==1){sumBarcode[files,]$sampleType<-gsub('\\.o.*','',gsub("([^_]*_[^_]*)*$", "\\1", Barcodes[files]))}
       else{sumBarcode[files,]$sampleType<-gsub('.*_|\\.o.*','',Barcodes[files])}
       #sumBarcode[files,]$sampleType<-gsub('\\.o.*','',gsub("([^_]*_[^_]*)*$", "\\1", Barcodes[files]))
       sumBarcode[files,]$Flowcell<-flowCells[flow] 
       #Number of total reads
       sumBarcode[files,][,5]<-splitLines('Number of total reads','\t')$X3
       
       #Flowcell name
       #gsub('aligned.','',flowCells[flow])
       
       #Total number of read bacodes
       sumBarcode[files,][,6]<-splitLines('Number of total read barcodes','\t')$X3
       
       #If no UMI information exists
       if (length(grep('No UMIs found',Sample))==0){
           #Number of unique BC+UMI
           sumBarcode[files,][,7]<-splitLines('Number of unique barcodes+UMI','\t')$X3
           #Number of unique BC
           sumBarcode[files,][,9]<-splitLines('Number of unique barcodes','\t')[2,]$X3
           #Check if the line after UMI length DOES NOT contains the number of unique BC+UMI
           if(!is.na(getLength('UMI lengths')$X2)) {sumBarcode[files,][,8]<-getLength('UMI lengths')$X2}
              else if (getLength('UMI lengths')$X1!=sumBarcode[files,][,7]){
                     sumBarcode[files,][,6] <- strsplit(Sample[grep(sumBarcode[files,][,5],Sample,fixed=T)][2], "\\s+")[[1]][2]
                     if (is.na(strsplit(Sample[grep(sumBarcode[files,][,7],Sample,fixed=T)][2], "\\s+")[[1]][2]) && length(Sample[grep(as.numeric(sumBarcode[files,][,7])-1,Sample,fixed=T)])==0){
                     sumBarcode[files,][,8]<-NA
                     }else sumBarcode[files,][,8] <- strsplit(Sample[grep(as.numeric(sumBarcode[files,][,7])-1,Sample,fixed=T)], "\\s+")[[1]][2]
           }else sumBarcode[files,][,8]<-getLength('UMI lengths')$X2
       }else if (length(splitLines('Number of unique barcodes','\t')$X3)==2){sumBarcode[files,][,8]<-NA
       }else sumBarcode[files,][,9]<-splitLines('Number of unique barcodes','\t')$X3
       
       
       #Find number of unique barcodes 
       if (getLength('Barcode lengths')$X1!=sumBarcode[files,][,9]){
           #This happens if there's empty barcodes in there
           sumBarcode[files,][,10] <- strsplit(Sample[grep(as.character(as.numeric(sumBarcode[files,][,9])-1),Sample,fixed=T)], "\\s+")[[1]][2]
       }else sumBarcode[files,][,10] <- getLength('Barcode lengths')$X2
       

        if (length(grep('Total PF tags',Sample))==1){ 
            #Integration sample
            sumBarcode[files,][,12]<-splitLines('Total PF tags','\t')$X3
            sumBarcode[files,][,13]<-splitLines('Number of tags passing all filters and having barcodes assigned','\t')$X2
            sumBarcode[files,][,15]<-splitLines('Total uniq sites','\t')$X3[1]
            sumBarcode[files,][,16]<-splitLines('Total uniq sites (within 5 bp, ignoring strand)','\t')$X3[1]
            sumBarcode[files,][17]<-getInsert(Sample)[1]
            sumBarcode[files,][18]<-getInsert(Sample)[2]
            
        
        }
            
    }   
    SummarizeFlowcell[[flow]]<-sumBarcode
    

 }
    
#names(SummarizeFlowcell)<-flowCells
SummarizeFlowcell<-do.call(rbind.data.frame, SummarizeFlowcell)

SummarizeFlowcell$Complexity<-0
for (i in 1:nrow(SummarizeFlowcell)){
    #Complexity
    if (is.na(SummarizeFlowcell[,7][i])) {SummarizeFlowcell$Complexity[i]<-1
    } else SummarizeFlowcell$Complexity[i]<-as.numeric(SummarizeFlowcell[,9][i])/as.numeric(SummarizeFlowcell[,7][i])
    cat(as.numeric(SummarizeFlowcell[,9][i])/as.numeric(SummarizeFlowcell[,7][i]),'\n')
    
    #Prop PF
    if(!is.na(SummarizeFlowcell[,12][i])) {SummarizeFlowcell[,14][i]<-as.numeric(SummarizeFlowcell[,13][i])/(as.numeric(SummarizeFlowcell[,5][i]))}
    #if(!is.na(SummarizeFlowcell[,12][i])) {SummarizeFlowcell[,14][i]<-as.numeric(SummarizeFlowcell[,13][i])/(as.numeric(SummarizeFlowcell[,12][i])+as.numeric(SummarizeFlowcell[,13][i]))}
    if(!is.na(SummarizeFlowcell[,12][i])) {SummarizeFlowcell[,19][i]<-as.numeric(SummarizeFlowcell[,18][i])/(as.numeric(SummarizeFlowcell[,17][i])+as.numeric(SummarizeFlowcell[,18][i]))}
}
    
for (i in 5:ncol(SummarizeFlowcell)){
SummarizeFlowcell[,i]<-format(as.numeric(SummarizeFlowcell[,i]),big.mark=",", trim=TRUE)
}  
    
write.table(SummarizeFlowcell,file=paste0('SummarisedFlowcells/SummarisedFlowcell_',Sys.Date( ),'.tsv'),row.names=F,sep='\t',quote=F)
write.table(SummarizeFlowcell,file=paste0('~/public_html/blog/Flowcells/SummarisedFlowcell_',Sys.Date( ),'.tsv'),row.names=F,sep='\t',quote=F)






