library(reshape)
library(ggplot2)
library(stringr)
library(RColorBrewer)

#list the flowcells based on aligned.
flowCells<-list.files('./',pattern='aligned')
#Function to read lines quicker

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

SummarizeFlowcell<-list()
for (flow in 1:length(flowCells)){# change to 1:length(flowCells) when all flowcells are done
    cat(flowCells[flow],'\n') 
    SummarizeFlowcell[[flow]]<-flowCells[flow]
    FlowSamples<-list.dirs(flowCells[flow])
    FlowSamples<-FlowSamples[-1]
    Barcodes<-list.files(flowCells[flow],pattern='Analyse') 
    Barcodes<-Barcodes[file.info(paste0(flowCells[flow],'/',Barcodes))$size>200]#Removes files with errors
    SumBarcodes <- list()
    sumBarcode<-as.data.frame(matrix(ncol=8,nrow=length(Barcodes))) 
    colnames(sumBarcode)<-c('Sample','Total reads','Total barcodes','BC+UMI','UMI length','Unique BC','BC length','Flowcell')
    allHistograms<-list()
    for (files in 1:length(Barcodes)){cat(Barcodes[files],'\n')
        #Read in Sample
        Sample<-readLines(paste0(flowCells[flow],'/',Barcodes[files]))
        #Remove leading whitespace
        Sample<-sub("^\\s+","",Sample)
        
        #Sample name
        sumBarcode[files,]$Sample<-gsub('Analyse_|\\..*','',Barcodes[files])
        
        #Number of total reads
        sumBarcode[files,][,2]<-splitLines('Number of total reads','\t')$X3
        
        #Flowcell name
        sumBarcode[files,][,8]<-gsub('aligned.','',flowCells[flow])
        
        #Total number of read bacodes
        sumBarcode[files,][,3]<-splitLines('Number of total read barcodes','\t')$X3
        
        #If no UMI information exists
        if (length(grep('No UMIs found',Sample))==0){
            #Number of unique BC+UMI
            sumBarcode[files,][,4]<-splitLines('Number of unique barcodes+UMI','\t')$X3
            #Number of unique BC
            sumBarcode[files,][,6]<-splitLines('Number of unique barcodes','\t')[2,]$X3
            #Check if the line after UMI length DOES NOT contains the number of unique BC+UMI
            if (getLength('UMI lengths')$X1!=sumBarcode[files,][,4]){
                sumBarcode[files,][,4] < - strsplit(Sample[grep(sumBarcode[files,][,4],Sample,fixed=T)][2], "\\s+")[[1]][2]
            }else sumBarcode[files,][,5]<-getLength('UMI lengths')$X2
        }else sumBarcode[files,][,6]<-splitLines('Number of unique barcodes','\t')$X3
       
        #Find number of unique barcodes 
        if (getLength('Barcode lengths')$X1!=sumBarcode[files,][,6]){
            #This happens if there's empty barcodes in there
            sumBarcode[files,][,7] <- strsplit(Sample[grep(as.character(as.numeric(sumBarcode[files,][,6])-1),Sample,fixed=T)], "\\s+")[[1]][2]
        }else sumBarcode[files,][,7] <- getLength('Barcode lengths')$X2
        

        
        ####
        #Create histogram for plotting 
        ####
        Bins<-c("1","2","3",'4','5','6','7','8','9','10+')
        Hist<-as.data.frame(matrix(ncol=4,nrow=20))
        colnames(Hist)<-c('Sample','Barcode','Group','Bin')
        Hist$Sample<-gsub('Analyse_|\\..*','',Barcodes[files])
        #If there's UMIs in the Sample
        if (length(grep('No UMIs found',Sample))==0){
            umiHist<-getUMIHist(Sample)
            #Check if there are alpha characters in vector. This means that no barcodes exists more than e.g. 5 times
            if (length(grep("^[[:alpha:]]*$",umiHist))>0){
                cutOff<-grep("^[[:alpha:]]*$",umiHist)[1]-1
                Hist[,2][1:cutOff]<-umiHist[1:cutOff]
                Hist[,3][c(1:10)]<-'BC_UMI'
                Hist[,4][c(1:10)]<-Bins
            }else {Hist[,2][1:10]<-umiHist[1:10]
                    Hist[,3][c(1:10)]<-'BC_UMI'
                    Hist[,4][c(1:10)]<-Bins}
        }
        if (length(grep('No UMIs found',Sample))==1){Hist[,2][1:10]<-NA
        Hist[,3][c(1:10)]<-'BC_UMI'
        Hist[,4][c(1:10)]<-Bins}
        #Remove the barcode+UMI line if it exists
        BCSample<-Sample[grep('Histogram of number of reads per barcode+UMI',Sample,fixed=T,invert=T)]
        bcHist<-getBCHist(BCSample)
        #Check if there are alpha characters in vector. This means that no barcodes exists more than e.g. 5 times
        if (length(grep("^[[:alpha:]]*$",bcHist))>0){
            cutOff<-grep("^[[:alpha:]]*$",bcHist)[1]-1
            Hist[,2][11:(10+cutOff)]<-bcHist[1:cutOff]
            Hist[,3][c(11:20)]<-'BC'
            Hist[,4][c(11:20)]<-Bins
        }else {Hist[,2][11:20]<-bcHist[1:10]
                Hist[,3][c(11:20)]<-'BC'
                Hist[,4][c(11:20)]<-Bins}
                
        allHistograms[[files]]<-Hist
   
    }   
    SummarizeFlowcell[[flow]]<-sumBarcode
  
  PlotHist<-do.call(rbind.data.frame, allHistograms)
  PlotHist$Barcode[is.na(PlotHist$Barcode)]<-1#Since we're plotting log10 we don't want negative values
  PlotHist$Barcode<- as.numeric(PlotHist$Barcode)
  PlotHist$Bin <- factor(PlotHist$Bin , levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10+"))
  
  p<-ggplot(PlotHist,aes(x=Bin,y=Barcode,fill=Group))+
      geom_bar(stat = "identity",position=position_dodge(), colour="black")+
      facet_wrap(~Sample, nrow = 1)+
      scale_fill_brewer(palette ='Set1')+
      theme_classic()+
      theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=45,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
      ylab('Number of barcodes')+
      xlab('Bins')+
      scale_y_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
      annotation_logticks(sides = "l",long = unit(0.1, "cm"))  +
      ggtitle(flowCells[flow])
  ggsave(p, file=paste0('Histogram_',gsub('aligned.','',flowCells[flow]),'.pdf'), width=3*length(Barcodes), height=4)  

    
    
}
#names(SummarizeFlowcell)<-flowCells
SummarizeFlowcell<-do.call(rbind.data.frame, SummarizeFlowcell)

SummarizeFlowcell$Complexity<-0
for (i in 1:nrow(SummarizeFlowcell)){
    if (is.na(SummarizeFlowcell[,4][i])) {SummarizeFlowcell$Complexity[i]<-1
    } else SummarizeFlowcell$Complexity[i]<-as.numeric(SummarizeFlowcell[,6][i])/as.numeric(SummarizeFlowcell[,4][i])
    cat(as.numeric(SummarizeFlowcell[,6][i])/as.numeric(SummarizeFlowcell[,4][i]),'\n')
}
write.table(do.call(rbind.data.frame, SummarizeFlowcell),file='SummarisedFlowcell.tsv',row.names=F,sep='\t',quote=F)
