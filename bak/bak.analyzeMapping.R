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
    colnames(sumBarcode)<-c('Flowcell','BSnumber','Sample','Type','Total reads','Total barcodes','BC+UMI','UMI length','Unique BC','BC length','Complexity','PF reads','QC mapped reads','Prop PF','Unique sites','Unique sites +-5bp','#BC 1 site','#BC 2+ sites','Prop 2+')
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
       sumBarcode[files,]$Sample<-gsub('\\.o.*','',Barcodes[files])
       sumBarcode[files,]$Flowcell<-flowCells[flow] 
       sumBarcode[files,]$Sample<-gsub('\\.o.*','',Barcodes[files])
       #Number of total reads
       sumBarcode[files,][,3]<-splitLines('Number of total reads','\t')$X3
       
       #Flowcell name
       #gsub('aligned.','',flowCells[flow])
       
       #Total number of read bacodes
       sumBarcode[files,][,4]<-splitLines('Number of total read barcodes','\t')$X3
       
       #If no UMI information exists
       if (length(grep('No UMIs found',Sample))==0){
           #Number of unique BC+UMI
           sumBarcode[files,][,5]<-splitLines('Number of unique barcodes+UMI','\t')$X3
           #Number of unique BC
           sumBarcode[files,][,7]<-splitLines('Number of unique barcodes','\t')[2,]$X3
           #Check if the line after UMI length DOES NOT contains the number of unique BC+UMI
           if(!is.na(getLength('UMI lengths')$X2)) {sumBarcode[files,][,6]<-getLength('UMI lengths')$X2}
              else if (getLength('UMI lengths')$X1!=sumBarcode[files,][,5]){
                     sumBarcode[files,][,6] <- strsplit(Sample[grep(sumBarcode[files,][,5],Sample,fixed=T)][2], "\\s+")[[1]][2]
                     if (is.na(strsplit(Sample[grep(sumBarcode[files,][,5],Sample,fixed=T)][2], "\\s+")[[1]][2]) && length(Sample[grep(as.numeric(sumBarcode[files,][,5])-1,Sample,fixed=T)])==0){
                     sumBarcode[files,][,6]<-NA
                     }else sumBarcode[files,][,6] <- strsplit(Sample[grep(as.numeric(sumBarcode[files,][,5])-1,Sample,fixed=T)], "\\s+")[[1]][2]
           }else sumBarcode[files,][,6]<-getLength('UMI lengths')$X2
       }else if (length(splitLines('Number of unique barcodes','\t')$X3)==2){sumBarcode[files,][,6]<-NA
       }else sumBarcode[files,][,7]<-splitLines('Number of unique barcodes','\t')$X3
       
       
       #Find number of unique barcodes 
       if (getLength('Barcode lengths')$X1!=sumBarcode[files,][,7]){
           #This happens if there's empty barcodes in there
           sumBarcode[files,][,8] <- strsplit(Sample[grep(as.character(as.numeric(sumBarcode[files,][,7])-1),Sample,fixed=T)], "\\s+")[[1]][2]
       }else sumBarcode[files,][,8] <- getLength('Barcode lengths')$X2
       

        if (length(grep('Total PF tags',Sample))==1){ 
            #Integration sample
            sumBarcode[files,][,10]<-splitLines('Total PF tags','\t')$X3
            sumBarcode[files,][,11]<-splitLines('Number of tags passing all filters and having barcodes assigned','\t')$X2
            sumBarcode[files,][,13]<-splitLines('Total uniq sites','\t')$X3[1]
            sumBarcode[files,][,14]<-splitLines('Total uniq sites (within 5 bp, ignoring strand)','\t')$X3[1]
            sumBarcode[files,][15]<-getInsert(Sample)[1]
            sumBarcode[files,][16]<-getInsert(Sample)[2]
            
        
            ###
            #Barplot of insertions per chromosome
            ###
            ChromInsert<-getChrom(Sample)
            ChromInsert<-ChromInsert[grep("chr",ChromInsert$X2),]
            ChromInsert$Sample<-gsub('\\.o.*','',Barcodes[files])
            ChromInsert$Flowcell<-flowCells[flow]
            allChromosome[[files]]<-ChromInsert
            
            
            ####
            #Barplot of reads per site
            ####
            ReadsSite<-getReadsPerSite(Sample)
            ReadsSite<-ReadsSite[grep("^[[:alpha:]]*$|[[:punct:]]",ReadsSite$X1,invert=T),]
            ReadsSite$Sample<-gsub('\\.o.*','',Barcodes[files])
            ReadsSite$Flowcell<-flowCells[flow]
            allReadsperSite[[files]]<-ReadsSite
            
            ####
            #Distance to refseq TSS
            #####
            #DistToTSS<-read.table(paste0(flowCells[flow],'/',gsub('\\.o.*','',gsub('merge-analyze-','',Barcodes[files])),'/DistToTSS.txt'))
            #colnames(DistToTSS) <- c("DistToTSS")
            #DistToTSS$Bins <- cut(DistToTSS$DistToTSS, breaks=c(-10^7, -10^6, -10^5, -10^4, -2500, 0, 5000, 10^4, 10^5, 10^6, 10^7), right=F, include.lowest=T, labels=c("-10 Mb", "-1 Mb", "-100 kb", "-10 kb", "-2.5 kb", "+5 kb", "+10 kb", "+100 kb", "+1 Mb", "+10 Mb"))
            #colnames(DistToTSS) <- c("DistToTSS")
            #DistToTSS$Sample<-gsub('\\.o.*','',Barcodes[files])
            #DistToTSS$Flowcell<-flowCells[flow]
            #allDistToTSS[[files]]<-DistToTSS
            
            ####
            #Genomic feature
            ####
            GeneLoc<-getGeneLoc(Sample)
            GeneLoc<-GeneLoc[grep("^[[:alpha:]]*$|[[:punct:]]|NA",GeneLoc$X1,invert=T),]
            GeneLoc<-GeneLoc[!is.na(GeneLoc$X1),]
            GeneLoc$Sample<-gsub('\\.o.*','',Barcodes[files])
            GeneLoc$Flowcell<-flowCells[flow]
            allGeneLoc[[files]]<-GeneLoc
        }
            
        ####
        #Create histogram of Barcodes  
        ####
        Bins<-c("1","2","3",'4','5','6','7','8','9','10+')
        Hist<-as.data.frame(matrix(ncol=4,nrow=20))
        colnames(Hist)<-c('Sample','Barcode','Group','Bin')
        Hist$Sample<-gsub('\\.o.*','',Barcodes[files])
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
    
    
    
    ####
    #PLOTTING
    #####
    #Barcodes
    allHistograms<-allHistograms[lapply(allHistograms,length)>0]
    PlotHist<-do.call(rbind.data.frame, allHistograms)
    if(length(allHistograms)>6) plotWidth <- length(allHistograms)/2 else plotWidth <-length(allHistograms)
    if(length(allHistograms)>6) plotRows<- 2 else plotRows <-1  
    
    PlotHist$Barcode[is.na(PlotHist$Barcode)]<-1#Since we're plotting log10 we don't want negative values
    PlotHist$Barcode<- as.numeric(PlotHist$Barcode)
    PlotHist$Bin <- factor(PlotHist$Bin , levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10+"))
    
    #p<-ggplot(PlotHist,aes(x=Bin,y=Barcode,fill=Group))+
    #    geom_bar(stat = "identity",position=position_dodge(), colour="black")+
    #    facet_wrap(~Sample, nrow = plotRows)+
    #    scale_fill_brewer(palette ='Set1')+
    #    theme_classic()+
    #    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=45,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
    #    ylab('Number of barcodes')+
    #    xlab('Bins')+
    #    scale_y_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
    #    annotation_logticks(sides = "l",long = unit(0.1, "cm"))  +
    #    ggtitle(flowCells[flow])
    #ggsave(p, file=paste0(gsub('aligned.','',flowCells[flow]),'_Barcodes','.pdf'), width=3*plotWidth, height=4) 
    
    #    allHistograms<-list()
    #allChromosome<-list()
    #allReadsperSite<-list()
    #allDistToTSS<-list()
    #allGeneLoc<-list()
    #
    if (length(allChromosome[lapply(allChromosome,length)>0])){ 
        ####
        #Insertion per chromosome
        #### 
        allChromosome<-allChromosome[lapply(allChromosome,length)>0]
        
        if(length(allChromosome)>6) plotWidth <- length(allChromosome)/2 else plotWidth <-length(allChromosome)
        if(length(allChromosome)>6) plotRows<- 2 else plotRows <-1  
        
        PlotChrom<-do.call(rbind.data.frame, allChromosome)   
        colnames(PlotChrom)<-c('Insertions','Chromosome','Sample','Flowcell')
        PlotChrom$Insertions<- as.numeric(PlotChrom$Insertions)
        PlotChrom$Chromosome <- factor(PlotChrom$Chromosome , levels = c(paste0('chr',1:22),'chrY','chrX','chrM',unique(PlotChrom[,2])[grep('_',unique(PlotChrom[,2]))]))
       
        #p<-ggplot(PlotChrom,aes(x=Chromosome,y=Insertions,fill=Insertions))+
        #    geom_bar(stat = "identity")+
        #    facet_wrap(~Sample, nrow = plotRows)+
        #    #scale_fill_brewer(palette ='PuBu')+
        #    theme_classic()+
        #    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=6,angle=45,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
        #    ylab('Number of insertions')+
        #    xlab('Chromosomes')+
        #    #scale_y_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
        #    #annotation_logticks(sides = "l",long = unit(0.1, "cm"))  +
        #    ggtitle(flowCells[flow])
        #ggsave(p, file=paste0(gsub('aligned.','',flowCells[flow]),'_Chromosomes','.pdf'), width=3*plotWidth, height=4) 
        
        ####
        #Reads per site
        #### 
        allReadsperSite<-allReadsperSite[lapply(allReadsperSite,length)>0]
        
        if(length(allReadsperSite)>6) plotWidth <- length(allReadsperSite)/2 else plotWidth <-length(allReadsperSite)
        if(length(allReadsperSite)>6) plotRows<- 2 else plotRows <-1  
        
        PlotReadsSite<-do.call(rbind.data.frame, allReadsperSite)   
        colnames(PlotReadsSite)<-c('Reads','Site','Sample','Flowcell')
        PlotReadsSite$Reads<- as.numeric(PlotReadsSite$Reads)
        PlotReadsSite$Site <- factor(PlotReadsSite$Site , levels = unique(PlotReadsSite$Site))
       
        #p<-ggplot(PlotReadsSite,aes(x=Site,y=Reads,fill=Reads))+
        #    geom_bar(stat = "identity")+
        #    facet_wrap(~Sample, nrow = plotRows)+
        #    #scale_fill_brewer(palette ='PuBu')+
        #    theme_classic()+
        #    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=45,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
        #    ylab('Barcodes')+
        #    xlab('Reads per barcode')+
        #    #scale_y_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
        #    #annotation_logticks(sides = "l",long = unit(0.1, "cm"))  +
        #    ggtitle(flowCells[flow])
        #ggsave(p, file=paste0(gsub('aligned.','',flowCells[flow]),'_ReadsPerBarcode','.pdf'), width=3*plotWidth, height=4) 
        
        ####
        #Distance to TSS
        #### 
        #allDistToTSS<-allDistToTSS[lapply(allDistToTSS,length)>0]
        #
        #if(length(allDistToTSS)>6) plotWidth <- length(allDistToTSS)/2 else plotWidth <-length(allDistToTSS)
        #if(length(allDistToTSS)>6) plotRows<- 2 else plotRows <-1  
        #
        #PlotDistToTSS<-do.call(rbind.data.frame, allDistToTSS)   
        #colnames(PlotDistToTSS)<-c('Distance','Bins','Sample','Flowcell')
        #PlotDistToTSS$Distance<- as.numeric(PlotDistToTSS$Distance)
#
        #
#
        ##split into minus and plus for grouping
        #PlotDistToTSS_min<-PlotDistToTSS[grep('-',PlotDistToTSS$Bins,fixed=T),]
        #PlotDistToTSS_min$Direction<-'-'
        #PlotDistToTSS_plus<-PlotDistToTSS[grep('+',PlotDistToTSS$Bins,fixed=T),]
        #PlotDistToTSS_plus$Direction<-'+'
        #PlotDistToTSS<-rbind(PlotDistToTSS_plus,PlotDistToTSS_min)
        #PlotDistToTSS$Bins<-gsub('-|\\+','',PlotDistToTSS$Bins)
        #PlotDistToTSS$Bins<-as.character(PlotDistToTSS$Bins)
        #PlotDistToTSS$Bins<-gsub(" ",'',PlotDistToTSS$Bins)
        #PlotDistToTSS$Bins <- factor(PlotDistToTSS$Bins , levels = c('2.5kb','5kb','10kb','100kb','1Mb','10Mb'))
        #PlotDistToTSS<-PlotDistToTSS[,-c(1)]
       
       
               
        ####
        #Density graph of distance to TSS 
        ####
		#<-ggplot(PlotDistToTSS,aes(abs(Distance),fill=Direction))+
		#	geom_density(alpha=0.3)+
		#	facet_wrap(~Sample, nrow = 1)+
		#	theme_classic()+
		#	scale_fill_brewer(palette='Set1')+
		#	ggtitle(flowCells[flow])+
		#	scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
        #	annotation_logticks(sides = "b",long = unit(0.2, "cm"))  
		#gsave(p, file=paste0('TSSdensity',gsub('aligned.','',flowCells[flow]),'.pdf'), width=3*length(Barcodes), height=4) 

               
        #####
        #Barchart graph of distance to TSS
        ######
        #library(plyr)
      	#PlotDistToTSS<-count(PlotDistToTSS, vars=c('Bins','Sample','Flowcell','Direction'))
       	#colnames(PlotDistToTSS)<-c('Bins','Sample','Flowcell','Direction','Frequency')
       #
        #p<-ggplot(PlotDistToTSS,aes(x=Bins,y=Frequency,fill=Direction))+
        #    geom_bar(stat = "identity",position=position_dodge(), colour="black")+
        #    facet_wrap(~Sample, nrow = plotRows)+
        #    scale_fill_brewer(palette ='Set1')+
        #    theme_classic()+
        #    theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=45,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
        #    ylab('Number of insertions')+
        #    xlab('Sites')+
        #    #scale_y_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
        #    #annotation_logticks(sides = "l",long = unit(0.1, "cm"))  +
        #    ggtitle(flowCells[flow])
        #ggsave(p, file=paste0(gsub('aligned.','',flowCells[flow]),'_DistanceToTSS','.pdf'), width=3*plotWidth, height=4) 
        
        
       #####
       #Barchart over genomic location
       #####
        allGeneLoc<-allGeneLoc[lapply(allGeneLoc,length)>0]
       
        if(length(allGeneLoc)>6) plotWidth <- length(allGeneLoc)/2 else plotWidth <-length(allGeneLoc)
        if(length(allGeneLoc)>6) plotRows<- 2 else plotRows <-1  
        PlotallGeneLoc<-do.call(rbind.data.frame, allGeneLoc)   
        colnames(PlotallGeneLoc)<-c('Number','Genomic_region','Sample','Flowcell')
       # PlotallGeneLoc$Number<-as.numeric(PlotallGeneLoc$Number)
       # PlotallGeneLoc[,2] <- factor(PlotallGeneLoc[,2]  , levels = c('TxStart-10kb','TxS-1000','5UTR','coding','intron','3UTR','intergenic'))
       # 
	#	p<-ggplot(PlotallGeneLoc,aes(x=Genomic_region,y=Number,fill=Genomic_region))+
       # 	geom_bar(stat = "identity")+
       # 	facet_wrap(~Sample, nrow = plotRows)+
       # 	scale_fill_brewer(palette ='Pastel1')+
       # 	theme_classic()+
       # 	theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=45,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
       # 	ylab('Number of insertions')+
       # 	xlab('Genomic locations')+
       # 	scale_y_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
       # 	annotation_logticks(sides = "l",long = unit(0.1, "cm"))  +
       # 	ggtitle(flowCells[flow])
    	#ggsave(p, file=paste0(gsub('aligned.','',flowCells[flow]),'_GenomicLoc','.pdf'), width=3*plotWidth, height=4) 
        
	}
 }
    
#names(SummarizeFlowcell)<-flowCells
SummarizeFlowcell<-do.call(rbind.data.frame, SummarizeFlowcell)

SummarizeFlowcell$Complexity<-0
for (i in 1:nrow(SummarizeFlowcell)){
    #Complexity
    if (is.na(SummarizeFlowcell[,5][i])) {SummarizeFlowcell$Complexity[i]<-1
    } else SummarizeFlowcell$Complexity[i]<-as.numeric(SummarizeFlowcell[,7][i])/as.numeric(SummarizeFlowcell[,5][i])
    cat(as.numeric(SummarizeFlowcell[,7][i])/as.numeric(SummarizeFlowcell[,5][i]),'\n')
    
    #Prop PF
    if(!is.na(SummarizeFlowcell[,10][i])) {SummarizeFlowcell[,12][i]<-as.numeric(SummarizeFlowcell[,11][i])/(as.numeric(SummarizeFlowcell[,10][i])+as.numeric(SummarizeFlowcell[,11][i]))}
    if(!is.na(SummarizeFlowcell[,10][i])) {SummarizeFlowcell[,17][i]<-as.numeric(SummarizeFlowcell[,16][i])/(as.numeric(SummarizeFlowcell[,15][i])+as.numeric(SummarizeFlowcell[,16][i]))}
}
    
write.table(SummarizeFlowcell,file=paste0('SummarisedFlowcell_',Sys.Date( ),'.tsv'),row.names=F,sep='\t',quote=F)






