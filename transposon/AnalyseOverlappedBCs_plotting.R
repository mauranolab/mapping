#!/bin/env Rscript

#Add variables from command line tsv files of DSnumbers and Institute
library("optparse")
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(scales)
library(dplyr)
library(reshape2)
#library(plyr)
library(data.table)
library(randomForest)
library(ggjoy)
basenames <- ls()
#source("~/.Rprofile")
option_list = list(
  make_option(c("-sample", "--SampleName"), type="character", default=NULL, 
              help="Sample name", metavar="character"),
	make_option(c("-o", "--OUTDIR"), type="character", default=NULL, 
              help="output file name", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$SampleName)){
  print_help(opt_parser)
  stop("Sample name needs to be provided.\n", call.=FALSE)
}

if (is.null(opt$OUTDIR)){
  stop("Output directory needs to be provided", call.=FALSE)
}


cat('Sample name: ', opt$SampleName, '\n')
cat('Output directory: ', opt$OUTDIR, '\n')

sampleColor <- data.frame(c('CMV','Bglobin','GGlo','GGlo-HS2','GGlo-Ins-HS2','Ins_fw-GGlo_fw-Ins-HS2','Ins_rv-GGlo-Ins_rv-HS2', 'Ins-GGlo-Ins-HS2'),c('#e41a1c','#999999','#FF7F00','#a65628','#377eb8','#4daf4a','#984EA3', '#fdc086'))
colnames(sampleColor) <- c('Sample', 'Color')

if (length(as.character(sampleColor[sampleColor$Sample%in%opt$SampleName,]$Color))==1) { 
       sampleColor <- as.character(sampleColor[sampleColor$Sample%in%opt$SampleName,]$Color)
       } else {sampleColor <- 'white'
}

#opt <- data.frame(SampleName='CMV', OUTDIR='/home/maagj01/public_html/blog/OverAllAnalysis/CMV')
#set colours for plotting
myColors <- brewer.pal(4, "Set1")
names(myColors) <- factor(c('RNA', 'iPCR', 'Plasmid', 'DNA'),levels=c('RNA', 'iPCR', 'Plasmid', 'DNA'))
colScale <- scale_fill_manual(name = "Sample", values = myColors)


filename <- paste0(opt$OUTDIR, "/AllBCs.txt")
data <- read.table(filename, header=T, stringsAsFactors=F)


#Normalisation step
normData<-data
normData[,"DNA"] <- normData[,"DNA"]/(sum(data[,"DNA"], na.rm=T)/1e6)
normData[,"RNA"] <- normData[,"RNA"]/(sum(data[,"RNA"], na.rm=T)/1e6)
normData[,"iPCR"] <- normData[,"iPCR"]/(sum(data[,"iPCR"], na.rm=T)/1e6)

#Combine Raw and Normalised counts
#TODO stack or make.group functions to replace this
datamelt <- melt(data)
datamelt[,"Process"] <- c('Raw')
normDatamelt <- melt(normData)
normDatamelt[,"Process"] <- c('Normalised')

combinedData <- rbind(datamelt,normDatamelt)
colnames(combinedData)[1:4] <- c('BC', 'Sample', 'Counts', 'Process')
combinedData<-combinedData[combinedData[,"Counts"]>=1,]

library(ggplot2)
pdf(file=paste0(opt$OUTDIR,"/graphs/Count_variability_normalised.pdf"))
ggplot(combinedData[grep('Normalised', combinedData[,"Process"]), ], aes(x= Sample, y= Counts, fill=Sample)) +
geom_boxplot(color="black", size=0.25, outlier.shape=NA)+ 
facet_wrap(~Process, scale='free')+
theme_classic()+
theme(axis.title.x= element_text(size=14), axis.text.x =element_text(size=12), axis.title.y= element_text(size=14), axis.text.y= element_text(size=12))+
colScale+
ggtitle(opt$SampleName," \n Variability of counts")+
xlab('Samples')+
ylab('Counts')+
coord_cartesian(ylim = c(0, as.numeric(1.5*summary(combinedData[grep('Normalised', combinedData[,4]),][,3])[5])))
#geom_blank(data = dummy,aes(Process,Counts))
dev.off()


filename <- paste0(opt$OUTDIR,"/AllBCs.txt")
data <- read.table(filename, header=T, stringsAsFactors=F)
data <- data[!is.na(data[,"DNA"]),]
data <- data[!is.na(data[,"RNA"]),]
data <- data[!is.na(data[,"iPCR"]),]
normData <- data
normData[,"DNA"]<-normData[,"DNA"]/(sum(data[,"DNA"], na.rm=T)/1e6)
normData[,"RNA"]<-normData[,"RNA"]/(sum(data[,"RNA"], na.rm=T)/1e6)
normData[,"iPCR"]<-normData[,"iPCR"]/(sum(data[,"iPCR"], na.rm=T)/1e6)

normData10 <- normData[normData[,"DNA"]>=1,]
normData10 <- normData10[normData10[,"RNA"]>=1,]

ggRD <- ggplot(normData10, aes(x = DNA, y = RNA)) +
geom_point(color="black", size=0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
colScale+
ggtitle(opt$SampleName," \n Variability of counts --- Over 10 normalised counts")+
xlab('DNA')+
ylab('RNA')+
geom_smooth(method = "lm", se = FALSE)+
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "lb", long = unit(0.2, "cm"))
pdf(file=paste0(opt$OUTDIR,"/graphs/Variability_Normalised_over10.pdf"), width=5, height=5)
ggRD
dev.off()

####
#iPCR specific
####
filename <- paste0(opt$OUTDIR,"/AllInsertions.annotated.bed")
cat('Reading ', filename,'\n')
dataiPCR <- read.table(filename, sep='\t', header=T, stringsAsFactors=F)
dataiPCR <- dataiPCR[!is.na(dataiPCR[,"DNA"]),]
dataiPCR <- dataiPCR[!is.na(dataiPCR[,"RNA"]),]
dataiPCR <- dataiPCR[!is.na(dataiPCR[,"iPCR"]),]


#Normalisation step
normdataiPCR <- dataiPCR
normdataiPCR[,"DNA"] <- normdataiPCR[,"DNA"]/(sum(data[,"DNA"], na.rm=T)/1e6)
normdataiPCR[,"RNA"] <- normdataiPCR[,"RNA"]/(sum(data[,"RNA"], na.rm=T)/1e6)
normdataiPCR[,"iPCR"] <- normdataiPCR[,"iPCR"]/(sum(data[,"iPCR"], na.rm=T)/1e6)
normdataiPCR.unfiltered <- normdataiPCR


#####
#Exploration of counts for DNA and RNA
#####
expExplore <- subset(dataiPCR,select =c(RNA, DNA,BC))
expExplore[,"normRNA"] <- expExplore[,"RNA"]/(sum(data[,"RNA"], na.rm=T)/1e6)
expExplore[,"normDNA"] <- expExplore[,"DNA"]/(sum(data[,"DNA"], na.rm=T)/1e6)
expExplore$expression <- log10(expExplore$normRNA/expExplore$normDNA+1)
cutoffiPCR<- 1
expExplore.cut <- subset(expExplore, expExplore[,"normDNA"] >=cutoffiPCR & expExplore[,"normRNA"] >=cutoffiPCR)

pdf(file=paste0(opt$OUTDIR,"/graphs/DNA_copies_ExpressionEffect.pdf"), width=7, height=5)
boxplot(expExplore.cut$expression,  #all
expExplore.cut[expExplore.cut$normDNA <   quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[2]],]$expression,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[2]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[3]],]$expression,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[3]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[4]],]$expression,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[4]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[5]],]$expression,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[5]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[6]],]$expression,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[6]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[7]],]$expression,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[7]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[8]],]$expression,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[8]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[9]],]$expression,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[9]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[10]],]$expression,
expExplore.cut[expExplore.cut$normDNA > quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[10]],]$expression, #over 3rd
outline=FALSE, notch=T, col=as.character(sampleColor),
names=c('All', c("<10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", ">90%")), las=2, frame=F,
xlab='', ylab='expression (log10(RNA/DNA+1)',
main = paste0(opt$SampleName,'\nFiltering for number of DNA copies \nnormalised expression'))

dev.off()

pdf(file=paste0(opt$OUTDIR,"/graphs/DNA_copies_ExpressionEffect_RawRNA.pdf"), width=7, height=5)
boxplot(expExplore.cut$RNA,  #all
expExplore.cut[expExplore.cut$normDNA <   quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[2]],]$RNA,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[2]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[3]],]$RNA,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[3]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[4]],]$RNA,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[4]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[5]],]$RNA,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[5]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[6]],]$RNA,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[6]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[7]],]$RNA,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[7]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[8]],]$RNA,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[8]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[9]],]$RNA,
expExplore.cut[expExplore.cut$normDNA >  quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[9]] & expExplore.cut$normDNA < quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[10]],]$RNA,
expExplore.cut[expExplore.cut$normDNA > quantile(expExplore.cut$normDNA, prob = seq(0, 1, length = 11), type = 5)[[10]],]$RNA, #over 3rd
outline=FALSE, notch=T, col=as.character(sampleColor),
names=c('All', c("<10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", ">90%")), las=2, frame=F,
xlab='', ylab='Raw RNA counts',
main = paste0(opt$SampleName,'\nFiltering for number of DNA copies\nraw RNA counts'))
dev.off()





GetSeqCount<- function(Sequence) {
str_count(Sequence,c('A','C','G','T'))}

expExplore.cut $A <- 0
expExplore.cut $C <- 0
expExplore.cut $G <- 0
expExplore.cut $T <- 0
for (i in 1:nrow(expExplore.cut )) {
#seqBias[[i]] <- data.frame(BC=expExplore.cut $V1[i], A=GetSeqCount(expExplore.cut $V1[i])[1],
# C=GetSeqCount(expExplore.cut $V1[i])[2],
# G=GetSeqCount(expExplore.cut $V1[i])[3],
# T=GetSeqCount(expExplore.cut $V1[i])[4])
expExplore.cut $A[i] <- GetSeqCount(expExplore.cut $BC[i])[1]
expExplore.cut $C[i] <- GetSeqCount(expExplore.cut $BC[i])[2]
expExplore.cut $G[i] <- GetSeqCount(expExplore.cut $BC[i])[3]
expExplore.cut $T[i] <- GetSeqCount(expExplore.cut $BC[i])[4]
 
}

expExplore.cut$GC_barcode <- (expExplore.cut$C+expExplore.cut$G)/16
expExplore.cut[,"PercentGC.bin_Barcode"] <- cut.pretty(expExplore.cut[,"GC_barcode"], breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")

#pdf(file=paste0(opt$OUTDIR,"/graphs/GC_barcode_ExpressionEffect.pdf"), width=8, height=8)
pdf(file=paste0("~/public_html/blog/2017Aug07/GC_barcode_ExpressionEffet_2",opt$SampleName,".pdf"), width=8, height=8)
ggplot(expExplore.cut, aes(x=log10(normRNA/normDNA), fill=PercentGC.bin_Barcode)) +
geom_histogram(colour='black')+
geom_vline(xintercept = 0,colour='red')+
theme_classic()+
theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12),plot.title=element_text(size=10))+
facet_wrap(~PercentGC.bin_Barcode, scale='free_y') +
#coord_cartesian(xlim = c(0, 4), ylim=c(0,4))+
xlab('log10(normRNA/normDNA)') +
ggtitle( opt$SampleName, '\nGC-bias in expression')+
scale_fill_manual(values=brewer.pal('Blues', n=9)) +
annotation_logticks(sides = "b", long = unit(0.2, "cm"))
dev.off()

#####
#Cutoff formalized reads 
#####
cutoffiPCR<- 1
normdataiPCR <- subset(normdataiPCR, normdataiPCR[,"DNA"] >=cutoffiPCR & normdataiPCR[,"RNA"] >=cutoffiPCR)


normdataiPCR[,"DistToTSS.bin"] <- cut(normdataiPCR[,"DistToTSS"], breaks=c(-10^7, -10^6, -10^5, -10^4, -2500, 0, 5000, 10^4, 10^5, 10^6, 10^7), right=F, include.lowest=T, labels=c("-10 Mb", "-1 Mb", "-100 kb", "-10 kb", "-2.5 kb", "+5 kb", "+10 kb", "+100 kb", "+1 Mb", "+10 Mb"))
normdataiPCR[,"DistToNearestDHS.bin"] <- cut(normdataiPCR[,"DistToNearestDHS"], breaks=c(-10^7, -10^6, -10^5, -10^4, -2500, 0, 5000, 10^4, 10^5, 10^6, 10^7), right=F, include.lowest=T, labels=c("-10 Mb", "-1 Mb", "-100 kb", "-10 kb", "-2.5 kb", "+5 kb", "+10 kb", "+100 kb", "+1 Mb", "+10 Mb"))
normdataiPCR[,"DistToCGI.bin"] <- cut(normdataiPCR[,"DistToCGI"], breaks=c(-10^7, -10^6, -10^5, -10^4, -2500, 0, 5000, 10^4, 10^5, 10^6, 10^7), right=F, include.lowest=T, labels=c("-10 Mb", "-1 Mb", "-100 kb", "-10 kb", "-2.5 kb", "+5 kb", "+10 kb", "+100 kb", "+1 Mb", "+10 Mb"))
normdataiPCR[,"PercentGC.bin"] <- cut.pretty(normdataiPCR[,"PercentGC"], breaks=c(0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")
normdataiPCR[,"HeLaDistToNearestDHS.bin"] <- cut(normdataiPCR[,"HeLaDistToNearestDHS"], breaks=c(-10^7, -10^6, -10^5, -10^4, -2500, 0, 5000, 10^4, 10^5, 10^6, 10^7), right=F, include.lowest=T, labels=c("-10 Mb", "-1 Mb", "-100 kb", "-10 kb", "-2.5 kb", "+5 kb", "+10 kb", "+100 kb", "+1 Mb", "+10 Mb"))
normdataiPCR[,"DistToGencodeGeneTPM.bin"] <- cut(normdataiPCR[,"DistToGencodeGeneTPM"], breaks=c(-10^7, -10^6, -10^5, -10^4, -2500, 0, 5000, 10^4, 10^5, 10^6, 10^7), right=F, include.lowest=T, labels=c("-10 Mb", "-1 Mb", "-100 kb", "-10 kb", "-2.5 kb", "+5 kb", "+10 kb", "+100 kb", "+1 Mb", "+10 Mb"))

dnaDeciles <- quantile(normdataiPCR$DNA, prob = seq(0, 1, length = 11), type = 5)

normdataiPCR$expression <- log10(normdataiPCR$RNA/normdataiPCR$DNA)
normdataiPCR$DNA_Deciles<- cut(normdataiPCR$DNA, 
       breaks=c(dnaDeciles[[1]], dnaDeciles[[2]], dnaDeciles[[3]], dnaDeciles[[4]], dnaDeciles[[5]], dnaDeciles[[6]], dnaDeciles[[7]], dnaDeciles[[8]], dnaDeciles[[9]], dnaDeciles[[10]], dnaDeciles[[11]]),
       right=T, include.lowest=T, labels=c("<10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", ">90%"))





#####
#Create track hub and bigwig of insertions
#####

normUCSC <- subset(normdataiPCR,select=c(chrom, chromStart, chromEnd, BC, strand, DNA, RNA))
#Normalise
normUCSC$Norm <- log10(normUCSC$RNA/normUCSC$DNA+1)
#Create Integrer for UCSC
normUCSC$Integrer <- round(normUCSC$RNA/normUCSC$DNA)
if (nrow(normUCSC[normUCSC$Integrer>=1000,])>0){
       normUCSC[normUCSC$Integrer>=1000,]$Integrer <- 1000
}
#Remove duplicated
normUCSC<-normUCSC[!duplicated(subset(normUCSC,select=c(chrom,chromStart,chromEnd))),]
write.table(subset(normUCSC,select=c(chrom, chromStart, chromEnd, BC, Integrer, strand)),file=paste0(opt$OUTDIR,"/AllInsertions.coords.Integrer.bed"), sep='\t', col.names=F, quote=F, row.names=F)

write.table(subset(normUCSC,select=c(chrom, chromStart, chromEnd, BC, Norm, strand)),file=paste0(opt$OUTDIR,"/AllInsertions.coords.bed"), sep='\t', col.names=F, quote=F, row.names=F)


#
#
#pdf(paste0(opt$OUTDIR,'/graphs/DNA_count_number.pdf'), height=4, width=7)
# normdataiPCR %>%
#       select(DNA_Deciles) %>%
#       count(DNA_Deciles) %>%
#       ggplot(aes(x=DNA_Deciles, y=n, fill=DNA_Deciles)) +
#       #geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=1, colour=sampleColor) +
#       geom_bar(stat='identity', fill=sampleColor, colour='black', alpha=0.6)+ 
#       theme_classic()+
#       theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=60,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
#       ggtitle('CMV - Number of insertions in each bin')+
#       ylab("Number")+
#      # annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#dev.off()
#
#
#cat('DistToNearesDHS for bins','\n')
#pdf(paste0(opt$OUTDIR,'/graphs/DistanceToDHS_Bin_DNAbins.pdf'), height=10, width=14)
#normdataiPCR %>% 
#       mutate(DistToNearestDHS.bin = gsub('-|\\+', '', DistToNearestDHS.bin)) %>%  
#       mutate(DistToNearestDHS.bin = factor(DistToNearestDHS.bin, levels=c('2.5 kb', '5 kb', '10 kb', '100 kb', '1 Mb' ,'10 Mb'))) %>%
#       ggplot(aes(x=DistToNearestDHS.bin, y =expression)) +
#       geom_boxplot(outlier.shape=NA,fill=sampleColor, alpha=0.6)+ 
#       theme_classic()+
#       geom_hline(yintercept=median(normdataiPCR$expression), lty=2, color='blue', size=1.5) +
#       facet_wrap(~DNA_Deciles, ncol=4) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName,"  \ndistToNearestDHS and RNA expression per DNA decile")+
#       xlab('distance to DNase')+
#       #scale_fill_manual(values=c('white', brewer.pal(name='Reds', n=9))) +
#       ylab('log10(RNA/DNA)')+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm"))
#dev.off()
#
#
#
#
#pdf(paste0(opt$OUTDIR,'/graphs/Expression_histogram_DNAbins.pdf'), height=10, width=14)
#ggplot(normdataiPCR, aes( x = expression)) +
#geom_histogram(color="black", binwidth=0.1, fill=sampleColor, alpha=0.6)+ 
#theme_classic()+
#theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#geom_vline(xintercept=median(normdataiPCR$expression), lty=2, color='blue', size=1.5) +
#ggtitle(opt$SampleName,"  \nExpression distribution per DNA decile")+
#facet_wrap(~DNA_Deciles, ncol=4)
#dev.off()
#
#
#
#pdf(paste0(opt$OUTDIR,'/graphs/ExpressionHistogram_RNA_DNAbins.pdf'), height=10, width=14)
#ggplot(normdataiPCR, aes( x = log10(RNA+1))) +
#geom_histogram(color="black", binwidth=0.1, fill=sampleColor, alpha=0.6)+ 
#theme_classic()+
#theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#geom_vline(xintercept=median(log10(normdataiPCR$RNA+1)), lty=2, color='blue', size=1.5) +
#ggtitle(opt$SampleName,"  \nRaw RNA distribution per DNA decile")+
#facet_wrap(~DNA_Deciles, ncol=4) 
#dev.off()
#
#pdf(paste0(opt$OUTDIR,'/graphs/ExpressionHistogram_RNA.pdf'), height=5, width=5)
#ggplot(normdataiPCR, aes( x = log10(RNA+1))) +
#geom_histogram(color="black", binwidth=0.1, fill=sampleColor, alpha=0.6)+ 
#theme_classic()+
#theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#geom_vline(xintercept=median(log10(normdataiPCR$RNA+1)), lty=2, color='blue', size=1.5) +
#ggtitle(opt$SampleName,"  \nRaw RNA distribution")
#dev.off()
#
#pdf(paste0(opt$OUTDIR,'/graphs/ExpressionHistogram_DNA.pdf'), height=5, width=5)
#ggplot(normdataiPCR, aes( x = log10(DNA+1))) +
#geom_histogram(color="black", binwidth=0.1, fill=sampleColor, alpha=0.6)+ 
#theme_classic()+
#theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#geom_vline(xintercept=median(log10(normdataiPCR$RNA+1)), lty=2, color='blue', size=1.5) +
#ggtitle(opt$SampleName,"  \nRaw DNA distribution")
#dev.off()
#
#
#pdf(paste0(opt$OUTDIR,'/graphs/Mappability_20kb.pdf'), height=5, width=5)
#hist(normdataiPCR$mappability20kb,
#n=100, 
#main=paste0(opt$SampleName, '\nmappability 20kb around insertion'), 
#col=sampleColor,
#xlab='Proportion mappability',
#xlim=c(0 ,1))
#dev.off()
#

#####
#Remove <20% DNA deciles, which inflates the expression (optional)
#####
#normdataiPCR <- normdataiPCR[normdataiPCR$DNA_Deciles!="<10%" & normdataiPCR$DNA_Deciles!="10-20%",]



#####
#Histogram of reporter activity distribution
#####
pdf(file=paste0(opt$OUTDIR,"/graphs/ReporterActivityDist.pdf"), width=5, height=5)
ggplot(normdataiPCR, aes(x = RNA/DNA)) +
       geom_histogram(color="black", size=0.25, fill=sampleColor, alpha=0.6)+ 
       theme_classic()+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       #colScale+
       #geom_smooth(method = "lm", se = FALSE)+
       ggtitle(opt$SampleName,"Reporter activity distribution")+
       xlab('log10(RNA/DNA)')+
       #ylab('RNA')+
       scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "b", long = unit(0.2, "cm"))
dev.off()



#####
#Scatterplot between DNA and RNA of insertions within 10kb of promoters
#####
ggRD <- ggplot(subset(normdataiPCR, abs(normdataiPCR[,"DistToTSS"]) <= 1e4), aes(x = DNA, y = RNA)) +
geom_point(color="black", size=0.25,)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
colScale+
geom_smooth(method = "lm", se = FALSE)+
ggtitle(opt$SampleName,"  \n Variability of counts +- 10kb within a promoter")+
xlab('DNA')+
ylab('RNA')+
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "lb", long = unit(0.2, "cm"))
pdf(file=paste0(opt$OUTDIR,"/graphs/Variability_Normalised_with10kb.pdf"), width=5, height=5)
ggRD
dev.off()


cat('DistToTSS','\n')
ggRD <- ggplot(normdataiPCR, aes(x = abs(DistToTSS), y = RNA/DNA)) +
geom_point(color="black", size=0.25,)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \n DistToTSS and RNA expression")+
xlab('distance to TSS')+
geom_smooth(method = "lm", se = FALSE)+
coord_cartesian(xlim = c(0, 1e6))+
ylab('log10(RNA/DNA)')+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))+
scale_x_continuous(label= fancy_scientific)
pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_DistToTSS.pdf"), width=5, height=5)
ggRD
dev.off()


pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_DistToTSS_expressedGenes.bin.pdf"), width=5, height=4)
cat('DistToGencodeGeneTPM>.bin=1','\n')
normdataiPCR %>% 
       mutate(DistToGencodeGeneTPM.bin = gsub('-|\\+', '', DistToGencodeGeneTPM.bin)) %>%  
       mutate(DistToGencodeGeneTPM.bin = factor(DistToGencodeGeneTPM.bin, levels=c('2.5 kb', '5 kb', '10 kb', '100 kb', '1 Mb' ,'10 Mb'))) %>%
       ggplot(aes(x =DistToGencodeGeneTPM.bin, y = RNA/DNA,fill=SampleName)) +
              geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
              theme_classic()+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
              ggtitle(opt$SampleName,"  \n distToNearestExpTSS and RNA expression")+
              xlab('distance to TSS')+
              ylab('log10(RNA/DNA)')+
              scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
              annotation_logticks(sides = "l",long = unit(0.2, "cm"))

dev.off()

cat('DistToGencodeGeneTPM>=1','\n')
ggExpresedGenes <- ggplot(normdataiPCR, aes(x = abs(DistToGencodeGeneTPM+1), y = RNA/DNA)) +
geom_point(color="black", size=0.25,)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName," \n DistToTSS of expressed genes (TPM>=1) and RNA expression")+
xlab('distance to TSS of expressed genes (TPM>=1)')+
coord_cartesian(xlim = c(0, 1e6))+
geom_smooth(method = "lm", se = FALSE)+
ylab('log10(RNA/DNA)')+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l", long = unit(0.2, "cm"))+
scale_x_continuous(label= fancy_scientific)


pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_DistToTSS_expressedGenes.pdf"), width=5, height=5)
ggExpresedGenes
dev.off()



cat('DistToTSSGencodeActive','\n')
ggDistGencode <- ggplot(normdataiPCR, aes(x = abs(DistToGencodeGene)+1, y = RNA/DNA)) +
geom_point( size=0.4,pch=19, colour=sampleColor)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName," \n DistToTSS and RNA expression")+
xlab('distance to TSS')+
geom_smooth(method = "lm")+
geom_smooth(method = "loess", linetype='dashed')+
ylab('log10(RNA/DNA)')+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "lb", long = unit(0.2, "cm"))

ggDistGencodeExpressed <- ggplot(normdataiPCR, aes(x = abs(DistToGencodeGeneTPM)+1, y = RNA/DNA)) +
geom_point( size=0.4,pch=19, colour=sampleColor)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName," \n DistToTSS (TPM>=1) and RNA expression")+
xlab('distance to TSS (TPM>=1)')+
geom_smooth(method = "lm")+
geom_smooth(method = "loess", linetype='dashed')+
ylab('log10(RNA/DNA)')+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "lb", long = unit(0.2, "cm"))

ggDistActivePromoter <- ggplot(normdataiPCR, aes(x = abs(DistToActivePromoter)+1, y = RNA/DNA)) +
geom_point( size=0.4,pch=19, colour=sampleColor)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName," \n DistToDHS-H3K27ac-H3K4me3")+
xlab('DistToDHS-H3K27ac-H3K4me3')+
geom_smooth(method = "lm")+
geom_smooth(method = "loess", linetype='dashed')+
ylab('log10(RNA/DNA)')+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "lb", long = unit(0.2, "cm"))

pdf(file=paste0(opt$OUTDIR,"/graphs/DistToTSS_Gencode_active.pdf"), width=12, height=4)
grid.arrange(ggDistGencode, ggDistGencodeExpressed, ggDistActivePromoter, ncol=3)
dev.off()


cat('DistToNearesDHS','\n')
ggDNaseAll <- ggplot(normdataiPCR, aes(x = abs(DistToNearestDHS)+1, y = RNA/DNA)) +
geom_point( size=0.4,pch=19, colour=sampleColor)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \ndistToNearestDHS and RNA expression")+
xlab('distance to DNase')+
ylab('log10(RNA/DNA)')+
geom_smooth(method = "lm", se = FALSE)+
geom_smooth(method = "loess", linetype='dashed')+
#coord_cartesian(xlim = c(0, 1e6))+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "lb",long = unit(0.2, "cm"))

pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_DistToNearestDHSAll.pdf"), width=5, height=5)
ggDNaseAll
dev.off()




cat('DistToDHS-H3K27ac-H3K4me1','\n')
ggDNaseActive <- ggplot(normdataiPCR, aes(x = abs(DistToActiveEnhancer)+1, y = RNA/DNA)) +
geom_point( size=0.4,pch=19, colour=sampleColor)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \ndistToDHS-H3K27ac-H3K4me1 and RNA expression")+
xlab('distToDHS-H3K27ac-H3K4me1')+
ylab('log10(RNA/DNA)')+
geom_smooth(method = "lm")+
geom_smooth(method = "loess", linetype='dashed')+
#coord_cartesian(xlim = c(0, 1e6))+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "lb",long = unit(0.2, "cm"))

pdf(file=paste0(opt$OUTDIR,"/graphs/DistToNearestActiveEnhancer.pdf"), width=5, height=5)
ggDNaseActive
dev.off()

pdf(file=paste0(opt$OUTDIR,"/graphs/DistToNearestActiveEnhancerVsNearestDHS.pdf"), width=8, height=4)
grid.arrange(ggDNaseAll, ggDNaseActive, ncol=2)
dev.off()


cat('DistToNearesDHS type','\n')
ggDI <- ggplot(normdataiPCR, aes(x = abs(DistToNearestDHS) ,y = RNA/DNA, colour=NearestDHStype)) +
geom_point(size=0.25,pch=21)+ 
theme_classic()+
facet_wrap(~NearestDHStype,ncol=3)+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(paste0(opt$SampleName," \ndistToNearestDHS and RNA expression, CTCF = ",nrow(normdataiPCR[normdataiPCR$NearestDHStype=='CTCF',]), ", Distal = ", nrow(normdataiPCR[normdataiPCR$NearestDHStype=='Distal',]),  ", Promoter = ",nrow(normdataiPCR[normdataiPCR$NearestDHStype=='Promoter',])))+
xlab('distance to DNase')+
ylab('log10(RNA/DNA)')+
coord_cartesian(xlim = c(0, 1e6))+
geom_smooth(method = "lm", se = FALSE)+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))+
scale_x_continuous(label= fancy_scientific)
pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_distToNearestDHS.pdf"),width=15,height=5)
ggDI
dev.off()



cat('DistToNearesDHS strand','\n')
ggDIstrand <- ggplot(normdataiPCR, aes(x = abs(DistToNearestDHS),y = RNA/DNA, colour=NearestDHStype)) +
geom_point(size=0.25,pch=21)+ 
theme_classic()+
facet_wrap(~NearestDHStype+NearestDHSstrand,ncol=9)+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y= element_text(size=14), axis.text.y= element_text(size=12))+
ggtitle(paste0(opt$SampleName,"  \ndistToNearestDHS and RNA expression, CTCF = ",nrow(normdataiPCR[normdataiPCR$NearestDHStype=='CTCF',]), ", Distal = ", nrow(normdataiPCR[normdataiPCR$NearestDHStype=='Distal',]),  ", Promoter = ",nrow(normdataiPCR[normdataiPCR$NearestDHStype=='Promoter',])))+
xlab('distance to DNase')+
ylab('log10(RNA/DNA)')+
coord_cartesian(xlim = c(0, 1e6))+
geom_smooth(method = "lm", se = FALSE)+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))+
scale_x_continuous(label= fancy_scientific)
pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_distToNearestDHS_Strand.pdf"), width=25, height=5)
ggDIstrand
dev.off()


cat('DistToDpn','\n')
ggIR <- ggplot(normdataiPCR, aes(x = abs(DistToDpn),y = RNA/DNA)) +
geom_point(color="black", size=0.25,)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \ndistDpn and RNA expression")+
xlab('distance to Dpn')+
ylab('log10(RNA/DNA)')+
coord_cartesian(xlim = c(0, 1e3))+
geom_smooth(method = "lm", se = FALSE)+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))+
scale_x_continuous(label= fancy_scientific)

pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_distDpn.pdf"), width=5, height=5)
ggIR
dev.off()

cat('DistToCTCF','\n')
ggCTCF <- ggplot(normdataiPCR, aes(x = abs(DistToNearestCTCF),y = RNA/DNA)) +
geom_point(color="black", size=0.25,)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14),axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \ndistToNearestDHS and RNA expression")+
xlab('distance to CTCF')+
ylab('log10(RNA/DNA)')+
coord_cartesian(xlim = c(0, 1e6))+
geom_smooth(method = "lm", se = FALSE)+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))+
scale_x_continuous(label= fancy_scientific)

pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_distCTCF.pdf"), width=5, height=5)
ggCTCF
dev.off()




pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised.pdf"), width=15, height=5)
grid.arrange(ggRD, ggDI, ggIR, ncol=3, nrow =1)
dev.off()



#####
#Correlation between DHS and histone marks
#####
library(corrplot)
histCorr <- subset(normdataiPCR, select=c(DistToNearestDHS, DistToTSS, DistToNearestCTCF, H3K27me3, H3K36me3, H3K4me1, H3K4me2, H3K4me3, H3K79me2, 
H3K9ac, H3K9me1, H3K9me3, H4K20me1, HCFC1, HDAC1, HDAC2, POLR2A, POLR2AphosphoS2, POLR3A))

colnames(histCorr)[1:3] <- c('DHS', 'TSS', 'CTCF')
histCorr <- abs(histCorr)
histCorrelation <- cor(histCorr)


pdf(file=paste0(opt$OUTDIR,"/graphs/HistoneCorrelation.pdf"), width=5, height=5)
corrplot(histCorrelation, type="upper", method="ellipse", tl.pos="d")
#corrplot(histCorrelation, type="lower", method="number", col="black", 
#         add=TRUE, diag=FALSE, tl.pos="n", cl.pos="n")
dev.off()


#Correlation with TFs

TFcorr <- subset(normdataiPCR, select=c(DistToNearestDHS, ARID1B, ARID3A, ARNT, ATF1, ATF2, ATF3, ATF7, 
BACH1, BCLAF1, BCOR,  BHLHE40, BMI1, BRD4,
 C11orf30, CBX1, CBX2, CBX3, CBX5, CBX8, 
CCNT2, CDC5L, CEBPB, CEBPD, CEBPZ, CHAMP1, CHD1, 
CHD4, CHD7, COPS2, CREB1, CREB3L1, CREBBP, CREM, 
CSDE1, CTBP1, CTCF, CTCFL, DDX20, DEAF1, DNMT1, 
DPF2, E2F1, E2F4, E2F6, EGR1, ELF1, ELF4, 
 EP300, ESRRA, ETS1, ETV6,  FOS, FOSL1, 
FOXA1, FOXK2, FOXM1, GABPA, GATA1, GATA2, GFI1B, 
 GTF2F1, GTF3C2, H2AFZ, HDGF, HES1, HLTF, HMBOX1, HMGN3, IKZF1, 
IRF1, IRF2, JUNB, JUN, JUND, KAT2B, KAT8, KDM1A, 
KDM4B, KDM5B, KLF16, KMT2B, L3MBTL2, LARP7, LEF1, 
MAFF, MAFK, MAX, MAZ, MBD2, MCM3, MCM5, MCM7, 
MEF2A, MEIS2, MGA, MIER1, MITF, MLLT1, MNT, MTA2, 
MTA3, MYC, MYNN, NBN, NCOA1, NCOR1, NELFE, NEUROD1, 
NFE2, NFE2L2, NFRKB, NFXL1, NFYA, NFYB, NKRF, NONO, 
NRF1, PHB2, PHF8, PKNOX1, PML,  POU5F1, PRPF4, RAD21, RAD51, RBBP5, 
RCOR1, REST, RFX1, RFX5, RING1, RNF2, RUNX1, SAP30, 
SETDB1, SIN3A, SIN3B, SIRT6, SIX5, SMAD1, SMAD2, 
SMAD5, SMARCA4, SMARCA5, SMARCB1, SMARCE1, SMC3, 
SOX6, SP1, SP2, SPI1, SREBF1, SRF, SRSF3, STAT1, 
STAT2, STAT5A, SUPT5H, TAF1, TAF7, TAL1, 
TARDBP, TBL1XR1, TBP, TCF12, TCF7, TCF7L2, TEAD4, 
TFDP1, THAP1, THRAP3, TOE1, TRIM24, TRIM25, TRIM28, 
UBTF, USF1,  WDR5, WHSC1, YBX1, YBX3, 
YY1, ZBED1, ZBTB11, ZBTB2, ZBTB33, ZBTB40, ZBTB7A, 
ZC3H11A, ZC3H8, ZEB2, ZHX1, ZKSCAN1, ZMIZ1, ZMYM3, 
ZNF24, ZNF263, ZNF274, ZNF316, ZNF318, ZNF384, ZNF407, 
ZSCAN29, ZZZ3))




colnames(TFcorr)[1] <- c('DHS')
TFcorr <- abs(TFcorr)

d <- dist(t(TFcorr), method = "euclidean")
fit <- hclust(d) 
tfCorrelation <- cor(TFcorr)
pdf(file=paste0(opt$OUTDIR,"/graphs/TfCorrelation.pdf"), width=5, height=5)
corrplot(tfCorrelation, type="upper", method="shade", tl.pos="d")
dev.off()

#####
#Surrounding CTCF sites direction and expression
####

normdataiPCR$surroundingCTCF <- paste0(normdataiPCR$CTCFdirection.up,normdataiPCR$CTCFdirection.dn)

surrCTCF <- normdataiPCR %>%
       mutate(surroundingCTCF = paste0(CTCFdirection.up,CTCFdirection.dn)) %>%
       filter(surroundingCTCF %in% c('--','-+','+-','++')) %>%
       mutate(expression = log10(RNA/DNA)) %>%
       ggplot(aes(x=surroundingCTCF, y=expression)) +
       geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
       theme_classic()+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \nSurrounding CTCF durection and RNA expression")+
       xlab('Upstream and downstream CTCF direction')+
       ylab('log10(RNA/DNA)')+
       #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm"))
       
pdf(file=paste0(opt$OUTDIR,"/graphs/surrounding_CTCFdirection.pdf"), width=5, height=5)
surrCTCF
dev.off()

###
#Bins 
#Using Dplyr to mutate bins instead of creating new variables
####

#DHS
cat('DistToDHS.bin','\n')
pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_DistToNearestDHS_Bins.pdf"), width=5, height=4)
normdataiPCR %>% mutate(DistToNearestDHS.bin = gsub('-|\\+', '', DistToNearestDHS.bin)) %>%  
mutate(DistToNearestDHS.bin = factor(DistToNearestDHS.bin, levels=c('2.5 kb', '5 kb', '10 kb', '100 kb', '1 Mb' ,'10 Mb'))) %>%
ggplot(aes(x =DistToNearestDHS.bin, y = RNA/DNA,fill=SampleName)) +
geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \ndistToNearestDHS and RNA expression")+
xlab('distance to K562 DHS')+
ylab('log10(RNA/DNA)')+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))
dev.off()



cat('HeLa DistToDHS.bin','\n')
pdf(file=paste0(opt$OUTDIR,"/graphs/HeLa_DistToNearestDHS_Bins.pdf"), width=5, height=4)
normdataiPCR %>% mutate(HeLaDistToNearestDHS.bin = gsub('-|\\+', '', HeLaDistToNearestDHS.bin)) %>%  
mutate(HeLaDistToNearestDHS.bin = factor(HeLaDistToNearestDHS.bin, levels=c('2.5 kb', '5 kb', '10 kb', '100 kb', '1 Mb' ,'10 Mb'))) %>%
ggplot(aes(x =HeLaDistToNearestDHS.bin, y = RNA/DNA,fill=SampleName)) +
geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nHeLa distToNearestDHS and RNA expression")+
xlab('distance to HeLa DHS')+
ylab('log10(RNA/DNA)')+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))
dev.off()


#TSS
cat('DistToTSS.bin','\n')
pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_DistToTSS_Bins.pdf"), width=5, height=4)
normdataiPCR %>% mutate(DistToTSS.bin = gsub('-|\\+', '', DistToTSS.bin)) %>%  
mutate(DistToTSS.bin = factor(DistToTSS.bin, levels=c('2.5 kb', '5 kb', '10 kb', '100 kb', '1 Mb' ,'10 Mb'))) %>%
ggplot(aes(x =DistToTSS.bin, y = RNA/DNA,fill=SampleName)) +
geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nDistToTSS and RNA expression")+
xlab('distance to TSS')+
ylab('log10(RNA/DNA)')+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))
dev.off()

#CGI
cat('DistToCGI.bin','\n')
pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_DistToCGI_Bins.pdf"), width=5, height=4)
normdataiPCR %>% mutate(DistToCGI.bin = gsub('-|\\+', '', DistToCGI.bin)) %>%  
mutate(DistToCGI.bin = factor(DistToCGI.bin, levels=c('2.5 kb', '5 kb', '10 kb', '100 kb', '1 Mb' ,'10 Mb'))) %>%
ggplot(aes(x =DistToCGI.bin, y = RNA/DNA,fill=SampleName)) +
geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nDistToCGI and RNA expression")+
xlab('distance to CGI')+
ylab('log10(RNA/DNA)')+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))
dev.off()

#Percentage CG within +- 75bp
cat('Percentage CG +- 75bp bin','\n')
pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_PercentGC_Bins.pdf"), width=5, height=4)
normdataiPCR %>%  mutate(PercentGC.bin = factor(PercentGC.bin, levels=c('0+','0.1+','0.2+','0.3+','0.4+','0.5+','0.6+','0.7+','0.8+','0.9+'))) %>%
ggplot(aes(x = PercentGC.bin, y = RNA/DNA,fill=SampleName)) +
geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nPercentGC.bin and RNA expression")+
xlab('Percentage CG +-75bp of insertion')+
ylab('log10(RNA/DNA)')+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))
dev.off()


####
#DNA instead of RNA
#####


#ggRD <- ggplot(normdataiPCR, aes(x = abs(DistToTSS), y = DNA)) +
#geom_point(color="black", size=0.25,)+ 
#theme_classic()+
#theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#ggtitle(opt$SampleName,"  \nDistToTSS and DNA counts")+
#xlab('distance to TSS')+
#ylab('log10(DNA)')+
#coord_cartesian(xlim = c(0, 1e5))+
#scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
#annotation_logticks(sides = "l",long = unit(0.2, "cm"))
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_DistToTSS_DNA.pdf"), width=5, height=5)
#ggRD
#dev.off()
#
#ggDI <- ggplot(normdataiPCR, aes(x = abs(DistToNearestDHS),y = DNA)) +
#geom_point(color="black", size=0.25,)+ 
#theme_classic()+
#theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#ggtitle(opt$SampleName,"  \ndistToNearestDHS and DNA counts")+
#xlab('distance to DNase')+
#ylab('log10(DNA)')+
#coord_cartesian(xlim = c(0, 1e5))+
#scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
#annotation_logticks(sides = "l",long = unit(0.2, "cm"))
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_distToNearestDHS_DNA.pdf"), width=5, height=5)
#ggDI
#dev.off()
#
#ggIR <- ggplot(normdataiPCR, aes(x = abs(DistToDpn), y = DNA)) +
#geom_point(color="black", size=0.25,)+ 
#theme_classic()+
#theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#ggtitle(opt$SampleName," \ndistDpn and DNA counts")+
#xlab('distance to Dpn')+
#ylab('log10(DNA)')+
#coord_cartesian(xlim = c(0, 1e3))+
#scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
#annotation_logticks(sides = "l",long = unit(0.2, "cm"))
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_distDpn_DNA.pdf"), width=5, height=5)
#ggIR
#dev.off()
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/iPCR_Normalised_DNA.pdf"), width=15, height=5)
#grid.arrange(ggRD, ggDI, ggIR, ncol=3, nrow =1)
#dev.off()
##



####
#Number of CTCF, DHS, CGI, TSS, Expressed genes within 100Kb
####
cat('CTCF +-100kb','\n')

ggCTCF <- ggplot(normdataiPCR, aes(x = nCTCFpeaks100kb, y = RNA/DNA)) +
geom_point(color = "black", size = 0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nmean signal and nCTCF")+
geom_smooth(method = "lm", se = FALSE)+
xlab('number of CTCF peaks')+
ylab('mean signal')+
#ylim(c(-2,2))+
coord_cartesian(xlim = c(0, 30))+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))

pdf(file=paste0(opt$OUTDIR,"/graphs/MeanSignal_CTCF.pdf"), width=5, height=5)
ggCTCF
dev.off()

cat('DHS +-100kb','\n')
ggDHS <- ggplot(normdataiPCR, aes(x = nDHS100kb, y = RNA/DNA)) +
geom_point(color = "black", size = 0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nmean signal and nDHS")+
geom_smooth(method = "lm", se = FALSE)+
xlab('number of DHS peaks')+
ylab('mean signal')+
#ylim(c(-2,2))+
coord_cartesian(xlim = c(0, 50))+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))


pdf(file=paste0(opt$OUTDIR,"/graphs/MeanSignal_DHS.pdf"), width=5, height=5)
ggDHS
dev.off()

cat('TSS +-100kb','\n')

ggTSS <- ggplot(normdataiPCR, aes(x = nTSS100kb, y = RNA/DNA)) +
geom_point(color = "black", size = 0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nmean signal and nTSS")+
geom_smooth(method = "lm", se = FALSE)+
xlab('number of gencode TSS')+
ylab('mean signal')+
#ylim(c(-2,2))+
coord_cartesian(xlim = c(0, 30))+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))


pdf(file=paste0(opt$OUTDIR,"/graphs/MeanSignal_TSS.pdf"), width=5, height=5)
ggTSS
dev.off()

cat('Expressed genes +-100kb','\n')
ggExpressed <- ggplot(normdataiPCR, aes(x = nExpressed100kb, y = RNA/DNA)) +
geom_point(color = "black", size = 0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nmean signal and nExpressed genes")+
geom_smooth(method = "lm", se = FALSE)+
xlab('number of expressed gencode')+
ylab('mean signal')+
#ylim(c(-2,2))+
coord_cartesian(xlim = c(0, 20))+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))


pdf(file=paste0(opt$OUTDIR,"/graphs/MeanSignal_Expressed.pdf"), width=5, height=5)
ggExpressed
dev.off()

cat('CGI +-100kb','\n')
ggCGI <- ggplot(normdataiPCR, aes(x = nCGI100kb, y = RNA/DNA)) +
geom_point(color = "black", size = 0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nmean signal and nCGI")+
geom_smooth(method = "lm", se = FALSE)+
xlab('number of CGI')+
ylab('mean signal')+
#ylim(c(-2,2))+
coord_cartesian(xlim = c(0, 40))+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))+



pdf(file=paste0(opt$OUTDIR,"/graphs/MeanSignal_CGI.pdf"), width=5, height=5)
ggCGI
dev.off()

pdf(file=paste0(opt$OUTDIR,"/graphs/100kb_regions.pdf"), width=15, height=3)
grid.arrange(ggDHS, ggCTCF, ggTSS, ggExpressed, ggCGI, ncol=5, nrow =1)
dev.off()



#####
#Analyse A/B compartments 
#####

cat('A/Bcomp','\n')
ABgroups <- normdataiPCR %>%
       mutate(ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
       mutate('expression' = log10(RNA/DNA)) %>%
       filter(!is.na(ABcomp), ABcomp!='NA') 

ABt.test <- wilcox.test(ABgroups[ABgroups$ABcomp=='A',]$expression, ABgroups[ABgroups$ABcomp=='B',]$expression)
ABpval <- ABt.test$p.value



Asize <- fread('~/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/K562_ABcomparments_A.bedGraph', sep='\t', header=F)
Bsize <- fread('~/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/K562_ABcomparments_B.bedGraph', sep='\t', header=F)
ABcompartments <- data.frame(ABcomp=c('A', 'B'), genomicSize=c(sum(Asize$V3-Asize$V2), sum(Bsize$V3-Bsize$V2)))
ABcompartments$genomicSize <- as.numeric(ABcompartments$genomicSize)
ABcompartments$prop <- ABcompartments$genomicSize/sum(ABcompartments$genomicSize)
ABcompartments$Expected <- ABcompartments$prop *nrow(normdataiPCR)
ABcompartments$Observed <- c(nrow(ABgroups[ABgroups$ABcomp=='A',]), nrow(ABgroups[ABgroups$ABcomp=='B',]))

ggABcompObsExp <- ABcompartments %>%
                            mutate(ObsExp=Observed/Expected) %>%
                            ggplot(aes(x = ABcomp, y = ObsExp, fill=ABcomp)) +
                            geom_bar(stat='identity',alpha=0.6, colour='black') +
                            theme_classic()+
                            scale_fill_brewer(palette='Set1') +
                            theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
                            ggtitle(opt$SampleName,"  \nA/B compartments --- Observed/Expected")+
                            geom_hline(yintercept=1, lty=2) +
                            xlab('Genomic compartment')+
                            ylab('Observed/Expected insertions')
                            

ggABcomp <- normdataiPCR %>%
       mutate(ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
       #mutate(expression = log10(RNA/DNA)) %>%
       filter(!is.na(ABcomp), ABcomp!='NA') %>%
       ggplot(aes(x = ABcomp, y = RNA/DNA, fill=ABcomp)) +
       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=21,aes(colour=ABcomp)) +
       geom_boxplot(outlier.shape=NA, alpha=0.6) +
       theme_classic()+
       scale_fill_brewer(palette='Set1') +
       scale_colour_brewer(palette='Set1') +
       annotate('text', label = paste0('p=', format(as.numeric(ABpval),big.mark=",", trim=TRUE, digits=2)), x = 1.5, y = 20,  color = 'black', size=3) +
       annotate('text', label = paste0(format(as.numeric(nrow(ABgroups[ABgroups$ABcomp=='A',])),big.mark=",", trim=TRUE)), x = 1, y = 100, color = 'black') +
       annotate('text', label = paste0(format(as.numeric(nrow(ABgroups[ABgroups$ABcomp=='B',])),big.mark=",", trim=TRUE)), x = 2, y = 100, color = 'black') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \nA/B compartments --- Expression")+
       xlab('Genomic compartment')+
       ylab('log10(RNA/DNA)')+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm"))
       

pdf(file=paste0(opt$OUTDIR,"/graphs/ABcomparment_obs.pdf"), width=7, height=4)
grid.arrange(ggABcomp, ggABcompObsExp, ncol=2)
dev.off()

pdf(file=paste0(opt$OUTDIR,"/graphs/ABcomparment.pdf"), width=4, height=4)
ggABcomp
dev.off()

#####
#Variation of expression in A/B compartments
#####
cat('A/Bcomp s.d./mean','\n')
ABsd <- data.frame('Compartment'=c('A','B'),
       'sd'=c(sd(ABgroups[ABgroups$ABcomp=='A',]$expression),sd(ABgroups[ABgroups$ABcomp=='B',]$expression)),
       'var'=c(var(ABgroups[ABgroups$ABcomp=='A',]$expression),var(ABgroups[ABgroups$ABcomp=='B',]$expression)),
       'mean'=c(mean(ABgroups[ABgroups$ABcomp=='A',]$expression), mean(ABgroups[ABgroups$ABcomp=='B',]$expression)))

ggABsd <- ABsd %>%
              melt() %>%
              ggplot(aes(x = variable, y = value, fill= Compartment)) +
              geom_bar(stat='identity', position = "dodge", colour='black') +
              theme_classic()+
              scale_fill_brewer(palette = 'Set1') +
              geom_hline(yintercept=0) +
              theme(axis.title.x = element_text(size = 14), axis.text.x = element_text(size=12), axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12))+
              ggtitle(opt$SampleName, 'Compartment variability')+
              coord_flip()

pdf(file=paste0(opt$OUTDIR,"/graphs/ABcomparment_variability.pdf"), width=6, height=2.5)
ggABsd
dev.off()




#####
#Distance to DHS in A/B compartments
#####
cat('ABcomp DistToDHS.bin','\n')
pdf(file=paste0(opt$OUTDIR,"/graphs/ABcomp_DistToNearestDHS_Bins.pdf"), width=8, height=4)
normdataiPCR %>% 
       mutate(ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
       filter(!is.na(ABcomp), ABcomp!='NA') %>%
       mutate(DistToNearestDHS.bin = gsub('-|\\+', '', DistToNearestDHS.bin)) %>%  
       mutate(DistToNearestDHS.bin = factor(DistToNearestDHS.bin, levels=c('2.5 kb', '5 kb', '10 kb', '100 kb', '1 Mb' ,'10 Mb'))) %>%
       ggplot(aes(x =DistToNearestDHS.bin, y = RNA/DNA,fill=ABcomp)) +
       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=19,aes(group=ABcomp), colour=sampleColor) +
       geom_boxplot(outlier.shape=NA,aes(fill=ABcomp),alpha=0.6)+ 
       theme_classic()+
       #facet_wrap(~ABcomp, ncol=2) +
       scale_fill_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \ndistToNearestDHS and RNA expression in A/B compartments")+
       xlab('distance to DNase')+
       ylab('log10(RNA/DNA)')+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm"))
dev.off()



#####
cat('ABcomp DistToDHS','\n')
pdf(file=paste0(opt$OUTDIR,"/graphs/ABcomp_DistToNearestDHS.pdf"), width=6, height=4)
normdataiPCR %>% 
       mutate(ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
       filter(!is.na(ABcomp), ABcomp!='NA') %>%
       #mutate(DistToNearestDHS.bin = gsub('-|\\+', '', DistToNearestDHS.bin)) %>%  
       #mutate(DistToNearestDHS.bin = factor(DistToNearestDHS.bin, levels=c('2.5 kb', '5 kb', '10 kb', '100 kb', '1 Mb' ,'10 Mb'))) %>%
       ggplot(aes(x =log10(abs(DistToNearestDHS)+1), y = RNA/DNA, fill=ABcomp, colour=ABcomp)) +
       geom_point(alpha=0.2, size=0.5, pch=21, aes(fill=ABcomp)) +
#       geom_boxplot(outlier.shape=NA,aes(fill=ABcomp),alpha=0.6)+ 
       theme_classic()+
       geom_smooth(method='loess')+
       #facet_wrap(~ABcomp, ncol=2) +
       scale_fill_brewer(palette='Set1') +
       scale_colour_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \ndistToNearestDHS and RNA expression in A/B compartments")+
       xlab('Distance to DHS')+
       ylab('log10(RNA/DNA)')+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "lb",long = unit(0.2, "cm"))
dev.off()


cat('ABcomp DistToDHS.bin by cell cycle','\n')
pdf(file=paste0(opt$OUTDIR,"/graphs/ABcomp_cellCycle_DistToNearestDHS_Bins.pdf"), width=20, height=4)
normdataiPCR %>% 
       mutate(ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
       filter(!is.na(ABcomp), ABcomp!='NA',
       !is.na(repliSeq), repliSeq!='NA') %>%
       mutate(DistToNearestDHS.bin = gsub('-|\\+', '', DistToNearestDHS.bin)) %>%  
       mutate(DistToNearestDHS.bin = factor(DistToNearestDHS.bin, levels=c('2.5 kb', '5 kb', '10 kb', '100 kb', '1 Mb' ,'10 Mb'))) %>%
       mutate(repliSeq = factor(repliSeq, levels=c('G1', 'S1', 'S2', 'S3', 'S4', 'G2'))) %>%
       ggplot(aes(x =DistToNearestDHS.bin, y = RNA/DNA,fill=ABcomp)) +
       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=19,aes(group=ABcomp), colour=sampleColor) +
       geom_boxplot(outlier.shape=NA,aes(fill=ABcomp),alpha=0.6)+ 
       theme_classic()+
       geom_hline(yintercept=1, linetype='dashed')+
       facet_wrap(~repliSeq, ncol=6) +
       scale_fill_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \ndistToNearestDHS and RNA expression in A/B compartments")+
       xlab('distance to DNase')+
       ylab('log10(RNA/DNA)')+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm"))
dev.off()



#####
#ChromHMM
#####
cat('chromHMM','\n')
ggchromHMM <- normdataiPCR %>%
       mutate(ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
       mutate(chromHMM = factor(chromHMM, levels=unique(normdataiPCR$chromHMM)[order(as.numeric(gsub('_.*','',unique(normdataiPCR$chromHMM))))])) %>%
       mutate(log10DNARNA = log10(RNA/DNA)) %>% 
       filter(!is.na(chromHMM), !is.na(ABcomp), ABcomp!='NA') %>%
       group_by(chromHMM,SampleName, ABcomp) %>%
       count() %>%
       ggplot(aes(x=chromHMM,y=n,fill=ABcomp)) +
       geom_bar(stat='identity',alpha=1,color='black')+ 
       theme_classic()+
       scale_fill_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=60,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
       ggtitle(opt$SampleName,' Insertions per chromHMM18 state --- number')+
       ylab("Insertions")

#cat('plot ChromHMMbar','\n')
#pdf(file=paste0(opt$OUTDIR,"/graphs/chromHMM_bar.pdf"), width=5, height=4)
#ggchromHMM
#dev.off()
ggchromHMMExpression <- normdataiPCR %>%
       mutate(ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
       mutate(log10DNARNA = log10(RNA/DNA)) %>% 
       filter(!is.na(chromHMM), !is.na(ABcomp), ABcomp!='NA') %>%
       mutate(chromHMM = factor(chromHMM, levels=unique(normdataiPCR$chromHMM)[order(as.numeric(gsub('_.*','',unique(normdataiPCR$chromHMM))))])) %>%
       ggplot(aes(x=chromHMM,y=log10DNARNA,fill=SampleName)) +
       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=1, colour=sampleColor) +
       geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
       theme_classic()+
       theme(axis.title.x = element_text(size=14),axis.text.x  = element_text(size=12,angle=60,hjust=1),axis.title.y=element_text(size=14),axis.text.y=element_text(size=12))+
       ggtitle(opt$SampleName,'Insertions per chromHMM18 state --- expression')+
       ylab("log10(RNA/DNA)")+
       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
       theme(legend.position="none")+
       coord_cartesian(ylim = c(-3, 3))

#cat('plot ChromHMMbox','\n')
#pdf(file=paste0(opt$OUTDIR,"/graphs/chromHMM_box.pdf"), width=5, height=4)
#ggchromHMM
#dev.off()

pdf(file=paste0(opt$OUTDIR,"/graphs/chromHMM.pdf"), width=10, height=4)
grid.arrange(ggchromHMM, ggchromHMMExpression, ncol=2)
dev.off()


######
#geneDesert, geneRegion, heteroChromatin
######
cat('Genomic division','\n')
ggGenomicDivision <- normdataiPCR %>%
       mutate(ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
       #mutate(expression = log10(RNA/DNA)) %>%
       filter(!is.na(genomicDivision), genomicDivision!='NA') %>%
       filter(!is.na(ABcomp), ABcomp!='NA') %>%
       filter(!is.na(genomicDivision), 
       genomicDivision!='NA',
       genomicDivision!="heteroChromatin") %>%
       mutate(genomicDivision = factor(genomicDivision, levels=c('geneRich', 'geneDesert'))) %>%
       ggplot(aes(x = genomicDivision, y = RNA/DNA, fill=genomicDivision, colour=genomicDivision)) +
       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=19, aes(colour=genomicDivision)) +
       geom_boxplot(outlier.shape=NA, alpha=0.6, aes(fill=genomicDivision), colour='black') +
       theme_classic()+
       #facet_wrap(~genomicDivision)+
       scale_fill_manual(values=c('#386cb0', '#fdc086')) +
       scale_colour_manual(values=c('#386cb0', '#fdc086')) +
       annotate('text', label = format(as.numeric(nrow(normdataiPCR[normdataiPCR$genomicDivision=='geneRich',])),big.mark=",", trim=TRUE), x = 1, y = 100, color = 'black') +
       annotate('text', label = format(as.numeric(nrow(normdataiPCR[normdataiPCR$genomicDivision=='geneDesert',])),big.mark=",", trim=TRUE), x = 2, y = 100, color = 'black') +
       annotate('text', label = paste0('p=',prettyNum(t.test(normdataiPCR[normdataiPCR$genomicDivision=='geneRich',]$expression, normdataiPCR[normdataiPCR$genomicDivision=='geneDesert',]$expression)$p.value, digits=2)), x = 1.5, y = 20, color = 'black', size=3) +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \nGenomic region")+
       xlab('Genomic compartment')+
       ylab('log10(RNA/DNA)')+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm"))
       

pdf(file=paste0(opt$OUTDIR,"/graphs/GenomicDivision.pdf"), width=4, height=4)
ggGenomicDivision
dev.off()



cat('Genomic division DistToDHS.bin','\n')
pdf(file=paste0(opt$OUTDIR,"/graphs/GenomicDivision_DistToNearestDHS_Bins.pdf"), width=8, height=4)
normdataiPCR %>% 
       mutate(ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
       filter(!is.na(ABcomp), ABcomp!='NA') %>%
       mutate(DistToNearestDHS.bin = gsub('-|\\+', '', DistToNearestDHS.bin)) %>%  
       mutate(DistToNearestDHS.bin = replace(DistToNearestDHS.bin , DistToNearestDHS.bin == '2.5 kb', '5 kb')) %>% 
       mutate(DistToNearestDHS.bin = factor(DistToNearestDHS.bin, levels=c('2.5 kb', '5 kb', '10 kb', '100 kb', '1 Mb' ,'10 Mb'))) %>%
       filter(!is.na(genomicDivision), 
       genomicDivision!='NA',
       genomicDivision!="heteroChromatin") %>%
       ggplot(aes(x =DistToNearestDHS.bin, y = RNA/DNA, group=DistToNearestDHS.bin)) +
       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=19, aes(fill=SampleName), colour=sampleColor) +
       geom_boxplot(outlier.shape=NA, fill=sampleColor, alpha=0.6)+ 
       theme_classic()+
       facet_wrap(~genomicDivision, ncol=3) +
       #scale_fill_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \ndistToNearestDHS and RNA expression in A/B compartments by genomics division")+
       xlab('distance to DNase')+
       ylab('log10(RNA/DNA)')+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm"))+
       theme(legend.position="none")
dev.off()

cat('Genomic division DistToDHS','\n')
ggDesertDHS <- normdataiPCR %>% 
       filter(!is.na(genomicDivision), 
       genomicDivision!='NA',
       genomicDivision!="heteroChromatin") %>%
       mutate(genomicDivision = factor(genomicDivision, levels=c('geneRich', 'geneDesert'))) %>%
       ggplot(aes(x =log10(abs(DistToNearestDHS)+1), y = RNA/DNA, fill=genomicDivision, colour=genomicDivision)) +
       geom_point(alpha=0.2, pch=21, size=0.4, aes(fill=genomicDivision)) +
       #geom_boxplot(outlier.shape=NA, fill=sampleColor, alpha=0.6)+ 
       theme_classic()+
       #facet_wrap(~genomicDivision, ncol=3) +
       geom_smooth(method ='loess') +
       scale_fill_manual(values=c('#386cb0', '#fdc086')) +
       scale_colour_manual(values=c('#386cb0', '#fdc086')) +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \ndistToNearestDHS and RNA expression by genomics division")+
       xlab('distance to DHS')+
       ylab('log10(RNA/DNA)')+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "lb",long = unit(0.2, "cm"))
       
pdf(file=paste0(opt$OUTDIR,"/graphs/GenomicDivision_DistToNearestDHS.pdf"), width=6, height=4)
ggDesertDHS
dev.off()


cat('ABcomp genomic division DistToDHS.bin per cell division bin','\n')
pdf(file=paste0(opt$OUTDIR,"/graphs/ABcomp_GenomicDivision_cellCycle_DistToNearestDHS_Bins.pdf"), width=20, height=12)
normdataiPCR %>% 
       mutate(ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
       filter(!is.na(ABcomp), ABcomp!='NA',
       !is.na(repliSeq), repliSeq!='NA',
       !is.na(genomicDivision), genomicDivision!='NA') %>%
       mutate(DistToNearestDHS.bin = gsub('-|\\+', '', DistToNearestDHS.bin)) %>%  
       mutate(DistToNearestDHS.bin = factor(DistToNearestDHS.bin, levels=c('2.5 kb', '5 kb', '10 kb', '100 kb', '1 Mb' ,'10 Mb'))) %>%
       mutate(repliSeq = factor(repliSeq, levels=c('G1', 'S1', 'S2', 'S3', 'S4', 'G2'))) %>%
       ggplot(aes(x =DistToNearestDHS.bin, y = RNA/DNA,fill=ABcomp)) +
       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=19,aes(group=ABcomp), colour=sampleColor) +
       geom_boxplot(outlier.shape=NA,aes(fill=ABcomp),alpha=0.6)+ 
       theme_classic()+
       geom_hline(yintercept=1, linetype='dashed')+
       facet_wrap(~genomicDivision+repliSeq, ncol=6) +
       scale_fill_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \ndistToNearestDHS and RNA expression in A/B compartments by genomics division")+
       xlab('distance to DNase')+
       ylab('log10(RNA/DNA)')+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm"))
dev.off()

#####
#Insertions in geneDeserts and their distance to the end
#####
cat('distance to gene desert border', '\n')
ggGeneDesertDist <- normdataiPCR %>%
       mutate(expression = log10(RNA/DNA),
       DistToGeneDesert=abs(DistToGeneDesert)) %>%
       filter(genomicDivision=='geneDesert') %>%
       mutate(geneDesertdistance.bin =cut(DistToGeneDesert, breaks=c(0, 5000, 10^4, 5*10^4, 10^5, 2*10^5, 5*10^5, 10^6, 10^7), right=F, include.lowest=F, labels=c("5 kb", "10 kb", "50 kb","100 kb", "200 kb", "500 kb", "1 Mb", "+1 Mb"))) %>%       
       ggplot(aes(x = geneDesertdistance.bin, y = expression, fill=opt$SampleName)) +
       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=19, colour=sampleColor) +
       geom_boxplot(outlier.shape=NA, width=.2, alpha=0.6, colour='black') +
       geom_violin( alpha=0.6, fill=sampleColor, colour='black') +
       theme_classic()+
       #annotate('text', label = paste0('Pvalue = ', format(as.numeric(ABpval),big.mark=",", trim=TRUE, digits=2)), x = 1.5, y = 3.5, color = 'black') +
       #annotate('text', label = paste0('nA = ', format(as.numeric(nrow(ABgroups[ABgroups$ABcomp=='A',])),big.mark=",", trim=TRUE)), x = 1, y = 3, color = 'black') +
       #annotate('text', label = paste0('nB = ', format(as.numeric(nrow(ABgroups[ABgroups$ABcomp=='B',])),big.mark=",", trim=TRUE)), x = 2, y = 3, color = 'black') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \nInsertions in gene deserts and how distance\nto the border influences expression")+
       xlab('Distance to geneDesert border')+
       ylab('log10(RNA/DNA)') +
       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
       theme(legend.position="none")

ggsave(file=paste0(opt$OUTDIR,"/graphs/GeneDesertDisttoBorder.pdf"), ggGeneDesertDist, width=6, height=4)


#####
#K562 gene desert distance to border
#####
ggk562GeneDesertDist <- normdataiPCR %>%
       mutate(expression = log10(RNA/DNA),
       DistTok562GeneDesert=abs(DistTok562GeneDesert)) %>%
       filter(K562geneDesert=='K562geneDesert') %>%
       mutate(K562geneDesertdistance.bin =cut(DistTok562GeneDesert, breaks=c(0, 5000, 10^4, 5*10^4, 10^5, 2*10^5, 5*10^5, 10^6, 10^7), right=F, include.lowest=F, labels=c("5 kb", "10 kb", "50 kb","100 kb", "200 kb", "500 kb", "1 Mb", "+1 Mb"))) %>%          
       ggplot(aes(x = K562geneDesertdistance.bin, y = expression, fill=opt$SampleName)) +
       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=19, colour=sampleColor) +
       geom_boxplot(outlier.shape=NA, width=.2, alpha=0.6, colour='black') +
       geom_violin( alpha=0.6, fill=sampleColor, colour='black') +
       theme_classic()+
       #annotate('text', label = paste0('Pvalue = ', format(as.numeric(ABpval),big.mark=",", trim=TRUE, digits=2)), x = 1.5, y = 3.5, color = 'black') +
       #annotate('text', label = paste0('nA = ', format(as.numeric(nrow(ABgroups[ABgroups$ABcomp=='A',])),big.mark=",", trim=TRUE)), x = 1, y = 3, color = 'black') +
       #annotate('text', label = paste0('nB = ', format(as.numeric(nrow(ABgroups[ABgroups$ABcomp=='B',])),big.mark=",", trim=TRUE)), x = 2, y = 3, color = 'black') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \nInsertions in K562 gene deserts and how distance\nto the border influences expression")+
       xlab('Distance to K562 geneDesert border')+
       ylab('log10(RNA/DNA)') +
       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
       theme(legend.position="none")

ggsave(file=paste0(opt$OUTDIR,"/graphs/GeneDesert_K562_DisttoBorder.pdf"), ggk562GeneDesertDist, width=6, height=4)

#
# test <- normdataiPCR %>%
#              mutate(expression = log10(RNA/DNA),
#              DistTok562GeneDesert=abs(DistTok562GeneDesert)) %>%
#              filter(K562geneDesert=='K562geneDesert') %>%
#              mutate(K562geneDesertdistance.bin =cut(DistTok562GeneDesert, breaks=c(0, 5000, 10^4, 5*10^4, 10^5, 2*10^5, 5*10^5, 10^6, 10^7), right=F, include.lowest=T, labels=c("+5 kb", "+10 kb", "+50 kb","+100 kb", "+200 kb", "+500 kb", "+1 Mb", "+10 Mb")))  %>% 
#              filter(K562geneDesert=="K562geneDesert",
#                     K562geneDesertdistance.bin =="+500 kb", 
#                     AB_compartment_value<0, 
#                     DNA_Deciles!='<10%', 
#                     DNA_Deciles!='10-20%')
#                     
#test <- test[rev(order(test$expression)),]
#
#              
######
#TF analysis
#######
cat('TF linear correlation','\n')
#TF <- read.table(paste0(opt$OUTDIR,"/",opt$SampleName,".TFbinding.tsv"), sep='\t', header=T, stringsAsFactors=F)
tfMarks <- readLines(paste0(opt$OUTDIR,"/TFmarks"))
TF <- normdataiPCR
TF$expression <- log10(TF$RNA/TF$DNA)
TF <- subset(TF,select=c('expression',tfMarks,'DistToNearestDHS'))

for (i in 2:ncol(TF)){
TF[,i] <- abs(TF[,i])

}
#TF$Expressed <- TF$expression>0

#TF <- subset(TF, select=-c(chrom,chromStart, chromEnd, BC,  strand, expression))

fit <- lm(expression ~ ., data=TF)
TFresults <- as.data.frame(summary(fit)$coefficient) # show results
TFresults$FDR <- p.adjust(TFresults[,"Pr(>|t|)"], method ="fdr", n = nrow(TFresults))

sigTFs <- rownames(TFresults[TFresults$FDR<0.05,])[-1]


f_plot_ROC <- function (outcome, risk, main="ROC Curve") {
  # This function plots the ROC Curve for a risk model
  # Example call
  # risk <- deathModelAll.glm$fitted.values
  # outcome <-deathModelAll.glm$y
  # f_plot_ROC(outcome,risk)
  # Version 20120707
  
  if (length(outcome) != length(risk)) stop ("outcome and risk must be same length")
  if (min(risk)<0.0 || max(risk) > 1.0 ) stop ("Not valid risk")
  if (length(table(outcome)) !=2  ) stop ("Not valid outcome")
  
  library(ROCR)
  preds <- prediction(as.numeric(risk),
                      as.numeric(outcome))
  perf <- performance(preds, "auc")
  rocC <- round(as.numeric(slot(perf,"y.values")),3)
  perf_plot <- performance(preds, "tpr", "fpr")
  plot(perf_plot, main=main, lwd=2, col="darkgreen")
  text(0.5, 0.5, paste("AUC = ",rocC), cex=2.2, col="darkred")
  return(rocC)
}




if (length(sigTFs)>1) {
       
       #pdf(paste0(opt$OUTDIR,"/graphs/TFsglm_ROC.pdf"))
       #f_plot_ROC(fit$y,fit$fitted.values) 
       #dev.off()
       
       #TF <- read.table(paste0(opt$OUTDIR,"/",opt$SampleName,".TFbinding.tsv"), sep='\t', header=T, stringsAsFactors=F)
       #TF <- subset(TF, select=-c(chrom,chromStart, chromEnd, BC,  strand, H3K4me3))

       sigTFs <- TF[,colnames(TF) %in% c(sigTFs,'expression')]
       
       
       
       ggTF <- sigTFs %>%
              melt(id=1) %>%
              mutate(DistanceBin= cut(value,  breaks=c(-10^7, -10^6, -10^5, -10^4, -2500, 0, 5000, 10^4, 10^5, 10^6, 10^7), right=F, include.lowest=T, labels=c("-10 Mb", "-1 Mb", "-100 kb", "-10 kb", "-2.5 kb", "+5 kb", "+10 kb", "+100 kb", "+1 Mb", "+10 Mb"))) %>%
              mutate(DistanceBin = gsub('-|\\+', '', DistanceBin)) %>%  
              mutate(DistanceBin = factor(DistanceBin, levels=c('2.5 kb', '5 kb', '10 kb', '100 kb', '1 Mb' ,'10 Mb'))) %>%
              filter(DistanceBin != "NA" ) %>%
              ggplot(aes(x = DistanceBin ,y =expression)) +
              geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
              theme_classic()+
              facet_wrap(~variable,ncol=7)+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
              ggtitle(paste0(opt$SampleName,"\nDistance to TFs and expression\nLinear regression FDR < 0.05"))+
              xlab('distance to TF mark')+
              ylab('log10(RNA/DNA)') +
              annotation_logticks(sides = "l",long = unit(0.2, "cm"))
       
       ggsave(paste0(opt$OUTDIR,"/graphs/TFslinearCorr.pdf"), ggTF, width=10)
       
      cat('Joy plot','\n')
       joyTFs <- subset(normdataiPCR,select=c(colnames(sigTFs)[-1],'AB_compartment_value','genomicDivision','DNA','RNA'))
       joyTFs$expression <- log10(joyTFs$RNA/joyTFs$DNA)
       joyTFs <- subset(joyTFs, select=-c(DNA, RNA))
       joyOrder <- subset(joyTFs,select=-c(AB_compartment_value, genomicDivision))
       joyOrder <- names(colMeans(joyOrder[joyOrder$expression>0,])[-length(colMeans(joyOrder[joyOrder$expression>0,]))][rev(order(colMeans(joyOrder[joyOrder$expression>0,])[-ncol(joyOrder)]))])

       ggsigTFsJoy <- joyTFs %>%
                     mutate(expression = ifelse(expression >0, 'Expressed', ifelse(expression <0, 'Nonexpressed',0))) %>%
                     mutate(AB_compartment_value = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
                     melt(id=c('expression','AB_compartment_value','genomicDivision')) %>%
                     filter(!is.na(genomicDivision), genomicDivision!='NA') %>%
                     filter(!is.na(AB_compartment_value), AB_compartment_value!='NA') %>%
                     mutate(variable = factor(variable, levels=joyOrder)) %>%
                     ggplot(aes(x=log10(abs(value)+1), y=variable, fill=expression, height=..density..)) +
                     geom_joy(scale=2,alpha=0.8) +
                    # facet_wrap(~AB_compartment_value+genomicDivision) +
                     scale_fill_brewer(palette = 'Set1') +
                     scale_y_discrete(expand=c(0.01, 0)) +
                     scale_x_continuous(expand=c(0, 0)) + theme_joy()+
                     ggtitle(opt$SampleName,"  \nDistance to TFs")+
                     xlab('log10(Distance +1)')+
                     ylab('Transcription factors')
       #ggsave(paste0(opt$OUTDIR,"/graphs/TFs_DistJoyTEST.pdf"), ggsigTFsJoy, width=15, height=20)

              ggsave(paste0(opt$OUTDIR,"/graphs/TFs_DistJoy.pdf"), ggsigTFsJoy, width=10, height=10)

       #####
       #Random forrest
       #####
       #cat('TF Random forrest','\n')
       #
       #
       #sigTFs$class <- TF$expression>0
       #
       #sigTFs <- subset(sigTFs,select=-expression)
       #
       #set.seed(415)
       #fit <- randomForest(as.factor(class) ~ . ,
       #                      data=sigTFs, 
       #                      importance=TRUE, 
       #                      ntree=2000)
       #
       #rfTFs <- melt(importance(fit)[,3:4])
       #
       #TForder<- rfTFs$X1[order(rfTFs[rfTFs$X2=='MeanDecreaseAccuracy',]$value)]
       #ggRF <- rfTFs %>%
       #              mutate(X2 = factor(X2, levels=c('MeanDecreaseAccuracy','MeanDecreaseGini'))) %>%
       #              mutate(X1 = factor(X1, levels=TForder)) %>%
       #              ggplot(aes(x = X1 ,y =value)) +
       #              geom_point(pch=21, fill=sampleColor,alpha=0.6, size=3)+ 
       #              theme_classic()+
       #              facet_wrap(~X2,ncol=2,scale='free_x')+
       #              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       #              ggtitle(opt$SampleName,paste0("\nRandom Forest analysis of significant TFs/Histone marks FDR < 0.05\nOOB error rate = ",round(tail(fit$err.rate[,1],1)*100,2),'%'))+
       #              xlab('Transcription factors/ Histones Marks')+
       #              ylab('Value') +
       #              coord_flip()
       #
       #ggsave(paste0(opt$OUTDIR,"/graphs/TF_RandomForest.pdf"), ggRF, width=7, height=8)

}



######
#Random forest for features
######
#####
#cat('Random forest of features','\n')
sumExpression <- summary(log10(normdataiPCR$RNA/normdataiPCR$DNA))
#
###DnaseContact <- read.table(paste0(opt$OUTDIR, "/DnaseContact.txt"),sep=' ', stringsAsFactors=F)
##pdf('~/public_html/blog/2017Jul17/CMV_exp.pdf')
##boxplot(log10(normdataiPCR$RNA/normdataiPCR$DNA), xlab='CMV', ylab='Expression',cex.axis=1.5, cex.lab=1.5, cex.main=1.5, main='CMV expression')
##abline(h=as.numeric(sumExpression[5]), lty=2, col='red', lwd=3)
##abline(h=as.numeric(sumExpression[2]), lty=2, col='blue', lwd=3)
##dev.off()
##
##pdf('~/public_html/blog/2017Jul17/CMV_exp_dist.pdf')
##hist(log10(normdataiPCR$RNA/normdataiPCR$DNA), n=50, xlab='Expression', ylab='Distribution',cex.axis=1.5, cex.lab=1.5,cex.main=1.5,  main='CMV expression distribution')
##abline(v=as.numeric(sumExpression[5]), lty=2, col='red',lwd = 3)
##abline(v=as.numeric(sumExpression[2]), lty=2, col='blue', lwd=3)
##dev.off()
#
#
#RFfeature <- subset(normdataiPCR,select=c(BC, DNA, RNA, DistToTSS, DistToNearestDHS, DistToNearestCTCF, DistToGencodeGeneTPM, DistToCGI, PercentGC, AB_compartment_value,
# chromHMM, GenicLocation.simple, nDHS100kb, nDHS10kb, nTSS100kb, repliSeq, nExpressed100kb, nCGI100kb, nCTCFpeaks100kb, DNasePksDens, NearestDHStype, CTCFPksDens, genomicDivision))
#  
#RFfeature.clean <- RFfeature %>%
#       filter(!is.na(repliSeq)) %>%
#       #mutate(expression = ifelse(log10(RNA/DNA) >, 'expressed', ifelse(log10(RNA/DNA) <0, 'nonexpressed',0)), 
#       mutate(Expression = ifelse(log10(RNA/DNA) >as.numeric(sumExpression[5]), 1, ifelse(log10(RNA/DNA) <(as.numeric(sumExpression[2])), 0,NA)), 
#       ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', NA)), 
#       chromHMM = factor(chromHMM),
#       repliSeq =factor(repliSeq),
#       #CGIovr = factor(CGIovr), 
#       NearestDHStype = factor(NearestDHStype),
#       genomicDivision = factor(genomicDivision), 
#       #DistToNearestDHS.up = abs(DistToNearestDHS.up), 
#       #DistToNearestDHS.dn = abs(DistToNearestDHS.dn),
#       #DistToNearestCTCF.up = abs(DistToNearestCTCF.up), 
#       #DistToNearestCTCF.dn = abs(DistToNearestCTCF.dn),
#       geneLoc = factor(GenicLocation.simple),
#       DistToTSS = abs(DistToTSS),
#       DistToNearestDHS = abs(DistToNearestDHS),
#       DistToCGI = abs(DistToCGI),
#       DistToNearestCTCF = abs(DistToNearestCTCF),
#       DistToGencodeGeneTPM = log10(abs(DistToGencodeGeneTPM)+1)) %>%
#       filter(!is.na(Expression)) %>%
#       select(DistToTSS, DistToNearestDHS, DistToNearestCTCF , DistToGencodeGeneTPM, DistToCGI, repliSeq, PercentGC, chromHMM, AB_compartment_value, geneLoc, Expression, nDHS10kb, nDHS100kb, nTSS100kb, 
#       nExpressed100kb, nCGI100kb, nCTCFpeaks100kb, DNasePksDens, NearestDHStype, CTCFPksDens, genomicDivision)
#       
#
#RFfeature.clean  <- subset(RFfeature.clean ,!is.na(AB_compartment_value) & !is.na(chromHMM) & !is.na(genomicDivision))
#
#RFfeature.clean$Expression <- as.factor(RFfeature.clean$Expression)
#set.seed(415)
#fit <- randomForest(Expression ~ . ,
#                      data=RFfeature.clean, 
#                      importance=TRUE, 
#                      ntree=2000)
#
#
##x <- ctree(expression ~ ., data=RFfeature.clean)
##plot(x, type="simple")
#
#
##pt <- party:::prettytree(cf@ensemble[[1]], names(cf@data@get("input"))) 
##pt 
##nt <- new("BinaryTree") 
##nt@tree <- pt 
##nt@data <- cf@data 
##nt@responses <- cf@responses 
##nt 
##plot(nt) 
#
#
#rfTFs <- melt(importance(fit)[,3:4])
#
#TForder<- rfTFs$X1[order(rfTFs[rfTFs$X2=='MeanDecreaseAccuracy',]$value)]
#ggRF <- rfTFs %>%
#              mutate(X2 = factor(X2, levels=c('MeanDecreaseAccuracy','MeanDecreaseGini'))) %>%
#              mutate(X1 = factor(X1, levels=TForder)) %>%
#              ggplot(aes(x = X1 ,y =value)) +
#              geom_point(pch=21, fill=sampleColor,alpha=0.6, size=3)+ 
#              theme_classic()+
#              facet_wrap(~X2,ncol=2,scale='free_x')+
#              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#              ggtitle(opt$SampleName,paste0("\nRandom Forest analysis of selected features\nOOB error rate = ",round(tail(fit$err.rate[,1],1)*100,2),'%'))+
#              xlab('Features')+
#              ylab('Value') +
#              coord_flip()
#
#ggsave(paste0(opt$OUTDIR,"/graphs/Features_RandomForest.pdf"), ggRF, width=6, height=8)
#
#cat('ROC plots of features','\n')
#featureROC.clean <- RFfeature %>%
#       filter(!is.na(repliSeq)) %>%
#       mutate(expression = ifelse(log10(RNA/DNA) >as.numeric(sumExpression[5]), 1, ifelse(log10(RNA/DNA) <(as.numeric(sumExpression[2])), 0,NA)), 
#       #mutate(expression = log10(RNA/DNA), 
#       #ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 0)), 
#       #chromHMM = as.numeric(gsub('_.*','',chromHMM)),
#       #geneLoc = as.numeric(as.factor(GenicLocation.simple)),
#       DistToTSS = abs(DistToTSS),
#       chromHMM = factor(chromHMM),
#       geneLoc = factor(GenicLocation.simple),
#       NearestDHStype = factor(NearestDHStype), 
#       repliSeq = factor(repliSeq),
#       genomicDivision = factor(genomicDivision),
#       DistToNearestDHS = abs(DistToNearestDHS),
#       DistToCGI = abs(DistToCGI),
#       DistToNearestCTCF = abs(DistToNearestCTCF),
#       DistToGencodeGeneTPM = log10(DistToGencodeGeneTPM+1)) %>%
#       select(DistToTSS, DistToNearestDHS, DistToNearestCTCF , DistToGencodeGeneTPM, DistToCGI, PercentGC, chromHMM, AB_compartment_value, geneLoc, expression, nDHS10kb, nDHS100kb, nTSS100kb, 
#       nExpressed100kb, nCGI100kb, nCTCFpeaks100kb, repliSeq, DNasePksDens, NearestDHStype, CTCFPksDens, genomicDivision)
#
#
#featureROC.clean <- featureROC.clean[!is.na(featureROC.clean$expression),]
#featureROC.clean <- subset(featureROC.clean, select=-c(DistToCGI, DistToGencodeGeneTPM, geneLoc, nCGI100kb, nExpressed100kb, nCTCFpeaks100kb, chromHMM, DistToTSS, CTCFPksDens, DNasePksDens, nTSS100kb))
#Expression <- featureROC.clean$expression
#featureROC.clean <- subset(featureROC.clean, select=-c(expression))
#
#
##featurePerm <- featureROC.clean[,sample(ncol(featureROC.clean))]
#featurePerm <-  featureROC.clean[,c("nDHS10kb", "nDHS100kb", "DistToNearestDHS", "AB_compartment_value", "NearestDHStype", "repliSeq", "DistToNearestCTCF", "PercentGC", "genomicDivision")]
#pdf(paste0(opt$OUTDIR,"/graphs/",opt$SampleName,"_ROC.pdf"), height=10, width=9)
#par(mfrow=c(3,4))
#for (i in length(colnames(featurePerm)):2) {cat(i, '\n')
#       featureSubset <- featurePerm[,c(1:i)]
#       
#       fit <- glm(formula=Expression ~. ,data=featureSubset, family='binomial')
#       
#       f_plot_ROC(fit$y,fit$fitted.values, main="")
#       text(x=0.5, y=0.2, labels= paste0(colnames(featureSubset),collapse='\n'))
#}
#fit <- glm(formula=Expression ~DistToNearestDHS ,data=featureROC.clean , family='binomial')
#f_plot_ROC(fit$y,fit$fitted.values, main="")
#text(x=0.5, y=0.2, labels= "DistToNearestDHS")
#
#fit <- glm(formula=Expression ~AB_compartment_value ,data=featureROC.clean , family='binomial')
#f_plot_ROC(fit$y,fit$fitted.values, main="")
#text(x=0.5, y=0.2, labels= "AB_compartment_value")
#
#fit <- glm(formula=Expression ~nDHS100kb ,data=featureROC.clean , family='binomial')
#f_plot_ROC(fit$y,fit$fitted.values, main="")
#text(x=0.5, y=0.2, labels= "nDHS100kb")
#
#fit <- glm(formula=Expression ~nDHS10kb ,data=featureROC.clean , family='binomial')
#f_plot_ROC(fit$y,fit$fitted.values, main="")
#text(x=0.5, y=0.2, labels= "nDHS10kb")
#
#
#title(paste0(opt$SampleName, "\nFeature to differentiate between Expressed (log10(RNA/DNA) > ", as.numeric(sumExpression[5]),") and Nonexpressed (log10(RNA/DNA) <- ",as.numeric(sumExpression[2]),")"),line=-3, cex=2, outer = TRUE)
#dev.off()
#

######
#geneDesert analysis
######
cat('Analyzing genomic context','\n')
genomicContext <- normdataiPCR
genomicContext$expression <- log10(genomicContext$RNA/genomicContext$DNA)

ggGenomicContext <- genomicContext %>%
                     ggplot(aes(x=genomicDivision, y=expression, fill=genomicDivision)) +
                     geom_boxplot(outlier.shape=NA) +
                     theme_classic()+
                     scale_fill_manual(values=c('#fdc086','#386cb0','#f0027f'))+
                     theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
                     annotate('text', label = paste0('n = ', format(nrow(genomicContext[genomicContext$genomicDivision=='geneDesert',]),big.mark=",", trim=TRUE, digits=2)), x = 1, y = 2, color = 'black') +
                     annotate('text', label = paste0('n = ', format(nrow(genomicContext[genomicContext$genomicDivision=='geneRich',]),big.mark=",", trim=TRUE)), x = 2, y = 2, color = 'black') +
                     annotate('text', label = paste0('n = ', format(nrow(genomicContext[genomicContext$genomicDivision=='heteroChromatin',]),big.mark=",", trim=TRUE)), x = 3, y = 2, color = 'black') +
                     ggtitle(opt$SampleName,paste0("Expression of insertions dependent on genomic context"))

pdf(paste0(opt$OUTDIR,"/graphs/GenomicContext.pdf"), width=5, height=4)
ggGenomicContext
dev.off()


ggGenomicContextNoExp <- genomicContext %>%
                     ggplot(aes(x=genomicDivisionNoExp, y=expression, fill=genomicDivisionNoExp)) +
                     geom_boxplot(outlier.shape=NA) +
                     theme_classic()+
                     scale_fill_manual(values=c('#fdc086','#386cb0'))+
                     theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
                     annotate('text', label = paste0('n = ', format(nrow(genomicContext[genomicContext$genomicDivisionNoExp=='geneDesert',]),big.mark=",", trim=TRUE, digits=2)), x = 1, y = 2, color = 'black') +
                     annotate('text', label = paste0('n = ', format(nrow(genomicContext[genomicContext$genomicDivisionNoExp=='geneRich',]),big.mark=",", trim=TRUE)), x = 2, y = 2, color = 'black') +
                     ggtitle(opt$SampleName,paste0("Expression of insertions dependent on genomic context"))

pdf(paste0(opt$OUTDIR,"/graphs/GenomicContextNoExp.pdf"), width=5, height=4)
ggGenomicContextNoExp
dev.off()





geneDesert <- subset(genomicContext, genomicDivision=='geneDesert')

cat('geneDesert','\n')
ggDesertJoy <- geneDesert %>%
       mutate(expression = ifelse(expression >0.3, 'Expressed', ifelse(expression <(-0.3), 'Nonexpressed',NA))) %>%
       select(POLR2AphosphoS2, POLR2A, DistToTSS, DistToNearestDHS, DistToNearestCTCF, DistToCGI, expression, AB_compartment_value) %>%
       mutate(AB_compartment_value = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
       melt(id=c('expression','AB_compartment_value')) %>%
       #filter(!is.na(genomicDivision), genomicDivision!='NA') %>%
       filter(!is.na(AB_compartment_value), AB_compartment_value!='NA', !is.na(expression)) %>%
      # mutate(variable = factor(variable, levels=joyOrder)) %>%
       ggplot(aes(x=log10(abs(value)+1), y=variable, fill=expression, height=..density..)) +
       geom_joy(scale=2, alpha=0.8) +
       facet_wrap(~AB_compartment_value) +
       scale_fill_brewer(palette = 'Set1') +
       scale_y_discrete(expand=c(0.01, 0)) +
       scale_x_continuous(expand=c(0, 0)) + 
       theme_joy() +
       ggtitle(opt$SampleName,"  \nDistance to features")+
       xlab('log10(Distance +1)')+
       ylab('Features')

pdf(paste0(opt$OUTDIR,"/graphs/DesertJoy.pdf"), width=8, height=8)
ggDesertJoy
dev.off()



  
geneDesert.clean <- geneDesert %>%
       #mutate(expression = ifelse(log10(RNA/DNA) >, 'expressed', ifelse(log10(RNA/DNA) <0, 'nonexpressed',0)), 
       mutate(expression = ifelse(expression >as.numeric(sumExpression[5]), 1, ifelse(expression <(as.numeric(sumExpression[2])), 0,NA)), 
       ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', NA)), 
       chromHMM = factor(chromHMM),
       CGIovr = factor(CGIovr), 
       NearestDHStype = factor(NearestDHStype),
       genomicDivision = factor(genomicDivision), 
       DistToNearestDHS.up = abs(DistToNearestDHS.up), 
       DistToNearestDHS.dn = abs(DistToNearestDHS.dn),
       DistToNearestCTCF.up = abs(DistToNearestCTCF.up), 
       DistToNearestCTCF.dn = abs(DistToNearestCTCF.dn),
       geneLoc = factor(GenicLocation.simple),
       DistToTSS = abs(DistToTSS),
       DistToNearestDHS = abs(DistToNearestDHS),
       DistToCGI = abs(DistToCGI),
       DistToNearestCTCF = abs(DistToNearestCTCF),
       DistToGencodeGeneTPM = log10(DistToGencodeGeneTPM+1)) %>%
       filter(!is.na(expression)) %>%
       select(DistToTSS, DistToNearestDHS, DistToNearestCTCF , DistToGencodeGeneTPM, DistToCGI, PercentGC, chromHMM, AB_compartment_value, geneLoc, expression, nDHS100kb, nTSS100kb, 
       nExpressed100kb, nCGI100kb, nCTCFpeaks100kb, DNasePksDens, NearestDHStype, CTCFPksDens, genomicDivision)
       

geneDesert.clean  <- subset(geneDesert.clean ,!is.na(AB_compartment_value) & !is.na(chromHMM))

#set.seed(415)
#fit <- randomForest(as.factor(expression) ~ . ,
#                      data=geneDesert.clean, 
#                      importance=TRUE, 
#                      ntree=2000)
#
dir.create(paste0(opt$OUTDIR,'/Sodaplot'))
if (nrow( subset(geneDesert,geneDesert$AB_compartment_value<0 & geneDesert$expression > as.numeric(sumExpression[5])))>0) {
       geneDesExp <- subset(geneDesert,geneDesert$AB_compartment_value<0 & geneDesert$expression > as.numeric(sumExpression[5]))
       write.table(geneDesExp[,c(1:6)], file=paste0(opt$OUTDIR,'/Sodaplot/geneDeserts.bed'), sep='\t', col.names=F, row.names=F, quote=F)
}



#####
#Replication phase
#####
cat('Replication phases' ,'\n')
repliSeqN <- normdataiPCR %>%
       filter(!is.na(repliSeq), repliSeq!='NA') %>%
       group_by(repliSeq) %>%
       count()

repSeq <- fread('~/public_html/blog/OverAllAnalysis/K562_geneExp_ABcomp/repliSeq/repliPhases.bed', sep='\t', header=F)

ggreplication <- repSeq %>% 
              mutate(length=V3-V2) %>% 
              select(length, V4) %>% 
              group_by( V4) %>%   
              summarise(sum = sum(length)/1e6) %>% 
              mutate(repliSeq=factor(V4, levels=c('G1','S1', 'S2','S3','S4','G2'))) %>%
              mutate(Expected = sum/ sum(sum)*nrow(normdataiPCR)) %>%
              left_join(repliSeqN, by='repliSeq') %>%
              mutate(ObsExp = n/Expected) %>%
              mutate(repliSeq=factor(V4, levels=c('G1','S1', 'S2','S3','S4','G2'))) %>%
              ggplot(aes(x=repliSeq, y=ObsExp)) + 
              geom_bar(stat='identity',fill=sampleColor,alpha=0.6, colour='black') + 
              theme_classic() +
              geom_hline(yintercept=1, lty=2, size=1)+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
              ggtitle(opt$SampleName,"  \nObserved vs expected insertion per replication phase")+
              xlab('Replication phase') + 
              ylab('Observed/Expected')

ggrepliSeq <- normdataiPCR %>%
       mutate(expression = log10(RNA/DNA)) %>%
       filter(!is.na(repliSeq), repliSeq!='NA') %>%
       mutate(repliSeq = factor(repliSeq, levels=c('G1', 'S1','S2','S3','S4','G2'))) %>%
       ggplot(aes(x = repliSeq, y = RNA/DNA, fill=SampleName)) +
       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=1, colour=sampleColor) +
       geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
       annotate('text', label = format(as.numeric(repliSeqN[repliSeqN$repliSeq=='G1',]$n),big.mark=",", trim=TRUE), x = 1, y = 100, color = 'black') +
       annotate('text', label = format(as.numeric(repliSeqN[repliSeqN$repliSeq=='S1',]$n),big.mark=",", trim=TRUE), x = 2, y = 100, color = 'black') +
       annotate('text', label = format(as.numeric(repliSeqN[repliSeqN$repliSeq=='S2',]$n),big.mark=",", trim=TRUE), x = 3, y = 100, color = 'black') +
       annotate('text', label = format(as.numeric(repliSeqN[repliSeqN$repliSeq=='S3',]$n),big.mark=",", trim=TRUE), x = 4, y = 100, color = 'black') +
       annotate('text', label = format(as.numeric(repliSeqN[repliSeqN$repliSeq=='S4',]$n),big.mark=",", trim=TRUE), x = 5, y = 100, color = 'black') +
       annotate('text', label = format(as.numeric(repliSeqN[repliSeqN$repliSeq=='G2',]$n),big.mark=",", trim=TRUE), x = 6, y = 100, color = 'black') +
       annotate('text', label = paste0('p=',prettyNum(t.test(normdataiPCR[normdataiPCR$repliSeq=='G1',]$expression, normdataiPCR[normdataiPCR$repliSeq=='S1',]$expression)$p.value, digits=2)), x = 1.5, y = 20, color = 'black', size=3) +
       annotate('text', label = paste0('p=',prettyNum(t.test(normdataiPCR[normdataiPCR$repliSeq=='S1',]$expression, normdataiPCR[normdataiPCR$repliSeq=='S2',]$expression)$p.value, digits=2)), x = 2.5, y = 20, color = 'black', size=3) +
       annotate('text', label = paste0('p=',prettyNum(t.test(normdataiPCR[normdataiPCR$repliSeq=='S2',]$expression, normdataiPCR[normdataiPCR$repliSeq=='S3',]$expression)$p.value, digits=2)), x = 3.5, y = 20, color = 'black', size=3) +
       annotate('text', label = paste0('p=',prettyNum(t.test(normdataiPCR[normdataiPCR$repliSeq=='S3',]$expression, normdataiPCR[normdataiPCR$repliSeq=='S4',]$expression)$p.value, digits=2)), x = 4.5, y = 20, color = 'black', size=3) +
       annotate('text', label = paste0('p=',prettyNum(t.test(normdataiPCR[normdataiPCR$repliSeq=='S4',]$expression, normdataiPCR[normdataiPCR$repliSeq=='G2',]$expression)$p.value, digits=2)), x = 5.5, y = 20, color = 'black', size=3) +
       theme_classic()+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \nReplication timing")+
       xlab('Replication phase')+
       ylab('log10(RNA/DNA)') +
       theme(legend.position="none")+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm"))
       
       
pdf(paste0(opt$OUTDIR,"/graphs/repliSeq.pdf"), width=10, height=4)
grid.arrange(ggrepliSeq, ggreplication, ncol=2)
dev.off()


ggrepliSeqColour <- normdataiPCR %>%
       mutate(expression = log10(RNA/DNA)) %>%
       filter(!is.na(repliSeq), repliSeq!='NA') %>%
       mutate(repliSeq = factor(repliSeq, levels=c('G1', 'S1','S2','S3','S4','G2'))) %>%
       ggplot(aes(x = repliSeq, y = RNA/DNA, fill=repliSeq, colour=repliSeq)) +
       geom_point(position=position_jitterdodge(dodge.width=0.75),aes(fill=repliSeq), alpha=0.2, pch=21) +
       geom_boxplot(outlier.shape=NA,alpha=0.6, colour='black')+ 
       annotate('text', label = format(as.numeric(repliSeqN[repliSeqN$repliSeq=='G1',]$n),big.mark=",", trim=TRUE), x = 1, y = 100, color = 'black') +
       annotate('text', label = format(as.numeric(repliSeqN[repliSeqN$repliSeq=='S1',]$n),big.mark=",", trim=TRUE), x = 2, y = 100, color = 'black') +
       annotate('text', label = format(as.numeric(repliSeqN[repliSeqN$repliSeq=='S2',]$n),big.mark=",", trim=TRUE), x = 3, y = 100, color = 'black') +
       annotate('text', label = format(as.numeric(repliSeqN[repliSeqN$repliSeq=='S3',]$n),big.mark=",", trim=TRUE), x = 4, y = 100, color = 'black') +
       annotate('text', label = format(as.numeric(repliSeqN[repliSeqN$repliSeq=='S4',]$n),big.mark=",", trim=TRUE), x = 5, y = 100, color = 'black') +
       annotate('text', label = format(as.numeric(repliSeqN[repliSeqN$repliSeq=='G2',]$n),big.mark=",", trim=TRUE), x = 6, y = 100, color = 'black') +
       annotate('text', label = paste0('p=',prettyNum(t.test(normdataiPCR[normdataiPCR$repliSeq=='G1',]$expression, normdataiPCR[normdataiPCR$repliSeq=='S1',]$expression)$p.value, digits=2)), x = 1.5, y = 20, color = 'black', size=3) +
       annotate('text', label = paste0('p=',prettyNum(t.test(normdataiPCR[normdataiPCR$repliSeq=='S1',]$expression, normdataiPCR[normdataiPCR$repliSeq=='S2',]$expression)$p.value, digits=2)), x = 2.5, y = 20, color = 'black', size=3) +
       annotate('text', label = paste0('p=',prettyNum(t.test(normdataiPCR[normdataiPCR$repliSeq=='S2',]$expression, normdataiPCR[normdataiPCR$repliSeq=='S3',]$expression)$p.value, digits=2)), x = 3.5, y = 20, color = 'black', size=3) +
       annotate('text', label = paste0('p=',prettyNum(t.test(normdataiPCR[normdataiPCR$repliSeq=='S3',]$expression, normdataiPCR[normdataiPCR$repliSeq=='S4',]$expression)$p.value, digits=2)), x = 4.5, y = 20, color = 'black', size=3) +
       annotate('text', label = paste0('p=',prettyNum(t.test(normdataiPCR[normdataiPCR$repliSeq=='S4',]$expression, normdataiPCR[normdataiPCR$repliSeq=='G2',]$expression)$p.value, digits=2)), x = 5.5, y = 20, color = 'black', size=3) +
       theme_classic()+
       scale_fill_brewer(palette='Set1') +
       scale_colour_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \nReplication timing")+
       xlab('Replication phase')+
       ylab('log10(RNA/DNA)') +
       theme(legend.position="none")+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm"))

       
pdf(paste0(opt$OUTDIR,"/graphs/repliSeq_colour.pdf"), width=6, height=4)
ggrepliSeqColour
dev.off()

#Repliseq and distance to DHS
cat('RepliSeq and distance to DHS', '\n')
ggrepliSeqDistDHS <- normdataiPCR %>%
       #mutate(expression = log10(RNA/DNA)) %>%
       filter(!is.na(repliSeq), repliSeq!='NA') %>%
       mutate(repliSeq = factor(repliSeq, levels=c('G1', 'S1','S2','S3','S4','G2'))) %>%
       ggplot(aes(x = log10(abs(DistToNearestDHS)+1), y = RNA/DNA, fill=repliSeq, colour=repliSeq)) +
       geom_point(alpha=0.2, size=0.4, pch=21, aes(fill=repliSeq)) +
       #geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
       theme_classic() +
       geom_smooth(method='loess') +
       scale_fill_brewer(palette='Set1') +
       scale_colour_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \nReplication timing and distance to DHS")+
       xlab('Distance to DHS')+
       ylab('log10(RNA/DNA)')+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "lb",long = unit(0.2, "cm"))

pdf(paste0(opt$OUTDIR,"/graphs/repliSeq_DistanceToDHS.pdf"), width=6, height=4)
ggrepliSeqDistDHS
dev.off()



####
#Dnase within 100kb
####
cat('number of DNase around +-100 kb', '\n')
ggDnaseContacts <- normdataiPCR %>%

       select(BC, DNA, RNA, SampleName, nDHS100kb) %>%
       mutate(bins = cut.pretty(nDHS100kb,breaks=c(0, 1, 2, 5, 10, 20, 50,100, 10000), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
       mutate(Expression = log10(RNA/DNA)) %>%
       ggplot(aes(x=bins, y=Expression, fill=SampleName)) +
       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=1, colour=sampleColor) +
       geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
       theme_classic()+
       #scale_fill_manual(values=sampleColor) +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \nNumber of DNase peaks surrounding\n +-100kb of insertion")+
       xlab('Number of DNase peaks (+-100kb)')+
       ylab('log10(RNA/DNA)')+
       coord_cartesian(ylim = c(-3, 3))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
       theme(legend.position="none")


pdf(file=paste0(opt$OUTDIR,"/graphs/DNasearound_insertion.pdf"), width=4, height=4)
ggDnaseContacts
dev.off()





ggNearestDHSdensity <- normdataiPCR %>%
       filter(!is.na(DNasePksDens)) %>%
       ggplot(aes(x=DNasePksDens, y=RNA/DNA, fill=SampleName)) +
       geom_point( pch=21, colour=sampleColor) +
       #geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
       geom_smooth(method = "lm") +
       theme_classic()+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName, 'Nearest DHS density')+
       xlab('Nearest DHS density')+
       ylab('log10(RNA/DNA)')+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "lb",long = unit(0.2, "cm"))


ggnearestDHS <- normdataiPCR %>%
       filter(!is.na(DistToNearestDHS)) %>%
       ggplot(aes(x=abs(DistToNearestDHS), y=RNA/DNA, fill=SampleName)) +
       geom_point( pch=21, colour=sampleColor) +
       #geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
       geom_smooth(method = "loess") +
       theme_classic()+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName, 'Distance to nearest DHS')+
       xlab('Distance to DHS')+
       ylab('log10(RNA/DNA)')+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "lb",long = unit(0.2, "cm"))


gg300kbDHS <- normdataiPCR %>%
       filter(!is.na(nDHS300kb)) %>%
       ggplot(aes(x=abs(nDHS300kb), y=RNA/DNA, fill=SampleName)) +
       geom_point( pch=21, colour=sampleColor) +
       #geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
       geom_smooth(method = "lm") +
       theme_classic()+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName, 'number of DHS +-300kb from insertion')+
       xlab('number of DHS +-300kb')+
       ylab('log10(RNA/DNA)')+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm"))


ggDistToNearestDHS <- normdataiPCR %>% mutate(DistToNearestDHS.bin = gsub('-|\\+', '', DistToNearestDHS.bin)) %>%  
mutate(DistToNearestDHS.bin = factor(DistToNearestDHS.bin, levels=c('2.5 kb', '5 kb', '10 kb', '100 kb', '1 Mb' ,'10 Mb'))) %>%
ggplot(aes(x =DistToNearestDHS.bin, y = RNA/DNA,fill=SampleName)) +
geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \ndistToNearestDHS and RNA expression")+
xlab('distance to K562 DHS')+
ylab('log10(RNA/DNA)')+
scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "l",long = unit(0.2, "cm"))


linearGenome <- data.frame(variable=c('DistToNearestDHS', 'nDHS10kb', 'nDHS100kb', 'nDHS300kb', 'DNasePksDens'),
       pval=c(signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~log10(abs(normdataiPCR$DistToNearestDHS)+1)))$coef[2,4],2), 
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$nDHS10kb)))$coef[2,4],2),
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$nDHS100kb)))$coef[2,4],2),
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$nDHS300kb)))$coef[2,4],2),
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~log10(abs(normdataiPCR$DNasePksDens)+1)))$coef[2,4],2)),
       slope=c(signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~log10(abs(normdataiPCR$DistToNearestDHS)+1)))$coef[[2]],2), 
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$nDHS10kb)))$coef[[2]],2),
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$nDHS100kb)))$coef[[2]],2),
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$nDHS300kb)))$coef[[2]],2),
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~log10(abs(normdataiPCR$DNasePksDens)+1)))$coef[[2]],2)),
       adjR2=c(signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~log10(abs(normdataiPCR$DistToNearestDHS)+1)))$adj.r.squared,2), 
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$nDHS10kb)))$adj.r.squared,2),
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$nDHS100kb)))$adj.r.squared,2),
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$nDHS300kb)))$adj.r.squared,2),
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~log10(abs(normdataiPCR$DNasePksDens)+1)))$adj.r.squared,2)),
       Intercept=c(signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~log10(abs(normdataiPCR$DistToNearestDHS)+1)))$coef[[1]],2), 
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$nDHS10kb)))$coef[[1]],2),
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$nDHS100kb)))$coef[[1]],2),
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$nDHS300kb)))$coef[[1]],2),
              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~log10(abs(normdataiPCR$DNasePksDens)+1)))$coef[[1]],2)),
              x=c(2, 2.5, 10, 30, 2))
              
ggLinearGenome <- normdataiPCR %>% 
       mutate(expression = log10(RNA/DNA),
                     DistToNearestDHS = log10(abs(DistToNearestDHS)+1),
                     TAD_armatus0.5_length = log10(TAD_armatus0.5_length),
                     DNasePksDens = log10(DNasePksDens)) %>%
              filter(!is.na(expression)) %>%
              select(expression, DistToNearestDHS, nDHS10kb, nDHS100kb, nDHS300kb, DNasePksDens) %>%
              melt(id='expression') %>%
       ggplot(aes(x=value, y=expression, colour=variable)) +
       geom_point( pch=21, colour=sampleColor) +
       facet_wrap(~variable, scale='free_x') +
       geom_smooth(method = "lm") +
       geom_smooth(method = "loess", linetype='dashed') +
       geom_text(data=linearGenome, aes(y=2, x=x, label=paste0('adjR2=',adjR2)), colour='black', size=5)+
       geom_text(data=linearGenome, aes(y=1.7, x=x, label=paste0('slope=',slope)), colour='black', size=5)+
       geom_text(data=linearGenome, aes(y=1.4, x=x, label=paste0('intercept=',Intercept)), colour='black', size=5)+
       theme_classic()+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName, 'Linear genome variables')+
       ylab('log10(RNA/DNA)')+
       annotation_logticks(sides = "l",long = unit(0.2, "cm"))


pdf(file=paste0(opt$OUTDIR,"/graphs/LinearGenome.pdf"), width=14, height=8)
ggLinearGenome
dev.off()



#####
#TAD analysis
#####

cat('TAD analysis - DHS densityin 0.2', '\n')
ggTADs0.2DHS <- normdataiPCR %>%
       mutate(Expression = log10(RNA/DNA)) %>%
       ggplot(aes(x=TAD_armatus0.2_DHS/TAD_armatus0.2_length, y=Expression, fill=SampleName)) +
       geom_point(alpha=0.6, pch=1, colour=sampleColor) +
       theme_classic()+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"TADs (armatus 0.2) and number of DHS correlated to expression")+
       xlab('DHS density in each armatus 0.2 TAD')+
       geom_smooth(method = "loess") +
       ylab('log10(RNA/DNA)')+
       coord_cartesian(ylim = c(-3, 3))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
       theme(legend.position="none")

pdf(file=paste0(opt$OUTDIR,"/graphs/TAD_0.2_expression.pdf"), width=6, height=4)
ggTADs0.2DHS 
dev.off()

cat('TAD analysis - DHS density in 0.5', '\n')
ggTADs0.5DHS <- normdataiPCR %>%
       mutate(Expression = log10(RNA/DNA)) %>%
       ggplot(aes(x=TAD_armatus0.5_DHS/TAD_armatus0.5_length, y=Expression, fill=SampleName)) +
       geom_point(alpha=0.6, pch=1, colour=sampleColor) +
       theme_classic()+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"TADs (armatus 0.5) and number of DHS correlated to expression")+
       xlab('DHS density in each armatus 0.5 TAD')+
       geom_smooth(method = "loess") +
       ylab('log10(RNA/DNA)')+
       coord_cartesian(ylim = c(-3, 3))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
       theme(legend.position="none")

pdf(file=paste0(opt$OUTDIR,"/graphs/TAD_0.5_expression.pdf"), width=6, height=4)
ggTADs0.5DHS 
dev.off()


cat('TAD analysis - Distance to DHS in TAD', '\n')
ggTADs0.5DHS <- normdataiPCR %>%
       mutate(Expression = log10(RNA/DNA)) %>%
       filter(!is.na(TAD_armatus0.5)) %>%
       mutate(DHStad = ifelse(TAD_armatus0.5 == TAD_nearestDHS, 'same', ifelse(TAD_armatus0.5 != TAD_nearestDHS , 'different', NA))) %>%
       mutate(DHStad = factor(DHStad, levels=c('same', 'different'))) %>%
       ggplot(aes(x=log10(abs(DistToNearestDHS)+1), y=Expression, fill=DHStad, colour=DHStad)) +
       geom_point(alpha=0.2, size=0.4, pch=21, aes(fill=DHStad)) +
       theme_classic()+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"Distance to DHS and if nearest DHS\n is within the same TAD as insertion or not ")+
       xlab('Distance to DHS')+
       scale_fill_brewer(palette='Set1') +
       scale_colour_brewer(palette='Set1') +
       geom_smooth(method = "loess") +
       ylab('log10(RNA/DNA)')+
       coord_cartesian(ylim = c(-3, 3))+
       annotation_logticks(sides = "lb",long = unit(0.2, "cm")) 

pdf(file=paste0(opt$OUTDIR,"/graphs/TAD_0.5_nearestDHS_expression.pdf"), width=6, height=4)
ggTADs0.5DHS 
dev.off()




cat('TAD analysis - DHS density in 0.8', '\n')
ggTADs0.8DHS <- normdataiPCR %>%
       mutate(Expression = log10(RNA/DNA)) %>%
       ggplot(aes(x=TAD_armatus0.8_DHS/TAD_armatus0.8_length, y=Expression, fill=SampleName)) +
       geom_point(alpha=0.6, pch=1, colour=sampleColor) +
       theme_classic()+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"TADs (armatus 0.8) and number of DHS correlated to expression")+
       xlab('DHS density in each armatus 0.8 TAD')+
       geom_smooth(method = "loess") +
       ylab('log10(RNA/DNA)')+
       coord_cartesian(ylim = c(-3, 3))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
       theme(legend.position="none")

pdf(file=paste0(opt$OUTDIR,"/graphs/TAD_0.8_expression.pdf"), width=6, height=4)
ggTADs0.8DHS 
dev.off()





######
#ROC plot from all features minus contacts
######
cat('ROC plots of features','\n')

featureROC <- subset(normdataiPCR,select=c(BC, DNA, RNA, DistToTSS, DistToNearestDHS, DistToNearestCTCF, DistToCGI, PercentGC, AB_compartment_value,
 chromHMM, GenicLocation.simple, nDHS100kb, nDHS10kb, nDHS300kb, nDHS500kb, nDHS1Mb, nDHS100kbDensity, nDHS300kbsum, nDHS500kbDensity, nDHS1MbDensity, nTSS100kb, 
 repliSeq, nExpressed100kb, nCGI100kb, nCTCFpeaks100kb, DNasePksDens, NearestDHStype, CTCFPksDens, genomicDivision, DistToGeneDesert, DistTok562GeneDesert, distToDHScluster))#,
# H2AFZ, H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me2, H3K4me3, H3K79me2, H3K9ac, H3K9me1, H3K9me3, H4K20me1, HDAC1, HDAC2, HDAC6))


quantile(normdataiPCR$expression, prob = seq(0, 1, length = 11), type = 5)[10]

featureROC.clean <- featureROC %>%
       filter(!is.na(repliSeq),
              !is.na(chromHMM),
              !is.na(NearestDHStype),
              !is.na(genomicDivision),
              !is.na(GenicLocation.simple),
              !is.na(AB_compartment_value)) %>%
       mutate(expression = ifelse(log10(RNA/DNA) >as.numeric(quantile(normdataiPCR$expression, prob = seq(0, 1, length = 11), type = 5)[8]), 1, ifelse(log10(RNA/DNA) <(as.numeric(summary(normdataiPCR$expression)[2])), 0,NA)), 
              DistToTSS = abs(DistToTSS),
              chromHMM = factor(chromHMM),
              geneLoc = factor(GenicLocation.simple),
              NearestDHStype = factor(NearestDHStype), 
              repliSeq = factor(repliSeq),
              genomicDivision = factor(genomicDivision),
              DistToNearestDHS = abs(DistToNearestDHS),
              DistToCGI = abs(DistToCGI),
              DistToNearestCTCF = abs(DistToNearestCTCF),
              distToDHScluster =abs(distToDHScluster))#,
              #H2AFZ= abs(H2AFZ),
              #H3K27ac = abs(H3K27ac), 
              #H3K27me3 = abs(H3K27me3), 
              #H3K36me3 = abs(H3K36me3), 
              #H3K4me1 =abs(H3K4me1), 
              #H3K4me2 = abs(H3K4me2), 
              #H3K4me3 = abs(H3K4me3),
              #H3K79me2 = abs(H3K79me2), 
              #H3K9ac =abs(H3K9ac), 
              #H3K9me1 =abs(H3K9me1), 
              #H3K9me3 =abs(H3K9me3), 
              #H4K20me1 =abs(H4K20me1), 
              #HDAC1 =abs(HDAC1), 
              #HDAC2 =abs(HDAC2), 
              #HDAC6 =abs(HDAC6))

featureROC.clean <- featureROC.clean[!is.na(featureROC.clean$expression),]

Expression <- featureROC.clean$expression
featureROC.clean <- subset(featureROC.clean, select=-c(expression, BC, DNA, RNA, GenicLocation.simple))


library(ROCR)
fit <- glm(formula=Expression ~. ,data=featureROC.clean, family='binomial')
#f_plot_ROC(fit$y,fit$fitted.values, main="")
preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
perf <- performance(preds, "auc")
rocC <- round(as.numeric(slot(perf, "y.values")), 3)
perf_plot <- performance(preds, "tpr", "fpr")

featureROC.clean.ROC <- data.frame(y=perf_plot@y.values[[1]] ,  x=perf_plot@x.values[[1]], SampleName=unique(normdataiPCR$SampleName))

AUC.features <- ggplot(featureROC.clean.ROC ,aes(x=x, y=y, colour=SampleName)) +
       geom_line(size = 2)+
       theme_classic()+
       ggtitle(opt$SampleName,"  \nAUC for bottom and top decile")+
       scale_colour_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       xlab('1-Specificity')+
       ylab('Sensitivity')+
       annotate('text', label = paste0('AUC = ',round(as.numeric(slot(perf, "y.values")), 3)), x = 0.5, y = 0.5,  color = 'black', size=10)+
       theme(legend.position="none")

cat('Plotting feature importance\n')
library(caret)
library(gridExtra)
library(cowplot)
#library(ggpubr)
impROC <- varImp(fit)
impROC$Feature <- rownames(impROC)
impROC$SampleName <- unique(normdataiPCR$SampleName)
impROC$Feature <- factor(impROC$Feature, levels= impROC[order(impROC$Overall),]$Feature)


Imp.features <- ggplot(impROC, aes(x=Feature, y=Overall, fill=SampleName)) +
       geom_bar(stat='identity', colour='black') +
       theme_classic() +
       scale_fill_brewer(palette='Set1') +
       ylab('Importance of features') +
       ggtitle(opt$SampleName,"  \nFeature importance")+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=9))+
       coord_flip()+
       theme(legend.position="none")



ggRepliSeqDistToDHS <- normdataiPCR %>%
       filter(!is.na(repliSeq)) %>%
       mutate(repliSeq = factor(repliSeq, levels=c('G1', 'S1', 'S2', 'S3', 'S4', 'G2'))) %>%
       ggplot(aes(x=repliSeq, y=abs(DistToNearestDHS)+1, fill=SampleName))+
       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=1, colour=sampleColor) +
       geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
       theme_classic()+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ggtitle(opt$SampleName,"  \nReplication timing")+
       xlab('Replication phase')+
       ylab('log10(distanceToDHS)') +
       theme(legend.position="none")+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm"))

ggABcompDistToDHS <- normdataiPCR %>%
       mutate(ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
       #mutate(expression = log10(RNA/DNA)) %>%
       filter(!is.na(ABcomp), ABcomp!='NA') %>%
       ggplot(aes(x = ABcomp, y = abs(DistToNearestDHS)+1, fill=ABcomp)) +
       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=21,aes(colour=ABcomp)) +
       geom_boxplot(outlier.shape=NA, alpha=0.6) +
       theme_classic()+
       scale_fill_brewer(palette='Set1') +
       scale_colour_brewer(palette='Set1') +
       ggtitle(opt$SampleName,"  \nA/B compartments ")+
       xlab('Genomic compartment')+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       ylab('log10(distanceToDHS)')+
       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
       theme(legend.position="none")



ggABcompNoLegend <- ggABcomp +  
                     theme(legend.position="none") +
                     ggtitle(opt$SampleName,"  \nA/B compartments ")


pdf(file=paste0(opt$OUTDIR,"/graphs/AUC_Features.pdf"), width=15, height=15)
grid.arrange(AUC.features, Imp.features, ggrepliSeq, ggABcompNoLegend,ggDistToNearestDHS, ggRepliSeqDistToDHS,ggABcompDistToDHS, ncol=3, nrow = 3, 
layout_matrix=rbind( c(1, 2, 2), 
                     c(3, 4, 5),
                     c(6, 7, NA)))
dev.off()



topDecile <- head(normdataiPCR[rev(order(normdataiPCR$expression)),],nrow(normdataiPCR)*0.1)
bottomDecile <- tail(normdataiPCR[rev(order(normdataiPCR$expression)),],nrow(normdataiPCR)*0.1)


cat('writing topDecile','\n')
write.table(subset(topDecile, select=c(chrom, chromStart, chromEnd, BC)), file=paste0(opt$OUTDIR,"/topDecile.bed"), sep='\t', quote=F, col.names=F, row.names=F)
cat('writing bottomDecile','\n')
write.table(subset(bottomDecile, select=c(chrom, chromStart, chromEnd, BC)), file=paste0(opt$OUTDIR,"/bottomDecile.bed"), sep='\t', quote=F, col.names=F, row.names=F)

#randomDecile <- read.table('/tmp/randomGeneticDivision.bed', sep='\t')
decilesGeneRegion <- rbind(data.frame(table(bottomDecile$genomicDivision),Decile='Bottom'),
       data.frame(table(topDecile$genomicDivision),Decile='Top'))#,
#       data.frame(table(randomDecile$V4),Decile='Random'))


ggdecilesGeneRegion <- decilesGeneRegion %>%
       ggplot(aes(x=Var1, y=Freq, fill=Decile)) +
       geom_bar(stat='identity',  position=position_dodge(), colour="black", alpha=0.6) +
       theme_classic() +
       scale_fill_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
       xlab('Genomic region') +
       ylab('Number of insertions') +
       ggtitle(opt$SampleName,"Genomic region per top and bottom expression decile")

pdf(file=paste0(opt$OUTDIR,"/graphs/decilesGeneRegion.pdf"), width=5, height=5)
ggdecilesGeneRegion
dev.off()

topDecile$Decile <- 'Top'
bottomDecile$Decile <- 'Bottom'
DecilesDistance <- rbind(topDecile, bottomDecile)

ggDecilesDistance <- DecilesDistance %>%  
       mutate(Decile = factor(Decile, levels=c('Top', 'Bottom'))) %>%
       ggplot(aes(x=Decile, y=log10(abs(DistToNearestDHS)+1), fill=Decile)) +
       geom_boxplot(outlier.shape=NA, alpha=0.6) +
       theme_classic() +
       scale_fill_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
       xlab('Expression decile') +
       ylab('log10(distToNearestDHS') +
       ggtitle(opt$SampleName,"Distance to DHS per top and bottom expression decile")

pdf(file=paste0(opt$OUTDIR,"/graphs/decilesDistToDHS.pdf"), width=5, height=5)
ggDecilesDistance
dev.off()

pdf(file=paste0(opt$OUTDIR,"/graphs/DecilesTopBottom.pdf"), width=10, height=10)
grid.arrange(ggdecilesGeneRegion, ggDecilesDistance, AUC.features ,Imp.features, ncol=2)
dev.off()


cat('\nDivide deciles into genomic regions\n')
GeneDesertDecile <- normdataiPCR[normdataiPCR$genomicDivision=='geneDesert',]
GeneDesertDecile <- GeneDesertDecile[!is.na(GeneDesertDecile$expression),]

topGeneDesertDecile <- head(GeneDesertDecile[rev(order(GeneDesertDecile$expression)),], nrow(GeneDesertDecile)*0.1)
topGeneDesertDecile$BC <- paste0(formatC(1:nrow(topGeneDesertDecile), width = 3, format = "d", flag = "0"), topGeneDesertDecile$BC)
bottomGeneDesertDecile <- tail(GeneDesertDecile[rev(order(GeneDesertDecile$expression)),],nrow(GeneDesertDecile)*0.1)
bottomGeneDesertDecile$BC <- paste0(formatC(1:nrow(bottomGeneDesertDecile), width = 3, format = "d", flag = "0"), bottomGeneDesertDecile$BC)

cat('writing gene Desert topDecile','\n')
write.table(subset(topGeneDesertDecile, select=c(chrom, chromStart, chromEnd, BC)), file=paste0(opt$OUTDIR,"/GeneDesertTopDecile.bed"), sep='\t', quote=F, col.names=F, row.names=F)
cat('writing gene Desert bottomDecile','\n')
write.table(subset(bottomGeneDesertDecile, select=c(chrom, chromStart, chromEnd, BC)), file=paste0(opt$OUTDIR,"/GeneDesertBottomDecile.bed"), sep='\t', quote=F, col.names=F, row.names=F)

topGeneDesertDecile$Decile <- 'Top'
bottomGeneDesertDecile$Decile <- 'Bottom'
DecilesGeneDesertDistance <- rbind(topGeneDesertDecile, bottomGeneDesertDecile)

ggDecilesGeneDesertDistance<- DecilesGeneDesertDistance %>%  
       mutate(Decile = factor(Decile, levels=c('Top', 'Bottom'))) %>%
       ggplot(aes(x=Decile, y=log10(abs(DistToNearestDHS)+1), fill=Decile)) +
       geom_boxplot(outlier.shape=NA, alpha=0.6) +
       theme_classic() +
       scale_fill_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
       xlab('Expression decile') +
       ylab('log10(distToNearestDHS') +
       ggtitle(opt$SampleName,"Distance to DHS per top and bottom expression decile\nof all insertions in Gene deserts")+
       annotation_logticks(sides = "lb",long = unit(0.2, "cm"))+
       ylim(0,7)



geneRichDecile <- normdataiPCR[normdataiPCR$genomicDivision=='geneRich',]
geneRichDecile <- geneRichDecile[!is.na(geneRichDecile$expression),]

topgeneRichDecile <- head(geneRichDecile[rev(order(geneRichDecile$expression)),], nrow(geneRichDecile)*0.1)
topgeneRichDecile$BC <- paste0(formatC(1:nrow(topgeneRichDecile), width = 3, format = "d", flag = "0"), topgeneRichDecile$BC)
bottomgeneRichDecile <- tail(geneRichDecile[rev(order(geneRichDecile$expression)),],nrow(geneRichDecile)*0.1)
bottomgeneRichDecile$BC <- paste0(formatC(1:nrow(bottomgeneRichDecile), width = 3, format = "d", flag = "0"), bottomgeneRichDecile$BC)

cat('writing gene Rich topDecile','\n')
write.table(subset(topgeneRichDecile, select=c(chrom, chromStart, chromEnd, BC)), file=paste0(opt$OUTDIR,"/geneRichTopDecile.bed"), sep='\t', quote=F, col.names=F, row.names=F)
cat('writing gene Rich bottomDecile','\n')
write.table(subset(bottomgeneRichDecile, select=c(chrom, chromStart, chromEnd, BC)), file=paste0(opt$OUTDIR,"/geneRichBottomDecile.bed"), sep='\t', quote=F, col.names=F, row.names=F)



topgeneRichDecile$Decile <- 'Top'
bottomgeneRichDecile$Decile <- 'Bottom'
DecilesgeneRichDistance <- rbind(topgeneRichDecile, bottomgeneRichDecile)

ggDecilesgeneRichDistance<- DecilesgeneRichDistance %>%  
       mutate(Decile = factor(Decile, levels=c('Top', 'Bottom'))) %>%
       ggplot(aes(x=Decile, y=log10(abs(DistToNearestDHS)+1), fill=Decile)) +
       geom_boxplot(outlier.shape=NA, alpha=0.6) +
       theme_classic() +
       scale_fill_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
       xlab('Expression decile') +
       ylab('log10(distToNearestDHS') +
       ggtitle(opt$SampleName,"Distance to DHS per top and bottom expression decile\nof all insertions in Gene rich ares")+
       annotation_logticks(sides = "lb",long = unit(0.2, "cm"))+
       ylim(0,7)
       

pdf(file=paste0(opt$OUTDIR,"/graphs/DecilesTopBottom_geneDesetRich.pdf"), width=10, height=5)
grid.arrange(ggDecilesGeneDesertDistance, ggDecilesgeneRichDistance, ncol=2)
dev.off()


bottomNearDHS <- normdataiPCR[normdataiPCR$expression<(-0.5) & abs(normdataiPCR$DistToNearestDHS)< 1e4,]
write.table(subset(bottomNearDHS, select=c(chrom, chromStart, chromEnd, BC)), file=paste0(opt$OUTDIR,"/bottomNearDHS.bed"), sep='\t', quote=F, col.names=F, row.names=F)








######
#Create AUROC curves of distance to DHS types
######
##Create a function to return colname and AUROC for each column in a data frame
cat('Creating AUROC for distance to DHS types\n')
ROCExpression <- function(data, BCcol) { 
       library(ROCR)
       library(dplyr)
       library(reshape2)
       data <- data[data[,BCcol] %in% normdataiPCR$BC,]
       data$expression <- normdataiPCR[match(data[,BCcol],normdataiPCR$BC),]$expression
       dataROC <- subset(data, select=-c(BCcol))
       ROCdf <- as.data.frame(matrix(ncol=(ncol(dataROC)-1),nrow=1))
       rownames(ROCdf) <- 'AUROC'
       for (i in 1:c(ncol(dataROC)-1)) {
              #cat(colnames(dataROC)[i],'\n')
              dataROCglm <- dataROC %>%
              mutate(expression = ifelse(expression>as.numeric(summary(dataROC$expression)[5]), 1, ifelse(expression<(as.numeric(summary(dataROC$expression)[2])), 0,NA)))
              
              dataROCglm <- subset(dataROCglm, select=c(colnames(dataROCglm)[i], 'expression'))
             
              dataROCglm <- dataROCglm[!is.na(dataROCglm$expression),]
              Expression <- dataROCglm$expression
              dataROCglm <- subset(dataROCglm, select=-c(expression))
              fit2 <- glm(formula=Expression ~. ,data=dataROCglm, family='binomial')
       
              preds <- prediction(as.numeric(fit2$fitted.values), as.numeric(fit2$y))
              perf <- performance(preds, "auc")
              rocC <- round(as.numeric(slot(perf, "y.values")), 3)
              #cat(rocC, '\n')
              colnames(ROCdf)[i] <- colnames(dataROC)[i]
              ROCdf[,i] <- rocC
              
       }
       return(ROCdf)
}
       

ROCExpressionDirection <- function(data, BCcol) { 
       library(ROCR)
       library(dplyr)
       library(reshape2)
       data <- data[data[,BCcol] %in% normdataiPCR$BC,]
       data$expression <- normdataiPCR[match(data[,BCcol],normdataiPCR$BC),]$expression
       dataROC <- subset(data, select=-c(BCcol))
       ROCdf <- as.data.frame(matrix(ncol=(ncol(dataROC)-1),nrow=2))
       rownames(ROCdf) <- c('AUROC', 'Estimate')
       for (i in 1:c(ncol(dataROC)-1)) {
              #cat(colnames(dataROC)[i],'\n')
              dataROCglm <- dataROC %>%
              mutate(expression = ifelse(expression>as.numeric(summary(dataROC$expression)[5]), 1, ifelse(expression<(as.numeric(summary(dataROC$expression)[2])), 0,NA)))
              
              dataROCglm <- subset(dataROCglm, select=c(colnames(dataROCglm)[i], 'expression'))
             
              dataROCglm <- dataROCglm[!is.na(dataROCglm$expression),]
              Expression <- dataROCglm$expression
              dataROCglm <- subset(dataROCglm, select=-c(expression))
              fit2 <- glm(formula=Expression ~. ,data=dataROCglm, family='binomial')
       
              preds <- prediction(as.numeric(fit2$fitted.values), as.numeric(fit2$y))
              perf <- performance(preds, "auc")
              rocC <- round(as.numeric(slot(perf, "y.values")), 3)
              #cat(rocC, '\n')
              colnames(ROCdf)[i] <- colnames(dataROC)[i]
              ROCdf[,i][1] <- rocC
              if (length(as.numeric(coef(summary(fit2))[,1][-1]))>1) {
                     ROCdf[,i][2] <-     mean(as.numeric(coef(summary(fit2))[,1][-1]))
              } else { 
                     ROCdf[,i][2] <- as.numeric(coef(summary(fit2))[,1][-1])
              }
       }
       return(ROCdf)
}
       


library(ROCR)
library(reshape2)
library(dplyr)
dhsTypeProm <- read.table(paste0(opt$OUTDIR,'/', opt$SampleName,'_DHSpromoter.bed'), sep='\t', stringsAsFactors=F)
dhsTypeProm$dhsType <- rep(paste0('Promoter_',formatC(1:5, width = 2, format = "d", flag = "0")), length(unique(dhsTypeProm$V1)))


dhsTypeCTCF <- read.table(paste0(opt$OUTDIR,'/', opt$SampleName,'_DHSctcf.bed'), sep='\t', stringsAsFactors=F)
dhsTypeCTCF$dhsType <- rep(paste0('CTCF_',formatC(1:5, width = 2, format = "d", flag = "0")), length(unique(dhsTypeCTCF$V1)))

dhsTypeDistal <- read.table(paste0(opt$OUTDIR,'/', opt$SampleName,'_DHSdistal.bed'), sep='\t', stringsAsFactors=F)
dhsTypeDistal$dhsType <- rep(paste0('Distal_',formatC(1:5, width = 2, format = "d", flag = "0")), length(unique(dhsTypeDistal$V1)))

dhsTypeall <- read.table(paste0(opt$OUTDIR,'/', opt$SampleName,'_DHSall.bed'), sep='\t', stringsAsFactors=F)
dhsTypeall$dhsType <- rep(paste0('allTogether_',formatC(1:5, width = 2, format = "d", flag = "0")), length(unique(dhsTypeall$V1)))

dhsTypenonCTCF <- read.table(paste0(opt$OUTDIR,'/', opt$SampleName,'_nonCTCF.bed'), sep='\t', stringsAsFactors=F)
dhsTypenonCTCF$dhsType <- rep(paste0('nonCTCF_',formatC(1:5, width = 2, format = "d", flag = "0")), length(unique(dhsTypenonCTCF$V1)))



dhsTypeROC <- rbind(dhsTypeProm, dhsTypeCTCF, dhsTypeDistal, dhsTypeall, dhsTypenonCTCF)


dhsTypeROC$V3 <- log10(abs(dhsTypeROC$V3)+1)
cat('\nCreate data frame with distance to dhsType\n')
distdhsTypeROC <- dcast(data = dhsTypeROC,formula = V1~dhsType,fun.aggregate = sum,value.var = "V3")

ROCdf <- ROCExpression(distdhsTypeROC, 1)
library(dplyr)
distdhsTypeROCglm <- distdhsTypeROC %>%
       filter(V1 %in% normdataiPCR$BC) %>%
       mutate(expression = normdataiPCR[match(V1, normdataiPCR$BC),]$expression) %>%
       mutate(expression = ifelse(expression>as.numeric(summary(expression)[5]), 1, ifelse(expression<(as.numeric(summary(expression)[2])), 0,NA))) %>%
       filter(!is.na(expression))
       
Expression <- distdhsTypeROCglm$expression
distdhsTypeROCsubset <- subset(distdhsTypeROCglm, select=-c(expression, V1))


fit <- glm(formula=Expression ~. ,data=subset(distdhsTypeROCsubset, select =c(1:5)), family='binomial')
preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
perf <- performance(preds, "auc")
rocC <- round(as.numeric(slot(perf, "y.values")), 3)
ROCdf[,26] <- rocC


fit <- glm(formula=Expression ~. ,data=subset(distdhsTypeROCsubset, select =c(6:10)), family='binomial')
preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
perf <- performance(preds, "auc")
rocC <- round(as.numeric(slot(perf, "y.values")), 3)
ROCdf[,27] <- rocC

fit <- glm(formula=Expression ~. ,data=subset(distdhsTypeROCsubset, select =c(11:15)), family='binomial')
preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
perf <- performance(preds, "auc")
rocC <- round(as.numeric(slot(perf, "y.values")), 3)
ROCdf[,28] <- rocC

fit <- glm(formula=Expression ~. ,data=subset(distdhsTypeROCsubset, select =c(15:20)), family='binomial')
preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
perf <- performance(preds, "auc")
rocC <- round(as.numeric(slot(perf, "y.values")), 3)
ROCdf[,29] <- rocC

fit <- glm(formula=Expression ~. ,data=subset(distdhsTypeROCsubset, select =c(21:25)), family='binomial')
preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
perf <- performance(preds, "auc")
rocC <- round(as.numeric(slot(perf, "y.values")), 3)
ROCdf[,30] <- rocC


fit <- glm(formula=Expression ~. ,data=subset(distdhsTypeROCsubset, select =c(6:15,21:25)), family='binomial')
preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
perf <- performance(preds, "auc")
rocC <- round(as.numeric(slot(perf, "y.values")), 3)
ROCdf[,31] <- rocC

colnames(ROCdf)[26] <- 'allTogether_nearest-05'
#colnames(ROCdf)[52] <- 'allTogether_nearest-10'
colnames(ROCdf)[27] <- 'CTCF_nearest-05'
#colnames(ROCdf)[54] <- 'CTCF_nearest-10'
colnames(ROCdf)[28] <- 'Distal_nearest-05'
#colnames(ROCdf)[56] <- 'Distal_nearest-10'
colnames(ROCdf)[29] <- 'nonCTCF_nearest-05'
#colnames(ROCdf)[58] <- 'nonCTCF_nearest-10'
colnames(ROCdf)[30] <- 'Promoter_nearest-05'
#colnames(ROCdf)[60] <- 'Promoter_nearest-10'
colnames(ROCdf)[31] <- 'allSeparately_nearest-05'
#colnames(ROCdf)[62] <- 'allSeparately_nearest-10'

NearestDHStype.ROC <- ROCdf
pdf(paste0(opt$OUTDIR,'/graphs/AUROC_DHStypes.pdf'), height=5, width=5)
ROCdf  %>%
       melt() %>%
       mutate(DHS_Type = gsub('_.*' ,'' , variable),
              variable=paste0('DHS_', gsub('.*_' ,'' , variable)),
              DHS_Type= factor(DHS_Type, levels=c('CTCF', 'Promoter', 'Distal' ,'nonCTCF', 'allTogether', 'allSeparately'))) %>%
       ggplot(aes(x=variable, y=value, fill=DHS_Type)) +
              geom_bar(stat='identity', position=position_dodge(), colour='black', alpha=0.6) +
              xlab('DHS neighbour') +
              ylab('AUROC') +
              coord_cartesian(ylim = c(0.5, 0.9)) +
              scale_fill_brewer(palette='Set1') +
              theme_classic() +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, ,hjust=0.95,vjust=0.95), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, subtitle='AUROC for nearest 5 DHS type individually\nas well as combined') #+
              #facet_zoom(y=SampleName=='CMV',zoom.size = 2, ylim=c(0,3))
       
dev.off()






cat('Binned distance around insertions\n')
library(data.table)
cat('nDHS\n')
linDist <- read.table(paste0(opt$OUTDIR,'/', opt$SampleName,'_linearDistanceDHS_bedtools.bed'), sep='\t', stringsAsFactors=F)

linDistSplit <- split.data.frame(linDist, linDist$V4)
#loci <- linDistSplit[[1]]
lociDistance <- lapply(linDistSplit, function(loci) { 
       loci$distance <- 'NA'
       if (nrow(loci)==10) {
             loci$distance <- c(seq(from=-5e4, to=-1e4, by=1e4), seq(from=10000, to=5e4, by=1e4))
             loci
       }
})

#lociDistanceDF <-  do.call(rbind, lociDistance)
lociDistanceDF <-  rbindlist(lociDistance)

lociDistance <- lociDistanceDF
lociDistance$distance <- abs(lociDistance$distance)
distLoci <- dcast(data = lociDistance,formula = V4~distance,fun.aggregate = sum,value.var = "V6")
rownames(distLoci) <- distLoci$V4
nDHSbins <-distLoci

cat('tagsDHSbins\n')
linDist <- read.table(paste0(opt$OUTDIR,'/', opt$SampleName,'_DistDensity_bedtools.bed'), sep='\t', stringsAsFactors=F)

linDistSplit <- split.data.frame(linDist, linDist$V4)
#loci <- linDistSplit[[1]]
lociDistance <- lapply(linDistSplit, function(loci) { 
       loci$distance <- 'NA'
       if (nrow(loci)==10) {
             loci$distance <- c(seq(from=-5e4, to=-1e4, by=1e4), seq(from=10000, to=5e4, by=1e4))
             loci
       }
})

lociDistanceDF <-  rbindlist(lociDistance)
lociDistance <- lociDistanceDF
lociDistance$V6 <- as.numeric(lociDistance$V6)
lociDistance$distance <- abs(lociDistance$distance)
lociDistance[lociDistance$V6=='NaN',]$V6 <- 0

distLoci <- dcast(data = lociDistance,formula = V4~distance,fun.aggregate = sum,value.var = "V6")
rownames(distLoci) <- distLoci$V4
tagsDHSbins <-distLoci

cat('maxTagsDHSbins\n')
linDist <- read.table(paste0(opt$OUTDIR,'/', opt$SampleName,'_DistDensityMax_bedtools.bed'), sep='\t', stringsAsFactors=F)

linDistSplit <- split.data.frame(linDist, linDist$V4)
#loci <- linDistSplit[[1]]
lociDistance <- lapply(linDistSplit, function(loci) { 
       loci$distance <- 'NA'
       if (nrow(loci)==10) {
             loci$distance <- c(seq(from=-5e4, to=-1e4, by=1e4), seq(from=10000, to=5e4, by=1e4))
             loci
       }
})

lociDistanceDF <-  rbindlist(lociDistance)
lociDistance <- lociDistanceDF
lociDistance$V6 <- as.numeric(lociDistance$V6)
lociDistance[lociDistance$V6=='NaN',]$V6 <- 0
lociDistance$distance <- abs(lociDistance$distance)
distLoci <- dcast(data = lociDistance,formula = V4~distance,fun.aggregate = sum,value.var = "V6")
rownames(distLoci) <- distLoci$V4
maxTagsDHSbins <-distLoci

cat('totalTagsbins\n')
linDist <- read.table(paste0(opt$OUTDIR,'/', opt$SampleName,'_10kbWindows_tagsFromBam.bed'), sep='\t', stringsAsFactors=F)

linDistSplit <- split.data.frame(linDist, linDist$V4)
#loci <- linDistSplit[[1]]
lociDistance <- lapply(linDistSplit, function(loci) { 
       loci$distance <- 'NA'
       if (nrow(loci)==10) {
             loci$distance <- c(seq(from=-5e4, to=-1e4, by=1e4), seq(from=10000, to=5e4, by=1e4))
             loci
       }
})

lociDistanceDF <-  rbindlist(lociDistance)
lociDistance <- lociDistanceDF

lociDistance$V6 <- as.numeric(lociDistance$V6)
lociDistance[lociDistance$V6=='NaN',]$V6 <- 0
lociDistance$distance <- abs(lociDistance$distance)
distLoci <- dcast(data = lociDistance,formula = V4~distance,fun.aggregate = sum,value.var = "V6")
rownames(distLoci) <- distLoci$V4
totalTagsbins <-distLoci



nDHSbins.subset <- subset(nDHSbins, select=c(1:6))
colnames(nDHSbins.subset) <- paste0('numberDHS_', colnames(nDHSbins.subset))

tagsDHSbins.subset <- subset(tagsDHSbins, select=c(1:6))
colnames(tagsDHSbins.subset) <- paste0('tagsDHS_', colnames(tagsDHSbins.subset))

maxTagsDHSbins.subset <- subset(maxTagsDHSbins, select=c(1:6))
colnames(maxTagsDHSbins.subset) <- paste0('maxTagsDHS_', colnames(maxTagsDHSbins.subset))

totalTagsbins.subset <- subset(totalTagsbins, select=c(1:6))
colnames(totalTagsbins.subset) <- paste0('totalTags_', colnames(totalTagsbins.subset))

DHSbinsAnalysis <- cbind(nDHSbins.subset, tagsDHSbins.subset, maxTagsDHSbins.subset, totalTagsbins.subset)
DHSbinsAnalysis <- subset(DHSbinsAnalysis, select=-c(tagsDHS_V4, maxTagsDHS_V4, totalTags_V4))
DHSbinsAnalysis <- as.data.frame(DHSbinsAnalysis)
#TODO
ROCdhsAnalysis <- ROCExpression(DHSbinsAnalysis, 1)
library(dplyr)
DHSbinsAnalysisglm <- DHSbinsAnalysis %>%
       filter(numberDHS_V4 %in% normdataiPCR$BC) %>%
       mutate(expression = normdataiPCR[match(numberDHS_V4, normdataiPCR$BC),]$expression) %>%
       mutate(expression = ifelse(expression>as.numeric(summary(expression)[5]), 1, ifelse(expression<(as.numeric(summary(expression)[2])), 0,NA))) %>%
       filter(!is.na(expression))
       
Expression <- DHSbinsAnalysisglm$expression
DHSbinsAnalysisglm <- subset(DHSbinsAnalysisglm, select=-c(numberDHS_V4))


fit <- glm(formula=Expression ~. ,data=subset(DHSbinsAnalysisglm, select =c(1:5)), family='binomial')
preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
perf <- performance(preds, "auc")
rocC <- round(as.numeric(slot(perf, "y.values")), 3)
ROCdhsAnalysis[,21] <- rocC

fit <- glm(formula=Expression ~. ,data=subset(DHSbinsAnalysisglm, select =c(6:10)), family='binomial')
preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
perf <- performance(preds, "auc")
rocC <- round(as.numeric(slot(perf, "y.values")), 3)
ROCdhsAnalysis[,22] <- rocC

fit <- glm(formula=Expression ~. ,data=subset(DHSbinsAnalysisglm, select =c(11:15)), family='binomial')
preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
perf <- performance(preds, "auc")
rocC <- round(as.numeric(slot(perf, "y.values")), 3)
ROCdhsAnalysis[,23] <- rocC

fit <- glm(formula=Expression ~. ,data=subset(DHSbinsAnalysisglm, select =c(16:20)), family='binomial')
preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
perf <- performance(preds, "auc")
rocC <- round(as.numeric(slot(perf, "y.values")), 3)
ROCdhsAnalysis[,24] <- rocC


colnames(ROCdhsAnalysis)[21] <- 'numberDHS_All'
colnames(ROCdhsAnalysis)[22] <- 'tagsDHS_All'
colnames(ROCdhsAnalysis)[23] <- 'maxTagsDHS_All'
colnames(ROCdhsAnalysis)[24] <- 'totalTags_All'


pdf(paste0(opt$OUTDIR,'/graphs/AUROC_BinnedDistance.pdf'), height=5, width=5)
ROCdhsAnalysis  %>%
       melt() %>%
       mutate(Type = gsub('_.*' ,'' , variable),
              variable=gsub('.*_' ,'' , variable),
              variable= factor(variable, levels=c('10000', '20000', '30000' ,'40000', '50000', '60000', '70000', '80000', '90000', '1e+05', 'All'))) %>%
       ggplot(aes(x=variable, y=value, fill=Type)) +
              geom_bar(stat='identity', position=position_dodge(), colour='black', alpha=0.6) +
              xlab('Binned distance from insertion') +
              ylab('AUROC') +
              coord_cartesian(ylim = c(0.5, 0.9)) +
              scale_fill_brewer(palette='Set1') +
              theme_classic() +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, ,hjust=0.95,vjust=0.95), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, subtitle='AUROC for number of DHS,\nand DNase accessibility surrounding the insertion') #+
              #facet_zoom(y=SampleName=='CMV',zoom.size = 2, ylim=c(0,3))
       
dev.off()





#####
#AUROC for each independent feature
#####
cat('AUROC for each individual feature\n')
dhsTypeall <- read.table(paste0(opt$OUTDIR,'/', opt$SampleName,'_DHSall.bed'), sep='\t', stringsAsFactors=F)
dhsTypeall$dhsType <- rep(paste0('allTogether_',formatC(1:5, width = 2, format = "d", flag = "0")), length(unique(dhsTypeall$V1)))

dhsTypeallExp <- dcast(data = dhsTypeall,formula = V1~dhsType,fun.aggregate = sum,value.var = "V3")

dhsTypeallExp$meanDist <- rowSums(dhsTypeallExp[,2:6])/5

dhsTypeallExp <- dhsTypeallExp[dhsTypeallExp$V1 %in% normdataiPCR$BC ,]
dhsTypeallExp$expression <- normdataiPCR[match(dhsTypeallExp$V1, normdataiPCR$BC),]$expression



totalTagsbins.subset$meanTags <- rowSums(totalTagsbins.subset[,2:6])/5
normdataiPCR$totalTags50kb <- totalTagsbins.subset[match(normdataiPCR$BC, totalTagsbins.subset$totalTags_V4) ,]$meanTags
normdataiPCR$meanDist5DHS <- dhsTypeallExp[match(normdataiPCR$BC, dhsTypeallExp$V1) ,]$meanDist

featureROC <- subset(normdataiPCR,select=c(BC, DNA, RNA, DistToTSS, DistToNearestDHS, DistToNearestCTCF, DistToCGI, PercentGC, AB_compartment_value,
 chromHMM, nDHS100kb, nDHS10kb, nDHS300kb, nDHS500kb, nDHS1Mb, nDHS100kbDensity, nDHS300kbsum, nDHS500kbDensity, nDHS1MbDensity, nTSS100kb, 
 repliSeq, nExpressed100kb, nCGI100kb, nCTCFpeaks100kb, DNasePksDens, NearestDHStype, CTCFPksDens, genomicDivision, DistToGeneDesert, DistTok562GeneDesert, meanDist5DHS, distToDHScluster, 
 H2AFZ, H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me2, H3K4me3, H3K79me2, H3K9ac, H3K9me1, H3K9me3, H4K20me1, HDAC1, HDAC2, HDAC6, totalTags50kb ))
 



featureROC.clean <- featureROC %>%
       filter(!is.na(repliSeq),
              !is.na(chromHMM),
              !is.na(NearestDHStype),
              !is.na(genomicDivision),
              !is.na(AB_compartment_value)) %>%
              mutate(expression = log10(RNA/DNA)) %>%
       mutate(expression = ifelse(log10(RNA/DNA) >as.numeric(summary(expression)[5]), 1, ifelse(log10(RNA/DNA) <(as.numeric(summary(expression)[2])), 0,NA)), 
              DistToTSS = abs(DistToTSS),
              chromHMM = factor(chromHMM),
              NearestDHStype = factor(NearestDHStype), 
              repliSeq = factor(repliSeq),
              genomicDivision = factor(genomicDivision),
              DistToNearestDHS = abs(DistToNearestDHS),
              DistToCGI = abs(DistToCGI),
              DistToNearestCTCF = abs(DistToNearestCTCF),
              distToDHScluster =abs(distToDHScluster),
              H2AFZ= abs(H2AFZ),
              H3K27ac = abs(H3K27ac), 
              H3K27me3 = abs(H3K27me3), 
              H3K36me3 = abs(H3K36me3), 
              H3K4me1 =abs(H3K4me1), 
              H3K4me2 = abs(H3K4me2), 
              H3K4me3 = abs(H3K4me3),
              H3K79me2 = abs(H3K79me2), 
              H3K9ac =abs(H3K9ac), 
              H3K9me1 =abs(H3K9me1), 
              H3K9me3 =abs(H3K9me3), 
              H4K20me1 =abs(H4K20me1), 
              HDAC1 =abs(HDAC1), 
              HDAC2 =abs(HDAC2), 
              HDAC6 =abs(HDAC6))

#featureROC.clean <- featureROC.clean[!is.na(featureROC.clean$expression),]

featureROClog <- featureROC.clean 
featureROClog <- featureROClog[!is.na(featureROClog$totalTags50kb),]
Expression <- featureROClog$expression
featureROClog <- subset(featureROClog, select=-c(expression))


logCol <- c("DistToTSS", "DistToNearestDHS", "DistToNearestCTCF", 
"DistToCGI",  
"nDHS100kbDensity", "nDHS300kbsum", "nDHS500kbDensity", "nDHS1MbDensity", 
"CTCFPksDens",
 "meanDist5DHS", "distToDHScluster", 
"H2AFZ", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me2", 
"H3K4me3", "H3K79me2", "H3K9ac", "H3K9me1", "H3K9me3", "H4K20me1", 
"HDAC1", "HDAC2", "HDAC6", "totalTags50kb", "DNA", "RNA")

featureROClog[,logCol] <- log10(abs(featureROClog[,logCol])+1)
#Expression <- featureROClog$expression
#featureROClog <- subset(featureROClog, select=-c(expression,  DNA, RNA, GenicLocation.simple))

ROC.featureROC <- ROCExpression(subset(featureROClog, DNA>=0 & RNA>=0, select=-c(DNA, RNA)), 1)
ROCExpressionDirection(featureROClog, 1)
names(ROC.featureROC) <- gsub('Density', 'Accessibility', names(ROC.featureROC)) 
names(ROC.featureROC) <- gsub('Dens', 'Accessibility', names(ROC.featureROC)) 
ROC.featureROC <- ROC.featureROC[,!names(ROC.featureROC) %in% c('DistToGeneDesert','DistTok562GeneDesert', 'geneLoc','nDHS300kbsum', 'meanDist5DHS', 'distToDHScluster')]
#ROC.featureROC <- cbind(ROC.featureROC ,NearestDHStype.ROC[,c(26, 27, 28, 29, 30, 31)])

ROC.featureROC.Order <- names(ROC.featureROC[,rev(order(ROC.featureROC) )] )

ROC.featureROC <- ROC.featureROC[,order(names(ROC.featureROC) )] 

ROC.featureROC  <- ROC.featureROC   %>%
       melt()
ROC.featureROC$Feature <- 'NA'
ROC.featureROC[grep("^n.",ROC.featureROC$variable),]$Feature <- 'NumberOfFeatures'
ROC.featureROC[grep("Access|totalTags",ROC.featureROC$variable),]$Feature <- 'Accessibility'
ROC.featureROC[grep("H3|H4|HDAC|H2AF",ROC.featureROC$variable),]$Feature <- 'Histones'
ROC.featureROC[grep("Dist|dist|Nearest|nearest",ROC.featureROC$variable),]$Feature <- 'DistanceToNearest' 
ROC.featureROC[grep("^n|Dist|dist|Nearest|nearest|Access|totalTags|H3|H4|HDAC|H2AF",ROC.featureROC$variable, invert=T),]$Feature <- 'Other' 

featureClassNumber <- normdataiPCR %>%
       filter(log10(DNA)>=0,
       log10(RNA)>=0) %>%
       mutate(expression = ifelse(log10(RNA/DNA) >as.numeric(summary(expression)[5]), 1, ifelse(log10(RNA/DNA) <(as.numeric(summary(expression)[2])), 0,NA))) %>%
       filter(!is.na(expression))
       
featureClassNumber <- table(featureClassNumber$expression)


cat('plotting AUROC_SelectedFeatures.pdf\n')
pdf(paste0(opt$OUTDIR,'/graphs/AUROC_SelectedFeatures.pdf'), height=8, width=10)
ROC.featureROC %>%
       mutate(SampleName=opt$SampleName,
              variable= factor(variable, levels=rev(ROC.featureROC.Order))) %>%
       ggplot(aes(x=variable, y=value, fill=Feature)) +
              geom_bar(stat='identity', position=position_dodge(), colour='black', alpha=0.6) +
              xlab('Features') +
              ylab('AUROC') +
              #facet_wrap(~Feature)+
              #coord_cartesian(ylim = c(0.5, 0.9)) +
              scale_fill_brewer(palette='Set1') +
              theme_classic() +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, subtitle=paste0('AUROC for each individual feature classifying top and bottom 25% insertions based on gene expression\n',
                      'top25% = ',as.numeric(featureClassNumber)[1], ', ', 'bottom25% = ',as.numeric(featureClassNumber)[2]))+
              coord_flip(ylim=c(0.5, 0.95))+
              geom_text(aes(label=value), hjust=-0.1)
             # theme(legend.position="none")
       
dev.off()


cat('plotting DNAcopiesVsActivity.pdf\n')
pdf(paste0(opt$OUTDIR,'/graphs/DNAcopiesVsActivity.pdf'), height=8, width=8)
scatter.smooth(log10(normdataiPCR$DNA), normdataiPCR$expression, col=sampleColor, xlab='Log10(DNA copies)', ylab='Reporter activity', main=paste0(opt$SampleName,'\nthe reporter activity as a function of number of DNA copies'), frame=F)
dev.off()
#####
#all features together
#####
#cat('DNase and Histone tags in promoter\n')
#tagsProm <- read.table(paste0(opt$OUTDIR, '/histoneDNaseTags.bed'), sep='\t', header=T, stringsAsFactors=F)
#tagsProm <- subset(tagsProm, select=-c(chrom, start, end))
#
#
#allFeatures <- featureROClog 
#allFeatures <- allFeatures %>%
#       left_join(tagsProm, by='BC') 
#allFeatures$Expression1 <- normdataiPCR[match(allFeatures$BC, normdataiPCR$BC),]$expression
#allFeatures.ROC <-  allFeatures %>%
#mutate(Expression1 = ifelse(Expression1>as.numeric(summary(Expression1)[5]), 1, ifelse(Expression1<(as.numeric(summary(Expression1)[2])), 0,NA))) %>%
#filter(!is.na(Expression1))
#
#fit <- glm(formula=Expression1 ~. ,data=subset(allFeatures.ROC, select=-c(BC)), family='binomial')
#preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
#perf <- performance(preds, "auc")
#allFeatures.ROC.results <- round(as.numeric(slot(perf, "y.values")), 3)
#

#ROCExpression(subset(allFeatures.ROC, select=-c(Expression1)), 1)
####
#Divided by gene desert and gene rich
#####
cat('AUROC for each individual feature.\nDivided by gene desert and gene rich area')


geneRich.ROC <- featureROClog[featureROClog$genomicDivision=='geneRich',]
ROC.geneRich <- ROCExpression(subset(geneRich.ROC, select=-c(genomicDivision)), 1)
#t(ROCExpressionDirection(subset(geneRich.ROC, select=-c(genomicDivision)), 1))
#t(ROCExpressionDirection(subset(geneDesert.ROC, select=-c(genomicDivision)), 1))
distdhsType <- distdhsTypeROC
distdhsType$Expression <- normdataiPCR[match(distdhsType$V1, normdataiPCR$BC),]$expression
distdhsType <- distdhsType[!is.na(distdhsType$Expression),]
dhsTypes <- list(c(2:6), c(7:11), c(12:16), c(17:21), c(22:26), c(7:16,22:26))

richDHStypeROC <- data.frame(matrix(nrow=1, ncol=6))
for (i in 1:length(dhsTypes)) {
       geneRichDistdhsType <- subset(distdhsType, select=c(1, dhsTypes[[i]], Expression)) %>%
       filter(V1 %in%geneRich.ROC$BC) %>%
       mutate(Expression = ifelse(Expression>as.numeric(summary(Expression)[5]), 1, ifelse(Expression<(as.numeric(summary(Expression)[2])), 0,NA))) %>%
       filter(!is.na(Expression))
       
       
       fit <- glm(formula=Expression ~. ,data=subset(geneRichDistdhsType, select=-c(V1)), family='binomial')
       preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
       perf <- performance(preds, "auc")
       richDHStypeROC[i] <-round(as.numeric(slot(perf, "y.values")), 3)
}

colnames(richDHStypeROC) <- c('allTogether_nearest5', 'CTCF_nearest5', 'Distal_neareast5', 'nonCTCF_nearest5', 'Promoter_nearest5', 'allSeparately_nearest5')
ROC.geneRich <- cbind(ROC.geneRich, richDHStypeROC)
classNumberRich <- table(geneRichDistdhsType$Expression)

names(ROC.geneRich) <- gsub('Density', 'Accessibility', names(ROC.geneRich)) 
names(ROC.geneRich) <- gsub('Dens', 'Accessibility', names(ROC.geneRich)) 
ROC.geneRich <- ROC.geneRich[,!names(ROC.geneRich) %in% c('DistToGeneDesert','DistTok562GeneDesert', 'geneLoc','nDHS300kbsum', 'meanDist5DHS' , 'distToDHScluster')]


ROC.geneRich <- ROC.geneRich %>%
                     melt()

ROC.geneRich$Feature <- 'NA'
ROC.geneRich[grep("^n.",ROC.geneRich$variable),]$Feature <- 'NumberOfFeatures'
ROC.geneRich[grep("Access|totalTags",ROC.geneRich$variable),]$Feature <- 'Accessibility'
ROC.geneRich[grep("H3|H4|HDAC|H2AF",ROC.geneRich$variable),]$Feature <- 'Histones'
ROC.geneRich[grep("Dist|dist|Nearest|nearest",ROC.geneRich$variable),]$Feature <- 'DistanceToNearest' 
ROC.geneRich[grep("^n|Dist|dist|Nearest|nearest|Access|totalTags|H3|H4|HDAC|H2AF",ROC.geneRich$variable, invert=T),]$Feature <- 'Other' 
ROC.geneRich.order <- ROC.geneRich[order(ROC.geneRich$value),]$variable

pdf(paste0(opt$OUTDIR,'/graphs/AUROC_SelectedFeatures_geneRich.pdf'), height=8, width=10)
ROC.geneRich %>%
       mutate(SampleName=opt$SampleName,
              variable= factor(variable, levels=ROC.geneRich.order)) %>%
       ggplot(aes(x=variable, y=value, fill=Feature)) +
              geom_bar(stat='identity', position=position_dodge(), colour='black', alpha=0.6) +
              xlab('Features') +
              ylab('AUROC') +
              #facet_wrap(~Feature)+
              #coord_cartesian(ylim = c(0.5, 0.9)) +
              scale_fill_brewer(palette='Set1') +
              theme_classic() +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(paste0(opt$SampleName, ' - geneRich'), subtitle=paste0('AUROC for each individual feature classifying top and bottom 25% insertions based on gene expression\n',
                      'top25% = ',as.numeric(classNumberRich)[1], ', ', 'bottom25% = ',as.numeric(classNumberRich)[2]))+
              coord_flip(ylim=c(0.5, 0.95))+
              geom_text(aes(label=value), hjust=-0.1)
             # theme(legend.position="none")
       
dev.off()



geneDesert.ROC <- featureROClog[featureROClog$genomicDivision=='geneDesert',]
ROC.geneDesert <- ROCExpression(subset(geneDesert.ROC, select=-c(genomicDivision)), 1)

desertDHStypeROC <- data.frame(matrix(nrow=1, ncol=6))
for (i in 1:length(dhsTypes)) {
       genedesertDistdhsType <- subset(distdhsType, select=c(1, dhsTypes[[i]], Expression)) %>%
       filter(V1 %in%geneDesert.ROC$BC) %>%
       mutate(Expression = ifelse(Expression>as.numeric(summary(Expression)[5]), 1, ifelse(Expression<(as.numeric(summary(Expression)[2])), 0,NA))) %>%
       filter(!is.na(Expression))
       
       
       fit <- glm(formula=Expression ~. ,data=subset(genedesertDistdhsType, select=-c(V1)), family='binomial')
       preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
       perf <- performance(preds, "auc")
       desertDHStypeROC[i] <-round(as.numeric(slot(perf, "y.values")), 3)
}
colnames(desertDHStypeROC) <- c('allTogether_nearest5', 'CTCF_nearest5', 'Distal_neareast5', 'nonCTCF_nearest5', 'Promoter_nearest5', 'allSeparately_nearest5')
ROC.geneDesert <- cbind(ROC.geneDesert, desertDHStypeROC)
classNumberDesert <- table(genedesertDistdhsType$Expression)
names(ROC.geneDesert) <- gsub('Density', 'Accessibility', names(ROC.geneDesert)) 
names(ROC.geneDesert) <- gsub('Dens', 'Accessibility', names(ROC.geneDesert)) 
ROC.geneDesert <- ROC.geneDesert[,!names(ROC.geneDesert) %in% c('DistToGeneDesert','DistTok562GeneDesert', 'geneLoc','nDHS300kbsum', 'meanDist5DHS' , 'distToDHScluster')]



ROC.geneDesert <- ROC.geneDesert %>%
                     melt()

ROC.geneDesert$Feature <- 'NA'
ROC.geneDesert[grep("^n.",ROC.geneDesert$variable),]$Feature <- 'NumberOfFeatures'
ROC.geneDesert[grep("Access|totalTags",ROC.geneDesert$variable),]$Feature <- 'Accessibility'
ROC.geneDesert[grep("H3|H4|HDAC|H2AF",ROC.geneDesert$variable),]$Feature <- 'Histones'
ROC.geneDesert[grep("Dist|dist|Nearest|nearest",ROC.geneDesert$variable),]$Feature <- 'DistanceToNearest' 
ROC.geneDesert[grep("^n|Dist|dist|Nearest|nearest|Access|totalTags|H3|H4|HDAC|H2AF",ROC.geneDesert$variable, invert=T),]$Feature <- 'Other' 
ROC.geneDesert.order <- ROC.geneDesert[order(ROC.geneDesert$value),]$variable

pdf(paste0(opt$OUTDIR,'/graphs/AUROC_SelectedFeatures_geneDesert.pdf'), height=8, width=10)
ROC.geneDesert %>%
       mutate(SampleName=opt$SampleName,
              variable= factor(variable, levels=ROC.geneRich.order)) %>%
       ggplot(aes(x=variable, y=value, fill=Feature)) +
              geom_bar(stat='identity', position=position_dodge(), colour='black', alpha=0.6) +
              xlab('Features') +
              ylab('AUROC') +
              #facet_wrap(~Feature)+
              #coord_cartesian(ylim = c(0.5, 0.9)) +
              scale_fill_brewer(palette='Set1') +
              theme_classic() +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(paste0(opt$SampleName, ' - geneDesert'), subtitle=paste0('AUROC for each individual feature classifying top and bottom 25% insertions based on gene expression\n',
              'top25% = ',as.numeric(classNumberDesert)[1], ', ', 'bottom25% = ',as.numeric(classNumberDesert)[2]))+
              coord_flip(ylim=c(0.5, 0.95))+
              geom_text(aes(label=value), hjust=-0.1)
             # theme(legend.position="none")
       
dev.off()


#####
#Compare expression for geneDesert and geneRich top and bottom 20%
#####
geneRichDistdhsType$Expression[geneRichDistdhsType$Expression==1] <-'top25'
geneRichDistdhsType$Expression[geneRichDistdhsType$Expression==0] <-'bottom25'
geneRichDistdhsType$geneRegion <- 'geneRich'
genedesertDistdhsType$Expression[genedesertDistdhsType$Expression==1] <-'top25'
genedesertDistdhsType$Expression[genedesertDistdhsType$Expression==0] <-'bottom25'
genedesertDistdhsType$geneRegion <- 'geneDesert'

desertRichExpression <- rbind(genedesertDistdhsType, geneRichDistdhsType) 
colnames(desertRichExpression)[17] <- 'expressionLevel'
desertRichExpression$Expression <- normdataiPCR[match(desertRichExpression$V1, normdataiPCR$BC),]$expression

pdf(paste0(opt$OUTDIR,'/graphs/AUROC_geneDesertRich_20per_expression.pdf'), height=5, width=6)
desertRichExpression %>%
       mutate(geneRegion = factor(geneRegion, levels=c('geneRich', 'geneDesert'))) %>%
       ggplot(aes(x=geneRegion, y=Expression, fill=expressionLevel)) +
       geom_boxplot()+
       theme_classic()+
       scale_fill_brewer(palette='Set2')+
       annotate("text", label=paste0('n = ',nrow(geneRichDistdhsType)), x=1, y=2)+
       annotate("text", label=paste0('n = ',nrow(genedesertDistdhsType)), x=2, y=2) +
       ggtitle(opt$SampleName, 'Expression of the top and bottom 25% of insertions\nin gene deserts or in gene rich ares') 

dev.off()



#######
#AUROC for gene rich ares and gene deserts based on ALL insertions top and bottom 20%
#######
cat('AUROC for each individual feature.\nDivided by gene desert and gene rich area\n')
cat('Creating AUROC for distance to DHS types\n')
ROCgenomicRegion <- function(data, BCcol, region) { 
       library(ROCR)
       library(dplyr)
       library(reshape2)
       data <- data[data[,BCcol] %in% normdataiPCR$BC,]
       data$expression <- normdataiPCR[match(data[,BCcol],normdataiPCR$BC),]$expression
       dataROC <- subset(data, select=-c(BCcol))
       ROCdf <- as.data.frame(matrix(ncol=(ncol(dataROC)-1),nrow=1))
       rownames(ROCdf) <- 'AUROC'
       for (i in 1:c(ncol(dataROC)-1)) {
              #cat(colnames(dataROC)[i],'\n')
              if (colnames(dataROC)[i] != 'genomicDivision') {
                     dataROCglm <- dataROC %>%
                     mutate(expression = ifelse(expression>as.numeric(summary(dataROC$expression)[5]), 1, ifelse(expression<(as.numeric(summary(dataROC$expression)[2])), 0,NA)))
                     
                     dataROCglm <- subset(dataROCglm, select=c(colnames(dataROCglm)[i], 'expression', 'genomicDivision'))
                    
                     dataROCglm <- dataROCglm[!is.na(dataROCglm$expression),]
                     dataROCglm <- dataROCglm[dataROCglm$genomicDivision==region,]
                     Expression <- dataROCglm$expression
                     dataROCglm <- subset(dataROCglm, select=-c(expression, genomicDivision))
                     fit2 <- glm(formula=Expression ~. ,data=dataROCglm, family='binomial')
              
                     preds <- prediction(as.numeric(fit2$fitted.values), as.numeric(fit2$y))
                     perf <- performance(preds, "auc")
                     rocC <- round(as.numeric(slot(perf, "y.values")), 3)
                     #cat(rocC, '\n')
                     colnames(ROCdf)[i] <- colnames(dataROC)[i]
                     ROCdf[,i] <- rocC
              }
              
       }
       return(ROCdf)
}

ROC.geneRich <- ROCgenomicRegion(featureROClog, 1, 'geneRich')
geneRich.ROC <- featureROClog[featureROClog$genomicDivision=='geneRich',]
distdhsType <- distdhsTypeROC
distdhsType$Expression <- normdataiPCR[match(distdhsType$V1, normdataiPCR$BC),]$expression
distdhsType <- distdhsType[!is.na(distdhsType$Expression),]
dhsTypes <- list(c(2:6), c(7:11), c(12:16), c(17:21), c(22:26), c(7:16,22:26))

richDHStypeROC <- data.frame(matrix(nrow=1, ncol=6))
for (i in 1:length(dhsTypes)) {
       geneRichDistdhsType <- subset(distdhsType, select=c(1, dhsTypes[[i]], Expression)) %>%
       mutate(Expression = ifelse(Expression>as.numeric(summary(Expression)[5]), 1, ifelse(Expression<(as.numeric(summary(Expression)[2])), 0,NA))) %>%
       filter(!is.na(Expression)) %>%
       filter(V1 %in%geneRich.ROC$BC)
       
       
       fit <- glm(formula=Expression ~. ,data=subset(geneRichDistdhsType, select=-c(V1)), family='binomial')
       preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
       perf <- performance(preds, "auc")
       richDHStypeROC[i] <-round(as.numeric(slot(perf, "y.values")), 3)
}

colnames(richDHStypeROC) <- c('allTogether_nearest5', 'CTCF_nearest5', 'Distal_neareast5', 'nonCTCF_nearest5', 'Promoter_nearest5', 'allSeparately_nearest5')
ROC.geneRich <- cbind(ROC.geneRich, richDHStypeROC)
classNumberRich <- table(geneRichDistdhsType$Expression)
names(ROC.geneRich) <- gsub('Density', 'Accessibility', names(ROC.geneRich)) 
names(ROC.geneRich) <- gsub('Dens', 'Accessibility', names(ROC.geneRich)) 
ROC.geneRich <- ROC.geneRich[,!names(ROC.geneRich) %in% c('DistToGeneDesert','DistTok562GeneDesert', 'geneLoc','nDHS300kbsum', 'meanDist5DHS' , 'distToDHScluster', 'V25')]


ROC.geneRich <- ROC.geneRich %>%
                     melt()

ROC.geneRich$Feature <- 'NA'
ROC.geneRich[grep("^n.",ROC.geneRich$variable),]$Feature <- 'NumberOfFeatures'
ROC.geneRich[grep("Access|totalTags",ROC.geneRich$variable),]$Feature <- 'Accessibility'
ROC.geneRich[grep("H3|H4|HDAC|H2AF",ROC.geneRich$variable),]$Feature <- 'Histones'
ROC.geneRich[grep("Dist|dist|Nearest|nearest",ROC.geneRich$variable),]$Feature <- 'DistanceToNearest' 
ROC.geneRich[grep("^n|Dist|dist|Nearest|nearest|Access|totalTags|H3|H4|HDAC|H2AF",ROC.geneRich$variable, invert=T),]$Feature <- 'Other' 
ROC.geneRich.order <- ROC.geneRich[order(ROC.geneRich$value),]$variable

pdf(paste0(opt$OUTDIR,'/graphs/AUROC_SelectedFeatures_geneRich_totalTopBottom.pdf'), height=8, width=10)
ROC.geneRich %>%
       mutate(SampleName=opt$SampleName,
              variable= factor(variable, levels=ROC.geneRich.order)) %>%
       ggplot(aes(x=variable, y=value, fill=Feature)) +
              geom_bar(stat='identity', position=position_dodge(), colour='black', alpha=0.6) +
              xlab('Features') +
              ylab('AUROC') +
              #facet_wrap(~Feature)+
              #coord_cartesian(ylim = c(0.5, 0.9)) +
              scale_fill_brewer(palette='Set1') +
              theme_classic() +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(paste0(opt$SampleName, ' - geneRich'), subtitle=paste0('AUROC for each individual feature classifying all top and bottom 25% insertions based on gene expression\n',
                      'top25% = ',as.numeric(classNumberRich)[1], ', ', 'bottom25% = ',as.numeric(classNumberRich)[2]))+
              coord_flip(ylim=c(0.5, 0.95))+
              geom_text(aes(label=value), hjust=-0.1)
             # theme(legend.position="none")
       
dev.off()



geneDesert.ROC <- featureROClog[featureROClog$genomicDivision=='geneDesert',]
ROC.geneDesert <- ROCgenomicRegion(featureROClog, 1, 'geneDesert')
desertDHStypeROC <- data.frame(matrix(nrow=1, ncol=6))
for (i in 1:length(dhsTypes)) {
       genedesertDistdhsType <- subset(distdhsType, select=c(1, dhsTypes[[i]], Expression)) %>%
       mutate(Expression = ifelse(Expression>as.numeric(summary(Expression)[5]), 1, ifelse(Expression<(as.numeric(summary(Expression)[2])), 0,NA))) %>%
       filter(!is.na(Expression)) %>%
       filter(V1 %in%geneDesert.ROC$BC)
       
       
       fit <- glm(formula=Expression ~. ,data=subset(genedesertDistdhsType, select=-c(V1)), family='binomial')
       preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
       perf <- performance(preds, "auc")
       desertDHStypeROC[i] <-round(as.numeric(slot(perf, "y.values")), 3)
}
colnames(desertDHStypeROC) <- c('allTogether_nearest5', 'CTCF_nearest5', 'Distal_neareast5', 'nonCTCF_nearest5', 'Promoter_nearest5', 'allSeparately_nearest5')
ROC.geneDesert <- cbind(ROC.geneDesert, desertDHStypeROC)
classNumberDesert <- table(genedesertDistdhsType$Expression)

names(ROC.geneDesert) <- gsub('Density', 'Accessibility', names(ROC.geneDesert)) 
names(ROC.geneDesert) <- gsub('Dens', 'Accessibility', names(ROC.geneDesert)) 
ROC.geneDesert <- ROC.geneDesert[,!names(ROC.geneDesert) %in% c('DistToGeneDesert','DistTok562GeneDesert', 'geneLoc','nDHS300kbsum', 'meanDist5DHS' , 'distToDHScluster', 'V25')]



ROC.geneDesert <- ROC.geneDesert %>%
                     melt()

ROC.geneDesert$Feature <- 'NA'
ROC.geneDesert[grep("^n.",ROC.geneDesert$variable),]$Feature <- 'NumberOfFeatures'
ROC.geneDesert[grep("Access|totalTags",ROC.geneDesert$variable),]$Feature <- 'Accessibility'
ROC.geneDesert[grep("H3|H4|HDAC|H2AF",ROC.geneDesert$variable),]$Feature <- 'Histones'
ROC.geneDesert[grep("Dist|dist|Nearest|nearest",ROC.geneDesert$variable),]$Feature <- 'DistanceToNearest' 
ROC.geneDesert[grep("^n|Dist|dist|Nearest|nearest|Access|totalTags|H3|H4|HDAC|H2AF",ROC.geneDesert$variable, invert=T),]$Feature <- 'Other' 
ROC.geneDesert.order <- ROC.geneDesert[order(ROC.geneDesert$value),]$variable

pdf(paste0(opt$OUTDIR,'/graphs/AUROC_SelectedFeatures_geneDesert_totalTopBottom.pdf'), height=8, width=10)
ROC.geneDesert %>%
       mutate(SampleName=opt$SampleName,
              variable= factor(variable, levels=ROC.geneRich.order)) %>%
       ggplot(aes(x=variable, y=value, fill=Feature)) +
              geom_bar(stat='identity', position=position_dodge(), colour='black', alpha=0.6) +
              xlab('Features') +
              ylab('AUROC') +
              #facet_wrap(~Feature)+
              #coord_cartesian(ylim = c(0.5, 0.9)) +
              scale_fill_brewer(palette='Set1') +
              theme_classic() +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(paste0(opt$SampleName, ' - geneDesert'), subtitle=paste0('AUROC for each individual feature classifying all top and bottom 25% insertions based on gene expression\n',
                      'top25% = ',as.numeric(classNumberDesert)[1], ', ', 'bottom25% = ',as.numeric(classNumberDesert)[2]))+
              coord_flip(ylim=c(0.5, 0.95))+
              geom_text(aes(label=value), hjust=-0.1)
             # theme(legend.position="none")
       
dev.off()


######
#Gene rich/ gene desert (only those two divisions)
######

geneRichDesertNoExpROC <- subset(normdataiPCR,select=c(BC, DNA, RNA, DistToTSS, DistToNearestDHS, DistToNearestCTCF, DistToCGI, PercentGC, AB_compartment_value,
 chromHMM, GenicLocation.simple, nDHS100kb, nDHS10kb, nDHS300kb, nDHS500kb, nDHS1Mb, nDHS100kbDensity, nDHS300kbsum, nDHS500kbDensity, nDHS1MbDensity, nTSS100kb, 
 repliSeq, nExpressed100kb, nCGI100kb, nCTCFpeaks100kb, DNasePksDens, NearestDHStype, CTCFPksDens, genomicDivision, genomicDivisionNoExp, DistToGeneDesert, DistTok562GeneDesert, meanDist5DHS, distToDHScluster, 
 H2AFZ, H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me2, H3K4me3, H3K79me2, H3K9ac, H3K9me1, H3K9me3, H4K20me1, HDAC1, HDAC2, HDAC6, totalTags50kb ))
 
geneRichDesertNoExpROC$genomicDivision <- geneRichDesertNoExpROC$genomicDivisionNoExp
geneRichDesertNoExpROC <- subset(geneRichDesertNoExpROC, select=-c(genomicDivisionNoExp))
geneRichDesertNoExpROC.clean <- geneRichDesertNoExpROC %>%
       filter(!is.na(repliSeq),
              !is.na(chromHMM),
              !is.na(NearestDHStype),
              !is.na(genomicDivision),
              !is.na(GenicLocation.simple),
              !is.na(AB_compartment_value)) %>%
       mutate(expression = ifelse(log10(RNA/DNA) >as.numeric(summary(normdataiPCR$expression)[5]), 1, ifelse(log10(RNA/DNA) <(as.numeric(summary(normdataiPCR$expression)[2])), 0,NA)), 
              DistToTSS = abs(DistToTSS),
              chromHMM = factor(chromHMM),
              geneLoc = factor(GenicLocation.simple),
              NearestDHStype = factor(NearestDHStype), 
              repliSeq = factor(repliSeq),
              genomicDivision = factor(genomicDivision),
              DistToNearestDHS = abs(DistToNearestDHS),
              DistToCGI = abs(DistToCGI),
              DistToNearestCTCF = abs(DistToNearestCTCF),
              distToDHScluster =abs(distToDHScluster),
              H2AFZ= abs(H2AFZ),
              H3K27ac = abs(H3K27ac), 
              H3K27me3 = abs(H3K27me3), 
              H3K36me3 = abs(H3K36me3), 
              H3K4me1 =abs(H3K4me1), 
              H3K4me2 = abs(H3K4me2), 
              H3K4me3 = abs(H3K4me3),
              H3K79me2 = abs(H3K79me2), 
              H3K9ac =abs(H3K9ac), 
              H3K9me1 =abs(H3K9me1), 
              H3K9me3 =abs(H3K9me3), 
              H4K20me1 =abs(H4K20me1), 
              HDAC1 =abs(HDAC1), 
              HDAC2 =abs(HDAC2), 
              HDAC6 =abs(HDAC6))

#geneRichDesertNoExpROC.clean <- geneRichDesertNoExpROC.clean[!is.na(geneRichDesertNoExpROC.clean$expression),]

geneRichDesertNoExpROClog <- geneRichDesertNoExpROC.clean 
geneRichDesertNoExpROClog <- geneRichDesertNoExpROClog[!is.na(geneRichDesertNoExpROClog$totalTags50kb),]
Expression <- geneRichDesertNoExpROClog$expression
geneRichDesertNoExpROClog <- subset(geneRichDesertNoExpROClog, select=-c(expression,  DNA, RNA, GenicLocation.simple))


logCol <- c("DistToTSS", "DistToNearestDHS", "DistToNearestCTCF", 
"DistToCGI",  
"nDHS100kbDensity", "nDHS300kbsum", "nDHS500kbDensity", "nDHS1MbDensity", 
"CTCFPksDens",
 "meanDist5DHS", "distToDHScluster", 
"H2AFZ", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me2", 
"H3K4me3", "H3K79me2", "H3K9ac", "H3K9me1", "H3K9me3", "H4K20me1", 
"HDAC1", "HDAC2", "HDAC6", "totalTags50kb")

geneRichDesertNoExpROClog[,logCol] <- log10(abs(geneRichDesertNoExpROClog[,logCol])+1)
ROC.geneRich <- ROCgenomicRegion(geneRichDesertNoExpROClog, 1, 'geneRich')
geneRich.ROC <- geneRichDesertNoExpROClog[geneRichDesertNoExpROClog$genomicDivision=='geneRich',]
distdhsType <- distdhsTypeROC
distdhsType$Expression <- normdataiPCR[match(distdhsType$V1, normdataiPCR$BC),]$expression
distdhsType <- distdhsType[!is.na(distdhsType$Expression),]
dhsTypes <- list(c(2:6), c(7:11), c(12:16), c(17:21), c(22:26), c(7:16,22:26))

richDHStypeROC <- data.frame(matrix(nrow=1, ncol=6))
for (i in 1:length(dhsTypes)) {
       geneRichDistdhsType <- subset(distdhsType, select=c(1, dhsTypes[[i]], Expression)) %>%
       mutate(Expression = ifelse(Expression>as.numeric(summary(Expression)[5]), 1, ifelse(Expression<(as.numeric(summary(Expression)[2])), 0,NA))) %>%
       filter(!is.na(Expression)) %>%
       filter(V1 %in%geneRich.ROC$BC)
       
       
       fit <- glm(formula=Expression ~. ,data=subset(geneRichDistdhsType, select=-c(V1)), family='binomial')
       preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
       perf <- performance(preds, "auc")
       richDHStypeROC[i] <-round(as.numeric(slot(perf, "y.values")), 3)
}

colnames(richDHStypeROC) <- c('allTogether_nearest5', 'CTCF_nearest5', 'Distal_neareast5', 'nonCTCF_nearest5', 'Promoter_nearest5', 'allSeparately_nearest5')
ROC.geneRich <- cbind(ROC.geneRich, richDHStypeROC)
classNumberRich <- table(geneRichDistdhsType$Expression)
names(ROC.geneRich) <- gsub('Density', 'Accessibility', names(ROC.geneRich)) 
names(ROC.geneRich) <- gsub('Dens', 'Accessibility', names(ROC.geneRich)) 
ROC.geneRich <- ROC.geneRich[,!names(ROC.geneRich) %in% c('DistToGeneDesert','DistTok562GeneDesert', 'geneLoc','nDHS300kbsum', 'meanDist5DHS' , 'distToDHScluster', 'V25')]


ROC.geneRich <- ROC.geneRich %>%
                     melt()

ROC.geneRich$Feature <- 'NA'
ROC.geneRich[grep("^n.",ROC.geneRich$variable),]$Feature <- 'NumberOfFeatures'
ROC.geneRich[grep("Access|totalTags",ROC.geneRich$variable),]$Feature <- 'Accessibility'
ROC.geneRich[grep("H3|H4|HDAC|H2AF",ROC.geneRich$variable),]$Feature <- 'Histones'
ROC.geneRich[grep("Dist|dist|Nearest|nearest",ROC.geneRich$variable),]$Feature <- 'DistanceToNearest' 
ROC.geneRich[grep("^n|Dist|dist|Nearest|nearest|Access|totalTags|H3|H4|HDAC|H2AF",ROC.geneRich$variable, invert=T),]$Feature <- 'Other' 
ROC.geneRich.order <- ROC.geneRich[order(ROC.geneRich$value),]$variable

pdf(paste0(opt$OUTDIR,'/graphs/AUROC_SelectedFeatures_geneRich_totalTopBottom_noExp.pdf'), height=8, width=10)
ROC.geneRich %>%
       mutate(SampleName=opt$SampleName,
              variable= factor(variable, levels=ROC.geneRich.order)) %>%
       ggplot(aes(x=variable, y=value, fill=Feature)) +
              geom_bar(stat='identity', position=position_dodge(), colour='black', alpha=0.6) +
              xlab('Features') +
              ylab('AUROC') +
              #facet_wrap(~Feature)+
              #coord_cartesian(ylim = c(0.5, 0.9)) +
              scale_fill_brewer(palette='Set1') +
              theme_classic() +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(paste0(opt$SampleName, ' - geneRich'), subtitle=paste0('AUROC for each individual feature classifying all top and bottom 25% insertions based on gene expression\n',
                      'top25% = ',as.numeric(classNumberRich)[1], ', ', 'bottom25% = ',as.numeric(classNumberRich)[2]))+
              coord_flip(ylim=c(0.5, 0.95))+
              geom_text(aes(label=value), hjust=-0.1)
             # theme(legend.position="none")
       
dev.off()



geneDesert.ROC <- geneRichDesertNoExpROClog[geneRichDesertNoExpROClog$genomicDivision=='geneDesert',]
ROC.geneDesert <- ROCgenomicRegion(geneRichDesertNoExpROClog, 1, 'geneDesert')
desertDHStypeROC <- data.frame(matrix(nrow=1, ncol=6))
for (i in 1:length(dhsTypes)) {
       genedesertDistdhsType <- subset(distdhsType, select=c(1, dhsTypes[[i]], Expression)) %>%
       mutate(Expression = ifelse(Expression>as.numeric(summary(Expression)[5]), 1, ifelse(Expression<(as.numeric(summary(Expression)[2])), 0,NA))) %>%
       filter(!is.na(Expression)) %>%
       filter(V1 %in%geneDesert.ROC$BC)
       
       
       fit <- glm(formula=Expression ~. ,data=subset(genedesertDistdhsType, select=-c(V1)), family='binomial')
       preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
       perf <- performance(preds, "auc")
       desertDHStypeROC[i] <-round(as.numeric(slot(perf, "y.values")), 3)
}
colnames(desertDHStypeROC) <- c('allTogether_nearest5', 'CTCF_nearest5', 'Distal_neareast5', 'nonCTCF_nearest5', 'Promoter_nearest5', 'allSeparately_nearest5')
ROC.geneDesert <- cbind(ROC.geneDesert, desertDHStypeROC)
classNumberDesert <- table(genedesertDistdhsType$Expression)

names(ROC.geneDesert) <- gsub('Density', 'Accessibility', names(ROC.geneDesert)) 
names(ROC.geneDesert) <- gsub('Dens', 'Accessibility', names(ROC.geneDesert)) 
ROC.geneDesert <- ROC.geneDesert[,!names(ROC.geneDesert) %in% c('DistToGeneDesert','DistTok562GeneDesert', 'geneLoc','nDHS300kbsum', 'meanDist5DHS' , 'distToDHScluster', 'V25')]



ROC.geneDesert <- ROC.geneDesert %>%
                     melt()

ROC.geneDesert$Feature <- 'NA'
ROC.geneDesert[grep("^n.",ROC.geneDesert$variable),]$Feature <- 'NumberOfFeatures'
ROC.geneDesert[grep("Access|totalTags",ROC.geneDesert$variable),]$Feature <- 'Accessibility'
ROC.geneDesert[grep("H3|H4|HDAC|H2AF",ROC.geneDesert$variable),]$Feature <- 'Histones'
ROC.geneDesert[grep("Dist|dist|Nearest|nearest",ROC.geneDesert$variable),]$Feature <- 'DistanceToNearest' 
ROC.geneDesert[grep("^n|Dist|dist|Nearest|nearest|Access|totalTags|H3|H4|HDAC|H2AF",ROC.geneDesert$variable, invert=T),]$Feature <- 'Other' 
ROC.geneDesert.order <- ROC.geneDesert[order(ROC.geneDesert$value),]$variable

pdf(paste0(opt$OUTDIR,'/graphs/AUROC_SelectedFeatures_geneDesert_totalTopBottom_noExp.pdf'), height=8, width=10)
ROC.geneDesert %>%
       mutate(SampleName=opt$SampleName,
              variable= factor(variable, levels=ROC.geneRich.order)) %>%
       ggplot(aes(x=variable, y=value, fill=Feature)) +
              geom_bar(stat='identity', position=position_dodge(), colour='black', alpha=0.6) +
              xlab('Features') +
              ylab('AUROC') +
              #facet_wrap(~Feature)+
              #coord_cartesian(ylim = c(0.5, 0.9)) +
              scale_fill_brewer(palette='Set1') +
              theme_classic() +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(paste0(opt$SampleName, ' - geneDesert'), subtitle=paste0('AUROC for each individual feature classifying all top and bottom 25% insertions based on gene expression\n',
                      'top25% = ',as.numeric(classNumberDesert)[1], ', ', 'bottom25% = ',as.numeric(classNumberDesert)[2]))+
              coord_flip(ylim=c(0.5, 0.95))+
              geom_text(aes(label=value), hjust=-0.1)
             # theme(legend.position="none")
       
dev.off()




#####
#AUROC for DHS clusters
#####
#cat('AUROC for DHS clusters\n')
#
#dhsClust10kb <- read.table(paste0(opt$OUTDIR,'/',opt$SampleName,'_DHShotCluster10kb.bed'),sep='\t', header=F, stringsAsFactors=F) 
#dhsClust10kb$V3 <- log10(abs(dhsClust10kb$V3)+1)
#dhsClust10kb$Cluster <- rep(paste0('Cluster-10kb_',formatC(1:5, width = 2, format = "d", flag = "0")), length(unique(dhsClust10kb$V1)))
##dhsClust10kbROC <- dcast(data = dhsClust10kb,formula = V1~Cluster,fun.aggregate = sum,value.var = "V3")
#
#dhsClust50kb <- read.table(paste0(opt$OUTDIR,'/',opt$SampleName,'_DHShotCluster50kb.bed'),sep='\t', header=F, stringsAsFactors=F) 
#dhsClust50kb$V3 <- log10(abs(dhsClust50kb$V3)+1)
#dhsClust50kb$Cluster <- rep(paste0('Cluster-50kb_',formatC(1:5, width = 2, format = "d", flag = "0")), length(unique(dhsClust50kb$V1)))
##dhsClust50kbROC <- dcast(data = dhsClust50kb,formula = V1~Cluster,fun.aggregate = sum,value.var = "V3")
#
#dhsClust100kb <- read.table(paste0(opt$OUTDIR,'/',opt$SampleName,'_DHShotCluster100kb.bed'),sep='\t', header=F, stringsAsFactors=F) 
#dhsClust100kb$V3 <- log10(abs(dhsClust100kb$V3)+1)
#dhsClust100kb$Cluster <- rep(paste0('Cluster-100kb_',formatC(1:5, width = 2, format = "d", flag = "0")), length(unique(dhsClust100kb$V1)))
##dhsClust100kbROC <- dcast(data = dhsClust100kb,formula = V1~Cluster,fun.aggregate = sum,value.var = "V3")
#
#
#dhsCluster <- rbind(dhsClust10kb, dhsClust50kb, dhsClust100kb)
#dhsCluster <- dcast(data = dhsCluster,formula = V1~Cluster,fun.aggregate = sum,value.var = "V3")
#
#dhsCluster.ROC <- ROCExpression(dhsCluster, 1)
#
#
#clustdhsTypes <- list(c(2:6), c(7:11), c(12:16))
#
#clustDHStypeROC <- data.frame(matrix(nrow=1, ncol=3))
#dhsCluster$Expression <- normdataiPCR[match(dhsCluster$V1, normdataiPCR$BC),]$expression
#dhsCluster <- dhsCluster[!is.na(dhsCluster$Expression),]
#for (i in 1:length(clustdhsTypes)) {
#       geneClustDistdhsType <- subset(dhsCluster, select=c(1, clustdhsTypes[[i]], Expression)) %>%
#       mutate(Expression = ifelse(Expression>as.numeric(summary(Expression)[5]), 1, ifelse(Expression<(as.numeric(summary(Expression)[2])), 0,NA))) %>%
#       filter(!is.na(Expression))
#       
#       
#       fit <- glm(formula=Expression ~. ,data=subset(geneClustDistdhsType, select=-c(V1)), family='binomial')
#       preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
#       perf <- performance(preds, "auc")
#       clustDHStypeROC[i] <-round(as.numeric(slot(perf, "y.values")), 3)
#}
#
#colnames(clustDHStypeROC) <- c('Cluster-100kb_nearest05', 'Cluster-10kb_nearest05', 'Cluster-50kb_nearest05')
#dhsCluster.ROC <- cbind(dhsCluster.ROC, clustDHStypeROC)
#
#
#pdf(paste0(opt$OUTDIR,'/graphs/AUROC_DHSclusters.pdf'), height=5, width=5)
#dhsCluster.ROC  %>%
#       melt() %>%
#       mutate(Cluster = gsub('_.*' ,'' , variable),
#              variable=paste0('Cluster_', gsub('.*_' ,'' , variable)),
#              Cluster= factor(Cluster, levels=c('Cluster-10kb', 'Cluster-50kb', 'Cluster-100kb'))) %>%
#       ggplot(aes(x=variable, y=value, fill=Cluster)) +
#              geom_bar(stat='identity', position=position_dodge(), colour='black', alpha=0.6) +
#              xlab('DHS cluster neighbour') +
#              ylab('AUROC') +
#              coord_cartesian(ylim = c(0.5, 0.9)) +
#              scale_fill_brewer(palette='Set1') +
#              theme_classic() +
#              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, ,hjust=0.95,vjust=0.95), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
#              ggtitle(opt$SampleName, subtitle='AUROC for nearest 5 DHS clusters individually\nas well as combined') #+
#              #facet_zoom(y=SampleName=='CMV',zoom.size = 2, ylim=c(0,3))
#       
#dev.off()
#
#
#
#####
#Plot mean distance to DHS vs nDHS100kb
#####
pdf(paste0(opt$OUTDIR,'/graphs/nDHS100kb_vs_meanDist5DHS.pdf'), height=5, width=12)
normdataiPCR %>%
       mutate(expression = ifelse(expression>as.numeric(summary(expression)[5]), 'top25', ifelse(expression<(as.numeric(summary(expression)[2])), 'bottom25','middle50'))) %>%
       ggplot(aes(x=log10(abs(meanDist5DHS)+1), y=nDHS100kb, fill=expression, colour=expression)) +
              geom_point(pch=21, colour='black')+ 
              theme_classic()+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              scale_fill_brewer(palette='Set1') +
              scale_colour_brewer(palette='Set1') +
              geom_smooth(method = "loess") +
              facet_wrap(~expression)+
              xlab('log10(mean distance to nearest 5 DHS+1)') +
              ylab('number of DHS within 100kb')
dev.off()

bottomExpDHS <- normdataiPCR %>% filter(log10(abs(meanDist5DHS)+1) < 4, nDHS100kb>10, expression <(-0.3))
bottomExpDHS$ExpDHS <- 'bottom'
topExpDHS <- normdataiPCR %>% filter(log10(abs(meanDist5DHS)+1) < 4, nDHS100kb>10, expression > 0.3)
topExpDHS$ExpDHS <- 'top'

ExpDHS <- rbind(bottomExpDHS, topExpDHS)

pdf(paste0(opt$OUTDIR,'/graphs/topBottom_DNAdeciles_DHS.pdf'), height=5, width=10)
ExpDHS %>%
       ggplot(aes(x=DNA_Deciles)) +
       geom_bar(stat="count", fill=sampleColor, colour='black')+
       facet_wrap(~ExpDHS, scale='free') +
       theme_classic()+
       ggtitle(opt$SampleName, 'Comparison of DNA deciles between insertions having mean distance to 5 DHS <10kb &\nhaving more than 10 DHS within 100kb & and an expression of >0.3 (top) and <(-0.3) (bottom)')+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) 
dev.off()

#
#ExpDHS %>%
#       ggplot(aes(x=ExpDHS, y=meanDist5DHS)) +
#       geom_boxplot()+
#       theme_classic()+
#       ggtitle(opt$SampleName, 'Comparison of DNA deciles between insertions having mean distance to 5 DHS <10kb &\nhaving more than 10 DHS within 100kb & and an expression of >0.3 (top) and <(-0.3) (bottom)')+
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) 
####
#Number of insertions overlapping DHS, and histone marks
####
histDHSpks <- c("DistToNearestDHS", "H3K27ac", 
"H3K27me3", "H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K79me2", 
"H3K9ac", "H3K9me1", "H3K9me3", "H4K20me1")

histDHSpks <- normdataiPCR[,histDHSpks]
#histDHSpks <- histDHSpks[apply(histDHSpks, 1, function(r) any(r==0)),]
histDHSpks$BC <- normdataiPCR[rownames(histDHSpks),]$BC
histDHSpks$Expression <- normdataiPCR[rownames(histDHSpks),]$expression
#histDHSpks <- histDHSpks %>% 
#              mutate(ExpressionLevel = ifelse(Expression>as.numeric(summary(Expression)[5]), 'top25', ifelse(Expression<(as.numeric(summary(Expression)[2])), 'bottom25','middle50')))
#histDHSpks <- histDHSpks[histDHSpks$ExpressionLevel !="middle50",]
#par(mfrow=c(3, 4))
#for (i in 1:(ncol(histDHSpks)-2)){
#       cat(colnames(histDHSpks)[i], '\n')
#       reg<-lm(histDHSpks$Expression ~ log10(abs(histDHSpks[,i])+1))
#       cat(summary(reg)$adj.r.squared, '\n')
# plot(log10(abs(histDHSpks[,i])+1), histDHSpks$Expression, xlab='log10(distance+1)', main=paste0(colnames(histDHSpks)[i], ', adj.R^2 = ', round(summary(reg)$adj.r.squared, 3)), frame=F)
# abline(reg, col='red')
# }
 
 
 
#####
#Selected distance features vs expression)
#####
cat('Selected distance features vs expression\n')
geneDesert.ROC$Expression <- normdataiPCR[match(geneDesert.ROC$BC, normdataiPCR$BC),]$expression
geneRich.ROC$Expression <- normdataiPCR[match(geneRich.ROC$BC, normdataiPCR$BC),]$expression

pdf(paste0(opt$OUTDIR,'/graphs/DHSHistoneDistance_part1.pdf'), height=15, width=8) 
par(mfrow=c(6,2), mar=c(2,2,1,1.1)) 
featureDistace <- c("DistToNearestDHS",
 "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1","H3K4me2") 
for (i in featureDistace){
        cat(i, '\n')
       reg<-lm(geneRich.ROC$Expression ~ geneRich.ROC[,i])
       cat(summary(reg)$adj.r.squared, '\n')
      #smoothScatter(geneRich.ROC[,i], geneRich.ROC$Expression, xlab='log10(distance+1)', col='darkgreen', main=paste0(i, ', adj.R^2 = ', round(summary(reg)$adj.r.squared, 3)), frame=F)
       scatter.smooth(geneRich.ROC[,i], geneRich.ROC$Expression, xlab='log10(distance+1)', col='darkgreen', main=paste0(i, ', adj.R^2 = ', round(summary(reg)$adj.r.squared, 3)), frame=F)
       abline(reg, col='red')
       
       reg<-lm(geneDesert.ROC$Expression ~ geneDesert.ROC[,i])
       cat(summary(reg)$adj.r.squared, '\n')
       scatter.smooth(geneDesert.ROC[,i], geneDesert.ROC$Expression, xlab='log10(distance+1)' , col='steelblue', ylim=c(-2, 2), main=paste0(i, ', adj.R^2 = ', round(summary(reg)$adj.r.squared, 3)), frame=F)
       abline(reg, col='red')
}
dev.off()

pdf(paste0(opt$OUTDIR,'/graphs/DHSHistoneDistance_part2.pdf'), height=15, width=8) 
par(mfrow=c(6,2), mar=c(2,2,1,1.1)) 
featureDistace <- c("H3K4me3", "H3K79me2", "H3K9ac", "H3K9me1", "H3K9me3", "H4K20me1")
for (i in featureDistace){
        cat(i, '\n')
       reg<-lm(geneRich.ROC$Expression ~ geneRich.ROC[,i])
       cat(summary(reg)$adj.r.squared, '\n')
       scatter.smooth(geneRich.ROC[,i], geneRich.ROC$Expression, xlab='log10(distance+1)', col='darkgreen', main=paste0(i, ', adj.R^2 = ', round(summary(reg)$adj.r.squared, 3)), frame=F)
       abline(reg, col='red')
       
       reg<-lm(geneDesert.ROC$Expression ~ geneDesert.ROC[,i])
       cat(summary(reg)$adj.r.squared, '\n')
       scatter.smooth(geneDesert.ROC[,i], geneDesert.ROC$Expression, xlab='log10(distance+1)' , col='steelblue', ylim=c(-2, 2), main=paste0(i, ', adj.R^2 = ', round(summary(reg)$adj.r.squared, 3)), frame=F)
       abline(reg, col='red')
}
dev.off()

########
#Correlation between insertions and neighbouring insertions
########

cat('Correlation between neighbouring insertions\n')
neighIns <- read.table(paste0(opt$OUTDIR,'/',opt$SampleName,'_neighbours.bed'), sep='\t' , header=F, stringsAsFactors=F)
neighIns <- neighIns[neighIns$V1 %in% normdataiPCR$BC,]
neighIns <- neighIns[neighIns$V1!=neighIns$V2 ,]

nonOverlapIns <- as.character(as.data.frame(table(neighIns$V1))[as.data.frame(table(neighIns$V1))$Freq==20,]$Var1)
neighIns <- neighIns[neighIns$V1 %in% nonOverlapIns,]
neighIns$neighbour <- rep(paste0('Neighbour_',formatC(1:20, width = 2, format = "d", flag = "0")), length(unique(neighIns$V1)))
neighIns$Expression <- normdataiPCR[match(neighIns$V1, normdataiPCR$BC),]$expression
neighIns$NeighbourExpression <- normdataiPCR[match(neighIns$V2, normdataiPCR$BC),]$expression

neighInsDcast <- dcast(data = neighIns,formula = V1~neighbour,fun.aggregate = sum,value.var = "NeighbourExpression")
neighInsDcast$Expression <- normdataiPCR[match(neighInsDcast$V1, normdataiPCR$BC),]$expression
 
library(Hmisc)
corrNeighbours <- as.data.frame(matrix(nrow=1, ncol=20))
neighInsDcast <- subset(neighInsDcast, select=-c(V1))
for (i in 1:(ncol(neighInsDcast)-1)) {
        corrNeighbours[1,][,i] <- as.numeric(rcorr(neighInsDcast$Expression, neighInsDcast[,i], type="pearson")$r[,1][2])
       }

colnames(corrNeighbours) <- 1:20

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions.pdf'), height=5, width=5) 
corrNeighbours %>%
       melt() %>%
       mutate(SampleName=opt$SampleName) %>%
       ggplot(aes(x=variable, y=value, fill=SampleName)) +    
              geom_bar(stat='identity', col='black', alpha=0.6)+
              theme_classic() +
              scale_fill_manual(values=sampleColor) +
              ylab('Pearson correlation') +
              xlab('Neighbouring insertions') +
              ylim(c(0, 0.5)) +
              ggtitle(opt$SampleName, 'Correlation between insertion and its nearest 20 neighbouring insertions')+
              theme(legend.position='none')
dev.off()


######
#Binned distance All vs All pairwise comparison
######
cat('Binned distance All vs All pairwise comparison\n')

distList <- list(c(0, 5e4), c(5e4, 1e5), c(1e5, 2e5), c(2e5, 5e5), c(5e5, 1e6), c(1e6, 2e6), c(2e6, 1e9)) 
binDistCorrAllvsAll <- data.frame(matrix(ncol=3,nrow=length(distList))) 
colnames(binDistCorrAllvsAll) <- c('Pearsons', 'Dist', 'Pvalue')

for (i in 1:length(distList)) {
       binCorr <- rcorr( neighIns[neighIns$V3 >= distList[[i]][1] & neighIns$V3 < distList[[i]][2],]$Expression  , neighIns[neighIns$V3 >= distList[[i]][1] & neighIns$V3 < distList[[i]][2],]$NeighbourExpression)  
       binDistCorrAllvsAll[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
       binDistCorrAllvsAll[,2][i] <- paste0(distList[[i]][1], ' <= x < ', distList[[i]][2])
       binDistCorrAllvsAll[,3][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
}

binDistCorrAllvsAll$Dist <- gsub('.* x', '', binDistCorrAllvsAll$Dist)

binDistCorrAllvsAll$Dist[nrow(binDistCorrAllvsAll)] <- '> 2e+06'
pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_BinnedDistance_AllvsAll.pdf'), height=5, width=4) 
binDistCorrAllvsAll %>%
       mutate(Dist = factor(Dist, levels=unique(binDistCorrAllvsAll$Dist)),
              SampleName = opt$SampleName) %>%
       ggplot(aes(x = Dist, y = Pearsons, fill = SampleName)) +
              geom_bar(stat='identity', color='black', alpha=0.6) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_fill_manual(values=sampleColor) +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All vs all pairwise comparison of neighbours (< 20)')+
              xlab('Distance bins')+
              ylim(0,0.5) +
              theme(legend.position = 'none')
dev.off()

######
cat('Binned distance All vs All pairwise comparison 100kb bins\n')

distList <- list(c(0,  1e5), c(1e5, 2e5), c(2e5, 3e5), c(3e5, 4e5), c(4e5, 5e5), c(5e5, 6e5), c(6e5, 7e5), c(7e5, 8e5), c(8e5, 9e5), c(9e5, 1e6)) 
binDistCorrAllvsAll <- data.frame(matrix(ncol=3,nrow=length(distList))) 
colnames(binDistCorrAllvsAll) <- c('Pearsons', 'Dist', 'Pvalue')

for (i in 1:length(distList)) {
       binCorr <- rcorr( neighIns[neighIns$V3 >= distList[[i]][1] & neighIns$V3 < distList[[i]][2],]$Expression  , neighIns[neighIns$V3 >= distList[[i]][1] & neighIns$V3 < distList[[i]][2],]$NeighbourExpression)  
       binDistCorrAllvsAll[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
       binDistCorrAllvsAll[,2][i] <- paste0(distList[[i]][1], ' <= x < ', distList[[i]][2])
       binDistCorrAllvsAll[,3][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
}

binDistCorrAllvsAll$Dist <- gsub('.* x', '', binDistCorrAllvsAll$Dist)

#binDistCorrAllvsAll$Dist[nrow(binDistCorrAllvsAll)] <- '> 2e+06'
pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_BinnedDistance_AllvsAll_100kbBins.pdf'), height=5, width=4) 
binDistCorrAllvsAll %>%
       mutate(Dist = factor(Dist, levels=unique(binDistCorrAllvsAll$Dist)),
              SampleName = opt$SampleName) %>%
       ggplot(aes(x = Dist, y = Pearsons, colour = SampleName, group=SampleName)) +
              geom_line(size=1) +
              geom_point(size=4, pch=20) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_fill_manual(values=sampleColor) +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All vs all pairwise comparison of neighbours (< 20)')+
              xlab('Distance bins')+
              ylim(0,0.5) +
              theme(legend.position = 'none')
dev.off()


######
#Neighbouring correlation based on distance and genomic location
######
library(Hmisc)

neighLoc <- neighIns
neighLoc$InsLoc <- normdataiPCR[match(neighLoc$V1, normdataiPCR$BC),]$genomicDivisionNoExp
neighLoc$NeighLoc <- normdataiPCR[match(neighLoc$V2, normdataiPCR$BC),]$genomicDivisionNoExp


pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertionsDistance.pdf'), height=5, width=6) 
neighLoc %>%
       mutate(neighbour = gsub('Neighbour_','', neighbour),
              SampleName=opt$SampleName) %>%
       ggplot(aes(x=neighbour, y=log10(V3+1), fill=SampleName)) +
       geom_boxplot(outlier.shape=NA, alpha=0.6)+
       #facet_wrap(~neighbour) +
       #geom_vline(aes(xintercept=med, group=neighbour), linetype='dashed', color='black')+
       theme_classic()+
       scale_fill_manual(values=sampleColor)+
       xlab('Neighbour') +
       ylab('Distance to neighbour')+
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
       #geom_text(data=unique(subset(neighLoc, select=c(neighbour, med, SampleName))), aes(label=paste0('med log10(dist) = ',round(med,1)), group=neighbour), y=2000, x=3)+
       theme(legend.position = "none")+
       ggtitle(opt$SampleName, 'Distance to neighbours')
dev.off()



cat('Correlate gene region insertions\n')


corrMatDesert <- data.frame(matrix(ncol=4, nrow=length(unique(neighLoc$neighbour))))
colnames(corrMatDesert) <-c('Pearson', 'Pvalue','Neighbour', 'geneRegion')

corrMatRich <- data.frame(matrix(ncol=4, nrow=length(unique(neighLoc$neighbour))))
colnames(corrMatRich) <-c('Pearson', 'Pvalue','Neighbour', 'geneRegion')
for (neighbour in 1:length(unique(neighLoc$neighbour))) {
       corrDesert <- rcorr( neighLoc[neighLoc$InsLoc=='geneDesert' & neighLoc$NeighLoc=='geneDesert' & neighLoc$neighbour==unique(neighLoc$neighbour)[neighbour],]$Expression, 
              neighLoc[neighLoc$InsLoc=='geneDesert' & neighLoc$NeighLoc=='geneDesert' & neighLoc$neighbour==unique(neighLoc$neighbour)[neighbour],]$NeighbourExpression)
       
       corrRich <- rcorr( neighLoc[neighLoc$InsLoc=='geneRich' & neighLoc$NeighLoc=='geneRich' & neighLoc$neighbour==unique(neighLoc$neighbour)[neighbour],]$Expression, 
              neighLoc[neighLoc$InsLoc=='geneRich' & neighLoc$NeighLoc=='geneRich' & neighLoc$neighbour==unique(neighLoc$neighbour)[neighbour],]$NeighbourExpression)
       
       
       corrMatDesert[,1][neighbour] <- signif(corrDesert$r[1,2], 2)
       corrMatDesert[,2][neighbour] <- signif(corrDesert$P[1,2], 2)
       corrMatDesert[,3][neighbour] <- unique(neighLoc$neighbour)[neighbour]
       corrMatDesert[,4][neighbour] <-'geneDesert'
       
       corrMatRich[,1][neighbour] <- signif(corrRich$r[1,2], 2)
       corrMatRich[,2][neighbour] <- signif(corrRich$P[1,2], 2)
       corrMatRich[,3][neighbour] <- unique(neighLoc$neighbour)[neighbour]
       corrMatRich[,4][neighbour]<-'geneRich'
}
corrRichDesert <- rbind(corrMatDesert, corrMatRich)

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertionsDistance_desertRich.pdf'), height=5, width=7) 
corrRichDesert %>%
       mutate(geneRegion = factor(geneRegion, levels=c('geneRich', 'geneDesert')),
              Neighbour=gsub('Neighbour_', '',Neighbour)) %>%
              ggplot(aes(x=Neighbour, y=Pearson, fill=geneRegion)) +
              geom_bar(stat='identity', position=position_dodge()) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_fill_manual(values=c("#beaed4","#fdc086")) +
              ggtitle('Correlation between neighbours in gene deserts and gene rich area')
dev.off()

####
#First neighbour
####
cat('Correlate Gene desert/genic neighbours\n')

CorrNeighDesert <- neighLoc[neighLoc$neighbour=='Neighbour_01',]
CorrNeighDesert <- as.data.frame(CorrNeighDesert)
CorrNeighDesert <- CorrNeighDesert[!is.na(CorrNeighDesert$InsLoc) & !is.na( CorrNeighDesert$NeighLoc),]
corrMatRich <- data.frame(matrix(ncol=4, nrow=3))
colnames(corrMatRich) <-c('Pearson', 'Pvalue', 'geneRegion')
corrDesert <- rcorr( CorrNeighDesert[CorrNeighDesert$InsLoc=='geneDesert' & CorrNeighDesert$NeighLoc=='geneDesert',]$Expression, 
       CorrNeighDesert[CorrNeighDesert$InsLoc=='geneDesert' & CorrNeighDesert$NeighLoc=='geneDesert',]$NeighbourExpression)

corrRich <- rcorr( CorrNeighDesert[CorrNeighDesert$InsLoc=='geneRich' & CorrNeighDesert$NeighLoc=='geneRich' ,]$Expression, 
       CorrNeighDesert[CorrNeighDesert$InsLoc=='geneRich' & CorrNeighDesert$NeighLoc=='geneRich' ,]$NeighbourExpression)

corrDesertRich  <- rcorr( CorrNeighDesert[CorrNeighDesert$InsLoc!=CorrNeighDesert$NeighLoc,]$Expression, 
       CorrNeighDesert[CorrNeighDesert$InsLoc!=CorrNeighDesert$NeighLoc,]$NeighbourExpression)



corrDesertNeigh <- rbind(data.frame(Pearsons=signif(corrDesert$r[1,2], 2),
              Pvalue= signif(corrDesert$P[1,2], 2),
              geneRegion = 'DesertDesert'),
       data.frame(Pearsons=signif(corrRich$r[1,2], 2),
              Pvalue= signif(corrRich$P[1,2], 2),
              geneRegion = 'RichRich'),
       data.frame(Pearsons=signif(corrDesertRich$r[1,2], 2),
              Pvalue= signif(corrDesertRich$P[1,2], 2),
              geneRegion = 'RichDesert'))


rects <- data.frame(xstart = c(0.3,2), xend = c(2,3.7), geneRegion=c('RichRich', 'DesertDesert'))


pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertionsDistance_desertRichRegion.pdf'), height=5, width=4) 
corrDesertNeigh %>%
       mutate(SampleName = opt$SampleName,
              geneRegion = factor(geneRegion, levels=c('RichRich', 'RichDesert', 'DesertDesert'))) %>%
        ggplot() +
              geom_bar(stat='identity', color='black', alpha=0.6, aes( x = geneRegion, y = Pearsons), fill=sampleColor) +
              geom_rect(data=rects, mapping=aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill=geneRegion), alpha=0.25)+
              #geom_bar(stat='identity', color='black', alpha=0.6, aes( x = Compartment, y = Pearsons), fill=sampleColor)+
              theme_classic() +
              scale_fill_manual(values=c("#fdc086","#beaed4")) +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'Correlation between nearest\nneighbour and genomic region')+
              theme(legend.position='none')+
              ylab('Pearsons correlation') +
              xlab('Genomic region of nearest neighbours')+
              ylim(0, 0.5)+
              annotate('text', label='geneRich', x='RichRich', y=0.4, size=6)+ 
              annotate('text', label='geneDesert', x='DesertDesert', y=0.4, size=6)
dev.off()


#####
#Distance between genomic regions
#####
cat('Distance between gene desert/genic insertions\n')

CorrNeighDesert$geneRegion <- paste0(CorrNeighDesert$InsLoc, CorrNeighDesert$NieghLoc)
CorrNeighDesert[CorrNeighDesert$InsLoc!=CorrNeighDesert$NeighLoc,]$geneRegion <- 'RichDesert'
CorrNeighDesert[CorrNeighDesert$InsLoc=='geneDesert' & CorrNeighDesert$NeighLoc=='geneDesert',]$geneRegion <- 'DesertDesert'
CorrNeighDesert[CorrNeighDesert$InsLoc=='geneRich' & CorrNeighDesert$NeighLoc=='geneRich',]$geneRegion <- 'RichRich'


pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertionsDistance_desertRichRegionDistance.pdf'), height=5, width=4) 
CorrNeighDesert %>%
       mutate(geneRegion = factor(geneRegion, levels=c('RichRich', 'RichDesert', 'DesertDesert'))) %>%
       ggplot() +
              geom_boxplot(outlier.shape=NA, alpha=0.6, fill=sampleColor, aes(x=geneRegion, y=log10(V3+1))) +
              geom_rect(data=rects, mapping=aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill=geneRegion), alpha=0.25)+
              theme_classic() +
              scale_fill_manual(values=c("#fdc086","#beaed4")) +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'Distance between neighbours\nand genomic region')+
              theme(legend.position='none')+
              ylab('Distance between neighbours') +
              xlab('Genomic region of nearest neighbours')+
              annotate('text', label='geneRich', x='RichRich', y=2.5, size=6)+ 
              annotate('text', label='geneDesert', x='DesertDesert', y=2.5, size=6)
              
dev.off()



#####
#Expression and distance to genomic feature
#####
cat('Distance between Gene desert/genic\n')

distDivision<- normdataiPCR
#colnames(distDivision) <- c('a', 'b', 'c', 'BC', '4', '5', 'genomicDivisionNoExp','genomicDivision_desertDist', 'genomicDivision_richDist')

distDivision <- distDivision[!is.na(distDivision$genomicDivisionNoExp),]
distDivision$genomicDivision_desertDist <- abs(distDivision$genomicDivision_desertDist) *-1
distDivision$genomicDivision_richDist <- abs(distDivision$genomicDivision_richDist)
distDivision$genomicDivision_richDist[distDivision$genomicDivision_richDist==0] <- distDivision$genomicDivision_desertDist[distDivision$genomicDivision_richDist==0]

distDivision <- subset(distDivision, select=c(BC, genomicDivisionNoExp, genomicDivision_richDist))
colnames(distDivision) <- c('BC', 'genomicRegion', 'Distance')
distDivision$Expression <- normdataiPCR[match(distDivision$BC, normdataiPCR$BC),]$expression

rects <- data.frame(xstart = c(-5e5,0), xend = c(0,5e5), genomicRegion=c('geneRich', 'geneDesert'))

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_desertRichRegionExpression.pdf'), height=5, width=5) 
distDivision %>%
       filter(Distance <5e5 & Distance> -5e5) %>%
       mutate(genomicRegion = factor(genomicRegion, levels=c('geneRich', 'geneDesert'))) %>%
       ggplot()+
              #geom_point()+
              geom_smooth(aes(x=Distance, y=Expression), method='loess', color=sampleColor) +
              geom_rect(data=rects, mapping=aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill=genomicRegion), alpha=0.25)+
              #geom_bar(stat='identity', color='black', alpha=0.6, aes( x = Compartment, y = Pearsons), fill=sampleColor)+
              theme_classic() +
              scale_fill_manual(values=c("#fdc086","#beaed4")) +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'Compartment expression and distance\nto genomic region boundary')+
              theme(legend.position='none')+
              ylab('Reporter activity') +
              xlab('Distance to genomic region boundary')+
              annotate('text', label='geneRich', x=-2.5e5, y=-0.05, size=6)+ 
              annotate('text', label='geneDesert', x=2.5e5, y=-0.05, size=6)
dev.off()




######
#Correlation within and outside TADs
######
cat('Correlate TADs\n')
####
#TODO create loop for different TADs
#CorrTAD <- neighLoc[neighLoc$neighbour=='Neighbour_01',]

CorrTAD <- neighLoc

TADlist <- list(c('TAD_armatus0.2', 'Armatus_5kb_'), c('TAD_armatus0.5', 'Armatus_5kb_'), c('TAD_armatus0.8', 'Armatus_5kb_'), c('Arrowhead', 'Arrowhead_'))

for (TAD in 1:length(TADlist)) {
       cat('analysing ',TADlist[[TAD]][1], '\n')
       CorrTAD$TAD <- as.numeric(gsub(TADlist[[TAD]][2],'',normdataiPCR[match(CorrTAD$V1, normdataiPCR$BC),][,TADlist[[TAD]][1]]))
       CorrTAD$neighTAD <-  as.numeric(gsub(TADlist[[TAD]][2],'',normdataiPCR[match(CorrTAD$V2, normdataiPCR$BC),][,TADlist[[TAD]][1]]))
       CorrTAD <- CorrTAD[abs(CorrTAD$TAD-CorrTAD$neighTAD) <=1,]
       
       CorrTADsame <- rcorr( CorrTAD[CorrTAD$TAD == CorrTAD$neighTAD,]$Expression, 
               CorrTAD[CorrTAD$TAD == CorrTAD$neighTAD,]$NeighbourExpression)
       
       CorrTADdiff <- rcorr( CorrTAD[CorrTAD$TAD !=CorrTAD$neighTAD,]$Expression, 
               CorrTAD[CorrTAD$TAD != CorrTAD$neighTAD,]$NeighbourExpression)
       

       CorrTADcorr <- rbind(data.frame(Pearsons=signif(CorrTADsame$r[1,2], 2),
                     Pvalue= signif(CorrTADsame$P[1,2], 2),
                     TAD = 'Same'),
              data.frame(Pearsons=signif(CorrTADdiff$r[1,2], 2),
                     Pvalue= signif(CorrTADdiff$P[1,2], 2),
                     TAD = 'Different'))
       
       
      
       ggNeighbouringInsertions_TAD <- CorrTADcorr %>%
              mutate(SampleName = opt$SampleName,
                     TAD = factor(TAD, levels=c('Same', 'Different'))) %>%
               ggplot() +
                     geom_bar(stat='identity', color='black', alpha=0.6, aes( x = TAD, y = Pearsons), fill=sampleColor) +
                     #geom_bar(stat='identity', color='black', alpha=0.6, aes( x = Compartment, y = Pearsons), fill=sampleColor)+
                     theme_classic() +
                     scale_fill_manual(values=sampleColor) +
                     theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
                     ggtitle(opt$SampleName, paste0('Correlation between nearest\nneighbour in same or different TAD\n',TADlist[[TAD]][1]))+
                     theme(legend.position='none')+
                     ylab('Pearsons correlation') +
                     xlab('TAD of nearest neighbour')+
                     ylim(0, 0.5)
        ggsave( paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_TAD_',TADlist[[TAD]][1],'.pdf'), plot=ggNeighbouringInsertions_TAD ,height=5, width=3) 
       
      #pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_TAD_Distance',TADlist[[TAD]][1],'.pdf'), height=6, width=5) 
      #boxplot(log10(CorrTAD[CorrTAD$TAD == CorrTAD$neighTAD,]$V3+1), log10(CorrTAD[CorrTAD$TAD != CorrTAD$neighTAD,]$V3+1), 
      #names=c('sameTAD', 'diffTAD'), frame=F, xlab='TAD fo neighbour', ylab='log10(Distance to neighbour+1)', main=paste0(opt$SampleName,'\nDistance between neighbours\nin same or different TAD\n',TADlist[[TAD]][1]), col=sampleColor)
      #dev.off()
       
       
       
       #####
       #Correlation for different distances between same and different TAD
       #####
       #CorrTADsame
       #CorrTADdiff
       
       distList <- list(c(0,  1e5), c(1e5, 2e5), c(2e5, 3e5), c(3e5, 4e5), c(4e5, 5e5), c(5e5, 6e5), c(6e5, 7e5), c(7e5, 8e5), c(8e5, 9e5), c(9e5, 1e6)) 
       binDistTADsame <- data.frame(matrix(ncol=4,nrow=length(distList))) 
       colnames(binDistTADsame) <- c('Pearsons', 'Dist', 'Pvalue', 'n')
       for (i in 1:length(distList)) {
              binCorr <- rcorr( CorrTAD[CorrTAD$V3 >= distList[[i]][1] & CorrTAD$V3 < distList[[i]][2] & CorrTAD$TAD==CorrTAD$neighTAD,]$Expression  , CorrTAD[CorrTAD$V3 >= distList[[i]][1] & CorrTAD$V3 < distList[[i]][2] & CorrTAD$TAD==CorrTAD$neighTAD,]$NeighbourExpression)  
              binDistTADsame[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistTADsame[,2][i] <- paste0(distList[[i]][1], ' <= x < ', distList[[i]][2])
              binDistTADsame[,3][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistTADsame[,4][i] <-  as.numeric(binCorr$n[,1][2])
       }
       binDistTADsame$TAD <- 'Same'
       
       binDistTADdiff <- data.frame(matrix(ncol=4,nrow=length(distList))) 
       colnames(binDistTADdiff) <- c('Pearsons', 'Dist', 'Pvalue', 'n')
       for (i in 1:length(distList)) {
              binCorr <- rcorr( CorrTAD[CorrTAD$V3 >= distList[[i]][1] & CorrTAD$V3 < distList[[i]][2] & CorrTAD$TAD!=CorrTAD$neighTAD,]$Expression  , CorrTAD[CorrTAD$V3 >= distList[[i]][1] & CorrTAD$V3 < distList[[i]][2] & CorrTAD$TAD!=CorrTAD$neighTAD,]$NeighbourExpression)  
              binDistTADdiff[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistTADdiff[,2][i] <- paste0(distList[[i]][1], ' <= x < ', distList[[i]][2])
              binDistTADdiff[,3][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistTADdiff[,4][i] <- as.numeric(binCorr$n[,1][2])
       }
       binDistTADdiff$TAD <- 'Different'
       bothTAD <- rbind(binDistTADsame, binDistTADdiff)
       
       
       library(ggrepel)
       
       ggTAD_Distance_bins <- bothTAD %>% 
              mutate(Dist = gsub('.*x ','', Dist), 
                     Dist = factor(Dist, levels=unique(Dist)),
                     SampleName = opt$SampleName) %>%
              ggplot(aes(x = Dist, y = Pearsons, group = TAD, colour =TAD)) +
                     geom_line(size=1) +
                      geom_point(size=4, pch=20) +
                     #geom_bar(stat='identity', color='black', alpha=0.6, position = position_dodge()) +
                     theme_classic() +
                     ylab('Pearsons correlation') +
                     scale_colour_brewer(palette='Set2') +
                     scale_fill_brewer(palette='Set2') +
                     #geom_text_repel(aes(label=Pvalue), point.padding = 0.5)+
                     theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
                     ggtitle(paste0(opt$SampleName, '\nAll vs all pairwise comparison of neighbours (< 20)\nand TADs of neighbours\n',TADlist[[TAD]][1]))+
                     xlab('Distance bins')+
                     ylim(-.20,0.5)
                     
       ggsave( paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_TAD_Distance_bins_',TADlist[[TAD]][1],'.pdf'), plot=ggTAD_Distance_bins ,height=5, width=7) 
       
}

#dir.create('~/public_html/blog/2017Nov20/CMVcorr')
#for (i in unique(neighLoc$neighbour)) {
#       #pdf(paste0(opt$OUTDIR,'/graphs/Neighbour1_distCorr.pdf'), height=8, width=12) 
#       pdf(paste0('~/public_html/blog/2017Nov20/CMVcorr/',i,'.pdf'), height=8, width=12) 
#       par(mfrow=c(2,2))
#       corrDesert <- rcorr( neighLoc[neighLoc$InsLoc=='geneDesert' & neighLoc$NeighLoc=='geneDesert' & neighLoc$neighbour==i,]$Expression, 
#              neighLoc[neighLoc$InsLoc=='geneDesert' & neighLoc$NeighLoc=='geneDesert' & neighLoc$neighbour==i,]$NeighbourExpression)
#       
#       corrRich <- rcorr( neighLoc[neighLoc$InsLoc=='geneRich' & neighLoc$NeighLoc=='geneRich' & neighLoc$neighbour==i,]$Expression, 
#              neighLoc[neighLoc$InsLoc=='geneRich' & neighLoc$NeighLoc=='geneRich' & neighLoc$neighbour==i,]$NeighbourExpression)
#       
#       pRich <- hist(log10(neighLoc[neighLoc$InsLoc=='geneRich' & neighLoc$NeighLoc=='geneRich' & neighLoc$neighbour==i,]$V3+1), n=100, col='#beaed4', xlab='Log10(distance to 1st neighbourig insertion+1)', main='geneRich', frame=F)
#       abline(v=mean(log10(neighLoc[neighLoc$InsLoc=='geneRich' & neighLoc$NeighLoc=='geneRich' & neighLoc$neighbour==i,]$V3+1)), lty=2, col='red', lwd=2)
#       text(paste0("n = ", nrow(neighLoc[neighLoc$InsLoc=='geneRich' & neighLoc$NeighLoc=='geneRich' & neighLoc$neighbour==i,]),
#       '\nmean = ',signif(mean(log10(neighLoc[neighLoc$InsLoc=='geneRich' & neighLoc$NeighLoc=='geneRich' & neighLoc$neighbour==i,]$V3+1)),2)), y=max(pRich$counts)*0.8, 
#       x=mean(log10(neighLoc[neighLoc$InsLoc=='geneRich' & neighLoc$NeighLoc=='geneRich' & neighLoc$neighbour==i,]$V3+1))*0.95)
#       
#       
#       pDesert <- hist(log10(neighLoc[neighLoc$InsLoc=='geneDesert' & neighLoc$NeighLoc=='geneDesert' & neighLoc$neighbour==i,]$V3+1), n=100, col='#fdc086', xlab='Log10(distance to 1st neighbourig insertion+1)', main='geneDesert', frame=F)
#       abline(v=mean(log10(neighLoc[neighLoc$InsLoc=='geneDesert' & neighLoc$NeighLoc=='geneDesert' & neighLoc$neighbour==i,]$V3+1)), lty=2, col='red', lwd=2)
#       text(paste0("n = ", nrow(neighLoc[neighLoc$InsLoc=='geneDesert' & neighLoc$NeighLoc=='geneDesert' & neighLoc$neighbour==i,]),
#       '\nmean = ',signif(mean(log10(neighLoc[neighLoc$InsLoc=='geneDesert' & neighLoc$NeighLoc=='geneDesert' & neighLoc$neighbour==i,]$V3+1)),2)), y=max(pDesert$counts)*0.8, 
#       x=mean(log10(neighLoc[neighLoc$InsLoc=='geneDesert' & neighLoc$NeighLoc=='geneDesert' & neighLoc$neighbour==i,]$V3+1))*0.95)
#       
#       scatter.smooth(neighLoc[neighLoc$InsLoc=='geneRich' & neighLoc$NeighLoc=='geneRich' & neighLoc$neighbour==i,]$Expression, 
#              neighLoc[neighLoc$InsLoc=='geneRich' & neighLoc$NeighLoc=='geneRich' & neighLoc$neighbour==i,]$NeighbourExpression, frame=F, xlab='Activity Insertion', ylab='Activity neighour 1', col='#beaed4', ylim=c(-2,2), xlim=c(-2,2))
#       text(paste0("Pearson's = ", signif(corrRich$r[1,2], 2)), y=1.5, x=-1.5)
#       
#       scatter.smooth(neighLoc[neighLoc$InsLoc=='geneDesert' & neighLoc$NeighLoc=='geneDesert' & neighLoc$neighbour==i,]$Expression, 
#       neighLoc[neighLoc$InsLoc=='geneDesert' & neighLoc$NeighLoc=='geneDesert' & neighLoc$neighbour==i,]$NeighbourExpression, frame=F, xlab='Activity Insertion', ylab='Activity neighour 1', col='#fdc086', ylim=c(-2,2), xlim=c(-2,2))
#       text(paste0("Pearson's = ", signif(corrDesert$r[1,2], 2)), y=1.5, x=-1.5)
#       dev.off()
#}

####
#Neighbour with CTCF in between
####

cat('Number of CTCF between neighbours\n')

cat('Number of CTCF between nearest neighbours\n')

neighCTCF <- read.table(paste0(opt$OUTDIR,'/Neighbour_CTCF.bed'), sep='\t' , header=T, stringsAsFactors=F)
neighCTCF <- neighCTCF[neighCTCF$BC %in% normdataiPCR$BC,]
neighCTCF <- neighCTCF[neighCTCF$neighBC %in% normdataiPCR$BC,]
neighCTCF$Expression <- normdataiPCR[match(neighCTCF$BC, normdataiPCR$BC),]$expression
neighCTCF$NeighbourExpression <- normdataiPCR[match(neighCTCF$neighBC, normdataiPCR$BC),]$expression

CTCFcorr <- data.frame(matrix(ncol=2, nrow=6))
colnames(CTCFcorr) <- c('Pearsons','nCTCF')
for (i in 1:5) {
       CTCFcorr[,1][i] <- as.numeric(rcorr(neighCTCF[neighCTCF$nCTCF==(i-1),]$Expression, neighCTCF[neighCTCF$nCTCF==(i-1),]$NeighbourExpression)$r[,1][2])
       CTCFcorr[,2][i] <- i-1
}
CTCFcorr[,1][nrow(CTCFcorr)] <- as.numeric(rcorr(neighCTCF[neighCTCF$nCTCF>=5,]$Expression, neighCTCF[neighCTCF$nCTCF>=5,]$NeighbourExpression)$r[,1][2])
CTCFcorr[,2][nrow(CTCFcorr)] <- '>=5'

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_CTCF.pdf'), height=5, width=4) 
CTCFcorr %>%
       mutate(SampleName=opt$SampleName,
       nCTCF = factor(nCTCF, levels=c(0, 1, 2, 3, 4, '>=5'))) %>%
       ggplot(aes(x=nCTCF, y=Pearsons, fill=SampleName))+
              geom_bar(stat='identity', col='black', alpha=0.6)+
              theme_classic() +
              scale_fill_manual(values=sampleColor) +
              ylab('Pearson correlation') +
              xlab('Number of CTCF between nearest neighbours') +
              ylim(c(0, 0.5)) +
              ggtitle(opt$SampleName, 'Correlation between nearest neighbours\nand number of CTCF sites in between')+
              theme(legend.position='none')
dev.off()

CTCFdist <- neighCTCF
CTCFdist[CTCFdist$nCTCF>=5,]$nCTCF <-'>=5'

cat('Number of CTCF between nearest neighbours - distance\n')


pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_CTCF_Distance.pdf'), height=5, width=4) 
CTCFdist %>%
   mutate(SampleName=opt$SampleName,
       nCTCF = factor(nCTCF, levels=c(0, 1, 2, 3, 4, '>=5'))) %>%
       ggplot(aes(x=nCTCF, y=log10(Distance+1), fill=SampleName))+
              geom_boxplot(outlier.shape=NA, alpha=0.6)+
              theme_classic() +
              scale_fill_manual(values=sampleColor) +
              ylab('Distance between neighbours') +
              xlab('Number of CTCF between nearest neighbours') +
              #ylim(c(0, 0.5)) +
              ggtitle(opt$SampleName, 'Distance between nearest neighbours and\nnumber of CTCF sites in between')+
              theme(legend.position='none')
dev.off()





#####
#Distance effect on correlation of nearest neighbour
#####
cat('Binned distance and correlation between nearest neighbours\n')

binDistCorr <- data.frame(matrix(ncol=3,nrow=6)) 
colnames(binDistCorr) <- c('Pearsons', 'Dist', 'Pvalue')
distList <- list(c(0, 1e4), c(1e4, 4e4), c(4e4, 1e5), c(1e5, 2e5), c(2e5, 5e5), c(5e5, 1e6)) 
for (i in 1:length(distList)) {
       binCorr <- rcorr( neighCTCF[neighCTCF$Distance >= distList[[i]][1] & neighCTCF$Distance < distList[[i]][2],]$Expression  , neighCTCF[neighCTCF$Distance >= distList[[i]][1] & neighCTCF$Distance < distList[[i]][2],]$NeighbourExpression)  
       binDistCorr[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
       binDistCorr[,2][i] <- paste0(distList[[i]][1], ' <= x < ', distList[[i]][2])
       binDistCorr[,3][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
}

binDistCorr$Dist = gsub('.*< ','< ',binDistCorr$Dist)
pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_BinnedDistance.pdf'), height=5, width=4) 
binDistCorr %>%
       mutate(Dist = factor(Dist, levels=binDistCorr$Dist),
              SampleName = opt$SampleName) %>%
       ggplot(aes(x = Dist, y = Pearsons, fill = SampleName)) +
              geom_bar(stat='identity', color='black', alpha=0.6) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_fill_manual(values=sampleColor) +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'Correlation between nearest\nneighbour and binned distance')+
              theme(legend.position='none')+
              xlab('Distance to nearest neighbour')+
              ylim(0,0.5)
dev.off()


####
#CTCF divide by distance and number of CTCF
####
cat('Binned distance and correlation between nearest neighbours and CTCF sites\n')


distList <- list(c(0, 5e4, 0), c(5e4,  1e5, 0), c(1e5, 2e5, 0), c(2e5, 5e5, 0),
                     c(0, 5e4, 1), c(5e4, 1e5, 1), c(1e5, 2e5, 1), c(2e5, 5e5, 1),
                     c(0, 5e4, 2), c(5e4, 1e5, 2), c(1e5, 2e5, 2), c(2e5, 5e5, 2),
                     c(0, 5e4, 3), c(5e4, 1e5, 3), c(1e5, 2e5, 3), c(2e5, 5e5, 3))
                     
                     
binDistCorr <- data.frame(matrix(ncol=4,nrow=length(distList))) 
colnames(binDistCorr) <- c('Pearsons', 'Dist', 'CTCF','Pvalue')
for (i in 1:length(distList)) { 
       if( distList[[i]][3] <3) {
              binCorr <- rcorr( neighCTCF[neighCTCF$Distance >= distList[[i]][1] & neighCTCF$Distance < distList[[i]][2] & neighCTCF$nCTCF == distList[[i]][3] ,]$Expression  , neighCTCF[neighCTCF$Distance >= distList[[i]][1] & neighCTCF$Distance < distList[[i]][2] & neighCTCF$nCTCF == distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorr[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorr[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorr[,3][i] <- distList[[i]][3]
              binDistCorr[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
       } else if(distList[[i]][3] ==3) {
              binCorr <- rcorr( neighCTCF[neighCTCF$Distance >= distList[[i]][1] & neighCTCF$Distance < distList[[i]][2] & neighCTCF$nCTCF >= distList[[i]][3] ,]$Expression  , neighCTCF[neighCTCF$Distance >= distList[[i]][1] & neighCTCF$Distance < distList[[i]][2] & neighCTCF$nCTCF >= distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorr[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorr[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorr[,3][i] <- paste0('>= ',distList[[i]][3])
              binDistCorr[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
       }
}


pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_BinnedDistance_CTCF.pdf'), height=5, width=4) 
binDistCorr %>%
       mutate(Dist = factor(Dist, levels=unique(binDistCorr$Dist)),
              SampleName = opt$SampleName,
              CTCF = factor(CTCF, levels=c(0, 1, 2, '>= 3'))) %>%
       ggplot(aes(x = Dist, y = Pearsons, fill = CTCF)) +
              geom_bar(stat='identity', color='black', alpha=0.6, position = position_dodge()) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_fill_brewer(palette='Set2') +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'Correlation between nearest\nneighbour and binned distance\nand number CTCF sites in between')+
              xlab('Distance to nearest neighbour')+
              ylim(0,0.5)
dev.off()



#####
#All pairwise correlation between distance and CTCF
#####
cat('All pairwise with CTCF sites\n')
library(ggrepel)
library(dplyr)
library(Hmisc)
allCTCF <- read.table(paste0(opt$OUTDIR, '/AllNeighbour_CTCF.bed'), sep='\t', header=F)
colnames(allCTCF) <- c('chrom', 'start', 'end', 'BC', 'neighBC', 'Distance', 'nCTCF')

allCTCF$Expression <- normdataiPCR[match(allCTCF$BC, normdataiPCR$BC),]$expression
allCTCF$NeighbourExpression <- normdataiPCR[match(allCTCF$neighBC, normdataiPCR$BC),]$expression

distList <- list(c(0, 5e4, 0), c(5e4, 1e5, 0), c(1e5, 2e5, 0), c(2e5, 3e5, 0), c(3e5, 4e5, 0), c(4e5, 5e5, 0), c(5e5, 1e6, 0),
              c(0, 5e4, 1), c(5e4, 1e5, 1), c(1e5, 2e5, 1), c(2e5, 3e5, 1), c(3e5, 4e5, 1), c(4e5, 5e5, 1), c(5e5, 1e6, 1))
binDistCorrAllvsAllCTCF <- data.frame(matrix(ncol=5,nrow=length(distList))) 
colnames(binDistCorrAllvsAllCTCF) <- c('Pearsons', 'Dist', 'CTCF', 'Pvalue', 'n')

for (i in 1:length(distList)) { 
       if( distList[[i]][3] <1 & nrow( allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF == distList[[i]][3] ,]) >=50) {
              binCorr <- rcorr( allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF == distList[[i]][3] ,]$Expression  , allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF == distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorrAllvsAllCTCF[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCTCF[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorrAllvsAllCTCF[,3][i] <- distList[[i]][3]
              binDistCorrAllvsAllCTCF[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCTCF[,5][i] <- as.numeric(binCorr$n[,1][2])
       } else if(distList[[i]][3] ==1 & nrow(allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF >= distList[[i]][3] ,])) {
              binCorr <- rcorr( allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF >= distList[[i]][3] ,]$Expression  , allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF >= distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorrAllvsAllCTCF[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCTCF[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorrAllvsAllCTCF[,3][i] <- paste0('>= ',distList[[i]][3])
              binDistCorrAllvsAllCTCF[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCTCF[,5][i] <- as.numeric(binCorr$n[,1][2])
       }
}

binDistCorrAllvsAllCTCF <- binDistCorrAllvsAllCTCF[!is.na(binDistCorrAllvsAllCTCF$Dist),]
binDistCorrAllvsAllCTCF <- binDistCorrAllvsAllCTCF[binDistCorrAllvsAllCTCF$Pvalue<=0.05,]

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_BinnedDistance_CTCF_allPairwise.pdf'), height=5, width=6) 
binDistCorrAllvsAllCTCF %>%
       mutate(Dist = factor(Dist, levels=unique(binDistCorrAllvsAllCTCF$Dist)),
              SampleName = opt$SampleName,
              CTCF = factor(CTCF, levels=c(0,  '>= 1'))) %>%
       ggplot(aes(x = Dist, y = Pearsons, group = CTCF, colour =CTCF)) +
              geom_line(size=1) +
               geom_point(size=4, pch=20) +
              #geom_bar(stat='identity', color='black', alpha=0.6, position = position_dodge()) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_colour_brewer(palette='Set2') +
              scale_fill_brewer(palette='Set2') +
              geom_text_repel(aes(label=Pvalue), point.padding = 0.5)+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All vs all pairwise comparison of neighbours (< 20)\nand number of CTCF sites in between')+
              xlab('Distance bins')+
              ylim(0,0.5)
dev.off()

######
#CTCF divided into multiple parts
##########
#All pairwise correlation between distance and CTCF
#####
cat('All pairwise with CTCF sites\n')
library(ggrepel)
library(dplyr)
library(Hmisc)
allCTCF <- read.table(paste0(opt$OUTDIR, '/AllNeighbour_CTCF.bed'), sep='\t', header=F)
colnames(allCTCF) <- c('chrom', 'start', 'end', 'BC', 'neighBC', 'Distance', 'nCTCF')

allCTCF$Expression <- normdataiPCR[match(allCTCF$BC, normdataiPCR$BC),]$expression
allCTCF$NeighbourExpression <- normdataiPCR[match(allCTCF$neighBC, normdataiPCR$BC),]$expression

distList <- list(c(0, 5e4, 0), c(5e4, 1e5, 0), c(1e5, 2e5, 0), c(2e5, 3e5, 0), c(3e5, 4e5, 0), c(4e5, 5e5, 0), c(5e5, 1e6, 0),
              c(0, 5e4, 1), c(5e4, 1e5, 1), c(1e5, 2e5, 1), c(2e5, 3e5, 1), c(3e5, 4e5, 1), c(4e5, 5e5, 1), c(5e5, 1e6, 1),
              c(0, 5e4, 2), c(5e4, 1e5, 2), c(1e5, 2e5, 2), c(2e5, 3e5, 2), c(3e5, 4e5, 2), c(4e5, 5e5, 2), c(5e5, 1e6, 2),
              c(0, 5e4, 3), c(5e4, 1e5, 3), c(1e5, 2e5, 3), c(2e5, 3e5, 3), c(3e5, 4e5, 3), c(4e5, 5e5, 3), c(5e5, 1e6, 3))
binDistCorrAllvsAllCTCFMultiple <- data.frame(matrix(ncol=5,nrow=length(distList))) 
colnames(binDistCorrAllvsAllCTCFMultiple) <- c('Pearsons', 'Dist', 'CTCF', 'Pvalue', 'n')

for (i in 1:length(distList)) { 
       if( distList[[i]][3] <3 & nrow( allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF == distList[[i]][3] ,]) >=50) {
              binCorr <- rcorr( allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF == distList[[i]][3] ,]$Expression  , allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF == distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorrAllvsAllCTCFMultiple[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorrAllvsAllCTCFMultiple[,3][i] <- distList[[i]][3]
              binDistCorrAllvsAllCTCFMultiple[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,5][i] <- as.numeric(binCorr$n[,1][2])
       } else if(distList[[i]][3] ==3 & nrow(allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF >= distList[[i]][3] ,])) {
              binCorr <- rcorr( allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF >= distList[[i]][3] ,]$Expression  , allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF >= distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorrAllvsAllCTCFMultiple[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorrAllvsAllCTCFMultiple[,3][i] <- paste0('>= ',distList[[i]][3])
              binDistCorrAllvsAllCTCFMultiple[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,5][i] <- as.numeric(binCorr$n[,1][2])
       }
}

binDistCorrAllvsAllCTCFMultiple <- binDistCorrAllvsAllCTCFMultiple[!is.na(binDistCorrAllvsAllCTCFMultiple$Dist),]
binDistCorrAllvsAllCTCFMultiple <- binDistCorrAllvsAllCTCFMultiple[binDistCorrAllvsAllCTCFMultiple$Pvalue<=0.05,]

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_BinnedDistance_CTCF_allPairwise_upTOfour.pdf'), height=5, width=6) 
binDistCorrAllvsAllCTCFMultiple %>%
       mutate(Dist = factor(Dist, levels=unique(binDistCorrAllvsAllCTCFMultiple$Dist)),
              SampleName = opt$SampleName) %>%
             # CTCF = factor(CTCF, levels=c(0,  '>= 1'))) %>%
       ggplot(aes(x = Dist, y = Pearsons, group = CTCF, colour =CTCF)) +
              geom_line(size=1) +
               geom_point(size=4, pch=20) +
              #geom_bar(stat='identity', color='black', alpha=0.6, position = position_dodge()) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_colour_brewer(palette='Set2') +
              scale_fill_brewer(palette='Set2') +
              #geom_text_repel(aes(label=Pvalue), point.padding = 0.5)+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All vs all pairwise comparison of neighbours (< 20)\nand number of CTCF sites in between')+
              xlab('Distance bins')+
              ylim(0,0.5)
dev.off()



##########
#All pairwise correlation between decile distance and CTCF
#####
cat('nCTCF and distance\n')
pdf('~/public_html/blog/2017Dec04/nCTCF_distance.pdf')
plot(log10(allCTCF$Distance+1), allCTCF$nCTCF)
dev.off()

distList <- list(c(0, 5e4, 0), c(5e4, 1e5, 0), c(1e5, 2e5, 0), c(2e5, 3e5, 0), c(3e5, 4e5, 0), c(4e5, 5e5, 0), c(5e5, 6e5, 0), c(6e5, 7e5, 0),  c(7e5, 8e5, 0), c(8e5, 9e5, 0),  c(9e5, 1e6, 0),
              c(0, 5e4, 1), c(5e4, 1e5, 1), c(1e5, 2e5, 1), c(2e5, 3e5, 1), c(3e5, 4e5, 1), c(4e5, 5e5, 1), c(5e5, 6e5, 1), c(6e5, 7e5, 1),  c(7e5, 8e5, 1), c(8e5, 9e5, 1),  c(9e5, 1e6, 1),
              c(0, 5e4, 2), c(5e4, 1e5, 2), c(1e5, 2e5, 2), c(2e5, 3e5, 2), c(3e5, 4e5, 2), c(4e5, 5e5, 2), c(5e5, 6e5, 2), c(6e5, 7e5, 2),  c(7e5, 8e5, 2), c(8e5, 9e5, 2),  c(9e5, 1e6, 2),
              c(0, 5e4, 3), c(5e4, 1e5, 3), c(1e5, 2e5, 3), c(2e5, 3e5, 3), c(3e5, 4e5, 3), c(4e5, 5e5, 3), c(5e5, 6e5, 3), c(6e5, 7e5, 3),  c(7e5, 8e5, 3), c(8e5, 9e5, 3),  c(9e5, 1e6, 3))
binDistCorrAllvsAllCTCFMultiple <- data.frame(matrix(ncol=5,nrow=length(distList))) 
colnames(binDistCorrAllvsAllCTCFMultiple) <- c('Pearsons', 'Dist', 'CTCF', 'Pvalue', 'n')

for (i in 1:length(distList)) { 
       if( distList[[i]][3] <3 & nrow( allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF == distList[[i]][3] ,]) >=50) {
              binCorr <- rcorr( allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF == distList[[i]][3] ,]$Expression  , allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF == distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorrAllvsAllCTCFMultiple[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorrAllvsAllCTCFMultiple[,3][i] <- distList[[i]][3]
              binDistCorrAllvsAllCTCFMultiple[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,5][i] <- as.numeric(binCorr$n[,1][2])
       } else if(distList[[i]][3] ==3 & nrow(allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF >= distList[[i]][3] ,])) {
              binCorr <- rcorr( allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF >= distList[[i]][3] ,]$Expression  , allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF >= distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorrAllvsAllCTCFMultiple[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorrAllvsAllCTCFMultiple[,3][i] <- paste0('>= ',distList[[i]][3])
              binDistCorrAllvsAllCTCFMultiple[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,5][i] <- as.numeric(binCorr$n[,1][2])
       }
}

binDistCorrAllvsAllCTCFMultiple <- binDistCorrAllvsAllCTCFMultiple[!is.na(binDistCorrAllvsAllCTCFMultiple$Dist),]
binDistCorrAllvsAllCTCFMultiple <- binDistCorrAllvsAllCTCFMultiple[binDistCorrAllvsAllCTCFMultiple$Pvalue<=0.05,]
cat('Neighbouring Insertions BinnedDistance CTCF allPairwise upTOfour\n')

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_BinnedDistance_CTCF_allPairwise_upTOfour.pdf'), height=5, width=6) 
binDistCorrAllvsAllCTCFMultiple %>%
       mutate(Dist = factor(Dist, levels=unique(binDistCorrAllvsAllCTCFMultiple$Dist)),
              SampleName = opt$SampleName) %>%
             # CTCF = factor(CTCF, levels=c(0,  '>= 1'))) %>%
       ggplot(aes(x = Dist, y = Pearsons, group = CTCF, colour =CTCF)) +
              geom_line(size=1) +
               geom_point(size=4, pch=20) +
              #geom_bar(stat='identity', color='black', alpha=0.6, position = position_dodge()) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_colour_brewer(palette='Set2') +
              scale_fill_brewer(palette='Set2') +
              #geom_text_repel(aes(label=Pvalue), point.padding = 0.5)+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All vs all pairwise comparison of neighbours (< 20)\nand number of CTCF sites in between')+
              xlab('Distance bins')+
              ylim(0,0.5)
dev.off()





######
#CTCF sites surrounding insertions
######
surrCTCF <- read.table(paste0(opt$OUTDIR,'/AllNeighbour_surroundingCTCF.bed'), sep='\t', header=T)
surrCTCF$InsExp <- normdataiPCR[match(surrCTCF$InsBC, normdataiPCR$BC),]$expression
surrCTCF$neighExp <- normdataiPCR[match(surrCTCF$neighBC, normdataiPCR$BC),]$expression

surrCTCF$InsUpCTCF <- abs(surrCTCF$InsUpCTCF)
surrCTCF$neighUpCTCF <- abs(surrCTCF$neighUpCTCF)

surrCTCFZero <- surrCTCF[surrCTCF$nCTCF==0,]
surrCTCFOne <- surrCTCF[surrCTCF$nCTCF>0,]

surrNear <- rcorr(surrCTCFZero[surrCTCFZero$InsUpCTCF<=1e4 & surrCTCFZero$neighDownCTCF<=1e4,]$InsExp, surrCTCFZero[surrCTCFZero$InsUpCTCF<=1e4 & surrCTCFZero$neighDownCTCF<=1e4,]$neighExp,)
surrFar <- rcorr(surrCTCFZero[surrCTCFZero$InsUpCTCF>1e4 & surrCTCFZero$neighDownCTCF>1e4,]$InsExp, surrCTCFZero[surrCTCFZero$InsUpCTCF>1e4 & surrCTCFZero$neighDownCTCF>1e4,]$neighExp,)
surrNearFar <- rcorr(surrCTCFZero[(surrCTCFZero$InsUpCTCF>1e4 & surrCTCFZero$neighDownCTCF<=1e4) | (surrCTCFZero$InsUpCTCF<=1e4 & surrCTCFZero$neighDownCTCF>1e4),]$InsExp, 
       surrCTCFZero[(surrCTCFZero$InsUpCTCF>1e4 & surrCTCFZero$neighDownCTCF<=1e4) | (surrCTCFZero$InsUpCTCF<=1e4 & surrCTCFZero$neighDownCTCF>1e4),]$neighExp,)
       
surrOneAll <- rcorr(surrCTCFOne$InsExp, surrCTCFOne$neighExp)

surrCTCFall <- data.frame(CTCFDist=c('both <=10kb', 'both >10kb', 'One <=10kb, other >10kb', 'Inbetween'), 
       Pearsons = c(signif(as.numeric(surrNear$r[,1][2]),3), signif(as.numeric(surrFar$r[,1][2]),3), signif(as.numeric(surrNearFar$r[,1][2]),3), signif(as.numeric(surrOneAll$r[,1][2]),3)),
       Pvalue = c(signif(as.numeric(surrNear$P[,1][2]),3), signif(as.numeric(surrFar$P[,1][2]),3), signif(as.numeric(surrNearFar$P[,1][2]),3), signif(as.numeric(surrOneAll$P[,1][2]),3)),
       n = c(nrow(surrCTCFZero[surrCTCFZero$InsUpCTCF<=1e4 & surrCTCFZero$neighDownCTCF<=1e4,]), 
              nrow(surrCTCFZero[surrCTCFZero$InsUpCTCF>1e4 & surrCTCFZero$neighDownCTCF>1e4,]), 
              nrow(surrCTCFZero[(surrCTCFZero$InsUpCTCF>1e4 & surrCTCFZero$neighDownCTCF<=1e4) | (surrCTCFZero$InsUpCTCF<=1e4 & surrCTCFZero$neighDownCTCF>1e4),]),
               nrow(surrCTCFOne)),
       nCTCF = c(0, 0, 0, '>=1'))

cat('Surrounding CTCF and Correlation\n')

pdf(paste0(opt$OUTDIR,'/graphs/SurroundingCTCFandCorrelation.pdf'), height=5, width=4) 

surrCTCFall %>%
       mutate(SampleName = opt$SampleName,
              CTCFDist = factor(CTCFDist, levels=c('both <=10kb',  'both >10kb', 'One <=10kb, other >10kb','Inbetween'))) %>%
       ggplot(aes(x = CTCFDist, y = Pearsons, fill = SampleName)) +
              geom_bar(stat='identity', color='black', alpha=0.6) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_fill_manual(values=sampleColor) +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All pariwise correlation between neighbours\nwithout CTCF site in between and\ndistance to CTCF')+
              xlab('Distance to CTCF sites')+
              ylim(0,0.5)+
              theme(legend.position='none')
dev.off()


####
#Redo with zero CTCF, 1 CTCF +- 10kb, 1 CTCF in between
####
cat('Redo with zero CTCF, 1 CTCF +- 10kb, 1 CTCF in between\n')

noCTCF <- rcorr(surrCTCF[surrCTCF$InsUpCTCF>1e4 & surrCTCF$neighDownCTCF>1e4 & surrCTCF$nCTCF==0,]$InsExp, surrCTCF[surrCTCF$InsUpCTCF>1e4 & surrCTCF$neighDownCTCF>1e4 & surrCTCF$nCTCF==0,]$neighExp,)
oneCTCF <- rcorr(surrCTCF[surrCTCF$InsUpCTCF<=1e4| 
                     surrCTCF$InsDownCTCF<=1e4|
                     surrCTCF$neighUpCTCF<=1e4 |
                     surrCTCF$neighDownCTC<=1e4 ,]$InsExp, 
              surrCTCF[surrCTCF$InsUpCTCF<=1e4| 
                     surrCTCF$InsDownCTCF<=1e4|
                     surrCTCF$neighUpCTCF<=1e4 |
                     surrCTCF$neighDownCTC<=1e4 ,]$neighExp,)
oneBetweenCTCF <- rcorr(surrCTCF[surrCTCF$nCTCF==1 ,]$InsExp, surrCTCF[(surrCTCF$nCTCF>=1) ,]$neighExp,)
#twoBetweenCTCF <- rcorr(surrCTCF[surrCTCF$nCTCF==2 ,]$InsExp, surrCTCF[(surrCTCF$nCTCF==2) ,]$neighExp,)
#threeBetweenCTCF <- rcorr(surrCTCF[surrCTCF$nCTCF==3 ,]$InsExp, surrCTCF[(surrCTCF$nCTCF==3) ,]$neighExp,)
#OverthreeBetweenCTCF <- rcorr(surrCTCF[surrCTCF$nCTCF>3 ,]$InsExp, surrCTCF[(surrCTCF$nCTCF>3) ,]$neighExp,)

surrCTCFall <- data.frame(CTCFDist=c('no CTCF', 'CTCF within +-10kb', '>=1 CTCF between insertions'), 
       Pearsons = c(signif(as.numeric(noCTCF$r[,1][2]),3), signif(as.numeric(oneCTCF$r[,1][2]),3), signif(as.numeric(oneBetweenCTCF$r[,1][2]),3)),
       Pvalue = c(signif(as.numeric(noCTCF$P[,1][2]),3), signif(as.numeric(oneCTCF$P[,1][2]),3), signif(as.numeric(oneBetweenCTCF$P[,1][2]),3)))

pdf(paste0(opt$OUTDIR,'/graphs/SurroundingCTCFandCorrelation_NoCTCF.pdf'), height=5, width=4) 
surrCTCFall %>%
       mutate(SampleName = opt$SampleName,
              CTCFDist = factor(CTCFDist, levels=c('no CTCF', 'CTCF within +-10kb', '>=1 CTCF between insertions'))) %>%
       ggplot(aes(x = CTCFDist, y = Pearsons, fill = SampleName)) +
              geom_bar(stat='identity', color='black', alpha=0.6) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_fill_manual(values=sampleColor) +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All pariwise correlation between neighbours\nwith or without CTCF sites present')+
              xlab('Distance to CTCF sites')+
              ylim(0,0.5)+
              theme(legend.position='none')
dev.off()

####
#Redo with zero CTCF, 1 CTCF +- 10kb, 1 CTCF in between
####
cat('Redo with zero CTCF, 1 CTCF +- 10kb, 1,2, 3 CTCF in between\n')

noCTCF <- rcorr(surrCTCF[surrCTCF$InsUpCTCF>1e4 & surrCTCF$neighDownCTCF>1e4 & surrCTCF$nCTCF==0,]$InsExp, surrCTCF[surrCTCF$InsUpCTCF>1e4 & surrCTCF$neighDownCTCF>1e4 & surrCTCF$nCTCF==0,]$neighExp,)
oneCTCF <- rcorr(surrCTCF[surrCTCF$InsUpCTCF<=1e4| 
                     surrCTCF$InsDownCTCF<=1e4|
                     surrCTCF$neighUpCTCF<=1e4 |
                     surrCTCF$neighDownCTC<=1e4 ,]$InsExp, 
              surrCTCF[surrCTCF$InsUpCTCF<=1e4| 
                     surrCTCF$InsDownCTCF<=1e4|
                     surrCTCF$neighUpCTCF<=1e4 |
                     surrCTCF$neighDownCTC<=1e4 ,]$neighExp,)
oneBetweenCTCF <- rcorr(surrCTCF[surrCTCF$nCTCF==1 ,]$InsExp, surrCTCF[(surrCTCF$nCTCF==1) ,]$neighExp,)
twoBetweenCTCF <- rcorr(surrCTCF[surrCTCF$nCTCF==2 ,]$InsExp, surrCTCF[(surrCTCF$nCTCF==2) ,]$neighExp,)
threeBetweenCTCF <- rcorr(surrCTCF[surrCTCF$nCTCF==3 ,]$InsExp, surrCTCF[(surrCTCF$nCTCF==3) ,]$neighExp,)
OverthreeBetweenCTCF <- rcorr(surrCTCF[surrCTCF$nCTCF>3 ,]$InsExp, surrCTCF[(surrCTCF$nCTCF>3) ,]$neighExp,)

surrCTCFall <- data.frame(CTCFDist=c('no CTCF', 'CTCF within +-10kb', '1 CTCF between insertions', '2 CTCF between insertions', '3 CTCF between insertions', '>=4 CTCF between insertions'), 
       Pearsons = c(signif(as.numeric(noCTCF$r[,1][2]),3), 
                     signif(as.numeric(oneCTCF$r[,1][2]),3), 
                     signif(as.numeric(oneBetweenCTCF$r[,1][2]),3), 
                     signif(as.numeric(twoBetweenCTCF$r[,1][2]),3), 
                     signif(as.numeric(threeBetweenCTCF$r[,1][2]),3), 
                     signif(as.numeric(OverthreeBetweenCTCF$r[,1][2]),3)),
       Pvalue = c(signif(as.numeric(noCTCF$P[,1][2]),3), 
                     signif(as.numeric(oneCTCF$P[,1][2]),3), 
                     signif(as.numeric(oneBetweenCTCF$P[,1][2]),3), 
                     signif(as.numeric(twoBetweenCTCF$P[,1][2]),3), 
                     signif(as.numeric(threeBetweenCTCF$P[,1][2]),3), 
                     signif(as.numeric(OverthreeBetweenCTCF$P[,1][2]),3)))

pdf(paste0(opt$OUTDIR,'/graphs/SurroundingCTCFandCorrelation_NoCTCF_multiple.pdf'), height=5, width=4) 
surrCTCFall %>%
       mutate(SampleName = opt$SampleName,
              CTCFDist = factor(CTCFDist, levels=c('no CTCF', 'CTCF within +-10kb', '1 CTCF between insertions', '2 CTCF between insertions', '3 CTCF between insertions', '>=4 CTCF between insertions'))) %>%
       ggplot(aes(x = CTCFDist, y = Pearsons, fill = SampleName)) +
              geom_bar(stat='identity', color='black', alpha=0.6) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_fill_manual(values=sampleColor) +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All pariwise correlation between neighbours\nwith or without CTCF sites present')+
              xlab('Distance to CTCF sites')+
              ylim(0,0.5)+
              theme(legend.position='none')
dev.off()


####
#No CTCF between for one +-10kb and only 1 CTCF between
####
cat('No CTCF between for one +-10kb and only 1 CTCF between\n')

noCTCF <- rcorr(surrCTCF[surrCTCF$InsUpCTCF>1e4 & surrCTCF$neighDownCTCF>1e4 & surrCTCF$nCTCF==0,]$InsExp, surrCTCF[surrCTCF$InsUpCTCF>1e4 & surrCTCF$neighDownCTCF>1e4 & surrCTCF$nCTCF==0,]$neighExp,)
oneCTCF <- rcorr(surrCTCF[(surrCTCF$InsUpCTCF<=1e4 & surrCTCF$nCTCF==0)| 
                     (surrCTCF$InsDownCTCF<=1e4 & surrCTCF$nCTCF==0) |
                     (surrCTCF$neighUpCTCF<=1e4 & surrCTCF$nCTCF==0) |
                     (surrCTCF$neighDownCTC<=1e4 & surrCTCF$nCTCF==0) ,]$InsExp, 
              surrCTCF[(surrCTCF$InsUpCTCF<=1e4 & surrCTCF$nCTCF==0)| 
                     (surrCTCF$InsDownCTCF<=1e4 & surrCTCF$nCTCF==0)|
                     (surrCTCF$neighUpCTCF<=1e4 & surrCTCF$nCTCF==0)|
                     (surrCTCF$neighDownCTC<=1e4 & surrCTCF$nCTCF==0),]$neighExp,)
oneBetweenCTCF <- rcorr(surrCTCF[surrCTCF$nCTCF==1 ,]$InsExp, surrCTCF[(surrCTCF$nCTCF==1) ,]$neighExp,)

surrCTCFall <- data.frame(CTCFDist=c('no CTCF', 'CTCF within +-10kb/No between', '1 CTCF between insertions'), 
       Pearsons = c(signif(as.numeric(noCTCF$r[,1][2]),3), signif(as.numeric(oneCTCF$r[,1][2]),3), signif(as.numeric(oneBetweenCTCF$r[,1][2]),3)),
       Pvalue = c(signif(as.numeric(noCTCF$P[,1][2]),3), signif(as.numeric(oneCTCF$P[,1][2]),3), signif(as.numeric(oneBetweenCTCF$P[,1][2]),3)))

pdf(paste0(opt$OUTDIR,'/graphs/SurroundingCTCFandCorrelation_NoCTCF_NoCTFbetween10kb_1CTCFbetween.pdf'), height=5, width=4) 
surrCTCFall %>%
       mutate(SampleName = opt$SampleName,
              CTCFDist = factor(CTCFDist, levels=c('no CTCF', 'CTCF within +-10kb/No between', '1 CTCF between insertions'))) %>%
       ggplot(aes(x = CTCFDist, y = Pearsons, fill = SampleName)) +
              geom_bar(stat='identity', color='black', alpha=0.6) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_fill_manual(values=sampleColor) +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All pariwise correlation between neighbours\nwith or without CTCF sites present')+
              xlab('Distance to CTCF sites')+
              ylim(0,0.5)+
              theme(legend.position='none')
dev.off()





####
#Divide into distance bins, numer of CTCFs and distance to nearest CTCF (3 lines)
####
cat('Divide into distance bins, numer of CTCFs and distance to nearest CTCF (3 lines)\n')

distList <- list(c(0, 5e4, 1), c(5e4, 1e5,  1), c(1e5, 2e5, 1), c(2e5, 3e5, 1), c(3e5, 4e5, 1), c(4e5, 5e5, 1), c(5e5, 6e6, 1),
              c(0, 5e4, 2), c(5e4, 1e5,  2), c(1e5, 2e5, 2), c(2e5, 3e5, 2), c(3e5, 4e5, 2), c(4e5, 5e5, 2), c(5e5, 6e6, 2),
              c(0, 5e4, 3), c(5e4, 1e5,  3), c(1e5, 2e5, 3), c(2e5, 3e5, 3), c(3e5, 4e5, 3), c(4e5, 5e5, 3), c(5e5, 6e6, 3))
binDistCorrAllvsAllCTCFsurrounding <- data.frame(matrix(ncol=5,nrow=length(distList))) 
colnames(binDistCorrAllvsAllCTCFsurrounding) <- c('Pearsons', 'Dist', 'CTCFdist', 'Pvalue', 'n')

for (i in 1:length(distList)) { 
       if( distList[[i]][3] ==1) {
              if(nrow(surrCTCFZero[surrCTCFZero$Distance >= distList[[i]][1] & surrCTCFZero$Distance < distList[[i]][2] & surrCTCFZero$InsUpCTCF <=1e4 & surrCTCFZero$neighDownCTCF <=1e4 ,])> 50) {
                     binCorr <- rcorr( surrCTCFZero[surrCTCFZero$Distance >= distList[[i]][1] & surrCTCFZero$Distance < distList[[i]][2] & surrCTCFZero$InsUpCTCF <=1e4 & surrCTCFZero$neighDownCTCF <=1e4 ,]$InsExp  , surrCTCFZero[surrCTCFZero$Distance >= distList[[i]][1] & surrCTCFZero$Distance < distList[[i]][2] & surrCTCFZero$InsUpCTCF <=1e4 & surrCTCFZero$neighDownCTCF <=1e4 ,]$neighExp)  
                     binDistCorrAllvsAllCTCFsurrounding[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
                     binDistCorrAllvsAllCTCFsurrounding[,2][i] <- paste0('< ', distList[[i]][2])
                     binDistCorrAllvsAllCTCFsurrounding[,3][i] <- '<=10kb'
                     binDistCorrAllvsAllCTCFsurrounding[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
                     binDistCorrAllvsAllCTCFsurrounding[,5][i] <- as.numeric(binCorr$n[,1][2])
              }
       } else        if( distList[[i]][3] ==2) {
              if(nrow(surrCTCFZero[surrCTCFZero$Distance >= distList[[i]][1] & surrCTCFZero$Distance < distList[[i]][2] & surrCTCFZero$InsUpCTCF >1e4 & surrCTCFZero$neighDownCTCF >1e4 ,])> 50) {
                     binCorr <- rcorr( surrCTCFZero[surrCTCFZero$Distance >= distList[[i]][1] & surrCTCFZero$Distance < distList[[i]][2] & surrCTCFZero$InsUpCTCF >1e4 & surrCTCFZero$neighDownCTCF >1e4 ,]$InsExp  , surrCTCFZero[surrCTCFZero$Distance >= distList[[i]][1] & surrCTCFZero$Distance < distList[[i]][2] & surrCTCFZero$InsUpCTCF >1e4 & surrCTCFZero$neighDownCTCF >1e4 ,]$neighExp)  
                     binDistCorrAllvsAllCTCFsurrounding[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
                     binDistCorrAllvsAllCTCFsurrounding[,2][i] <- paste0('< ', distList[[i]][2])
                     binDistCorrAllvsAllCTCFsurrounding[,3][i] <- '>10kb'
                     binDistCorrAllvsAllCTCFsurrounding[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
                     binDistCorrAllvsAllCTCFsurrounding[,5][i] <- as.numeric(binCorr$n[,1][2])
              }
       } else        if( distList[[i]][3] ==3) {
              if(nrow(surrCTCFZero[(surrCTCFZero$Distance >= distList[[i]][1] & surrCTCFZero$Distance < distList[[i]][2] & surrCTCFZero$InsUpCTCF <=1e4 & surrCTCFZero$neighDownCTCF >1e4) | (surrCTCFZero$Distance >= distList[[i]][1] & surrCTCFZero$Distance < distList[[i]][2] & surrCTCFZero$InsUpCTCF >1e4 & surrCTCFZero$neighDownCTCF <=1e4),])> 50) {
                     binCorr <- rcorr( surrCTCFZero[(surrCTCFZero$Distance >= distList[[i]][1] & surrCTCFZero$Distance < distList[[i]][2] & surrCTCFZero$InsUpCTCF <=1e4 & surrCTCFZero$neighDownCTCF >1e4) | (surrCTCFZero$Distance >= distList[[i]][1] & surrCTCFZero$Distance < distList[[i]][2] & surrCTCFZero$InsUpCTCF >1e4 & surrCTCFZero$neighDownCTCF <=1e4),]$InsExp  , surrCTCFZero[(surrCTCFZero$Distance >= distList[[i]][1] & surrCTCFZero$Distance < distList[[i]][2] & surrCTCFZero$InsUpCTCF <=1e4 & surrCTCFZero$neighDownCTCF >1e4) | (surrCTCFZero$Distance >= distList[[i]][1] & surrCTCFZero$Distance < distList[[i]][2] & surrCTCFZero$InsUpCTCF >1e4 & surrCTCFZero$neighDownCTCF <=1e4),]$neighExp)  
                     binDistCorrAllvsAllCTCFsurrounding[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
                     binDistCorrAllvsAllCTCFsurrounding[,2][i] <- paste0('< ', distList[[i]][2])
                     binDistCorrAllvsAllCTCFsurrounding[,3][i] <- '>10kb and <=10kb'
                     binDistCorrAllvsAllCTCFsurrounding[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
                     binDistCorrAllvsAllCTCFsurrounding[,5][i] <- as.numeric(binCorr$n[,1][2])
              }
       }
       
}
binDistCorrAllvsAllCTCFsurrounding$nCTCF <-0

binDistCorrAllvsAllCTCFsurroundingOne <- data.frame(matrix(ncol=5,nrow=length(distList))) 
colnames(binDistCorrAllvsAllCTCFsurroundingOne) <- c('Pearsons', 'Dist', 'CTCFdist', 'Pvalue', 'n')

for (i in 1:length(distList)) { 
       if( distList[[i]][3] ==1) {
              if(nrow(surrCTCFOne[surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF <=1e4 & surrCTCFOne$neighDownCTCF <=1e4 ,])> 50) {
                     binCorr <- rcorr( surrCTCFOne[(surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF <=1e4 & surrCTCFOne$neighDownCTCF <=1e4) | 
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF <=1e4 & surrCTCFOne$neighDownCTCF <=1e4) |
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsDownCTCF <=1e4 & surrCTCFOne$neighUpCTCF <=1e4) |
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF <=1e4 & surrCTCFOne$neighDownCTCF <=1e4),]$InsExp  , 
                                          surrCTCFOne[(surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF <=1e4 & surrCTCFOne$neighDownCTCF <=1e4) | 
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF <=1e4 & surrCTCFOne$neighDownCTCF <=1e4) |
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsDownCTCF <=1e4 & surrCTCFOne$neighUpCTCF <=1e4) |
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF <=1e4 & surrCTCFOne$neighDownCTCF <=1e4),]$neighExp)  
                     binDistCorrAllvsAllCTCFsurroundingOne[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
                     binDistCorrAllvsAllCTCFsurroundingOne[,2][i] <- paste0('< ', distList[[i]][2])
                     binDistCorrAllvsAllCTCFsurroundingOne[,3][i] <- '<=10kb'
                     binDistCorrAllvsAllCTCFsurroundingOne[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
                     binDistCorrAllvsAllCTCFsurroundingOne[,5][i] <- as.numeric(binCorr$n[,1][2])
              }
       } else        if( distList[[i]][3] ==2) {
              if(nrow(surrCTCFOne[surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF >1e4 & surrCTCFOne$neighDownCTCF >1e4 ,])> 50) {
                     binCorr <- rcorr( surrCTCFOne[(surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF >1e4 & surrCTCFOne$neighDownCTCF >1e4) | 
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF >1e4 & surrCTCFOne$neighDownCTCF >1e4) |
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsDownCTCF >1e4 & surrCTCFOne$neighUpCTCF >1e4) |
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF >1e4 & surrCTCFOne$neighDownCTCF >1e4),]$InsExp  , 
                                          surrCTCFOne[(surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF >1e4 & surrCTCFOne$neighDownCTCF >1e4) | 
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF >1e4 & surrCTCFOne$neighDownCTCF >1e4) |
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsDownCTCF >1e4 & surrCTCFOne$neighUpCTCF >1e4) |
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF >1e4 & surrCTCFOne$neighDownCTCF >1e4),]$neighExp)  
                     binDistCorrAllvsAllCTCFsurroundingOne[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
                     binDistCorrAllvsAllCTCFsurroundingOne[,2][i] <- paste0('< ', distList[[i]][2])
                     binDistCorrAllvsAllCTCFsurroundingOne[,3][i] <- '>10kb'
                     binDistCorrAllvsAllCTCFsurroundingOne[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
                     binDistCorrAllvsAllCTCFsurroundingOne[,5][i] <- as.numeric(binCorr$n[,1][2])
              }
       } else        if( distList[[i]][3] ==3) {
              if(nrow(surrCTCFOne[(surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF <=1e4 & surrCTCFOne$neighDownCTCF >1e4) | (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF >1e4 & surrCTCFOne$neighDownCTCF <=1e4),])> 50) {
                     binCorr <- rcorr( surrCTCFOne[(surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF <=1e4 & surrCTCFOne$neighDownCTCF >1e4) | 
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF >1e4 & surrCTCFOne$neighDownCTCF <=1e4) |
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsDownCTCF >1e4 & surrCTCFOne$neighUpCTCF <=1e4) |
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF <=1e4 & surrCTCFOne$neighDownCTCF >1e4),]$InsExp,
                                                 surrCTCFOne[(surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF <=1e4 & surrCTCFOne$neighDownCTCF >1e4) | 
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF >1e4 & surrCTCFOne$neighDownCTCF <=1e4) |
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsDownCTCF >1e4 & surrCTCFOne$neighUpCTCF <=1e4) |
                                                                      (surrCTCFOne$Distance >= distList[[i]][1] & surrCTCFOne$Distance < distList[[i]][2] & surrCTCFOne$InsUpCTCF <=1e4 & surrCTCFOne$neighDownCTCF >1e4),]$neighExp)
                     binDistCorrAllvsAllCTCFsurroundingOne[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
                     binDistCorrAllvsAllCTCFsurroundingOne[,2][i] <- paste0('< ', distList[[i]][2])
                     binDistCorrAllvsAllCTCFsurroundingOne[,3][i] <- '>10kb and <=10kb'
                     binDistCorrAllvsAllCTCFsurroundingOne[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
                     binDistCorrAllvsAllCTCFsurroundingOne[,5][i] <- as.numeric(binCorr$n[,1][2])
              }
       }
       
}
binDistCorrAllvsAllCTCFsurroundingOne$nCTCF <- '>=1'

library(ggrepel)
pdf(paste0(opt$OUTDIR,'/graphs/SurroundingCTCFandCorrelation_distanceToCTCF.pdf'), height=5, width=10) 
rbind(binDistCorrAllvsAllCTCFsurrounding, binDistCorrAllvsAllCTCFsurroundingOne) %>%
       filter(!is.na(CTCFdist)) %>%
       mutate(Dist = factor(Dist, levels=unique(binDistCorrAllvsAllCTCFsurroundingOne$Dist)),
              SampleName = opt$SampleName,
              nCTCF = factor(nCTCF, levels=c(0,  '>=1'))) %>%
       ggplot(aes(x = Dist, y = Pearsons, group = CTCFdist, colour =CTCFdist)) +
              geom_line(size=1) +
               geom_point(size=4, pch=20) +
              #geom_bar(stat='identity', color='black', alpha=0.6, position = position_dodge()) +
              theme_classic() +
              ylab('Pearsons correlation') +
              facet_wrap(~nCTCF) +
              scale_colour_brewer(palette='Set2') +
              scale_fill_brewer(palette='Set2') +
             # geom_text_repel(aes(label=Pvalue), point.padding = 0.5)+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All vs all pairwise comparison of neighbours (< 20)\nand number of CTCF sites in between')+
              xlab('Distance bins')+
              ylim(0,0.5)
dev.off()



#####
#CTCF in between limited to 500kb and quantiles
#####
cat('CTCF in between limited to 500kb and >=3\n')


distList <- list(c(0, 1e5, 0), c(1e5, 2e5, 0), c(2e5, 3e5, 0), c(3e5, 4e5, 0), c(4e5, 5e5, 0),
              c(0, 1e5, 1), c(1e5, 2e5, 1), c(2e5, 3e5, 1), c(3e5, 4e5, 1), c(4e5, 5e5, 1), 
              c(0, 1e5, 2), c(1e5, 2e5, 2), c(2e5, 3e5, 2), c(3e5, 4e5, 2), c(4e5, 5e5, 2), 
              c(0,1e5, 3), c(1e5, 2e5, 3), c(2e5, 3e5, 3), c(3e5, 4e5, 3), c(4e5, 5e5, 3))
binDistCorrAllvsAllCTCFMultiple <- data.frame(matrix(ncol=5,nrow=length(distList))) 
colnames(binDistCorrAllvsAllCTCFMultiple) <- c('Pearsons', 'Dist', 'CTCF', 'Pvalue', 'n')

for (i in 1:length(distList)) { 
       if( distList[[i]][3] <3 & nrow( allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF == distList[[i]][3] ,]) >=50) {
              binCorr <- rcorr( allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF == distList[[i]][3] ,]$Expression  , allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF == distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorrAllvsAllCTCFMultiple[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorrAllvsAllCTCFMultiple[,3][i] <- distList[[i]][3]
              binDistCorrAllvsAllCTCFMultiple[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,5][i] <- as.numeric(binCorr$n[,1][2])
       } else if(distList[[i]][3] ==3 & nrow(allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF >= distList[[i]][3] ,])>50) {
              binCorr <- rcorr( allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF >= distList[[i]][3] ,]$Expression  , allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF >= distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorrAllvsAllCTCFMultiple[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorrAllvsAllCTCFMultiple[,3][i] <- paste0('>= ',distList[[i]][3])
              binDistCorrAllvsAllCTCFMultiple[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,5][i] <- as.numeric(binCorr$n[,1][2])
       }
}

binDistCorrAllvsAllCTCFMultiple <- binDistCorrAllvsAllCTCFMultiple[!is.na(binDistCorrAllvsAllCTCFMultiple$Dist),]
#binDistCorrAllvsAllCTCFMultiple <- binDistCorrAllvsAllCTCFMultiple[binDistCorrAllvsAllCTCFMultiple$Pvalue<=0.05,]

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_BinnedDistance_CTCF_allPairwise_upTOfour_500bk.pdf'), height=5, width=6) 
binDistCorrAllvsAllCTCFMultiple %>%
       mutate(Dist = factor(Dist, levels=unique(binDistCorrAllvsAllCTCFMultiple$Dist)),
              SampleName = opt$SampleName) %>%
             # CTCF = factor(CTCF, levels=c(0,  '>= 1'))) %>%
       ggplot(aes(x = Dist, y = Pearsons, group = CTCF, colour =CTCF)) +
              geom_line(size=1) +
               geom_point(size=4, pch=20) +
              #geom_bar(stat='identity', color='black', alpha=0.6, position = position_dodge()) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_colour_brewer(palette='Set2') +
              scale_fill_brewer(palette='Set2') +
              #geom_text_repel(aes(label=Pvalue), point.padding = 0.5)+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All vs all pairwise comparison of neighbours (< 20)\nand number of CTCF sites in between')+
              xlab('Distance bins')+
              ylim(0,0.5)
dev.off()




cat('CTCF in between limited to 500kb and >=4 CTCF, quantiles\n')

CTCF500kb <- allCTCF[allCTCF$Distance<=5e5,]
CTCF500kb$Decile <- ntile(CTCF500kb$Distance, 10)


distDecile <- quantile(CTCF500kb$Distance, prob = seq(0, 1, length = 11), type = 5)
distDecile <- rep(distDecile[2:11], 4)
distList <- list(c(1, 0), c( 2, 0), c( 3, 0), c( 4, 0), c( 5, 0), c( 6, 0), c( 7, 0), c( 8, 0), c( 9, 0), c( 10, 0),
              c(1, 1), c( 2, 1), c( 3, 1), c( 4, 1), c( 5, 1), c( 6, 1), c( 7, 1), c( 8, 1), c( 9, 1), c( 10, 1),
              c(1, 2), c( 2, 2), c( 3, 2), c( 4, 2), c( 5, 2), c( 6, 2), c( 7, 2), c( 8, 2), c( 9, 2), c( 10, 2),
              c(1, 3), c( 2, 3), c( 3, 3), c( 4, 3), c( 5, 3), c( 6, 3), c( 7, 3), c( 8, 3), c( 9, 3), c( 10, 3))

binDistCorrAllvsAllCTCF500kb<- data.frame(matrix(ncol=6,nrow=length(distList))) 
colnames(binDistCorrAllvsAllCTCF500kb) <- c('Pearsons', 'Decile', 'Dist', 'CTCF', 'Pvalue', 'n')


for (i in 1:length(distList)) {
       if(distList[[i]][2]<3 & nrow(CTCF500kb[CTCF500kb$Decile==distList[[i]][1] & CTCF500kb$nCTCF==distList[[i]][2],])>50 ) {
              binCorr <- rcorr(CTCF500kb[CTCF500kb$Decile==distList[[i]][1] & CTCF500kb$nCTCF==distList[[i]][2],]$Expression, CTCF500kb[CTCF500kb$Decile==distList[[i]][1] & CTCF500kb$nCTCF==distList[[i]][2],]$NeighbourExpression)
              binDistCorrAllvsAllCTCF500kb[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCTCF500kb[,2][i] <- names(distDecile)[i]
              binDistCorrAllvsAllCTCF500kb[,3][i] <-  signif(as.numeric(distDecile)[i], 2)
              binDistCorrAllvsAllCTCF500kb[,4][i] <- distList[[i]][2]
              binDistCorrAllvsAllCTCF500kb[,5][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCTCF500kb[,6][i] <- as.numeric(binCorr$n[,1][2])
       } else if(distList[[i]][2]==3 & nrow(CTCF500kb[CTCF500kb$Decile==distList[[i]][1] & CTCF500kb$nCTCF>=distList[[i]][2],])>50 ) {
              binCorr <- rcorr(CTCF500kb[CTCF500kb$Decile==distList[[i]][1] & CTCF500kb$nCTCF>=distList[[i]][2],]$Expression, CTCF500kb[CTCF500kb$Decile==distList[[i]][1] & CTCF500kb$nCTCF>=distList[[i]][2],]$NeighbourExpression)
              binDistCorrAllvsAllCTCF500kb[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCTCF500kb[,2][i] <- names(distDecile)[i]
              binDistCorrAllvsAllCTCF500kb[,3][i] <-  signif(as.numeric(distDecile)[i], 2)
              binDistCorrAllvsAllCTCF500kb[,4][i] <- paste0('> ', distList[[i]][2])
              binDistCorrAllvsAllCTCF500kb[,5][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCTCF500kb[,6][i] <- as.numeric(binCorr$n[,1][2])
       }
}
       
binDistCorrAllvsAllCTCF500kb <- binDistCorrAllvsAllCTCF500kb[!is.na(binDistCorrAllvsAllCTCF500kb$Dist),]

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_BinnedDistance_CTCF_allPairwise_500bk_upToFour_quantiles.pdf'), height=5, width=6) 
binDistCorrAllvsAllCTCF500kb %>%
       mutate(Dist = factor(Dist, levels=unique(binDistCorrAllvsAllCTCF500kb$Dist)),
              SampleName = opt$SampleName) %>%
             # CTCF = factor(CTCF, levels=c(0,  '>= 1'))) %>%
       ggplot(aes(x = Dist, y = Pearsons, group = CTCF, colour =CTCF)) +
              geom_line(size=1) +
               geom_point(size=4, pch=20) +
              #geom_bar(stat='identity', color='black', alpha=0.6, position = position_dodge()) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_colour_brewer(palette='Set2') +
              scale_fill_brewer(palette='Set2') +
              #geom_text_repel(aes(label=Pvalue), point.padding = 0.5)+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All vs all pairwise comparison of neighbours (< 20)\nand number of CTCF sites in between')+
              xlab('Distance bins')+
              ylim(0,0.55)
dev.off()



cat('CTCF in between limited to 500kb and >=1 CTCF\n')
distList <- list(c(0, 1e5, 0), c(1e5, 2e5, 0), c(2e5, 3e5, 0), c(3e5, 4e5, 0), c(4e5, 5e5, 0),
              c(0, 1e5, 1), c(1e5, 2e5, 1), c(2e5, 3e5, 1), c(3e5, 4e5, 1), c(4e5, 5e5, 1))
binDistCorrAllvsAllCTCFMultiple <- data.frame(matrix(ncol=5,nrow=length(distList))) 
colnames(binDistCorrAllvsAllCTCFMultiple) <- c('Pearsons', 'Dist', 'CTCF', 'Pvalue', 'n')

for (i in 1:length(distList)) { 
       if( distList[[i]][3] <1 & nrow( allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF == distList[[i]][3] ,]) >=50) {
              binCorr <- rcorr( allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF == distList[[i]][3] ,]$Expression  , allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF == distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorrAllvsAllCTCFMultiple[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorrAllvsAllCTCFMultiple[,3][i] <- distList[[i]][3]
              binDistCorrAllvsAllCTCFMultiple[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,5][i] <- as.numeric(binCorr$n[,1][2])
       } else if(distList[[i]][3] ==1 & nrow(allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF >= distList[[i]][3] ,])>50) {
              binCorr <- rcorr( allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF >= distList[[i]][3] ,]$Expression  , allCTCF[allCTCF$Distance >= distList[[i]][1] & allCTCF$Distance < distList[[i]][2] & allCTCF$nCTCF >= distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorrAllvsAllCTCFMultiple[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorrAllvsAllCTCFMultiple[,3][i] <- paste0('>= ',distList[[i]][3])
              binDistCorrAllvsAllCTCFMultiple[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCTCFMultiple[,5][i] <- as.numeric(binCorr$n[,1][2])
       }
}

binDistCorrAllvsAllCTCFMultiple <- binDistCorrAllvsAllCTCFMultiple[!is.na(binDistCorrAllvsAllCTCFMultiple$Dist),]
#binDistCorrAllvsAllCTCFMultiple <- binDistCorrAllvsAllCTCFMultiple[binDistCorrAllvsAllCTCFMultiple$Pvalue<=0.05,]

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_BinnedDistance_CTCF_allPairwise_500bk_.pdf'), height=5, width=6) 
binDistCorrAllvsAllCTCFMultiple %>%
       mutate(Dist = factor(Dist, levels=unique(binDistCorrAllvsAllCTCFMultiple$Dist)),
              SampleName = opt$SampleName) %>%
             # CTCF = factor(CTCF, levels=c(0,  '>= 1'))) %>%
       ggplot(aes(x = Dist, y = Pearsons, group = CTCF, colour =CTCF)) +
              geom_line(size=1) +
               geom_point(size=4, pch=20) +
              #geom_bar(stat='identity', color='black', alpha=0.6, position = position_dodge()) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_colour_brewer(palette='Set2') +
              scale_fill_brewer(palette='Set2') +
              #geom_text_repel(aes(label=Pvalue), point.padding = 0.5)+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All vs all pairwise comparison of neighbours (< 20)\nand number of CTCF sites in between')+
              xlab('Distance bins')+
              ylim(0,0.5)
dev.off()


#####
#CTCF and quantiles up to 500kb
#####

cat('CTCF in between limited to 500kb and >=1 CTCF, quantiles\n')

CTCF500kb <- allCTCF[allCTCF$Distance<=5e5,]
CTCF500kb$Decile <- ntile(CTCF500kb$Distance, 10)


distDecile <- quantile(CTCF500kb$Distance, prob = seq(0, 1, length = 11), type = 5)
distDecile <- rep(distDecile[2:11], 2)
distList <- list(c(1, 0), c( 2, 0), c( 3, 0), c( 4, 0), c( 5, 0), c( 6, 0), c( 7, 0), c( 8, 0), c( 9, 0), c( 10, 0),
              c(1, 1), c( 2, 1), c( 3, 1), c( 4, 1), c( 5, 1), c( 6, 1), c( 7, 1), c( 8, 1), c( 9, 1), c( 10, 1))

binDistCorrAllvsAllCTCF500kb<- data.frame(matrix(ncol=6,nrow=length(distList))) 
colnames(binDistCorrAllvsAllCTCF500kb) <- c('Pearsons', 'Decile', 'Dist', 'CTCF', 'Pvalue', 'n')


for (i in 1:length(distList)) {
       if(distList[[i]][2]==0 & nrow(CTCF500kb[CTCF500kb$Decile==distList[[i]][1] & CTCF500kb$nCTCF==distList[[i]][2],])>50 ) {
              binCorr <- rcorr(CTCF500kb[CTCF500kb$Decile==distList[[i]][1] & CTCF500kb$nCTCF==distList[[i]][2],]$Expression, CTCF500kb[CTCF500kb$Decile==distList[[i]][1] & CTCF500kb$nCTCF==distList[[i]][2],]$NeighbourExpression)
              binDistCorrAllvsAllCTCF500kb[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCTCF500kb[,2][i] <- names(distDecile)[i]
              binDistCorrAllvsAllCTCF500kb[,3][i] <-  signif(as.numeric(distDecile)[i], 2)
              binDistCorrAllvsAllCTCF500kb[,4][i] <- distList[[i]][2]
              binDistCorrAllvsAllCTCF500kb[,5][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCTCF500kb[,6][i] <- as.numeric(binCorr$n[,1][2])
       } else if(distList[[i]][2]==1 & nrow(CTCF500kb[CTCF500kb$Decile==distList[[i]][1] & CTCF500kb$nCTCF>=distList[[i]][2],])>50 ) {
              binCorr <- rcorr(CTCF500kb[CTCF500kb$Decile==distList[[i]][1] & CTCF500kb$nCTCF>=distList[[i]][2],]$Expression, CTCF500kb[CTCF500kb$Decile==distList[[i]][1] & CTCF500kb$nCTCF>=distList[[i]][2],]$NeighbourExpression)
              binDistCorrAllvsAllCTCF500kb[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCTCF500kb[,2][i] <- names(distDecile)[i]
              binDistCorrAllvsAllCTCF500kb[,3][i] <-  signif(as.numeric(distDecile)[i], 2)
              binDistCorrAllvsAllCTCF500kb[,4][i] <- paste0('> ', distList[[i]][2])
              binDistCorrAllvsAllCTCF500kb[,5][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCTCF500kb[,6][i] <- as.numeric(binCorr$n[,1][2])
       }
}
       
binDistCorrAllvsAllCTCF500kb <- binDistCorrAllvsAllCTCF500kb[!is.na(binDistCorrAllvsAllCTCF500kb$Dist),]

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_BinnedDistance_CTCF_allPairwise_500bk_quantiles.pdf'), height=5, width=6) 
binDistCorrAllvsAllCTCF500kb %>%
       mutate(Dist = factor(Dist, levels=unique(binDistCorrAllvsAllCTCF500kb$Dist)),
              SampleName = opt$SampleName) %>%
             # CTCF = factor(CTCF, levels=c(0,  '>= 1'))) %>%
       ggplot(aes(x = Dist, y = Pearsons, group = CTCF, colour =CTCF)) +
              geom_line(size=1) +
               geom_point(size=4, pch=20) +
              #geom_bar(stat='identity', color='black', alpha=0.6, position = position_dodge()) +
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_colour_brewer(palette='Set2') +
              scale_fill_brewer(palette='Set2') +
              #geom_text_repel(aes(label=Pvalue), point.padding = 0.5)+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All vs all pairwise comparison of neighbours (< 20)\nand number of CTCF sites in between')+
              xlab('Distance bins')+
              ylim(0,0.5)
dev.off()



#####
#CTCF MCV of 1 CTCF between insertions
#####
cat('CTCF MCV of 1 CTCF between insertions\n')

mcvCTCF <- read.table(paste0(opt$OUTDIR,'/AllNeighbour_CTCF_MCV.bed'), stringsAsFactors=F)
colnames(mcvCTCF) <- c('chrom', 'chromStart', 'chromEnd', 'insBC', 'neighBC', 'Distance', 'MCV', 'nCTCF')
mcvCTCF$insExp <- normdataiPCR[match(mcvCTCF$insBC, normdataiPCR$BC),]$expression
mcvCTCF$neighExp <- normdataiPCR[match(mcvCTCF$neighBC, normdataiPCR$BC),]$expression

mcvCTCF$CTCFspecificity <- NA
mcvCTCF$CTCFspecificity[mcvCTCF$MCV==0] <- 'NoCTCF'
mcvCTCF$CTCFspecificity[mcvCTCF$MCV==1] <- 'Specific'
mcvCTCF$CTCFspecificity[mcvCTCF$MCV>1 & mcvCTCF$MCV<32] <- 'nonSpecific'
mcvCTCF$CTCFspecificity[mcvCTCF$MCV>=32] <- 'Ubiquitous (>90%)'

mcvCTCF <- mcvCTCF[mcvCTCF$Distance <=1e5,]

binNoCTCF <- rcorr(mcvCTCF[mcvCTCF$CTCFspecificity=='NoCTCF',]$insExp, mcvCTCF[mcvCTCF$CTCFspecificity=='NoCTCF',]$neighExp)
binSpecCTCF <- rcorr(mcvCTCF[mcvCTCF$CTCFspecificity=='Specific',]$insExp, mcvCTCF[mcvCTCF$CTCFspecificity=='Specific',]$neighExp)
binNonSpecCTCF <- rcorr(mcvCTCF[mcvCTCF$CTCFspecificity=='nonSpecific',]$insExp, mcvCTCF[mcvCTCF$CTCFspecificity=='nonSpecific',]$neighExp)
binUbiCTCF <- rcorr(mcvCTCF[mcvCTCF$CTCFspecificity=='Ubiquitous (>90%)',]$insExp, mcvCTCF[mcvCTCF$CTCFspecificity=='Ubiquitous (>90%)',]$neighExp)

CTCFscpecificity <- data.frame(CTCFSpecificity = c('NoCTCF', 'K562Specific', 'nonSpecific', 'Ubiquitous (>90%)'), 
              Pearsons =c(signif(as.numeric(binNoCTCF$r[1,2]),2), 
                            signif(as.numeric(binSpecCTCF$r[1,2]),2),
                            signif(as.numeric(binNonSpecCTCF$r[1,2]),2), 
                            signif(as.numeric(binUbiCTCF$r[1,2]),2)))


pdf(paste0(opt$OUTDIR,'/graphs/Neighbour_CTCFspecificity.pdf'), height=5, width=5) 
CTCFscpecificity %>%
       mutate(SampleName = opt$SampleName,
              CTCFSpecificity = factor(CTCFSpecificity, levels=c('NoCTCF', 'K562Specific', 'nonSpecific', 'Ubiquitous (>90%)'))) %>%
       ggplot(aes(x=CTCFSpecificity, y =Pearsons, fill = SampleName)) +
              geom_bar(stat='identity', colour='black', alpha=0.6)+
              theme_classic() +
              ylab('Pearsons correlation') +
              scale_fill_manual(values=sampleColor) +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12), legend.position='none') +
              ggtitle(opt$SampleName, 'MCV score of CTCF site between insertions\n(if only one exists) within 100kb')+
              xlab('CTCF specificity between insertions')+
              ylim(0,0.5)
dev.off()
#####
#A/B compartments
#####
cat('Correlation between A/B compartments\n')

neighAB <- neighCTCF
neighAB$InsComp <- normdataiPCR[match(neighAB$BC, normdataiPCR$BC),]$ABcomp
neighAB$NeighComp <- normdataiPCR[match(neighAB$neighBC, normdataiPCR$BC),]$ABcomp
neighAB <- neighAB[!is.na(neighAB$InsComp) & !is.na(neighAB$NeighComp),]

binCompCorr <- data.frame(matrix(ncol=3,nrow=4)) 
colnames(binCompCorr) <- c('Pearsons', 'Comp', 'Pvalue')

binCorrAA <- rcorr( neighAB[neighAB$InsComp == 'A' & neighAB$NeighComp == 'A',]$Expression  , neighAB[neighAB$InsComp  == 'A' & neighAB$NeighComp == 'A',]$NeighbourExpression) 
binCorrBB <- rcorr( neighAB[neighAB$InsComp == 'B' & neighAB$NeighComp == 'B',]$Expression  , neighAB[neighAB$InsComp  == 'B' & neighAB$NeighComp == 'B',]$NeighbourExpression)   
binCorrAB <- rcorr(neighAB[neighAB$InsComp == 'A' & neighAB$NeighComp == 'B' | neighAB$InsComp == 'B' & neighAB$NeighComp == 'A',]$Expression , 
                     neighAB[neighAB$InsComp == 'A' & neighAB$NeighComp == 'B' | neighAB$InsComp == 'B' & neighAB$NeighComp == 'A',]$NeighbourExpression)
                     
CompCorr <- rbind(data.frame(Pearsons=signif(as.numeric(binCorrAA$r[,1][2]),3), Compartment= 'AA', Pvalue = signif(as.numeric(binCorrAA$P[,1][2]),3), n=as.numeric(binCorrAA$n[1,][1])),
       data.frame(Pearsons=signif(as.numeric(binCorrAB$r[,1][2]),3), Compartment= 'AB', Pvalue = signif(as.numeric(binCorrAB$P[,1][2]),3), n=as.numeric(binCorrAB$n[1,][1])),
       data.frame(Pearsons=signif(as.numeric(binCorrBB$r[,1][2]),3), Compartment= 'BB', Pvalue = signif(as.numeric(binCorrBB$P[,1][2]),3), n=as.numeric(binCorrBB$n[1,][1])))


rects <- data.frame(xstart = c(0.3,2), xend = c(2,3.7), Compartment=c('AA', 'BB'))

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_Compartments.pdf'), height=5, width=4) 
CompCorr %>%
       mutate(Compartment = factor(Compartment, levels=c('AA', 'AB', 'BB')),
              SampleName= opt$SampleName) %>%
              ggplot() +
              geom_bar(stat='identity', color='black', alpha=0.6, aes( x = Compartment, y = Pearsons), fill=sampleColor) +
              geom_rect(data=rects, mapping=aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill=Compartment), alpha=0.25)+
              #geom_bar(stat='identity', color='black', alpha=0.6, aes( x = Compartment, y = Pearsons), fill=sampleColor)+
              theme_classic() +
              scale_fill_brewer(palette='Set1') +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'Correlation between nearest\nneighbour and A/B compartments')+
              theme(legend.position='none')+
              ylab('Pearsons correlation') +
              xlab('Compartments of nearest neighbours')+
              ylim(0, 0.5)+
              annotate('text', label='A', x='AA', y=0.4, size=10)+ 
              annotate('text', label='B', x='BB', y=0.4, size=10)
dev.off()


#####
#Distance between compartment neighbours
#####
cat('Distance between A/B comaprtments\n')

neighAB$Compartment <-''
neighAB[neighAB$InsComp == 'A' & neighAB$NeighComp == 'A',]$Compartment <-'AA'
neighAB[neighAB$InsComp == 'B' & neighAB$NeighComp == 'B',]$Compartment <-'BB'
neighAB[neighAB$InsComp != neighAB$NeighComp,]$Compartment <-'AB'
rects <- data.frame(xstart = c(0.3,2), xend = c(2,3.7), Compartment=c('AA', 'BB'))

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_Compartments_Distance.pdf'), height=5, width=4) 
neighAB %>%
       mutate(Compartment = factor(Compartment, levels=c('AA', 'AB', 'BB')),
              SampleName= opt$SampleName) %>%
              ggplot() +
              geom_boxplot(outlier.shape=NA, alpha=0.6, aes( x = Compartment, y = log10(Distance+1)), fill=sampleColor) +
              geom_rect(data=rects, mapping=aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill=Compartment), alpha=0.25)+
              #geom_bar(stat='identity', color='black', alpha=0.6, aes( x = Compartment, y = Pearsons), fill=sampleColor)+
              theme_classic() +
              scale_fill_brewer(palette='Set1') +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'Distance between neighbour\nand A/B compartments')+
              theme(legend.position='none')+
              ylab('Distante to neighbour') +
              xlab('Compartments of nearest neighbours')+
              annotate('text', label='A', x='AA', y=2.5, size=10)+ 
              annotate('text', label='B', x='BB', y=2.5, size=10)
dev.off()

#####
#Distance between compartments and expression
#####
cat('Distance between A/B comaprtments and expression\n')

distComp <- normdataiPCR
distComp <- distComp[!is.na(distComp$ABcomp),]
distComp$ABcomp_Adist <- abs(distComp$ABcomp_Adist)
distComp$ABcomp_Bdist <- abs(distComp$ABcomp_Bdist)*-1
distComp$ABcomp_Adist[distComp$ABcomp_Adist==0] <- distComp$ABcomp_Bdist[distComp$ABcomp_Adist==0]

distComp <- subset(distComp, select=c(BC, ABcomp, ABcomp_Adist))
colnames(distComp) <- c('BC', 'Compartment', 'Distance')
distComp$Expression <- normdataiPCR[match(distComp$BC, normdataiPCR$BC),]$expression

rects <- data.frame(xstart = c(-5e5,0), xend = c(0,5e5), Compartment=c('A', 'B'))

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_CompartmentsExpression.pdf'), height=5, width=5) 
distComp %>%
       filter(Distance <5e5 & Distance> -5e5) %>%
       ggplot()+
              #geom_point()+
              geom_smooth(aes(x=Distance, y=Expression), method='loess', color=sampleColor) +
              geom_rect(data=rects, mapping=aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill=Compartment), alpha=0.25)+
              #geom_bar(stat='identity', color='black', alpha=0.6, aes( x = Compartment, y = Pearsons), fill=sampleColor)+
              theme_classic() +
              scale_fill_brewer(palette='Set1') +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'Compartment expression and distance\nto compartment boundary')+
              theme(legend.position='none')+
              ylab('Reporter activity') +
              xlab('Distance to compartment boundary')+
              annotate('text', label='A', x=-2.5e5, y=0, size=10)+ 
              annotate('text', label='B', x=2.5e5, y=0, size=10)
dev.off()



######
#A/B compartments for all 20 neighbouring insertions
######
AB20 <- neighLoc


AB20$insComp <- normdataiPCR[match(AB20$V1, normdataiPCR$BC),]$ABcompNumber
AB20$neighComp <- normdataiPCR[match(AB20$V2, normdataiPCR$BC),]$ABcompNumber

AB20 <- AB20[!is.na(AB20$insComp) & !is.na(AB20$neighComp),]
AB20$insCompVal <- as.numeric(gsub('A_|B_', '', AB20$insComp))
AB20$neighCompVal <- as.numeric(gsub('A_|B_', '', AB20$neighComp))

AB20$insCompType <- gsub('_.*', '', AB20$insComp)
AB20$neighCompType <- gsub('_.*', '', AB20$neighComp)

AB20 <- AB20[AB20$insCompVal - AB20$neighCompVal <=1,]

AB20corr <- data.frame(Compartment=c('AA', 'BB', 'AB'),
              Pearsons=c(as.numeric(rcorr(AB20[AB20$insComp == AB20$neighComp & AB20$insCompType=='A' ,]$Expression, AB20[AB20$insComp == AB20$neighComp & AB20$insCompType=='A' ,]$NeighbourExpression)$r[,1][2]),
                     as.numeric(rcorr(AB20[AB20$insComp == AB20$neighComp & AB20$insCompType=='B' ,]$Expression, AB20[AB20$insComp == AB20$neighComp & AB20$insCompType=='B' ,]$NeighbourExpression)$r[,1][2]),
                     as.numeric(rcorr(AB20[AB20$insComp != AB20$neighComp  ,]$Expression, AB20[AB20$insComp != AB20$neighComp  ,]$NeighbourExpression)$r[,1][2])))



rects <- data.frame(xstart = c(0.3,2), xend = c(2,3.7), Compartment=c('AA', 'BB'))

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_Compartments_allvsall.pdf'), height=5, width=5) 
AB20corr %>%
       mutate(Compartment = factor(Compartment, levels=c('AA', 'AB', 'BB')),
              SampleName= opt$SampleName) %>%
              ggplot() +
              geom_bar(stat='identity', color='black', alpha=0.6, aes( x = Compartment, y = Pearsons), fill=sampleColor) +
              geom_rect(data=rects, mapping=aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill=Compartment), alpha=0.25)+
              #geom_bar(stat='identity', color='black', alpha=0.6, aes( x = Compartment, y = Pearsons), fill=sampleColor)+
              theme_classic() +
              scale_fill_brewer(palette='Set1') +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All vs all correlation betweennearest\nneighbours and A/B compartments')+
              theme(legend.position='none')+
              ylab('Pearsons correlation') +
              geom_hline(yintercept=0) +
              xlab('Compartments of nearest neighbours')+
              ylim(-0.02, 0.5)+
              annotate('text', label='A', x='AA', y=0.4, size=10)+ 
              annotate('text', label='B', x='BB', y=0.4, size=10)
dev.off()




distList <- list(c(0, 1e5), c(1e5, 2e5), c(2e5, 3e5), c(3e5, 4e5), c(4e5, 5e5), c(5e5, 6e5), c(6e5, 7e5),  c(7e5, 8e5), c(8e5, 9e5),  c(9e5, 1e6))
binDistCorrAllvsAllAB20BB<- data.frame(matrix(ncol=5,nrow=length(distList))) 
colnames(binDistCorrAllvsAllAB20BB) <- c('Pearsons', 'Dist', 'Compartment', 'Pvalue', 'n')


for (i in 1:length(distList)) {
       binCorr <- rcorr(AB20[AB20$insComp == AB20$neighComp & AB20$insCompType=='B' & AB20$V3 >= distList[[i]][1] & AB20$V3 < distList[[i]][2]  ,]$Expression, 
              AB20[AB20$insComp == AB20$neighComp & AB20$insCompType=='B' & AB20$V3 >= distList[[i]][1] & AB20$V3 < distList[[i]][2]  ,]$NeighbourExpression)
      binDistCorrAllvsAllAB20BB[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
      binDistCorrAllvsAllAB20BB[,2][i] <- paste0('< ', distList[[i]][2])
      binDistCorrAllvsAllAB20BB[,3][i] <- "BB"
      binDistCorrAllvsAllAB20BB[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
      binDistCorrAllvsAllAB20BB[,5][i] <- as.numeric(binCorr$n[,1][2])
}

distList <- list(c(0,  1e5), c(1e5, 2e5), c(2e5, 3e5), c(3e5, 4e5), c(4e5, 5e5), c(5e5, 6e5), c(6e5, 7e5),  c(7e5, 8e5), c(8e5, 9e5),  c(9e5, 1e6))
binDistCorrAllvsAllAB20AA<- data.frame(matrix(ncol=5,nrow=length(distList))) 
colnames(binDistCorrAllvsAllAB20AA) <- c('Pearsons', 'Dist', 'Compartment', 'Pvalue', 'n')


for (i in 1:length(distList)) {
       binCorr <- rcorr(AB20[AB20$insComp == AB20$neighComp & AB20$insCompType=='A' & AB20$V3 >= distList[[i]][1] & AB20$V3 < distList[[i]][2]  ,]$Expression, 
              AB20[AB20$insComp == AB20$neighComp & AB20$insCompType=='A' & AB20$V3 >= distList[[i]][1] & AB20$V3 < distList[[i]][2]  ,]$NeighbourExpression)
       binDistCorrAllvsAllAB20AA[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
       binDistCorrAllvsAllAB20AA[,2][i] <- paste0('< ', distList[[i]][2])
       binDistCorrAllvsAllAB20AA[,3][i] <- "AA"
       binDistCorrAllvsAllAB20AA[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
       binDistCorrAllvsAllAB20AA[,5][i] <- as.numeric(binCorr$n[,1][2])
}


distList <- list(c(0,  1e5), c(1e5, 2e5), c(2e5, 3e5), c(3e5, 4e5), c(4e5, 5e5), c(5e5, 6e5), c(6e5, 7e5),  c(7e5, 8e5), c(8e5, 9e5),  c(9e5, 1e6))
binDistCorrAllvsAllAB20AB<- data.frame(matrix(ncol=5,nrow=length(distList))) 
colnames(binDistCorrAllvsAllAB20AB) <- c('Pearsons', 'Dist', 'Compartment', 'Pvalue', 'n')


for (i in 1:length(distList)) {
       binCorr <- rcorr(AB20[AB20$insComp != AB20$neighComp  & AB20$V3 >= distList[[i]][1] & AB20$V3 < distList[[i]][2]  ,]$Expression, 
              AB20[AB20$insComp != AB20$neighComp & AB20$V3 >= distList[[i]][1] & AB20$V3 < distList[[i]][2]  ,]$NeighbourExpression)
       binDistCorrAllvsAllAB20AB[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
       binDistCorrAllvsAllAB20AB[,2][i] <- paste0('< ', distList[[i]][2])
       binDistCorrAllvsAllAB20AB[,3][i] <- "AB"
       binDistCorrAllvsAllAB20AB[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
       binDistCorrAllvsAllAB20AB[,5][i] <- as.numeric(binCorr$n[,1][2])
}

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_BinnedDistance_ABcompartment_allPairwis.pdf'), height=5, width=7) 
rbind(binDistCorrAllvsAllAB20BB, binDistCorrAllvsAllAB20AA, binDistCorrAllvsAllAB20AB) %>%
              filter(!is.na(Dist)) %>%
       mutate(Dist = factor(Dist, levels=unique(binDistCorrAllvsAllAB20BB$Dist)),
              SampleName = opt$SampleName,
              Compartment = factor(Compartment, levels=c('AA', 'BB', 'AB'))) %>%
       ggplot(aes(x = Dist, y = Pearsons, group = Compartment, colour =Compartment)) +
              geom_line(size=1) +
               geom_point(size=4, pch=20) +
              #geom_bar(stat='identity', color='black', alpha=0.6, position = position_dodge()) +
              theme_classic() +
              ylab('Pearsons correlation') +
#              facet_wrap(~nDHS) +
              scale_colour_brewer(palette='Set1') +
              scale_fill_brewer(palette='Set1') +
             # geom_text_repel(aes(label=Pvalue), point.padding = 0.5)+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All vs all pairwise comparison of neighbours (< 20)\nand AB compartment')+
              xlab('Distance bins')+
              ylim(-0.1,0.5)
dev.off()
#####
#nDHS in between insertions
####
cat('Distance to neighbours and DHS between\n')
neighCTCF20 <- read.table(paste0(opt$OUTDIR,'/Neighbour_20neigh_CTCF.bed'), sep='\t', header=T, stringsAsFactors=F)
neighCTCF20$Expression <- normdataiPCR[match(neighCTCF20$BC, normdataiPCR$BC),]$expression
neighCTCF20$NeighbourExpression <- normdataiPCR[match(neighCTCF20$neighBC, normdataiPCR$BC),]$expression

distList <- list(c(0, 1e4, 0), c(1e4, 5e4, 0), c(5e4, 1e5, 0), c(1e5, 2e5, 0), c(2e5, 3e5, 0), c(3e5, 4e5, 0), c(4e5, 5e5, 0), c(5e5, 1e6, 0),
               c(0, 1e4, 1), c(1e4, 5e4, 1), c(5e4, 1e5, 1), c(1e5, 2e5, 1), c(2e5, 3e5, 1), c(3e5, 4e5, 1), c(4e5, 5e5, 1), c(5e5, 1e6, 1))
binDistCorrAllvsAllDHS <- data.frame(matrix(ncol=5,nrow=length(distList))) 
colnames(binDistCorrAllvsAllDHS) <- c('Pearsons', 'Dist', 'nDHS', 'Pvalue', 'n')

for (i in 1:length(distList)) { 
       if( distList[[i]][3] <1 & nrow(neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nDHS == distList[[i]][3] ,])>50) {
              binCorr <- rcorr( neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nDHS == distList[[i]][3] ,]$Expression  , neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nDHS == distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorrAllvsAllDHS[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllDHS[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorrAllvsAllDHS[,3][i] <- distList[[i]][3]
              binDistCorrAllvsAllDHS[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllDHS[,5][i] <- nrow(neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nDHS == distList[[i]][3] ,])
       } else if(distList[[i]][3] ==1 & nrow(neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nDHS >= distList[[i]][3] ,]) >50) {
              binCorr <- rcorr( neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nDHS >= distList[[i]][3] ,]$Expression  , neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nDHS >= distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorrAllvsAllDHS[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllDHS[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorrAllvsAllDHS[,3][i] <- paste0('>= ',distList[[i]][3])
              binDistCorrAllvsAllDHS[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllDHS[,5][i] <- nrow(neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nDHS >= distList[[i]][3] ,])
       }
}

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_BinnedDistance_DHS_allPairwis.pdf'), height=5, width=5) 
binDistCorrAllvsAllDHS %>%
       filter(!is.na(Dist)) %>%
       mutate(Dist = factor(Dist, levels=c('< 10000', '< 50000', '< 1e+05', '< 2e+05', '< 3e+05', '< 4e+05', '< 5e+05', '< 1e+06')),
              SampleName = opt$SampleName,
              nDHS = factor(nDHS, levels=c(0,  '>= 1'))) %>%
       ggplot(aes(x = Dist, y = Pearsons, group = nDHS, colour =nDHS)) +
              geom_line(size=1) +
               geom_point(size=4, pch=20) +
              #geom_bar(stat='identity', color='black', alpha=0.6, position = position_dodge()) +
              theme_classic() +
              ylab('Pearsons correlation') +
#              facet_wrap(~nDHS) +
              scale_colour_brewer(palette='Set2') +
              scale_fill_brewer(palette='Set2') +
             # geom_text_repel(aes(label=Pvalue), point.padding = 0.5)+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All vs all pairwise comparison of neighbours (< 20)\nand number of DHS in between')+
              xlab('Distance bins')+
              ylim(0,0.5)
dev.off()

#####
#nCGI in between insertions
####
cat('Distance to neighbours and CGI between\n')

distList <- list(c(0, 1e4, 0), c(1e4, 5e4, 0), c(5e4, 1e5, 0), c(1e5, 2e5, 0), c(2e5, 3e5, 0), c(3e5, 4e5, 0), c(4e5, 5e5, 0), c(5e5, 1e6, 0),
               c(0, 1e4, 1), c(1e4, 5e4, 1), c(5e4, 1e5, 1), c(1e5, 2e5, 1), c(2e5, 3e5, 1), c(3e5, 4e5, 1), c(4e5, 5e5, 1), c(5e5, 1e6, 1))
binDistCorrAllvsAllCGI <- data.frame(matrix(ncol=5,nrow=length(distList))) 
colnames(binDistCorrAllvsAllCGI) <- c('Pearsons', 'Dist', 'nCGI', 'Pvalue', 'n')

for (i in 1:length(distList)) { 
       if( distList[[i]][3] <1 & nrow(neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nCGI == distList[[i]][3] ,])>50) {
              binCorr <- rcorr( neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nCGI == distList[[i]][3] ,]$Expression  , neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nCGI == distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorrAllvsAllCGI[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCGI[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorrAllvsAllCGI[,3][i] <- distList[[i]][3]
              binDistCorrAllvsAllCGI[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCGI[,5][i] <- nrow(neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nCGI == distList[[i]][3] ,])
       } else if(distList[[i]][3] ==1 & nrow(neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nCGI >= distList[[i]][3] ,]) >50) {
              binCorr <- rcorr( neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nCGI >= distList[[i]][3] ,]$Expression  , neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nCGI >= distList[[i]][3] ,]$NeighbourExpression)  
              binDistCorrAllvsAllCGI[,1][i] <- signif(as.numeric(binCorr$r[,1][2]),3)
              binDistCorrAllvsAllCGI[,2][i] <- paste0('< ', distList[[i]][2])
              binDistCorrAllvsAllCGI[,3][i] <- paste0('>= ',distList[[i]][3])
              binDistCorrAllvsAllCGI[,4][i] <- signif(as.numeric(binCorr$P[,1][2]),3)
              binDistCorrAllvsAllCGI[,5][i] <- nrow(neighCTCF20[neighCTCF20$Distance >= distList[[i]][1] & neighCTCF20$Distance < distList[[i]][2] & neighCTCF20$nCGI >= distList[[i]][3] ,])
       }
}

pdf(paste0(opt$OUTDIR,'/graphs/NeighbouringInsertions_BinnedDistance_CGI_allPairwise.pdf'), height=5, width=5) 
binDistCorrAllvsAllCGI %>%
       filter(!is.na(Dist)) %>%
       mutate(Dist = factor(Dist, levels=c('< 10000', '< 50000', '< 1e+05', '< 2e+05', '< 3e+05', '< 4e+05', '< 5e+05', '< 1e+06')),
              SampleName = opt$SampleName,
              nCGI = factor(nCGI, levels=c(0,  '>= 1'))) %>%
       ggplot(aes(x = Dist, y = Pearsons, group = nCGI, colour =nCGI)) +
              geom_line(size=1) +
               geom_point(size=4, pch=20) +
              #geom_bar(stat='identity', color='black', alpha=0.6, position = position_dodge()) +
              theme_classic() +
              ylab('Pearsons correlation') +
#              facet_wrap(~nCGI) +
              scale_colour_brewer(palette='Set2') +
              scale_fill_brewer(palette='Set2') +
             # geom_text_repel(aes(label=Pvalue), point.padding = 0.5)+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
              ggtitle(opt$SampleName, 'All vs all pairwise comparison of neighbours (< 20)\nand number of CGI in between')+
              xlab('Distance bins')+
              ylim(0,0.5)
dev.off()

cat('Reporter activity - nearest CTCF or nearest DHS\n')
pdf(paste0(opt$OUTDIR,'/graphs/distDHSorCTCFnearest.pdf'), height=7, width=5) 
boxplot(normdataiPCR[abs(normdataiPCR$DistToNearestCTCF) > abs(normdataiPCR$DistToNearestDHS),]$expression,
       normdataiPCR[abs(normdataiPCR$DistToNearestCTCF) < abs(normdataiPCR$DistToNearestDHS),]$expression, 
       frame=F, names=c('distCTCF>distDHS', 'distCTCF<distDHS'), ylab='Reporter activity', col='#fb6a4a', outline=FALSE, ylim=c(-2, 2),
       main=paste0(opt$SampleName,'\nDistance to CTCF and DHS\nand reporter activity')) 
text(x=c(1, 2, 1.5), y=c(1.5, 1.5, 1.8), labels=c(nrow(normdataiPCR[abs(normdataiPCR$DistToNearestCTCF) > abs(normdataiPCR$DistToNearestDHS),]),
                                   nrow(normdataiPCR[abs(normdataiPCR$DistToNearestCTCF) < abs(normdataiPCR$DistToNearestDHS),]),
                                   paste0('pvalue = ',signif(wilcox.test(normdataiPCR[abs(normdataiPCR$DistToNearestCTCF) > abs(normdataiPCR$DistToNearestDHS),]$expression,
                                                               normdataiPCR[abs(normdataiPCR$DistToNearestCTCF) < abs(normdataiPCR$DistToNearestDHS),]$expression)$p.value,3))))
dev.off()
cat('Reporter activity - nearest DHS type\n')

pdf(paste0(opt$OUTDIR,'/graphs/nearestDHStype_Activity.pdf'), height=7, width=7) 

boxplot( normdataiPCR$expression~ normdataiPCR$NearestDHStype, frame=F, ylab='Reporter activity', col= brewer.pal(n=3, name='Accent'), outline=FALSE, ylim=c(-2, 2),
       main=paste0(opt$SampleName, '\nNearest DHS type and reporter activity'))
text(y=c(1.5, 1.8), x=c(1.5, 2.5), labels=c(paste0('pvalue = ',signif(wilcox.test(normdataiPCR[normdataiPCR$NearestDHStype=='CTCF',]$expression, normdataiPCR[normdataiPCR$NearestDHStype=='Distal',]$expression)$p.value,3)),
                                          paste0('pvalue = ',signif(wilcox.test(normdataiPCR[normdataiPCR$NearestDHStype=='Promoter',]$expression, normdataiPCR[normdataiPCR$NearestDHStype=='Distal',]$expression)$p.value,3))))
dev.off()




########
#Correlation between insertions and K562 gencode genes
########
library(dplyr)
library(Hmisc)
correlateInsertionGencode <- function(data, distCol, distInterval, reporterExpCol, corrExpCol, FeatureCol) {
       listColumn <- list()
       for (colFactor in levels(factor(data[,FeatureCol]))) {
              dataFactor <- data[data[,FeatureCol] == colFactor,]
              #cat(nrow(dataFactor[abs(dataFactor[,distCol]) >= distInterval[1] & abs(dataFactor[,distCol]) < distInterval[2],]), '\n')
              if (nrow(dataFactor[abs(dataFactor[,distCol]) >= distInterval[1] & abs(dataFactor[,distCol]) < distInterval[2],]) >4) {
                     binCorr <- rcorr( dataFactor[abs(dataFactor[,distCol]) >= distInterval[1] & abs(dataFactor[,distCol]) < distInterval[2],][,reporterExpCol]  , 
                     dataFactor[abs(dataFactor[,distCol]) >= distInterval[1] & abs(dataFactor[,distCol]) < distInterval[2],][,corrExpCol])
                     listColumn[[paste0(colFactor)]] <- data.frame(Pearsons =   signif(as.numeric(binCorr$r[,1][2]),3),
                                   Dist =  paste0('< ', distInterval[2]),
                                   Pvalue = signif(as.numeric(binCorr$P[,1][2]),3),
                                   #n = nrow(dataFactor[abs(dataFactor[,distCol]) >= distInterval[1] & abs(dataFactor[,distCol]) < distInterval[2],]),
                                   n = as.numeric(binCorr$n[,1][2]),
                                   Factor = colFactor)
              }
       }
       df <- do.call(rbind, listColumn)
       return(df)
}

distList <- list(c(0, 1e3), c(1e3, 1e4), c(1e4, 5e4), c(5e4, 1e5), c(1e5, 2e5), c(2e5, 3e5), c(3e5, 4e5), c(4e5, 5e5), c(5e5, 1e6)) 

corrGene <- normdataiPCR
#
#distToGencode <- subset(corrGene , select=c(DistToTSS, DistToGencodeTSS_protLinc))
#distToGencode %>% 
#       melt() %>%
#       ggplot(aes( x = variable, y=log10(abs(value)+1), fill=variable)) +
#              geom_boxplot()


corrGene$DistToGencodeTSSDecile <- ntile(abs(corrGene$DistToGencodeTSS_protLinc), n=10)

decileDist <- list()
for (i in 1:10) {
       Decile <- corrGene[corrGene$DistToGencodeTSSDecile==i,]
       decileDist[[i]] <- c(min(abs(Decile$DistToGencodeTSS)), max(abs(Decile$DistToGencodeTSS)))
}

#decileDist <- list(c(0, 1e5))
corrGeneLoc <- list()
for (i in distList) {
       Distance <- as.character(i[2])
       corrGeneLoc[[Distance]] <- correlateInsertionGencode(data = normdataiPCR, 
                                                               distCol = "DistToGencodeTSS_protLinc", 
                                                               distInterval = i, 
                                                               reporterExpCol = "expression", 
                                                               corrExpCol = "ToGencodeTSSTPM_protLinc", 
                                                               FeatureCol="SampleName")
}
corrGeneLoc <- do.call(rbind, corrGeneLoc)

pdf(paste0(opt$OUTDIR,'/graphs/genicCorrelation_dist.pdf'), height=5, width=7) 
corrGeneLoc %>%
       filter(Pvalue <0.05) %>%
       ggplot(aes( x = Dist, y = Pearsons, colour = Factor , fill = Factor, group = Factor)) +
              geom_line()+
              geom_point(pch=21) +
              theme_classic() +
              scale_fill_manual(values = sampleColor) +
              scale_colour_manual(values = sampleColor) +
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(angle=60, hjust=1, size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))# +
              #ylim(0, 1)
dev.off()
#
#locList <- list()
#for (locIns in unique(normdataiPCR$GenicLocation)) {
#       corrLoc <- rcorr(normdataiPCR[normdataiPCR$GenicLocation==locIns,]$expression, normdataiPCR[normdataiPCR$GenicLocation==locIns,]$ToGencodeTSSTPM)
#       locList[[locIns]] <- data.frame(Factor = locIns, 
#                                   Pearsons =   signif(as.numeric(corrLoc$r[,1][2]),3),
#                                   Pvalue = signif(as.numeric(corrLoc$P[,1][2]),3),
#                                   #n = nrow(dataFactor[abs(dataFactor[,distCol]) >= distInterval[1] & abs(dataFactor[,distCol]) < distInterval[2],]),
#                                   n = as.numeric(corrLoc$n[,1][2])
#                                   )
#}
#locList <- do.call(rbind, locList)
#
#pdf(paste0(opt$OUTDIR,'/graphs/genicCorrelation.pdf'), height=5, width=7) 
#locList %>%
#       mutate(GenicLocation = factor(Factor, levels=levels(Factor))) %>%
#       ggplot(aes(x=GenicLocation, y=Pearsons))+
#              geom_bar(stat = 'identity', colour='black', fill=sampleColor, size=0.01)+
#              theme_classic() +
#              #scale_fill_manual(values=c('#f0f0f0', '#4292c6')) +
#              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, hjust=1), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#              xlab('Intertion location') +
#              ylab('Pearsons')+
#              ggtitle(opt$SampleName, 'Correlation between of insertions and nearest gene\nbased on genic location')+
#              geom_text(aes(x=GenicLocation, y=Pearsons, vjust=-1.5, label=paste0("P=",Pvalue)))+
#              geom_text(aes(x=GenicLocation, y=Pearsons, vjust=-0.3, label=paste0("n=",n)))+
#              ylim(0,1)
#dev.off()
#


##decileDist <- list(c(0, 1e5))
#corrGeneLoc <- list()
#for (i in distList) {
#       Distance <- as.character(i[2])
#       corrGeneLoc[[Distance]] <- correlateInsertionGencode(data = normdataiPCR, 
#                                                               distCol = "DistToGencodeTSS", 
#                                                               distInterval = i, 
#                                                               reporterExpCol = "expression", 
#                                                               corrExpCol = "ToGencodeTSSTPM", 
#                                                               FeatureCol="GenicLocation")
#}
#corrGeneLoc <- do.call(rbind, corrGeneLoc)
#
#corrGeneLoc <- corrGeneLoc[grep('intron|intergenic|coding', corrGeneLoc$Factor),]
#
#pdf(paste0('~/public_html/blog/20180305/genicCorrelation_genicDist.pdf'), height=5, width=7) 
#corrGeneLoc %>%
#       ggplot(aes( x = Dist, y = Pearsons, colour = Factor , fill = Factor, group = Factor)) +
#              geom_line()+
#              geom_point(pch=21) +
#              theme_classic() +
#              scale_fill_brewer(palette = "Set1") +
#              scale_colour_brewer(palette = "Set1") +
#              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(angle=60, hjust=1, size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
#              ylim(-0.5, 1) +
#              ggtitle(opt$SampleName, 'Correlation between of insertions and nearest gene\nbased on genic location and distance')
#dev.off()

#cat('Cohesin and CTCF loops\n')
#cohesin <- read.table('/tmp/AllInsertions.coords.bed', sep='\t' ,stringsAsFactors=F)
#surrCTCF$InsLoop <- cohesin[match(surrCTCF$InsBC, cohesin$V4),]$V7
#surrCTCF$neighLoop <- cohesin[match(surrCTCF$neighBC, cohesin$V4),]$V7
#
#cohesin <- surrCTCF[!is.na(surrCTCF$InsLoop) & !is.na(surrCTCF$neighLoop),]
#cohesin$InsLoopNum <- as.numeric(gsub('CohesinLoop_','', cohesin$InsLoop))
#cohesin$neighLoopNum <- as.numeric(gsub('CohesinLoop_','', cohesin$neighLoop))
#
#cohesin <- cohesin[abs(cohesin$InsLoopNum-cohesin$neighLoopNum)<=1,]
#
#sameLoop <- rcorr(cohesin[abs(cohesin$InsLoopNum-cohesin$neighLoopNum)==0,]$InsExp, cohesin[abs(cohesin$InsLoopNum-cohesin$neighLoopNum)==0,]$neighExp)
#neighLoop <- rcorr(cohesin[abs(cohesin$InsLoopNum-cohesin$neighLoopNum)==1,]$InsExp, cohesin[abs(cohesin$InsLoopNum-cohesin$neighLoopNum)==1,]$neighExp)
#rcorr(cohesin[abs(cohesin$InsLoopNum-cohesin$neighLoopNum)==3,]$InsExp, cohesin[abs(cohesin$InsLoopNum-cohesin$neighLoopNum)==3,]$neighExp)
#
#loop <- data.frame(matrix(ncol=4, nrow=6))
#colnames(loop) <- c('Pearson', 'diffLoop', 'Pvalue', 'n')
#for (i in 1:6) {
#       bincorr <- rcorr(cohesin[abs(cohesin$InsLoopNum-cohesin$neighLoopNum)==(i-1),]$InsExp, cohesin[abs(cohesin$InsLoopNum-cohesin$neighLoopNum)==(i-1),]$neighExp)
#       loop[i,][1] <- as.numeric(bincorr$r[1,][2])
#       loop[i,][2] <- (i-1)
#       loop[i,][3] <- as.numeric(bincorr$P[1,][2])
#       loop[i,][4] <- nrow(cohesin[abs(cohesin$InsLoopNum-cohesin$neighLoopNum)==(i-1),])
#}
#
#pdf(paste0(opt$OUTDIR,'/graphs/Loop_correlation.pdf'), height=5, width=5) 
#loop %>% 
#       mutate(SampleName = opt$SampleName,
#              diffLoop = factor(diffLoop)) %>%
#       ggplot(aes(x=diffLoop, y=Pearson, fill=SampleName)) +
#              geom_bar(stat='identity', colour='black', alpha=0.6) +
#              scale_fill_manual(values=sampleColor) +
#              theme_classic() +
#              xlab('Cohesin+CTCF between insertions') +
#              ggtitle(opt$SampleName, 'Insertion in same or different loop\nSMAD21 + SMC3 + CTCF')+
#              theme(legend.position='none')
#dev.off()
#


#cat('Correlation of insertions in contact\n')
#
#contact <- read.table(paste0('/tmp/',opt$SampleName,'_InsContacts.tsv'), sep='\t', stringsAsFactors=F)
#contact$exp <- normdataiPCR[match(contact$V4, normdataiPCR$BC),]$expression
#contact$Neighexp <- normdataiPCR[match(contact$V5, normdataiPCR$BC),]$expression
#contact <- contact[!is.na(contact$exp) & !is.na(contact$Neighexp) & contact$V4 != contact$V5,]
#rcorr(contact$exp, contact$Neighexp)
#
#
#####
#10 MB windows for insertions TODO with more insertions
######
#library(Hmisc)
#library(gtools)
#MBwindow <- read.table('/tmp/CMV_1MBinsertion.tsv', sep='\t')
#
#MBwindowSplit <- split.data.frame(MBwindow, factor(MBwindow$V1))
#
#test <- lapply(MBwindowSplit, function(MBwin) {
#              if(nrow(MBwin)>4) {
#                     MBwin2<-data.frame(x=MBwin$V2[combinations(n=length(MBwin$V2), r=2, v=)[,1]], y=MBwin$V2[combinations(n=length(MBwin$V2), r=2, v=)[,2]])
#                     MBwin2$xExp <- normdataiPCR[match(MBwin2$x, normdataiPCR$BC),]$expression
#                     MBwin2$yExp <- normdataiPCR[match(MBwin2$y, normdataiPCR$BC),]$expression
#                     MBwin2corr <- rcorr(MBwin2$xExp , MBwin2$yExp )
#                     data.frame(MBpart = as.character(unique(MBwin[[1]])),
#                                                 Pearson = as.numeric(MBwin2corr$r[1,][2]),
#                                                 Pvalue = as.numeric(MBwin2corr$P[1,][2]))
#              }
#})
#MBlist<-as.data.frame(rbind_list(test))
#
#Methylation 
#####
cat('Methylation and activity\n')

meth <- subset(normdataiPCR, select=c(CGIovr, DistToCGI, CGImethylation, InsMethylation2kb, InsCpG2kb, PercentGC, expression))
meth$Quartile <- ntile(meth$InsMethylation2kb,4 )  
pdf(paste0('~/public_html/blog/2017Nov27/',opt$SampleName,'_2kbMethylation.pdf'))
boxplot( meth$expression~meth$Quartile, col='red', frame=F, xlab='Methylation quartiles', ylab='Reporter activity', main=paste0(opt$SampleName,'\nMean CpG-methylation +-2kb of insertion'), ylim=c(-2,2.1), outline=F)
text(labels=c(paste0('mean = ',signif(mean(meth[meth$Quartile==1,]$InsMethylation2kb),2)), paste0('mean = ',signif(mean(meth[meth$Quartile==2,]$InsMethylation2kb),2)), paste0('mean = ',signif(mean(meth[meth$Quartile==3,]$InsMethylation2kb),2)), paste0('mean = ',signif(mean(meth[meth$Quartile==4,]$InsMethylation2kb),2))),
       x=c(1, 2, 3, 4), y =c(1.5, 1.5,1.5 ,1.5)) 
dev.off()

cat('CGI Methylation and activity\n')

meth$CGImethQuartile <- ntile(meth$CGImethylation, 4)
pdf(paste0('~/public_html/blog/2017Nov27/',opt$SampleName,'_CGIdistance_andMethylation.pdf'), height=5, width=7)
meth %>%
       mutate(DistToCGI.bin =cut(abs(DistToCGI), breaks=c(0, 5e4, 1e5, 2e5, 3e5, 4e5, 5e5, 1e6), right=F, include.lowest=F, labels=c("< 50 kb", "< 100 kb", "< 200 kb","< 300 kb", "< 400 kb", "< 500 kb", "< 1 Mb")),
              CGImethQuartile = as.factor(CGImethQuartile)) %>% 
       filter(!is.na(DistToCGI.bin)) %>%
       ggplot(aes(x=DistToCGI.bin, y=expression)) +
              geom_boxplot(aes(fill=CGImethQuartile), outlier.shape=NA)+
              scale_fill_brewer(palette='Set2')+
              theme_classic()+
              ggtitle(opt$SampleName,'Distance to nearest CGI and methylation of the CGI')
dev.off()

######
#MCV of, and K562 specific, DHS
######
cat('MCV of, and K562 specific, DHS\n')

cat('mean MCV of nearest DHS for top and bottom 25%\n')
library(reshape2)
mcv <- read.table(paste0(opt$OUTDIR, '/', opt$SampleName, '_k562_mcv_k20.bed'), header=F)

mcv$DHS  <- rep(paste0('DHS_',formatC(1:20, width = 2, format = "d", flag = "0")), length(unique(mcv$V4)))
mcv <- dcast(data = mcv,formula = V4~DHS,fun.aggregate = sum,value.var = "V5")
  
mcv$expression <- normdataiPCR[match(mcv$V4, normdataiPCR$BC),]$expression

mcv2<- mcv %>%
       mutate(Expression = ifelse(expression>as.numeric(summary(expression)[5]), 1, ifelse(expression<(as.numeric(summary(expression)[2])), 0,NA))) %>%
       filter(!is.na(Expression))
rownames(mcv2) <- mcv2$V4
mcv2 <- subset(mcv2, select=-c(V4))

mcv2$meanMCV <- as.numeric(rowSums(mcv2)/20)

pdf(paste0(opt$OUTDIR,'/graphs/ExpDHS_meanMCV.pdf'))
boxplot(mcv2$meanMCV ~ mcv2$Expression, outline=F, names=c('nonExpressed', 'expressed'), ylab=c('mean MCV of 20 nearest DHS'), col=sampleColor,
main=paste0(opt$SampleName,'\nmean MCV of nearest 20 DHS\n top vs bottom 25% of reporter activity'), frame=F)
text(x=1.5, y=10, labels=paste0('pvalue = ', signif( t.test(mcv2[mcv2$Expression==1,]$meanMCV, mcv2[mcv2$Expression==0,]$meanMCV)$p.value,3)))
dev.off()



cat('Distance to nearest K562/Heme specific DHS\n')
K562specDHS <- subset(normdataiPCR, select=c(BC, DistToNearestDHS, DNasePksDens, nDHS100kb, MCVK562NonSpecDHS, DistToNearestK562SpecDHS, nK562SpecDHS, DistToNearestK562NonSpecDHS, nK562NonSpecDHS, expression))
#K562specDHS <-read.table('/tmp/AllInsertions.coords.bed')
#colnames(K562specDHS) <- c("V1", "V2", "V3", "V4", "V5", "V6", "DistToNearestK562SpecDHS", "nK562SpecDHS", "DistToNearestK562NonSpecDHS", "nK562NonSpecDHS", 'MCVK562NonSpecDHS')
#K562specDHS$expression <- normdataiPCR[match(K562specDHS$V4, normdataiPCR$BC),]$expression
#K562specDHS$DistToNearestDHS <- normdataiPCR[match(K562specDHS$V4, normdataiPCR$BC),]$DistToNearestDHS
#K562specDHS$nDHS100kb <- normdataiPCR[match(K562specDHS$V4, normdataiPCR$BC),]$nDHS100kb
#K562specDHS$DNasePksDens <- normdataiPCR[match(K562specDHS$V4, normdataiPCR$BC),]$DNasePksDens
cat('Dist To K562 specific DHS vs nonspecific DHS\n')
pdf(paste0(opt$OUTDIR,'/graphs/DistToK562DHSvsDHS.pdf'))
plot(log10(abs(K562specDHS$DistToNearestK562NonSpecDHS)+1), log10(abs(K562specDHS$DistToNearestK562SpecDHS)+1), frame=F, col=sampleColor, 
       xlab='Distance to nearest non specific DHS', ylab='Distance to nearest K562 specific DHS',
       main=paste0(opt$SampleName,'\nDistance to nearest DHS vs nearest K562 specific DHS'))
dev.off()

K562specDHS1E5 <- K562specDHS[abs(K562specDHS$DistToNearestK562SpecDHS)<1e5,]
K562specDHS1E5$Decile <- ntile(abs(K562specDHS1E5$DistToNearestK562SpecDHS), n=4)
K562specDHS1E5$SpecNearer <- abs(K562specDHS1E5$DistToNearestK562SpecDHS) < abs(K562specDHS1E5$DistToNearestK562NonSpecDHS)
K562specDHS1E5[K562specDHS1E5$SpecNearer==TRUE,]$SpecNearer <- 'Specific'
K562specDHS1E5[K562specDHS1E5$SpecNearer==FALSE,]$SpecNearer <- 'nonSpecific'


cat('Dist To K562 specific DHS vs nonspecific DHS expression boxplot\n')

pdf(paste0(opt$OUTDIR,'/graphs/DistToK562DHSvsDHS_expressionboxplot.pdf'), width=9)
K562specDHS1E5 %>%
       mutate(SampleName = opt$SampleName) %>%
       ggplot(aes(x=as.factor(Decile), y=expression, fill=SpecNearer)) +
              geom_boxplot(outlier.shape=NA, notch=T)+
              theme_classic() +
              ylab('Reporter activity') +
              xlab('Distance quartile up to 100kb') +
              scale_fill_brewer(palette='Pastel2') +
              ggtitle(opt$SampleName,'Expression and distance to DHS and if the nearest DHS is K562 specific or not')+
              annotate('text', x=c(1, 1 ,1), y=c(1.7, 1.3, 1.5), label =c(paste0('pvalue = ',signif(t.test(K562specDHS1E5[K562specDHS1E5$Decile==1 & K562specDHS1E5$SpecNearer=='Specific',]$expression, K562specDHS1E5[K562specDHS1E5$Decile==1 & K562specDHS1E5$SpecNearer=='nonSpecific',]$expression)$p.value,2)),
                                                                      paste0('nSpecific = ',nrow(K562specDHS1E5[K562specDHS1E5$Decile==1 & K562specDHS1E5$SpecNearer=='Specific',])),
                                                                      paste0('nNonSpecific = ', nrow(K562specDHS1E5[K562specDHS1E5$Decile==1 & K562specDHS1E5$SpecNearer=='nonSpecific',])))) +
              annotate('text', x=c(2, 2 ,2), y=c(1.7, 1.3, 1.5), label =c(paste0('pvalue = ',signif(t.test(K562specDHS1E5[K562specDHS1E5$Decile==2 & K562specDHS1E5$SpecNearer=='Specific',]$expression, K562specDHS1E5[K562specDHS1E5$Decile==2 & K562specDHS1E5$SpecNearer=='nonSpecific',]$expression)$p.value,2)),
                                                                      paste0('nSpecific = ',nrow(K562specDHS1E5[K562specDHS1E5$Decile==2 & K562specDHS1E5$SpecNearer=='Specific',])),
                                                                      paste0('nNonSpecific = ', nrow(K562specDHS1E5[K562specDHS1E5$Decile==2 & K562specDHS1E5$SpecNearer=='nonSpecific',])))) +
              annotate('text', x=c(3, 3 ,3), y=c(1.7, 1.3, 1.5), label =c(paste0('pvalue = ',signif(t.test(K562specDHS1E5[K562specDHS1E5$Decile==3 & K562specDHS1E5$SpecNearer=='Specific',]$expression, K562specDHS1E5[K562specDHS1E5$Decile==3 & K562specDHS1E5$SpecNearer=='nonSpecific',]$expression)$p.value,2)),
                                                                      paste0('nSpecific = ',nrow(K562specDHS1E5[K562specDHS1E5$Decile==3 & K562specDHS1E5$SpecNearer=='Specific',])),
                                                                      paste0('nNonSpecific = ', nrow(K562specDHS1E5[K562specDHS1E5$Decile==3 & K562specDHS1E5$SpecNearer=='nonSpecific',])))) +
              annotate('text', x=c(4, 4 ,4), y=c(1.7, 1.3, 1.5), label =c(paste0('pvalue = ',signif(t.test(K562specDHS1E5[K562specDHS1E5$Decile==4 & K562specDHS1E5$SpecNearer=='Specific',]$expression, K562specDHS1E5[K562specDHS1E5$Decile==4 & K562specDHS1E5$SpecNearer=='nonSpecific',]$expression)$p.value,2)),
                                                                      paste0('nSpecific = ',nrow(K562specDHS1E5[K562specDHS1E5$Decile==4 & K562specDHS1E5$SpecNearer=='Specific',])),
                                                                      paste0('nNonSpecific = ', nrow(K562specDHS1E5[K562specDHS1E5$Decile==4 & K562specDHS1E5$SpecNearer=='nonSpecific',])))) +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))
dev.off()

cat('Accessibility of nearest K562 specific DHS vs nonspecific DHS\n')
pdf(paste0(opt$OUTDIR,'/graphs/DensityofNearestK562DHSvsDHS_expressionboxplot.pdf'), height=5, width=5)
K562specDHS1E5 %>%
         mutate(SampleName = opt$SampleName) %>%
       ggplot(aes(x=as.factor(SpecNearer), y=log10(DNasePksDens+1), fill=SampleName)) +
              geom_boxplot(outlier.shape=NA, notch=T)+
              theme_classic() +
              ylab('DHS peak accessibility') +
              xlab('Specificity of DHS') +
              scale_fill_manual(values=sampleColor) +
              ylim(0, 2)+
              ggtitle(opt$SampleName,'Accessibility of nearest DHS peak')+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12), legend.position='none')+
              annotate('text', x=1.5, y=2.0, label = paste0('pvalue = ', signif(t.test(log10(K562specDHS1E5[K562specDHS1E5$SpecNearer=='Specific',]$DNasePksDens+1), log10(K562specDHS1E5[K562specDHS1E5$SpecNearer=='nonSpecific',]$DNasePksDens+1))$p.value,2)))
dev.off()


cat('Number of K562 specific DHS vs nonspecific DHS +- 100kb and expression\n')
pdf(paste0(opt$OUTDIR,'/graphs/n100kb562DHSvsDHS_expressionboxplot.pdf'), width=9)
K562specDHS1E5$moreDHSisSpec <-  K562specDHS1E5$nK562SpecDHS > K562specDHS1E5$nK562NonSpecDHS
K562specDHS1E5$nDHSqaurtile <- ntile(K562specDHS1E5$nDHS100kb, n=4)
K562specDHS1E5  %>%
       mutate(SampleName = opt$SampleName) %>%
       ggplot(aes(x=as.factor(nDHSqaurtile), y=expression, fill=moreDHSisSpec)) +
              geom_boxplot(outlier.shape=NA, notch=T)+
              theme_classic() +
              ylab('Reporter activity') +
              xlab('Quartile of number of DHS +- 100kb') +
              scale_fill_brewer(palette='Pastel2') +
              ggtitle(opt$SampleName,'Expression and number of DHS +- 100kb and if theres more K562 specific or nonspecific DHS')+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) 
dev.off()

cat('Number of K562 specific DHS vs nonspecific DHS +- 100kb\n')
pdf(paste0(opt$OUTDIR,'/graphs/n100kbK562DHSvsDHS.pdf'))
smoothScatter(y=K562specDHS$nK562SpecDHS, x=K562specDHS$nK562NonSpecDHS, colramp = colorRampPalette(c("white", sampleColor)),
              ylab='number of K562 specific DHS', xlab='number of K562 nonspecific DHS',
              main=paste0(opt$SampleName,'\nNumber of surrounding DHS +- 100kb and specificity of DHS'))
dev.off()






#cat('Divide nearest DHS into K562 specific, medium, and nonspecific \n')
#K562specDHS1E5$specificity <- NA
#K562specDHS1E5[K562specDHS1E5$MCVK562NonSpecDHS==0,]$specificity <- 'Specific'
#K562specDHS1E5[K562specDHS1E5$MCVK562NonSpecDHS>0  & K562specDHS1E5$MCVK562NonSpecDHS<15 ,]$specificity <- 'mediumSpecific'
#K562specDHS1E5[K562specDHS1E5$MCVK562NonSpecDHS>=15 ,]$specificity <- 'nonSpecific'
#
#K562specDHS1E5 %>%
#       mutate(SampleName = opt$SampleName,
#       specificity = factor(specificity, levels=c('nonSpecific', 'mediumSpecific', 'Specific'))) %>%
#       ggplot(aes(x=as.factor(Decile), y=expression, fill=specificity)) +
#              geom_boxplot(outlier.shape=NA, notch=T)+
#              theme_classic() +
#              ylab('Reporter activity') +
#              xlab('Distance quartile up to 100kb') +
#              scale_fill_brewer(palette='Pastel2')
#
#




######
#ChIP-enrichment for top and bottom 25% insertions
######
library(reshape2)
library(pheatmap)
library(gplots)
library(dplyr)

if (which(list.files(opt$OUTDIR) == paste0(opt$SampleName, "_ChIPhighExpressed.tsv")) >0 ) {
       cat('ChIP-enrichment for top and bottom 25% insertions\n')
       library(dplyr)
       library(ggrepel)
       library(ggforce)
       ChIPhigh <- read.table(paste0(opt$OUTDIR,'/', opt$SampleName, "_ChIPhighExpressed.tsv"))
       ChIPhigh$Enrichment <- ChIPhigh$V4/ChIPhigh$V5
       
       ChIPlow <- read.table(paste0(opt$OUTDIR,'/', opt$SampleName, "_ChIPlowExpressed.tsv"))
       ChIPlow$Enrichment <- ChIPlow$V4/ChIPlow$V5

       
       ChIP <- rbind(ChIPhigh, ChIPlow)
       ChIP$V1 <- gsub('K562-', '', ChIP$V1)
       ChIP$V1 <- gsub('-.*', '', ChIP$V1)
       
       
       
       
       
       for (res in unique(ChIP$V2)) {
              cat(res ,'\n')
              ChIPres <- ChIP[ChIP$V2==res,]
              HighChIP <- ChIPres[ChIPres$V6=='High' & ChIPres$V4>=100,]
              FactorOrder <- HighChIP[rev(order(HighChIP$Enrichment)),]$V1
              ChIPres <- ChIPres[ChIPres$V1 %in% HighChIP$V1,]
              
              ggOrdChIPres<- ChIPres %>% 
                     #filter(V4>=100) %>%
                     mutate(V1 = factor(V1, levels=FactorOrder),
                            Expression = V6) %>%
                     ggplot(aes(x=V1, y=Enrichment, fill=Expression, group=Expression)) +
                           geom_point(pch=21, size=3)+
                           theme_classic()+
                           ggtitle(paste0('Enrichment of ChIP-seq marks +-',res/2,'bp around High and low expressed insertions\ncompared to random')) +
                           scale_fill_brewer(palette='Set1')+
                           #geom_text_repel(data=ChIP[ChIP$Enrichment>=2.5 & ChIP$V4>=100,], aes(x=V1, y=Enrichment, label = V1))+
                            #geom_text_repel(data=ChIP[ChIP$Enrichment>=3 & ChIP$V4>=100,], aes(x=V1, y=Enrichment, label = V1))+
                            geom_hline(yintercept=1, col='red', linetype='dashed')+
                            xlab('ChIP-seq marks')#+
                            #facet_zoom(y = Enrichment>=2.5)
              ggsave(ggOrdChIPres, file=paste0(opt$OUTDIR,'/graphs/ChIPenrichment_',res/2,'ordered.pdf'), height=6, width=12)
       
              
              namedRes <-  head(HighChIP[rev(order(HighChIP$Enrichment)),], nrow(HighChIP[rev(order(HighChIP$Enrichment)),])*0.20)$V1
              ChIPnamedRes <- ChIPres[ChIPres$V1 %in%namedRes,]
              
              ggNameChIPres<-  ChIPnamedRes[ChIPnamedRes$V6=='High',] %>% 
                     #filter(V4>=100) %>%
                     mutate(V1 = factor(V1, levels=rev(namedRes)),
                            Expression = V6) %>% 
                     ggplot(aes(x=V1, y=Enrichment, fill=Expression, group=V1)) +
                            geom_segment(aes(x = V1, y = 0, xend = V1, yend = Enrichment), color = "grey50") +
                            geom_point(pch=21, size=4)+
                            theme_classic()+
                            ylim(0, 4) +
                            ggtitle(paste0('Enrichment of ChIP marks +-',paste0(res/2000, 'kb'),'\naround High expressed insertions\ncompared to random') )+
                            scale_fill_brewer(palette='Set1')+
                            theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12), legend.position='none')+
                            #geom_text_repel(data=ChIP[ChIP$Enrichment>=2.5 & ChIP$V4>=100,], aes(x=V1, y=Enrichment, label = V1))+
                            #geom_label_repel(data=ChIP[ChIP$V1=='DNase' & ChIP$Enrichment>=2.2 ,], aes(x=V1, y=Enrichment, label = V1), fill='white')+
                            geom_hline(yintercept=1, col='red', linetype='dashed')+
                            xlab('ChIP-seq marks')+
                            coord_flip()
         
              ggsave(ggNameChIPres, file=paste0(opt$OUTDIR,'/graphs/ChIPenrichment_',res/2,'namedHigh.pdf'), height=8, width=5)
              
              
              namedRes <-  tail(HighChIP[rev(order(HighChIP$Enrichment)),], nrow(HighChIP[order(HighChIP$Enrichment),])*0.20)$V1
              ChIPnamedRes <- ChIPres[ChIPres$V1 %in%namedRes,]
              
                ggNameChIPres<-   ChIPnamedRes[ChIPnamedRes$V6=='High',] %>% 
                     #filter(V4>=100) %>%
                     mutate(V1 = factor(V1, levels=rev(namedRes)),
                            Expression = V6) %>%
                     ggplot(aes(x=V1, y=Enrichment, fill=Expression, group=V1)) +
                            geom_segment(aes(x = V1, y = 0, xend = V1, yend = Enrichment), color = "grey50") +
                            geom_point(pch=21, size=4)+
                            theme_classic()+
                            ggtitle(paste0('Enrichment of ChIP marks +-',paste0(res/2000, 'kb'),'\naround High expressed insertions\ncompared to random') )+
                            theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12), legend.position='none')+
                            scale_fill_brewer(palette='Set1')+
                            ylim(0, 4) +
                            #geom_text_repel(data=ChIP[ChIP$Enrichment>=2.5 & ChIP$V4>=100,], aes(x=V1, y=Enrichment, label = V1))+
                            #geom_label_repel(data=ChIP[ChIP$V1=='DNase' & ChIP$Enrichment>=2.2 ,], aes(x=V1, y=Enrichment, label = V1), fill='white')+
                            geom_hline(yintercept=1, col='red', linetype='dashed')+
                            xlab('ChIP-seq marks')+
                            coord_flip()
         
              ggsave(ggNameChIPres, file=paste0(opt$OUTDIR,'/graphs/ChIPenrichment_',res/2,'namedLow.pdf'), height=8, width=5)
       
       }
       
       
       ChIPhighDF <- ChIPhigh
       ChIPhighDF$V1 <- gsub('K562-', '', ChIPhighDF$V1)
       ChIPhighDF$V1 <- gsub('-.*', '', ChIPhighDF$V1)
       ChIPhighDF <- ChIPhighDF[ChIPhighDF$V4 >=100,]
       ChIPhighDF$V2 <- ChIPhighDF$V2/2
       
       ChIPhighDF <- dcast(data = ChIPhighDF,formula = V1~V2,fun.aggregate = sum,value.var = "Enrichment")
       rownames(ChIPhighDF) <- ChIPhighDF$V1
       ChIPhighDF <- subset(ChIPhighDF, select=-c(V1))
       colnames(ChIPhighDF) <- c('1kb', '2.5kb', '10kb', '25kb', '50kb', '100kb')
       ChIPhighDF <- ChIPhighDF[rev(order(ChIPhighDF[,"10kb"])),]
       ChIPhighDF[,"1kb"][ChIPhighDF[,"1kb"]==0] <- NA
       #ChIPhighDF <- ChIPhighDF[!(apply(ChIPhighDF, 1, function(y) any(y == 0))),]
       #pdf(paste0(opt$OUTDIR,'/graphs/ChIPenrichment_Comparison.pdf'), height=7, width=6)
       #heatmap.2(as.matrix(head(ChIPhighDF, 50)),
       #       Rowv=F, Colv=F, dendrogram="none", 
       #       col=colorRampPalette(brewer.pal(n = 9, name = "Reds"))(100), 
       #      sepwidth=c(0.01,0.01),
       #       sepcolor="lightgray",
       #       colsep=1:ncol(head(ChIPhighDF, 50)),
       #       rowsep=1:nrow(head(ChIPhighDF, 50)),
       #       linecol=F,
       #       trace='none',
       #       density.info="none",
       #       keysize=1,
       #       key.title=NA,
       #       key.xlab="Enrichment",
       #       key.ylab=NA,
       #       lmat = rbind(4:3,2:1),
       #       lhei = c(0.8,4),
       #       lwid = c(1.5,4),
       #       main=paste0(opt$SampleName,'\nEnrichment of ChIP-seq \npeaks around the insertions'))
       #dev.off()
       
       ChIPhighDF$ChIP <- rownames(ChIPhighDF)
       ChIPhighDF <- ChIPhighDF[rev(order(ChIPhighDF[,"10kb"])),]
       ChIPorder <- head(ChIPhighDF[rev(order(ChIPhighDF[,"10kb"])),]$ChIP, 50)
       
      ggChIPHeat <- head(ChIPhighDF, 50) %>%
              melt() %>%
              mutate(Bins = factor(variable, levels=c('1kb', '2.5kb','10kb', '25kb', '50kb', '100kb')),
                     ChIP = factor(ChIP, levels=rev(ChIPorder)),
                     SampleName=opt$SampleName,
                     Enrichment = value) %>%
              ggplot(aes(x=Bins,y=ChIP)) +
              geom_tile(aes(fill = Enrichment), colour = "black") + 
              scale_fill_gradient(low = "white",high = sampleColor)+
              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12), legend.position='none')+
              theme_minimal()+
              ggtitle(opt$SampleName,'Enrichment of ChIP-seq \npeaks around the insertions')
              
       ggsave(ggChIPHeat,  file=paste0(opt$OUTDIR,'/graphs/ChIPenrichment_Comparison.pdf'), height=7, width=6)
       
       ggChIPFreqMax <- as.data.frame(table(colnames(ChIPhighDF)[apply(ChIPhighDF,1,which.max)])) %>%
              mutate(Bins = factor(Var1, levels=c('1kb', '2.5kb','10kb', '25kb', '50kb', '100kb')),
                     SampleName=opt$SampleName) %>%
              ggplot(aes(x=Bins, y=Freq, fill=SampleName)) +
                     geom_bar(stat='identity', colour='black', alpha=0.6) +
                     scale_fill_manual(values=sampleColor)+
                     theme_classic()+
                     theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12), legend.position='none')+
                     ggtitle(opt$SampleName, 'Which bin around the insertions\nhas highest enrichment per ChIP-mark')
       ggsave(ggChIPFreqMax, file=paste0(opt$OUTDIR,'/graphs/ChIPenrichment_Comparison_FreqMax.pdf'), height=5, width=4)
       
       ggChIPFreqIns <- ChIPhigh %>%
              mutate( Bins = paste0(V2/2000, 'kb'),
                     Bins = factor(Bins, levels=c('1kb', '2.5kb','10kb', '25kb', '50kb', '100kb')),
                     ChIP = gsub('K562-', '', V1),
                     ChIP = gsub('-.*', '', ChIP),
                     SampleName=opt$SampleName) %>%
              ggplot(aes(x=Bins, y=V4, fill=SampleName)) +
                     geom_boxplot(outlier.shape=NA, alpha=0.6) +
                     scale_fill_manual(values=sampleColor)+
                     theme_classic()+
                     ylab('Number of insertions with  ChIP-mark') +
                     theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12), legend.position='none')+
                     ggtitle(opt$SampleName, paste0('How many insertions has a specific\nChIP-mark around them\nn = ',unique(ChIPhigh$V3)))
       ggsave(ggChIPFreqIns, file=paste0(opt$OUTDIR,'/graphs/ChIPenrichment_Comparison_FreqChIP.pdf'), height=5, width=4)
              

              
       
       #pheatmap(head(ChIPhighDF, 50), cluster_cols=F, cluster_rows=F, main='Enr)
}




#boxplot(K562specDHS$expression~K562specDHS$K562specDist
#####
#Contribution to variability and featuers between neighours 
#####
#neighFeat <- neighAB
#neighFeat$InsLoc <- normdataiPCR[match(neighFeat$BC, normdataiPCR$BC),]$genomicDivisionNoExp
#neighFeat$NeighLoc <- normdataiPCR[match(neighFeat$neighBC, normdataiPCR$BC),]$genomicDivisionNoExp
#
#neighFeat$geneRegion <-''
#neighFeat[neighFeat$InsLoc!=neighFeat$NeighLoc,]$geneRegion <- 'RichDesert'
#neighFeat[neighFeat$InsLoc=='geneDesert' & neighFeat$NeighLoc=='geneDesert',]$geneRegion <- 'DesertDesert'
#neighFeat[neighFeat$InsLoc=='geneRich' & neighFeat$NeighLoc=='geneRich',]$geneRegion <- 'RichRich'
#
#
#neighFeat$Variability <- apply(subset(neighFeat, select=c(Expression,NeighbourExpression)) , 1, function(x) var(x))
#neighFeat$VarianceDecile <- ntile(neighFeat$Variability, 10) 
#neighFeatGLM <- neighFeat %>%
#       mutate(Compartment = as.factor(Compartment),
#              geneRegion = as.factor(geneRegion))
#
#neighFeatGLM <- neighFeatGLM[neighFeatGLM$VarianceDecile==1 | neighFeatGLM$VarianceDecile==10,]
#neighFeatGLM[neighFeatGLM$VarianceDecile==10,]$VarianceDecile <- 0
#neighFeatGLM <- subset(neighFeatGLM, select=c(Distance, nCTCF, nDHS, nCGI, Compartment, geneRegion, VarianceDecile))
#
#neighFeatGLM$VarianceDecile <- as.factor(neighFeatGLM$VarianceDecile)
#
#neighROC <- data.frame(matrix(ncol=length(colnames(neighFeatGLM))-1, nrow=1))
#for (i in 1:(length(colnames(neighFeatGLM))-1)) {
#
#       fitNeigh <-subset(neighFeatGLM, select=c(colnames(neighFeatGLM)[i], 'VarianceDecile'))
#       fit <- glm(formula=VarianceDecile ~. ,data=fitNeigh, family='binomial')
#       preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
#       perf <- performance(preds, "auc")
#       neighROC[i] <-round(as.numeric(slot(perf, "y.values")), 3)
#}
#
#fit <-glm(formula=VarianceDecile ~. ,data=neighFeatGLM, family='binomial')
# preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
#       perf <- performance(preds, "auc")
#       round(as.numeric(slot(perf, "y.values")), 3)
#
####
#GAM
#####
library(mgcv)
library(ROCR)
library(gridExtra)
cat('GLM for features\n')

featureGAM <- subset(normdataiPCR,select=c(BC, DNA, RNA, DistToTSS, DistToNearestDHS, DistToNearestCTCF, DistToCGI, PercentGC, AB_compartment_value,
 chromHMM, nDHS100kb, nDHS10kb, nDHS300kb, nDHS500kb, nDHS1Mb, nDHS100kbDensity, nDHS300kbsum, nDHS500kbDensity, nDHS1MbDensity, nTSS100kb, 
 repliSeq, nExpressed100kb, nCGI100kb, nCTCFpeaks100kb, DNasePksDens, NearestDHStype, CTCFPksDens, genomicDivision, DistToGeneDesert, DistTok562GeneDesert, meanDist5DHS, distToDHScluster, 
 H2AFZ, H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me2, H3K4me3, H3K79me2, H3K9ac, H3K9me1, H3K9me3, H4K20me1, HDAC1, HDAC2, HDAC6, totalTags50kb, DistToNearestDHShot2, expression))
 
featureGAM2 <- cbind(featureGAM,distdhsTypeROC[match(featureGAM$BC, distdhsTypeROC$V1),])


featureGAM.clean <- featureGAM2 %>%
       filter(!is.na(repliSeq),
              !is.na(chromHMM),
              !is.na(NearestDHStype),
              !is.na(genomicDivision),
              !is.na(AB_compartment_value),
              DNA>=1) %>%
       mutate(Expression = ifelse(log10(RNA/DNA) >as.numeric(summary(expression)[5]), 1, ifelse(log10(RNA/DNA) <(as.numeric(summary(expression)[2])), 0,NA)), 
              DistToTSS = abs(DistToTSS),
              chromHMM = factor(chromHMM),
              NearestDHStype = factor(NearestDHStype),
              DistToNearestDHShot2 = abs(DistToNearestDHShot2),  
              repliSeq = factor(repliSeq),
              genomicDivision = factor(genomicDivision),
              DistToNearestDHS = abs(DistToNearestDHS),
              DistToCGI = abs(DistToCGI),
              DistToNearestCTCF = abs(DistToNearestCTCF),
              distToDHScluster =abs(distToDHScluster),
              H2AFZ= abs(H2AFZ),
              H3K27ac = abs(H3K27ac), 
              H3K27me3 = abs(H3K27me3), 
              H3K36me3 = abs(H3K36me3), 
              H3K4me1 =abs(H3K4me1), 
              H3K4me2 = abs(H3K4me2), 
              H3K4me3 = abs(H3K4me3),
              H3K79me2 = abs(H3K79me2), 
              H3K9ac =abs(H3K9ac), 
              H3K9me1 =abs(H3K9me1), 
              H3K9me3 =abs(H3K9me3), 
              H4K20me1 =abs(H4K20me1), 
              HDAC1 =abs(HDAC1), 
              HDAC2 =abs(HDAC2), 
              HDAC6 =abs(HDAC6)) %>%
       filter(!is.na(Expression))


GLMfeat <- glm(Expression ~ #log10(DistToCGI+1) + 
              log10(DistToNearestDHS+1)+ 
              nDHS100kb,
            #log10(nDHS100kbDensity+1) + 
            #AB_compartment_value +
            #log10(meanDist5DHS+1) +
            #log10(H3K9ac+1) + 
            #log10(H3K27ac+1) +
            #log10(H3K79me2+1) + 
            #log10(H3K4me2+1) +
            #log10(H2AFZ+1) + log10(H3K36me3+1) + log10(H3K36me3+1) + log10(H3K4me1+1) + log10(H3K9me1+1) + log10(H3K9me3+1) +
            #log10(H4K20me1+1) + log10(HDAC1+1) + log10(HDAC2+1) + log10(HDAC6+1) + repliSeq, 
              data = featureGAM.clean, family=binomial)
preds <- prediction(GLMfeat$fitted.values, GLMfeat$y)
perf <- performance(preds, "auc")
rocC <- round(as.numeric(slot(perf, "y.values")), 3)
rocC

GLM_plot <- performance(preds, "tpr", "fpr")


random_order_Expressed <- sample(which(featureGAM.clean$Expression == 1))
random_order_nonExpressed <- sample(which(featureGAM.clean$Expression == 0))


inTrain <- c(random_order_Expressed[1:2500], random_order_nonExpressed[1:2500])
inValidate <- c(random_order_Expressed[2501:3036], random_order_nonExpressed[2501:3038])


training <- featureGAM.clean[inTrain,]
validate <- featureGAM.clean[inValidate,]


GLMtrain <- glm(Expression ~ #log10(DistToCGI+1) + 
              log10(DistToNearestDHS+1)+ 
              nDHS100kb,
            #log10(nDHS100kbDensity+1) + 
            #AB_compartment_value +
            #log10(meanDist5DHS+1) +
            #log10(H3K9ac+1) + 
            #log10(H3K27ac+1) +
            #log10(H3K79me2+1) + 
            #log10(H3K4me2+1) +
            #log10(H2AFZ+1) + log10(H3K36me3+1) + log10(H3K36me3+1) + log10(H3K4me1+1) + log10(H3K9me1+1) + log10(H3K9me3+1) +
            #log10(H4K20me1+1) + log10(HDAC1+1) + log10(HDAC2+1) + log10(HDAC6+1) + repliSeq, 
              data = training, family=binomial)

preds <- prediction(GLMtrain$fitted.values, GLMtrain$y)
perf <- performance(preds, "auc")
rocCTrain <- round(as.numeric(slot(perf, "y.values")), 3)
rocCTrain
GLM_trainplot <- performance(preds, "tpr", "fpr")

glmPred <- predict(GLMtrain, validate, type='response')

preds <- prediction(as.numeric(glmPred), validate$Expression)
GLM_Valplot <- performance(preds, "tpr", "fpr")
perf <- performance(preds, "auc")
rocCVal <- round(as.numeric(slot(perf, "y.values")), 3)
rocCVal


validate$predExp <- predict(GLMtrain, validate)

cat('GLM for features plotting\n')
GLM_plot.ROC <- data.frame(y=GLM_plot@y.values[[1]] ,  x=GLM_plot@x.values[[1]], SampleName=unique(normdataiPCR$SampleName))
cat('Plot AUROC for GLM features\n')
AUC.features <- ggplot(GLM_plot.ROC ,aes(x=x, y=y, colour=SampleName)) +
       geom_line(size = 2)+
       theme_classic()+
       ggtitle(opt$SampleName,"  \nAUC for bottom and top quartile\nDistToDHS & nDHS100kb")+
       scale_colour_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       xlab('1-Specificity')+
       ylab('Sensitivity')+
       annotate('text', label = paste0('AUC = ',round(as.numeric(slot(perf, "y.values")), 3)), x = 0.5, y = 0.5,  color = 'black', size=10)+
       theme(legend.position="none")

pdf(paste0(opt$OUTDIR,'/graphs/GLM_features.pdf'), height=5, width=5) 
AUC.features
dev.off()
cat('K562 GLM\n')
K562gencode <- read.table('~/public_html/blog/2017Nov06/K562Gencode/AllInsertions.annotated.bed', sep='\t', header=T)
K562gencode$Expression <- K562gencode$value>=1
K562gencode$Expression[K562gencode$Expression=='TRUE'] <- 1
K562gencode$Expression[K562gencode$Expression=='FALSE'] <- 0



K562 <- K562gencode %>%
       filter(!is.na(repliSeq),
              !is.na(chromHMM),
              !is.na(NearestDHStype),
              !is.na(genomicDivision),
              !is.na(AB_compartment_value)) %>%
       mutate(DistToTSS = abs(DistToTSS),
              chromHMM = factor(chromHMM),
              NearestDHStype = factor(NearestDHStype),
              repliSeq = factor(repliSeq),
              genomicDivision = factor(genomicDivision),
              DistToNearestDHS = abs(DistToNearestDHS),
              DistToCGI = abs(DistToCGI),
              DistToNearestCTCF = abs(DistToNearestCTCF),
              distToDHScluster =abs(distToDHScluster),
              H2AFZ= abs(H2AFZ),
              H3K27ac = abs(H3K27ac), 
              H3K27me3 = abs(H3K27me3), 
              H3K36me3 = abs(H3K36me3), 
              H3K4me1 =abs(H3K4me1), 
              H3K4me2 = abs(H3K4me2), 
              H3K4me3 = abs(H3K4me3),
              H3K79me2 = abs(H3K79me2), 
              H3K9ac =abs(H3K9ac), 
              H3K9me1 =abs(H3K9me1), 
              H3K9me3 =abs(H3K9me3), 
              H4K20me1 =abs(H4K20me1), 
              HDAC1 =abs(HDAC1), 
              HDAC2 =abs(HDAC2), 
              HDAC6 =abs(HDAC6)) 

glmK562Pred <- predict(GLMtrain, K562, type='response')
preds <- prediction(as.numeric(glmK562Pred), K562$Expression)
GLM_K562plot <- performance(preds, "tpr", "fpr")
perf <- performance(preds, "auc")
rocCK562 <- round(as.numeric(slot(perf, "y.values")), 3)
rocCK562

K562$predExp <- predict(GLMtrain, K562)
GLM_All <- rbind( data.frame(y=GLM_trainplot@y.values[[1]], x=GLM_trainplot@x.values[[1]], SampleName='Training set'),
       data.frame(y=GLM_Valplot@y.values[[1]], x=GLM_Valplot@x.values[[1]], SampleName='Validation set'),
       data.frame(y=GLM_K562plot@y.values[[1]], x=GLM_K562plot@x.values[[1]], SampleName='K562 Gencode predicted'))
cat('Training, validating, and K562 AUROC GLM\n')
pdf(paste0(opt$OUTDIR,'/graphs/GLM_distDHSnDHSpredictExp.pdf'), height=5, width=7) 
ggplot(GLM_All ,aes(x=x, y=y, colour=SampleName)) +
       geom_line(size = 1)+
       theme_classic()+
       #ggtitle(opt$SampleName,"  \nAUC for bottom and top quartile\nDistToDHS & nDHS100kb")+
       scale_colour_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       xlab('1-Specificity')+
       ylab('Sensitivity')+
       annotate('text', label = paste0('AUC Train= ',round(rocCTrain, 3)), x = 0.5, y = 0.7,  color = '#E41A1C', size=6)+
       annotate('text', label = paste0('AUC Val= ',round(rocCVal, 3)), x = 0.5, y = 0.6,  color = '#377EB8', size=6)+
       annotate('text', label = paste0('AUC K562= ',round(rocCK562, 3)), x = 0.5, y = 0.5,  color = '#4DAF4A', size=6)
dev.off()
       #theme(legend.position="none")
#K562rand <- rbind(K562[K562$Expression=='TRUE',], K562[sample(nrow(K562[K562$Expression=='FALSE',]), nrow(K562[K562$Expression=='FALSE',])),])
#K562 <- K562[K562$value>0,]
GLMK562 <- glm(Expression ~ #log10(DistToCGI+1) + 
              log10(DistToNearestDHS+1)+ 
              nDHS100kb,
             #log10(nDHS100kbDensity+1) + 
             #AB_compartment_value +
             #log10(meanDist5DHS+1) +
             #log10(H3K9ac+1) + 
             #log10(H3K27ac+1) +
             #log10(H3K79me2+1) + 
             #log10(H3K4me2+1) + 
             #log10(H2AFZ+1) + log10(H3K36me3+1) + log10(H3K36me3+1) + log10(H3K4me1+1) + log10(H3K9me1+1) + log10(H3K9me3+1) +
             #log10(H4K20me1+1) + log10(HDAC1+1) + log10(HDAC2+1) + log10(HDAC6+1) + repliSeq, 
              data = K562, family=binomial)
preds <- prediction(GLMK562$fitted.values, GLMK562$y)
perf <- performance(preds, "auc")
rocC <- round(as.numeric(slot(perf, "y.values")), 3)
rocC

K562_plot <- performance(preds, "tpr", "fpr")
cat('Plot AUROC for GLM features for K562\n')

K562_plot.ROC <- data.frame(y=K562_plot@y.values[[1]] ,  x=K562_plot@x.values[[1]], SampleName='K562')

K562.features <- ggplot(K562_plot.ROC ,aes(x=x, y=y, colour=SampleName)) +
       geom_line(size = 2)+
       theme_classic()+
       ggtitle("K562  \nAUC for expressed (>=1 TPM) Gencode genes\nDistToDHS & nDHS100kb")+
       scale_colour_manual(values='#2ca25f') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       xlab('1-Specificity')+
       ylab('Sensitivity')+
       annotate('text', label = paste0('AUC = ',round(as.numeric(slot(perf, "y.values")), 3)), x = 0.5, y = 0.5,  color = 'black', size=10)+
       theme(legend.position="none")

pdf(paste0(opt$OUTDIR,'/graphs/K562_features.pdf'), height=5, width=5) 
K562.features
dev.off()



####
#Predicted expression based on GLM CMV model
####
featureLM<- featureGAM2 %>%
       filter(!is.na(repliSeq),
              !is.na(chromHMM),
              !is.na(NearestDHStype),
              !is.na(genomicDivision),
              !is.na(AB_compartment_value),
              DNA>=10) %>%
              mutate(DistToTSS = abs(DistToTSS),
              chromHMM = factor(chromHMM),
              NearestDHStype = factor(NearestDHStype),
              DistToNearestDHShot2 = abs(DistToNearestDHShot2),  
              repliSeq = factor(repliSeq),
              genomicDivision = factor(genomicDivision),
              DistToNearestDHS = abs(DistToNearestDHS),
              DistToCGI = abs(DistToCGI),
              DistToNearestCTCF = abs(DistToNearestCTCF),
              distToDHScluster =abs(distToDHScluster),
              H2AFZ= abs(H2AFZ),
              H3K27ac = abs(H3K27ac), 
              H3K27me3 = abs(H3K27me3), 
              H3K36me3 = abs(H3K36me3), 
              H3K4me1 =abs(H3K4me1), 
              H3K4me2 = abs(H3K4me2), 
              H3K4me3 = abs(H3K4me3),
              H3K79me2 = abs(H3K79me2), 
              H3K9ac =abs(H3K9ac), 
              H3K9me1 =abs(H3K9me1), 
              H3K9me3 =abs(H3K9me3), 
              H4K20me1 =abs(H4K20me1), 
              HDAC1 =abs(HDAC1), 
              HDAC2 =abs(HDAC2), 
              HDAC6 =abs(HDAC6))

LMfeat <- glm(expression ~ #log10(DistToCGI+1) + 
              log10(DistToNearestDHS+1)+ 
              nDHS100kb ,
            #log10(nDHS100kbDensity+1) + 
            #AB_compartment_value +
            ##log10(meanDist5DHS+1) +
            #log10(H3K9ac+1) + 
            #log10(H3K27ac+1) +
            #log10(H3K79me2+1) + 
            #log10(H3K4me2+1) +
            #log10(H2AFZ+1) + log10(H3K36me3+1) + log10(H3K36me3+1) + log10(H3K4me1+1) + log10(H3K9me1+1) + log10(H3K9me3+1) +
            #log10(H4K20me1+1) + log10(HDAC1+1) + log10(HDAC2+1) + log10(HDAC6+1) + repliSeq, 
             data = featureLM)

K562$predExp <- predict(LMfeat, K562)
#pdf(paste0('~/public_html/blog/2017Nov27/',opt$SampleName,'_predExpK562.pdf'))
#scatter.smooth(K562$predExp, log10(K562$value+1), col='red', frame=F, xlab='Predicted expression', ylab='log10(RPKM+1)', main='CMV\nPredicted expression vs gencode expression')
#abline(h=1, lty=2)
#abline(v=0.5, lty=2)
#text(labels='High predicted expression\nand high expression', x=0.8, y=4)
#dev.off()
#
#K562_expandPred <- K562[K562$value>=10 & K562$predExp>=0.5,]
#
#if (nrow(K562_expandPred)>0) {
#       write.table(subset(K562_expandPred, select=c(chrom, start,end,BC, value, strand)), file=paste0(opt$OUTDIR,'/K562PredictedExpressed.bed'), sep='\t', col.names=F, row.names=F, quote=F)
#}
#
#
#pdf(paste0('~/public_html/blog/2017Nov27/',opt$SampleName,'_predExpK562_repliSeq.pdf'))
#barplot(table(K562_expandPred$repliSeq), xlab='Replication phase', ylab='Number of genes', main='CMV\nPredicted expressed and expressed genes \nand replication phase')
#dev.off()
#
#
#pdf(paste0('~/public_html/blog/2017Nov27/',opt$SampleName,'_predExpK562_AB.pdf'))
#barplot(table(K562_expandPred$ABcomp), xlab='A/B compartment', ylab='Number of genes', main='CMV\nPredicted expressed and expressed genes \nand A/B compartment')
#dev.off()

######
#GAM model
#######
cat('Plot AUROC for GAM features\n')

GAMfeat <- gam(Expression ~ s(log10(DistToNearestDHS+1)) + 
       #s(log10(DistToCGI+1)) + 
       s(nDHS100kb),
       #s(log10(DistToNearestDHShot2+1)) +
       #s(log10(nDHS100kbDensity+1))+
     #s(AB_compartment_value) +
     
     #(CTCF_01) +
     #(CTCF_02) +
     #(CTCF_03) +
     #(CTCF_04) +
     #(CTCF_05) +
     #
     #(Distal_01) +
     #(Distal_02) +
     #(Distal_03) +
     #(Distal_04) +
     #(Distal_05) +
     #
     #(Promoter_01) +
     #(Promoter_02) +
     #(Promoter_03) +
     #(Promoter_04) +
     #(Promoter_05) +
      #s(log10(meanDist5DHS+1))+
      ##s(log10(H3K9ac+1)) +
      #s(log10(H3K27ac+1)) +
      #s(log10(H3K79me2+1)) +
      #s(log10(H3K4me2+1)) +
      #s(log10(H2AFZ+1))+
      #s(log10(H3K36me3+1))+
      #s(log10(H3K4me1+1))+
      ##s(log10(H3K9me1+1))+
      #s(log10(H3K9me3+1))+ 
      #s(log10(H4K20me1+1))+ 
      ##s(log10(HDAC1+1))+ 
      #s(log10(HDAC2+1))+ 
      #s(log10(HDAC6+1))+
      #s(repliSeq, bs='re') ,
       #s(chromHMM, bs='re') +
       #s(NearestDHStype, bs='re'),
       data = featureGAM.clean, 
       method="REML", 
       family=binomial, select=T ) 


preds <- prediction(GAMfeat$fitted.values, GAMfeat$y)
perf <- performance(preds, "auc")
rocC <- round(as.numeric(slot(perf, "y.values")), 3)
rocC

GAM_plot <- performance(preds, "tpr", "fpr")

GAM_plot.ROC <- data.frame(y=GAM_plot@y.values[[1]] ,  x=GAM_plot@x.values[[1]], SampleName=unique(normdataiPCR$SampleName))

AUC.features <- ggplot(GAM_plot.ROC ,aes(x=x, y=y, colour=SampleName)) +
       geom_line(size = 2)+
       theme_classic()+
       ggtitle(opt$SampleName,"  \nAUC for bottom and top quartile")+
       scale_colour_brewer(palette='Set1') +
       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
       xlab('1-Specificity')+
       ylab('Sensitivity')+
       annotate('text', label = paste0('AUC = ',round(as.numeric(slot(perf, "y.values")), 3)), x = 0.5, y = 0.5,  color = 'black', size=10)+
       theme(legend.position="none")

pdf(paste0(opt$OUTDIR,'/graphs/GAM_features.pdf'), height=5, width=5) 
AUC.features
dev.off()

nums <- sapply(featureGAM.clean, is.numeric)

numFeatuers <- featureGAM.clean[,nums]


logCol <- c("DistToTSS", "DistToNearestDHS", "DistToNearestCTCF", 
"DistToCGI",  
"nDHS100kbDensity", "nDHS300kbsum", "nDHS500kbDensity", "nDHS1MbDensity", 
"CTCFPksDens",
 "meanDist5DHS", "distToDHScluster", 
"H2AFZ", "H3K27ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me2", 
"H3K4me3", "H3K79me2", "H3K9ac", "H3K9me1", "H3K9me3", "H4K20me1", 
"HDAC1", "HDAC2", "HDAC6", "totalTags50kb")

numFeatuers[,logCol] <- log10(abs(numFeatuers[,logCol])+1)
numFeatuers <- subset(numFeatuers, select=-c(DNA, RNA, DistToGeneDesert, DistTok562GeneDesert))
library(gridExtra)
ggList <- list()
for (i in 1:(ncol(numFeatuers)-2)) {
       plotData <- numFeatuers[,c('DistToNearestDHS',colnames(numFeatuers)[i], 'Expression')]
       colnames(plotData) <-c('DistToNearestDHS','feature', 'Expression')
       ggList[[i]] <- ggplot(plotData, aes(x=DistToNearestDHS, y=feature, fill=as.factor(Expression))) +
                     geom_point(pch=21) +
                     theme_classic() +
                     scale_fill_brewer(palette='Set1')+
                     #facet_wrap(~as.factor(Expression))+
                     #geom_density2d(aes(colour=as.factor(Expression)), size=1) +
                     scale_colour_brewer(palette='Set1') + 
                     stat_ellipse(aes(colour=as.factor(Expression)), size=2)+
                     xlab('Distance to DHS') +
                     ylab(colnames(numFeatuers)[i])+
                     ggtitle(colnames(numFeatuers)[i])+
                     theme(legend.position='none')
}

#ggDHSfeatures <-do.call("grid.arrange", c(ggList, ncol=6))
#
#pdf(paste0(opt$OUTDIR,'/graphs/DHSdistanceFeatures.pdf'), height=10, width=15) 
#grid.draw(ggDHSfeatures)
#dev.off()
#
#####
#Expression features
#####
ggList <- list()
for (i in 1:(ncol(numFeatuers)-2)) {
       plotData <- numFeatuers[,c('expression',colnames(numFeatuers)[i], 'Expression')]
       colnames(plotData) <-c('expression','feature', 'Expression')
       ggList[[i]] <- ggplot(plotData, aes(y=feature, x=as.factor(Expression), fill=as.factor(Expression))) +
                     #geom_point(pch=21) +
                     geom_boxplot(outlier.shape=NA) +
                     theme_classic() +
                     scale_fill_brewer(palette='Set1')+
                     #facet_wrap(~as.factor(Expression))+
                     #geom_density2d(aes(colour=as.factor(Expression)), size=1) +
                     annotate('text', label=paste0('p = ',signif(t.test(plotData$feature[plotData$Expression==1], plotData$feature[plotData$Expression==0])$p.val,1)), x=1.5, y=max(plotData$feature)*0.95)+
                     scale_colour_brewer(palette='Set1') + 
                     stat_ellipse(aes(colour=as.factor(Expression)), size=2)+
                     xlab('Reporter activity') +
                     ylab(colnames(numFeatuers)[i])+
                     ggtitle(colnames(numFeatuers)[i])+
                     theme(legend.position='none')
}








#ggDHSfeatures <-do.call("grid.arrange", c(ggList, ncol=8))
#
#pdf(paste0(opt$OUTDIR,'/graphs/ExpressioneFeatures.pdf'), height=10, width=15) 
#grid.draw(ggDHSfeatures)
#dev.off()
#dhsTypeMCV <- rea
#dhsTypeMCV <- read.table(paste0(opt$OUTDIR,'/', opt$SampleName,'_MCV.bed'), sep='\t', stringsAsFactors=F)
#dhsTypeMCV$MCV <- rep(paste0('MCV_',formatC(1:20, width = 2, format = "d", flag = "0")), length(unique(dhsTypeMCV$V1)))
#dhsTypeMCVROC <- dcast(data = dhsTypeMCV,formula = V1~MCV,fun.aggregate = sum,value.var = "V3")
#dhsTypeMCVROC <- dhsTypeMCVROC[dhsTypeMCVROC$V1 %in% normdataiPCR$BC,]
#rocMCV <- ROCExpression(dhsTypeMCVROC, 1)
#dhsTypeMCVROCglm <- dhsTypeMCVROC %>%
#       filter(V1 %in% normdataiPCR$BC) %>%
#       mutate(Expression = normdataiPCR[match(V1, normdataiPCR$BC),]$expression) %>%
#       mutate(Expression = ifelse(Expression>as.numeric(summary(Expression)[5]), 1, ifelse(Expression<(as.numeric(summary(Expression)[2])), 0,NA))) %>%
#       filter(!is.na(Expression))
#fit <- glm(formula=Expression ~. ,data=dhsTypeMCVROCglm[,-1], family='binomial')
#preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
#perf <- performance(preds, "auc")
#rocC <- round(as.numeric(slot(perf, "y.values")), 3)
#rocMCV[,21] <-rocC
#colnames(rocMCV)[21] <- 'Nearest_20'
#
#pdf(paste0(opt$OUTDIR,'/graphs/AUROC_MCVnearest20.pdf'), height=5, width=10)
#rocMCV %>%
#       melt() %>%
#       mutate(SampleName= opt$SampleName) %>%
#       ggplot(aes(x=variable, y=value, fill=SampleName)) +
#       geom_bar(stat='identity', position=position_dodge(), colour='black', alpha=0.6) +
#       xlab('DHS neighbour') +
#       ylab('AUROC') +
#       coord_cartesian(ylim = c(0.5, 0.9)) +
#       scale_fill_brewer(palette='Set1') +
#       theme_classic() +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, ,hjust=0.95,vjust=0.95), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
#       ggtitle(opt$SampleName, subtitle='AUROC of MCV for nearest 5 DHS type individually\nas well as combined') #+
#dev.off()
#
#
#cat('MCV of binned distance\n')
#MCVdist <- read.table(paste0(opt$OUTDIR,'/', opt$SampleName,'_Dist_MCV.bed'), sep='\t', stringsAsFactors=F)
#
#MCVdistSplit <- split.data.frame(MCVdist, MCVdist$V5)
##loci <- MCVdistSplit[[1]]
#lociDistance <- lapply(MCVdistSplit, function(loci) { 
#       loci$distance <- 'NA'
#       if (nrow(loci)==10) {
#             loci$distance <- c(seq(from=-5e4, to=-1e4, by=1e4), seq(from=10000, to=5e4, by=1e4))
#             loci
#       }
#})
#
#lociDistanceDF <-  rbindlist(lociDistance)
#lociDistance <- lociDistanceDF
#lociDistance$V6 <- as.numeric(lociDistance$V6)
#lociDistance$distance <- abs(lociDistance$distance)
#lociDistance[lociDistance$V6=='NaN',]$V6 <- 0
#
#MCVdistLoci <- as.data.frame(dcast(data = lociDistance,formula = V5~distance,fun.aggregate = sum,value.var = "V6"))
#rocMCVbin <- ROCExpression(MCVdistLoci, 1)
#
#
#rocMCVbinglm <- MCVdistLoci %>%
#       filter(V5 %in% normdataiPCR$BC) %>%
#       mutate(Expression = normdataiPCR[match(V5, normdataiPCR$BC),]$expression) %>%
#       mutate(Expression = ifelse(Expression>as.numeric(summary(Expression)[5]), 1, ifelse(Expression<(as.numeric(summary(Expression)[2])), 0,NA))) %>%
#       filter(!is.na(Expression))
#fit <- glm(formula=Expression ~. ,data=rocMCVbinglm[,-1], family='binomial')
#preds <- prediction(as.numeric(fit$fitted.values), as.numeric(fit$y))
#perf <- performance(preds, "auc")
#rocC <- round(as.numeric(slot(perf, "y.values")), 3)
#rocMCVbin[,6] <-rocC
#colnames(rocMCVbin)[6] <- 'First_50kb'
#
#pdf(paste0(opt$OUTDIR,'/graphs/AUROC_MCVbinned100kb.pdf'), height=5, width=5)
#rocMCVbin %>%
#       melt() %>%
#       mutate(SampleName= opt$SampleName) %>%
#       ggplot(aes(x=variable, y=value, fill=SampleName)) +
#       geom_bar(stat='identity', position=position_dodge(), colour='black', alpha=0.6) +
#       xlab('Binned distance from insertions') +
#       ylab('AUROC') +
#       coord_cartesian(ylim = c(0.5, 0.9)) +
#       scale_fill_brewer(palette='Set1') +
#       theme_classic() +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12, angle=60, ,hjust=0.95,vjust=0.95), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12)) +
#       ggtitle(opt$SampleName, subtitle='AUROC of max MCV for surrounding DHS') #+
#dev.off()




#cat('NonCTCF of binned distance\n')
#nonCTCFdist <- read.table(paste0(opt$OUTDIR,'/', opt$SampleName,'_Dist_nDHSnonCTCF.bed'), sep='\t', stringsAsFactors=F)
#
#nonCTCFdistSplit <- split.data.frame(nonCTCFdist, nonCTCFdist$V5)
##loci <- nonCTCFdistSplit[[1]]
#lociDistance <- lapply(nonCTCFdistSplit, function(loci) { 
#       loci$distance <- 'NA'
#       if (nrow(loci)==10) {
#             loci$distance <- c(seq(from=-5e4, to=-1e4, by=1e4), seq(from=10000, to=5e4, by=1e4))
#             loci
#       }
#})
#
#lociDistanceDF <-  rbindlist(lociDistance)
#lociDistance <- lociDistanceDF
#lociDistance$V6 <- as.numeric(lociDistance$V6)
#lociDistance$distance <- abs(lociDistance$distance)
#lociDistance[lociDistance$V6=='NaN',]$V6 <- 0
#
#nonCTCFdistLoci <- as.data.frame(dcast(data = lociDistance,formula = V4~distance,fun.aggregate = sum,value.var = "V6"))
#rocnonCTCFbin <- ROCExpression(nonCTCFdistLoci, 1)
#



#RocFeature.Plot <- ggdraw() + 
#       draw_plot(AUC.features , x = 0, y = 0, width = .33, height = 1) +
#       draw_plot(Imp.features, x = .33, y = 0, width = .66, height = 1) 
#RocFeature.title <-ggdraw() +draw_label("ROC and feature importance",fontface='bold')
#plot_grid(RocFeature.title , RocFeature.Plot, ncol=1, rel_heights=c(0.1, 1))
#

#dev.off()

######
#Contacts from contact bed file
######
#####
#Find the bump
#####
#library(MESS)
#library(data.table) 
#library(parallel)
#contacts <-fread(paste0(opt$OUTDIR,'/',opt$SampleName,'_1Mb_contact_loess_annotated.bed'), header=T)
#
#
#
#
#
#bumpContacts <- mysplit.data.frame(contacts, contacts$BC)
#
#
#chromList <-  mysplit.data.frame(contacts, contacts$chrom)
#
#chromOutput <- lapply(chromList, function(chrom) {
#              loci <- unique(data.frame(distance=chrom$distance,
#                            lowess=chrom$lowess,
#                            lowessPercentage=chrom$lowessPercentage))
#                     loci <- loci[order(loci$distance),]
#                     loci
#              })
#
#
#AUClist <-list()
#for (i in 1:length(bumpContacts)) {cat(unique(bumpContacts[[i]]$BC),'\n')
#       
#
#       
#       lociZero <- data.frame(matrix(ncol=13, nrow=401))
#                     colnames(lociZero) <- c("chrom", "start", "end", "BC", "normContact", "distance", "percentageOfMax", 
#                                                 "lowess", "ObseExp", "lowessPercentage", "ObsExpPercentage", 
#                                                 "DHSdensity", "CTCFdensity")
#                     lociZero$BC <- unique(bumpContacts[[i]]$BC)
#                     lociZero$chrom <- unique(bumpContacts[[i]]$chrom)
#                     lociZero$normContact <- 0
#                     lociZero$percentageOfMax <- 0
#                     lociZero$ObseExp <- 0
#                     lociZero$ObsExpPercentage <- 0
#                     lociZero$distance <- seq(from=-1e6, to=1e6, by=5000)
#                     lociZero$DHSdensity <- 0
#                     lociZero$CTCFdensity <- 0
#                     lociZero$lowess <- chromOutput[[unique(bumpContacts[[i]]$chrom)]]$lowess
#                     lociZero$lowessPercentage <- chromOutput[[unique(bumpContacts[[i]]$chrom)]]$lowessPercentage
#                     
#                     
#       bumpContacts[[i]] <- rbind(bumpContacts[[i]],lociZero[!lociZero$distance %in% bumpContacts[[i]]$distance,])
#       bumpContacts[[i]] <- bumpContacts[[i]][order(bumpContacts[[i]]$distance),]
#                                
#       AUClowewss <- auc(x = subset(bumpContacts[[1]], distance > 1e5| distance < -1e5)$distance, y = subset(bumpContacts[[1]], distance > 1e5| distance < -1e5)$lowessPercentage)
#       
#       bumpContacts[[i]]$smoothSpline <- predict(smooth.spline(bumpContacts[[i]]$distance, y=bumpContacts[[i]]$percentageOfMax))$y
#
#       bumpMinusCenter <- subset(bumpContacts[[i]], distance > 1e5 | distance < -1e5)
#       
#
#       
#       
#
#      if (any(head(bumpMinusCenter[rev(order(bumpMinusCenter$smoothSpline)),],10)$DHSdensity>0 & head(bumpMinusCenter[rev(order(bumpMinusCenter$smoothSpline)),],10)$ObsExpPercentage>5)) { top25DHS <- TRUE 
#      } else {top25DHS <- FALSE}
#
#       
#      
#       aucBC <- auc(x = bumpMinusCenter$distance, y = bumpMinusCenter$smoothSpline)
#      
#      
#       
#       AUClist[[i]] <- data.frame(BC=unique(bumpContacts[[i]]$BC),
#              AUC= aucBC,
#              AUCvsAUSlowess = aucBC/AUClowewss,
#              topDHS = top25DHS)
#              #nrow(head(bumpContacts[[i]][bumpContacts[[i]]$distance > 1e5 | bumpContacts[[i]]$distance < -1e5,][rev(order(bumpContacts[[i]][bumpContacts[[i]]$distance > 1e5 | bumpContacts[[i]]$distance < -1e5,]$smoothSpline)),],20)[head(bumpContacts[[i]][bumpContacts[[i]]$distance > 1e5 | bumpContacts[[i]]$distance < -1e5,][rev(order(bumpContacts[[i]][bumpContacts[[i]]$distance > 1e5 | bumpContacts[[i]]$distance < -1e5,]$smoothSpline)),],20)$DHSdensity > 0,]))
#
#}
#
#aucBC <- as.data.frame(rbindlist(AUClist))
#aucBC <- aucBC[rev(order(aucBC$AUCvsAUSlowess)),]
#
#aucBC <- aucBC[aucBC$topDHS==TRUE,]
#pdf(paste0(opt$OUTDIR,"/graphs/BumpHunting_smoothed.pdf"), width=20, height=10, )
##aucBC <- aucBC[order(aucBC$AUCvsAUSlowess),]
#
#par(mfrow=c(5,10) , mar=c(1,2,1,1))
#for (i in 1:50) {
#
#       #BCcontact <- contacts[contacts$BC == aucBC$BC[i],]
#       BCcontact <- bumpContacts[[aucBC$BC[i]]]
#       BCcontact$smoothSpline <- predict(smooth.spline(BCcontact$distance, y=BCcontact$percentageOfMax, spar=0.40))$y
#       plot(BCcontact$distance, BCcontact$smoothSpline, type='l', axes=FALSE, frame=F, main=aucBC$BC[i])
#       Axis(side=1, labels=FALSE)
#       Axis(side=2, labels=TRUE)
#       lines(BCcontact$distance, BCcontact$lowessPercentage, col='blue', lty='dashed')
#       lines(BCcontact[BCcontact$ObsExpPercentage >5 & BCcontact$DHSdensity > 1 ,]$distance, BCcontact[BCcontact$ObsExpPercentage >5 & BCcontact$DHSdensity > 1 ,]$smoothSpline, type='p', col='red', pch=19)
#}
#
#dev.off()
#
#
#
#pdf(paste0(opt$OUTDIR,"/graphs/BumpHunting_smoothed_ObsVsExp_annotated.pdf"), width=20, height=10, )
##aucBC <- aucBC[order(aucBC$AUCvsAUSlowess),]
#
#par(mfrow=c(5,10) , mar=c(1,2,1,1))
#for (i in 1:50) {
#
#       #BCcontact <- contacts[contacts$BC == aucBC$BC[i],]
#       BCcontact <- bumpContacts[[aucBC$BC[i]]]
#       BCcontact$smoothSpline <- predict(smooth.spline(BCcontact$distance, y=BCcontact$ObsExpPercentage, spar=0.40))$y
#       plot(BCcontact$distance, BCcontact$smoothSpline, type='l', axes=FALSE, frame=F, main=aucBC$BC[i])
#       Axis(side=1, labels=FALSE)
#       Axis(side=2, labels=TRUE)
#       abline(v=1e5, lty='dashed')
#       abline(v=-1e5, lty='dashed')
#       #lines(BCcontact$distance, BCcontact$lowessPercentage, col='blue', axes=FALSE, frame=F, lty='dashed')
#       lines(BCcontact[BCcontact$ObsExpPercentage >5 & BCcontact$DHSdensity > 1 ,]$distance, BCcontact[BCcontact$ObsExpPercentage >5 & BCcontact$DHSdensity > 1 ,]$smoothSpline, type='p', col='springgreen2', pch=19)
#       lines(BCcontact[BCcontact$ObsExpPercentage >5 & BCcontact$CTCFdensity > 1 ,]$distance, BCcontact[BCcontact$ObsExpPercentage >5 & BCcontact$CTCFdensity > 1 ,]$smoothSpline, type='p', col='red', pch=19)
##       lines(BCcontact[BCcontact$ObsExpPercentage >5 & BCcontact$CTCFdensity > 1 & BCcontact$DHSdensity > 1,]$distance, BCcontact[BCcontact$ObsExpPercentage >5 & BCcontact$CTCFdensity > 1 & BCcontact$DHSdensity > 1,]$smoothSpline, type='p', col='red', pch=19)
#       
#       if (nrow(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$DHSdensity > 1 ,]) > 0) {
#              lines(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$DHSdensity > 1 ,]$distance, BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$DHSdensity > 1 ,]$smoothSpline, type='p', col='springgreen2', pch=19) }
#       if (nrow(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$CTCFdensity > 1 ,]) > 0 ) {
#              lines(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$CTCFdensity > 1 ,]$distance, BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$CTCFdensity > 1 ,]$smoothSpline, type='p', col='red', pch=19)}
#      # if (nrow(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$CTCFdensity > 1 & BCcontact$DHSdensity > 1,]) > 0) {
#      #        lines(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$DHSdensity > 1 & BCcontact$CTCFdensity > 1 ,]$distance, BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$CTCFdensity > 1 & BCcontact$DHSdensity > 1,]$smoothSpline, type='p', col='red', pch=19)}
#}
#
#dev.off()
#
#writePeaks <- aucBC[aucBC$AUCvsAUSlowess>1.5,]
#writePeaks$rank <- formatC(1:nrow(writePeaks), width = 3, format = "d", flag = "0")
#write.table(writePeaks, file=paste0(opt$OUTDIR,'/',opt$SampleName,'_topbumpsAucAucLowess1.5.tsv'), sep='\t' , quote=F, col.names=F, row.names=F)
#
#AUC_CTCF <- list() 
#for (i in 1:nrow(aucBC)) {
#       BCcontact <- contacts[contacts$BC == aucBC$BC[i],]
#       AUC_CTCF[[i]] <- data.frame(BC=unique(BCcontact$BC),
#                                   AUCvsAUCloss = aucBC$AUCvsAUSlowess[i],
#                                   DHS=nrow(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$DHSdensity > 1 ,]),
#                                   CTCF=nrow(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$CTCFdensity > 1 & BCcontact$DHSdensity > 1,]),
#                                   DHS_CTCF=nrow(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$CTCFdensity > 1 & BCcontact$DHSdensity > 1,]) )
#}
#AUC_CTCF <- rbindlist(AUC_CTCF)
#
##
##pdf('~/public_html/blog/2017Aug28/ProprtionDHS.pdf', height=5, width=5)
##barplot(c(sum(AUC_CTCF[AUC_CTCF$AUCvsAUCloss>1.5,]$CTCF)/nrow(AUC_CTCF[AUC_CTCF$AUCvsAUCloss>1.5,]), sum(AUC_CTCF[AUC_CTCF$AUCvsAUCloss<1.5,]$CTCF)/nrow(AUC_CTCF[AUC_CTCF$AUCvsAUCloss<1.5,]),46385/205899),
## names=c('FCauc>1.5','FCauc<1.5', 'All 15kb windows'),  ylab='Proportion windows having a CTCF' , main='Proportion of insertion windows +-5kb \n(i.e. 15kb around insertion) having a CTCF')
## dev.off()
## 
##pdf('~/public_html/blog/2017Aug28/Hist_FCauc.pdf', height=5, width=6) 
##hist(AUC_CTCF$AUCvsAUCloss, n=100, xlab='FCauc', main='Distribution of FCauc for filtered insertions')
##dev.off()
##
#
##pdf('~/public_html/blog/2017Aug28/FCauc_vsExpression.pdf', height=5, width=6)
##AUC_CTCF %>%
##       left_join(normdataiPCR, by='BC') %>%
##       select(BC, AUCvsAUCloss, RNA, DNA) %>%
##       filter( ! is.na(RNA)) %>%
##       ggplot(aes(x=AUCvsAUCloss, y=log10(RNA/DNA))) +
##       geom_point() +
##       theme_classic() +
##       geom_smooth(method='loess')
##dev.off()
##
#
#
#
#####
##Find the peaks
#####
#BCpeaks <- as.data.frame(rbindlist(AUClist))
#BCpeaks <- BCpeaks[rev(order(BCpeaks$AUCvsAUSlowess)),]
#
#peakWindow <- seq(from=-10, to=10, by=1)
#pdf("~/public_html/blog/2017Sep05/Peaks_diffs.pdf", width=20, height=10 )
#par(mfrow=c(5,10) , mar=c(1,2,1,1))
#for (i in 1:50) {
#
#       #BCcontact <- contacts[contacts$BC == aucBC$BC[i],]
#       BCcontact <- bumpContacts[[BCpeaks$BC[i]]]
#       BCcontact$smoothSpline <- predict(smooth.spline(BCcontact$distance, y=BCcontact$ObsExpPercentage, spar=0.3))$y
#       Peaks.pred <- which(diff(sign(diff(BCcontact$smoothSpline)))==-2)+1
#       #Peaks.pred <-c((which(diff(sign(diff(BCcontact$smoothSpline)))==-2)+1)-2, (which(diff(sign(diff(BCcontact$smoothSpline)))==-2)+1)-1,(which(diff(sign(diff(BCcontact$smoothSpline)))==-2)+1), (which(diff(sign(diff(BCcontact$smoothSpline)))==-2)+1)+1, (which(diff(sign(diff(BCcontact$smoothSpline)))==-2)+1)+2)
#       Peak.width.list <- list()
#              
#              for (j in 1:length(peakWindow)) {
#              peakWidth <- peakWindow[j]
#              Peak.width.list[[j]] <- Peaks.pred+peakWidth 
#              Peak.width.list[[j]][Peak.width.list[[j]]<0] <- 0
#              Peak.width.list[[j]] <- unique(Peak.width.list[[j       ]])}
#       Peaks.contact <- BCcontact[unlist(Peak.width.list),]
#       
#       
#       #Peaks.contact <- BCcontact[ Peaks.pred ,]
#      # Peaks.contact <- Peaks.contact[Peaks.contact$smoothSpline >mean(BCcontact$smoothSpline),]
#       plot(BCcontact$distance, BCcontact$smoothSpline, type='l', axes=FALSE, frame=F, main=BCpeaks$BC[i], ylim=c(0, max(BCcontact$smoothSpline)))
#       Axis(side=1, labels=FALSE)
#       Axis(side=2, labels=TRUE)
#       abline(v=1e5, lty='dashed')
#       abline(v=-1e5, lty='dashed')
#       #abline(h=mean(BCcontact$smoothSpline), lty='dashed')
#       lines(Peaks.contact[Peaks.contact$ObsExpPercentage >5 & Peaks.contact$DHSdensity > 1 ,]$distance, Peaks.contact[Peaks.contact$ObsExpPercentage >5 & Peaks.contact$DHSdensity > 1 ,]$smoothSpline, type='p', col='springgreen2', pch=19)
#       lines(Peaks.contact[Peaks.contact$ObsExpPercentage >5 & Peaks.contact$CTCFdensity > 1 ,]$distance, Peaks.contact[Peaks.contact$ObsExpPercentage >5 & Peaks.contact$CTCFdensity > 1 ,]$smoothSpline, type='p', col='red', pch=19)
#       
#       if (nrow(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$DHSdensity > 1 ,]) > 0) {
#              lines(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$DHSdensity > 1 ,]$distance, BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$DHSdensity > 1 ,]$smoothSpline, type='p', col='springgreen2', pch=19) }
#       if (nrow(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$CTCFdensity > 1 ,]) > 0 ) {
#              lines(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$CTCFdensity > 1 ,]$distance, BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$CTCFdensity > 1 ,]$smoothSpline, type='p', col='red', pch=19)}
#     
#     
#       lines(Peaks.contact$distance, rep(0.02*max(BCcontact$smoothSpline),nrow(Peaks.contact)), type='h', col='red', pch=1, cex=2)
#       
#       
#       text(x=0, y=0.90* max(BCcontact$smoothSpline), label=paste0('npeaks = ',length(which(diff(sign(diff(BCcontact$smoothSpline)))==-2)+1)), cex=1.5)
#}
#dev.off()
#
#
#####
##Find boundaries i.e. where one side of the loci is more in contact than the other
#####
#AUCBoundarylist <-list()
#for (i in 1:length(bumpContacts)) {cat(unique(bumpContacts[[i]]$BC),'\n')
#       
#          lociZero <- data.frame(matrix(ncol=13, nrow=401))
#                     colnames(lociZero) <- c("chrom", "start", "end", "BC", "normContact", "distance", "percentageOfMax", 
#                                                 "lowess", "ObseExp", "lowessPercentage", "ObsExpPercentage", 
#                                                 "DHSdensity", "CTCFdensity")
#                     lociZero$BC <- unique(bumpContacts[[i]]$BC)
#                     lociZero$normContact <- 0
#                     lociZero$percentageOfMax <- 0
#                     lociZero$ObseExp <- 0
#                     lociZero$ObsExpPercentage <- 0
#                     lociZero$distance <- seq(from=-1e6, to=1e6, by=5000)
#                     lociZero$DHSdensity <- 0
#                     lociZero$CTCFdensity <- 0
#                     lociZero$lowess <- chromOutput[[unique(bumpContacts[[i]]$chrom)]]$lowess
#                     lociZero$lowessPercentage <- chromOutput[[unique(bumpContacts[[i]]$chrom)]]$lowessPercentage
#                     
#                     
#       bumpContacts[[i]] <- rbind(bumpContacts[[i]],lociZero[!lociZero$distance %in% bumpContacts[[i]]$distance,])
#       bumpContacts[[i]] <- bumpContacts[[i]][order(bumpContacts[[i]]$distance),]
#
#       bumpContacts[[i]]$smoothSpline <- predict(smooth.spline(bumpContacts[[i]]$distance, y=bumpContacts[[i]]$ObsExpPercentage))$y
#
#       bumpMinusCenter <- subset(bumpContacts[[i]], distance > 1e5 | distance < -1e5)
#       
#       aucBCboth <- auc(x = bumpMinusCenter$distance, y = bumpMinusCenter$smoothSpline)
#       aucBCLeft <- auc(x = subset(bumpMinusCenter, distance<0)$distance, y = subset(bumpMinusCenter, distance<0)$smoothSpline)
#       aucBCRight <- auc(x = subset(bumpMinusCenter, distance>0)$distance, y = subset(bumpMinusCenter, distance>0)$smoothSpline)
#      
#       
#       AUCBoundarylist[[i]] <- data.frame(BC=unique(bumpContacts[[i]]$BC),
#              AUC= aucBCboth,
#              AUCvsAUSlowess = aucBCboth/AUClowewss,
#              AUCLeft = aucBCLeft,
#              AUCright = aucBCRight,
#              AUCratioLR = aucBCLeft/aucBCRight,
#              AUCratioNorm = abs(aucBCLeft/aucBCRight-1))
#}
#
#aucBoundary<- as.data.frame(rbindlist(AUCBoundarylist))
#
#aucBoundary <- aucBoundary[rev(order(aucBoundary$AUCratioNorm)),]
#
#
#par(mfrow=c(5,10) , mar=c(1,2,1,1))
#for (i in 1:50) {
#
#       #BCcontact <- contacts[contacts$BC == aucBC$BC[i],]
#       BCcontact <- bumpContacts[[aucBoundary$BC[i]]]
#       BCcontact$smoothSpline <- predict(smooth.spline(BCcontact$distance, y=BCcontact$ObsExpPercentage, spar=0.40))$y
#       plot(BCcontact$distance, BCcontact$smoothSpline, type='l', axes=FALSE, frame=F, main=aucBoundary$BC[i])
#       Axis(side=1, labels=FALSE)
#       Axis(side=2, labels=TRUE)
#        abline(v=1e5, lty='dashed')
#       abline(v=-1e5, lty='dashed')
#       #lines(BCcontact$distance, BCcontact$lowessPercentage, col='blue', lty='dashed')
#       lines(BCcontact[BCcontact$ObsExpPercentage >5 & BCcontact$DHSdensity > 1 ,]$distance, BCcontact[BCcontact$ObsExpPercentage >5 & BCcontact$DHSdensity > 1 ,]$smoothSpline, type='p', col='red', pch=19)
#}
#
#
######
##How peaks correlates with expression
#######
#peakWindow <- seq(from=-10, to=10, by=1)
#
#peakList <-list()
#for (i in 1:nrow(BCpeaks)) {
#       BCcontact <- bumpContacts[[BCpeaks$BC[i]]]
#       BCcontact$smoothSpline <- smooth.spline(BCcontact$distance, y=BCcontact$ObsExpPercentage, spar=0.1)$y
#       Peaks.pred <- which(diff(sign(diff(BCcontact$smoothSpline)))==-2)+1
#       BCpeaks$nPeaks[i] <- nrow(BCcontact[Peaks.pred ,][BCcontact[Peaks.pred ,]$smoothSpline>=2,])
#       
#       peakList[[i]] <- BCcontact[Peaks.pred ,][BCcontact[Peaks.pred ,]$smoothSpline>=2,]
#       
#       Peak.width.list <- list()
#       for (j in 1:length(peakWindow)) {
#              peakWidth <- peakWindow[j]
#              Peak.width.list[[j]] <- Peaks.pred+peakWidth 
#              Peak.width.list[[j]][Peak.width.list[[j]]<0] <- 0
#              Peak.width.list[[j]] <- unique(Peak.width.list[[j]])
#              }
#       Peaks.contact <- BCcontact[unlist(Peak.width.list),]
#       
#       
#       BCpeaks$nDHS[i] <- nrow(Peaks.contact[Peaks.contact$DHSdensity>0,])
#       BCpeaks$nCTCF[i] <- nrow(Peaks.contact[Peaks.contact$CTCFdensity>0,])
#       BCpeaks$DHSaccessibility[i]  <- sum(Peaks.contact$DHSdensity, na.rm=T)
#       BCpeaks$sumObsExp[i] <- sum(Peaks.contact$ObsExpPercentage, na.rm=T)
#       
#       if (nrow(Peaks.contact[Peaks.contact$DHSdensity>0,]) >0) {
#              BCpeaks$sumObsExpDHS[i] <- sum(Peaks.contact[Peaks.contact$DHSdensity>0,]$ObsExpPercentage , na.rm=T) 
#       } else { BCpeaks$sumObsExpDHS[i] <- 0 }
#}
#
#BCpeaks <- BCpeaks[rev(order(BCpeaks$nPeaks)),]
#
#writePeaks <- BCpeaks[c(1:50, 6001:6050, 12903:12952),]
#writePeaks$rank <- formatC(1:nrow(writePeaks), width = 3, format = "d", flag = "0")
#
#write.table(writePeaks, file=paste0(opt$OUTDIR,'/',opt$SampleName,'_HighMedLowPeaks.tsv'), sep='\t' , quote=F, col.names=F, row.names=F)
#
#
#######
##Compare peaks finding for different smooth options
#######
##High number of peaks
#
#SparTest <- c( 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
#peakWindow <- seq(from=-10, to=10, by=1)
#
#BCpeaksNoZero <- BCpeaks[BCpeaks$nPeaks>0,]
#
#
#set.seed(12345)
#for (i in 1:10) {
#       pdf(paste0("~/public_html/blog/2017Sep05/Peaks_HighMedLow_part",i,".pdf"), width=20, height=10 )
#       par(mfrow=c(6,10) , mar=c(1,2,1,1))
#       highPeaks <- sample(1:1000, 2)
#       medPeaks <- sample(4001:6000, 2)
#       lowPeaks <- sample(10000:12000, 2)
#       for (i in c(highPeaks, medPeaks, lowPeaks)) {
#              
#              
#              
#              BCcontact <- bumpContacts[[BCpeaksNoZero$BC[i]]]
#              #BCcontact$smoothSpline <- predict(smooth.spline(BCcontact$distance, y=BCcontact$percentageOfMax, spar=0.10))$y
#              plot(BCcontact$distance, BCcontact$percentageOfMax, type='l', axes=FALSE, frame=F, main='Percentage of max')
#              Axis(side=1, labels=FALSE)
#              Axis(side=2, labels=TRUE)
#              lines(BCcontact$distance, BCcontact$lowessPercentage, col='blue', lty='dashed')
#              
#              BCcontact <- bumpContacts[[BCpeaksNoZero$BC[i]]]
#            
#
#              Peaks.pred <- which(diff(sign(diff(BCcontact$ObsExpPercentage)))==-2)+1
#              Peaks.contact <- BCcontact[Peaks.pred ,][BCcontact[Peaks.pred ,]$ObsExpPercentage>=2,]
#               
#              plot(BCcontact$distance, BCcontact$ObsExpPercentage, type='l', axes=FALSE, frame=F, main=BCpeaksNoZero$BC[i], ylim=c(0, max(BCcontact$ObsExpPercentage)))
#                     Axis(side=1, labels=FALSE)
#                     Axis(side=2, labels=TRUE)
#                     abline(v=1e5, lty='dashed')
#                     abline(v=-1e5, lty='dashed')
#                     abline(h=2, lty='dashed')
#                     #lines(BCcontact$distance, BCcontact$lowessPercentage, col='blue', lty='dashed')
#                    # lines(BCcontact[BCcontact$ObsExpPercentage >5 & BCcontact$DHSdensity > 1 ,]$distance, BCcontact[BCcontact$ObsExpPercentage >5 & BCcontact$DHSdensity > 1 ,]$smoothSpline, type='p', col='red', pch=19)
#                     lines(Peaks.contact$distance, rep(0.02*max(BCcontact$ObsExpPercentage),nrow(Peaks.contact)), type='h', col='red', pch=1, cex=2)
#                     text(x=0, y=0.95* max(BCcontact$ObsExpPercentage), label=paste0('Raw Obs/Exp'), cex=1.5)
#                     text(x=0, y=0.85* max(BCcontact$ObsExpPercentage), label=paste0('npeaks = ', nrow(Peaks.contact)), cex=1.5)
#               
#              for (k in SparTest) { 
#       
#                     BCcontact$smoothSpline <- smooth.spline(BCcontact$distance, y=BCcontact$ObsExpPercentage, spar=k)$y
#                     
#                     
#                     Peaks.pred <- which(diff(sign(diff(BCcontact$smoothSpline)))==-2)+1
#                     cat(length(Peaks.pred),'\n')
#              
#
#                     Peaks.contact <- BCcontact[Peaks.pred ,][BCcontact[Peaks.pred ,]$smoothSpline>=2,]
#                     
#                     plot(BCcontact$distance, BCcontact$smoothSpline, type='l', axes=FALSE, frame=F, main=BCpeaksNoZero$BC[i], ylim=c(0, max(BCcontact$smoothSpline)))
#                     Axis(side=1, labels=FALSE)
#                     Axis(side=2, labels=TRUE)
#                     abline(v=1e5, lty='dashed')
#                     abline(v=-1e5, lty='dashed')
#                     abline(h=2, lty='dashed')
#                     #lines(BCcontact$distance, BCcontact$lowessPercentage, col='blue', lty='dashed')
#                    # lines(BCcontact[BCcontact$ObsExpPercentage >5 & BCcontact$DHSdensity > 1 ,]$distance, BCcontact[BCcontact$ObsExpPercentage >5 & BCcontact$DHSdensity > 1 ,]$smoothSpline, type='p', col='red', pch=19)
#                     lines(Peaks.contact$distance, rep(0.02*max(BCcontact$smoothSpline),nrow(Peaks.contact)), type='h', col='red', pch=1, cex=2)
#                     text(x=0, y=0.95* max(BCcontact$smoothSpline), label=paste0('spar = ',k), cex=1.5)
#                     text(x=0, y=0.85* max(BCcontact$smoothSpline), label=paste0('npeaks = ', nrow(Peaks.contact)), cex=1.5)
#              }
#       }
#       dev.off()
#}
#
#BCpeaksSpar <- BCpeaks
#pdf(paste0("~/public_html/blog/2017Sep05/SparNumberPeaks.pdf"), width=20, height=4 )
#par(mfrow=c(1,10))
#       for (i in 1:nrow(BCpeaks)) {
#       BCcontact <- bumpContacts[[BCpeaksSpar$BC[i]]]       
#       Peaks.pred <- which(diff(sign(diff(BCcontact$ObsExpPercentage)))==-2)+1
#       BCpeaksSpar$nPeaks[i] <- nrow(BCcontact[Peaks.pred ,][BCcontact[Peaks.pred ,]$smoothSpline>=2,])
#}
#hist(BCpeaksSpar$nPeaks, n=100, xlab='Number of peaks', ylab='Frequency', main=paste0('Raw Obs/Exp'), xlim=c(0, 129))
#
#BCpeaksSpar <- BCpeaks
#for (k in SparTest) {
#       for (i in 1:nrow(BCpeaks)) {
#              BCcontact <- bumpContacts[[BCpeaksSpar$BC[i]]]
#       
#                     BCcontact$smoothSpline <- smooth.spline(BCcontact$distance, y=BCcontact$ObsExpPercentage, spar=k)$y
#                     Peaks.pred <- which(diff(sign(diff(BCcontact$smoothSpline)))==-2)+1
#                     BCpeaksSpar$nPeaks[i] <- nrow(BCcontact[Peaks.pred ,][BCcontact[Peaks.pred ,]$smoothSpline>=2,])
#              }
#       hist(BCpeaksSpar$nPeaks, n=100, xlab='Number of peaks', ylab='Frequency', main=paste0('Spar = ',k), xlim=c(0, 45))
#       
#}
#dev.off()
#
#BCpeaksSpar <- BCpeaksSpar[rev(order(BCpeaksSpar$nPeaks)),]
#
#
#
#
#
#ggBCpeaks <- BCpeaks %>%
#       left_join(normdataiPCR, by='BC') %>%
#       mutate(Expression = log10(RNA/DNA),
#              allPeaksAccessibility = log10(DHSaccessibility+1),
#              dhsPeaksObsExp = log10(sumObsExpDHS+1),
#              allPeaksObsExp = log10(sumObsExp+1),
#              numerOfPeaks = nPeaks,
#              numberOfDHS = nDHS,
#              numberOfCTCF = nCTCF) %>%
#       select(numerOfPeaks, numberOfDHS, numberOfCTCF, Expression, allPeaksAccessibility, dhsPeaksObsExp, allPeaksObsExp)  %>%
#       melt(id='Expression')
#
#library(plyr)
#lm_eqn <- function(df){
#    m <- lm(value ~ Expression, df);
#    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
#         list(a = format(coef(m)[1], digits = 2), 
#              b = format(coef(m)[2], digits = 2), 
#             r2 = format(summary(m)$r.squared, digits = 3)))
#    as.character(as.expression(eq));                 
#}
#
#eq <- ddply(ggBCpeaks,.(variable),lm_eqn)
#
#eq$x <- c(6, 30, 25, 2, 2, 3)
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/threeDGenome.pdf"), width=14, height=8)
#ggplot(ggBCpeaks, aes(x=value, y=Expression, colour=variable)) +
#       geom_point( pch=21, colour=sampleColor) +
#       facet_wrap(~variable, scale='free_x') +
#       theme_classic() +
#       geom_smooth(method = "lm") +
#       geom_smooth(method = "loess", linetype='dashed') +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ylab('log10(RNA/DNA)')+
#       geom_text(data=eq,aes(x = x, y = 2,label=V1), parse = TRUE, inherit.aes=FALSE) 
#       annotation_logticks(sides = "l",long = unit(0.2, "cm"))
#dev.off()

####Find peaks from the following package with 0.20 of peak width
#library(pracma)
#pdf("~/public_html/blog/2017Sep05/FindPeaks.pdf", width=20, height=10 )
#par(mfrow=c(5,10) , mar=c(1,2,1,1))
#for (i in 1:50) {
#       BCcontact <- bumpContacts[[aucBC$BC[i]]]
#       BCcontact$smoothSpline <- predict(smooth.spline(BCcontact$distance, y=BCcontact$ObsExpPercentage, spar=0.40))$y
#       
#       BCpeaks <- findpeaks(BCcontact$smoothSpline)
#       BCpeaks <- as.data.frame(BCpeaks[BCpeaks[,1]>mean(BCcontact$smoothSpline),])
#       BCpeaks$length <- BCpeaks[,4]-BCpeaks[,3]
#       Peak.width <- data.frame(start=round(BCpeaks[,3]+(BCpeaks$length*.20)),end=round(BCpeaks[,4]-(BCpeaks$length*.20)))
#       
#       Peak.width.list <- list()
#              for (i in 1:nrow(Peak.width)) {Peak.width.list[[i]] <- Peak.width$start[i]:Peak.width$end[i]}
#       Peaks.contact <- BCcontact[unlist(Peak.width.list),]
#       
#       plot(BCcontact$distance, BCcontact$smoothSpline, type='l', axes=FALSE, frame=F, main=aucBC$BC[i], ylim=c(0, max(BCcontact$smoothSpline)))
#              Axis(side=1, labels=FALSE)
#              Axis(side=2, labels=TRUE)
#              abline(v=1e5, lty='dashed')
#              abline(v=-1e5, lty='dashed')
#              abline(h=mean(BCcontact$smoothSpline), lty='dashed')
#              
#              #abline(v=BCcontact$distance[BCpeaks$V2], lty='longdash', lwd=0.01, col='blue')
#              lines(BCcontact$distance[BCpeaks$V2], rep(0.04*max(BCcontact$smoothSpline),nrow(BCpeaks)), type='h', col='red', pch=1, cex=2, lwd=2)
#              segments(x0=BCcontact$distance[round(BCpeaks[,3]+(BCpeaks$length*.20))],y0=rep(0,nrow(BCpeaks)),x1= BCcontact$distance[round(BCpeaks[,4]-(BCpeaks$length*.20))], y1=rep(0,nrow(BCpeaks)), col='red', cex=2, lwd=2)######
#              
#              lines(Peaks.contact[Peaks.contact$ObsExpPercentage >5 & Peaks.contact$DHSdensity > 1 ,]$distance, Peaks.contact[Peaks.contact$ObsExpPercentage >5 & Peaks.contact$DHSdensity > 1 ,]$smoothSpline, type='p', col='springgreen2', pch=19)
#              lines(Peaks.contact[Peaks.contact$ObsExpPercentage >5 & Peaks.contact$CTCFdensity > 1 ,]$distance, Peaks.contact[Peaks.contact$ObsExpPercentage >5 & Peaks.contact$CTCFdensity > 1 ,]$smoothSpline, type='p', col='red', pch=19)
#              
#              
#              if (nrow(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$DHSdensity > 1 ,]) > 0) {
#                     lines(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$DHSdensity > 1 ,]$distance, BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$DHSdensity > 1 ,]$smoothSpline, type='p', col='springgreen2', pch=19) }
#              if (nrow(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$CTCFdensity > 1 ,]) > 0 ) {
#                     lines(BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$CTCFdensity > 1 ,]$distance, BCcontact[BCcontact$distance <= 5e3 & BCcontact$distance >= -5e3 & BCcontact$CTCFdensity > 1 ,]$smoothSpline, type='p', col='red', pch=19)}
#     
#     
#       
#}
#dev.off()
#


####
#running window
#####
#library(caTools)
#plot(BCcontact$distance, BCcontact$ObsExpPercentage, type='l', col='red'){
#par(mfrow=c(2, 4))
#for (i in c(1,2, 3, 5, 10, 15, 20, 40)) {
#       plot(BCcontact$distance, runmean(BCcontact$ObsExpPercentage, i), type='l')
#       
#}
#CTCF contacts
######
######
#CTCFContact <- read.table(paste0(opt$OUTDIR, "/Contacts/CTCF_Contact.txt"),sep=' ', stringsAsFactors=F)
#
##CTCFContact <- read.table('/tmp/CMV_CTCFContact.txt',sep=' ', stringsAsFactors=F)
#
#ggCTCFContacts <- normdataiPCR %>%
#       filter(BC %in% CTCFContact$V1) %>%
#       mutate(CTCFContact = CTCFContact[match(BC, CTCFContact$V1),]$V2) %>%
#       select(BC, DNA, RNA, SampleName, CTCFContact) %>%
#       mutate(bins = cut.pretty(CTCFContact,breaks=c(0, 1, 2, 5, 10, 20,100), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
#       mutate(Expression = log10(RNA/DNA)) %>%
#       ggplot(aes(x=bins, y=Expression, fill=SampleName)) +
#       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=1, colour=sampleColor) +
#       geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       theme_classic()+
#       #scale_fill_manual(values=sampleColor) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName,"  \nNumber of CTCF peaks in contact\nwith insertion (5kb resolution)")+
#       xlab('Number of CTCF peaks in contact (HiC)')+
#       ylab('log10(RNA/DNA)')+
#       coord_cartesian(ylim = c(-3, 3))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/CTCFContacts.pdf"), width=4, height=4)
#ggCTCFContacts
#dev.off()
#
#
#ggCTCFContactsAB <- normdataiPCR %>%
#       filter(BC %in% CTCFContact$V1) %>%
#       mutate(CTCFContact = CTCFContact[match(BC, CTCFContact$V1),]$V2) %>%
#       select(BC, DNA, RNA, SampleName, CTCFContact, AB_compartment_value) %>%
#       mutate(ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
#       filter(!is.na(ABcomp), ABcomp!='NA')  %>%
#       mutate(bins = cut.pretty(CTCFContact,breaks=c(0, 1, 2, 5, 10, 20,100), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
#       mutate(Expression = log10(RNA/DNA)) %>%
#       ggplot(aes(x=bins, y=Expression, fill=ABcomp)) +
#       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=1, colour=sampleColor) +
#       geom_boxplot(outlier.shape=NA,alpha=0.6)+ 
#       theme_classic()+
#       scale_fill_brewer(palette='Set1') +
#       facet_wrap(~ABcomp) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(paste0(opt$SampleName," A/B comp"),"  \nNumber of CTCF peaks in contact\nwith insertion (5kb resolution)")+
#       xlab('Number of CTCF peaks in contact (HiC)')+
#       ylab('log10(RNA/DNA)')+
#       coord_cartesian(ylim = c(-3, 3))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/CTCFContacts_ABcomp.pdf"), width=6, height=4)
#ggCTCFContactsAB
#dev.off()
#



####
#Creating regions between insertions
####
#TODO column names instead of numbers 
#percentin
iPCRregion <- subset(normdataiPCR,select=c(chrom, chromStart, chromEnd, DNA, RNA))
#Beginning of next insertion
iPCRregion$End <- c(iPCRregion$chromStart[-1], 'Test')
iPCRregion$EndChr <- c(iPCRregion$chrom[-1], 'Test')

#Signal for left and right insertion
iPCRregion$leftSignal <- c(log10(iPCRregion$RNA/iPCRregion$DNA))
iPCRregion$rightSignal <- c(log10(iPCRregion$RNA[-1]/iPCRregion$DNA[-1]), 'Test')

#Remove lines where the chromosomes are different
#This should remove the last name
iPCRregion <- iPCRregion[iPCRregion$chrom==iPCRregion$EndChr,]

iPCRregion <- subset(iPCRregion, select=c(chrom, chromStart, End, leftSignal, rightSignal))
#Remove last line
#TODO double check if we want to keep last name
iPCRregion <- iPCRregion[-nrow(iPCRregion),] 
colnames(iPCRregion) <- c('chrom','start','end','leftSignal','rightSignal')
iPCRregion <- iPCRregion[iPCRregion$start!=iPCRregion$end,]

write.table(iPCRregion,file=paste0(opt$OUTDIR,"/AllRegions.tsv"),sep='\t',col.names=T,quote=F,row.names=F)

cat("\nWriting output files\n")
cat("Saving -- "); print(date())
save(list=setdiff(ls(), c("basenames", basenames)), file=paste0(opt$OUTDIR,"/all.RData"), compress="bzip2")

cat("\nDone!!!\n")
print(date())


#############
#############
#Old contact code
#############
#############
######
#Contact maps
######
#cat('ContactMaps','\n')
#Contacts <- read.table(paste0(opt$OUTDIR,'/OBSvsExpcontactMap.tsv'),sep='\t', header=F, stringsAsFactors=F)
#ContactsCentered <- split(Contacts, factor(Contacts$V5))
#Contacts<- lapply(ContactsCentered, function(x) {
#       x$distance <- seq(from=300000, to=-300000, by=-5000)
#       x
#       })
#Contacts<-do.call(rbind.data.frame, Contacts)
#colnames(Contacts) <- c('chrom', 'start', 'end', 'contacts','BC','distance')
#Contacts[is.na(Contacts$contacts),]$contacts <- 0
##
##clustContacts <- subset(Contacts,select=c(BC, distance, contacts))
##clustContacts <- dcast(clustContacts, distance ~ BC)
##rownames(clustContacts) <- clustContacts$distance
##clustContacts <- subset(clustContacts, selec=-c(distance))
##Contacts$BC <- factor(Contacts$BC, levels=names(colSums(clustContacts>=1)[order(colSums(clustContacts>=1))])) 
#
#Contacts2<- Contacts
#Contacts2[log10(Contacts2$contacts)>2,]$contacts <- 2 
#Contacts2[log10(Contacts2$contacts+1) <= summary(log10(Contacts2$contacts+1))[5],]$contacts <- 0
##
#######
###Low and High maps
########
#
#High <-  tail(normdataiPCR[order(log10(normdataiPCR$RNA/normdataiPCR$DNA)),]$BC,1000)
#Low <- head(normdataiPCR[order(log10(normdataiPCR$RNA/normdataiPCR$DNA)),]$BC,1000)
#
#
#ggContactHigh <- ggplot(Contacts2[Contacts2$BC%in%High,], aes(x=distance, y=BC, fill=log10(contacts+1))) +
#geom_tile(aes(fill = log10(contacts+1))) + 
#scale_fill_gradient(low = "white",high = "red")+
#theme_classic()+
#scale_x_continuous(label= fancy_scientific,limits=c(-3e5,3e5))+
#ggtitle(paste0(opt$SampleName,' Normalised contact map'))+
#theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=0))
#
#test <- subset(normdataiPCR, select=c(BC, DNA, RNA))
#
#test$Expression <- log10(test$RNA/test$DNA+1)
#test$BC <- factor(test$BC, levels=test$BC[order(test$Expression)]) 
#ggTestHigh <- ggplot(test[test$BC%in%High,], aes(x=Expression, y=BC)) +
#       geom_point() +
#       theme_classic()+
#       ggtitle(paste0(opt$SampleName,' Normalised contact map'))+
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=0))+
#       coord_cartesian(xlim = c(0, 3))
#
#
#pdf(paste0(opt$OUTDIR,"/graphs/4C_contactMatrix_NEW_High.pdf"), width=8, height=10)
##ggMarginal(ggContact + theme_gray())
#grid.arrange(ggContactHigh, ggTestHigh, ncol = 4, layout_matrix = cbind(1, 1,  1, 2))
#dev.off()
#
#
#ggContactLow <- ggplot(Contacts2[Contacts2$BC%in%Low,], aes(x=distance, y=BC, fill=log10(contacts+1))) +
#geom_tile(aes(fill = log10(contacts+1))) + 
#scale_fill_gradient(low = "white",high = "red")+
#theme_classic()+
#scale_x_continuous(label= fancy_scientific,limits=c(-3e5,3e5))+
#ggtitle(paste0(opt$SampleName,' Normalised contact map'))+
#theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=0))
#
#test <- subset(normdataiPCR, select=c(BC, DNA, RNA))
#
#test$Expression <- log10(test$RNA/test$DNA+1)
#test$BC <- factor(test$BC, levels=test$BC[order(test$Expression)]) 
#ggTestLow <- ggplot(test[test$BC%in%Low,], aes(x=Expression, y=BC)) +
#       geom_point() +
#       theme_classic()+
#       ggtitle(paste0(opt$SampleName,' Normalised contact map'))+
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=0))+
#       #ylim(c(-2,2))+
#       coord_cartesian(xlim = c(0, 3))
#
#
#pdf(paste0(opt$OUTDIR,"/graphs/4C_contactMatrix_NEW_Low.pdf"), width=8, height=10)
#grid.arrange(ggContactLow, ggTestLow, ncol = 4, layout_matrix = cbind(1, 1,  1, 2))
#dev.off()
#
#
#
#####
##Number of obs/exp >2 contacts per loci for A/B compartments
#####
#cat('Number of Observed/Expected contacts >2 per +- 300kb around insertion','\n') 
#
#sigContacts <- Contacts %>% 
#       filter(contacts>2) %>% 
#       group_by(BC) %>%
#       count() 
#
#normdataiPCR$ObsExpover2<- sigContacts[match(normdataiPCR$BC,sigContacts$BC),]$n
#
#normdataiPCR$ObsExpover2[is.na(normdataiPCR$ObsExpover2)] <- 0
#
#sigContacts.test <- wilcox.test(normdataiPCR[normdataiPCR$AB_compartment_value>0,]$ObsExpover2, ABgroups[normdataiPCR$AB_compartment_value<0,]$ObsExpover2)
#sigContactspval <- sigContacts.test$p.value
#
#ggObsExpover2 <- normdataiPCR %>%
#       mutate(ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
#       #mutate(log10DNARNA = log10(RNA/DNA)) %>%
#       filter(!is.na(ABcomp), ABcomp!='NA') %>%
#       ggplot(aes(x = ABcomp, y = ObsExpover2, fill=ABcomp)) +
#       geom_boxplot(alpha=0.6, outlier.shape=NA) +
#       theme_classic()+
#       scale_fill_brewer(palette='Set1') +
#       annotate('text', label = paste0('Pvalue = ', format(as.numeric(sigContactspval),big.mark=",", trim=TRUE, digits=2)), x = 1.5, y = 40, color = 'black') +
#       annotate('text', label = paste0('nA = ', format(as.numeric(nrow(ABgroups[ABgroups$ABcomp=='A',])),big.mark=",", trim=TRUE)), x = 1, y = 35, color = 'black') +
#       annotate('text', label = paste0('nB = ', format(as.numeric(nrow(ABgroups[ABgroups$ABcomp=='B',])),big.mark=",", trim=TRUE)), x = 2, y = 35, color = 'black') +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName,"  \nA/B compartments")+
#       xlab('Genomic compartment')+
#       coord_cartesian(ylim = c(0, 45))+
#       ylab('Number of Obs/Exp>2 contacts\nwithin +-300kb')
##pdf(file=paste0("~/public_html/blog/2017Jul03/test.pdf"), width=4, height=4)
#pdf(file=paste0(opt$OUTDIR,"/graphs/AB_Contacts_ObsvsExpover2.pdf"), width=5, height=4)
#ggObsExpover2
#dev.off()
#
#

#####
#Compare DNAse in contact and in surrounding area of A/B comp insertions
#####
#DnaseContact <- read.table(paste0(opt$OUTDIR, "/DnaseContact.txt"),sep=' ', stringsAsFactors=F)
#
#ggDnaseContactsAB <- normdataiPCR %>%
#       select(BC, DNA, RNA, SampleName, contactDNase, AB_compartment_value) %>%
#       mutate(ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
#       filter(!is.na(ABcomp), ABcomp!='NA')  %>%
#       mutate(bins = cut.pretty(contactDNase,breaks=c(0, 1, 2, 5, 10, 20,100), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
#       mutate(Expression = log10(RNA/DNA)) %>%
#       ggplot(aes(x=bins, y=Expression, fill=ABcomp)) +
#       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=1, colour=sampleColor) +
#       geom_boxplot(outlier.shape=NA,alpha=0.6)+ 
#       theme_classic()+
#       scale_fill_brewer(palette='Set1') +
#       facet_wrap(~ABcomp) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(paste0(opt$SampleName," A/B comp"),"  \nNumber of DNase peaks in contact\nwith insertion (5kb resolution)")+
#       xlab('Number of DNase peaks in contact (HiC)')+
#       ylab('log10(RNA/DNA)')+
#       coord_cartesian(ylim = c(-3, 3))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/DNaseContacts_ABcomp.pdf"), width=6, height=4)
#ggDnaseContactsAB
#dev.off()
#
##DnaseContactAll <- read.table(paste0(opt$OUTDIR, "/DnaseContact_all.txt"),sep=' ', stringsAsFactors=F)
#
#ggDnaseContactsAB_All <- normdataiPCR %>%
#       select(BC, DNA, RNA, SampleName, nDHS100kb, AB_compartment_value) %>%
#       mutate(ABcomp = ifelse(AB_compartment_value > 0, 'A', ifelse(AB_compartment_value < 0, 'B', 'NA'))) %>%
#       filter(!is.na(ABcomp), ABcomp!='NA')  %>%
#       mutate(bins = cut.pretty(nDHS100kb,breaks=c(0, 1, 2, 5, 10, 20,100), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
#       mutate(Expression = log10(RNA/DNA)) %>%
#       ggplot(aes(x=bins, y=Expression, fill=ABcomp)) +
#       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=1, colour=sampleColor) +
#       geom_boxplot(outlier.shape=NA,alpha=0.6)+ 
#       theme_classic()+
#       facet_wrap(~ABcomp) +
#       scale_fill_brewer(palette='Set1') +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(paste0(opt$SampleName," A/B comp"),"  \nNumber of DNase peaks surrounding\n +-100kb of insertion")+
#       xlab('Number of DNase peaks (+-100kb)')+
#       ylab('log10(RNA/DNA)')+
#       coord_cartesian(ylim = c(-3, 3))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/DNaseContacts_ABcomp_all.pdf"), width=6, height=4)
#ggDnaseContactsAB_All
#dev.off()
#
#
######
##Contacts filtered by hingest DHS density or strongest contact
######
#ggMaxDHSinContact <- normdataiPCR %>%
#       mutate(maxDHSdensity=replace(maxDHSdensity, is.na(maxDHSdensity), 0)) %>%
#       #filter(!is.na(maxDHSdensity)) %>%
#       ggplot(aes(x=maxDHSdensity, y=RNA/DNA, fill=SampleName)) +
#       geom_point( pch=21, colour=sampleColor) +
#       #geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       geom_smooth(method = "lm") +
#       theme_classic()+
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName, 'Max DHS density in contact with insertion')+
#       xlab('Max DHS density in contact')+
#       ylab('log10(RNA/DNA)')+
#       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
#       scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
#       annotation_logticks(sides = "lb",long = unit(0.2, "cm"))
#
#

######
##How number of DNase peaks in contact with inertions affects expression 
######
##DnaseContact <- read.table(paste0(opt$OUTDIR, "/DnaseContact.txt"),sep=' ', stringsAsFactors=F)
#
##DnaseContact <- read.table('/tmp/CMV_dnaseContact.txt',sep=' ', stringsAsFactors=F)
#cat('DNase contacts (counts all DHS in windows)','\n')
#ggDnaseContacts <- normdataiPCR %>%
#       select(BC, DNA, RNA, SampleName, contactDNase) %>%
#       mutate(bins = cut.pretty(contactDNase,breaks=c(0, 1, 2, 5, 10, 20,100), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
#       mutate(Expression = log10(RNA/DNA)) %>%
#       ggplot(aes(x=bins, y=Expression, fill=SampleName)) +
#       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=1, colour=sampleColor) +
#       geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       theme_classic()+
#       #scale_fill_manual(values=sampleColor) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName,"  \nNumber of DNase peaks in contact\nwith insertion (5kb resolution)")+
#       xlab('Number of DNase peaks in contact (HiC)')+
#       ylab('log10(RNA/DNA)')+
#       coord_cartesian(ylim = c(-3, 3))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/DNaseContacts.pdf"), width=5, height=4)
#ggDnaseContacts
#dev.off()
#
#
#cat('DNase contacts scatterplot (counts all DHS in windows)','\n')
#ggDnaseContactsScatter <- normdataiPCR %>%
#       select(BC, DNA, RNA, SampleName, contactDNase) %>%
#       #mutate(bins = cut.pretty(contactDNase,breaks=c(0, 1, 2, 5, 10, 20,100), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
#       mutate(Expression = log10(RNA/DNA)) %>%
#       ggplot(aes(x=contactDNase, y=Expression, fill=SampleName)) +
#       geom_point( pch=21, colour=sampleColor) +
#       #geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       geom_smooth(method = "loess") +
#       theme_classic()+
#       #scale_fill_manual(values=sampleColor) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName,"  \nNumber of DNase peaks in contact\nwith insertion (5kb resolution)")+
#       xlab('Number of K562 DHS in contact')+
#       ylab('log10(RNA/DNA)')+
#       coord_cartesian(ylim = c(-3, 3))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/DNaseContacts_scatter.pdf"), width=4, height=4)
#ggDnaseContactsScatter
#dev.off()
#
#
##DnaseContact <- read.table('/tmp/CMV_dnaseContact.txt',sep=' ', stringsAsFactors=F)
#cat('DNase contacts (only look if theres one DHS peak in each window) ','\n')
#ggDnaseContacts <- normdataiPCR %>%
#       select(BC, DNA, RNA, SampleName, contactDNase_window) %>%
#       mutate(bins = cut.pretty(contactDNase_window,breaks=c(0, 1, 2, 5, 10, 20,100), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
#       mutate(Expression = log10(RNA/DNA)) %>%
#       ggplot(aes(x=bins, y=Expression, fill=SampleName)) +
#       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=1, colour=sampleColor) +
#       geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       theme_classic()+
#       #scale_fill_manual(values=sampleColor) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName,"  \nNumber of windows with a DNase peaks in contact\nwith insertion (5kb resolution)")+
#       xlab('Number of windows cointaining a DHS in contact')+
#       ylab('log10(RNA/DNA)')+
#       coord_cartesian(ylim = c(-3, 3))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/DNaseContacts_windows.pdf"), width=5, height=4)
#ggDnaseContacts
#dev.off()
#
#
#cat('DNase contacts scatterplot (only look if theres one DHS peak in each window)','\n')
#ggDnaseContactsScatter <- normdataiPCR %>%
#       select(BC, DNA, RNA, SampleName, contactDNase_window) %>%
#       #mutate(bins = cut.pretty(contactDNase,breaks=c(0, 1, 2, 5, 10, 20,100), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
#       mutate(Expression = log10(RNA/DNA)) %>%
#       ggplot(aes(x=contactDNase_window, y=Expression, fill=SampleName)) +
#       geom_point( pch=21, colour=sampleColor) +
#       #geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       geom_smooth(method = "lm") +
#       theme_classic()+
#       #scale_fill_manual(values=sampleColor) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName,"  \nNumber of windows with a DNase peaks in contact\nwith insertion (5kb resolution)")+
#       xlab('Number of windows cointaining a DHS in contact')+
#       ylab('log10(RNA/DNA)')+
#       coord_cartesian(ylim = c(-3, 3))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/DNaseContacts_scatter_windows.pdf"), width=5, height=4)
#ggDnaseContactsScatter
#dev.off()
#
#
#cat('All contacts','\n')
#ggDnaseContacts <- normdataiPCR %>%
#       select(BC, DNA, RNA, SampleName, ObsExpover2) %>%
#       mutate(bins = cut.pretty(ObsExpover2,breaks=c(0, 1, 2, 5, 10, 20,100), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
#       mutate(Expression = log10(RNA/DNA)) %>%
#       ggplot(aes(x=bins, y=Expression, fill=SampleName)) +
#       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=1, colour=sampleColor) +
#       geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       theme_classic()+
#       #scale_fill_manual(values=sampleColor) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName,"  \nNumber of 5kb windows in contact with insertion")+
#       xlab('Number of 5kb windows in contact')+
#       ylab('log10(RNA/DNA)')+
#       coord_cartesian(ylim = c(-3, 3))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/Contacts_all.pdf"), width=5, height=4)
#ggDnaseContacts
#dev.off()
#
#cat('All contacts scatterplot','\n')
#ggContactsScatter <- normdataiPCR %>%
#       select(BC, DNA, RNA, SampleName, ObsExpover2) %>%
#       #mutate(bins = cut.pretty(contactDNase,breaks=c(0, 1, 2, 5, 10, 20,100), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
#       mutate(Expression = log10(RNA/DNA)) %>%
#       ggplot(aes(x=ObsExpover2, y=Expression, fill=SampleName)) +
#       geom_point( pch=21, colour=sampleColor) +
#       #geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       geom_smooth(method = "loess") +
#       theme_classic()+
#       #scale_fill_manual(values=sampleColor) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName,"  \nNumber of 5kb windows in contact with insertion")+
#       xlab('Number of 5kb windows in contact')+
#       ylab('log10(RNA/DNA)')+
#       coord_cartesian(ylim = c(-3, 3))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/Contacts_all_scatter.pdf"), width=5, height=4)
#ggContactsScatter
#dev.off()
#
#
##DnaseContact <- read.table('/tmp/CMV_dnaseContact.txt',sep=' ', stringsAsFactors=F)
#cat('HeLa DNase contacts','\n')
#ggDnaseContactsHeLa <- normdataiPCR %>%
#       select(BC, DNA, RNA, SampleName, HeLacontactDNase) %>%
#       mutate(bins = cut.pretty(HeLacontactDNase,breaks=c(0, 1, 2, 5, 10, 20,100), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
#       mutate(Expression = log10(RNA/DNA)) %>%
#       ggplot(aes(x=bins, y=Expression, fill=SampleName)) +
#       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=1, colour=sampleColor) +
#       geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       theme_classic()+
#       #scale_fill_manual(values=sampleColor) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName,"  \nNumber of HeLa DNase peaks in contact\nwith insertion (5kb resolution)")+
#       xlab('Number of HeLa DHS in contact')+
#       ylab('log10(RNA/DNA)')+
#       coord_cartesian(ylim = c(-3, 3))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/HeLa_DNaseContacts.pdf"), width=5, height=4)
#ggDnaseContactsHeLa
#dev.off()
#
#
#
#cat('HeLa DNase contacts scatterplot','\n')
#ggDnaseContactsScatterHeLa <- normdataiPCR %>%
#       select(BC, DNA, RNA, SampleName, HeLacontactDNase) %>%
#       #mutate(bins = cut.pretty(contactDNase,breaks=c(0, 1, 2, 5, 10, 20,100), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
#       mutate(Expression = log10(RNA/DNA)) %>%
#       ggplot(aes(x=HeLacontactDNase, y=Expression, fill=SampleName)) +
#       geom_point( pch=21, colour=sampleColor) +
#       #geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       geom_smooth(method = "loess") +
#       theme_classic()+
#       #scale_fill_manual(values=sampleColor) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName,"  \nNumber of DNase peaks in contact\nwith insertion (5kb resolution)")+
#       xlab('Number of HeLa DHS in contact')+
#       ylab('log10(RNA/DNA)')+
#       coord_cartesian(ylim = c(-3, 3))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/HeLa_DNaseContacts_scatter.pdf"), width=5, height=4)
#ggDnaseContactsScatterHeLa
#dev.off()
#
#
#cat('HeLa DNase contacts windows','\n')
#ggDnaseContactsHeLa <- normdataiPCR %>%
#       select(BC, DNA, RNA, SampleName, HeLacontactDNase_window) %>%
#       mutate(bins = cut.pretty(HeLacontactDNase_window,breaks=c(0, 1, 2, 5, 10, 20,100), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
#       mutate(Expression = log10(RNA/DNA)) %>%
#       ggplot(aes(x=bins, y=Expression, fill=SampleName)) +
#       geom_point(position=position_jitterdodge(dodge.width=0.75),alpha=0.2, pch=1, colour=sampleColor) +
#       geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       theme_classic()+
#       #scale_fill_manual(values=sampleColor) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName,"  \nNumber of HeLa DNase peaks in contact\nwith insertion (5kb resolution)")+
#       xlab('Number of HeLa windows cointaining a DHS in contact')+
#       ylab('log10(RNA/DNA)')+
#       coord_cartesian(ylim = c(-3, 3))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/HeLa_DNaseContacts_window.pdf"), width=5, height=4)
#ggDnaseContactsHeLa
#dev.off()
#
#
#
#cat('HeLa DNase contacts scatterplot windows','\n')
#ggDnaseContactsScatterHeLa <- normdataiPCR %>%
#       select(BC, DNA, RNA, SampleName, HeLacontactDNase_window) %>%
#       #mutate(bins = cut.pretty(contactDNase,breaks=c(0, 1, 2, 5, 10, 20,100), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
#       mutate(Expression = log10(RNA/DNA)) %>%
#       ggplot(aes(x=HeLacontactDNase_window, y=Expression, fill=SampleName)) +
#       geom_point( pch=21, colour=sampleColor) +
#       #geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       geom_smooth(method = "lm") +
#       theme_classic()+
#       #scale_fill_manual(values=sampleColor) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName,"  \nNumber of HeLa DNase peaks in contact\nwith insertion (5kb resolution)")+
#       xlab('Number of HeLa windows cointaining a DHS in contact')+
#       ylab('log10(RNA/DNA)')+
#       coord_cartesian(ylim = c(-3, 3))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/HeLa_DNaseContacts_scatter_window.pdf"), width=5, height=4)
#ggDnaseContactsScatterHeLa
#dev.off()
#
#
#cat('HeLa DNase contacts scatterplot windows without windows containing K562 DHS','\n')
#ggDnaseContactsScatterHeLa <- normdataiPCR %>%
#       select(BC, DNA, RNA, SampleName, HeLacontactDNase_window_noK562) %>%
#       #mutate(bins = cut.pretty(contactDNase,breaks=c(0, 1, 2, 5, 10, 20,100), include.lowest=T, right=F, include.na.separate=T, pretty.labels="left")) %>%
#       mutate(Expression = log10(RNA/DNA)) %>%
#       ggplot(aes(x=HeLacontactDNase_window_noK562, y=Expression, fill=SampleName)) +
#       geom_point( pch=21, colour=sampleColor) +
#       #geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       geom_smooth(method = "lm") +
#       theme_classic()+
#       #scale_fill_manual(values=sampleColor) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName,"  \nNumber of HeLa DNase peaks in contact\nwith insertion (5kb resolution)\nminus K562 DHS containing windows")+
#       xlab('Number of HeLa windows cointaining a DHS in contact')+
#       ylab('log10(RNA/DNA)')+
#       coord_cartesian(ylim = c(-3, 3))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm")) +
#       theme(legend.position="none")
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/HeLa_DNaseContacts_scatter_window_noK562.pdf"), width=5, height=4)
#ggDnaseContactsScatterHeLa
#dev.off()
#
#
#
#ggProportionDHScontacts <- normdataiPCR %>%
#       ggplot(aes(x=contactDNase_window/ObsExpover2, fill=SampleName)) +
#       geom_histogram(binwidth=0.05, fill=sampleColor, colour='black', alpha=0.6) +
#       theme_classic()+
#       #scale_fill_manual(values=sampleColor) +
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName,"  \nProportion of 5kb windows +-300kb around insertions\nthat are in contact with DHS")+
#       xlab('Proportion of windows in contact with DHS')+
#       theme(legend.position="none")
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/Contacts_proportionDHS.pdf"), width=5, height=4)
#ggProportionDHScontacts
#dev.off()
#
#
#
#maxContactsDHS <- normdataiPCR %>%
#       select(chrom, chromStart, chromEnd, BC, maxContactContact) %>%
#       filter(!is.na(maxContactContact)) %>% 
#       arrange(desc(maxContactContact)) %>%
#       head(100) %>%
#       mutate(chromStart = chromStart-3e5,
#              chromEnd= chromEnd+3e5)
#       
#write.table(maxContactsDHS, file='~/public_html/blog/2017Aug21/maxContactsDHS.bed', sep='\t', quote=F, col.names=F, row.names=F)
#
#ggMaxContactWithDHS <- normdataiPCR %>%
#       mutate(maxContactDHSdensity=replace(maxContactDHSdensity, is.na(maxContactDHSdensity), 0)) %>%
#       #filter(!is.na(maxContactDHSdensity)) %>%
#       ggplot(aes(x=maxContactDHSdensity, y=RNA/DNA, fill=SampleName)) +
#       geom_point( pch=21, colour=sampleColor) +
#       #geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       geom_smooth(method = "loess") +
#       theme_classic()+
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName, 'DHS density of strongest contact')+
#       xlab('DHS density')+
#       ylab('log10(RNA/DNA)')+
#       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
#       scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
#       annotation_logticks(sides = "lb",long = unit(0.2, "cm"))
#
#ggContactMCV <- normdataiPCR %>%
#       filter(!is.na(maxContactMCV)) %>%
#       ggplot(aes(x=maxContactMCV, y=RNA/DNA, fill=SampleName)) +
#       geom_point( pch=21, colour=sampleColor) +
#       #geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       geom_smooth(method = "loess") +
#       theme_classic()+
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName, 'MCV of strongest connected DHS')+
#       xlab('MCV')+
#       ylab('log10(RNA/DNA)')+
#       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm"))
#       
#pdf(file=paste0(opt$OUTDIR,"/graphs/contactsMaxContactMCV.pdf"), width=5, height=4)
#ggContactMCV
#dev.off()
#
#
#ggContactDHSsum <- normdataiPCR %>%
#       mutate(maxContactDHSsum=replace(maxContactDHSsum, is.na(maxContactDHSsum), 0)) %>%
#       filter(!is.na(maxContactDHSsum)) %>%
#       ggplot(aes(x=log10(maxContactDHSsum), y=RNA/DNA, fill=SampleName)) +
#       geom_point( pch=21, colour=sampleColor) +
#       #geom_boxplot(outlier.shape=NA,fill=sampleColor,alpha=0.6)+ 
#       geom_smooth(method = "loess") +
#       theme_classic()+
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName, 'MCV of strongest connected DHS')+
#       xlab('MCV')+
#       ylab('log10(RNA/DNA)')+
#       scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm"))
#
#

####
#3d genome
####
#contactDNase_window
#maxDHScontact
#maxDHSdensity
#maxContactContact
#maxContactDHSdensity
#contactDNase
##TAD_armatus0.5_DHS 
#
#
#threeDGenome <- data.frame(variable=c('contactDNase_window',  'maxDHSdensity', 'contactDNase', 'TAD_armatus0.5_DHS', 'maxContactDHSdensity'),
#       pval=c(signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$contactDNase_window)))$coef[2,4],2), 
#              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$maxDHSdensity)))$coef[2,4],2),
#              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$contactDNase)))$coef[2,4],2),
#              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$TAD_armatus0.5_DHS)))$coef[2,4],2),
#              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$maxContactDHSdensity)))$coef[2,4],2)),
#       slope=c(signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$contactDNase_window)))$coef[[2]],2), 
#              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$maxDHSdensity)))$coef[[2]],2),
#              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$contactDNase)))$coef[[2]],2),
#              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$TAD_armatus0.5_DHS)))$coef[[2]],2),
#              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$maxContactDHSdensity)))$coef[[2]],2)),
#       adjR2=c(signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$contactDNase_window)))$adj.r.squared,2), 
#              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$maxDHSdensity)))$adj.r.squared,2),
#              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$contactDNase)))$adj.r.squared,2),
#              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$TAD_armatus0.5_DHS)))$adj.r.squared,2),
#              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$maxContactDHSdensity)))$adj.r.squared,2)),
#       Intercept=c(signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$contactDNase_window)))$coef[[1]],2), 
#              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$maxDHSdensity)))$coef[[1]],2),
#       signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$contactDNase)))$coef[[1]],2),
#              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$TAD_armatus0.5_DHS)))$coef[[1]],2),
#              signif(summary(lm(log10(normdataiPCR$RNA/normdataiPCR$DNA)~abs(normdataiPCR$maxContactDHSdensity)))$coef[[1]],2)),
#              x=c(10, 1, 25, 100, 1))
#
#
#ggthreeDGenome <- normdataiPCR %>% 
#              mutate(expression = log10(RNA/DNA),
#                     maxDHScontact = log10(maxDHScontact),
#                     maxDHSdensity = log10(maxDHSdensity),
#                     maxContactContact = log10(maxContactContact),
#                     maxContactDHSdensity = log10(maxContactDHSdensity), 
#                     maxDHScontact=replace(maxDHScontact, is.na(maxDHScontact), 0),
#                     maxDHSdensity=replace(maxDHSdensity, is.na(maxDHSdensity), 0),
#                     maxContactDHSdensity=replace(maxContactDHSdensity, is.na(maxContactDHSdensity), 0),
#                     contactDNase=replace(contactDNase, is.na(contactDNase), 0),
#                     TAD_armatus0.5_DHS=replace(TAD_armatus0.5_DHS, is.na(TAD_armatus0.5_DHS), 0)) %>%
#              filter(!is.na(expression)) %>%
#              select(expression, contactDNase_window, maxDHSdensity, maxContactDHSdensity, contactDNase, TAD_armatus0.5_DHS )%>%
#              melt(id='expression') %>%
#       ggplot(aes(x=value, y=expression, colour=variable)) +
#       geom_point( pch=21, colour=sampleColor) +
#       facet_wrap(~variable, scale='free_x') +
#       geom_smooth(method = "lm") +
#       geom_smooth(method = "loess", linetype='dashed') +
#       geom_text(data=threeDGenome, aes(y=2, x=x, label=paste0('adjR2=',adjR2)), colour='black', size=5)+
#       geom_text(data=threeDGenome, aes(y=1.7, x=x, label=paste0('slope=',slope)), colour='black', size=5)+
#       geom_text(data=threeDGenome, aes(y=1.4, x=x, label=paste0('intercept=',Intercept)), colour='black', size=5)+
#       theme_classic()+
#       theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#       ggtitle(opt$SampleName, '3d genome variables')+
#       ylab('log10(RNA/DNA)')+
#       annotation_logticks(sides = "l",long = unit(0.2, "cm"))
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/threeDGenome.pdf"), width=14, height=8)
#ggthreeDGenome
#dev.off()
#
#
#
#
#
#linearGenomeLM <- normdataiPCR %>%
#                     mutate(expression = log10(RNA/DNA),
#                     DistToNearestDHS = log10(abs(DistToNearestDHS)+1),
#                     TAD_armatus0.5_length = log10(TAD_armatus0.5_length)) %>%
#                     mutate(expression = ifelse(log10(RNA/DNA) >as.numeric(sumExpression[5]), 1, ifelse(log10(RNA/DNA) <(as.numeric(sumExpression[2])), 0,NA))) %>%
#                     filter(!is.na(expression)) %>%
#                     select(expression, DistToNearestDHS, nDHS10kb, nDHS100kb, nDHS300kb, DNasePksDens)
#                     
#Expression <- linearGenomeLM$expression
#
#linearGenomeLM <- subset(linearGenomeLM, select=-c(expression))
#fitlin <- glm(formula=Expression ~. ,data=linearGenomeLM, family='binomial')
#
#f_plot_ROC(fitlin$y,fitlin$fitted.values, main="")
#
#
#
#threeDGenomeLM <- normdataiPCR %>% 
#              mutate(expression = log10(RNA/DNA),
#                     maxDHScontact = log10(maxDHScontact),
#                     maxDHSdensity = log10(maxDHSdensity),
#                     maxContactContact = log10(maxContactContact),
#                     maxContactDHSdensity = log10(maxContactDHSdensity), 
#                     maxDHScontact=replace(maxDHScontact, is.na(maxDHScontact), 0),
#                     maxDHSdensity=replace(maxDHSdensity, is.na(maxDHSdensity), 0),
#                     maxContactDHSdensity=replace(maxContactDHSdensity, is.na(maxContactDHSdensity), 0),
#                     contactDNase=replace(contactDNase, is.na(contactDNase), 0),
#                     TAD_armatus0.5_DHS=replace(TAD_armatus0.5_DHS, is.na(TAD_armatus0.5_DHS), 0)) %>%
#              mutate(expression = ifelse(log10(RNA/DNA) >as.numeric(sumExpression[5]), 1, ifelse(log10(RNA/DNA) <(as.numeric(sumExpression[2])), 0,NA))) %>%
#              filter(!is.na(expression)) %>%
#              select(expression, contactDNase_window, maxDHSdensity, maxContactDHSdensity, contactDNase, TAD_armatus0.5_DHS )
#              
#Expression3D <- threeDGenomeLM$expression
#threeDGenomeLM <- subset(threeDGenomeLM, select=-c(expression))
#fit3d <- glm(formula=Expression ~. ,data=threeDGenomeLM, family='binomial')
#
#f_plot_ROC(fit3d$y,fit3d$fitted.values, main="")
#
#
#fitboth <- glm(formula=Expression ~. ,data=cbind(threeDGenomeLM,linearGenomeLM), family='binomial')
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/Linear_vs_3D.pdf"), width=12, height=4)
#par(mfrow=c(1,3))
#f_plot_ROC(fitlin$y,fitlin$fitted.values, main="Linear features")
#f_plot_ROC(fit3d$y,fit3d$fitted.values, main="3D features")
#f_plot_ROC(fitboth$y,fitboth$fitted.values, main="Combined")
#dev.off()
##

#ObsExp.filter <- 2
#percentageOfMax.filter <- 1
#DHSdensity.filter <- 1
#
#contacts <- contacts[contacts$BC %in% normdataiPCR$BC,]
#allDHS <- contacts %>%
#       group_by(BC) %>%
#       filter(percentageOfMax >= percentageOfMax.filter,
#       ObsExpPercentage >= ObsExp.filter ,
#       DHSdensity >=DHSdensity.filter,
#       distance<=1e6 & distance>=-1e6) %>%
#       mutate(DHSdensity = ifelse(DHSdensity > 0, 1, ifelse(DHSdensity == 0, 0, 0))) %>%
#       select(BC, DHSdensity) %>% 
#       summarize(n=sum(DHSdensity)) %>%
#       rename('contactsDHS_1Mb'=n)
#       
#
#allDHS <- contacts %>%
#       group_by(BC) %>%
#       filter(percentageOfMax >= percentageOfMax.filter,
#       ObsExpPercentage >= ObsExp.filter ,
#       DHSdensity >=DHSdensity.filter,
#       distance<=5e5 & distance>=-5e5) %>%
#       mutate(DHSdensity = ifelse(DHSdensity > 0, 1, ifelse(DHSdensity == 0, 0, 0))) %>%
#       select(BC, DHSdensity) %>%
#       summarize(n=sum(DHSdensity)) %>%
#       rename('contactsDHS_500kb'=n) %>%
#       left_join(allDHS, by='BC')
#       
#       
#allDHS <- contacts %>%
#       group_by(BC) %>%
#       filter(percentageOfMax >= percentageOfMax.filter,
#       ObsExpPercentage >= ObsExp.filter ,
#       DHSdensity >=DHSdensity.filter,
#       distance<=1e5 & distance>=-1e5) %>%
#       mutate(DHSdensity = ifelse(DHSdensity > 0, 1, ifelse(DHSdensity == 0, 0, 0))) %>%
#       select(BC, DHSdensity) %>%
#       summarize(n=sum(DHSdensity)) %>%
#       rename('contactsDHS_100kb'=n) %>%
#       left_join(allDHS, by='BC')
#       
#       
#allDHS <- contacts %>%
#       group_by(BC) %>%
#       filter(percentageOfMax >= percentageOfMax.filter,
#       ObsExpPercentage >= ObsExp.filter ,
#       DHSdensity >=DHSdensity.filter,
#       distance<=1e6 & distance>=-1e6) %>%
#       mutate(CTCFdensity = ifelse(CTCFdensity > 0, 1, ifelse(CTCFdensity == 0, 0, 0))) %>%
#       select(BC, CTCFdensity) %>%
#       summarize(n=sum(CTCFdensity)) %>%
#       rename('contactsCTCF_1Mb'=n) %>%
#       left_join(allDHS, by='BC')
#
#
#allDHS <- contacts %>%
#       group_by(BC) %>%
#       filter(percentageOfMax >= percentageOfMax.filter,
#       ObsExpPercentage >= ObsExp.filter ,
#       DHSdensity >=DHSdensity.filter,
#       distance<=5e5 & distance>=-5e5) %>%
#       mutate(CTCFdensity = ifelse(CTCFdensity > 0, 1, ifelse(CTCFdensity == 0, 0, 0))) %>%
#       select(BC, CTCFdensity) %>%
#       summarize(n=sum(CTCFdensity)) %>%
#       rename('contactsCTCF_500kb'=n) %>%
#       left_join(allDHS, by='BC')
#       
#allDHS <- contacts %>%
#       group_by(BC) %>%
#       filter(percentageOfMax >= percentageOfMax.filter,
#       ObsExpPercentage >= ObsExp.filter ,
#       DHSdensity >=DHSdensity.filter,
#       distance<=1e5 & distance>=-1e5) %>%
#       mutate(CTCFdensity = ifelse(CTCFdensity > 0, 1, ifelse(CTCFdensity == 0, 0, 0))) %>%
#       select(BC, CTCFdensity) %>%
#       summarize(n=sum(CTCFdensity)) %>%
#       rename('contactsCTCF_100kb'=n) %>%
#       left_join(allDHS, by='BC')
#
##
#allDHS <- contacts %>%
#       group_by(BC) %>%
#       filter(percentageOfMax >1,
#       ObsExpPercentage >= ObsExp.filter ,
#       distance<=1e6 & distance>= -1e6) %>%
#       select(BC, DHSdensity) %>%
#      # mutate(DHSdensity = ifelse(DHSdensity > 0, 1, ifelse(DHSdensity == 0, 0, 0))) %>%
#       summarize(n=sum(DHSdensity)) %>%
#       rename('contactsDHSdensity_1Mb'=n) %>%
#       left_join(allDHS, by='BC')
#
#allDHS <- contacts %>%
#       group_by(BC) %>%
#       filter(percentageOfMax >1,
#       ObsExpPercentage >= ObsExp.filter ,
#       distance<=5e5 & distance>=-5e5) %>%
#       select(BC, DHSdensity) %>%
#       summarize(n=sum(DHSdensity)) %>%
#       rename('contactsDHSdensity_500kb'=n) %>%
#       left_join(allDHS, by='BC')
#
#
#allDHS <- contacts %>%
#       group_by(BC) %>%
#       filter(percentageOfMax >1,
#       ObsExpPercentage >= ObsExp.filter ,
#       distance<=1e5 & distance>=-1e5) %>%
#       select(BC, DHSdensity) %>%
#       summarize(n=sum(DHSdensity)) %>%
#       rename('contactsDHSdensity_100kb'=n) %>%
#       left_join(allDHS, by='BC')
#
#
#
#
#
#
#allDHS <- contacts %>%
#       group_by(BC) %>%
#       filter(distance<=1e6 & distance>=-1e6) %>%
#       mutate(DHSdensity = ifelse(DHSdensity > 0, 1, ifelse(DHSdensity == 0, 0, 0))) %>%
#       select(BC, DHSdensity) %>% 
#       summarize(n=sum(DHSdensity)) %>%
#       rename('allDHS_1Mb'=n) %>%
#       left_join(allDHS, by='BC')
#       
#
#allDHS <- contacts %>%
#       group_by(BC) %>%
#       filter(distance<=5e5 & distance>=-5e5) %>%
#       mutate(DHSdensity = ifelse(DHSdensity > 0, 1, ifelse(DHSdensity == 0, 0, 0))) %>%
#       select(BC, DHSdensity) %>%
#       summarize(n=sum(DHSdensity)) %>%
#       rename('allDHS_500kb'=n) %>%
#       left_join(allDHS, by='BC')
#       
#       
#allDHS <- contacts %>%
#       group_by(BC) %>%
#       filter(distance<=1e5 & distance>=-1e5) %>%
#       mutate(DHSdensity = ifelse(DHSdensity > 0, 1, ifelse(DHSdensity == 0, 0, 0))) %>%
#       select(BC, DHSdensity) %>%
#       summarize(n=sum(DHSdensity)) %>%
#       rename('allDHS_100kb'=n) %>%
#       left_join(allDHS, by='BC')
#       
#
#
#allDHS <- contacts %>%
#       group_by(BC) %>%
#       filter(distance<=1e6 & distance>= -1e6) %>%
#       select(BC, DHSdensity) %>%
#      # mutate(DHSdensity = ifelse(DHSdensity > 0, 1, ifelse(DHSdensity == 0, 0, 0))) %>%
#       summarize(n=sum(DHSdensity)) %>%
#       rename('allDHSdensity_1Mb'=n) %>%
#       left_join(allDHS, by='BC')
#
#allDHS <- contacts %>%
#       group_by(BC) %>%
#       filter(distance<=5e5 & distance>=-5e5) %>%
#       select(BC, DHSdensity) %>%
#       summarize(n=sum(DHSdensity)) %>%
#       rename('allDHSdensity_500kb'=n) %>%
#       left_join(allDHS, by='BC')
#
#
#allDHS <- contacts %>%
#       group_by(BC) %>%
#       filter(distance<=1e5 & distance>=-1e5) %>%
#       select(BC, DHSdensity) %>%
#       summarize(n=sum(DHSdensity)) %>%
#       rename('allDHSdensity_100kb'=n) %>%
#       left_join(allDHS, by='BC')
#
#
#
#threeDGenome <- normdataiPCR %>%
#       left_join(allDHS, by='BC') 
#
#threeD <- data.frame(variable=c('allDHS_1Mb', 'allDHS_500kb', 'allDHS_100kb', 'contactsDHS_1Mb', 'contactsDHS_500kb', 'contactsDHS_100kb'),
#       pval=c(signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$allDHS_1Mb)))$coef[2,4],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$allDHS_500kb)))$coef[2,4],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$allDHS_100kb)))$coef[2,4],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$contactsDHS_1Mb)))$coef[2,4],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$contactsDHS_500kb)))$coef[2,4],2),  
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$contactsDHS_100kb)))$coef[2,4],2)),
#       slope=c(signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$allDHS_1Mb)))$coef[[2]],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$allDHS_500kb)))$coef[[2]],2), 
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$allDHS_100kb)))$coef[[2]],2), 
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$contactsDHS_1Mb)))$coef[[2]],2), 
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$contactsDHS_500kb)))$coef[[2]],2),  
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$contactsDHS_100kb)))$coef[[2]],2)),
#       adjR2=c(signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$allDHS_1Mb)))$adj.r.squared,2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$allDHS_500kb)))$adj.r.squared,2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$allDHS_100kb)))$adj.r.squared,2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$contactsDHS_1Mb)))$adj.r.squared,2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$contactsDHS_500kb)))$adj.r.squared,2), 
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$contactsDHS_100kb)))$adj.r.squared,2)),
#       Intercept=c(signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$allDHS_1Mb)))$coef[[1]],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$allDHS_500kb)))$coef[[1]],2), 
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$allDHS_100kb)))$coef[[1]],2),  
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$contactsDHS_1Mb)))$coef[[1]],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$contactsDHS_500kb)))$coef[[1]],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~abs(threeDGenome$contactsDHS_100kb)))$coef[[1]],2)),
#       x=c(100,50,15,50,50,10))
#       
#       
#       
#ggLinearVs3D <- normdataiPCR %>%
#       left_join(allDHS, by='BC') %>%
#       mutate(expression=log10(RNA/DNA)) %>%
#       select(expression, allDHS_1Mb, allDHS_500kb, allDHS_100kb, contactsDHS_1Mb, contactsDHS_500kb, contactsDHS_100kb) %>%
#       melt(id='expression') %>%
#       ggplot(aes(x=value, y=expression, colour=variable)) +
#              geom_point( pch=21, colour=sampleColor) +
#              facet_wrap(~variable, scale='free_x') +
#              geom_smooth(method = "lm") +
#              geom_smooth(method = "loess", linetype='dashed') +
#              geom_text(data=threeD, aes(y=2, x=x, label=paste0('adjR2=',adjR2)), colour='black', size=5)+
#              geom_text(data=threeD, aes(y=1.7, x=x, label=paste0('slope=',slope)), colour='black', size=5)+
#              geom_text(data=threeD, aes(y=1.4, x=x, label=paste0('intercept=',Intercept)), colour='black', size=5)+
#              theme_classic()+
#              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#              ggtitle(opt$SampleName, 'number of DHS surrounding the insertion (top)\nnumber of 5kb windows with a DHS in contact with the insertion (bottom)')+
#              ylab('log10(RNA/DNA)')+
#              annotation_logticks(sides = "l",long = unit(0.2, "cm"))
#       
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/LinearVs3D_scatter.pdf"), width=14, height=8)
#ggLinearVs3D
#dev.off()
#
#
#
#
#threeDensity <- data.frame(variable=c('allDHSdensity_1Mb', 'allDHSdensity_500kb', 'allDHSdensity_100kb', 'contactsDHSdensity_1Mb', 'contactsDHSdensity_500kb', 'contactsDHSdensity_100kb'),
#       pval=c(signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$allDHSdensity_1Mb)+1)))$coef[2,4],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$allDHSdensity_500kb)+1)))$coef[2,4],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$allDHSdensity_100kb)+1)))$coef[2,4],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$contactsDHSdensity_1Mb)+1)))$coef[2,4],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$contactsDHSdensity_500kb)+1)))$coef[2,4],2),  
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$contactsDHSdensity_100kb)+1)))$coef[2,4],2)),
#       slope=c(signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$allDHSdensity_1Mb)+1)))$coef[[2]],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$allDHSdensity_500kb)+1)))$coef[[2]],2), 
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$allDHSdensity_100kb)+1)))$coef[[2]],2), 
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$contactsDHSdensity_1Mb)+1)))$coef[[2]],2), 
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$contactsDHSdensity_500kb)+1)))$coef[[2]],2),  
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$contactsDHSdensity_100kb)+1)))$coef[[2]],2)),
#       adjR2=c(signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$allDHSdensity_1Mb)+1)))$adj.r.squared,2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$allDHSdensity_500kb)+1)))$adj.r.squared,2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$allDHSdensity_100kb)+1)))$adj.r.squared,2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$contactsDHSdensity_1Mb)+1)))$adj.r.squared,2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$contactsDHSdensity_500kb)+1)))$adj.r.squared,2), 
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$contactsDHSdensity_100kb)+1)))$adj.r.squared,2)),
#       Intercept=c(signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$allDHSdensity_1Mb)+1)))$coef[[1]],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$allDHSdensity_500kb)+1)))$coef[[1]],2), 
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$allDHSdensity_100kb)+1)))$coef[[1]],2),  
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$contactsDHSdensity_1Mb)+1)))$coef[[1]],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$contactsDHSdensity_500kb)+1)))$coef[[1]],2),
#              signif(summary(lm(log10(threeDGenome$RNA/threeDGenome$DNA)~log10(abs(threeDGenome$contactsDHSdensity_100kb)+1)))$coef[[1]],2)),
#       x=c(2,2,2,2,2,2))
#       
#       
#       
#ggLinearVs3DDensity <- normdataiPCR %>%
#       left_join(allDHS, by='BC') %>%
#       mutate(expression=log10(RNA/DNA)) %>%
#       select(expression, allDHSdensity_1Mb, allDHSdensity_500kb, allDHSdensity_100kb, contactsDHSdensity_1Mb, contactsDHSdensity_500kb, contactsDHSdensity_100kb) %>%
#       melt(id='expression') %>%
#       mutate(value=log10(value+1)) %>%
#       ggplot(aes(x=value, y=expression, colour=variable)) +
#              geom_point( pch=21, colour=sampleColor) +
#              facet_wrap(~variable, scale='free_x') +
#              geom_smooth(method = "lm") +
#              geom_smooth(method = "loess", linetype='dashed') +
#              geom_text(data=threeDensity, aes(y=2, x=x, label=paste0('adjR2=',adjR2)), colour='black', size=5)+
#              geom_text(data=threeDensity, aes(y=1.7, x=x, label=paste0('slope=',slope)), colour='black', size=5)+
#              geom_text(data=threeDensity, aes(y=1.4, x=x, label=paste0('intercept=',Intercept)), colour='black', size=5)+
#              theme_classic()+
#              theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#              ggtitle(opt$SampleName, 'Accessibility of all DHS surrounding the insertion (top)\nAccessibility of all 5kb windows with a DHS in contact with the insertion (bottom)')+
#              ylab('log10(RNA/DNA)')+
#              xlab('log10(total accessibility+1)')+
#              annotation_logticks(sides = "l",long = unit(0.2, "cm"))
#       
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/LinearVs3DDensity_scatter.pdf"), width=14, height=8)
#ggLinearVs3DDensity
#dev.off()
#
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/LinearVs3D_nDSH_vscontactsDHS.pdf"), width=14, height=4)
#par(mfrow=c(1,3))
#plot(threeDGenome$allDHS_1Mb, threeDGenome$contactsDHS_1Mb, xlab='Number of DHS surrounding the insertion', ylab='Number of 5kb windows with DHS in contact', main='1 Mb', frame=F, xlim=c(0,250), ylim=c(0,250))
#plot(threeDGenome$allDHS_500kb, threeDGenome$contactsDHS_500kb, xlab='Number of DHS surrounding the insertion', ylab='Number of 5kb windows with DHS in contact', main='500 kb', frame=F, xlim=c(0,150), ylim=c(0,150))
#plot(threeDGenome$allDHS_100kb, threeDGenome$contactsDHS_100kb, xlab='Number of DHS surrounding the insertion', ylab='Number of 5kb windows with DHS in contact', main='100 kb', frame=F, xlim=c(0,40), ylim=c(0,40))
#dev.off()
#
#
#pdf(file=paste0(opt$OUTDIR,"/graphs/LinearVs3D_nDSH_vscontactsDHSDensity.pdf"), width=14, height=4)
#par(mfrow=c(1,3))
#plot(log10(threeDGenome$allDHSdensity_1Mb+1),log10(threeDGenome$contactsDHSdensity_1Mb+1), xlab='Accessibility of surrounding windows  ', ylab='Accessibility of windows in contact', main='1 Mb', frame=F, xlim=c(0,4.5), ylim=c(0,4.5))
#plot(log10(threeDGenome$allDHSdensity_500kb+1),log10(threeDGenome$contactsDHSdensity_500kb+1), xlab='Accessibility of surrounding windows  ', ylab='Accessibility of windows in contact', main='500 kb', frame=F, xlim=c(0,4.5), ylim=c(0,4.5))
#plot(log10(threeDGenome$allDHSdensity_100kb+1),log10(threeDGenome$contactsDHSdensity_100kb+1), xlab='Accessibility of surrounding windows  ', ylab='Accessibility of windows in contact', main='100 kb', frame=F, xlim=c(0,4.5), ylim=c(0,4.5))
#dev.off()
#