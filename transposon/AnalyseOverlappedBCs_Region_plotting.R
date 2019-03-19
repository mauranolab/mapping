#!/bin/env Rscript

#Add variables from command line tsv files of DSnumbers and Institute
library("optparse")
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(scales)
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





####
#Regions 
####
filename <- paste0(opt$OUTDIR, "/AllRegions.annotated.bed")
cat('Reading ', filename,'\n')
dataRegions <- read(filename, header=T, stringsAsFactors=F)


#Correlation of signal
ggSignal <- ggplot(dataRegions, aes(x = rightSignal, y = leftSignal)) +
geom_point(color = "black", size = 0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nleft and right signal correlation")+
geom_smooth(method = "lm", se = FALSE)+
xlab('Right signal')+
ylab('Left signal')+
coord_cartesian(ylim = c(-2,2))+
annotation_logticks(sides = "lb",long = unit(0.2, "cm"))

pdf(file=paste0(opt$OUTDIR,"/Regiongraphs/Signal_Correlation.pdf"), width=5, height=5)
ggSignal
dev.off()



dataRegions$meanSignal <- (dataRegions$leftSignal+dataRegions$rightSignal)/2

ggCTCF <- ggplot(dataRegions, aes(x = nCTCFpeaks, y = meanSignal)) +
geom_point(color = "black", size = 0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nmean signal and nCTCF")+
geom_smooth(method = "lm", se = FALSE)+
xlab('number of CTCF peaks')+
ylab('mean signal')+
coord_cartesian(ylim = c(-2, 2),xlim=c(0,100))

pdf(file=paste0(opt$OUTDIR,"/Regiongraphs/MeanSignal_CTCF.pdf"), width=5, height=5)
ggCTCF
dev.off()

ggDHS <- ggplot(dataRegions, aes(x = nDHS, y = meanSignal)) +
geom_point(color = "black", size = 0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
#ggtitle(opt$SampleName,"  \nmean signal and nDHS")+
geom_smooth(method = "lm", se = FALSE)+
xlab('number of DHS peaks')+
ylab('mean signal')+
coord_cartesian(ylim = c(-2, 2),xlim=c(0,100))

pdf(file=paste0(opt$OUTDIR,"/Regiongraphs/MeanSignal_DHS.pdf"), width=5, height=5)
ggDHS
dev.off()

ggTSS <- ggplot(dataRegions, aes(x = nGencodeGenes, y = meanSignal)) +
geom_point(color = "black", size = 0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nmean signal and nTSS")+
geom_smooth(method = "lm", se = FALSE)+
xlab('number of gencode TSS')+
ylab('mean signal')+
coord_cartesian(ylim = c(-2, 2),xlim=c(0,100))

pdf(file=paste0(opt$OUTDIR,"/Regiongraphs/MeanSignal_TSS.pdf"), width=5, height=5)
ggTSS
dev.off()

ggExpressed <- ggplot(dataRegions, aes(x = nExpressedGenes, y = meanSignal)) +
geom_point(color = "black", size = 0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nmean signal and nExpressed genes")+
geom_smooth(method = "lm", se = FALSE)+
xlab('number of expressed gencode')+
ylab('mean signal')+
coord_cartesian(ylim = c(-2, 2),xlim=c(0,50))

pdf(file=paste0(opt$OUTDIR,"/Regiongraphs/MeanSignal_Expressed.pdf"), width=5, height=5)
ggExpressed
dev.off()

ggCGI <- ggplot(dataRegions, aes(x = nCGI, y = meanSignal)) +
geom_point(color = "black", size = 0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nmean signal and nCGI")+
geom_smooth(method = "lm", se = FALSE)+
xlab('number of CGI')+
ylab('mean signal')+
coord_cartesian(ylim = c(-2, 2),xlim=c(0,100))

pdf(file=paste0(opt$OUTDIR,"/Regiongraphs/MeanSignal_CGI.pdf"), width=5, height=5)
ggCGI
dev.off()


ggDpn <- ggplot(dataRegions, aes(x = nDpn, y = meanSignal)) +
geom_point(color = "black", size = 0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nmean signal and nDpn")+
geom_smooth(method = "lm", se = FALSE)+
xlab('number of Dpn')+
ylab('mean signal')+
coord_cartesian(ylim = c(-2, 2))+
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "b",long = unit(0.2, "cm"))

pdf(file=paste0(opt$OUTDIR,"/Regiongraphs/MeanSignal_Dpn.pdf"), width=5, height=5)
ggDpn
dev.off()


ggRegionLength <- ggplot(dataRegions, aes(x = regionWidth, y = meanSignal)) +
geom_point(color = "black", size = 0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=14), axis.text.x  = element_text(size=12), axis.title.y = element_text(size=14), axis.text.y = element_text(size=12))+
ggtitle(opt$SampleName,"  \nmean signal and Region width")+
geom_smooth(method = "lm", se = FALSE)+
xlab('regionWidth')+
ylab('mean signal')+
coord_cartesian(ylim = c(-2, 2))+
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "b",long = unit(0.2, "cm"))

pdf(file=paste0(opt$OUTDIR,"/Regiongraphs/MeanSignal_RegionLength.pdf"), width=5, height=5)
ggRegionLength
dev.off()


###
#Combine all plots
###
pdf(file=paste0(opt$OUTDIR,"/Regiongraphs/Regions_all.pdf"), width=23, height=3)
grid.arrange(ggDHS, ggCTCF, ggTSS, ggExpressed, ggCGI, ggDpn, ggRegionLength, ncol=7, nrow =1)
dev.off()

cat("\nWriting output files\n")
cat("Saving -- "); print(date())
save(list=setdiff(ls(), c("basenames", basenames)), file=paste0(opt$OUTDIR,"/allRegions.RData"), compress="bzip2")

cat("\nDone!!!\n")
print(date())



