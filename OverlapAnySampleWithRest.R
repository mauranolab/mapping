#!/bin/env Rscript

#Add variables from command line tsv files of DSnumbers and Institute
suppressWarnings(suppressMessages(library("optparse")))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(scales)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(pheatmap)))
suppressWarnings(suppressMessages(library(gtools)))
options(bitmapType="cairo") 

basenames <- ls()
#source("~/.Rprofile")
option_list = list(
       make_option(c("-s", "--samples"), type="character", default=NULL, 
              help="Comma seperated list of Samples to compare agains each other", metavar="character"),
       make_option(c("-t", "--threshold"), type="numeric", default=1, 
              help="Comma seperated list of SamplesB", metavar="character"),
       make_option(c("-f", "--File"), action="store_true", type="character", default=FALSE, 
              help="Samples are in a new line separated file", metavar="character"),
       make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output folder for all comparisons", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$samples)){
       print_help(opt_parser)
       stop("Samples A must be provided\n", call.=FALSE)
}


if (is.null(opt$threshold)){
       threshold <- 1
} else { 
       threshold <- opt$threshold
}


if (is.null(opt$output)){
  stop("Output file needs to be provided", call.=FALSE)
}


#threshold = 1

#samples <- c('aligned.FCHFNJJBGX5/BS00734B-BJM_20171011_GGlo~A1fw~HS2_1X_50k_iPCR,aligned.FCHFNJJBGX5/BS00769A-RDL_20180108_K562_T0049_A1f~GGlo~A1f~HS2_4C_iPCR,aligned.FCHFNJJBGX5/BS00772A-RDL_20171211_K562_T0055_LTR7_d1sort_iPCR,aligned.FCHFNJJBGX5/BS00748A-RDL_20171211_K562_Hsp68_day5_GFPHigh_RNA,aligned.FCHFNJJBGX5/BS00745B-RDL_20171211_K562_HSP68_day1_RNA,aligned.FCHFNJJBGX5/BS00742B-BJM_20171011_GGlo_1_2_50k_iPCR,aligned.FCHFNJJBGX5/BS00770A-RDL_20171106_K562_T0053_Hsp68_Day10_iPCR,aligned.FCHFNJJBGX5/BS00739B-BJM_20171011_A1fw~GGlo~A1fw~HS2_2X_50k_iPCR,aligned.FCHFNJJBGX5/BS00733B-BJM_20171011_GGlo_1X_50k_iPCR,aligned.FCHFNJJBGX5/BS00736B-BJM_20171011_BGlo_1X_50k_iPCR,aligned.FCHFNJJBGX5/BS00768A-RDL_20180108_K562_T0046_GGlo~HS2_4C_iPCR,aligned.FCHFNJJBGX5/BS00767A-RDL_20180108_K562_T0045_GGlo_4C_iPCR,aligned.FCHFNJJBGX5/BS00738B-BJM_20171011_GGlo~A1fw~HS2_2X_50k_iPCR,aligned.FCHFNJJBGX5/BS00740B-BJM_20171011_BGlo_2X_50k_iPCR,aligned.FCHFNJJBGX5/BS00745A-RDL_20171211_K562_HSP68_day1_RNA,aligned.FCHFNJJBGX5/BS00737B-BJM_20171011_GGlo_2X_50k_iPCR,aligned.FCHFNJJBGX5/BS00774A-RDL_20171106_K562_T0053_Hsp68_Day10_DNA,aligned.FCHFNJJBGX5/BS00735B-BJM_20171011_A1fw~GGlo~A1fw~HS2_1X_50k_iPCR,aligned.FCHFNJJBGX5/BS00775A-RDL_20171211_K562_T0056_CMV_d1sort_DNA,aligned.FCHFNJJBGX5/BS00773A-RDL_20171211_K562_T0056_CMV_d1sort_iPCR,aligned.FCHFNJJBGX5/BS00741B-BJM_20171011_GGlo_1_10_50k_iPCR,aligned.FCHFNJJBGX5/BS00746B-RDL_20171211_K562_CMV_day5_GFPHigh_RNA,aligned.FCHFNJJBGX5/BS00657D-RDL_20171026_K562_Hsp68_Day10_RNA,aligned.FCHFNJJBGX5/BS00771A-RDL_20171106_K562_T0054_BGlo_Day10_iPCR')

if (!is.null(opt$samples)){
       if (opt$File == TRUE) {
              samples <- readLines(opt$samples)
       } else {       
              samples <-strsplit(opt$samples,split=',')[[1]]
       }
       cat(samples[1], '\n')

}


dir.create(dirname(opt$output))


#For now remove all duplicated sample names
samples <- samples[!duplicated(gsub('.*/', '', samples))]
samples <-  gsub('_060517_', '_20170605_', samples)
#####
#For each sample in A compare with all samples in B
#####

#BCpairs <- as.data.frame(combinations(n = length(samples), r = 2, v = samples, repeats.allowed = FALSE))
BCpairs <- expand.grid(samples, samples, stringsAsFactors = FALSE)
colnames(BCpairs) <- c('V1', 'V2')
BCpairs$V1 <- as.character(BCpairs$V1)
BCpairs$V2 <- as.character(BCpairs$V2)
#BCpairs <- BCpairs[BCpairs$V1!=BCpairs$V2,]

usefulBCs<-data.frame(matrix(ncol=12,nrow=nrow(BCpairs)))
colnames(usefulBCs)<-c('A','B','#A','#B','Intersect','Union','Intersect/Union',paste0('#A_',threshold),paste0('#B_',threshold),paste0('Intersect_',threshold),paste0('Union_',threshold),paste0('Intersect_',threshold,'/','Union_',threshold))



for (i in 1:nrow(BCpairs)){

       cat(BCpairs$V1[i],BCpairs$V2[i],'\n')
       if (length(grep('iPCR', BCpairs$V1[i])) >0) {
              samA_BC<-fread(paste0("/home/maagj01/scratch/transposon/",BCpairs$V1[i],'/',gsub('.*/','',BCpairs$V1[i]),'.barcodes.coords.bed'),header=F,stringsAsFactors=F)
              samA_BC <- subset(samA_BC, select=c(V4, V5))
              colnames(samA_BC) <- c('V1', 'V2')
              # To avoid filtering iPCR samples, set barcode count to threshold + 1. We already filter for 2 Barcode for mapped iPCR
              samA_BC$V2 <- threshold+1
       } else {
              samA_BC<-fread(paste0("/home/maagj01/scratch/transposon/",BCpairs$V1[i],'/',gsub('.*/','',BCpairs$V1[i]),'.barcode.counts.txt'),header=F,stringsAsFactors=F)
       }
       
       

       if (length(grep('iPCR', BCpairs$V2[i])) >0) {
              samB_BC<-fread(paste0("/home/maagj01/scratch/transposon/",BCpairs$V2[i],'/',gsub('.*/','',BCpairs$V2[i]),'.barcodes.coords.bed'),header=F,stringsAsFactors=F)
              samB_BC <- subset(samB_BC, select=c(V4, V5))
              colnames(samB_BC) <- c('V1', 'V2')
              # To avoid filtering iPCR samples, set barcode count to threshold + 1. We already filter for 2 Barcode for mapped iPCR
              samB_BC$V2 <- threshold+1
       } else {
              samB_BC<-fread(paste0("/home/maagj01/scratch/transposon/",BCpairs$V2[i],'/',gsub('.*/','',BCpairs$V2[i]),'.barcode.counts.txt'),header=F,stringsAsFactors=F)
       }


       samA_samB<-intersect(samA_BC$V1,samB_BC$V1)
       usamA_samB<-union(samA_BC$V1,samB_BC$V1)

       samA_BC10<-samA_BC[samA_BC$V2>=threshold,]
       samB_BC10<-samB_BC[samB_BC$V2>=threshold,]
       
       samA_samB_10<-intersect(samA_BC10$V1,samB_BC10$V1)
       usamA_samB_10<-union(samA_BC10$V1,samB_BC10$V1)
              #Populate data frame
       usefulBCs[,1][i]<-gsub('_','-',gsub('.*/','',BCpairs$V1[i]))
       usefulBCs[,2][i]<-gsub('_','-',gsub('.*/','',BCpairs$V2[i]))
       usefulBCs[,3][i]<-nrow(samA_BC)
       usefulBCs[,4][i]<-nrow(samB_BC)
       usefulBCs[,5][i] <- length(samA_samB)
       usefulBCs[,6][i] <- length(usamA_samB)
       usefulBCs[,7][i] <- round(length(samA_samB)/length(usamA_samB),2)
       
       usefulBCs[,8][i]<-nrow(samA_BC10)
       usefulBCs[,9][i]<-nrow(samB_BC10)

       usefulBCs[,10][i]<-length(samA_samB_10)
       usefulBCs[,11][i]<-length(usamA_samB_10)
       usefulBCs[,12][i] <- round(length(samA_samB_10)/length(usamA_samB_10),2)
}



write.table(file=paste0(opt$output, '.tsv'), usefulBCs, quote=F, sep='\t', col.names=T, row.names=F)



cont <- usefulBCs

colnames(cont) <- paste0('V', 1:ncol(cont))

colnames(cont)
#Keep dates for Raven
cont$name1 <- gsub('T0053-|T0052-|T0054-|T0055-|T0056-|T0058-|T0059-|T0060-|T0061-|-2xgDNA', '', paste0( gsub('.*-K562-|.*-K562~|.*-BJM-20171011-', '',cont$V1), '-', gsub('-.*', '', cont$V1), '-', gsub('.*-', '', strtrim(cont$V1, width=21))))
cont$name2 <- gsub('T0053-|T0052-|T0054-|T0055-|T0056-|T0058-|T0059-|T0060-|T0061-|-2xgDNA', '', paste0( gsub('.*-K562-|.*-K562~|.*-BJM-20171011-', '',cont$V2), '-', gsub('-.*', '', cont$V2), '-', gsub('.*-', '', strtrim(cont$V2, width=21))))

cont$name1  <- gsub('K$', '20170605', cont$name1)
cont$name2 <- gsub('K$', '20170605', cont$name2)

cont$name1[grep('NEM', cont$V1)] <- gsub('Gglo', 'Nich', cont$name1[grep('NEM', cont$V1)])
cont$name2[grep('NEM', cont$V2)] <- gsub('Gglo', 'Nich', cont$name2[grep('NEM', cont$V2)])
cont$name1<- gsub('A1fw', 'A1f', cont$name1)
cont$name2<- gsub('A1fw', 'A1f', cont$name2)
cont$name1<- gsub('GGlo-DNA-', 'GGlo-1-DNA-', cont$name1)
cont$name2<- gsub('GGlo-DNA-', 'GGlo-1-DNA-', cont$name2)

cont$name1 <- gsub('GFP', '',  gsub('sort', '', gsub('Day|day', 'd', cont$name1)))
cont$name2 <- gsub('GFP', '',  gsub('sort', '', gsub('Day|day', 'd', cont$name2)))


#cont$name1<- gsub('~', '-', cont$name1)
#cont$name2<- gsub('~', '-', cont$name2)
#

contDF <- dcast(data = cont,formula = name1~name2,fun.aggregate = mean,value.var = "V12")


rownames(contDF) <-contDF$name1

contDF <- contDF[,-1]

#contDF[contDF<100] <-0
 
# 
#pheatmap(log10(contDF+1),  color = colorRampPalette(brewer.pal(n = 9, name ="Reds"))(100), main= "Log10(Shared barcodes)", fontsize_row = 8, fontsize_col = 8)
#
#pheatmap(log10(contDF[,colnames(contDF) %in% rownames(contDF)]+1),  color = colorRampPalette(brewer.pal(n = 9, name ="Reds"))(100), main= "Log10(Shared barcodes)", fontsize_row = 9, fontsize_col = 9)
#


#contDFSubset <- contDF[,colnames(contDF) %in% rownames(contDF)]
#contDFSubset <- contDFSubset[rownames(contDFSubset) %in% colnames(contDFSubset),]
contDFSubset <- contDF[,rownames(contDF)[match(colnames(contDF), rownames(contDF))]]



mat_col <- data.frame(Type = gsub('.*-', '', gsub("-BS.*", '', rownames(contDFSubset))),
                     Date = str_extract(samples, "2018[0-9]{4}|2017[0-9]{4}|2016[0-9]{4}|2019[0-9]{4}"),
                     #Sample = gsub('Nich.*', 'Nich', 
                     #              gsub('HSP68', 'Hsp68', 
                     #              gsub('Gglo', 'GGlo', 
                     #              gsub('A1fw', 'A1f',
                     #              gsub('-DNA.*|-RNA.*|-iPCR.*', '', rownames(contDFSubset)))))),
                     Transfection = str_extract(rownames(contDFSubset), "T[0-9]{4}"))

#mat_col$Sample <- gsub('HSP68', 'Hsp68', gsub('Gglo', 'GGlo', gsub('A1fw', 'A1f',gsub('-DNA.*|-RNA.*|-iPCR.*', '', rownames(contDFSubset)))))

#mat_col$Sample[grep('Nich', mat_col$Sample)] <- 'Nich'

rownames(contDFSubset) <- gsub('-20.*', '', rownames(contDFSubset) )
colnames(contDFSubset) <- gsub('-20.*', '', colnames(contDFSubset) )
rownames(mat_col) <- rownames(contDFSubset)
mat_col$Type <- as.factor(mat_col$Type)

mat_colors <- list(Type = brewer.pal( 8, "Set1")[1:length(levels(mat_col$Type))],
                     Date= c(brewer.pal( 8, "Accent"),brewer.pal( 8, "Dark2"))[1:length(levels(mat_col$Date))],
                     #Sample = c(brewer.pal( 9, "Set1"), brewer.pal(9, "Set3"), brewer.pal(8, "Set2"), brewer.pal(5, "Dark2"))[1:length(levels(mat_col$Sample))],
                     Transfection = c(brewer.pal( 8, "Pastel1"),brewer.pal( 8, "Pastel2"))[1:length(levels(mat_col$Transfection))])
names(mat_colors$Type) <- levels(mat_col$Type)
names(mat_colors$Date) <-  levels(mat_col$Date)
names(mat_colors$Sample) <-  levels(mat_col$Sample)
names(mat_colors$Transfection) <-  levels(mat_col$Transfection)
#as.factor(names(mat_colors$Type))


#pdf('~/public_html/blog/2017Dec04/CMV_featureHeaptmap.pdf', height=8, width=12, onefile=F) 
#pheatmap(log10(test+1), scale='column', annotation_row=mat_col, annotation_colors= mat_colors)
#dev.off()


contDFSubset[contDFSubset>0.3] <- 0.3


plotWidth = nrow(contDFSubset)*0.4
pdf(paste0(opt$output, '.pdf'), width=plotWidth+1, height=plotWidth, onefile=F)
pheatmap(contDFSubset,  
       color = colorRampPalette(brewer.pal(n = 9, name ="Reds"))(100), 
       main= "Contamination between samples\nJaccard Index", 
       fontsize_row = 8, fontsize_col = 8, 
       cluster_rows = FALSE, cluster_cols = FALSE ,
       annotation_row=mat_col, annotation_col=mat_col, annotation_colors= mat_colors)
dev.off()
#
##
##pdf('~/public_html/blog/20180326/OverlapContamination_clustered.pdf', width=12, height=12, onefile=F)
##pheatmap(log10(contDFSubset+1),  
##       color = colorRampPalette(brewer.pal(n = 9, name ="Reds"))(100), 
#       main= "Contamination between samples\nlog10(Shared barcodes)\n(only shows barcode overlap >100)", 
#       fontsize_row = 10, fontsize_col = 10, 
#       cluster_rows = T, cluster_cols = T ,
#       annotation_row=mat_col, annotation_col=mat_col, annotation_colors= mat_colors)
#dev.off()
#


#contLower <- as.matrix(contLower)
#colnames(contLower) <- colnames(contDFSubset)
#rownames(contLower) <- rownames(contDFSubset)
#
#contLower[contLower>0.3] <- 0.3
#
#
#pdf(paste0(opt$output, '.tsv'), width=12, height=12, onefile=F)
#pheatmap(contLower,  
#       #color = brewer.pal(n = 9, name ="Reds"), 
#       color = colorRampPalette(brewer.pal(n = 9, name ="Reds"))(100),
#       main= "Contamination between samples\nOverlap/Union", 
#       fontsize_row = 10, fontsize_col = 10, 
#       cluster_rows = FALSE, cluster_cols = FALSE ,
#       annotation_row=mat_col, annotation_col=mat_col, annotation_colors= mat_colors)
#dev.off()
