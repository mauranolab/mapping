#!/bin/env Rscript

#TODO change ABC in column headers to DNA, RNA, iPCR and make separate sample column more like https://cascade.isg.med.nyu.edu/~mauram01/blog/tag_transposon.shtml#1535134189
#BUGBUG why does it replace hyphens in sample names with underscores?

suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(RColorBrewer)))
suppressWarnings(suppressMessages(library(gridExtra)))
suppressWarnings(suppressMessages(library(scales)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(gtools)))
options(bitmapType="cairo") 

basenames <- ls()
#source("~/.Rprofile")
option_list = list(
	make_option(c("-a", "--samplesA"), type="character", default=NULL, 
		help="Comma separated list of SamplesA", metavar="character"),
	make_option(c("-b", "--samplesB"), type="character", default=NULL, 
		help="Comma separated list of SamplesB", metavar="character"),
	make_option(c("-c", "--samplesC"), type="character", default=NULL, 
		help="Comma separated list of SamplesC (assumed to be iPCR)", metavar="character"),
	make_option(c("-d", "--samplesC2"), type="character", default=NULL, 
		help="Comma separated list of SamplesC (assumed to be iPCR)", metavar="character"),
	make_option(c("-t", "--threshold"), type="numeric", default=1, 
		help="thresholdRNADNA", metavar="character"),
	make_option(c("-o", "--output"), type="character", default=NULL, 
		help="output folder for all comparisons", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$samplesA) & is.null(opt$samplesC2)){
	print_help(opt_parser)
	stop("Samples A must be provided\n", call.=FALSE)
}

if (!is.null(opt$samplesB)){
	message("Samples B are provided\n")
}

if (!is.null(opt$samplesC)){
	message("Samples C are provided \n")
}


if (!is.null(opt$samplesC2)){
	message("Samples C2 are provided \n")
}


if (is.null(opt$output)){
	stop("Output file needs to be provided", call.=FALSE)
}

if (is.null(opt$threshold)){
	thresholdRNADNA <- 10
} else { 
	thresholdRNADNA <- opt$threshold
}
thresholdiPCR = 1

if (!is.null(opt$samplesA)){
	samplesA <-strsplit(opt$samplesA, split=',')[[1]]
	message(samplesA[1], '\n')
}


#
#samplesA <- c("aligned.FCHVTCCBGX2/BS00370A-MSH_060517_K562~CMV~2roundLinear_Rep1_RNA", 'aligned.FCHVTCCBGX2/BS00372A-MSH_060517_K562~CMV~1roundLinear_Rep1_RNA', 'aligned.FCHVTCCBGX2/BS00372A-MSH_060517_K562~CMV~1roundLinear_Rep1_RNA')
#samplesB <- c('aligned.FCHVTCCBGX2/BS00371A-MSH_060517_K562~CMV~2roundLinear_Rep2_RNA', 'aligned.FCHVTCCBGX2/BS00373A-MSH_060517_K562~CMV~1roundLinear_Rep2_RNA', "aligned.FCHVTCCBGX2/BS00375A-MSH_060517_K562~GGloHS2~2roundLinear_Rep2_RNA")
#samplesC <- c('aligned.Merged/BS00244-K562_CMVBC4_iPCR_Merged', 'aligned.Merged/BS00244-K562_CMVBC4_iPCR_Merged', 'aligned.Merged/BS00244-K562_CMVBC4_iPCR_Merged')


#samplesA <- c(samplesA, samplesB)
dir.create(dirname(opt$output))
#samplesA <- c("aligned.FCH7HTKBGX3/BS00490A-MSH_K562~A1~GGlo~HS2~A1~1Linear~NoExo_Rep1_DNA, aligned.FCH7HTKBGX3/BS00491A-MSH_K562~A1~GGlo~HS2~A1~1Linear~NoExo_Rep2_DNA, aligned.FCH7HTKBGX3/BS00492A-MSH_K562~A1~GGlo~HS2~A1~1Linear~ExoI_Rep1_DNA, aligned.FCH7HTKBGX3/BS00493A-MSH_K562~A1~GGlo~HS2~A1~1Linear~ExoI_Rep2_DNA, aligned.Merged/BS00467_468A--MSH-K562~fw_GGlo_fw_OldProtocol-DNA_Merged, aligned.Merged/BS00217-K562d11_A1fw_GGlo_A1fw_HS2_BC4_1000ng_25Ksort_GTG_DNA_Merged")
#samplesA <-strsplit(samplesA, split=',')[[1]]


loadBCfile <- function(filename, type) {
	if(type=="iPCR") {
		suffix <- ".barcodes.coords.bed"
		col.names <- c("chrom", "chromStart", "chromEnd", "bc", "count", "strand")
	} else if (type %in% c("DNA", "RNA")) {
		suffix <- ".barcode.counts.UMI.corrected.txt"
		col.names <- c("bc", "count")
	} else {
		stop("Do not recognize file type", type)
	}
	data <- read(paste0(filename, '/', basename(filename), suffix), header=F)
	colnames(data) <- col.names
	return(data)
}

#TODO no arg handling
rnadnapair <- function(BCpairs, typeA, typeB, threshA, threshB) {
	message("\nComparing ", typeA, " vs. ", typeB, ". Thresholds=", threshA, "/", threshB)
	
	usefulBCs <- data.frame(matrix(ncol=6, nrow=nrow(BCpairs)))
	colnames(usefulBCs) <- c("Sample", paste0("#A_", typeA, "_", threshA), paste0("#B_", typeB, "_", threshB), paste0("#Intersect_A", threshA, "_B", threshB), paste0("#Union_A", threshA, "_B", threshB), paste0("Intersect/Union_A", threshA, "_B", threshB))
	
	for (i in 1:nrow(BCpairs)){
		sampleA <- basename(BCpairs$A[i])
		sampleB <- basename(BCpairs$B[i])
		
		#TODO need better sample name than gsub with typeA
		sample <- gsub(paste0('_', typeA), '', unlist(strsplit(sampleA, split='-'))[2])
		bsA <- unlist(strsplit(sampleA, split='-'))[1]
		bsB <- unlist(strsplit(sampleB, split='-'))[1]
		
		cat(sample, bsA, bsB, '\n')
		
		sampleA_BC <- loadBCfile(BCpairs$A[i], type=typeA)
		sampleB_BC <- loadBCfile(BCpairs$B[i], type=typeB)
		
		sampleA_BC <- sampleA_BC[sampleA_BC$count >= threshA,]
		sampleB_BC <- sampleB_BC[sampleB_BC$count >= threshB,]
		
		sampleA_sampleB <- intersect(sampleA_BC$bc, sampleB_BC$bc)
		usampleA_sampleB <- union(sampleA_BC$bc, sampleB_BC$bc)
		
		usefulBCs[,1][i] <- sample
		#TODO split out BS numbers for each sample
		usefulBCs[,2][i] <- format(as.numeric(nrow(sampleA_BC)), big.mark=",", trim=TRUE)
		usefulBCs[,3][i] <- format(as.numeric(nrow(sampleB_BC)), big.mark=",", trim=TRUE)
		usefulBCs[,4][i] <- format(as.numeric(length(sampleA_sampleB)), big.mark=",", trim=TRUE)
		usefulBCs[,5][i] <- format(as.numeric(length(usampleA_sampleB)), big.mark=",", trim=TRUE)
		usefulBCs[,6][i] <- format(as.numeric(round(length(sampleA_sampleB)/length(usampleA_sampleB), 2)), big.mark=",", trim=TRUE)
	}
	return(usefulBCs)
}


#####
#If only A exists. Compare all A
#####
if (!is.null(opt$samplesA) & is.null(opt$samplesB) & is.null(opt$samplesC)) {
	BCpairs <- as.data.frame(combinations(n = length(samplesA), r = 2, v = samplesA, repeats.allowed = FALSE))
	colnames(BCpairs) <- c("A", "B")
	BCpairs$A <- as.character(BCpairs$A)
	BCpairs$B <- as.character(BCpairs$B)
	BCpairs <- BCpairs[BCpairs$A != BCpairs$B,]
	usefulBCs <- rnadnapair(BCpairs, "DNA", "DNA", 0, 0)
	usefulBCs <- merge(usefulBCs, rnadnapair(BCpairs, "DNA", "DNA", thresholdRNADNA, thresholdRNADNA), by="Sample")
	write.table(usefulBCs[order(usefulBCs[,1]),], file=opt$output, row.names=F, sep='\t', quote=F)
}


#####
#If A and B exist
#####
if (!is.null(opt$samplesA) & !is.null(opt$samplesB) & is.null(opt$samplesC)) {
	samplesB <- strsplit(opt$samplesB, split=',')[[1]]
	message(samplesB[1], '\n')
	BCpairs <- data.frame("A"=samplesA, "B"=samplesB, stringsAsFactors=F)
	usefulBCs <- rnadnapair(BCpairs, "DNA", "RNA", 0, 0)
	usefulBCs <- merge(usefulBCs, rnadnapair(BCpairs, "DNA", "RNA", thresholdRNADNA, thresholdRNADNA), by="Sample")
	write.table(usefulBCs[order(usefulBCs[,1]),], file=opt$output, row.names=F, sep='\t', quote=F)
}



######
#If A, B, and C exists
######
if (!is.null(opt$samplesA) & !is.null(opt$samplesB) & !is.null(opt$samplesC)) {
	samplesB <- strsplit(opt$samplesB, split=',')[[1]]
	message(samplesB[1], '\n')
	samplesC <- strsplit(opt$samplesC, split=',')[[1]]
	BCpairs <- data.frame("A"=samplesA, "B"=samplesB, "C"=samplesC, stringsAsFactors=F)


#TODO too many columns
#usefulBCs <- rnadnapair(BCpairs, "DNA", "RNA", thresholdRNADNA, thresholdRNADNA)
#usefulBCs <- rnadnapair(BCpairs, "DNA", "iPCR", thresholdRNADNA, thresholdiPCR)
#usefulBCs <- rnadnapair(BCpairs, "RNA", "iPCR", thresholdRNADNA, thresholdiPCR)
#missing 3way

###
	usefulBCs <- data.frame(matrix(ncol=13, nrow=nrow(BCpairs)))
	colnames(usefulBCs) <- c('A', 'B', 'C', '#A', '#B', '#C_mapped', paste0('#A_', thresholdRNADNA), paste0('#B_', thresholdRNADNA), paste0('#C_', thresholdiPCR), paste0('#Intersect_A', thresholdRNADNA, '_B', thresholdRNADNA), paste0('#Intersect_A', thresholdRNADNA, '_C', thresholdiPCR), paste0('#Intersect_B', thresholdRNADNA, '_C', thresholdiPCR), paste0('#Intersect_A', thresholdRNADNA, '_B', thresholdRNADNA, '_C', thresholdiPCR))
	
	for (i in 1:nrow(BCpairs)) {
		samA_BC <- read(paste0(BCpairs$A[i], '/', gsub('.*/', '', BCpairs$A[i]), '.barcode.counts.UMI.corrected.txt'), header=F)
		colnames(samA_BC) <- c("bc", "count")
		
		samB_BC <- read(paste0(BCpairs$B[i], '/', gsub('.*/', '', BCpairs$B[i]), '.barcode.counts.UMI.corrected.txt'), header=F)
		colnames(samB_BC) <- c("bc", "count")
		
		samC <- read(paste0(BCpairs$C[i], '/', gsub('.*/', '', BCpairs$C[i]), '.barcodes.coords.bed'), header=F)
		colnames(samC) <- c("chrom", "chromStart", "chromEnd", "bc", "count", "strand")
		
		samC_BC1 <- samC[samC$count >= thresholdiPCR,]
		samA_BC10 <- samA_BC[samA_BC$count >= thresholdRNADNA,]
		samB_BC10 <- samB_BC[samB_BC$count >= thresholdRNADNA,]
		
		samA_BC10 <- samA_BC[samA_BC$count >= thresholdRNADNA,]
		samB_BC10 <- samB_BC[samB_BC$count >= thresholdRNADNA,]
		
		samA_samB_10 <- intersect(samA_BC10$bc, samB_BC10$bc)
		usamA_samB_10 <- union(samA_BC10$bc, samB_BC10$bc)
		
		
		samA10_samC2 <- intersect(samA_BC10$bc, samC_BC1$bc)
		samB10_samC2 <- intersect(samC_BC1$bc, samB_BC10$bc)
		
		samB10_samC2_samA10 <- intersect(samB10_samC2, samA_BC10$bc)
		#Populate data frame
		usefulBCs[,1][i] <- gsub('_', '-', gsub('.*/', '', BCpairs$A[i]))
		usefulBCs[,2][i] <- gsub('_', '-', gsub('.*/', '', BCpairs$B[i]))
		usefulBCs[,3][i] <- gsub('_', '-', gsub('.*/', '', BCpairs$C[i]))
		usefulBCs[,4][i] <- nrow(samA_BC)
		usefulBCs[,5][i] <- nrow(samB_BC)
		usefulBCs[,6][i] <- nrow(samC)
		usefulBCs[,7][i] <- nrow(samA_BC10)
		usefulBCs[,8][i] <- nrow(samB_BC10)
		usefulBCs[,9][i] <- nrow(samC_BC1)
		usefulBCs[,10][i] <- length(samA_samB_10)
		usefulBCs[,11][i] <- length(samA10_samC2)
		usefulBCs[,12][i] <- length(samB10_samC2)
		usefulBCs[,13][i] <- length(samB10_samC2_samA10)
	}
	for (i in c(4:13)) {
		usefulBCs[,i] <- format(as.numeric(usefulBCs[,i]), big.mark=",", trim=TRUE)
	}
	
	usefulBCs <- rnadnapair(BCpairs, "DNA", "RNA", 0, 0)
###
	
	
	write.table(usefulBCs[order(usefulBCs[,1]),], file=opt$output, row.names=F, sep='\t', quote=F)
}



#####
#If A and C exists
#####
if (!is.null(opt$samplesA) & is.null(opt$samplesB) & !is.null(opt$samplesC)) {
	samplesC <- strsplit(opt$samplesC, split=',')[[1]]
	message(samplesC[1], '\n')
	BCpairs <- data.frame("A"=samplesA, "B"=samplesC, stringsAsFactors=F)
	usefulBCs <- rnadnapair(BCpairs, "DNA", "iPCR", 0, 0)
	usefulBCs <- merge(usefulBCs, rnadnapair(BCpairs, "DNA", "iPCR", thresholdRNADNA, thresholdiPCR), by="Sample")
	write.table(usefulBCs[order(usefulBCs[,1]),], file=opt$output, row.names=F, sep='\t', quote=F)
}


#####
#Compare 2 iPCR samples
#####
if (is.null(opt$samplesC2) & !is.null(opt$samplesC) & !is.null(opt$samplesC2)) {
	samplesC <- strsplit(opt$samplesC, split=',')[[1]]
	samplesC2 <- strsplit(opt$samplesC2, split=',')[[1]]
	message(samplesC[1], '\n')
	BCpairs <- data.frame("C"=samplesC, "C2"=samplesC2, stringsAsFactors=F)
	usefulBCs <- rnadnapair(BCpairs, "iPCR", "iPCR", thresholdiPCR, thresholdiPCR)
	write.table(usefulBCs[order(usefulBCs[,1]),], file=opt$output, row.names=F, sep='\t', quote=F)
}


