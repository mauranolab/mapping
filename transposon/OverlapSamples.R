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
suppressWarnings(suppressMessages(library(data.table)))
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

if (is.null(opt$samplesB)){
	print_help(opt_parser)
	cat("Samples B must be provided\n")
}

if (!is.null(opt$samplesC)){
	cat("Samples C are provided \n")
}


if (!is.null(opt$samplesC2)){
	cat("Samples C2 are provided \n")
}


if (is.null(opt$output)){
	stop("Output file needs to be provided", call.=FALSE)
}

if (is.null(opt$threshold)){
	thresholdRNADNA <- 1
} else { 
	thresholdRNADNA <- opt$threshold
}
thresholdIPCR = 1

if (!is.null(opt$samplesA)){
	samplesA <-strsplit(opt$samplesA, split=',')[[1]]
	cat(samplesA[1], '\n')
}


#
#samplesA <- c("aligned.FCHVTCCBGX2/BS00370A-MSH_060517_K562~CMV~2roundLinear_Rep1_RNA", 'aligned.FCHVTCCBGX2/BS00372A-MSH_060517_K562~CMV~1roundLinear_Rep1_RNA', 'aligned.FCHVTCCBGX2/BS00372A-MSH_060517_K562~CMV~1roundLinear_Rep1_RNA')
#samplesB <- c('aligned.FCHVTCCBGX2/BS00371A-MSH_060517_K562~CMV~2roundLinear_Rep2_RNA', 'aligned.FCHVTCCBGX2/BS00373A-MSH_060517_K562~CMV~1roundLinear_Rep2_RNA', "aligned.FCHVTCCBGX2/BS00375A-MSH_060517_K562~GGloHS2~2roundLinear_Rep2_RNA")
#samplesC <- c('aligned.Merged/BS00244-K562_CMVBC4_iPCR_Merged', 'aligned.Merged/BS00244-K562_CMVBC4_iPCR_Merged', 'aligned.Merged/BS00244-K562_CMVBC4_iPCR_Merged')


#samplesA <- c(samplesA, samplesB)
dir.create(dirname(opt$output))
#samplesA <- c("aligned.FCH7HTKBGX3/BS00490A-MSH_K562~A1~GGlo~HS2~A1~1Linear~NoExo_Rep1_DNA, aligned.FCH7HTKBGX3/BS00491A-MSH_K562~A1~GGlo~HS2~A1~1Linear~NoExo_Rep2_DNA, aligned.FCH7HTKBGX3/BS00492A-MSH_K562~A1~GGlo~HS2~A1~1Linear~ExoI_Rep1_DNA, aligned.FCH7HTKBGX3/BS00493A-MSH_K562~A1~GGlo~HS2~A1~1Linear~ExoI_Rep2_DNA, aligned.Merged/BS00467_468A--MSH-K562~fw_GGlo_fw_OldProtocol-DNA_Merged, aligned.Merged/BS00217-K562d11_A1fw_GGlo_A1fw_HS2_BC4_1000ng_25Ksort_GTG_DNA_Merged")
#samplesA <-strsplit(samplesA, split=',')[[1]]


#####
#If only A exists. Compare all A
#####
if (!is.null(opt$samplesA) & is.null(opt$samplesC) & is.null(opt$samplesB)) {
	BCpairs <- as.data.frame(combinations(n = length(samplesA), r = 2, v = samplesA, repeats.allowed = FALSE))
	BCpairs$A <- as.character(BCpairs$A)
	BCpairs$B <- as.character(BCpairs$B)
	BCpairs <- BCpairs[BCpairs$A != BCpairs$B,]
	
	usefulBCs <- data.frame(matrix(ncol=12, nrow=nrow(BCpairs)))
	colnames(usefulBCs) <- c('A', 'B', '#A', '#B', '#Intersect', '#Union', 'Intersect/Union', paste0('#A_', thresholdRNADNA), paste0('#B_', thresholdRNADNA), paste0('#Intersect_', thresholdRNADNA), paste0('#Union_', thresholdRNADNA), paste0('Intersect_', thresholdRNADNA, '/', 'Union_', thresholdRNADNA))
	
	for (i in 1:nrow(BCpairs)){
		cat(BCpairs$A[i], BCpairs$B[i], '\n')
		
		samA_BC <- fread(paste0(BCpairs$A[i], '/', gsub('.*/', '', BCpairs$A[i]), '.barcode.counts.UMI.corrected.txt'), header=F, stringsAsFactors=F)
		colnames(samA_BC) <- c("bc", "count")
		
		samB_BC <- fread(paste0(BCpairs$B[i], '/', gsub('.*/', '', BCpairs$B[i]), '.barcode.counts.UMI.corrected.txt'), header=F, stringsAsFactors=F)
		colnames(samB_BC) <- c("bc", "count")
		
		
		samA_samB <- intersect(samA_BC$bc, samB_BC$bc)
		usamA_samB <- union(samA_BC$bc, samB_BC$bc)
		
		samA_BC10 <- samA_BC[samA_BC$count >= thresholdRNADNA,]
		samB_BC10 <- samB_BC[samB_BC$count >= thresholdRNADNA,]
		
		samA_samB_10 <- intersect(samA_BC10$bc, samB_BC10$bc)
		usamA_samB_10 <- union(samA_BC10$bc, samB_BC10$bc)
		#Populate data frame
		usefulBCs[,1][i] <- gsub('_', '-', gsub('.*/', '', BCpairs$A[i]))
		usefulBCs[,2][i] <- gsub('_', '-', gsub('.*/', '', BCpairs$B[i]))
		usefulBCs[,3][i] <- nrow(samA_BC)
		usefulBCs[,4][i] <- nrow(samB_BC)
		usefulBCs[,5][i] <- length(samA_samB)
		usefulBCs[,6][i] <- length(usamA_samB)
		usefulBCs[,7][i] <- round(length(samA_samB)/length(usamA_samB), 2)
		usefulBCs[,8][i] <- nrow(samA_BC10)
		usefulBCs[,9][i] <- nrow(samB_BC10)
		usefulBCs[,10][i] <- length(samA_samB_10)
		usefulBCs[,11][i] <- length(usamA_samB_10)
		usefulBCs[,12][i] <- round(length(samA_samB_10)/length(usamA_samB_10), 2)
	}
	for (i in c(3:6, 8:11)){
		usefulBCs[,i] <- format(as.numeric(usefulBCs[,i]), big.mark=",", trim=TRUE)
	}
	for (i in 1:ncol(usefulBCs)) {
		if (diff(as.numeric(factor(usefulBCs[,i])))[-1]==0) {
			usefulBCs[c(FALSE, diff(as.numeric(factor(usefulBCs[,i])))==0),][,i] <- '"'
		}
	}
	write.table(usefulBCs[order(usefulBCs[,1]),], file=opt$output, row.names=F, sep='\t', quote=F)
}



#####
#If A and B exists
#####
if (!is.null(opt$samplesA) & !is.null(opt$samplesB) & is.null(opt$samplesC)) {
	samplesB <- strsplit(opt$samplesB, split=',')[[1]]
	cat(samplesB[1], '\n')
	BCpairs <- data.frame("A"=samplesA, "B"=samplesB, stringsAsFactors=F)
	usefulBCs <- data.frame(matrix(ncol=12, nrow=nrow(BCpairs)))
	colnames(usefulBCs) <- c('A', 'B', '#A', '#B', '#Intersect', '#Union', 'Intersect/Union', paste0('#A_', thresholdRNADNA), paste0('#B_', thresholdRNADNA), paste0('#Intersect_', thresholdRNADNA), paste0('#Union_', thresholdRNADNA), paste0('Intersect_', thresholdRNADNA, '/', 'Union_', thresholdRNADNA))
	
	for (i in 1:nrow(BCpairs)){
		cat(BCpairs$A[i], BCpairs$B[i], '\n')
		
		samA_BC <- fread(paste0(BCpairs$A[i], '/', gsub('.*/', '', BCpairs$A[i]), '.barcode.counts.UMI.corrected.txt'), header=F, stringsAsFactors=F)
		colnames(samA_BC) <- c("bc", "count")
		
		samB_BC <- fread(paste0(BCpairs$B[i], '/', gsub('.*/', '', BCpairs$B[i]), '.barcode.counts.UMI.corrected.txt'), header=F, stringsAsFactors=F)
		colnames(samB_BC) <- c("bc", "count")
		
		
		samA_samB <- intersect(samA_BC$bc, samB_BC$bc)
		usamA_samB <- union(samA_BC$bc, samB_BC$bc)
		
		samA_BC10 <- samA_BC[samA_BC$count >= thresholdRNADNA,]
		samB_BC10 <- samB_BC[samB_BC$count >= thresholdRNADNA,]
		
		samA_samB_10 <- intersect(samA_BC10$bc, samB_BC10$bc)
		usamA_samB_10 <- union(samA_BC10$bc, samB_BC10$bc)
		#Populate data frame
		usefulBCs[,1][i] <- gsub('_', '-', gsub('.*/', '', BCpairs$A[i]))
		usefulBCs[,2][i] <- gsub('_', '-', gsub('.*/', '', BCpairs$B[i]))
		usefulBCs[,3][i] <- nrow(samA_BC)
		usefulBCs[,4][i] <- nrow(samB_BC)
		usefulBCs[,5][i] <- length(samA_samB)
		usefulBCs[,6][i] <- length(usamA_samB)
		usefulBCs[,7][i] <- round(length(samA_samB)/length(usamA_samB), 2)
		usefulBCs[,8][i] <- nrow(samA_BC10)
		usefulBCs[,9][i] <- nrow(samB_BC10)
		usefulBCs[,10][i] <- length(samA_samB_10)
		usefulBCs[,11][i] <- length(usamA_samB_10)
		usefulBCs[,12][i] <- round(length(samA_samB_10)/length(usamA_samB_10), 2)
	}
	for (i in c(3:6, 8:11)){
		usefulBCs[,i] <- format(as.numeric(usefulBCs[,i]), big.mark=",", trim=TRUE)
	}
	write.table(usefulBCs[order(usefulBCs[,1]),], file=opt$output, row.names=F, sep='\t', quote=F)
}



######
#If A, B, and C exists
######
if (!is.null(opt$samplesA) & !is.null(opt$samplesB) & !is.null(opt$samplesC)) {
	samplesB <- strsplit(opt$samplesB, split=',')[[1]]
	cat(samplesB[1], '\n')
	samplesC <- strsplit(opt$samplesC, split=',')[[1]]
	BCpairs <- data.frame("A"=samplesA, "B"=samplesB, "C"=samplesC, stringsAsFactors=F)
	usefulBCs <- data.frame(matrix(ncol=13, nrow=nrow(BCpairs)))
	colnames(usefulBCs) <- c('A', 'B', 'C', '#A', '#B', '#C_mapped', paste0('#A_', thresholdRNADNA), paste0('#B_', thresholdRNADNA), paste0('#C_', thresholdIPCR), paste0('#Intersect_A', thresholdRNADNA, '_B', thresholdRNADNA), paste0('#Intersect_A', thresholdRNADNA, '_C', thresholdIPCR), paste0('#Intersect_B', thresholdRNADNA, '_C', thresholdIPCR), paste0('#Intersect_A', thresholdRNADNA, '_B', thresholdRNADNA, '_C', thresholdIPCR))
	
	for (i in 1:nrow(BCpairs)) {
		samA_BC <- fread(paste0(BCpairs$A[i], '/', gsub('.*/', '', BCpairs$A[i]), '.barcode.counts.UMI.corrected.txt'), header=F, stringsAsFactors=F)
		colnames(samA_BC) <- c("bc", "count")
		
		samB_BC <- fread(paste0(BCpairs$B[i], '/', gsub('.*/', '', BCpairs$B[i]), '.barcode.counts.UMI.corrected.txt'), header=F, stringsAsFactors=F)
		colnames(samB_BC) <- c("bc", "count")
		
		samC <- fread(paste0(BCpairs$C[i], '/', gsub('.*/', '', BCpairs$C[i]), '.barcodes.coords.bed'), header=F, stringsAsFactors=F)
		colnames(samC) <- c("chrom", "chromStart", "chromEnd", "bc", "count", "strand")
		
		samC_BC1 <- samC[samC$count >= thresholdIPCR,]
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
	write.table(usefulBCs[order(usefulBCs[,1]),], file=opt$output, row.names=F, sep='\t', quote=F)
}



#####
#If A and C exists
#####
if (!is.null(opt$samplesA) & is.null(opt$samplesB) & !is.null(opt$samplesC)) {
	samplesB <- strsplit(opt$samplesC, split=',')[[1]]
	cat(samplesB[1], '\n')
	BCpairs <- data.frame("A"=samplesA, "B"=samplesB, stringsAsFactors=F)
	usefulBCs <- data.frame(matrix(ncol=12, nrow=nrow(BCpairs)))
	colnames(usefulBCs) <- c('A', 'B', '#A', '#B', 'Intersect', 'Union', 'Intersect/Union', paste0('#A_', thresholdRNADNA), paste0('#B_', thresholdRNADNA), paste0('Intersect_', thresholdRNADNA), paste0('Union_', thresholdRNADNA), paste0('Intersect_', thresholdRNADNA, '/', 'Union_', thresholdRNADNA))
	
	for (i in 1:nrow(BCpairs)) {
		cat(BCpairs$A[i], BCpairs$B[i], '\n')
		
		samA_BC <- fread(paste0(BCpairs$A[i], '/', gsub('.*/', '', BCpairs$A[i]), '.barcode.counts.UMI.corrected.txt'), header=F, stringsAsFactors=F)
		colnames(samA_BC) <- c("bc", "count")
		
		samC_BC <- fread(paste0(BCpairs$B[i], '/', gsub('.*/', '', BCpairs$B[i]), '.barcodes.coords.bed'), header=F, stringsAsFactors=F)
		colnames(samC_BC) <- c("chrom", "chromStart", "chromEnd", "bc", "count", "strand")
		
		
		samA_samC <- intersect(samA_BC$bc, samC_BC$bc)
		usamA_samC <- union(samA_BC$bc, samC_BC$bc)
		
		samA_BC10 <- samA_BC[samA_BC$count >= thresholdRNADNA,]
		#samC_BC10 <- samC_BC[samC_BC$count >= thresholdRNADNA,]
		
		samA_samC_10 <- intersect(samA_BC10$bc, samC_BC$bc)
		usamA_samC_10 <- union(samA_BC10$bc, samC_BC$bc)
		#Populate data frame
		usefulBCs[,1][i] <- gsub('_', '-', gsub('.*/', '', BCpairs$A[i]))
		usefulBCs[,2][i] <- gsub('_', '-', gsub('.*/', '', BCpairs$B[i]))
		usefulBCs[,3][i] <- nrow(samA_BC)
		usefulBCs[,4][i] <- nrow(samC_BC)
		usefulBCs[,5][i] <- length(samA_samC)
		usefulBCs[,6][i] <- length(usamA_samC)
		usefulBCs[,7][i] <- round(length(samA_samC)/length(usamA_samC), 2)
		usefulBCs[,8][i] <- nrow(samA_BC10)
		usefulBCs[,9][i] <- nrow(samC_BC)
		usefulBCs[,10][i] <- length(samA_samC_10)
		usefulBCs[,11][i] <- length(usamA_samC_10)
		usefulBCs[,12][i] <- round(length(samA_samC_10)/length(usamA_samC_10), 2)
	}
	for (i in c(3:6, 8:11)){
		usefulBCs[,i] <- format(as.numeric(usefulBCs[,i]), big.mark=",", trim=TRUE)
	}
	write.table(usefulBCs[order(usefulBCs[,1]),], file=opt$output, row.names=F, sep='\t', quote=F)
}


#####
#Compare 2 iPCR samples
#####
if (is.null(opt$samplesA) & !is.null(opt$samplesC) & !is.null(opt$samplesC2)) {
	samplesB <- strsplit(opt$samplesC, split=',')[[1]]
	samplesA <- strsplit(opt$samplesC2, split=',')[[1]]
	cat(samplesB[1], '\n')
	BCpairs <- data.frame("A"=samplesA, "B"=samplesB, stringsAsFactors=F)
	usefulBCs <- data.frame(matrix(ncol=7, nrow=nrow(BCpairs)))
	colnames(usefulBCs) <- c('A', 'B', '#A', '#B', 'Intersect', 'Union', 'Intersect/Union')
	
	for (i in 1:nrow(BCpairs)){
		cat(BCpairs$A[i], BCpairs$B[i], '\n')
		
		samC_BC <- fread(paste0(BCpairs$A[i], '/', gsub('.*/', '', BCpairs$A[i]), '.barcodes.coords.bed'), header=F, stringsAsFactors=F)
		colnames(samC_BC) <- c("chrom", "chromStart", "chromEnd", "bc", "count", "strand")
		
		samC2_BC <- fread(paste0(BCpairs$B[i], '/', gsub('.*/', '', BCpairs$B[i]), '.barcodes.coords.bed'), header=F, stringsAsFactors=F)
		colnames(samC2_BC) <- c("chrom", "chromStart", "chromEnd", "bc", "count", "strand")
		
		samA_samB <- intersect(samC_BC$bc, samC2_BC$bc)
		usamA_samB <- union(samC_BC$bc, samC2_BC$bc)
		
		#samC_BC10 <- samC_BC[samC_BC$count >= thresholdRNADNA,]
		#samC2_BC10 <- samC2_BC[samC2_BC$count >= thresholdRNADNA,]
		
		samA_samB_10 <- intersect(samC_BC$bc, samC2_BC$bc)
		usamA_samB_10 <- union(samC_BC$bc, samC2_BC$bc)
		#Populate data frame
		usefulBCs[,1][i] <- gsub('_', '-', gsub('.*/', '', BCpairs$A[i]))
		usefulBCs[,2][i] <- gsub('_', '-', gsub('.*/', '', BCpairs$B[i]))
		usefulBCs[,3][i] <- nrow(samC_BC)
		usefulBCs[,4][i] <- nrow(samC2_BC)
		usefulBCs[,5][i] <- length(samA_samB)
		usefulBCs[,6][i] <- length(usamA_samB)
		usefulBCs[,7][i] <- round(length(samA_samB)/length(usamA_samB), 2)
	}
	for (i in c(3:6)){
		usefulBCs[,i] <- format(as.numeric(usefulBCs[,i]), big.mark=",", trim=TRUE)
	}
	write.table(usefulBCs[order(usefulBCs[,1]),], file=opt$output, row.names=F, sep='\t', quote=F)
}


