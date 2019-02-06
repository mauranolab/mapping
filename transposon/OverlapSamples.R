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
option_list = list(
	make_option(c("-a", "--samplesA"), type="character", default=NULL, 
		help="Comma separated list of SamplesA", metavar="character"),
	make_option(c("-b", "--samplesB"), type="character", default=NULL, 
		help="Comma separated list of SamplesB", metavar="character"),
	make_option(c("-c", "--samplesC"), type="character", default=NULL, 
		help="Comma separated list of SamplesC (assumed to be iPCR)", metavar="character"),
	make_option(c("--typeA"), type="character", default=NULL, 
		help="Sample type for sample A (for thresholding and inferring filename): DNA, RNA, iPCR", metavar="character"),
	make_option(c("--typeB"), type="character", default=NULL, 
		help="Sample type for sample B (for thresholding and inferring filename): DNA, RNA, iPCR", metavar="character"),
	make_option(c("--typeC"), type="character", default=NULL, 
		help="Sample type for sample C (for thresholding and inferring filename): DNA, RNA, iPCR", metavar="character"),
	make_option(c("-o", "--output"), type="character", default=NULL, 
		help="output folder for all comparisons", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if(is.null(opt$samplesA)){
	print_help(opt_parser)
	stop("Samples A must be provided", call.=FALSE)
}
if(is.null(opt$typeA)){
	print_help(opt_parser)
	stop("Samples A type must be provided", call.=FALSE)
}
samplesA <- strsplit(opt$samplesA, split=',')[[1]]
message(samplesA[1], '\n')

if(!is.null(opt$samplesB)){
	message("Samples B are provided")
	
	if(is.null(opt$typeB)){
		print_help(opt_parser)
		stop("Samples B type must be provided", call.=FALSE)
	}
	
	samplesB <- strsplit(opt$samplesB, split=',')[[1]]
	message(samplesB[1], '\n')
}

if(!is.null(opt$samplesC)){
	message("Samples C are provided")
	if(is.null(opt$typeC)){
		print_help(opt_parser)
		stop("Samples C type must be provided", call.=FALSE)
	}
	if(is.null(opt$samplesB)){
		stop("Samples B must be provided", call.=FALSE)
	}
	samplesC <- strsplit(opt$samplesC, split=',')[[1]]
	message(samplesC[1], '\n')
}
#TODO verify --type options


if(is.null(opt$output)){
	stop("Output file needs to be provided", call.=FALSE)
}



#
#samplesA <- c("aligned.FCHVTCCBGX2/BS00370A-MSH_060517_K562~CMV~2roundLinear_Rep1_RNA", 'aligned.FCHVTCCBGX2/BS00372A-MSH_060517_K562~CMV~1roundLinear_Rep1_RNA', 'aligned.FCHVTCCBGX2/BS00372A-MSH_060517_K562~CMV~1roundLinear_Rep1_RNA')
#samplesB <- c('aligned.FCHVTCCBGX2/BS00371A-MSH_060517_K562~CMV~2roundLinear_Rep2_RNA', 'aligned.FCHVTCCBGX2/BS00373A-MSH_060517_K562~CMV~1roundLinear_Rep2_RNA', "aligned.FCHVTCCBGX2/BS00375A-MSH_060517_K562~GGloHS2~2roundLinear_Rep2_RNA")
#samplesC <- c('aligned.Merged/BS00244-K562_CMVBC4_iPCR_Merged', 'aligned.Merged/BS00244-K562_CMVBC4_iPCR_Merged', 'aligned.Merged/BS00244-K562_CMVBC4_iPCR_Merged')


#samplesA <- c(samplesA, samplesB)
#samplesA <- c("aligned.FCH7HTKBGX3/BS00490A-MSH_K562~A1~GGlo~HS2~A1~1Linear~NoExo_Rep1_DNA, aligned.FCH7HTKBGX3/BS00491A-MSH_K562~A1~GGlo~HS2~A1~1Linear~NoExo_Rep2_DNA, aligned.FCH7HTKBGX3/BS00492A-MSH_K562~A1~GGlo~HS2~A1~1Linear~ExoI_Rep1_DNA, aligned.FCH7HTKBGX3/BS00493A-MSH_K562~A1~GGlo~HS2~A1~1Linear~ExoI_Rep2_DNA, aligned.Merged/BS00467_468A--MSH-K562~fw_GGlo_fw_OldProtocol-DNA_Merged, aligned.Merged/BS00217-K562d11_A1fw_GGlo_A1fw_HS2_BC4_1000ng_25Ksort_GTG_DNA_Merged")
#samplesA <-strsplit(samplesA, split=',')[[1]]


loadBCfile <- function(filename, type, thresh) {
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
	
	data <- data[data$count >= thresh,]
	
	return(data)
}

thresholdRNADNA <- 10
thresholdiPCR = 1

getThreshold <- function(type, applyThreshold) {
	if(applyThreshold) {
		if(type=="iPCR") {
			thresh=thresholdiPCR
		} else {
			thresh=thresholdRNADNA
		}
	} else {
		thresh=0
	}
	return(thresh)
}

#TODO no arg handling
comparePair <- function(BCpairs, typeA, typeB, applyThreshold=FALSE) {
	threshA = getThreshold(typeA, applyThreshold)
	threshB = getThreshold(typeB, applyThreshold)
	
	message("\nComparing ", typeA, " vs. ", typeB, ". Thresholds=", threshA, "/", threshB)
	
	outputCols <- c("Sample", "BS_A", "BS_B", paste0("#A_", typeA, "_", threshA), paste0("#B_", typeB, "_", threshB), paste0("#Intersect_A", threshA, "_B", threshB), paste0("#Union_A", threshA, "_B", threshB), paste0("Intersect/Union_A", threshA, "_B", threshB))
	usefulBCs <- data.frame(matrix(ncol=length(outputCols), nrow=nrow(BCpairs)))
	colnames(usefulBCs) <- outputCols
	
	for (i in 1:nrow(BCpairs)){
		sampleA <- basename(BCpairs$A[i])
		sampleB <- basename(BCpairs$B[i])
		
		#TODO need better sample name than gsub with typeA
		sample <- gsub(paste0('_', typeA), '', unlist(strsplit(sampleA, split='-'))[2])
		bsA <- unlist(strsplit(sampleA, split='-'))[1]
		bsB <- unlist(strsplit(sampleB, split='-'))[1]
		
		message(paste(sample, bsA, bsB, sep=" "))
		
		sampleA_BC <- loadBCfile(BCpairs$A[i], type=typeA, threshA)
		sampleB_BC <- loadBCfile(BCpairs$B[i], type=typeB, threshB)
		
		sampleA_sampleB <- intersect(sampleA_BC$bc, sampleB_BC$bc)
		usampleA_sampleB <- union(sampleA_BC$bc, sampleB_BC$bc)
		
		usefulBCs[,1][i] <- sample
		usefulBCs[,2][i] <- bsA
		usefulBCs[,3][i] <- bsB
		usefulBCs[,4][i] <- format(as.numeric(nrow(sampleA_BC)), big.mark=",", trim=TRUE)
		usefulBCs[,5][i] <- format(as.numeric(nrow(sampleB_BC)), big.mark=",", trim=TRUE)
		usefulBCs[,6][i] <- format(as.numeric(length(sampleA_sampleB)), big.mark=",", trim=TRUE)
		usefulBCs[,7][i] <- format(as.numeric(length(usampleA_sampleB)), big.mark=",", trim=TRUE)
		usefulBCs[,8][i] <- format(as.numeric(round(length(sampleA_sampleB)/length(usampleA_sampleB), 2)), big.mark=",", trim=TRUE)
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
	usefulBCs <- comparePair(BCpairs, opt$typeA, opt$typeA, applyThreshold=F)
	usefulBCs <- merge(usefulBCs, comparePair(BCpairs, opt$typeA, opt$typeA, applyThreshold=T), by=c("Sample", "BS_A", "BS_B"))
	write.table(usefulBCs[order(usefulBCs[,1]),], file=opt$output, row.names=F, sep='\t', quote=F)
}


#####
#If A and B exist
#####
if (!is.null(opt$samplesA) & !is.null(opt$samplesB) & is.null(opt$samplesC)) {
	BCpairs <- data.frame("A"=samplesA, "B"=samplesB, stringsAsFactors=F)
	usefulBCs <- comparePair(BCpairs, opt$typeA, opt$typeB, applyThreshold=F)
	usefulBCs <- merge(usefulBCs, comparePair(BCpairs, opt$typeA, opt$typeB, applyThreshold=T), by=c("Sample", "BS_A", "BS_B"))
	write.table(usefulBCs[order(usefulBCs[,1]),], file=opt$output, row.names=F, sep='\t', quote=F)
}



######
#If A, B, and C exists
######
if (!is.null(opt$samplesA) & !is.null(opt$samplesB) & !is.null(opt$samplesC)) {
	BCpairs <- data.frame("A"=samplesA, "B"=samplesB, "C"=samplesC, stringsAsFactors=F)
	
	typeA = opt$typeA
	typeB = opt$typeB
	typeC = opt$typeC
	
	threshA = getThreshold(typeA, applyThreshold=TRUE)
	threshB = getThreshold(typeB, applyThreshold=TRUE)
	threshC = getThreshold(typeC, applyThreshold=TRUE)
	
	outputCols <- c("Sample", "BS_A", "BS_B", "BS_C", paste0('#A_', threshA), paste0('#B_', threshB), paste0('#C_', threshC), paste0('#Intersect_A', threshA, '_B', threshB), paste0('#Intersect_A', threshA, '_C', threshC), paste0('#Intersect_B', threshB, '_C', threshC), paste0('#Intersect_A', threshA, '_B', threshB, '_C', threshC), paste0('#Union_A', threshA, '_B', threshB, '_C', threshC))
	usefulBCs <- data.frame(matrix(ncol=length(outputCols), nrow=nrow(BCpairs)))
	colnames(usefulBCs) <- outputCols
	
	message("\nComparing ", typeA, " vs. ", typeB, " vs. ", typeC, ". Thresholds=", threshA, "/", threshB, "/", threshC)
	
	for (i in 1:nrow(BCpairs)) {
		sampleA <- basename(BCpairs$A[i])
		sampleB <- basename(BCpairs$B[i])
		sampleC <- basename(BCpairs$C[i])
		
		#TODO need better sample name than gsub with typeA
		sample <- gsub(paste0('_', typeA), '', unlist(strsplit(sampleA, split='-'))[2])
		bsA <- unlist(strsplit(sampleA, split='-'))[1]
		bsB <- unlist(strsplit(sampleB, split='-'))[1]
		bsC <- unlist(strsplit(sampleC, split='-'))[1]
		
		message(paste(sample, bsA, bsB, bsC, sep=" "))
		
		sampleA_BC <- loadBCfile(BCpairs$A[i], type=typeA, thresh=threshA)
		sampleB_BC <- loadBCfile(BCpairs$B[i], type=typeB, thresh=threshB)
		sampleC_BC <- loadBCfile(BCpairs$C[i], type=typeC, thresh=threshC)
		
		sampleA_sampleB <- intersect(sampleA_BC$bc, sampleB_BC$bc)
		sampleA_sampleC <- intersect(sampleA_BC$bc, sampleC_BC$bc)
		sampleB_sampleC <- intersect(sampleC_BC$bc, sampleB_BC$bc)
		sampleB_sampleC_sampleA <- intersect(sampleB_sampleC, sampleA_BC$bc)
		usampleB_sampleC_sampleA <- union(sampleB_sampleC, sampleA_BC$bc)
		
		usefulBCs[,1][i] <- sample
		usefulBCs[,2][i] <- bsA
		usefulBCs[,3][i] <- bsB
		usefulBCs[,4][i] <- bsC
		usefulBCs[,5][i] <- format(as.numeric(nrow(sampleA_BC)), big.mark=",", trim=TRUE)
		usefulBCs[,6][i] <- format(as.numeric(nrow(sampleB_BC)), big.mark=",", trim=TRUE)
		usefulBCs[,7][i] <- format(as.numeric(nrow(sampleC_BC)), big.mark=",", trim=TRUE)
		usefulBCs[,8][i] <- format(as.numeric(length(sampleA_sampleB)), big.mark=",", trim=TRUE)
		usefulBCs[,9][i] <- format(as.numeric(length(sampleA_sampleC)), big.mark=",", trim=TRUE)
		usefulBCs[,10][i] <- format(as.numeric(length(sampleB_sampleC)), big.mark=",", trim=TRUE)
		usefulBCs[,11][i] <- format(as.numeric(length(sampleB_sampleC_sampleA)), big.mark=",", trim=TRUE)
		usefulBCs[,12][i] <- format(as.numeric(length(usampleB_sampleC_sampleA)), big.mark=",", trim=TRUE)
	}
	write.table(usefulBCs[order(usefulBCs[,1]),], file=opt$output, row.names=F, sep='\t', quote=F)
}
