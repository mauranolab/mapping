#!/bin/env Rscript
source("/vol/mauranolab/transposon/src/maagj01_profile.R")

library(tableHTML)
library(dplyr)

if(length(commandArgs(TRUE) == 1)) {
	outdir <- commandArgs(TRUE)[1]
} else {
	stop("ERROR: Missing output directory")
}


Barcodes <- list.files('./', pattern='^BS', include.dirs = FALSE)
Barcodes <- Barcodes[file.info(Barcodes)$size>200]#Removes files with errors
Barcodes <- Barcodes[grep('\\.o[0-9]', Barcodes)]
Barcodes <- Barcodes[file.info(Barcodes)$size>200]#Removes files with errors
sumBarcode <- as.data.frame(matrix(ncol=19, nrow=length(Barcodes)))
colnames(sumBarcode) <- c('Flowcell', 'BSnumber', 'Name', 'sampleType', 'Total reads', 'Total barcodes', 'BC+UMI', 'UMI length', 'Unique BC', 'BC length', 'Complexity', 'PF reads', 'QC mapped reads', 'Prop PF', 'Unique sites', 'Unique sites +-5bp', '#BC 1 site', '#BC 2+ sites', 'Prop 2+')

for (files in 1:length(Barcodes)) { cat(Barcodes[files], '\n')
	Sample <- readLines(Barcodes[files])
	Sample <- sub("^\\s+", "", Sample)
	if (tail(Sample, 2)[1]=="Done!!!") {
	sumBarcode[files,]$BSnumber <- gsub('-.*', '', Barcodes[files])
	sumBarcode[files,]$Name <- gsub('_RNA$|_DNA$|_iPCR$', '', gsub('\\.o.*|.*-', '', Barcodes[files]))
	if(length(grep('Merged', Barcodes[files]))==1){ sumBarcode[files,]$sampleType <- gsub('.*_|\\.o.*', '', gsub('_Merged', '', Barcodes[files])) }
	else { sumBarcode[files,]$sampleType <- gsub('.*_|\\.o.*', '', Barcodes[files]) }
	sumBarcode[files,]$Flowcell <- gsub('.*/', '', getwd())
	sumBarcode[files,][,5] <- splitLines('Number of total reads', '\t')$X3
	sumBarcode[files,][,6] <- splitLines('Number of total read barcodes', '\t')$X3
	if (length(grep('No UMIs found', Sample))==0){
		sumBarcode[files,][,7] <- splitLines('Number of unique barcodes+UMI', '\t')$X3
		sumBarcode[files,][,9] <- splitLines('Number of unique barcodes', '\t')[2,]$X3
		if(!is.na(getLength('UMI lengths')$X2)) { sumBarcode[files,][,8] <- getLength('UMI lengths')$X2 }
		else if (getLength('UMI lengths')$X1!=sumBarcode[files,][,7]){
			sumBarcode[files,][,6] <- strsplit(Sample[grep(sumBarcode[files,][,5], Sample, fixed=T)][2], "\\s+")[[1]][2]
			if (is.na(strsplit(Sample[grep(sumBarcode[files,][,7], Sample, fixed=T)][2], "\\s+")[[1]][2]) && length(Sample[grep(as.numeric(sumBarcode[files,][,7])-1, Sample, fixed=T)])==0){
			sumBarcode[files,][,8] <- NA
			} else sumBarcode[files,][,8] <- strsplit(Sample[grep(as.numeric(sumBarcode[files,][,7])-1, Sample, fixed=T)], "\\s+")[[1]][2]
		} else sumBarcode[files,][,8] <- getLength('UMI lengths')$X2
	} else if (length(splitLines('Number of unique barcodes', '\t')$X3)==2){ sumBarcode[files,][,8] <- NA
	} else sumBarcode[files,][,9] <- splitLines('Number of unique barcodes', '\t')$X3
	if (getLength('Barcode lengths')$X1!=sumBarcode[files,][,9]){
		sumBarcode[files,][,10] <- strsplit(Sample[grep(as.character(as.numeric(sumBarcode[files,][,9])-1), Sample, fixed=T)], "\\s+")[[1]][2]
	} else sumBarcode[files,][,10] <- getLength('Barcode lengths')$X2
	
	 if (length(grep('Total PF tags', Sample))==1){ 
		 sumBarcode[files,][,12] <- splitLines('Total PF tags', '\t')$X3
		 sumBarcode[files,][,13] <- splitLines('Number of tags passing all filters and having barcodes assigned', '\t')$X2
		 sumBarcode[files,][,15] <- splitLines('Total uniq sites', '\t')$X3[1]
		 sumBarcode[files,][,16] <- splitLines('Total uniq sites (within 5 bp, ignoring strand)', '\t')$X3[1]
		 sumBarcode[files,][17] <- getInsert(Sample)[1]
		 sumBarcode[files,][18] <- getInsert(Sample)[2]
	 
	 }
	} else {
		sumBarcode[files,]$Flowcell <- gsub('.*/', '', getwd())
		sumBarcode[files,]$BSnumber <- gsub('-.*', '', Barcodes[files])
		sumBarcode[files,]$Name <- gsub('_RNA$|_DNA$|_iPCR$', '', gsub('\\.o.*|.*-', '', Barcodes[files]))
		sumBarcode[files,]$sampleType <- gsub('.*_|\\.o.*', '', Barcodes[files])
		sumBarcode[files,][,5:18] <- NA
	 }
}

SummarizeFlowcell <- sumBarcode
SummarizeFlowcell$Complexity <- 0
for (i in 1:nrow(SummarizeFlowcell)){
	#Complexity
	if (is.na(SummarizeFlowcell[,7][i])) { SummarizeFlowcell$Complexity[i] <- 1
	} else SummarizeFlowcell$Complexity[i] <- as.numeric(SummarizeFlowcell[,9][i])/as.numeric(SummarizeFlowcell[,7][i])
	cat(as.numeric(SummarizeFlowcell[,9][i])/as.numeric(SummarizeFlowcell[,7][i]),'\n')
	
	#Prop PF
	if(!is.na(SummarizeFlowcell[,12][i])) { SummarizeFlowcell[,14][i] <- as.numeric(SummarizeFlowcell[,13][i])/(as.numeric(SummarizeFlowcell[,12][i])+as.numeric(SummarizeFlowcell[,13][i])) }
	if(!is.na(SummarizeFlowcell[,12][i])) { SummarizeFlowcell[,19][i] <- as.numeric(SummarizeFlowcell[,18][i])/(as.numeric(SummarizeFlowcell[,17][i])+as.numeric(SummarizeFlowcell[,18][i])) }
}



for (i in 5:ncol(SummarizeFlowcell)) {
    SummarizeFlowcell[,i] <- format(as.numeric(SummarizeFlowcell[,i]), big.mark=",", trim=TRUE)
}  
library(tableHTML)
SumHTMLtable <- tableHTML(SummarizeFlowcell) %>%  add_css_row(css = list('background-color', 'lightblue'), rows = odd(1:nrow(SummarizeFlowcell)))
write_tableHTML(SumHTMLtable, file = paste0(outdir, "/FlowcellSummary/index.html"))
write.table(SummarizeFlowcell, file = paste0(outdir, "/FlowcellSummary/summaryflowcell.tsv"), sep='\t', quote=F, col.names=T, row.names=F)
