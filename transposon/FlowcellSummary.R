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
outputCols <- c("Flowcell", "BS", "Name", "sampleType", "Total reads", "Total barcodes", "BC+UMI", "UMI length", "Unique BC", "BC length", "Complexity", "PF reads", "QC mapped reads", "Prop PF", "Unique sites", "Unique sites +-5bp", "#BC 1 site", "#BC 2+ sites", "Prop 2+")
sumBarcode <- as.data.frame(matrix(ncol=length(outputCols), nrow=length(Barcodes)))
colnames(sumBarcode) <- outputCols

for (files in 1:length(Barcodes)) {
	cat(Barcodes[files], '\n')
	Sample <- readLines(Barcodes[files])
	Sample <- sub("^\\s+", "", Sample)
	if(tail(Sample, 2)[1]=="Done!!!") {
		sumBarcode[files, "BS"] <- gsub('-.*', '', Barcodes[files])
		sumBarcode[files, "Name"] <- gsub('_RNA$|_DNA$|_iPCR$', '', gsub('\\.o.*|.*-', '', Barcodes[files]))
		if(grepl('Merged', Barcodes[files])) {
			sumBarcode[files,]$sampleType <- gsub('.*_|\\.o.*', '', gsub('_Merged', '', Barcodes[files]))
		} else {
			sumBarcode[files,]$sampleType <- gsub('.*_|\\.o.*', '', Barcodes[files])
		}
		sumBarcode[files, "Flowcell"] <- basename(getwd())
		sumBarcode[files, "Total reads"] <- as.numeric(splitLines('Number of total reads', '\t')$X3)
		sumBarcode[files, "Total barcodes"] <- as.numeric(splitLines('Number of total read barcodes', '\t')$X3)
		sumBarcode[files, "Unique BC"] <- as.numeric(tail(splitLines('Number of unique barcodes', '\t'),1)$X3)
		if(getLength('Barcode lengths')$X1!=sumBarcode[files, "Unique BC"]) {
			sumBarcode[files, "BC length"] <- strsplit(Sample[grep(as.character(as.numeric(sumBarcode[files, "Unique BC"])-1), Sample, fixed=T)], "\\s+")[[1]][2]
		} else {
			sumBarcode[files, "BC length"] <- getLength('Barcode lengths')$X2
		}
		if(!any(grepl('No UMIs found', Sample))) {
			sumBarcode[files, "BC+UMI"] <- as.numeric(splitLines('Number of unique barcodes+UMI', '\t')$X3)
			if(!is.na(getLength('UMI lengths')$X2)) {
				sumBarcode[files, "UMI length"] <- getLength('UMI lengths')$X2
			} else if(getLength('UMI lengths')$X1!=sumBarcode[files, "BC+UMI"]) {
				sumBarcode[files, "Total barcodes"] <- strsplit(Sample[grep(sumBarcode[files, "Total reads"], Sample, fixed=T)][2], "\\s+")[[1]][2]
				if(is.na(strsplit(Sample[grep(sumBarcode[files, "BC+UMI"], Sample, fixed=T)][2], "\\s+")[[1]][2]) && length(Sample[grep(sumBarcode[files, "BC+UMI"]-1, Sample, fixed=T)])==0) {
					sumBarcode[files, "UMI length"] <- NA
				} else {
					sumBarcode[files, "UMI length"] <- strsplit(Sample[grep(sumBarcode[files, "BC+UMI"]-1, Sample, fixed=T)], "\\s+")[[1]][2]
				}
			} else {
				sumBarcode[files, "UMI length"] <- getLength('UMI lengths')$X2
			}
			sumBarcode[files, "UMI length"] <- as.numeric(sumBarcode[files, "UMI length"])
			sumBarcode[files, "Complexity"] <- sumBarcode[files, "Unique BC"] / sumBarcode[files, "BC+UMI"]
		}
		
		if(any(grepl('Total PF tags', Sample))) {
			sumBarcode[files, "PF reads"] <- as.numeric(splitLines('Total PF tags', '\t')$X3)
			sumBarcode[files, "QC mapped reads"] <- as.numeric(splitLines('Number of tags passing all filters and having barcodes assigned', '\t')$X2)
			sumBarcode[files, "Prop PF"] <- sumBarcode[files, "QC mapped reads"] / (sumBarcode[files, "PF reads"] + sumBarcode[files, "QC mapped reads"])
			sumBarcode[files, "Unique sites"] <- as.numeric(splitLines('Total uniq sites', '\t')$X3[1])
			sumBarcode[files, "Unique sites +-5bp"] <- as.numeric(splitLines('Total uniq sites (within 5 bp, ignoring strand)', '\t')$X3[1])
			sumBarcode[files, "#BC 1 site"] <- as.numeric(getInsert(Sample)[1])
			sumBarcode[files, "#BC 2+ sites"] <- as.numeric(getInsert(Sample)[2])
			sumBarcode[files, "Prop 2+"] <- sumBarcode[files, "#BC 2+ sites"] / (sumBarcode[files, "#BC 1 site" ] + sumBarcode[files, "#BC 2+ sites"])
		}
	}
}


#Format for printing
for (i in 5:ncol(sumBarcode)) {
	sumBarcode[,i] <- format(sumBarcode[,i], big.mark=",", trim=TRUE)
}


SumHTMLtable <- tableHTML(sumBarcode) %>%  add_css_row(css = list('background-color', 'lightblue'), rows = odd(1:nrow(sumBarcode)))
write_tableHTML(SumHTMLtable, file = paste0(outdir, "/FlowcellSummary/index.html"))
write.table(sumBarcode, file = paste0(outdir, "/FlowcellSummary/summaryflowcell.tsv"), sep='\t', quote=F, col.names=T, row.names=F)
