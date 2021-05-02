#!/bin/env Rscript
library(tableHTML)
library(dplyr)


###Arguments
if(length(commandArgs(TRUE) == 1)) {
	outdir <- commandArgs(TRUE)[1]
} else {
	stop("ERROR: Missing output directory")
}


###Define functions for splitting the lines
splitLines <- function(Term,sep) {
	data.frame(do.call('rbind', strsplit(as.character(Sample[grep(Term,Sample,fixed=T)]),sep,fixed=TRUE)),stringsAsFactors=F)
}

#As above but +1 line
getLength <- function(Term) {
	data.frame(do.call('rbind', strsplit(as.character(Sample[grep(Term,Sample,fixed=T)+1]),' ',fixed=TRUE)),stringsAsFactors=F)
}

#getUMIHist <- function(Sample) {
#	data.frame(do.call('rbind', strsplit(as.character(Sample[grep('Histogram of number of reads per barcode+UMI',Sample,fixed=T)+1:10]),' ',fixed=TRUE)),stringsAsFactors=F)$X1
#}
#
#getBCHist<- function(Sample) {
#	data.frame(do.call('rbind', strsplit(as.character(Sample[grep('Histogram of number of reads per barcode',Sample,fixed=T)+1:10]),' ',fixed=TRUE)),stringsAsFactors=F)$X1
#}

getInsert <- function(Sample) {
	data.frame(do.call('rbind', strsplit(as.character(Sample[grep('Histogram of number of insertion sites per barcode',Sample,fixed=T)+1:2]),' ',fixed=TRUE)),stringsAsFactors=F)$X1
}

#getChrom <- function(Sample) {
#	data.frame(do.call('rbind', strsplit(as.character(Sample[grep('Integration sites by chrom',Sample,fixed=T)+1:50]),' ',fixed=TRUE)),stringsAsFactors=F)[,c(1:2)]
#}

getReadsPerSite <- function(Sample) {
	data.frame(do.call('rbind', strsplit(as.character(Sample[grep('Histogram of number of insertion sites per barcode',Sample,fixed=T)+4:15]),' ',fixed=TRUE)),stringsAsFactors=F)[,c(1:2)]
}

#getGeneLoc <- function(Sample) {
#	data.frame(do.call('rbind', strsplit(as.character(Sample[grep('doing GenicLocation',Sample,fixed=T)+1:50]),' ',fixed=TRUE)),stringsAsFactors=F)[,c(1:2)]
#}


###
analysisFiles <- list.files('./', pattern='\\.o[0-9]+$', include.dirs = FALSE)
analysisFiles <- analysisFiles[grep('^submit\\.', analysisFiles, invert=T)]
analysisFiles <- analysisFiles[file.info(analysisFiles)$size>200]#Removes i with errors
outputCols <- c("Flowcell", "BS", "Name", "Type", "Total read pairs", "Total barcodes", "BC+UMI", "UMI length", "Unique BC", "BC length", "Complexity", "Mapped reads", "Prop pSB/pTR", "Mapped+BC reads", "Read lengths", "Raw unique sites", "Unique sites", "Prop Sites 2+ BCs", "#BC 1 site", "#BC 2+ sites", "Prop BC 2+ sites", "Prop Reads Chimeric")
data <- as.data.frame(matrix(ncol=length(outputCols), nrow=length(analysisFiles)))
colnames(data) <- outputCols

for (i in 1:length(analysisFiles)) {
	cat(analysisFiles[i], '\n')
	Sample <- readLines(analysisFiles[i])
	Sample <- sub("^\\s+", "", Sample)
	data[i, "Flowcell"] <- basename(getwd())
	data[i, "Name"] <- gsub('_(10x)?RNA$|_DNA$|_iPCR$', '', gsub('\\.o.*|.*-', '', analysisFiles[i]))
	data[i, "BS"] <- gsub('-.*', '', analysisFiles[i])
	if(grepl('Merged', analysisFiles[i])) {
		data[i, "Type"] <- gsub('.*_|\\.o.*', '', gsub('_Merged', '', analysisFiles[i]))
	} else {
		data[i, "Type"] <- gsub('.*_|\\.o.*', '', analysisFiles[i])
	}
	#NB there are some exit routes that print Done!!! but don't print data below
	if(tail(Sample, 2)[1]=="Done!!!") {
		if(!any(grepl('WARNING: no barcodes left, so exiting', Sample))) {
			data[i, "Total read pairs"] <- as.numeric(splitLines('Number of total reads', '\t')$X3)
			data[i, "Total barcodes"] <- as.numeric(splitLines('Number of total read barcodes', '\t')$X3)
			#Note tail() takes "Number of unique barcodes passing minimum read cutoff" when present; "Number of unique barcodes" if not
			#Note we grep for "Number of unique barcodes" but actually use tail to get "Number of unique barcodes passing minimum read cutoff" for recent versions of pipeline -- TODO just get the latter directly
			data[i, "Unique BC"] <- as.numeric(tail(splitLines('Number of unique barcodes passing', '\t'),1)$X3)
			data[i, "BC length"] <- getLength('Barcode lengths')$X2
		}
		if(!any(grepl('No UMIs found', Sample))) {
			data[i, "BC+UMI"] <- as.numeric(splitLines('Number of unique barcodes+UMI', '\t')$X3)
			if(!is.na(getLength('UMI lengths')$X2)) {
				data[i, "UMI length"] <- getLength('UMI lengths')$X2
			} else if(getLength('UMI lengths')$X1!=data[i, "BC+UMI"]) {
				data[i, "Total barcodes"] <- strsplit(Sample[grep(data[i, "Total read pairs"], Sample, fixed=T)][2], "\\s+")[[1]][2]
				if(is.na(strsplit(Sample[grep(data[i, "BC+UMI"], Sample, fixed=T)][2], "\\s+")[[1]][2]) && length(Sample[grep(data[i, "BC+UMI"]-1, Sample, fixed=T)])==0) {
					data[i, "UMI length"] <- NA
				} else {
					data[i, "UMI length"] <- strsplit(Sample[grep(data[i, "BC+UMI"]-1, Sample, fixed=T)], "\\s+")[[1]][2]
				}
			} else {
				data[i, "UMI length"] <- getLength('UMI lengths')$X2
			}
			data[i, "UMI length"] <- as.numeric(data[i, "UMI length"])
			data[i, "Complexity"] <- signif(data[i, "BC+UMI"] / data[i, "Total barcodes"], 2)
		}
		
		if(data[i, "Type"] == "iPCR") {
			data[i, "Mapped reads"] <- as.numeric(splitLines('Total PF reads', '\t')$X3)
			data[i, "Prop pSB/pTR"] <- signif(as.numeric(splitLines('Total reads mapping to pSB or pTR', '\t')$X3) / data[i, "Mapped reads"], 2)
			data[i, "Mapped+BC reads"] <- as.numeric(splitLines('Number of reads passing all filters and having barcodes assigned', '\t')$X3)
			data[i, "Read lengths"] <- gsub(" ", "", gsub(" \\([0-9]+\\)", "", splitLines('Read lengths', '\t')$X3))
			data[i, "Read lengths"] <- paste(min(unlist(strsplit(data[i, "Read lengths"], ","))), max(unlist(strsplit(data[i, "Read lengths"], ","))), sep="-")
			data[i, "Raw unique sites"] <- as.numeric(splitLines('Number of unique insertion sites before collapsing nearby ones', '\t')$X3[1])
			#Note we use tail to get "Number of unique insertion sites" rather than the "before collapsing nearby ones" one
			data[i, "Unique sites"] <- as.numeric(tail(splitLines('Number of unique insertion sites', '\t'), 1)$X3[1])
			#BUGBUG value in numerator depends on #BCs at site, e.g. 3 BCs is more than 2
			data[i, "Prop Sites 2+ BCs"] <- signif(1 - data[i, "Unique sites"] / as.numeric(splitLines('Number of BC+insertions after minimum read cutoff', '\t')$X3[1]), 2)
			data[i, "#BC 1 site"] <- as.numeric(getInsert(Sample)[1])
			data[i, "#BC 2+ sites"] <- as.numeric(getInsert(Sample)[2])
			data[i, "Prop BC 2+ sites"] <- signif(data[i, "#BC 2+ sites"] / (data[i, "#BC 1 site" ] + data[i, "#BC 2+ sites"]), 2)
			data[i, "Prop. BC+insertions chimeric"] <- signif(1 - as.numeric(splitLines('Number of BC+insertions passing minPropReadsAtSite cutoff', '\t')$X3[1]) / as.numeric(splitLines('Number of BC+insertions after collapsing nearby ones', '\t')$X3[1]), 2)
		} else if(data[i, "Type"] == "10xRNA") {
			data[i, "Prop Reads Chimeric"] <- signif(1 - as.numeric(splitLines('Number of reads passing minPropReadsForBC filter', '\t')$X3[1]) / as.numeric(splitLines('Number of reads passing minPropReadsForBC filter', '\t')$X4[1]), 2)
			print(t(splitLines('Number of reads passing minPropReadsForBC filter', '\t')))
		}
	}
}


#Save a tsv without commas in numeric columns
write.table(data, file = paste0(outdir, "/FlowcellSummary/summaryflowcell.tsv"), sep='\t', quote=F, col.names=T, row.names=F)


#Format for printing
for (i in 1:ncol(data)) {
	if(is.numeric(data[,i])) {
		data[,i] <- format(data[,i], big.mark=",", trim=TRUE)
	}
}


SumHTMLtable <- tableHTML(data, rownames=F, footer='N.B.: Unique BC and integration site counts require 10+ reads for DNA/RNA, and 2+ reads for iPCR') %>% add_css_row(css = list('background-color', 'lightblue'), rows=odd(1:nrow(data)))
write_tableHTML(SumHTMLtable, file = paste0(outdir, "/FlowcellSummary/index.html"))
