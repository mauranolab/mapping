###From ~maagj01/.Rprofile
fancy_scientific <- function(l) { 
	# turn in to character string in scientific notation 
	l <- format(l, scientific = TRUE) 
	# quote the part before the exponent to keep all the digits 
	l <- gsub("^(.*)e", "'\\1'e", l) 
	# turn the 'e+' into plotmath format 
	l <- gsub("e", "%*%10^", l) 
	l <- gsub("\\+",'',l)
	# return this as an expression 
	parse(text=l) 
} 


#Define functions for splitting the lines
splitLines <- function(Term,sep) {
	data.frame(do.call('rbind', strsplit(as.character(Sample[grep(Term,Sample,fixed=T)]),sep,fixed=TRUE)),stringsAsFactors=F)
}

#As above but +1 line
getLength <- function(Term) {
	data.frame(do.call('rbind', strsplit(as.character(Sample[grep(Term,Sample,fixed=T)+1]),' ',fixed=TRUE)),stringsAsFactors=F)
}

getUMIHist <- function(Sample) {
	data.frame(do.call('rbind', strsplit(as.character(Sample[grep('Histogram of number of reads per barcode+UMI',Sample,fixed=T)+1:10]),' ',fixed=TRUE)),stringsAsFactors=F)$X1
}

getBCHist<- function(Sample) {
	data.frame(do.call('rbind', strsplit(as.character(Sample[grep('Histogram of number of reads per barcode',Sample,fixed=T)+1:10]),' ',fixed=TRUE)),stringsAsFactors=F)$X1
}

getInsert <- function(Sample) {
	data.frame(do.call('rbind', strsplit(as.character(Sample[grep('Histogram of number of insertion sites per barcode',Sample,fixed=T)+1:2]),' ',fixed=TRUE)),stringsAsFactors=F)$X1
}

getChrom <- function(Sample) {
	data.frame(do.call('rbind', strsplit(as.character(Sample[grep('Integration sites by chrom',Sample,fixed=T)+1:50]),' ',fixed=TRUE)),stringsAsFactors=F)[,c(1:2)]
}

getReadsPerSite <- function(Sample) {
	data.frame(do.call('rbind', strsplit(as.character(Sample[grep('Histogram of number of insertion sites per barcode',Sample,fixed=T)+4:15]),' ',fixed=TRUE)),stringsAsFactors=F)[,c(1:2)]
}

getGeneLoc <- function(Sample) {
	data.frame(do.call('rbind', strsplit(as.character(Sample[grep('doing GenicLocation',Sample,fixed=T)+1:50]),' ',fixed=TRUE)),stringsAsFactors=F)[,c(1:2)]
}
