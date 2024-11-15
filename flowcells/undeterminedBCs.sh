#!/bin/bash
#BUGBUG I think head below causes nonzero exit code from zcat, preventing using pipefail
#set -eu -o pipefail
set -u

fc=$1

echo "Postprocessing ${fc}"


echo
date
echo "Undetermined BCs"

zcat Undetermined_S0_R1_001.fastq.gz | awk -F ":" 'BEGIN {OFS="\t"} NR % 4 ==1 {split($10, bc, "+"); print bc[1], bc[2]}' | head -10000000 | sort -S 6G -g | uniq -c | sort -k1,1g |
awk -v nreads=10000000 'BEGIN {OFS="\t"; print "#reads", "Prop._of_undet_reads", "BC1", "BC2"} $1/nreads>0.001 {print $1, $1/nreads, $2, $3}' | 
tee undetermined_barcodes.txt |
#Remove unlikely barcodes probably from sequencing errors
awk -F "\t" 'BEGIN {OFS="\t"} $3!~/^[GN]+$/ && $4!~/^[GN]+$/'

#
##to debug undetermined barcodes
#zcat Undetermined_S0_R1_001.fastq.gz | awk -F ":" 'BEGIN {OFS="\t"} NR % 4 ==1 {split($10, bc, "+");  bc1=bc[1]; bc2=bc[2]} NR % 4 ==2 && bc2=="TAAGGC"' |m


echo
echo "Possible index hopping in undetermined reads"
samplesheetfile="/gpfs/data/isg_sequencing/data/${fc}/SampleSheet.csv"
Rscript --quiet --no-save - ${samplesheetfile} << 'EOF'
#read() -- smarter wrapper around read.table
#BUGBUG who knows why i'm getting "zcat: stdout: Broken pipe" on first read.table
#TODO seems I need to use gzcat on the mac and zcat on linux
read <- function(filename, nrows=100, header=F, col.classes=NULL, stringsAsFactors=F, sep="\t", comment.char="", quote="", ...) {
	#BUGBUG should print example line of file upon failure
	#TODO can't use a pipe or specify which columns to use
	#BUGBUG when the first column is all characters, R tries to use it as row names. But this goofs up my hardcoded colClasses. row.names=NULL doesn't work in the function, though it works interactively.
	
	filename <- Sys.glob(filename)
	if(length(filename) == 0 || !file.exists(filename)) {
		stop("ERROR read() -- ", filename, " does not exist!\n")
	} else if(length(filename) > 1) {
		stop("ERROR read() -- Multiple files were matched: ",  paste(filename, collapse=", "), "!\n")
	}
	#TODO seems I need to use gzcat on the mac and zcat on linux
	#readCmd <- paste("gzcat -S \"\" -f \"", filename, "\"", sep="")
	readCmd <- paste("zcat -f \"", filename, "\"", sep="")
#	cat("reading with:", readCmd, "\n")

	#Let UNIX tell us how much memory to allocate
	wcpipe <- pipe(paste(readCmd, " | wc -l | awk '{print $1}' ", sep=""))
	num.lines <- as.integer(readLines(wcpipe))
	close(wcpipe)
	
	if(num.lines==0) {
		message("WARNING read() -- input file had zero lines")
		return(NULL)
	} else {
		if(is.null(col.classes)) {
			#Let R estimate the column classes so we can hardcode it later
			#Hijack nrows parameter for the number of rows upon which to estimate classes (the real number of rows is used for the final read.table)
		
			shortcopy <- read.table(pipe(readCmd), sep=sep, comment.char = comment.char, quote = quote, strip.white = TRUE, stringsAsFactors = F, nrows=nrows, header=header, ...)
			guessed.classes <- lapply(shortcopy, typeof)
		} else {
			guessed.classes <- col.classes
		}
		
		return(read.table(pipe(readCmd), sep=sep, comment.char = comment.char, quote = quote, strip.white = TRUE, stringsAsFactors = stringsAsFactors, nrows=num.lines, colClasses=guessed.classes, header=header, ...))
	}
}

samplesheetfile <- commandArgs(TRUE)[1]

duplicaterows <- function(x, copies=2) {
	x.rle <- rle(sort(x))
	dupls <- x.rle$values[x.rle$lengths >= copies]
	return(x %in% dupls)
}


undetermined <- read("undetermined_barcodes.txt", header=T)


#TODO not lane aware
data <- read.table(samplesheetfile, header=T, skip=19, sep=",", comment.char = "", quote = "", strip.white = TRUE, stringsAsFactors = F)


undetermined$BC1.index <- sapply(undetermined$BC1, FUN=function(x) {paste(which(x==data$index), collapse=",")})
undetermined$BC2.index <- sapply(undetermined$BC2, FUN=function(x) {paste(which(x==data$index2), collapse=",")})

print.data.frame(subset(undetermined, BC1.index!="" & BC2.index!=""), row.names=F)
EOF


echo
echo "Done!!!"
date
