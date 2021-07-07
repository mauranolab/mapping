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
samplesheetfile="/vol/mauranolab/flowcells/data/${fc}/SampleSheet.csv"
Rscript --quiet --no-save - ${samplesheetfile} << 'EOF'
samplesheetfile <- commandArgs(TRUE)[1]

duplicaterows <- function(x, copies=2) {
	x.rle <- rle(x)
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
