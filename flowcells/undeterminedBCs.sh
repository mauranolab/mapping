#!/bin/bash
#set -eu -o pipefail
set -u

echo
date
echo "Undetermined BCs"

#BUGBUG I think head causes nonzero exit code from zcat
zcat Undetermined_S0_R1_001.fastq.gz | awk -F ":" 'BEGIN {OFS="\t"} NR % 4 ==1 {split($10, bc, "+"); print bc[1], bc[2]}' | head -10000000 | sort -S 6G -g | uniq -c | sort -k1,1g |
awk -v nreads=10000000 'BEGIN {OFS="\t"; print "#reads", "Prop._of_undet_reads", "BC1", "BC2"} $1/nreads>0.001 {print $1, $1/nreads, $2, $3}' | 
tee unexplained_barcodes.txt |
#Remove unlikely barcodes probably from sequencing errors
awk -F "\t" 'BEGIN {OFS="\t"} $3!~/^[GN]+$/ && $4!~/^[GN]+$/'

#
##to debug unexplained barcodes
#zcat Undetermined_S0_R1_001.fastq.gz | awk -F ":" 'BEGIN {OFS="\t"} NR % 4 ==1 {split($10, bc, "+");  bc1=bc[1]; bc2=bc[2]} NR % 4 ==2 && bc2=="TAAGGC"' |m


echo
echo "Done!!!"
date
