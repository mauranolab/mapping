#!/bin/bash
set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'


name=$1


sampleOutdir=$name
mappedgenome=hg38_noalt
chromsizes="/vol/isg/annotation/fasta/${mappedgenome}/${mappedgenome}.chrom.sizes"



analyzedReads=`unstarch --list ${sampleOutdir}/${name}.${mappedgenome}.reads.starch | perl -pe 's/ +\|/\t/g;' | awk -F "\t" 'BEGIN {OFS="\t"; sum=0} NR==1 {for(i=1; i<=NF; i++) {if($i=="chr") {chromCol=i} else if($i=="uncompressedLineCount") {uncompressedLineCountCol=i} } } NR>1 && $chromCol!="chrM" {sum+=$uncompressedLineCountCol} END {print sum}'`


echo "Doing sample ${name}, ${analyzedReads} reads"
date

echo
echo "Making density track"

#Make density track of number of overlapping reads per 150-bp window
#BUGBUG double counts fragments where both reads are in window
cat ${chromsizes} | 
awk -F "\t" 'BEGIN {OFS="\t"} $1!="chrM" && $1!="chrEBV"' |
awk '{OFS="\t"; $3=$2; $2=0; print}' | sort-bed - | cut -f1,3 | awk -v step=20 -v binwidth=150 'BEGIN {OFS="\t"} {for(i=0; i<=$2-binwidth; i+=step) {print $1, i, i+binwidth, "."} }' | 
#--faster is ok since we are dealing with bins and read starts
bedmap --faster --delim "\t" --bp-ovr 1 --echo --count - ${sampleOutdir}/${name}.${mappedgenome}.reads.starch |
#resize intervals down from full bin width to step size
#Intervals then conform to Richard's convention that the counts are reported in 20bp windows including reads +/-75 from the center of that window
awk -v step=20 -v binwidth=150 -F "\t" 'BEGIN {OFS="\t"} {offset=(binwidth-step)/2 ; $2+=offset; $3-=offset; print}' |
#Normalizes the density to 1M reads
awk -v analyzedReads=${analyzedReads} -F "\t" 'BEGIN {OFS="\t"} {$5=$5/analyzedReads*1000000; print}' |
#tee $TMPDIR/${name}.${mappedgenome}.density.bed |
#Remember wig is 1-indexed
#NB assumes span == step
#bedgraph would be simpler but produces 5x larger file. Would variablestep result in smaller .bw file?
awk -F "\t" 'lastChrom!=$1 || $2 != lastChromEnd || $3-$2 != curspan {curstep=$3-$2; curspan=curstep; print "fixedStep chrom=" $1 " start=" $2+1 " step=" curstep " span=" curspan} {lastChrom=$1; lastChromStart=$2; lastChromEnd=$3; print $5}' > $TMPDIR/${name}.${mappedgenome}.density.wig

#awk -F "\t" 'BEGIN {OFS="\t"} $5!=0' $TMPDIR/${name}.${mappedgenome}.density.bed | starch - > ${sampleOutdir}/${name}.${mappedgenome}.density.starch

#Kent tools can't use STDIN
wigToBigWig $TMPDIR/${name}.${mappedgenome}.density.wig ${chromsizes} ${sampleOutdir}/${name}.${mappedgenome}.density.bw

echo
echo "Done!"
date
