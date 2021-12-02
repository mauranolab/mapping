#!/bin/bash
set -e -o pipefail

src=$( dirname "${BASH_SOURCE[0]}" )

chromsizes=/vol/isg/annotation/fasta/hg38/hg38.chrom.sizes


getcolor () {
    local trackcolor
    
    trackcolor="51,128,195" #blue
    if [[ "$1" =~ A1fw ]] || [[ "$1" =~ C1fw ]] || [[ "$1" =~ A1rv ]] || [[ "$1" =~ C1rv ]] || [[ "$1" =~ DoubleIns ]]; then
        trackcolor="255,153,0" #orange
    elif [[ "$1" =~ HS2 ]]; then
        trackcolor="238,54,36" #red
    fi
    
    echo ${trackcolor}
}


OUTBASE=$1

PREFIX=`basename $OUTBASE`
#TMPDIR=/tmp

echo
echo "Analyzing ${PREFIX}"


echo
echo "Processing data in R"
${src}/analyzeCombinedSignalWithInsertions.R ${PREFIX} ${OUTBASE}


echo
echo "Generating UCSC tracks"
trackcolor=$(getcolor $PREFIX)
UCSCbaseURL="https://cascade.isg.med.nyu.edu/~mauram01/transposon/${OUTBASE}"
numUCSCsites=`cat ${OUTBASE}.activity.bed | wc -l`
echo "Num of sites in browser track: ${numUCSCsites}"


echo
echo "Generating insertion density UCSC track"
#NB only considers insertions passing all filters, including minimum DNA count

cat ${chromsizes} | 
awk -F "\t" 'BEGIN {OFS="\t"} $1!="chrM" && $1!="chrEBV"' |
awk -F "\t" '{OFS="\t"; print $1, 0, $2}' | sort-bed - | cut -f1,3 | awk -v step=10000 -v binwidth=10000 'BEGIN {OFS="\t"} {for(i=0; i<=$2-binwidth; i+=step) {print $1, i, i+binwidth, "."} }' | 
#--faster is ok since we are dealing with bins and read starts
bedmap --faster --delim "\t" --bp-ovr 1 --echo --count - ${OUTBASE}.activity.bed | perl -pe 's/NAN$/0/g;' |
#normalize to insertions per 10^6
awk -v totalinsertions=${numUCSCsites} -F "\t" 'BEGIN {OFS="\t"} {$NF=$NF/totalinsertions*10^6; print}' |
#resize intervals down from full bin width to step size
#Intervals then conform to Richard's convention that the counts are reported in 20bp windows including reads +/-75 from the center of that window
awk -v step=10000 -v binwidth=10000 -F "\t" 'BEGIN {OFS="\t"} {offset=(binwidth-step)/2 ; $2+=offset; $3-=offset; print}' |
tee $TMPDIR/${PREFIX}.insertion.density.bed |
#Remember wig is 1-indexed
#NB assumes span == step
#bedgraph would be simpler but produces 5x larger file. Would variablestep result in smaller .bw file?
awk -F "\t" 'lastChrom!=$1 || $2 != lastChromEnd || $3-$2 != curspan {curstep=$3-$2; curspan=curstep; print "fixedStep chrom=" $1 " start=" $2+1 " step=" curstep " span=" curspan} {lastChrom=$1; lastChromStart=$2; lastChromEnd=$3; print $5}' > $TMPDIR/${PREFIX}.insertion.density.wig

cat $TMPDIR/${PREFIX}.insertion.density.bed | starch - > ${OUTBASE}.insertion.density.starch

#Kent tools can't use STDIN
wigToBigWig $TMPDIR/${PREFIX}.insertion.density.wig ${chromsizes} ${OUTBASE}.insertion.density.bw


echo "track name=${PREFIX}-insertiondens description=\"${PREFIX} normalized insertion density (20 kb windows), ${numUCSCsites} sites\" maxHeightPixels=15 color=$trackcolor viewLimits=0:20 autoScale=off visibility=full db=hg38 type=bigWig bigDataUrl=${UCSCbaseURL}.insertion.density.bw"


echo
echo "Generating activity UCSC track"
cut -f1-3,5 ${OUTBASE}.activity.bed > $TMPDIR/${PREFIX}.activity.bedGraph
bedGraphToBigWig $TMPDIR/${PREFIX}.activity.bedGraph ${chromsizes} ${OUTBASE}.activity.bw

echo "track name=${PREFIX}-activity description=\"${PREFIX} activity (log(RNA/DNA)), ${numUCSCsites} sites\" maxHeightPixels=30 color=$trackcolor viewLimits=0:4 autoScale=off visibility=full db=hg38 type=bigWig bigDataUrl=${UCSCbaseURL}.activity.bw"
sort-bed ${OUTBASE}.activity.bed | starch - > ${OUTBASE}.activity.starch
rm -f ${OUTBASE}.activity.bed


echo
echo "Done!!!"
date

