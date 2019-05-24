#!/bin/bash
set -e -o pipefail

src=/vol/mauranolab/mapped/src/transposon/

chromsizes=/vol/isg/annotation/fasta/hg38/hg38.chrom.sizes


getcolor () {
    local trackcolor
    
    trackcolor="51,128,195" #blue
    if [[ "$1" =~ A1fw ]] || [[ "$1" =~ C1fw ]] || [[ "$1" =~ A1rv ]] || [[ "$1" =~ C1rv ]]; then
        trackcolor="120,88,165" #purple
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
echo "Generating UCSC track"
trackcolor=$(getcolor $PREFIX)


numUCSCsites=`cat ${OUTBASE}.zscore.bed | wc -l`
echo "Num of sites in browser track: ${numUCSCsites}"

cut -f1-3,5 ${OUTBASE}.zscore.bed > $TMPDIR/${PREFIX}.zscore.bedGraph

bedGraphToBigWig $TMPDIR/${PREFIX}.zscore.bedGraph ${chromsizes} ${OUTBASE}.zscore.bw
UCSCbaseURL="https://mauranolab@cascade.isg.med.nyu.edu/~mauram01/transposon/${OUTBASE}"
echo "track name=${PREFIX}-activity description=\"${PREFIX} activity (zscore of log(RNA/DNA) scaled onto [0,1]), ${numUCSCsites} sites\" maxHeightPixels=30 color=$trackcolor viewLimits=0:1 autoScale=off visibility=full db=hg38 type=bigWig bigDataUrl=${UCSCbaseURL}.zscore.bw"


echo
echo "Generating smoothed activity UCSC track"


cat ${chromsizes} | 
awk -F "\t" 'BEGIN {OFS="\t"} $1!="chrM" && $1!="chrEBV"' |
awk -F "\t" '{OFS="\t"; print $1, 0, $2}' | sort-bed - | cut -f1,3 | awk -v step=10000 -v binwidth=10000 'BEGIN {OFS="\t"} {for(i=0; i<=$2-binwidth; i+=step) {print $1, i, i+binwidth, "."} }' | 
#--faster is ok since we are dealing with bins and read starts
bedmap --faster --delim "\t" --bp-ovr 1 --echo --mean - ${OUTBASE}.zscore.bed | 
#Leave NAN for wigToBigWig for now
#perl -pe 's/NAN$/NA/g;' |
#resize intervals down from full bin width to step size
#Intervals then conform to Richard's convention that the counts are reported in 20bp windows including reads +/-75 from the center of that window
awk -v step=10000 -v binwidth=10000 -F "\t" 'BEGIN {OFS="\t"} {offset=(binwidth-step)/2 ; $2+=offset; $3-=offset; print}' |

tee $TMPDIR/${PREFIX}.zscore.density.bed |
#Remember wig is 1-indexed
#NB assumes span == step
#bedgraph would be simpler but produces 5x larger file. Would variablestep result in smaller .bw file?
awk -F "\t" 'lastChrom!=$1 || $2 != lastChromEnd || $3-$2 != curspan {curstep=$3-$2; curspan=curstep; print "fixedStep chrom=" $1 " start=" $2+1 " step=" curstep " span=" curspan} {lastChrom=$1; lastChromStart=$2; lastChromEnd=$3; print $5}' > $TMPDIR/${PREFIX}.zscore.density.wig

awk -F "\t" 'BEGIN {OFS="\t"} $5!=0' $TMPDIR/${PREFIX}.zscore.density.bed | starch - > ${OUTBASE}.zscore.density.starch


#Kent tools can't use STDIN
wigToBigWig $TMPDIR/${PREFIX}.zscore.density.wig ${chromsizes} ${OUTBASE}.zscore.density.bw

echo "track name=${PREFIX}-activitydens description=\"${PREFIX} activity density (smoothed zscore of log(RNA/DNA) scaled onto [0,1]), ${numUCSCsites} sites\" maxHeightPixels=30 color=$trackcolor viewLimits=0:1 autoScale=off visibility=full db=hg38 type=bigWig bigDataUrl=${UCSCbaseURL}.zscore.density.bw"


echo
echo "Done!!!"
date

