#!/bin/bash

set -e -o pipefail

module load bedops

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'


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

chromsizes=/vol/isg/annotation/fasta/hg38/hg38.chrom.sizes


echo "Starting analysis" 
date 
#first input should be DNA (sample.barcode.counts.txt file)
DNA=$1
#Second input should be RNA (sample.barcode.counts.txt file)
RNA=$2
#Third input should be iPCR (sample.barcodes.coords.bed file)
iPCR=$3
#Fourth input should be output location + prefix
OUTBASE=$4
PREFIX=`basename $OUTBASE`
#TMPDIR=/tmp

dnafile=`find ${DNA} -name "*.barcode.counts.UMI_corrected.txt"`
if [ -z "${dnafile}" ]; then
    #Fall back on non-UMI counts
    dnafile=`find ${DNA} -name "*.barcode.counts.txt"`
fi
rnafile=`find ${RNA} -name "*.barcode.counts.UMI_corrected.txt"`
if [ -z "${rnafile}" ]; then
    #Fall back on non-UMI counts
    rnafile=`find ${RNA} -name "*.barcode.counts.txt"`
fi
ipcrfile=`find ${iPCR} -name "*.barcodes.coords.bed"`

echo
echo "Analyzing ${PREFIX}"
echo "Input DNA ${dnafile}"
echo "Input RNA ${rnafile}"
echo "Input iPCR ${ipcrfile}"
echo

trackcolor=$(getcolor $PREFIX)

#Join all barcodes together and include NAs for each sample 
join -e NA -a1 -a2 ${dnafile} -o 0 1.2 2.2 ${rnafile} |
join -e NA -a1 -a2 -1 4 <(sort -k4,4 ${ipcrfile}) -o 1.1 1.2 1.3 0 1.6 2.2 2.3 1.5 -2 1 - | awk -F " " 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' > $TMPDIR/AllBCs.txt

header="chrom\tchromStart\tchromEnd\tBC\tstrand\tDNA\tRNA\tiPCR"
echo -e $header | cat - $TMPDIR/AllBCs.txt | awk -F "\t" 'BEGIN {OFS="\t"} {print $4, $6, $7, $8}' > ${OUTBASE}.AllBCs.txt


echo
echo -n "Number of starting sites: "
cat ${OUTBASE}.AllBCs.txt | wc -l


OUTPUT=${OUTBASE}.txt

#Remove iPCR insertions with more than one location
cat ${ipcrfile} | awk -F "\t" 'BEGIN {OFS="\t"} {print $4}' | sort  | uniq -c | sort -nk1 | awk '$1==1 {print $2}' > $TMPDIR/${PREFIX}.singleIns.txt

echo
echo "Creating bed file for all BCs with a single integration site"
header="chrom\tchromStart\tchromEnd\tBC\tvalue\tstrand\tDNA\tRNA\tiPCR"
cat $TMPDIR/AllBCs.txt | awk  -F "\t" 'BEGIN {OFS="\t"} $1!="NA" && $1!="chrY" && $1!="chrM" {print $1, $2, $3, $4, 0, $5, $6, $7, $8}' |
#Retains barcodes from $OUTPUT that are present in $TMPDIR/${PREFIX}.singleIns.txt
awk -F "\t" 'BEGIN {OFS="\t"} NR==FNR{a[$1];next} ($4) in a' $TMPDIR/${PREFIX}.singleIns.txt - |
sort-bed - > $OUTPUT.new 
mv $OUTPUT.new $OUTPUT


echo
echo -n "Number of remaining sites: "
cat $OUTPUT | wc -l



echo
echo "Adding annotation"
#Sample name
echo "doing sample name"
header="$header\tName"
awk -F "\t" -v sampleName="${PREFIX}" 'BEGIN {OFS="\t"} {print $0, sampleName}' $OUTPUT > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


#Distance to TSS
echo 'doing distance to nearest TSS'
header="$header\tDistToTSS\tNearestRefseqGeneTSS\tNearestRefseqGeneStrand"
closest-features --delim '\t' --closest --dist --no-ref $OUTPUT /vol/isg/annotation/bed/hg38/refseq_gene/refGene.CombinedTxStarts.bed | 
awk -F "\t" 'BEGIN {OFS="\t"} {print $NF, $(NF-3), $(NF-1)}' | paste $OUTPUT - > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


#Distance to DNase
echo 'doing distance to nearest DHS '
header="$header\tDistToNearestDHS"
sort-bed $OUTPUT | closest-features --closest --dist --no-ref - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.pks.bed |
awk -F'|' 'BEGIN {OFS="\t"} {print $2}' | paste $OUTPUT - > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo 'doing MCV of nearest DHS'
header="$header\tMCVNearestDHS"
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "K562-DS9764"}' /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.pks.bed |
bedops -u /home/mauram01/scratch/hybridmice/dnase/lineagemcv/all.lineage.dhs.unmerged.starch - |
bedmap  --delim '\t' --skip-unmapped --echo --fraction-ref 1.0 --count  /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.pks.bed - | 
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, NR, $4}' > $TMPDIR/DHS_MCV.bed
closest-features --delim '\t' --dist --closest $OUTPUT $TMPDIR/DHS_MCV.bed | awk '{print $(NF-1)}' | paste $OUTPUT - > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


#Print to file
echo -e $header | cat - $OUTPUT > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo
echo "Processing data in R"
/vol/mauranolab/mapped/src/transposon/combineSignalWithInsertions.R ${PREFIX} ${OUTBASE}


echo
echo "Generating UCSC track"
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

