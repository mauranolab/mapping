#!/bin/bash
set -e -o pipefail

src=/vol/mauranolab/mapped/src/transposon

module load bedops

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'


echo "Combining signal with insertions" 
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

if [[ "${DNA}" == "none" ]]; then
    dnafile=/dev/null
    #All entries will be NA
    rnafile=`find ${RNA} -name "*.clones.counts.filtered.txt"`
    echo "Preprocessing 10xRNA files $rnafile"
    cat ${rnafile} | mlr --tsv --headerless-csv-output put '$count=$count/$nCells' then cut -f BC,count then sort -f BC > $TMPDIR/$PREFIX.counts.txt
    rnafile="$TMPDIR/$PREFIX.counts.txt"
else
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
fi
ipcrfile=`find ${iPCR} -name "*.barcodes.coords.bed"`

echo
echo "Analyzing ${PREFIX}"
echo "Input DNA ${dnafile}"
echo "Input RNA ${rnafile}"
echo "Input iPCR ${ipcrfile}"
echo

#Join all barcodes together and include NAs for each sample 
join -e NA -a1 -a2 ${dnafile} -o 0 1.2 2.2 ${rnafile} |
join -e NA -a1 -a2 -1 4 -2 1 -o 1.1 1.2 1.3 0 1.6 2.2 2.3 1.5 <(sort -k4,4 ${ipcrfile}) - | awk -F " " 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9}' > $TMPDIR/AllBCs.txt

header="chrom\tchromStart\tchromEnd\tBC\tstrand\tDNA\tRNA\tiPCR"
echo -e $header | cat - $TMPDIR/AllBCs.txt | awk -F "\t" 'BEGIN {OFS="\t"} {print $4, $6, $7, $8}' > ${OUTBASE}.AllBCs.txt


echo
echo -n "Number of starting sites: "
cat ${OUTBASE}.AllBCs.txt | wc -l


OUTPUT=${OUTBASE}.mapped.txt

echo
echo "Creating bed file for all BCs with a single integration site"
header="#chrom\tchromStart\tchromEnd\tBC\tvalue\tstrand\tDNA\tRNA\tiPCR"
cat $TMPDIR/AllBCs.txt | awk  -F "\t" 'BEGIN {OFS="\t"} $1!="NA" && $1!="chrY" && $1!="chrM" {print $1, $2, $3, $4, 0, $5, $6, $7, $8}' |
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


echo 'doing distance to nearest TSS'
header="$header\tDistToTSS\tNearestRefseqGeneTSS\tNearestRefseqGeneStrand"
closest-features --delim '\t' --closest --dist --no-ref $OUTPUT /vol/isg/annotation/bed/hg38/refseq_gene/refGene.CombinedTxStarts.bed | 
awk -F "\t" 'BEGIN {OFS="\t"} {print $NF, $(NF-3), $(NF-1)}' | paste $OUTPUT - > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo 'doing distance to nearest DHS '
header="$header\tDistToNearestDHS"
sort-bed $OUTPUT | closest-features --closest --dist --no-ref - /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch |
awk -F'|' 'BEGIN {OFS="\t"} {print $2}' | paste $OUTPUT - > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo 'doing MCV of nearest DHS'
header="$header\tMCVNearestDHS"
unstarch /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "K562-DS9764"}' |
bedops -u /home/mauram01/scratch/hybridmice/dnase/lineagemcv/all.lineage.dhs.unmerged.starch - |
bedmap  --delim '\t' --skip-unmapped --echo --fraction-ref 1.0 --count  /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch - | 
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, NR, $4}' > $TMPDIR/DHS_MCV.bed

closest-features --delim '\t' --dist --closest $OUTPUT $TMPDIR/DHS_MCV.bed | awk '{print $(NF-1)}' | paste $OUTPUT - > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo 'doing distance to nearest CTCF site'
header="$header\tDistToNearestCTCF"
sort-bed $OUTPUT | closest-features --closest --dist --no-ref - /vol/isg/encode/CTCF_maurano_cell_reports_2015/ctcf.bycelltype.K562.hg38.starch |
awk -F'|' 'BEGIN {OFS="\t"} {print $2}' | paste $OUTPUT - > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo 'doing distance to nearest DHS (excluding CTCF sites)'
bedops -n -1 /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch /vol/isg/encode/CTCF_maurano_cell_reports_2015/ctcf.bycelltype.K562.hg38.starch > $TMPDIR/K562-DS9764.noctcf.hg38_noalt.fdr0.01.pks.bed


header="$header\tDistToNearestDHSnoCTCF"
sort-bed $OUTPUT | closest-features --closest --dist --no-ref - $TMPDIR/K562-DS9764.noctcf.hg38_noalt.fdr0.01.pks.bed |
awk -F'|' 'BEGIN {OFS="\t"} {print $2}' | paste $OUTPUT - > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo 'doing number of DHSs within +/-100 kbp'
header="$header\tNumDHS100kb"
bedmap --delim '\t' --range 100000 --echo --count $OUTPUT /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo 'doing number of CTCF sites within +/-10 kbp'
header="$header\tNumCTCF10kb"
bedmap --delim '\t' --range 10000 --echo --count $OUTPUT /vol/isg/encode/CTCF_maurano_cell_reports_2015/ctcf.bycelltype.K562.hg38.starch > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo 'doing number of CTCF sites within +/-50 kbp'
header="$header\tNumCTCF50kb"
bedmap --delim '\t' --range 50000 --echo --count $OUTPUT /vol/isg/encode/CTCF_maurano_cell_reports_2015/ctcf.bycelltype.K562.hg38.starch > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo 'doing number of CTCF sites within +/-100 kbp'
header="$header\tNumCTCF100kb"
bedmap --delim '\t' --range 100000 --echo --count $OUTPUT /vol/isg/encode/CTCF_maurano_cell_reports_2015/ctcf.bycelltype.K562.hg38.starch > $OUTPUT.new
mv $OUTPUT.new $OUTPUT



#Print to file
echo -e $header | cat - $OUTPUT > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


${src}/analyzeCombinedSignalWithInsertions.sh ${OUTBASE}


echo
echo "Done!!!"
date

