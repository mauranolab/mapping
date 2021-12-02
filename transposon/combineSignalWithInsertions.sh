#!/bin/bash
set -e -o pipefail

src=$( dirname "${BASH_SOURCE[0]}" )

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


echo 'doing number of DHSs within +/-5 kbp'
header="$header\tNumDHS5kb"
bedmap --delim '\t' --range 5000 --echo --count $OUTPUT /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo 'doing number of DHSs within +/-10 kbp'
header="$header\tNumDHS10kb"
bedmap --delim '\t' --range 10000 --echo --count $OUTPUT /vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch > $OUTPUT.new
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


echo 'doing local insertion density'
header="$header\tInsDens1Mb"
bedmap --delim '\t' --range 1000000 --echo --count $OUTPUT ${ipcrfile} > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo 'doing local insertion density'
header="$header\tInsDens"
bedmap --delim '\t' --range 100000 --echo --count $OUTPUT ${ipcrfile} > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo 'doing TADid'
header="$header\tTADid"
bedmap --delim '\t' --bp-ovr 1 --echo --echo-map-id $OUTPUT /home/maagj01/scratch/transposon/Analysis/K562_reference_epigenome/Hi-C/K562/raodomains_5kb_KR/hg38/TADs/Id/Armatus_5kb_gamma.0.5.0.bed > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo 'doing TSS density'
header="$header\tNumTSS100kb"
bedmap --delim '\t' --range 100000 --echo --count $OUTPUT /vol/isg/annotation/bed/hg38/refseq_gene/refGene.CombinedTxStarts.bed > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo "doing Overlap with CGI"
header="$header\tCGIovr"
awk -F "\t" 'BEGIN {OFS="\t"} {$4="TRUE"; print}' /vol/isg/annotation/bed/hg38/cpg_islands/cpgIslands.bed |
bedmap --faster --delim '\t' --echo --echo-map-id $OUTPUT - | awk -F "\t" 'BEGIN {OFS="\t"} $NF~/TRUE/ {$NF="TRUE"} $NF=="" {$NF="FALSE"} {print}' > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo "doing G+C content +/-75bp"
header="$header\tPercentGC"
awk -F "\t" -v widen=75 'BEGIN {OFS="\t"} {$2=$2 > widen ? $2-widen : 0; $3+=widen; print;}' $OUTPUT | ~/bin/bed2fasta.pl - /vol/isg/annotation/fasta/hg38 2>/dev/null |
grep -v -e "^>" | tr '[a-z]' '[A-Z]' | perl -ne 'chomp; print length($_) != 0 ? ($_ =~ tr/[gcGC]//) / length($_) : "NA"; print "\n"' | paste $OUTPUT - > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo "doing GenicLocation.simple"
header="$header\tGenicLocation.simple"
gcat /vol/isg/annotation/bed/hg38/refseq_paint/refGene_hg38.bed6.starch | grep -v TxS-1000 | grep -v TxStart-10kb > $TMPDIR/refGene_hg38.simple.bed6
bedmap --delim '\t' --echo-map-id $OUTPUT $TMPDIR/refGene_hg38.simple.bed6 | awk -F "\t" 'BEGIN {OFS="\t"} {col=NF} $col~/coding/ {$col="coding"} $col~/promoter/ {$col="promoter"} $col~/3.UTR/ {$col="3UTR"} $col~/5.UTR/ {$col="5UTR"} $col~/intron/ {$col="intron"} $col=="" || $col~/3.proximal/ {$col="intergenic"} {print;}' | paste $OUTPUT - > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


echo "doing repeatmasker density"
#/vol/isg/annotation/bed/hg38/repeat_masker/DNA?.bed /vol/isg/annotation/bed/hg38/repeat_masker/SINE?.bed /vol/isg/annotation/bed/hg38/repeat_masker/LTR?.bed /vol/isg/annotation/bed/hg38/repeat_masker/RC?.bed
rmsk=" /vol/isg/annotation/bed/hg38/repeat_masker/DNA.bed /vol/isg/annotation/bed/hg38/repeat_masker/LINE.bed /vol/isg/annotation/bed/hg38/repeat_masker/Low_complexity.bed  /vol/isg/annotation/bed/hg38/repeat_masker/LTR.bed  /vol/isg/annotation/bed/hg38/repeat_masker/RC.bed /vol/isg/annotation/bed/hg38/repeat_masker/Retroposon.bed /vol/isg/annotation/bed/hg38/repeat_masker/RNA.bed /vol/isg/annotation/bed/hg38/repeat_masker/rRNA.bed /vol/isg/annotation/bed/hg38/repeat_masker/Satellite.bed /vol/isg/annotation/bed/hg38/repeat_masker/scRNA.bed /vol/isg/annotation/bed/hg38/repeat_masker/Simple_repeat.bed /vol/isg/annotation/bed/hg38/repeat_masker/SINE.bed /vol/isg/annotation/bed/hg38/repeat_masker/snRNA.bed /vol/isg/annotation/bed/hg38/repeat_masker/srpRNA.bed /vol/isg/annotation/bed/hg38/repeat_masker/tRNA.bed /vol/isg/annotation/bed/hg38/repeat_masker/Unknown.bed"
for rmskf in ${rmsk}; do
    base=`basename ${rmskf} .bed`
    echo "${base}"
    header="$header\trmsk.${base}"
    bedmap --delim '\t' --range 100000 --echo --count $OUTPUT ${rmskf} > $OUTPUT.new
    mv $OUTPUT.new $OUTPUT
done


echo "doing Overlap with LAD"
header="$header\tLADovr"
awk -F "\t" 'BEGIN {OFS="\t"} {$4="TRUE"; print}' /vol/mauranolab/publicdata/4DN/4DNFIT7Q8TTV.bed |
bedmap --faster --delim '\t' --echo --echo-map-id $OUTPUT - | awk -F "\t" 'BEGIN {OFS="\t"} $NF~/TRUE/ {$NF="TRUE"} $NF=="" {$NF="FALSE"} {print}' > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


#Print to file
echo -e $header | cat - $OUTPUT > $OUTPUT.new
mv $OUTPUT.new $OUTPUT


${src}/analyzeCombinedSignalWithInsertions.sh ${OUTBASE}


echo
echo "Done!!!"
date

