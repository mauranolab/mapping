#!/bin/bash
set -eu -o pipefail


mappedgenome=$1
bait=$2

shift
shift
#A list of sample directories output by mapping pipeline
dirs=$@


case "${mappedgenome}" in
hg38_full)
    annotationgenome="hg38"
    mappableFile="/vol/isg/annotation/bed/${annotationgenome}/mappability/${annotationgenome}.K36.mappable_only.starch"
    ;;
mm10)
    annotationgenome="mm10"
    mappableFile="/vol/isg/annotation/bed/${annotationgenome}/mappability/${annotationgenome}.K36.mappable_only.starch"
    ;;
cegsvectors*)
    annotationgenome="cegsvectors"
    mappableFile="/dev/null"
    ;;
rn6)
    annotationgenome="rn6"
    mappableFile="/dev/null"
    ;;
*)
    echo "ERROR: Don't recognize genome ${mappedgenome}";
    exit 1;;
esac


case "${bait}" in
HPRT1_bait)
    grep -w "HPRT1_assembly" /vol/cegs/sequences/hg38/HPRT1/HPRT1_assembly.bed | bedops -m - > $TMPDIR/target.bed
    ;;
RnHoxa_bait)
    grep -w "RnHoxa_assembly" /vol/cegs/sequences/rn6/RnHoxa/RnHoxa_assembly.bed | bedops -m - > $TMPDIR/target.bed
    ;;
LP087_bait)
    #TODO not counting backbone since we don't support multiple chromosomes
    cat /vol/cegs/sequences/cegsvectors/vectors.incells.chrom.sizes | awk -F "\t" 'BEGIN {OFS="\t"} $1=="LP087" {print $1, 0, $2}' > $TMPDIR/target.bed
    ;;
*)
    echo "ERROR: Don't recognize bait ${bait}";
    exit 2;;
esac


#speed things up for local analyses
if [ `cut -f1 $TMPDIR/target.bed | uniq | wc -l` -ne 1 ]; then
    echo "ERROR: can't handle bait spanning >1 chrom"
    exit 3
fi
chrom=`cut -f1 $TMPDIR/target.bed | uniq`

echo "Analyzing ${bait} bait mapped to ${mappedgenome}"


#Summary of on-target coverage
#NB genomecov is over all chrom.sizes -- right way?

echo -e "Sample\tBait\tNum_sequenced_reads\tNonredundant_reads_analyzed\tNonredundant_reads_on_target\tDuplicate_reads_on_target\tcov.wholegenome\tcov.mean\tcov.median\tcov.5th.pctile\tcov.sd"
for covfile in `find ${dirs} -name "*.${mappedgenome}.coverage.binned.starch"`; do
    base=`basename ${covfile} .${mappedgenome}.coverage.binned.starch`
    echo -n -e "${base}\t${bait}\t"
    
    readsfile=`echo ${covfile} | perl -pe 's/\.coverage.binned.starch$/.reads.starch/g;'`
    bamfile=`echo ${covfile} | perl -pe 's/\.coverage.binned.starch$/.bam/g;'`
    
    analysisFile=`dirname ${covfile} | xargs -I {} find {} -name "analysis.*.${mappedgenome}.o*"`
    sequencedReads=`awk -F "\t" 'BEGIN {OFS="\t"} $1=="Num_sequenced_reads" {print $2; exit}' ${analysisFile}`
    echo -n -e "${sequencedReads}\t"
    
    #Nonredundant_reads_analyzed
    samtools view -F 1536 -c ${bamfile} | perl -pe 's/\n/\t/g;'
    
    #Nonredundant_reads_on_target, Duplicate_reads_on_target
    bedops --chrom ${chrom} -e -1 ${readsfile} $TMPDIR/target.bed | awk -F "\t" 'BEGIN {OFS="\t"; nonredundantReads=0; dupReads=0} {if(and($7, 1024)) {dupReads+=1} else {nonredundantReads+=1}} END {print nonredundantReads,dupReads}' | perl -pe 's/\n/\t/g;'
    
    #cov.wholegenome
    unstarch ${covfile} | awk -F "\t" 'BEGIN {OFS="\t"; sum=0; count=0} {sum+=$5; count+=1} END {if(count==0) {print "NA"} else {print sum/count}}' | perl -pe 's/\n/\t/g;'
    
    #cov.mappablegenome
    #unstarch ${covfile} | bedops -e -1 - ${mappableFile} | awk -F "\t" 'BEGIN {OFS="\t"; sum=0; count=0} {sum+=$5; count+=1} END {if(count==0) {print "NA"} else {print sum/count}}' | perl -pe 's/\n/\t/g;'
    
    #cov.mean, cov.sd, cov.5th.pctile, cov.median
    #bedmap --chrom ${chrom} --faster --delim "\t" --bp-ovr 1 --prec 2 --mean --stdev --kth 0.05 --median $TMPDIR/target.bed ${covfile}
    #Use mlr to handle discontinuous target regions
    bedmap --chrom ${chrom} --faster --delim "\n" --multidelim "\n" --bp-ovr 1 --prec 3 --echo-map-score $TMPDIR/target.bed ${covfile} | mlr --ofmt %.1f --tsv --implicit-csv-header --headerless-csv-output stats1 -a mean,p50,p5,stddev -f 1
done



echo
echo "Done!"
date
