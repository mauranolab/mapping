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

echo -e "Sample\tNum_sequenced_reads\tNonredundant_reads_analyzed\tNonredundant_reads_on_target\tcov.wholegenome\tcov.mappablegenome\tcov.mean\tcov.sd\tcov.5th.pctile\tcov.median"
for covfile in `find ${dirs} -name "*.${mappedgenome}.coverage.binned.starch"`; do
    base=`basename ${covfile} .${mappedgenome}.coverage.binned.starch`
    echo -n -e "${base}\t"
    
    readsfile=`echo ${covfile} | perl -pe 's/\.coverage.binned.starch$/.reads.starch/g;'`
    bamfile=`echo ${covfile} | perl -pe 's/\.coverage.binned.starch$/.bam/g;'`
    
    analysisFile=`dirname ${covfile} | xargs -I {} find {} -name "analysis.*.${mappedgenome}.o*"`
    sequencedReads=`awk -F "\t" 'BEGIN {OFS="\t"} $1=="Num_sequenced_reads" {print $2; exit}' ${analysisFile}`
    echo -n -e "${sequencedReads}\t"
    
    #Nonredundant_reads_analyzed
    samtools view -F 1536 -c ${bamfile} | perl -pe 's/\n/\t/g;'
    
    #Nonredundant_reads_on_target
    bedops --chrom ${chrom} -e -1 ${readsfile} $TMPDIR/target.bed | awk -F "\t" 'BEGIN {OFS="\t"} !and($7, 1024)' | wc -l | perl -pe 's/\n/\t/g;'
    
    #cov.wholegenome
    unstarch ${covfile} | awk -F "\t" 'BEGIN {OFS="\t"; sum=0; count=0} {sum+=$5; count+=1} END {if(count==0) {print "NA"} else {print sum/count}}' | perl -pe 's/\n/\t/g;'
    
    #cov.mappablegenome
    unstarch ${covfile} | bedops -e -1 - ${mappableFile} | awk -F "\t" 'BEGIN {OFS="\t"; sum=0; count=0} {sum+=$5; count+=1} END {if(count==0) {print "NA"} else {print sum/count}}' | perl -pe 's/\n/\t/g;'
    
    #cov.mean, cov.sd, cov.5th.pctile, cov.median
    #bedmap --chrom ${chrom} --faster --delim "\t" --bp-ovr 1 --prec 2 --mean --stdev --kth 0.05 --median $TMPDIR/target.bed ${covfile}
    #Use mlr to handle discontinuous target regions
    bedmap --chrom ${chrom} --faster --delim "\n" --multidelim "\n" --bp-ovr 1 --prec 3 --echo-map-score $TMPDIR/target.bed ${covfile} | mlr  --tsv --implicit-csv-header --headerless-csv-output stats1 -a mean,stddev,p5,p50 -f 1
done



echo
echo "Done!"
date