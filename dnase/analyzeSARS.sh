#!/bin/bash

if [ "$#" -ge 1 ]; then
    dir="$1"
else
    dir="."
fi

echo
echo "Processing ${dir}"

sampletype=`basename ${dir}`

#Not perfect
fc=`pwd | xargs basename`

echo -e "FC\tSample_Type\tSample\tBS\tnumAnalyzedViralReads\tnumAnalyzedReads\tMean_Viral_Coverage\tNum_bp_20x" > ${dir}/viralcoverage.txt
for bamfile in `find ${dir} -mindepth 2 -name "*.bam" -not -name "*R1R2*" -not -path "*trash*"`; do
    sample=`basename ${bamfile} .hg38_full_wuhCor1.bam`
    sampleid=`echo ${sample} | cut -d "-" -f1`
    bs=`echo ${sample} | cut -d "-" -f2`
    echo -n -e "${fc}\t${sampletype}\t${sampleid}\t${bs}\t"
    
    samtools view -F 1024 -c ${bamfile} NC_045512v2 | perl -pe 's/\n/\t/g;'
    
    samtools view -F 1024 -c ${bamfile} | perl -pe 's/\n/\t/g;'
    
    coveragefile=`echo ${bamfile} | perl -pe 's/.bam/.coverage.starch/g;'`
    awk -F "\t" '{OFS="\t"; print $1, 0, $2}' /vol/isg/annotation/fasta/wuhCor1/wuhCor1.chrom.sizes | bedmap --prec 0 --wmean - ${coveragefile} | perl -pe 's/NAN$/0/g;' | perl -pe 's/\n/\t/g;'
    
    unstarch NC_045512v2 ${coveragefile} | awk '$5>=20' | awk -F "\t" 'BEGIN {sum=0} {sum+=$3-$2} END {print sum}'
done >> ${dir}/viralcoverage.txt


echo
echo -e "\nDone!"
