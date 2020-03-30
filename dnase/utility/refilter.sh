#!/bin/bash
set -e -o pipefail

sample=$1
src=$2

permittedMismatches=3
minMAPQ=20

NSLOTS=1

echo
echo "Processing $sample bam"
date

samtools view -h bak.reMakeTracks/$sample/$sample.bam |
awk -F "\t" 'length($10) >= 27 || $1 ~ /^@/' |
samtools sort -@ $NSLOTS -m 2000M -O bam -T $TMPDIR/${sample}.sortbyname -l 1 -n - |
$src/filter_reads.py --max_mismatches $permittedMismatches --min_mapq $minMAPQ - - |
samtools sort -@ $NSLOTS -m 2000M -O bam -T $TMPDIR/${sample}.sort -l 9 - > $sample/$sample.bam


echo
echo "Indexing"
date
samtools index $sample/$sample.bam


echo
echo "samtools flagstat"
date
echo "Before"
samtools flagstat bak.reMakeTracks/$sample/$sample.bam

echo
echo "After"
samtools flagstat $sample/$sample.bam

echo
echo "Done!"
date
