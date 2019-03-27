#!/bin/bash
set -eu -o pipefail

src=/vol/mauranolab/transposon/src

#TODO inconsistently applied
minReadCutoff=$1
sample=$2


OUTDIR=${sample}

if [ ! -s "${OUTDIR}/${sample}.barcodes.preFilter.txt.gz" ]; then
    echo "analyzeBCcounts.sh ERROR: barcode input file ${OUTDIR}/${sample}.barcodes.preFilter.txt.gz does not exist!"
    exit 1
fi


echo "Analyzing data for ${sample} (minReadCutoff=${minReadCutoff})"
zcat ${OUTDIR}/${sample}.barcodes.preFilter.txt.gz | ${src}/removeOverrepresentedBCs.py --col 1 --umicol 3 --freq 0.01 -o ${TMPDIR}/${sample}.barcodes.txt -


if [ `cat ${TMPDIR}/${sample}.barcodes.txt | cut -f1 | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | wc -l` -eq 0 ]; then
    echo "analyzeBCcounts.sh WARNING: no barcodes left, so exiting!"
    echo "Done!!!"
    exit 0
fi


pigz -p ${NSLOTS} -c -9 ${TMPDIR}/${sample}.barcodes.txt > ${OUTDIR}/${sample}.barcodes.txt.gz

###Overall summary statistics
echo
echo -n -e "${sample}\tNumber of total reads\t"
zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | wc -l
echo -n -e "${sample}\tNumber of total read barcodes\t"
zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | awk -F "\t" '$1!="" {count+=1} END {print count}'


zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | cut -f1 | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | shuf -n 1000000 | awk -F "\t" 'BEGIN {OFS="\t"} {print ">id-" NR; print}' |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} barcodes" --stacks-per-line 100 > $TMPDIR/${sample}.barcodes.eps
convert $TMPDIR/${sample}.barcodes.eps ${OUTDIR}/${sample}.barcodes.png

##UMIs weblogo
minUMILength=`zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | awk -F "\t" 'BEGIN {OFS="\t"} $1!="" {print length($3)}' | uniq | sort -n | awk 'NR==1'`
echo
echo -e "${sample}\tMinimum UMI length\t${minUMILength}"
if [ "${minUMILength}" -gt 0 ]; then
    echo "Generating UMI weblogo"
    #Need to trim in case not all the same length
    zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | cut -f3 | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | shuf -n 1000000 | awk -v trim=${minUMILength} -F "\t" 'BEGIN {OFS="\t"} {print substr($0, 0, trim)}' | awk -F "\t" 'BEGIN {OFS="\t"} {print ">id-" NR; print}' |
    weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} UMI" --stacks-per-line 100 > $TMPDIR/${sample}.UMIs.eps
    convert $TMPDIR/${sample}.UMIs.eps ${OUTDIR}/${sample}.UMIs.png
fi


if [ "${minUMILength}" -ge 5 ]; then
    echo
    echo "Analyzing UMIs"
    
    echo
    echo "Barcode counts (before UMI)"
    zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz |
    awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | cut -f1,3 |
    sort -k1,1 -k2,2 | uniq -c | sort -k2,2 |
    #Format: BC, UMI, n
    awk 'BEGIN {OFS="\t"} {print $2, $3, $1}' |
    pigz -p ${NSLOTS} -c -9 > ${OUTDIR}/${sample}.barcode.counts.withUMI.txt.gz
    zcat -f ${OUTDIR}/${sample}.barcode.counts.withUMI.txt.gz | cut -f1 | sort | uniq -c | awk 'BEGIN {OFS="\t"} {print $2, $1}' > ${OUTDIR}/${sample}.barcode.counts.UMI_corrected.txt
    echo -n -e "${sample}\tNumber of unique barcodes+UMI\t"
    #TODO move to extract.py
    zcat -f ${OUTDIR}/${sample}.barcode.counts.withUMI.txt.gz | wc -l
    
    echo
    echo "Histogram of number of reads per barcode+UMI"
    zcat -f ${OUTDIR}/${sample}.barcode.counts.withUMI.txt.gz | cut -f3 | awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g
    
    echo
    echo "UMI lengths"
    #No more empty BCs are present so no need to check
    zcat -f ${OUTDIR}/${sample}.barcode.counts.withUMI.txt.gz | cut -f2 |
    awk '{print length($0)}' | sort -g | uniq -c | sort -k2,2 | awk '$2!=0'
    
    #TODO why bother with the bcfile variable since it's the same in both cases?
    bcfile="${OUTDIR}/${sample}.barcodes.txt.gz"
else
    echo
    echo "No UMIs found"
    
    bcfile="${OUTDIR}/${sample}.barcodes.txt.gz"
fi

echo
echo "Barcode counts"
#Empty BCs haven't been filtered yet for non-UMI data
zcat -f ${bcfile} | awk -F "\t" '$1!=""' | cut -f1 | sort | uniq -c | awk 'BEGIN {OFS="\t"} {print $2, $1}' > ${OUTDIR}/${sample}.barcode.counts.txt

echo -n -e "${sample}\tNumber of unique barcodes\t"
#TODO move to extract.py
cat ${OUTDIR}/${sample}.barcode.counts.txt | wc -l

echo -n -e "${sample}\tNumber of unique barcodes passing minimum read cutoff\t"
#TODO move to extract.py
cat ${OUTDIR}/${sample}.barcode.counts.txt | awk -v minReadCutoff=${minReadCutoff} -F "\t" 'BEGIN {OFS="\t"} $2>=minReadCutoff' | wc -l

echo -n -e "${sample}\tNumber of analyzed reads\t"
cat ${OUTDIR}/${sample}.barcode.counts.txt | awk -v minReadCutoff=${minReadCutoff} -F "\t" 'BEGIN {OFS="\t"; sum=0; total=0} $2>=minReadCutoff {sum+=$2} {total+=$2} END {print sum}'


echo
echo "Histogram of number of reads per barcode"
cat ${OUTDIR}/${sample}.barcode.counts.txt | cut -f2 | awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g


echo
echo "Barcode lengths"
cut -f1 ${OUTDIR}/${sample}.barcode.counts.txt |
awk '{print length($0)}' | sort -g | uniq -c | sort -k2,2 | awk '$2!=0'


echo
echo "Saturation curves"
numlines=`zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | wc -l`

for i in `echo {10,1000,10000,50000,100000,250000,500000,1000000,2000000,3000000,4000000,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000,50000000,55000000,60000000,70000000,80000000,90000000,100000000,110000000,120000000} | tr ' ' '\n' | gawk -v subset=$numlines '{if ($1<=subset) print}'`; do 
saturation=$(zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | shuf -n $i | awk -F "\t" '{print $1}' | sort | uniq | wc -l); echo $i  $saturation 'minreads1' ;done > ${OUTDIR}/${sample}.Saturation_minReads1.txt

if [[ `zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | awk -F "\t" 'BEGIN {OFS="\t"} length($3) > 4 {found=1} END {print found}'` == 1 ]]; then
    for i in `echo {10,1000,10000,50000,100000,250000,500000,1000000,2000000,3000000,4000000,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000,50000000,55000000,60000000,70000000,80000000,90000000,100000000,110000000,120000000} | tr ' ' '\n' | gawk -v subset=$numlines '{if ($1<=subset) print}'`; do 
    saturation=$(zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | shuf -n $i | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $3}' | sort | uniq | awk -F "\t" '{print $1}' | sort | uniq -c | awk '{if ($1>=10) print $2}' | wc -l); echo $i $saturation 'minreads10' ;done > ${OUTDIR}/${sample}.Saturation_minReads10.txt
else
    for i in `echo {10,1000,10000,50000,100000,250000,500000,1000000,2000000,3000000,4000000,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000,50000000,55000000,60000000,70000000,80000000,90000000,100000000,110000000,120000000} | tr ' ' '\n' | gawk -v subset=$numlines '{if ($1<=subset) print}'`; do 
    saturation=$(zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | shuf -n $i | awk -F "\t" '{print $1}' | sort | uniq -c | awk '{if ($1>=10) print $2}' | wc -l); echo $i $saturation 'minreads10' ;done > ${OUTDIR}/${sample}.Saturation_minReads10.txt
fi

echo
echo "Done!!!"
date
