#!/bin/bash
set -eu -o pipefail

src=/vol/mauranolab/transposon/src

sample=$1


#TODO inconsistently applied
minReadCutoff=10

OUTDIR=${sample}

if [ ! -s "${OUTDIR}/${sample}.barcodes.preFilter.txt" ]; then
    echo "analyzeBCcounts.sh ERROR: barcode input file ${OUTDIR}/${sample}.barcodes.preFilter.txt does not exist!"
    exit 1
fi


echo "Analyzing data for ${sample} (minReadCutoff=${minReadCutoff})"
${src}/removeOverrepresentedBCs.py --col 1 --freq 0.01 -o ${OUTDIR}/${sample}.barcodes.txt ${OUTDIR}/${sample}.barcodes.preFilter.txt


if [ `cat ${OUTDIR}/${sample}.barcodes.txt | cut -f1 | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | wc -l` -eq 0 ]; then
    echo "analyzeBCcounts.sh WARNING: no barcodes left, so exiting!"
    echo "Done!!!"
    exit 0
fi


###Overall summary statistics
cat ${OUTDIR}/${sample}.barcodes.txt | cut -f1 | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | shuf -n 1000000 | awk -F "\t" 'BEGIN {OFS="\t"} {print ">id-" NR; print}' |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} barcodes" --stacks-per-line 100 > $TMPDIR/${sample}.barcodes.eps
convert $TMPDIR/${sample}.barcodes.eps ${OUTDIR}/${sample}.barcodes.png

##UMIs weblogo
if [[ `awk -F "\t" 'BEGIN {OFS="\t"} {print $3}' ${OUTDIR}/${sample}.barcodes.txt | awk '$1!=""' | awk '{print length}' | sort | uniq | wc -l` == 1 ]]; then
    cat ${OUTDIR}/${sample}.barcodes.txt | awk -F "\t" '{print $3}' | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | shuf -n 1000000 | awk -F "\t" 'BEGIN {OFS="\t"} {print ">id-" NR; print}' |
    weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} UMI" --stacks-per-line 100 > $TMPDIR/${sample}.UMIs.eps
    convert $TMPDIR/${sample}.UMIs.eps ${OUTDIR}/${sample}.UMIs.png
fi

echo
echo -n -e "${sample}\tNumber of total reads\t"
cat ${OUTDIR}/${sample}.barcodes.txt | wc -l
echo -n -e "${sample}\tNumber of total read barcodes\t"
awk -F "\t" '$1!="" {count+=1} END {print count}' ${OUTDIR}/${sample}.barcodes.txt


minUMIlength=5
if [[ `awk -v minUMIlength=${minUMIlength} -F "\t" 'BEGIN {OFS="\t"} length($3) >= minUMIlength {found=1} END {print found}' ${OUTDIR}/${sample}.barcodes.txt` == 1 ]]; then
    echo
    echo "Analyzing UMIs"
    
    echo
    echo "Barcode counts (before UMI)"
    cat ${OUTDIR}/${sample}.barcodes.txt | 
    awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | cut -f1,3 |
    #BUGBUG for non-UMI samples counts reads/barcode 
    sort -k1,1 -k2,2 | uniq -c | sort -k2,2 |
    #Format: BC, UMI, n
    awk 'BEGIN {OFS="\t"} {print $2, $3, $1}' > ${OUTDIR}/${sample}.barcode.counts.withUMI.txt
    awk '{print $1}' ${OUTDIR}/${sample}.barcode.counts.withUMI.txt | sort | uniq -c | awk 'BEGIN {OFS="\t"} {print $2, $1}' >  ${OUTDIR}/${sample}.barcode.counts.UMI.corrected.txt
    echo -n -e "${sample}\tNumber of unique barcodes+UMI\t"
    #TODO move to extract.py
    cat ${OUTDIR}/${sample}.barcode.counts.withUMI.txt | wc -l
    
    echo
    echo "Histogram of number of reads per barcode+UMI"
    cat ${OUTDIR}/${sample}.barcode.counts.withUMI.txt | cut -f3 | awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g
    
    echo
    echo "UMI lengths"
    cut -f2 ${OUTDIR}/${sample}.barcode.counts.withUMI.txt | 
    awk '{print length($0)}' | sort -g | uniq -c | sort -k2,2 | awk '$2!=0'
    #EROR 2017Feb16 bcfile="${OUTDIR}/${sample}.barcode.counts.withUMI.txt"
    bcfile="${OUTDIR}/${sample}.barcodes.txt"
else
    echo
    echo "No UMIs found"
    
    bcfile="${OUTDIR}/${sample}.barcodes.txt"
fi

echo
echo "Barcode counts"
cat ${bcfile} | cut -f1 |
sort -k1,1 | uniq -c | sort -k2,2 | awk '$2!=0' | awk 'BEGIN {OFS="\t"} {print $2, $1}' | awk -F "\t" '$1!=""' > ${OUTDIR}/${sample}.barcode.counts.txt

echo -n -e "${sample}\tNumber of unique barcodes\t"
#TODO move to extract.py
cat ${OUTDIR}/${sample}.barcode.counts.txt | wc -l

echo -n -e "${sample}\tNumber of unique barcodes with 10+ reads\t"
#TODO move to extract.py
cat ${OUTDIR}/${sample}.barcode.counts.txt | awk -v minReadCutoff=${minReadCutoff} -F "\t" 'BEGIN {OFS="\t"} $2>minReadCutoff' | wc -l

echo
echo "Histogram of number of reads per barcode"
cat ${OUTDIR}/${sample}.barcode.counts.txt | cut -f2 | awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g


echo
echo "Barcode lengths"
cut -f1 ${OUTDIR}/${sample}.barcode.counts.txt | 
awk '{print length($0)}' | sort -g | uniq -c | sort -k2,2 | awk '$2!=0'


echo
echo "Saturation curves"
numlines=`cat ${OUTDIR}/${sample}.barcodes.txt | wc -l`

for i in `echo {10,1000,10000,50000,100000,250000,500000,1000000,2000000,3000000,4000000,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000,50000000,55000000,60000000,70000000,80000000,90000000,100000000,110000000,120000000} | tr ' ' '\n' | gawk -v subset=$numlines '{if ($1<=subset) print}'`; do 
saturation=$(shuf -n $i ${OUTDIR}/${sample}.barcodes.txt | awk -F "\t" '{print $1}' | sort | uniq | wc -l); echo $i  $saturation 'minreads1' ;done > ${OUTDIR}/${sample}.Saturation_minReads1.txt

if [[ `awk -F "\t" 'BEGIN {OFS="\t"} length($3) > 4 {found=1} END {print found}' ${OUTDIR}/${sample}.barcodes.txt` == 1 ]]; then
    for i in `echo {10,1000,10000,50000,100000,250000,500000,1000000,2000000,3000000,4000000,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000,50000000,55000000,60000000,70000000,80000000,90000000,100000000,110000000,120000000} | tr ' ' '\n' | gawk -v subset=$numlines '{if ($1<=subset) print}'`; do 
    saturation=$(shuf -n $i ${OUTDIR}/${sample}.barcodes.txt | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $3}' | sort | uniq | awk -F "\t" '{print $1}' | sort | uniq -c | awk '{if ($1>=10) print $2}' | wc -l); echo $i $saturation 'minreads10' ;done > ${OUTDIR}/${sample}.Saturation_minReads10.txt
else
    for i in `echo {10,1000,10000,50000,100000,250000,500000,1000000,2000000,3000000,4000000,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000,50000000,55000000,60000000,70000000,80000000,90000000,100000000,110000000,120000000} | tr ' ' '\n' | gawk -v subset=$numlines '{if ($1<=subset) print}'`; do 
    saturation=$(shuf -n $i ${OUTDIR}/${sample}.barcodes.txt | awk -F "\t" '{print $1}' | sort | uniq -c | awk '{if ($1>=10) print $2}' | wc -l); echo $i $saturation 'minreads10' ;done > ${OUTDIR}/${sample}.Saturation_minReads10.txt
fi

echo
echo "Done!!!"
date
