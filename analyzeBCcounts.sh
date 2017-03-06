#!/bin/bash
set -e # -o pipefail

sample=$1

OUTDIR=$sample

#One input file:
#$OUTDIR/$sample.barcodes.preFilter.txt


echo "Analyzing data for $sample"

~/scratch/transposon/src/removeOverrepresentedBCs.py  --col 1 --freq 0.01 -o $OUTDIR/$sample.barcodes.txt $OUTDIR/$sample.barcodes.preFilter.txt 


###Overall summary statistics
shuf -n 1000000 $OUTDIR/$sample.barcodes.txt |cut -f1 | awk -F "\t" 'BEGIN {OFS="\t"} $1!="" {print ">id-" NR; print $1}' |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} barcodes" --stacks-per-line 100 > $TMPDIR/${sample}.barcodes.eps
convert $TMPDIR/${sample}.barcodes.eps $OUTDIR/${sample}.barcodes.png

###UMIs weblogo
if [[ `awk -F "\t" 'BEGIN {OFS="\t"} length($3) >= 1 {found=1} END {print found}' $OUTDIR/$sample.barcodes.txt` == 1 ]]; then
    shuf -n 1000000 $OUTDIR/$sample.barcodes.txt| awk -F'\t' '{print $3}'| awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | 
    awk -F "\t" 'BEGIN {OFS="\t"} {print ">id-" NR; print}' |
    weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${i} UMI" --stacks-per-line 100 > $TMPDIR/${sample}.UMIs.eps
    convert $TMPDIR/${sample}.UMIs.eps $OUTDIR/${sample}.UMIs.png
fi

echo
echo -n -e "$sample\tNumber of total reads\t"
cat $OUTDIR/$sample.barcodes.txt | wc -l
echo -n -e "$sample\tNumber of total read barcodes\t"
awk -F "\t" '$1!="" {count+=1} END {print count}' $OUTDIR/$sample.barcodes.txt


if [[ `awk -F "\t" 'BEGIN {OFS="\t"} length($3) > 4 {found=1} END {print found}' $OUTDIR/$sample.barcodes.txt` == 1 ]]; then
       echo
       echo "Analyzing UMIs"
       
       echo
       echo "Barcode counts (before UMI)"
       cat $OUTDIR/$sample.barcodes.txt | 
       awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | cut -f1,3 |
       #BUGBUG for non-UMI samples counts reads/barcode 
       sort -k1,1 -k2,2 | uniq -c | sort -k2,2 |
       #Format: BC, UMI, n
       awk 'BEGIN {OFS="\t"} {print $2, $3, $1}' > $OUTDIR/$sample.barcode.counts.withUMI.txt

       echo
       echo -n -e "$sample\tNumber of unique barcodes+UMI\t"
       #TODO move to extract.py
       cat $OUTDIR/$sample.barcode.counts.withUMI.txt | wc -l

       echo
       echo "Histogram of number of reads per barcode+UMI"
       cat $OUTDIR/$sample.barcode.counts.withUMI.txt | cut -f3 | awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g

       echo
       echo "UMI lengths"
       cut -f2 $OUTDIR/$sample.barcode.counts.withUMI.txt | 
       awk '{print length($0)}' | sort -g | uniq -c | sort -k2,2|awk '$2!=0'
       #EROR 2017Feb16 bcfile="$OUTDIR/$sample.barcode.counts.withUMI.txt"
       bcfile="$OUTDIR/$sample.barcodes.txt"
else
       echo
       echo "No UMIs found"
       
       bcfile="$OUTDIR/$sample.barcodes.txt"
fi

echo "Barcode counts"
cat $bcfile | cut -f1 |
sort -k1,1 | uniq -c | sort -k2,2 | awk '$2!=0'| awk 'BEGIN {OFS="\t"} {print $2, $1}'| awk -F'\t' '$1!=""' > $OUTDIR/$sample.barcode.counts.txt

echo
echo -n -e "$sample\tNumber of unique barcodes\t"
#TODO move to extract.py
cat $OUTDIR/$sample.barcode.counts.txt | wc -l

echo
echo "Histogram of number of reads per barcode"
cat $OUTDIR/$sample.barcode.counts.txt | cut -f2 | awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g


echo
echo "Barcode lengths"
cut -f1 $OUTDIR/$sample.barcode.counts.txt | 
awk '{print length($0)}' | sort -g | uniq -c | sort -k2,2|awk '$2!=0'


echo
echo "Done!!!"
date
