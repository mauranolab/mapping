#!/bin/bash
set -e


NSLOTS=1


sample=$1
f1=$2/*_R1_*.fastq.gz
f2=$2/*_R2_*.fastq.gz

amplicon=$2

#clip before anything (i.e. for NNNN used to create diversity)
R1trim=$3
R2trim=$4

#BUGBUG WTF
#R2primerlen=18
#R2trim=`echo $R2trim + $R2primerlen | bc -l -q`

bclen=$5

bcread=$6

shift 6
basedir=$@
echo "Looking for files in $basedir"
f1=`find $basedir/ -name "*_R1_*.fastq.gz"`
f2=`find $basedir/ -name "*_R2_*.fastq.gz"`



#trim after alignment (i.e. for sequencing errors)
trim=0

minNumOccurrences=1


OUTDIR=$sample
mkdir -p $OUTDIR #note -p also makes $PREFIX here, required at the end


echo "Running on $HOSTNAME. Using $OUTDIR as tmp"
date


zcat -f $f1 | awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2 {print ">id-" NR; print}' | head -500000 | tail -100000 |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} R1 raw sequence" --stacks-per-line 100 > $TMPDIR/${sample}.R1.raw.eps
convert $TMPDIR/${sample}.R1.raw.eps $OUTDIR/${sample}.R1.raw.png

zcat -f $f2 | awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2 {print ">id-" NR; print}' | head -500000 | tail -100000 |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} R2 raw sequence" --stacks-per-line 100 > $TMPDIR/${sample}.R2.raw.eps
convert $TMPDIR/${sample}.R2.raw.eps $OUTDIR/${sample}.R2.raw.png


if [[ "$bcread" == "R1" ]]; then
       echo "Trimming $R1trim bp from R1 files $f1"
       zcat -f $f1 | awk -v trim=$R1trim '{if(NR % 4==2 || NR % 4==0) {print substr($0, trim+1)} else {print}}' | gzip -c > $OUTDIR/$sample.R1.fastq.gz
       echo -n "# tags R1: "
       zcat $OUTDIR/$sample.R1.fastq.gz | awk 'END {print NR/4}'


elif [[ "$bcread" == "R2" ]]; then
       echo "Trimming $R2trim bp from R2 files $f2"
       zcat -f $f2 | awk -v trim=$R2trim '{if(NR % 4==2 || NR % 4==0) {print substr($0, trim+1)} else {print}}' | gzip -c > $OUTDIR/$sample.R2.fastq.gz
       echo -n "# tags R2: "
       zcat $OUTDIR/$sample.R2.fastq.gz | awk 'END {print NR/4}'
else
       echo "Can't understand read $bcread"
       exit 1
fi

echo "FASTQC"
fastqc --outdir $OUTDIR $OUTDIR/$sample.$bcread.fastq.gz


echo
echo "Trimmomatic"
#Don't care about second read
#TODO verify that second read matches plasmid sequence?
trimmomaticBaseOpts="-threads $NSLOTS -trimlog $OUTDIR/${sample}.trim.log.txt"
trimmomaticSteps="TOPHRED33 AVGQUAL:30"

java org.usadellab.trimmomatic.TrimmomaticSE $trimmomaticBaseOpts $OUTDIR/$sample.$bcread.fastq.gz $OUTDIR/$sample.trimmed.$bcread.fastq.gz $trimmomaticSteps


echo
echo "Extracting barcodes from $bcread"
date
zcat -f $OUTDIR/$sample.trimmed.$bcread.fastq.gz | ~mauram01/scratch/transposon/src/extractBarcode.py - --referenceSeq $amplicon --bclen $bclen > $TMPDIR/$sample.barcodes.txt > $OUTDIR/$sample.barcodes.txt
# --align
#Format: BC, Readname
date

cut -f1 $OUTDIR/$sample.barcodes.txt | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' |
awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2 {print ">id-" NR; print}' | head -500000 | tail -100000 |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} barcodes" --stacks-per-line 100 > $TMPDIR/${sample}.barcodes.eps
convert $TMPDIR/${sample}.barcodes.eps $OUTDIR/${sample}.barcodes.png

cut -f1 $OUTDIR/$sample.barcodes.txt | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' |
sort -g | uniq -c | sort -k2,2 | awk 'BEGIN {OFS="\t"} {print $2, $1}' > $OUTDIR/$sample.barcode.counts.txt


echo
echo "Barcode lengths"
cut -f1 $OUTDIR/$sample.barcode.counts.txt | awk '{print length($0)}' | sort -g | uniq -c | sort -k2,2


echo
echo "Number of unique barcodes"
#TODO move to extract.py
cat $OUTDIR/$sample.barcode.counts.txt | wc -l


echo
echo "Histogram of number of reads per barcode"
cat $OUTDIR/$sample.barcode.counts.txt | cut -f2 | awk -v cutoff=20 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g


echo "Done!!!"
date
