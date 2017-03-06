#!/bin/bash
#BUGBUG weblogo breaking pipefail
set -e # -o pipefail

module load trimmomatic/0.33 weblogo/3.5.0 ImageMagick picard/1.140 FastQC/0.11.4 samtools/1.3.1 bwa/0.7.7 
# the function "round()" was taken from 
# https://stempell.com/2009/08/rechnen-in-bash/

floor()
{
       echo $(printf %.$2f $(echo "scale=0;$1/1" | bc))
};


sample=$1

amplicon=$2

#clip before anything (i.e. for NNNN used to create diversity)
R1trim=$3
R2trim=$4
bclen=$5

shift 5
basedir=$@
echo "Looking for files in $basedir"
f1=`find $basedir/ -name "*_R1_*.fastq.gz"`
f2=`find $basedir/ -name "*_R2_*.fastq.gz"`


OUTDIR=${sample}
mkdir -p $OUTDIR #note -p also makes $PREFIX here, required at the end


export NSLOTS=1


echo "Processing sample ${sample}"
echo "Running on $HOSTNAME. Output to $OUTDIR"
date


###First process all reads together
echo "Weblogo of raw reads"
zcat -f $f1 | awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2 {print ">id-" NR; print}' | head -500000 | tail -100000 |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} R1 raw sequence" --stacks-per-line 100 > $TMPDIR/${sample}.R1.raw.eps
convert $TMPDIR/${sample}.R1.raw.eps $OUTDIR/${sample}.R1.raw.png

zcat -f $f2 | awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2 {print ">id-" NR; print}' | head -500000 | tail -100000 |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} R2 raw sequence" --stacks-per-line 100 > $TMPDIR/${sample}.R2.raw.eps
convert $TMPDIR/${sample}.R2.raw.eps $OUTDIR/${sample}.R2.raw.png


echo "Trimming/extracting UMI from R1 files $f1 and R2 files $f2"
#zcat -f $f1 | awk -v trim=$R1trim '{if(NR % 4==2 || NR % 4==0) {print substr($0, trim+1)} else if (NR % 4==1) {print $1} else {print}}' | gzip -c > $OUTDIR/${sample}.R1.fastq.gz
#zcat -f $f2 | awk -v trim=$R2trim '{if(NR % 4==2 || NR % 4==0) {print substr($0, trim+1)} else if (NR % 4==1) {print $1} else {print}}' | gzip -c > $OUTDIR/${sample}.R2.fastq.gz

if [[ "$R1trim" > "0" ]]; then
       echo "Trimming $R1trim bp from R1"
       bc1pattern="--bc-pattern "`printf 'N%.0s' $(seq 1 $R1trim)`
else
       bc1pattern="--bc-pattern X"
fi
if [[ "$R2trim" > "0" ]]; then
       echo "Trimming $R2trim bp from R2"
       bc2pattern=" --bc-pattern2 "`printf 'N%.0s' $(seq 1 $R2trim)`
else
       bc2pattern=" --bc-pattern2 X"
fi
echo "umi_tools extract $bc1pattern $bc2pattern"

zcat -f $f2 > $TMPDIR/${sample}.R2.fastq
#/home/mauram01/.local/bin/umi_tools extract --help
zcat -f $f1 | /home/maagj01/.local/bin/umi_tools extract $bc1pattern $bc2pattern --read2-in=$TMPDIR/${sample}.R2.fastq --read2-out=$TMPDIR/${sample}.R2.out.fastq -v 0 --log=$TMPDIR/${sample}.umi.log |
gzip -9 -c > $OUTDIR/${sample}.R1.fastq.gz
gzip -9 -c $TMPDIR/${sample}.umi.log > $OUTDIR/${sample}.umi.log.gz

echo -n "# tags R1: "
zcat -f $OUTDIR/${sample}.R1.fastq.gz | awk 'END {print NR/4}'
echo -n "# tags R2: "
zcat -f $TMPDIR/${sample}.R2.out.fastq | awk 'END {print NR/4}'

#Weblogo after UMI extraction
echo "Weblogo of raw UMI"
zcat -f $OUTDIR/${sample}.R1.fastq.gz| awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2 {print ">id-" NR; print}' | head -500000 | tail -100000 |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} R1 UMI sequence" --stacks-per-line 100 > $TMPDIR/${sample}.R1.raw.UMI.eps
convert $TMPDIR/${sample}.R1.raw.UMI.eps $OUTDIR/${sample}.R1.raw.UMI.png




echo "FASTQC"
fastqc --outdir $OUTDIR $OUTDIR/${sample}.R1.fastq.gz
echo
fastqc --outdir $OUTDIR $TMPDIR/${sample}.R2.out.fastq






#Trim primer before mapping to genome
R2primerlen=18
echo "Trimming $R2primerlen bp primer from R2"
cat $TMPDIR/${sample}.R2.out.fastq | 
#Remove UMI at the same time
awk -F "\t" 'BEGIN {OFS="\t"} {if(NR % 4==1 ) {split($0, name, "_"); print name[1]} else {print}}' |
awk -v trim=$R2primerlen '{if(NR % 4==2 || NR % 4==0) {print substr($0, trim+1)} else if (NR % 4==1) {print $1} else {print}}' | gzip -9 -c > $OUTDIR/${sample}.R2.fastq.gz



#####
#Remove all poly-G reads with >75% G
#####
echo 'Removing polyG reads'
cp $OUTDIR/${sample}.R1.fastq.gz $TMPDIR/${sample}.R1.fastq.gz
cp $OUTDIR/${sample}.R2.fastq.gz $TMPDIR/${sample}.R2.fastq.gz
echo 'Gunzipping fastq'
gunzip -f $TMPDIR/${sample}.R1.fastq.gz
gunzip -f $TMPDIR/${sample}.R2.fastq.gz
echo "Filtering out polyG reads"
python /home/maagj01/scratch/transposon/src/filterPolyGseq.py --inputR1 $TMPDIR/${sample}.R1.fastq --inputR2 $TMPDIR/${sample}.R2.fastq --maxPolyG 75 -oR1 $OUTDIR/${sample}.R1.fastq -oR2 $OUTDIR/${sample}.R2.fastq
echo 'Gzipping'
gzip -f $OUTDIR/${sample}.R1.fastq
gzip -f $OUTDIR/${sample}.R2.fastq 




echo
echo "Trimmomatic"
#NB redundant right now
#
#TODO I think these match
#illuminaAdapters=/cm/shared/apps/trimmomatic/0.33/adapters/TruSeq3-PE-2.fa
#seedmis=2
##Pretty much anything below 10 works
#PEthresh=5
#SEthresh=5
#mintrim=1
#keepReverseReads=true
trimmomaticBaseOpts="-threads $NSLOTS -trimlog $OUTDIR/${sample}.trim.log.txt"
trimmomaticSteps="TOPHRED33"
#AVGQUAL:25 TRAILING:20 MINLEN:27
#was: ILLUMINACLIP:$illuminaAdapters:$seedmis:$PEthresh:$SEthresh:$mintrim:$keepReverseReads AVGQUAL:35 TRAILING:35

java org.usadellab.trimmomatic.TrimmomaticPE $trimmomaticBaseOpts $OUTDIR/${sample}.R1.fastq.gz $OUTDIR/${sample}.R2.fastq.gz $OUTDIR/${sample}.trimmed.R1.fastq.gz $OUTDIR/${sample}.trimmed.R1.unpaired.fastq.gz $OUTDIR/${sample}.trimmed.R2.fastq.gz $OUTDIR/${sample}.trimmed.R2.unpaired.fastq.gz $trimmomaticSteps


numlines=`zcat $OUTDIR/${sample}.trimmed.R1.fastq.gz | wc -l`
chunksize=2000000
numjobs=`echo "$numlines / $chunksize" | bc -l -q`
numjobs=$(floor $numjobs)
numjobs=`echo "$numjobs + 1" | bc -l -q`
#BUGBUG verify must be int. what to do with last chunk?
echo "$numlines lines to process in chunks of $chunksize"

echo
echo "Submitting $numjobs jobs"
qsub -S /bin/bash -q testq -t 1-${numjobs} -terse -j y -N map.${sample} -o ${sample} -b y "~/scratch/transposon/src/mapIntegrations.sh ${sample} $amplicon $bclen $chunksize" | perl -pe 's/[^\d].+$//g;' > sgeid.map.${sample}

echo "Will merge $numjobs files"
bcfiles=`seq 1 $numjobs | xargs -L 1 -I {} echo -n "${sample}/${sample}.{}.barcodes.txt "`
echo -e "Will merge barcode files: $bcfiles\n"
bamfiles=`seq 1 $numjobs | xargs -L 1 -I {} echo -n "${sample}/${sample}.{}.bam "`
echo -e "Will merge bamfiles files: $bamfiles\n"
cat <<EOF | qsub -S /bin/bash -terse -hold_jid `cat sgeid.map.${sample}` -j y -N ${sample} -b y | perl -pe 's/[^\d].+$//g;' # > sgeid.merge.${sample}
set -e -o pipefail
echo Merging barcodes
cat $bcfiles > $OUTDIR/$sample.barcodes.txt
#rm -f $bcfiles

echo Merging bam files
samtools merge -f -l 9 $OUTDIR/$sample.bam $bamfiles
#TODO sort?
#samtools index $OUTDIR/$sample.bam
#rm -f $bamfiles

bash ~/scratch/transposon/src/analyzeIntegrations.sh ${sample}
EOF
rm -f sgeid.map.${sample}


echo "Done!!!"
date
