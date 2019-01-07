#!/bin/bash
set -euo pipefail

module load trimmomatic/0.38
module load weblogo/3.5.0
module load ImageMagick
module load picard/1.140
module load FastQC/0.11.4
module load cutadapt/1.9.1

# the function "round()" was taken from 
# https://stempell.com/2009/08/rechnen-in-bash/
floor()
{
    set +u
    echo $(printf %.$2f $(echo "scale=0;$1/1" | bc))
    set -u
};

src=/vol/mauranolab/transposon/src


###Parse command line args
if [ "$#" -lt 9 ]; then
    echo "Wrong number of arguments"
    exit 1
fi

sample=$1
BCreadSeq=$2
R1trim=$3
R2trim=$4
bclen=$5
bcread=$6
plasmidSeq=$7
extractBCargs=$8

shift 8
basedir=$@
echo "Looking for files in ${basedir}"
f1=`find ${basedir}/ -maxdepth 1 -name "*_R1_*.fastq.gz"`
f2=`find ${basedir}/ -maxdepth 1 -name "*_R2_*.fastq.gz"`


OUTDIR=${sample}
mkdir -p $OUTDIR #note -p also makes $PREFIX here, required at the end


#export NSLOTS=1


echo "Processing sample ${sample}"
echo "Running on $HOSTNAME. Output to $OUTDIR"
date


###First process all reads together
echo "Weblogo of raw reads"
zcat -f $f1 | awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2' | shuf -n 1000000 | awk '{print ">id-" NR; print}' |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} R1 raw sequence" --stacks-per-line 100 > $TMPDIR/${sample}.R1.raw.eps
convert $TMPDIR/${sample}.R1.raw.eps $OUTDIR/${sample}.R1.raw.png

zcat -f $f2 | awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2' | shuf -n 1000000 | awk '{print ">id-" NR; print}' |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} R2 raw sequence" --stacks-per-line 100 > $TMPDIR/${sample}.R2.raw.eps
convert $TMPDIR/${sample}.R2.raw.eps $OUTDIR/${sample}.R2.raw.png


echo "Trimming/extracting UMI from R1 files $f1 and R2 files $f2"
#Get UMI from both reads even if only one is used for further analysis
if [[ "$R1trim" > "0" ]]; then
    echo "Trimming $R1trim bp from R1"
    bc1pattern=`printf 'N%.0s' $(seq 1 $R1trim)`
else
    bc1pattern="X"
fi

#BUGBUG really?
#Only extract UMI from R2 for RNA samples where bcread is R1
if [[ "$R2trim" > "0" && "$bcread" == "R1" ]]; then
    echo "Trimming $R2trim bp from R2"
    bc2pattern=`printf 'N%.0s' $(seq 1 $R2trim)`
    zcat -f $f2 > $TMPDIR/${sample}.R2.fastq
else
    bc2pattern="X"
    echo "Removing $R2trim bp from R2"
    cutadapt --quiet -u $R2trim $f2 > $TMPDIR/${sample}.R2.fastq
fi

if [[ "$bc1pattern" != "" ]]; then
    bc1pattern="--bc-pattern $bc1pattern"
fi
if [[ "$bc2pattern" != "" ]]; then
    bc2pattern=" --bc-pattern2 $bc2pattern"
fi
echo "umi_tools extract $bc1pattern $bc2pattern"

#BUGBUG really?
#zcat -f $f2 > $TMPDIR/${sample}.R2.fastq
#/home/mauram01/.local/bin/umi_tools extract --help
#BUGBUG umi_tools installed in python3.5 module missing matplotlib??? also .local by default is not group-readable
zcat -f $f1 | /home/mauram01/.local/bin/umi_tools extract $bc1pattern $bc2pattern --read2-in=$TMPDIR/${sample}.R2.fastq --read2-out=$TMPDIR/${sample}.R2.out.fastq -v 0 --log=$TMPDIR/${sample}.umi.log |
gzip -1 -c > $TMPDIR/${sample}.R1.fastq.gz
gzip -1 -c $TMPDIR/${sample}.R2.out.fastq > $TMPDIR/${sample}.R2.fastq.gz
gzip -9 -c $TMPDIR/${sample}.umi.log > $OUTDIR/${sample}.umi.log.gz


echo
echo "Filtering out reads with >75% G content"
$src/filterNextSeqReadsForPolyG.py --inputfileR1 $TMPDIR/${sample}.R1.fastq.gz --inputfileR2 $TMPDIR/${sample}.R2.fastq.gz --maxPolyG 75 --outputfileR1 $OUTDIR/${sample}.R1.fastq.gz --outputfileR2 $OUTDIR/${sample}.R2.fastq.gz


echo
echo -n "# reads R1: "
zcat -f $OUTDIR/${sample}.R1.fastq.gz | awk 'END {print NR/4}'
echo -n "# reads R2: "
zcat -f $OUTDIR/${sample}.R2.fastq.gz | awk 'END {print NR/4}'


echo
echo "FASTQC"
fastqc --outdir $OUTDIR $OUTDIR/${sample}.R1.fastq.gz
echo
fastqc --outdir $OUTDIR $OUTDIR/${sample}.R2.fastq.gz


echo
echo "Weblogo of raw UMI"
#BUGBUG why hardcoded to R2?
#Use tail to run through to end of file so zcat doesn't throw an error code
UMIlength=`zcat -f $OUTDIR/${sample}.R2.fastq.gz | awk 'BEGIN {OFS="\t"} NR % 4 == 1 {split($1, readname, "_"); print length(readname[2])}' | tail -1`
if [[ "${UMIlength}" > "0" ]]; then
    echo "Making weblogo of UMI"
    zcat -f $OUTDIR/${sample}.R2.fastq.gz | awk 'BEGIN {OFS="\t"} NR % 4 == 1 {split($1, readname, "_"); print readname[2]}' | shuf -n 1000000 | awk -F "\t" 'BEGIN {OFS="\t"} {print ">id-" NR; print}' | weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} R1 UMI sequence" --stacks-per-line 100 > $TMPDIR/${sample}.R1.raw.UMI.eps
    convert $TMPDIR/${sample}.R1.raw.UMI.eps $OUTDIR/${sample}.R1.raw.UMI.png
fi


echo
echo "Trimmomatic"
case "$bcread" in
R1)
    R1file="BC";
    R2file="plasmid";;
R2)
    R1file="plasmid";
    R2file="BC";;
*)
    echo "bcread $bcread must be either R1 or R2";
    exit 1;;
esac

trimmomaticBaseOpts="-threads $NSLOTS"
trimmomaticSteps="CROP:68 TOPHRED33"
java org.usadellab.trimmomatic.TrimmomaticPE $trimmomaticBaseOpts $OUTDIR/${sample}.R1.fastq.gz $OUTDIR/${sample}.R2.fastq.gz $OUTDIR/${sample}.trimmed.${R1file}.fastq.gz $OUTDIR/${sample}.trimmed.${R1file}.unpaired.fastq.gz $OUTDIR/${sample}.trimmed.${R2file}.fastq.gz $OUTDIR/${sample}.${R2file}.unpaired.fastq.gz $trimmomaticSteps


###Weblogo of processed reads
echo "Weblogo of processed reads"
zcat -f $OUTDIR/${sample}.trimmed.BC.fastq.gz | awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2' | shuf -n 1000000 | awk '{print ">id-" NR; print}' |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} BC processed sequence" --stacks-per-line 100 > $TMPDIR/${sample}.BC.processed.eps
convert $TMPDIR/${sample}.BC.processed.eps $OUTDIR/${sample}.BC.processed.png

zcat -f $OUTDIR/${sample}.trimmed.plasmid.fastq.gz | awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2' | shuf -n 1000000 | awk '{print ">id-" NR; print}' |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} plasmid processed sequence" --stacks-per-line 100 > $TMPDIR/${sample}.plasmid.processed.eps
convert $TMPDIR/${sample}.plasmid.processed.eps $OUTDIR/${sample}.plasmid.processed.png


###Finally submit jobs
numlines=`zcat $OUTDIR/${sample}.trimmed.BC.fastq.gz | wc -l`
chunksize=2000000 #Split fastq into 500,000 reads for deduplication (500,000 x 4)
numjobs=`echo "$numlines / $chunksize" | bc -l -q`
numjobs=$(floor $numjobs)
numjobs=`echo "$numjobs + 1" | bc -l -q`
#BUGBUG verify must be int. what to do with last chunk?
echo "$numlines lines to process in chunks of $chunksize"

echo
echo "Submitting $numjobs jobs"
qsub -S /bin/bash -t 1-${numjobs} -terse -j y --qos=full -N extract.${sample} -o ${sample} -b y "$src/extractBCcounts.sh ${sample} $BCreadSeq $bclen $chunksize $plasmidSeq $extractBCargs" | perl -pe 's/[^\d].+$//g;' > sgeid.map.${sample}
echo "Will merge $numjobs files"
bcfiles=`seq 1 $numjobs | xargs -L 1 -I {} echo -n "${sample}/${sample}.{}.barcodes.txt "`
echo -e "Will merge barcode files: $bcfiles\n"
cat <<EOF | qsub -S /bin/bash -terse -hold_jid `cat sgeid.map.${sample}` -j y --qos=full -N ${sample} -b y | perl -pe 's/[^\d].+$//g;' # > sgeid.merge.${sample}
set -e -o pipefail
echo "Merging barcodes"
cat $bcfiles > $OUTDIR/$sample.barcodes.preFilter.txt
#rm -f $bcfiles

$src/analyzeBCcounts.sh ${sample}
EOF

rm -f sgeid.${sample}


echo "Done!!!"
date
