#!/bin/bash
#BUGBUG weblogo breaking pipefail
set -e # -o pipefail
module load trimmomatic/0.33 weblogo/3.5.0 ImageMagick picard/1.140 FastQC/0.11.4
# the function "round()" was taken from 
# https://stempell.com/2009/08/rechnen-in-bash/

floor()
{
       echo $(printf %.$2f $(echo "scale=0;$1/1" | bc))
};

sample=$1
#amplicon=$2
bcread=$2
OUTDIR=${sample}
numlines=`zcat $OUTDIR/${sample}.trimmed.$bcread.fastq.gz | wc -l`
chunksize=2000000 #Split fastq into 500,000 reads for deduplication (500,000 x 4)
numjobs=`echo "$numlines / $chunksize" | bc -l -q`
numjobs=$(floor $numjobs)
numjobs=$(floor $numjobs)



echo "Will merge $numjobs files"
bcfiles=`seq 1 $numjobs | xargs -L 1 -I {} echo -n "${sample}/${sample}.{}.barcodes.deduped.txt "`
echo -e "Will merge barcode files: $bcfiles\n"

cat <<EOF | qsub -S /bin/bash -terse  -j y -N ${sample} -b y | perl -pe 's/[^\d].+$//g;' # > sgeid.merge.${sample}
set -e
echo Merging barcodes
cat $bcfiles > $OUTDIR/$sample.barcodes.deduped.txt
#rm -f $bcfiles

bash ~/scratch/transposon/src/analyzeBCcounts_DEDUP.sh ${sample}
EOF
