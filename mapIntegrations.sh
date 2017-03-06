#!/bin/bash
set -e -o pipefail
module load trimmomatic/0.33 weblogo/3.5.0 ImageMagick picard/1.140 FastQC/0.11.4 samtools/1.3.1 bwa/0.7.7 

NSLOTS=1


sample=$1
amplicon=$2
bclen=$3
chunksize=$4


jobid=${SGE_TASK_ID}
#jobid=1


echo "Analyzing barcodes"
#Will read SGE_TASK_ID independently
~/scratch/transposon/src/extractBCcounts.sh $sample $amplicon $bclen R1 $chunksize


#Now set up for the iPCR-specific parts
OUTDIR=$sample
f2=$OUTDIR/$sample.trimmed.R2.fastq.gz
sample="${sample}.$jobid"
#BUGBUG @RG wrong below--includes jobids
#TODO
DS="BS00000A"


firstline=`echo "$chunksize * ($jobid - 1) + 1" | bc -l -q`
lastline=`echo "$chunksize * $jobid" | bc -l -q`
echo "Running on $HOSTNAME. Output to $OUTDIR. Jobid=$jobid (lines ${firstline}-${lastline})"



#TODO trim primer off read



echo
echo "Aligning reads"
permittedMismatches=2
curStrain="hg38"
bwaAlnOpts="-n $permittedMismatches -l 32 $userAlnOptions -t $NSLOTS -Y"

echo "Will map to strains $curStrain"


#NB am losing about 15" to load index when submit multiple jobs
bwaIndexBase=/vol/isg/annotation/bwaIndex

echo
echo "Mapping to genome for strain $curStrain"
case "$curStrain" in
hg19)
       bwaIndex=$bwaIndexBase/hg19all/hg19all;;
hg38)
       bwaIndex=$bwaIndexBase/hg38all/hg38all;;
mm10)
       bwaIndex=$bwaIndexBase/mm10all/mm10all;;
*)
       echo "Don't recognize strain $curStrain";
       exit 3;;
esac

zcat $f2 | awk -v firstline=$firstline -v lastline=$lastline 'NR>=firstline && NR<=lastline' | gzip -c > $TMPDIR/$sample.trimmed.R2.fastq.gz


date
echo "bwa aln $bwaAlnOpts $bwaIndex ..."
bwa aln $bwaAlnOpts $bwaIndex $TMPDIR/$sample.trimmed.R2.fastq.gz > $OUTDIR/$sample.R2.sai


date
DS_nosuffix=`echo $DS | perl -pe 's/[A-Z]$//g;'`
bwaExtractOpts="-n 3 -r @RG\\tID:${sample}\\tLB:$DS\\tSM:${DS_nosuffix}"
extractcmd="samse $bwaExtractOpts $bwaIndex $OUTDIR/$sample.R2.sai $TMPDIR/$sample.trimmed.R2.fastq.gz"
echo "Extracting"
echo -e "extractcmd=bwa $extractcmd | (...)"
bwa $extractcmd |
#No need to sort SE data
#samtools sort -@ $NSLOTS -O bam -T $OUTDIR/${sample}.sortbyname -l 1 -n - |
python /vol/mauranolab/maagj01/transposon/src/filter_reads.py --max_mismatches $permittedMismatches - - |
samtools view -@ NSLOTS -1 - > $OUTDIR/$sample.bam


echo "Done!!!"
date
