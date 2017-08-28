#!/bin/bash
set -e -o pipefail
module load trimmomatic/0.33 weblogo/3.5.0 ImageMagick picard/1.140 FastQC/0.11.4 samtools/1.3.1 bwa/0.7.7 

NSLOTS=1

src=/home/maagj01/scratch/transposon/src

sample=$1
BCreadSeq=$2
bclen=$3
chunksize=$4
plasmidSeq=$5
OUTDIR=$sample
jobid=${SGE_TASK_ID}
#jobid=1
#TODO 

#TMPDIR=$OUTDIR/bamintermediate
#mkdir -p $TMPDIR
#echo "using $TMPDIR as TMPDIR"

echo "Analyzing barcodes"
#Will read SGE_TASK_ID independently
$src/extractBCcounts.sh $sample $BCreadSeq $bclen $chunksize $plasmidSeq 


#TODO f2 is the read containing the primer sequence and the genomic regions. 
f2=$OUTDIR/${sample}.plasmid.fastq.gz
sample="${sample}.$jobid"
#BUGBUG @RG wrong below--includes jobids
#TODO 
DS="BS00000A"


firstline=`echo "$chunksize * ($jobid - 1) + 1" | bc -l -q`
lastline=`echo "$chunksize * $jobid" | bc -l -q`
echo "Running on $HOSTNAME. Output to $OUTDIR. Jobid=$jobid (lines ${firstline}-${lastline})"




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


#Now set up for the iPCR-specific parts
##Trim primer before mapping to genome
R2primerlen=18
echo "Trimming $R2primerlen bp primer from R2"
zcat $f2 | awk -v firstline=$firstline -v lastline=$lastline 'NR>=firstline && NR<=lastline' | 
awk -F "\t" 'BEGIN {OFS="\t"} {if(NR % 4==1 ) {split($0, name, "_"); print name[1]} else {print}}' |
awk -v trim=$R2primerlen '{if(NR % 4==2 || NR % 4==0) {print substr($0, trim+1)} else if (NR % 4==1) {print $1} else {print}}' | gzip -9 -c > $TMPDIR/${sample}.genome.fastq.gz


date
echo "bwa aln $bwaAlnOpts $bwaIndex ..."
bwa aln $bwaAlnOpts $bwaIndex $TMPDIR/$sample.genome.fastq.gz > $OUTDIR/$sample.genome.sai


date
DS_nosuffix=`echo $DS | perl -pe 's/[A-Z]$//g;'`
bwaExtractOpts="-n 3 -r @RG\\tID:${sample}\\tLB:$DS\\tSM:${DS_nosuffix}"
extractcmd="samse $bwaExtractOpts $bwaIndex $OUTDIR/$sample.genome.sai $TMPDIR/$sample.genome.fastq.gz"
echo "Extracting"
echo -e "extractcmd=bwa $extractcmd | (...)"
bwa $extractcmd |
#No need to sort SE data
#samtools sort -@ $NSLOTS -O bam -T $OUTDIR/${sample}.sortbyname -l 1 -n - |
$src/filter_reads.py --max_mismatches $permittedMismatches - - |
samtools view -@ NSLOTS -1 - > $OUTDIR/$sample.bam


echo "Done!!!"
date
