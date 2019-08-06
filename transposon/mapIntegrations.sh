#!/bin/bash
set -eu -o pipefail

#NSLOTS=1

src=/vol/mauranolab/mapped/src/transposon


###Parse command line args
if [ "$#" -ne 6 ]; then
    echo "Wrong number of arguments"
    exit 1
fi

sample=$1
BCreadSeq=$2
bclen=$3
chunksize=$4
plasmidSeq=$5
extractBCargs=$6

OUTDIR=${sample}
jobid=${SGE_TASK_ID}
#jobid=1

#TMPDIR=$OUTDIR/bamintermediate
#mkdir -p $TMPDIR
#echo "using $TMPDIR as TMPDIR"

echo "Analyzing barcodes"
#Will read SGE_TASK_ID independently
${src}/extractBCcounts.sh ${sample} $BCreadSeq $bclen $chunksize $plasmidSeq $extractBCargs


#TODO f2 is the read containing the primer sequence and the genomic regions. 
f2=$OUTDIR/${sample}.plasmid.fastq.gz
sample="${sample}.$jobid"
#BUGBUG @RG wrong below--includes jobids
#TODO 
DS="BS00000A"


firstline=`echo "$chunksize * ($jobid - 1) + 1" | bc -l -q`
lastline=`echo "$chunksize * $jobid" | bc -l -q`
echo "Running on $HOSTNAME. Output to $OUTDIR"
echo "Jobid=$jobid (lines ${firstline}-${lastline})"



echo
echo "Aligning reads"
userAlnOptions=""
permittedMismatches=2
curGenome="hg38_noalt"
bwaAlnOpts="-n ${permittedMismatches} -l 32 ${userAlnOptions} -t ${NSLOTS} -Y"
minMAPQ=10

echo
echo "Mapping to reference ${curGenome}"
bwaIndexBase=/vol/isg/annotation/bwaIndex
case "${curGenome}" in
hg38_noalt)
    bwaIndex=${bwaIndexBase}/hg38_noalt_transposon/hg38_noalt_transposon;;
mm10)
    bwaIndex=${bwaIndexBase}/mm10_no_alt_analysis_set/mm10_no_alt_analysis_set;;
*)
    echo "Don't recognize genome ${curGenome}";
    exit 3;;
esac


#Now set up for the iPCR-specific parts
##Trim primer before mapping to genome
#BUGBUG hardcoded primer length
R2primerlen=18
echo "Trimming $R2primerlen bp primer from R2"
zcat -f $f2 | awk -v firstline=$firstline -v lastline=$lastline 'NR>=firstline && NR<=lastline' |
#awk -F "\t" 'BEGIN {OFS="\t"} {if(NR % 4==1 ) {split($0, name, "_"); print name[1]} else {print}}' |
awk -v trim=$R2primerlen '{if(NR % 4==2 || NR % 4==0) {print substr($0, trim+1)} else if (NR % 4==1) {print $1} else {print}}' | pigz -p ${NSLOTS} -c -1 > $TMPDIR/${sample}.genome.fastq.gz


date
echo "bwa aln ${bwaAlnOpts} ${bwaIndex} ..."
bwa aln ${bwaAlnOpts} ${bwaIndex} $TMPDIR/${sample}.genome.fastq.gz > $TMPDIR/${sample}.genome.sai


echo
date
DS_nosuffix=`echo $DS | perl -pe 's/[A-Z]$//g;'`
bwaExtractOpts="-n 3 -r @RG\\tID:${sample}\\tLB:$DS\\tSM:${DS_nosuffix}"
extractcmd="samse ${bwaExtractOpts} ${bwaIndex} $TMPDIR/${sample}.genome.sai $TMPDIR/${sample}.genome.fastq.gz"
echo "Extracting"
echo -e "extractcmd=bwa ${extractcmd} | (...)"
bwa ${extractcmd} |
#No need to sort SE data
#samtools sort -@ ${NSLOTS} -O bam -T $OUTDIR/${sample}.sortbyname -l 1 -n - |
#TODO UMIs don't seem to be in format filter_reads.py expects so they are not getting passed into bam
${src}/../dnase/filter_reads.py --reqFullyAligned --failUnwantedRefs --unwanted_refs_list "hap|random|^chrUn_|_alt$|scaffold|^C\d+|^pSB$|^pTR$" --max_mismatches ${permittedMismatches} --min_mapq ${minMAPQ} - - |
samtools sort -@ $NSLOTS -m 1750M -O bam -T $TMPDIR/${sample}.sortbyname -l 1 > $OUTDIR/${sample}.bam


echo
echo "Done!!!"
date
