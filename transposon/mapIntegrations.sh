#!/bin/bash
set -eu -o pipefail

#NSLOTS=1

src=$( dirname "${BASH_SOURCE[0]}" )


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

echo
echo "Analyzing barcodes"
#Will read SGE_TASK_ID independently
${src}/extractBCcounts.sh ${sample} $BCreadSeq $bclen $chunksize $plasmidSeq $extractBCargs


#TODO f2 is the read containing the primer sequence and the genomic regions.
f1=$OUTDIR/${sample}.BC.fastq.gz
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
permittedMismatches="0.1"
curGenome="hg38_noalt"
DS_nosuffix=`echo $DS | perl -pe 's/[A-Z]$//g;'`
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


###Now set up for the iPCR-specific parts
# Configurable parameters
minLengthForPEmapping=24
cutadaptR1args=""
cutadaptR2args=""
#Trim primer before mapping to genome
#R1 will get trimmed if PE mapping is selected
#BUGBUG hardcoded primer lengths. FIXED BY using sequence template as base length
R1primerlen=$(expr $(echo -n ${BCreadSeq} | wc -c) - 4)
R2primerlen=$(echo -n ${plasmidSeq} | wc -c)
minR1Len=$(expr ${R1primerlen} + ${minLengthForPEmapping})
# R1 read adapters
#BUGBUG move this so it isn't hard coded
altDpnSeq="GATCTTTGTCCAAACTCATCGAGCTCGG"
R2PrimerSeq=$(echo $plasmidSeq | tr "[ATGC]" "[TACG]" | rev)
# R2 read adapters
firstDpnRevSeq=$(echo ${BCreadSeq} | awk '{ print substr($0, 0, length($0)-4) }' | tr "[ATGCB]" "[TACGN]" | rev)
altDpnRevSeq=$(echo ${altDpnSeq} | tr "[ATGC]" "[TACG]" | rev)

echo
date
## Require BC read have R1primerlen + 20 bp length to run paired mapping
if [ "$(zcat -f $f1 | head -n 4000 | awk 'NR % 4 == 2 { sum += length($0) }; END { print 4 * sum/NR }')" -le $minR1Len ]; then

    echo "Trimming $R2primerlen bp primer from R2"
    zcat -f $f2 | 
    awk -v firstline=$firstline -v lastline=$lastline 'NR>=firstline && NR<=lastline' |
    cutadapt -j $NSLOTS -Z -o $TMPDIR/${sample}.genome.fastq.gz -u $R2primerlen ${cutadaptR2args} -

    # Single-end mapping using bwa aln
    bwaAlnOpts="-n ${permittedMismatches} -l 32 ${userAlnOptions} -t ${NSLOTS} -Y"
    
    echo "bwa aln ${bwaAlnOpts} ${bwaIndex} ..."
    bwa aln ${bwaAlnOpts} ${bwaIndex} $TMPDIR/${sample}.genome.fastq.gz > $TMPDIR/${sample}.genome.sai
    
    echo
    bwaAlnExtractOpts="-n 3 -r @RG\\tID:${sample}\\tLB:$DS\\tSM:${DS_nosuffix}"
    extractcmd="samse ${bwaAlnExtractOpts} ${bwaIndex} $TMPDIR/${sample}.genome.sai $TMPDIR/${sample}.genome.fastq.gz"
else
    # Paired-end mapping using bwa mem
    echo "Adaptive trimming of R1"
    zcat -f $f1 |
    awk -v firstline=$firstline -v lastline=$lastline 'NR>=firstline && NR<=lastline' |
    cutadapt -Z -j $NSLOTS -o $TMPDIR/${sample}.R1.fastq.gz -u ${R1primerlen} -a ${R2PrimerSeq} -g X${altDpnSeq} ${cutadaptR1args} -
    
    echo "Trimming $R2primerlen bp primer from R2"
    zcat -f $f2 | 
    awk -v firstline=$firstline -v lastline=$lastline 'NR>=firstline && NR<=lastline' |
    cutadapt -Z -j $NSLOTS -u 18 -a ${firstDpnRevSeq} -a ${altDpnRevSeq} ${cutadaptR2args} -o $TMPDIR/${sample}.genome.fastq.gz -
    
    bwaMemOptions="-Y -K 100000000"
    extractcmd="mem ${bwaMemOptions} -t ${NSLOTS} -R @RG\\tID:${sample}\\tLB:$DS\\tSM:${DS_nosuffix} ${bwaIndex} $TMPDIR/${sample}.R1.fastq.gz $TMPDIR/${sample}.genome.fastq.gz"
fi

echo "Extracting"
echo -e "extractcmd=bwa ${extractcmd} | (...)"
bwa ${extractcmd} |
#SE data technically doesn't need a sort
samtools sort -@ ${NSLOTS} -O bam -T $OUTDIR/${sample}.sortbyname -l 1 -n - |
#TODO UMIs don't seem to be in format filter_reads.py expects so they are not getting passed into bam
${src}/../dnase/filter_reads.py --reqFullyAligned --failUnwantedRefs --unwanted_refs_list "hap|random|^chrUn_|_alt$|scaffold|^C\d+|^pSB$|^pTR$" --max_mismatches ${permittedMismatches} --min_mapq ${minMAPQ} --min_insert_size 0 --max_insert_size 1000 --maxPermittedTrailingOverrun 2 - - |
samtools sort -@ $NSLOTS -m 1750M -O bam -T $TMPDIR/${sample}.sortbyname -l 1 > $OUTDIR/${sample}.bam


echo
echo "Done!!!"
date
