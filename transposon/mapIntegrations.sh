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
#curGenome="mm10"
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
#Trim primer before mapping to genome
#R1 will get trimmed if PE mapping is selected
R1primerlen=$(expr $(echo -n ${BCreadSeq} | wc -c) - 4)
R2primerlen=$(echo -n ${plasmidSeq} | wc -c)


echo
date

echo "Trimming $R2primerlen bp primer from R2"
zcat -f $f2 |
awk -v firstline=$firstline -v lastline=$lastline 'NR>=firstline && NR<=lastline' |
cutadapt -j $NSLOTS -u ${R2primerlen} - -o - |
#hardcoded fix for LP305-derived libraries
#cutadapt -j $NSLOTS -o - -g XTTATG - |
gzip -1 -c > $TMPDIR/${sample}.genome.fastq.gz

## Require BC read have 24 bp left after primer for PE mapping
minLengthForPEmapping=24
if [ "$(zcat -f $f1 | head -n 4000 | awk 'NR % 4 == 2 { sum += length($0) }; END { print 4 * sum/NR }')" -le $(expr $(echo -n ${BCreadSeq} | wc -c) + ${minLengthForPEmapping}) ]; then
    echo "Single-end mapping using bwa aln"
    
    bwaAlnOpts="-n ${permittedMismatches} -l 32 ${userAlnOptions} -t ${NSLOTS} -Y"
    echo "bwa aln ${bwaAlnOpts} ${bwaIndex} ..."
    bwa aln ${bwaAlnOpts} ${bwaIndex} $TMPDIR/${sample}.genome.fastq.gz > $TMPDIR/${sample}.genome.sai
    
    bwaAlnExtractOpts="-n 3 -r @RG\\tID:${sample}\\tLB:$DS\\tSM:${DS_nosuffix}"
    extractcmd="samse ${bwaAlnExtractOpts} ${bwaIndex} $TMPDIR/${sample}.genome.sai $TMPDIR/${sample}.genome.fastq.gz"
else
    echo "Paired-end mapping using bwa mem"
    
    #BUGBUG hard coded secondary trimming
    altDpnSeq="GATCTTTGTCCAAACTCATCGAGCTCGG"
    firstDpnRevSeq=$(echo ${BCreadSeq} | awk '{ print substr($0, 0, length($0)-4) }' | tr "[ATGCB]" "[TACGN]" | rev)
    altDpnRevSeq=$(echo ${altDpnSeq} | tr "[ATGC]" "[TACG]" | rev)
    
    echo "Adaptive trimming of R1"
    #Trim readthrough to the other read primer for long reads / short sites
    R2PrimerSeq=$(echo $plasmidSeq | tr "[ATGC]" "[TACG]" | rev)
    zcat -f $f1 |
    awk -v firstline=$firstline -v lastline=$lastline 'NR>=firstline && NR<=lastline' |
    cutadapt -j $NSLOTS  -u ${R1primerlen} -a ${R2PrimerSeq} -g X${altDpnSeq} - -Z -o $TMPDIR/${sample}.R1.fastq.gz
    
    
    echo
    echo "Trimming alternate Dpn site on R2"
    # R2 read adapters
    cutadapt -j $NSLOTS -a ${firstDpnRevSeq} -a ${altDpnRevSeq} $TMPDIR/${sample}.genome.fastq.gz -Z -o $TMPDIR/${sample}.genome.fastq.gz.new && mv $TMPDIR/${sample}.genome.fastq.gz.new $TMPDIR/${sample}.genome.fastq.gz
    
    bwaMemOptions="-Y -K 100000000"
    extractcmd="mem ${bwaMemOptions} -t ${NSLOTS} -R @RG\\tID:${sample}\\tLB:$DS\\tSM:${DS_nosuffix} ${bwaIndex} $TMPDIR/${sample}.R1.fastq.gz $TMPDIR/${sample}.genome.fastq.gz"
fi

echo
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
