#!/bin/bash
set -euo pipefail

#Limit thread usage by python processes using OPENBLAS (esp. scipy). Set here and will be inherited by spawned jobs
#https://stackoverflow.com/questions/51256738/multiple-instances-of-python-running-simultaneously-limited-to-35
export OPENBLAS_NUM_THREADS=1

module load miller
module load pigz
module load weblogo/3.5.0
module load python/3.8.1
module load ImageMagick
module load picard/1.140
module load FastQC/0.11.4
module load cutadapt/1.9.1
module load samtools/1.9
module load bwa/0.7.15
module load bedops/2.4.35
module load bedtools/2.25
module load ucsckentutils/12152017

#NB python packages swalign and leven are required to be installed


# the function "round()" was taken from 
# https://stempell.com/2009/08/rechnen-in-bash/
floor()
{
    set +u
    echo $(printf %.$2f $(echo "scale=0;$1/1" | bc))
    set -u
};

src=/vol/mauranolab/mapped/src/transposon

#Hardcoded right now rather than as parameters like dnase pipeline
runPreprocess=1
runExtract=1
runMerge=1

###Parse command line args
if [ "$#" -lt 9 ]; then
    echo "ERROR submit: Wrong number of arguments"
    exit 1
fi

sample=$1
sampleType=$2
R1trim=$3
R2trim=$4
bcread=$5
bclen=$6
BCreadSeq=$7
plasmidSeq=$8
extractBCargs=$9

shift 9
basedir=$@
echo "Looking for files in ${basedir}"
f1=`find ${basedir}/ -maxdepth 1 -name "*_R1_*.fastq.gz" | sort`
f2=`find ${basedir}/ -maxdepth 1 -name "*_R2_*.fastq.gz" | sort`


if [[ "${sampleType}" != "DNA" ]] && [[ "${sampleType}" != "10xRNA" ]] && [[ "${sampleType}" != "RNA" ]] && [[ "${sampleType}" != "iPCR" ]]; then
    echo "ERROR submit: unknown sample type ${sampleType}"
    exit 2
fi


OUTDIR=${sample}
mkdir -p $OUTDIR


#export NSLOTS=1


echo "Processing ${sampleType} sample ${sample}"
echo "Running on $HOSTNAME. Output to $OUTDIR"
date


#TODO could refactor as separate preprocess job but need to know the number of total lines in trimmed fastq. Would take:
#floor()
#sample=$1
#R1trim=$2
#R2trim=$3
#bcread=$4
#src=$5
if [ ${runPreprocess} -eq 1 ]; then
    ###First process all reads together
    echo
    echo "Extracting BCs from ${bcread}"
    
    case "${bcread}" in
    R1)
        R1file="BC";
        R2file="plasmid";;
    R2)
        R1file="plasmid";
        R2file="BC";;
    *)
        echo "ERROR submit: bcread ${bcread} must be either R1 or R2";
        exit 3;;
    esac
    
    
    echo "Weblogo of raw reads"
    zcat -f ${f1} | awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2' | shuf -n 1000000 | awk '{print ">id-" NR; print}' |
    weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} R1 raw sequence" --stacks-per-line 100 > $TMPDIR/${sample}.R1.raw.eps
    convert $TMPDIR/${sample}.R1.raw.eps $OUTDIR/${sample}.R1.raw.png
    
    zcat -f ${f2} | awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2' | shuf -n 1000000 | awk '{print ">id-" NR; print}' |
    weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} R2 raw sequence" --stacks-per-line 100 > $TMPDIR/${sample}.R2.raw.eps
    convert $TMPDIR/${sample}.R2.raw.eps $OUTDIR/${sample}.R2.raw.png
    
    
    echo
    echo "Filtering out reads with >75% G content"
    date
    zcat -f ${f1} | gzip -1 > $TMPDIR/${sample}.R1.fastq.gz
    zcat -f ${f2} | gzip -1 > $TMPDIR/${sample}.R2.fastq.gz
    ${src}/../dnase/filterNextSeqReadsForPolyG.py --inputfileR1 $TMPDIR/${sample}.R1.fastq.gz --inputfileR2 $TMPDIR/${sample}.R2.fastq.gz --maxPolyG 75 --outputfileR1 $TMPDIR/${sample}.filtered.R1.fastq.gz --outputfileR2 $TMPDIR/${sample}.filtered.R2.fastq.gz
    
    
    echo
    echo -n "# reads R1: "
    zcat -f $TMPDIR/${sample}.filtered.R1.fastq.gz | awk 'END {print NR/4}'
    echo -n "# reads R2: "
    zcat -f $TMPDIR/${sample}.filtered.R2.fastq.gz | awk 'END {print NR/4}'
    
    
    echo
    echo "FASTQC"
    date
    fastqc --outdir $OUTDIR $TMPDIR/${sample}.filtered.R1.fastq.gz
    echo
    fastqc --outdir $OUTDIR $TMPDIR/${sample}.filtered.R2.fastq.gz
    
    
    ##BCread:
    #RNA, iPCR - R1
    #DNA, RNA 10x - R2
    
    ##UMI (additional Ns may be present from PCR primers but do not represent a true UMI)
    #RNA: UMI is preserved on both reads
    #RNA 10x: UMI is in the on-bead RT primer (plasmid read); Ns on RNA (R1) are not true UMI
    #DNA: UMI is only preserved in the plasmid read 
    #iPCR: No UMI
    #ChIP: UMI is only preserved in the plasmid read (4N for all)
    
    
    echo
    echo "Trimming/extracting UMI from R1 files ${f1} and R2 files ${f2}"
    date
    
    if [[ "${sampleType}" == "10xRNA" ]]; then
        echo "Extracting cell BC and UMI for 10x data"
        bc1pattern="^(?P<cell_1>.{16})(?P<umi_1>.{12})"
    else
        echo "Trimming ${R1trim} bp from R1"
        bc1pattern="^(?P<umi_1>.{$R1trim})"
    fi
    
    #Only extract UMI from R2 for RNA or iPCR samples (bcread is R1), throw it away for DNA samples
    if [[ "${R2trim}" -gt 0 && "${bcread}" == "R1" ]]; then
        echo "Trimming ${R2trim} bp from R2"
        bc2pattern="^(?P<umi_1>.{$R2trim})"
    else
        bc2pattern="^(?P<discard_1>.{$R2trim})"
        echo "Removing ${R2trim} bp from R2"
    fi
    
    
    echo "umi_tools extract --bc-pattern \"${bc1pattern}\" --bc-pattern2 \"${bc2pattern}\""
    #TODO this is pretty slow for big 10xRNA samples
    #TODO --quality-filter-threshold for UMI filtering?
    #umi_tools extract --help
    umi_tools extract --quality-encoding=phred33 --quality-filter-threshold=30 --compresslevel=9 -v 0 --log2stderr --extract-method regex --bc-pattern "${bc1pattern}" --bc-pattern2 "${bc2pattern}" -I $TMPDIR/${sample}.filtered.R1.fastq.gz -S $OUTDIR/${sample}.${R1file}.fastq.gz --read2-in=$TMPDIR/${sample}.filtered.R2.fastq.gz --read2-out=$OUTDIR/${sample}.${R2file}.fastq.gz
    
    
    echo
    echo "Weblogo of raw UMI"
    date
    #Use tail to run through to end of file so zcat doesn't throw an error code
    #UMI gets put on both fastq, so just look at BC file
    UMIlength=`zcat -f $OUTDIR/${sample}.BC.fastq.gz | awk 'BEGIN {OFS="\t"} NR % 4 == 1 {split($1, readname, "_"); print length(readname[2])}' | tail -1`
    if [[ "${UMIlength}" -gt 0 ]]; then
        echo "Making weblogo of UMI"
        zcat -f $OUTDIR/${sample}.BC.fastq.gz | awk 'BEGIN {OFS="\t"} NR % 4 == 1 {split($1, readname, "_"); print readname[2]}' | shuf -n 1000000 | awk -F "\t" 'BEGIN {OFS="\t"} {print ">id-" NR; print}' | weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} UMI sequence" --stacks-per-line 100 > $TMPDIR/${sample}.raw.UMI.eps
        convert $TMPDIR/${sample}.raw.UMI.eps $OUTDIR/${sample}.raw.UMI.png
    fi
    
    echo
    echo "Weblogo of raw cellBC"
    date
    #Use tail to run through to end of file so zcat doesn't throw an error code
    #cellBC gets put on both fastq, so just look at BC file
    cellBClength=`zcat -f $OUTDIR/${sample}.BC.fastq.gz | awk 'BEGIN {OFS="\t"} NR % 4 == 1 {split($1, readname, "_"); print length(readname[3])}' | tail -1`
    if [[ "${cellBClength}" -gt 0 ]]; then
        echo "Making weblogo of cellBC"
        zcat -f $OUTDIR/${sample}.BC.fastq.gz | awk 'BEGIN {OFS="\t"} NR % 4 == 1 {split($1, readname, "_"); print readname[3]}' | shuf -n 1000000 | awk -F "\t" 'BEGIN {OFS="\t"} {print ">id-" NR; print}' | weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} R1 cellBC sequence" --stacks-per-line 100 > $TMPDIR/${sample}.raw.cellBC.eps
        convert $TMPDIR/${sample}.raw.cellBC.eps $OUTDIR/${sample}.raw.cellBC.png
    fi
    
    
    ###Weblogo of processed reads
    echo "Weblogo of processed reads"
    date
    zcat -f $OUTDIR/${sample}.BC.fastq.gz | awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2' | shuf -n 1000000 | awk '{print ">id-" NR; print}' |
    weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} BC processed sequence" --stacks-per-line 100 > $TMPDIR/${sample}.BC.processed.eps
    convert $TMPDIR/${sample}.BC.processed.eps $OUTDIR/${sample}.BC.processed.png
    
    #BUGBUG needs to be manually disabled if nothing left on plasmid read (e.g. 28bp sequencing of 10xRNA R1)
    if [ 1 -eq 1 ]; then
        zcat -f $OUTDIR/${sample}.plasmid.fastq.gz | awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2' | shuf -n 1000000 | awk '{print ">id-" NR; print}' |
        weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} plasmid processed sequence" --stacks-per-line 100 > $TMPDIR/${sample}.plasmid.processed.eps
        convert $TMPDIR/${sample}.plasmid.processed.eps $OUTDIR/${sample}.plasmid.processed.png
    fi
fi


###Finally submit jobs
numlines=`zcat -f $OUTDIR/${sample}.BC.fastq.gz | wc -l`
chunksize=2000000 #Split fastq into 500,000 reads for deduplication (500,000 x 4)
numjobs=`echo "${numlines} / ${chunksize}" | bc -l -q`
numjobs=$(floor ${numjobs})
numjobs=`echo "${numjobs} + 1" | bc -l -q`
#BUGBUG verify must be int. what to do with last chunk?
echo
echo "${numlines} lines to process in chunks of ${chunksize}"
date

if [ ${runExtract} -eq 1 ]; then
    if [[ ${sampleType} == "iPCR" ]]; then
        qsub -S /bin/bash -t 1-${numjobs} -terse -j y -N mapintegrations.${sample} -o ${sample} -b y "${src}/mapIntegrations.sh ${sample} ${BCreadSeq} ${bclen} ${chunksize} ${plasmidSeq} ${extractBCargs}" | perl -pe 's/[^\d].+$//g;' > sgeid.extract.${sample}
    else
        qsub -S /bin/bash -t 1-${numjobs} -terse -j y -N extract.${sample} -o ${sample} -b y "${src}/extractBCcounts.sh ${sample} ${BCreadSeq} ${bclen} ${chunksize} ${plasmidSeq} ${extractBCargs}" | perl -pe 's/[^\d].+$//g;' > sgeid.extract.${sample}
    fi
fi


if [ ${runExtract} -eq 1 ]; then
    mergeHold="-hold_jid `cat sgeid.extract.${sample}`"
else
    mergeHold=""
fi

if [ ${runMerge} -eq 1 ]; then
    echo
    echo "Submitting merge"
    
    echo "Will merge ${numjobs} files"
    bcfiles=`seq 1 ${numjobs} | xargs -L 1 -I {} echo -n "${sample}/${sample}.{}.barcodes.txt.gz "`
    echo -e "Will merge barcode files: ${bcfiles}\n"
    if [[ ${sampleType} == "iPCR" ]]; then
        bamfiles=`seq 1 ${numjobs} | xargs -L 1 -I {} echo -n "${sample}/${sample}.{}.bam "`
        echo -e "Will merge bamfiles files: ${bamfiles}\n"
    else
        bamfiles=""
    fi
    
    cat <<EOF | qsub -S /bin/bash -j y -b y -N merge.${sample} -o ${OUTDIR} -terse ${mergeHold} | perl -pe 's/[^\d].+$//g;' > sgeid.merge.${sample}
    set -eu -o pipefail
    echo "Merging barcodes"
    zcat -f ${bcfiles} | pigz -p ${NSLOTS} -c -9 > $OUTDIR/${sample}.barcodes.preFilter.txt.gz
    rm -f ${bcfiles}
    
    if [[ ${sampleType} == "iPCR" ]]; then
        echo "Merging bam files"
        if [[ `echo ${bamfiles} | wc | awk '{print $2}'` -gt 1 ]]; then 
            samtools merge -f -l 9 $OUTDIR/${sample}.bam ${bamfiles}
        else 
            cp ${bamfiles} $OUTDIR/${sample}.bam
        fi
        samtools index $OUTDIR/${sample}.bam
        rm -f ${bamfiles}
    fi
    
    echo "Done"
EOF
    rm -f sgeid.extract.${sample}
fi


if [ ${runMerge} -eq 1 ]; then
    analysisHold="-hold_jid `cat sgeid.merge.${sample}`"
else
    analysisHold=""
fi

echo
echo "Submitting analysis"
if [[ ${sampleType} == "iPCR" ]]; then
    minReadCutoff=2
else
    minReadCutoff=10
fi

#-o ${OUTDIR} breaks Jesper's Flowcell_Info.sh
cat <<EOF | qsub -S /bin/bash -j y -b y -N ${sample} -terse ${analysisHold} # | perl -pe 's/[^\d].+$//g;' > sgeid.analysis
set -eu -o pipefail
${src}/analyzeBCcounts.sh ${minReadCutoff} ${sample}

if [[ ${sampleType} == "iPCR" ]]; then
    ${src}/analyzeIntegrations.sh ${sample}
fi
EOF
rm -f sgeid.merge.${sample}


echo "Done!!!"
date
