#!/bin/bash
set -euo pipefail

#TODO merge with submit.sh a la dnase pipeline -- main issue is whether to remove files or not


#Limit thread usage by python processes using OPENBLAS (esp. scipy). Set here and will be inherited by spawned jobs
#https://stackoverflow.com/questions/51256738/multiple-instances-of-python-running-simultaneously-limited-to-35
export OPENBLAS_NUM_THREADS=1

module load python/3.8.1
module load weblogo/3.5.0
module load ImageMagick
module load bedtools
module load ucsckentutils/379
module load miller
module load pigz
module load samtools/1.14

src=$( dirname "${BASH_SOURCE[0]}" )

#Hardcoded right now rather than as parameters like dnase pipeline
runMerge=1

###Parse command line args
sample=$1
analyzeBCargs="--freq 0.01"
shift 1

while getopts ':b:' opt; do
    case $opt in
    b) analyzeBCargs="$OPTARG";;
    ?) echo "ERROR submit: unrecognize flag '-$OPTARG'"
       exit 2;;
    *) echo "ERROR submit: impossible!"
       exit 2;;
    esac
done
shift $((OPTIND-1))

indivsamples=$@


echo "Processing sample ${sample}"
echo "Running on $HOSTNAME. Using $TMPDIR as tmp"
date


OUTDIR=${sample}
mkdir -p ${OUTDIR}


#export NSLOTS=1


if [[ ${indivsamples} =~ _iPCR ]]; then
    sampleType="iPCR"
    minReadCutoff=2
    bamfiles=`echo ${indivsamples} | perl -pe 's/\/+$//g;' -e 's/\/+( +)/\1/g;' | awk 'BEGIN {ORS=" "} {for(i=1; i<=NF; i++) {split($i, path, "/"); print $i "/" path[length(path)] ".bam" }}'`
    echo -e "Will merge bam files: ${bamfiles}\n"
else
    #As long as it's not iPCR
    sampleType="DNA"
    minReadCutoff=10
    bamfiles=""
fi


if [ ${runMerge} -eq 1 ]; then
    bcfiles=`echo ${indivsamples} | perl -pe 's/\/+$//g;' -e 's/\/+( +)/\1/g;' | awk 'BEGIN {ORS=" "} {for(i=1; i<=NF; i++) {split($i, path, "/"); print $i "/" path[length(path)] ".barcodes.preFilter.txt.gz" }}'`
    echo -e "Will merge barcode files: ${bcfiles}\n"
    
    cat <<EOF | qsub -S /bin/bash -j y -b y -N merge.${sample} -o ${sample} -terse > sgeid.merge.${sample}
    set -eu -o pipefail
    echo "Merging barcodes from files: ${bcfiles}"
    date
    
    zcat -f ${bcfiles} | pigz -p ${NSLOTS} -c -9 > ${OUTDIR}/${sample}.barcodes.preFilter.txt.gz
    if [[ ${sampleType} == "iPCR" ]]; then
        echo "Merging bam files"
        date
        if [[ `echo ${bamfiles} | wc | awk '{print $2}'` -gt 1 ]]; then 
            samtools merge -f -l 9 $OUTDIR/${sample}.bam ${bamfiles}
        else 
            cp ${bamfiles} $OUTDIR/${sample}.bam
        fi
        samtools index $OUTDIR/${sample}.bam
    fi
    
    echo
    echo Done
    date
EOF
fi

if [ ${runMerge} -eq 1 ]; then
    analysisHold="-hold_jid `cat sgeid.merge.${sample}`"
else
    analysisHold=""
fi


#-o ${sample} breaks Jesper's Flowcell_Info.sh
cat <<EOF | qsub -S /bin/bash -j y -b y -N ${sample} -terse ${analysisHold} | perl -pe 's/[^\d].+$//g;'  #> sgeid.analysis
set -eu -o pipefail

${src}/analyzeBCcounts.sh ${minReadCutoff} ${sample} ${analyzeBCargs}
if [[ ${sampleType} == "iPCR" ]]; then
    ${src}/analyzeIntegrations.sh ${sample}
fi
EOF
rm -f sgeid.merge.${sample}

