#!/bin/bash
set -euo pipefail

#Based on /home/maagj01/scratch/transposon/submitMerged.sh

#Limit thread usage by python processes using OPENBLAS (esp. scipy). Set here and will be inherited by spawned jobs
#https://stackoverflow.com/questions/51256738/multiple-instances-of-python-running-simultaneously-limited-to-35
export OPENBLAS_NUM_THREADS=1

module load pigz
module load samtools/1.9

src=/vol/mauranolab/transposon/src

sample=$1
shift
indivsamples=$@

bcfiles=`echo ${indivsamples} | perl -pe 's/\/+$//g;' -e 's/\/+( +)/\1/g;' | awk 'BEGIN {ORS=" "} {for(i=1; i<=NF; i++) {split($i, path, "/"); print $i "/" path[length(path)] ".barcodes.preFilter.txt.gz" }}'`
echo -e "Will merge barcode files: ${bcfiles}\n"

if [[ ${indivsamples} =~ _iPCR ]]; then
    mergeIntegrations=1
    minReadCutoff=2
    bamfiles=`echo ${indivsamples} | perl -pe 's/\/+$//g;' -e 's/\/+( +)/\1/g;' | awk 'BEGIN {ORS=" "} {for(i=1; i<=NF; i++) {split($i, path, "/"); print $i "/" path[length(path)] ".bam" }}'`
    echo -e "Will merge bam files: ${bamfiles}\n"
else
    mergeIntegrations=0
    minReadCutoff=10
    bamfiles=""
fi

qsub -S /bin/bash -j y -N ${sample} -b y <<EOF
set -eu -o pipefail
mkdir -p ${sample}
echo "Merging barcodes"
zcat ${bcfiles} | pigz -p ${NSLOTS} -9 > ${sample}/${sample}.barcodes.preFilter.txt.gz
${src}/analyzeBCcounts.sh ${minReadCutoff} ${sample}

if [ ${mergeIntegrations} -eq 1 ]; then
    echo
    echo "Merging integration sites"
    samtools merge -f -l 9 $sample/${sample}.bam ${bamfiles}
    ${src}/analyzeIntegrations.sh ${sample}
fi
EOF
