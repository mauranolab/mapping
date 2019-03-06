#!/bin/bash
set -euo pipefail

#Based on /home/maagj01/scratch/transposon/submitMerged.sh

src=/vol/mauranolab/transposon/src

sample=$1
shift
indivsamples=$@

bcfiles=`echo ${indivsamples} | perl -pe 's/\/+$//g;' -e 's/\/+( +)/\1/g;' | awk 'BEGIN {ORS=" "} {for(i=1; i<=NF; i++) {split($i, path, "/"); print $i "/" path[length(path)] ".barcodes.preFilter.txt" }}'`
echo -e "Will merge barcode files: ${bcfiles}\n"

if [[ ${indivsamples} =~ _iPCR ]]; then
    mergeIntegrations=1
    bamfiles=`echo ${indivsamples} | perl -pe 's/\/+$//g;' -e 's/\/+( +)/\1/g;' | awk 'BEGIN {ORS=" "} {for(i=1; i<=NF; i++) {split($i, path, "/"); print $i "/" path[length(path)] ".bam" }}'`
    echo -e "Will merge bam files: ${bamfiles}\n"
else
    mergeIntegrations=0
fi

qsub -S /bin/bash -j y --qos=full -N ${sample} -b y <<EOF
set -eu -o pipefail
mkdir -p ${sample}
echo "Merging barcodes"
cat ${bcfiles} > ${sample}/${sample}.barcodes.preFilter.txt
${src}/analyzeBCcounts.sh ${sample}

if [ ${mergeIntegrations} -eq 1 ]; then
    echo
    echo "Merging integration sites"
    samtools merge -f -l 9 $sample/${sample}.bam ${bamfiles}
    ${src}/analyzeIntegrations.sh ${sample}
fi
EOF
