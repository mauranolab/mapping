#!/bin/bash
set -euo pipefail

#Based on /home/maagj01/scratch/transposon/submitMerged.sh

src=/vol/mauranolab/transposon/src

sample=$1
shift
indivsamples=$@

bcfiles=`echo $indivsamples | perl -pe 's/\/+$//g;' -e 's/\/+( +)/\1/g;' | awk 'BEGIN {ORS=" "} {for(i=1; i<=NF; i++) {split($i, path, "/"); print $i "/" path[length(path)] ".barcodes.preFilter.txt" }}'`
echo -e "Will merge barcode files: ${bcfiles}\n"
qsub -S /bin/bash -j y --qos=full -N ${sample} -b y <<EOF
set -eu -o pipefail
mkdir -p ${sample}
echo "Merging barcodes"
cat ${bcfiles} > ${sample}/${sample}.barcodes.preFilter.txt
${src}/analyzeBCcounts.sh ${sample}
EOF
