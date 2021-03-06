#!/bin/bash
set -e -o pipefail

src=$( dirname "${BASH_SOURCE[0]}" )

OUTBASE=$1
shift
outbases=$@


echo
echo "Merging AllBCs"
#Just used for "Overlap between libraries by read depth" plots
bcfiles=""
for f in ${outbases}; do
    bcfiles="${bcfiles} ${f}.AllBCs.txt"
done
echo ${bcfiles}

mlr --tsv cat ${bcfiles} > ${OUTBASE}.AllBCs.txt


echo
echo "Merging mapped.txt"
mappedfiles=""
for f in ${outbases}; do
    mappedfiles="${mappedfiles} ${f}.mapped.txt"
    head -1 ${f}.mapped.txt > ${OUTBASE}.mapped.txt
done
echo ${mappedfiles}

tail -q -n +2 ${mappedfiles} | sort-bed - >> ${OUTBASE}.mapped.txt


${src}/analyzeCombinedSignalWithInsertions.sh ${OUTBASE}


echo
echo "Done!!!"
date

