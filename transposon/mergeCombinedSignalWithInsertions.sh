#!/bin/bash
set -e -o pipefail

src=/vol/mauranolab/mapped/src/transposon/

OUTBASE=$1
shift
outbases=$@


echo
echo "Merging AllBCs"
bcfiles=""
for f in ${outbases}; do
    bcfiles="${bcfiles} ${f}.AllBCs.txt"
done
echo ${bcfiles}

mlr --tsv cat ${bcfiles} > ${OUTBASE}.AllBCs.txt


echo
echo "Merging txt"
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

