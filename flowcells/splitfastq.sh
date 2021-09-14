#!/bin/bash
set -eu -o pipefail

#Based on stampipes/scripts/fastq/splitfastq.bash
#Requires sequences/qualities to be on single line


module add pigz


SAMPLE_NAME=$1
FASTQ=$2
SUFFIX=$3
CHUNK_SIZE=$4

if [ "$#" -eq 5 ]; then
    FASTQ_BAK=$5
else
    FASTQ_BAK=""
fi

CHUNK_SIZE="${CHUNK_SIZE:-16000000}"
LINE_CHUNK=$((${CHUNK_SIZE} * 4))


mkdir -p ${TMPDIR}/splitfastq.${SAMPLE_NAME}
echo "Splitting FASTQ files into ${CHUNK_SIZE} reads using ${TMPDIR}/splitfastq.${SAMPLE_NAME}"
date

echo "Splitting ${SUFFIX} ${FASTQ}"
pigz -dc ${FASTQ} | split -l ${LINE_CHUNK} -d -a 3 - "${TMPDIR}/splitfastq.${SAMPLE_NAME}/${SAMPLE_NAME}_${SUFFIX}_"
date


SPLIT_COUNT=`ls ${TMPDIR}/splitfastq.${SAMPLE_NAME}/${SAMPLE_NAME}_${SUFFIX}_??? | wc -l`
echo "Split into ${SPLIT_COUNT} files"


orig_dir=`dirname ${FASTQ}`
for RAW_FILE in `find ${TMPDIR}/splitfastq.${SAMPLE_NAME} -name "${SAMPLE_NAME}_${SUFFIX}_???"`; do
    FASTQ_FILENAME=`basename ${RAW_FILE}`
    if [ -e "${orig_dir}/${FASTQ_FILENAME}.fastq.gz" ]; then
        echo "ERROR: ${orig_dir}/${FASTQ_FILENAME}.fastq.gz already exists!"
        exit 3
    fi
    echo "Compressing ${RAW_FILE}"
    date
    mv ${RAW_FILE} ${RAW_FILE}.fastq
    # -11 (zopfli) saves ~6% beyond -9 in one test at about 10x runtime
    pigz -p ${NSLOTS} -c -9 ${RAW_FILE}.fastq > ${orig_dir}/${FASTQ_FILENAME}.fastq.gz
    rm -f ${RAW_FILE}.fastq
done


if [ -z "${FASTQ_BAK}" ]; then
    echo "Removing original file"
    rm -f ${FASTQ}
else
    echo "Backing up file to ${FASTQ_BAK}"
    mkdir -p ${FASTQ_BAK}
    mv ${FASTQ} ${FASTQ_BAK}
fi


echo "Done compressing"
date

