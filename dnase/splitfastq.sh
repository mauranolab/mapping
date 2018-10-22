#!/bin/bash
set -eu -o pipefail

#Based on stampipes/scripts/fastq/splitfastq.bash

SAMPLE_NAME=$1
R1_FASTQ=$2
R2_FASTQ=$3
FASTQ_BAK=$4
CHUNK_SIZE=$5

CHUNK_SIZE="${CHUNK_SIZE:-16000000}"
LINE_CHUNK=$((${CHUNK_SIZE} * 4))


mkdir -p ${TMPDIR}/splitfastq.${SAMPLE_NAME}
echo "Splitting FASTQ files into ${CHUNK_SIZE} reads using ${TMPDIR}/splitfastq.${SAMPLE_NAME}"
date

echo "Splitting R1 ${R1_FASTQ}"
zcat ${R1_FASTQ} | split -l ${LINE_CHUNK} -d -a 3 - "${TMPDIR}/splitfastq.${SAMPLE_NAME}/${SAMPLE_NAME}_R1_"
date

#BUGBUG doesn't work with SE fastq file
if [ true ]; then
    echo "Splitting R2 ${R2_FASTQ}"
    zcat ${R2_FASTQ} | split -l ${LINE_CHUNK} -d -a 3 - "${TMPDIR}/splitfastq.${SAMPLE_NAME}/${SAMPLE_NAME}_R2_"
    date
fi

SPLIT_COUNT=`ls ${TMPDIR}/splitfastq.${SAMPLE_NAME}/${SAMPLE_NAME}_R1_??? | wc -l`

if [ "$SPLIT_COUNT" -eq 1 ]; then
    # We only have one file split
    echo "Files do not exceed ${CHUNK_SIZE}; not touching original files"
    exit 1
else
    mkdir -p ${FASTQ_BAK}

    mv ${R1_FASTQ} ${FASTQ_BAK}
    mv ${R2_FASTQ} ${FASTQ_BAK}
    orig_dir=`dirname ${R1_FASTQ}`
    orig_dir2=`dirname ${R2_FASTQ}`
    if [[ "${orig_dir}" != "${orig_dir2}" ]]; then
        echo "ERROR: ${R1_FASTQ} in different location than ${R2_FASTQ}"
        exit 2
    fi

    for RAW_FILE in `find ${TMPDIR}/splitfastq.${SAMPLE_NAME} -name "${SAMPLE_NAME}_R?_???"`; do
        FASTQ_FILENAME=`basename ${RAW_FILE}`
        if [ -e "${orig_dir}/${FASTQ_FILENAME}.fastq.gz" ]; then
            echo "ERROR: ${orig_dir}/${FASTQ_FILENAME}.fastq.gz already exists!"
            exit 3
        fi
        echo "Compressing ${RAW_FILE}"
        date
        mv ${RAW_FILE} ${RAW_FILE}.fastq
        bgzip -@ $NSLOTS -c ${RAW_FILE}.fastq > ${orig_dir}/${FASTQ_FILENAME}.fastq.gz
    done
fi

echo "Done compressing"
date
