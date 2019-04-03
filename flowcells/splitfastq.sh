#!/bin/bash
set -eu -o pipefail

#Based on stampipes/scripts/fastq/splitfastq.bash
#Requires sequences/qualities to be on single line
#NB doesn't work with SE fastq files


module add pigz


SAMPLE_NAME=$1
R1_FASTQ=$2
R2_FASTQ=$3
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

echo "Splitting R1 ${R1_FASTQ}"
pigz -dc ${R1_FASTQ} | split -l ${LINE_CHUNK} -d -a 3 - "${TMPDIR}/splitfastq.${SAMPLE_NAME}/${SAMPLE_NAME}_R1_"
date

#BUGBUG doesn't work with SE fastq file
if [ true ]; then
    echo "Splitting R2 ${R2_FASTQ}"
    pigz -dc ${R2_FASTQ} | split -l ${LINE_CHUNK} -d -a 3 - "${TMPDIR}/splitfastq.${SAMPLE_NAME}/${SAMPLE_NAME}_R2_"
    date
fi

SPLIT_COUNT=`ls ${TMPDIR}/splitfastq.${SAMPLE_NAME}/${SAMPLE_NAME}_R1_??? | wc -l`

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
    # -11 (zopfli) saves ~6% beyond -9 in one test at about 10x runtime
    pigz -p ${NSLOTS} -c -9 ${RAW_FILE}.fastq > ${orig_dir}/${FASTQ_FILENAME}.fastq.gz
    rm -f ${RAW_FILE}.fastq
done


if [ -z "${FASTQ_BAK}" ]; then
    echo "Removing original files"
    rm -f ${R1_FASTQ} ${R2_FASTQ}
else
    echo "Backing up files to ${FASTQ_BAK}"
    mkdir -p ${FASTQ_BAK}
    mv ${R1_FASTQ} ${FASTQ_BAK}
    mv ${R2_FASTQ} ${FASTQ_BAK}
fi


echo "Done compressing"
date

