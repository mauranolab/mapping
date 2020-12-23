#!/bin/bash
set -eu -o pipefail

######################################################################
# This script identifies reads in a bam file where one of the mates overlaps exactly once with regions defined in a bed file
#
# Outputs single bed12 file: ${outputBed}
######################################################################

bam=${1}           # This bam file needs to have an associated bai file.
regionsOfInterestFile=${2}
exclude_flags=${3}
outputBed=${4}
src=${5}


###################
# Preserve the header.
samtools view -H ${bam} > ${TMPDIR}/header.sam

# Get alignments where one read is in region of interest
samtools view -F ${exclude_flags} -L ${regionsOfInterestFile} ${bam} | cut -f1 |
#Get the alignments with only one read in region of interest
sort | uniq -c | awk '$1==1 {print $2}' > ${TMPDIR}/zone1_readnames.txt

num_lines=$(wc -l < ${TMPDIR}/zone1_readnames.txt)
echo "[bedOverlapPE] ${num_lines} eligible alignments"
if [ "${num_lines}" = "0" ]; then
    touch ${outputBed}
    echo "[bedOverlapPE] No eligible reads; quitting successfully"
    exit 0
fi

# Get all alignments where one of the reads is in region of interest; one alignment will be outside
${src}/subsetBAM.py --exclude_flags ${exclude_flags} --include_readnames ${TMPDIR}/zone1_readnames.txt ${bam} ${TMPDIR}/subsetBAM_output.bam

# Need to sort & index here so we can filter by region below
samtools sort ${TMPDIR}/subsetBAM_output.bam > ${TMPDIR}/subsetBAM_output_sorted.bam
samtools index ${TMPDIR}/subsetBAM_output_sorted.bam

# Get only the reads inside region of interest
samtools view  -L ${regionsOfInterestFile} ${TMPDIR}/subsetBAM_output_sorted.bam | tee ${TMPDIR}/zone1reads_noHdr.sam | cat ${TMPDIR}/header.sam - | samtools view -h -b | samtools sort -n > ${TMPDIR}/bamfile_zone1reads.bam

# Get the reads which are outside region of interest
samtools view ${TMPDIR}/subsetBAM_output_sorted.bam | grep -v -F -f ${TMPDIR}/zone1reads_noHdr.sam | cat ${TMPDIR}/header.sam - | samtools view -h -b | samtools sort -n > ${TMPDIR}/bamfile_Zone2reads.bam

# Now call bamintersect.py to build the bed12 file.
${src}/bamintersect.py --ReqFullyAligned ${TMPDIR}/bamfile_Zone2reads.bam ${TMPDIR}/bamfile_zone1reads.bam --src ${src} --bedout ${TMPDIR}/bamintersect_out.bed12

cat ${TMPDIR}/bamintersect_out.bed12 |
#In the case that we return reads on the same chromosome regularize output so that leftmost starting coord appears as bam1
awk -F "\t" 'BEGIN {OFS="\t"} {if($1==$7 && $8<$2) {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6} else {print}}' |
sort-bed - > ${outputBed}

echo "[bedOverlapPE] Done!"
date
