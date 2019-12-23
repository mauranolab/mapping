#!/bin/bash
######################################################################
set -eu -o pipefail
######################################################################
# This bam file needs to have an associated bai file:
bamfile=$1
sampleOutdir=$2
sample_name=$3
INTERMEDIATEDIR=$4
src=$5
ReqFullyAligned=$6

# Exclude read with flags: unmapped, failed QC flag, read is duplicate, or read is sup alignment.
exclude_flags=3596

######################################################################
echo "Generating HA table"

# Get the header.
samtools view -H ${bamfile} > "${TMPDIR}/header.sam"

# Get reads from both HAs, where both mates are mapped.  Then get the read IDs with only one read in a HA (via cut & awk).
samtools view -F ${exclude_flags} -L ${INTERMEDIATEDIR}/HA_coords.bed ${bamfile} | cut -f1 | sort | uniq -c | awk '$1==1 {print $2}' > "${TMPDIR}/Reads.txt"
# Reads.txt is all the reads with mapped mates, and only one read in a HA.

# Get the records for the read IDs with only one end in a HA.
# Note that we will get two lines for each such read ID: one in the HA, and one outside the HA.
${src}/subsetBAM.py  \
    --exclude_flags ${exclude_flags} \
    --readNames "${TMPDIR}/Reads.txt" \
    --bamFile_in ${bamfile} \
    --bamFile_out "${TMPDIR}/subsetBAM_output.bam"

# Need to sort & index here so we can filter by region below.
samtools sort "${TMPDIR}/subsetBAM_output.bam" > "${TMPDIR}/subsetBAM_output_sorted.bam"
samtools index "${TMPDIR}/subsetBAM_output_sorted.bam"

# Get only the reads inside the HAs.
samtools view -L ${INTERMEDIATEDIR}/HA_coords.bed "${TMPDIR}/subsetBAM_output_sorted.bam" | tee "${TMPDIR}/HAreads_noHdr.sam" | cat ${TMPDIR}/header.sam - | samtools view -h -b | samtools sort -n > "${TMPDIR}/bamfile_HAreads.bam"

# Get the reads which are not inside the HAs.
samtools view "${TMPDIR}/subsetBAM_output_sorted.bam" | grep -v -F -f "${TMPDIR}/HAreads_noHdr.sam" | cat ${TMPDIR}/header.sam - | samtools view -h -b | samtools sort -n > ${TMPDIR}/bamfile_mates.bam

# Now call bamintersect.py to build the bed12 file.
${src}/bamintersect.py --src ${src} --bam1 "${TMPDIR}/bamfile_HAreads.bam" --bam2 "${TMPDIR}/bamfile_mates.bam" --outdir ${sampleOutdir} --make_bed ${ReqFullyAligned}
mv "${sampleOutdir}/dsgrep_out.csv" "${sampleOutdir}/${sample_name}.HA.bed"

