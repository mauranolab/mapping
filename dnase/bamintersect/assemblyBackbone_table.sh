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
echo "Generating assemblyBackbone table"

# Get the header.
samtools view -H ${bamfile} > "${TMPDIR}/header.sam"

# For now, assume the assembly is first.
# Need some rules to determine when to run this script.
assembly=$(samtools idxstats ${bamfile} | head -n 1 | cut -f1)
backbone=${assembly}_backbone
backboneCheck=$(samtools idxstats ${bamfile} | tail -n +2 | head -n 1 | cut -f1)
if [ "${backbone}" != "${backboneCheck}" ]; then
    echo "[assemblyBackbone_table.sh] ERROR: Problem with assignment of backbone names: ${assembly} ${backbone} ${backboneCheck}"
    echo ""
    echo "samtools idxstats ${bamfile} :"
    samtools idxstats ${bamfile}
    # exit 1

    # For now:
    touch "${sampleOutdir}/${sample_name}.assemblyBackbone.bed"
    exit 0
fi
echo assembly is:  ${assembly}
echo backbone is: ${backbone}

# Get reads from both chromosomes, where both mates are mapped.  Then get the read IDs with only one read in a chromosome (via cut & awk).
samtools view -F ${exclude_flags} ${bamfile} ${assembly} | cut -f1 | sort | uniq -c | awk '$1==1 {print $2}' > "${TMPDIR}/Reads.txt"
# Reads.txt is all the reads with mapped mates, and only one read in the assembly chromosome.
# The others must be in the backbone, since there are only two chromosomes allowed.

num_lines=$(wc -l < "${TMPDIR}/Reads.txt")
echo num_lines is: ${num_lines}
if [ "${num_lines}" = "0" ]; then
    touch "${sampleOutdir}/${sample_name}.assemblyBackbone.bed"
    echo "Leaving assemblyBackbone_table.sh due to no eligible reads."
    exit 0
fi

# Get the records for the read IDs with only one end in the assembly.
# Note that we will get two lines for each such read ID: one in the assembly, and one in the backbone.
${src}/subsetBAM.py  \
    --exclude_flags ${exclude_flags} \
    --readNames "${TMPDIR}/Reads.txt" \
    --bamFile_in ${bamfile} \
    --bamFile_out "${TMPDIR}/subsetBAM_output.bam"

# Need to sort & index here so we can filter by region below.
samtools sort "${TMPDIR}/subsetBAM_output.bam" > "${TMPDIR}/subsetBAM_output_sorted.bam"
samtools index "${TMPDIR}/subsetBAM_output_sorted.bam"

# Get only the reads inside the assembly.
samtools view "${TMPDIR}/subsetBAM_output_sorted.bam" ${assembly} | tee "${TMPDIR}/HAreads_noHdr.sam" | cat ${TMPDIR}/header.sam - | samtools view -h -b | samtools sort -n > "${TMPDIR}/bamfile_HAreads.bam"

# Get the reads which are not inside the HAs.
samtools view "${TMPDIR}/subsetBAM_output_sorted.bam" | grep -v -F -f "${TMPDIR}/HAreads_noHdr.sam" | cat ${TMPDIR}/header.sam - | samtools view -h -b | samtools sort -n > ${TMPDIR}/bamfile_mates.bam

# Now call bamintersect.py to build the bed12 file.
${src}/bamintersect.py --src ${src} --bam1 "${TMPDIR}/bamfile_HAreads.bam" --bam2 "${TMPDIR}/bamfile_mates.bam" --outdir ${sampleOutdir} --make_bed ${ReqFullyAligned}
mv "${sampleOutdir}/dsgrep_out.csv" "${sampleOutdir}/${sample_name}.assemblyBackbone.bed"

