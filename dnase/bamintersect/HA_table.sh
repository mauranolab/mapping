#!/bin/bash
######################################################################
set -eu -o pipefail
######################################################################
# This bam file needs to have an associated bai file:
bamfile=$1

# 3076=read is unmapped, or read is PCR, or read is sup alignment.
exclude_flags=$2

sampleOutdir=$3
sample_name=$4
HA5p=$5
HA3p=$6
src=$7

######################################################################
mkdir -p ${sampleOutdir}
TMPDIR=`mktemp -d`   # TMPDIR has no trailing slash
######################################################################
# Convert HA coordinates from bed to positional:

IFS=$':' read chrom coords <<< "${HA5p}"
IFS=$'-' read coord1 coord2 <<< "${coords}"
echo -e "${chrom}\t${coord1}\t${coord2}" > "${TMPDIR}/HA_BEDfile.bed"

IFS=$':' read chrom coords <<< "${HA3p}"
IFS=$'-' read coord1 coord2 <<< "${coords}"
echo -e "${chrom}\t${coord1}\t${coord2}" >> "${TMPDIR}/HA_BEDfile.bed"

######################################################################
# Get the header.
samtools view -H ${bamfile} > "${TMPDIR}/header.sam"

# Gets reads from both HAs, where both mates are mapped:
let "exclude_flags_plus_unmapped = ${exclude_flags} | 8" # 8 means mate is unmapped.
samtools view -F ${exclude_flags_plus_unmapped} -L "${TMPDIR}/HA_BEDfile.bed" ${bamfile} > "${TMPDIR}/all_HA_reads.sam"

# Get the read IDs with only one read in a HA:
cut -f1 "${TMPDIR}/all_HA_reads.sam" | sort | uniq -c | sed 's/^ *//g' | grep "^1 " | sed 's/^1 *//g' > "${TMPDIR}/Reads.txt"

# Reads.txt is all the reads with mapped mates, and only one read in a HA.

# Get the records for the read IDs with only one end in a HA.
# Note that we will get two lines for each such read ID: one in the HA, and one outside the HA.
/vol/mauranolab/cadlej01/git_repositories/mapping/dnase/bamintersect/subsetBAM.py  \
    --flags ${exclude_flags_plus_unmapped} \
    --readNames "${TMPDIR}/Reads.txt" \
    --bamFile_in ${bamfile} \
    --bamFile_out "${TMPDIR}/subsetBAM_output.bam"

# Need to do this so we can filter by region later.
samtools sort "${TMPDIR}/subsetBAM_output.bam" > "${TMPDIR}/subsetBAM_output_sorted.bam"
samtools index "${TMPDIR}/subsetBAM_output_sorted.bam"

# Get all the eligible reads (both inside and outside the HAs).
samtools view "${TMPDIR}/subsetBAM_output_sorted.bam" > "${TMPDIR}/all.sam"

# Get only the reads inside the HAs.
samtools view -L "${TMPDIR}/HA_BEDfile.bed" "${TMPDIR}/subsetBAM_output_sorted.bam" | sort > "${TMPDIR}/HA_noHdr.sam"

# Build sam files with headers for the reads inside the HAs.
cp "${TMPDIR}/header.sam" "${TMPDIR}/HA.sam"
cat "${TMPDIR}/HA_noHdr.sam" >> "${TMPDIR}/HA.sam"

# Get the reads which are not inside the HAs, and build sam files with headers.
cp "${TMPDIR}/header.sam" "${TMPDIR}/mates.sam"
grep -v -F -f "${TMPDIR}/HA_noHdr.sam" "${TMPDIR}/all.sam" | sort >> "${TMPDIR}/mates.sam"

# Now call bamintersect.py to build the bed12 file.
samtools view -b "${TMPDIR}/HA.sam" > "${TMPDIR}/bamfile_HA.bam"
samtools view -b "${TMPDIR}/mates.sam" > "${TMPDIR}/bamfile_mates.bam"
${src}/bamintersect.py --src ${src} --bam1 "${TMPDIR}/bamfile_HA.bam" --bam2 "${TMPDIR}/bamfile_mates.bam" --outdir ${sampleOutdir} --make_csv --ReqFullyAligned
mv "${sampleOutdir}/dsgrep_out.csv" "${sampleOutdir}/${sample_name}_HA.bed"

