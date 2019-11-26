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

# Get the header.
samtools view -H ${bamfile} > "${TMPDIR}/header.sam"

# Gets reads from both HAs, where both mates are mapped:
let "exclude_flags_plus_unmapped = ${exclude_flags} | 8" # 8 means mate is unmapped.
samtools view -F ${exclude_flags_plus_unmapped} ${bamfile} ${HA5p} > "${TMPDIR}/all_p5HA_reads.sam"
samtools view -F ${exclude_flags_plus_unmapped} ${bamfile} ${HA3p}  > "${TMPDIR}/all_p3HA_reads.sam"

# Get the read IDs with only one read in a HA:
cut -f1 "${TMPDIR}/all_p5HA_reads.sam" | sort | uniq -c | sed 's/^ *//g' | grep "^1 " | sed 's/^1 *//g' > "${TMPDIR}/p5Reads.txt"
cut -f1 "${TMPDIR}/all_p3HA_reads.sam" | sort | uniq -c | sed 's/^ *//g' | grep "^1 " | sed 's/^1 *//g' > "${TMPDIR}/p3Reads.txt"

# p?Reads.txt is all the reads with mapped mates, and only one read in a HA.

# Get the records for the read IDs with only one end in a HA.
# Note that we will get two lines for each such read ID: one in the HA, and one outside the HA.
/vol/mauranolab/cadlej01/git_repositories/mapping/dnase/bamintersect/subsetBAM.py  \
    --flags ${exclude_flags_plus_unmapped} \
    --readNames "${TMPDIR}/p5Reads.txt" \
    --bamFile_in ${bamfile} \
    --bamFile_out "${TMPDIR}/subsetBAM_output_p5.bam"

/vol/mauranolab/cadlej01/git_repositories/mapping/dnase/bamintersect/subsetBAM.py  \
    --flags ${exclude_flags_plus_unmapped} \
    --readNames "${TMPDIR}/p3Reads.txt" \
    --bamFile_in ${bamfile} \
    --bamFile_out "${TMPDIR}/subsetBAM_output_p3.bam"

# Need to do this so we can filter by region later.
samtools sort "${TMPDIR}/subsetBAM_output_p5.bam" > "${TMPDIR}/subsetBAM_output_p5sorted.bam"
samtools index "${TMPDIR}/subsetBAM_output_p5sorted.bam"
samtools sort "${TMPDIR}/subsetBAM_output_p3.bam" > "${TMPDIR}/subsetBAM_output_p3sorted.bam"
samtools index "${TMPDIR}/subsetBAM_output_p3sorted.bam"

# Get all the eligible reads (both inside and outside the HAs).
samtools view "${TMPDIR}/subsetBAM_output_p5sorted.bam" > "${TMPDIR}/allp5.sam"
samtools view "${TMPDIR}/subsetBAM_output_p3sorted.bam" > "${TMPDIR}/allp3.sam"

# Get the reads inside the HAs.
samtools view "${TMPDIR}/subsetBAM_output_p5sorted.bam" ${HA5p} | sort > "${TMPDIR}/p5HA_noHdr.sam"
samtools view "${TMPDIR}/subsetBAM_output_p3sorted.bam" ${HA3p} | sort > "${TMPDIR}/p3HA_noHdr.sam"

# Build sam files with headers for the reads inside the HAs.
cp "${TMPDIR}/header.sam" "${TMPDIR}/p5HA.sam"
cat "${TMPDIR}/p5HA_noHdr.sam" >> "${TMPDIR}/p5HA.sam"
cp "${TMPDIR}/header.sam" "${TMPDIR}/p3HA.sam"
cat "${TMPDIR}/p3HA_noHdr.sam" >> "${TMPDIR}/p3HA.sam"

# Get the reads which are not inside the HAs, and build sam files with headers.
cp "${TMPDIR}/header.sam" "${TMPDIR}/p5mates.sam"
grep -v -F -f "${TMPDIR}/p5HA_noHdr.sam" "${TMPDIR}/allp5.sam" | sort >> "${TMPDIR}/p5mates.sam"
cp "${TMPDIR}/header.sam" "${TMPDIR}/p3mates.sam"
grep -v -F -f "${TMPDIR}/p3HA_noHdr.sam" "${TMPDIR}/allp3.sam" | sort >> "${TMPDIR}/p3mates.sam"


# Now call python to build the bed12 file
samtools view -b "${TMPDIR}/p5HA.sam" > "${TMPDIR}/bamfile_HA.bam"
samtools view -b "${TMPDIR}/p5mates.sam" > "${TMPDIR}/bamfile_mates.bam"
${src}/bamintersect.py --src ${src} --bam1 "${TMPDIR}/bamfile_HA.bam" --bam2 "${TMPDIR}/bamfile_mates.bam" --outdir ${TMPDIR} --max_mismatches 0 --make_csv --ReqFullyAligned
mv "${TMPDIR}/dsgrep_out.csv" "${sampleOutdir}/${sample_name}_p5HA_output.bed12"

samtools view -b "${TMPDIR}/p3HA.sam" > "${TMPDIR}/bamfile_HA.bam"
samtools view -b "${TMPDIR}/p3mates.sam" > "${TMPDIR}/bamfile_mates.bam"
${src}/bamintersect.py --src ${src} --bam1 "${TMPDIR}/bamfile_HA.bam" --bam2 "${TMPDIR}/bamfile_mates.bam" --outdir ${TMPDIR} --max_mismatches 0 --make_csv --ReqFullyAligned
mv "${TMPDIR}/dsgrep_out.csv" "${sampleOutdir}/${sample_name}_p3HA_output.bed12"

