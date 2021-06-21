#!/bin/bash
#
# delly variant caller.
#
#######################################################################################################################
# From Original - Needed for delly (could be put in user's .bashrc)
# Should no longer be needed, added to delly module load
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cm/local/apps/boost/1.68.0.a/lib64
#export LD_LIBRARY_PATH

module load delly/0.8.7
module load python/3.8.1 
module load bcftools/1.10.2
module load htslib/1.9       # bgzip
module load samtools/1.10    # tabix

bam_input=$1
genome=$2
BSID=$3
flowCell=$4

#######################################################################################################################

bam_input_BASE="${bam_input%.bam}"                               # remove .bam
bam_input_sampleName="${bam_input_BASE##*/}.${flowCell}"         # retain only <sampleName>.<genome>.<flowCell>    Added flowCell. Not the full path.

# Make a place for the temp bam files with qcfail turned off.
bam_output_qcOK="${TMPDIR}/QC_OK_bamfile.bam"
bam_output_sorted="${TMPDIR}/QC_OK_bamfile_sorted.bam"

# Turn off all qcfail flags.
python changeFlags.py ${bam_input} ${bam_output_qcOK}

# Assume $baminput was already sorted...
# samtools sort -o ${bam_output_sorted} ${bam_output_qcOK}
mv ${bam_output_qcOK} ${bam_output_sorted}

samtools index ${bam_output_sorted}
# rm ${bam_output_qcOK}

# Look for variants (DEL INS DUP INV TRA).
delly call -t ALL \
    -o "${TMPDIR}/${bam_input_sampleName}.bcf" \
    -g "/vol/isg/annotation/fasta/${genome}/${genome}.fa.gz" \
    ${bam_output_sorted}

rm "${bam_output_sorted}" "${bam_output_sorted}.bai"

#mv "${TMPDIR}/${bam_input_sampleName}.bcf" "/home/grivam01/scripts/cegs/delly/tmp/test.bcf"

# Create a vcf file, and set it up for viewing in vcfTabix format.
# Only accept variants that pass the FILTER test.
outdir=/vol/mauranolab/grivam01/projects/delly_batch/browser_files
bcftools filter -i 'FILTER="PASS" & PE>10' "${TMPDIR}/${bam_input_sampleName}.bcf" > "${outdir}/${bam_input_sampleName}.vcf"
# bcftools view "${TMPDIR}/${bam_input_sampleName}.bcf" > "${outdir}/${bam_input_sampleName}.vcf"    # This way does not apply the FILTER test.
bgzip -f "${outdir}/${bam_input_sampleName}.vcf"
tabix -f -p vcf "${outdir}/${bam_input_sampleName}.vcf.gz"

# Save the custom track code.
if [ ! "${genome}" = "mm10" ]; then
    # deal with rn6 and t2t later
    genome=hg38
fi

echo -e "track type=vcfTabix name=\"${BSID}_${flowCell}_vcf\" description=\"${bam_input_sampleName}_vcf\" bigDataUrl=https://mauranolab:chromatin@cascade.isg.med.nyu.edu/~grivam01/analysis_tracks/delly_batch_browser_files/${bam_input_sampleName}.vcf.gz visibility=pack maxWindowToDraw=50000000 db=${genome}" > "${outdir}/${bam_input_sampleName}.customTrackCode.txt"
# echo -e "browser position chr3:34728000-34741000" >> "${outdir}/${bam_input_sampleName}.customTrackCode.txt"

