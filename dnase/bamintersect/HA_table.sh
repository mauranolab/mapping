#!/bin/bash
######################################################################
set -eu -o pipefail
######################################################################
# This script is run on a single bam file in order to find read pairs
# in which one read is in "Zone 1" and the other is in "Zone 2".
# 
# The parameter "runHA" is used to indicate how the zones are defined:
#     runHA=HA indicates "Way 1" Zone 1 is within a Homology Arm, and Zone 2 is outside the Homology Arms.
#     runHA=AB indicates "Way 2" Zone 1 is within a cegsvector assembly, and zone 2 is within the associated cegsvector assembly backbone
######################################################################
runHA=$1             # can equal HA ("Homology Arms") or AB ("Assembly Backbone")
bamfile=$2           # This bam file needs to have an associated bai file.
sampleOutdir=$3
sample_name=$4
INTERMEDIATEDIR=$5
src=$6
exclude_flags=$7
ReqFullyAligned=$8

######################################################################

if [ ${runHA} = "AB" ]; then
    echo "[HA_table] Generating assemblyBackbone table"
    
    # Assume there are only 2 chromosomes: the assembly and the assembly backbone, with the assembly first.
    read -r assembly assembly_readlength all_other <<< $(samtools idxstats ${bamfile})
    
    # Make sure the naming conventions are being followed, since some manual effort goes into making these assembly files.
    # We require only two chromosomes: an [assembly] and an [assembly]_backbone
    # If there is a misspelling, extra chromosomes, or if they are out of order, an error happens.
    backbone=$(samtools idxstats ${bamfile} | tail -n +2 | head -n 1 | cut -f1)  # This should be the backbone chromosome.
    numChroms=$(samtools idxstats ${bamfile} | awk -F "\t" 'BEGIN {OFS="\t"} $1!="*"' | wc -l)
    if [[ "${numChroms}" != 2 ]]; then
        echo "[HA_table] ${assembly} has wrong number of chromosomes, quitting successfully."
        touch "${sampleOutdir}/${sample_name}.assemblyBackbone.bed"
        exit 0
    elif [ "${backbone}" != "${assembly}_backbone" ]; then
        # The second chromosome is not named [assembly]_backbone
        echo "[HA_table] WARNING: Assembly has two chromosomes but the second does not match the expected backbone name, quitting successfully: ${assembly} ${backbone}"
        echo ""
        echo "[HA_table] samtools idxstats ${bamfile} :"
        samtools idxstats ${bamfile}
        touch "${sampleOutdir}/${sample_name}.assemblyBackbone.bed"
        exit 0
    else
        # The naming conventions are being followed.
        echo "[HA_table] assembly is: ${assembly}"
        echo "[HA_table] backbone is: ${backbone}"
        
        # Define Zone 1. We want it to be the entire assembly.
        echo -e "${assembly}\t0\t${assembly_readlength}" > ${TMPDIR}/zone1.bed
        outputBed="${sampleOutdir}/${sample_name}.assemblyBackbone.bed"
    fi
elif [ ${runHA} = "HA" ]; then
    echo "[HA_table] Generating HA table"
    # Define Zone 1. We want it to be the HAs (both of them).
    cp ${INTERMEDIATEDIR}/HA_coords.bed ${TMPDIR}/zone1.bed
    outputBed="${sampleOutdir}/${sample_name}.HA.bed"
else
    echo "[HA_table] Bad runHA parameter in consolidated_HA_AB_table.sh"
    exit 1
fi

# Preserve the header.
samtools view -H ${bamfile} > "${TMPDIR}/header.sam"

# Get reads from Zone 1, where the mates are mapped to anywhere.  Then get the read IDs with only one read in Zone 1 (via cut & awk).
samtools view -F ${exclude_flags} -L ${TMPDIR}/zone1.bed ${bamfile} | cut -f1 | sort | uniq -c | awk '$1==1 {print $2}' > "${TMPDIR}/Reads.txt"
# Reads.txt is all the reads with mapped mates, and only one read in Zone 1.
# The other reads must be somewhere outside the HAs (Way 1) or inside the backbone (Way 2).

if [ ${runHA} = "AB" ]; then
    # This code block is run only for runHA="AB" ("AssemblyBackbone").
    num_lines=$(wc -l < "${TMPDIR}/Reads.txt")
    echo "[HA_table] num_lines is: ${num_lines}"
    if [ "${num_lines}" = "0" ]; then
        touch "${sampleOutdir}/${sample_name}.assemblyBackbone.bed"
        echo "[HA_table] Leaving due to no eligible reads."
        exit 0
    fi
fi

# Get the bam record lines for the read IDs with only one read in Zone 1.
# Note that we will get two lines for each such read ID: one in Zone 1, and one in Zone 2.
${src}/subsetBAM.py  \
    --exclude_flags ${exclude_flags} \
    --readNames "${TMPDIR}/Reads.txt" \
    --bamFile_in ${bamfile} \
    --bamFile_out "${TMPDIR}/subsetBAM_output.bam"

# Need to sort & index here so we can filter by region below.
samtools sort "${TMPDIR}/subsetBAM_output.bam" > "${TMPDIR}/subsetBAM_output_sorted.bam"
samtools index "${TMPDIR}/subsetBAM_output_sorted.bam"

# Get only the reads inside Zone 1.
samtools view  -L ${TMPDIR}/zone1.bed "${TMPDIR}/subsetBAM_output_sorted.bam" | tee "${TMPDIR}/Zone1reads_noHdr.sam" | cat ${TMPDIR}/header.sam - | samtools view -h -b | samtools sort -n > "${TMPDIR}/bamfile_Zone1reads.bam"

# Get the reads which are in Zone 2.
samtools view "${TMPDIR}/subsetBAM_output_sorted.bam" | grep -v -F -f "${TMPDIR}/Zone1reads_noHdr.sam" | cat ${TMPDIR}/header.sam - | samtools view -h -b | samtools sort -n > ${TMPDIR}/bamfile_Zone2reads.bam

# Now call bamintersect.py to build the bed12 file.
${src}/bamintersect.py --src ${src} --bam1 "${TMPDIR}/bamfile_Zone1reads.bam" --bam2 "${TMPDIR}/bamfile_Zone2reads.bam" --outdir ${sampleOutdir} --make_bed ${ReqFullyAligned}

sort-bed ${sampleOutdir}/dsgrep_out.csv > ${outputBed}
rm -f ${sampleOutdir}/dsgrep_out.csv

echo "[HA_table] Done]"
date
