#!/bin/bash
######################################################################
set -eu -o pipefail
######################################################################
# This script is run on a single bam file in order to find read pairs
# in which one read is in "Zone 1" and the other is in "Zone 2". The
# zones can be defined in two ways:
# 
# Way 1:  Zone 1 is within a Homology Arm, and Zone 2 is outside the Homology Arms.
# Way 2:  Zone 1 is within a cegsvector assembly, and zone 2 is within the associated
#         cegsvector assembly backbone.
#
# The parameter "runHA" is used to indicate which way is being utilized:
#     runHA=HA indicates "Way 1"
#     runHA=AB indicates "Way 2"
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
    echo "Generating assemblyBackbone table"

    # Assume there are only 2 chromosomes: the assembly and the assembly backbone, with the assembly first.
    read -r assembly assembly_readlength all_other <<< $(samtools idxstats ${bamfile})

    # Make sure the naming conventions are being followed, since some manual effort goes into making these assembly files.
    # We are expecting only two chromosomes: an [assembly] and an [assembly]_backbone
    # If there is a misspelling, extra chromosomes, or if they are out of order, an error happens.
    backbone=${assembly}_backbone
    backboneCheck=$(samtools idxstats ${bamfile} | tail -n +2 | head -n 1 | cut -f1)  # This should be the backbone chromosome.
    numChroms=$(wc -l <(samtools idxstats ${bamfile}) | cut -d ' ' -f1)   # Actually 1 bigger than the number of chroms. "*" is always the last line.

    if [ "${numChroms}" -gt 3 ]; then
        # There are three or more chromosomes.
        echo "${assembly} has too many chromosomes."
        touch "${sampleOutdir}/${sample_name}.assemblyBackbone.bed"
        exit 0
    elif [[ "${backboneCheck}" == "*" ]] && [[ "${numChroms}" == "2" ]]; then
        # There is no backbone, only a single [assembly]. As an example, this happens for pSpCas9.
        echo "${assembly} does not have a backbone chromosome."
        touch "${sampleOutdir}/${sample_name}.assemblyBackbone.bed"
        exit 0
    elif [ "${backbone}" != "${backboneCheck}" ]; then
        # We are not getting [assembly] and [assembly]_backbone
        echo "[assemblyBackbone_table.sh] ERROR: Problem with assignment of backbone names: ${assembly} ${backbone} ${backboneCheck}"
        echo ""
        echo "samtools idxstats ${bamfile} :"
        samtools idxstats ${bamfile}
        touch "${sampleOutdir}/${sample_name}.assemblyBackbone.bed"
        exit 0
    else
        # The naming conventions are being followed.
        echo assembly is: ${assembly}
        echo backbone is: ${backbone}

        # Define Zone 1. We want it to be the entire assembly.
        echo -e "${assembly}\t0\t${assembly_readlength}" > ${TMPDIR}/zone1.bed
        outputBed="${sampleOutdir}/${sample_name}.assemblyBackbone.bed"
    fi
elif [ ${runHA} = "HA" ]; then
    echo "Generating HA table"
    # Define Zone 1. We want it to be the HAs (both of them).
    cp ${INTERMEDIATEDIR}/HA_coords.bed ${TMPDIR}/zone1.bed
    outputBed="${sampleOutdir}/${sample_name}.HA.bed"
else
    echo "Bad runHA parameter in consolidated_HA_AB_table.sh"
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
    echo num_lines is: ${num_lines}
    if [ "${num_lines}" = "0" ]; then
        touch "${sampleOutdir}/${sample_name}.assemblyBackbone.bed"
        echo "Leaving consolidated_HA_AB_table.sh due to no eligible reads."
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

