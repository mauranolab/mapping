#!/bin/bash
set -eu -o pipefail

## Variables are passed in via sbatch export.
##########################################################################################################
## This function implements verbose output for debugging purposes.
debug_fa() {
    if [ ${verbose} = "True" ]; then
        echo "${1}"
    fi
}
##########################################################################################################
## Merge the bam_intersect output files:

echo "Starting merge_bamintersect.sh"

# First check for the "no reads to merge" problem. Also deal with the usual pipefail issue via "true".
dir_names=($(ls -d ${INTERMEDIATEDIR}/bamintersectPyOut/*/ 2> /dev/null || true))      ## These have a trailing /
numElements="${#dir_names[@]}"
if [ "${numElements}" = "0" ]; then
    # Don't generate an error. Create some appropriately empty output files.
    echo -e "chromosome-bam2\tStart_Pos\tEnd_Pos\tWidth\tNearest_Gene\tPost-filter_Reads\tchromosome-bam1\tSample" > "${sampleOutdir}/${sample_name}.counts.txt"
    echo -e "chromosome-bam2\tStart_Pos\tEnd_Pos\tWidth\tNearest_Gene\tPost-filter_Reads\tchromosome-bam1\tSample" > "${sampleOutdir}/${sample_name}.informative.counts.txt"
    echo "" > "${sampleOutdir}/${sample_name}.bed"
    echo "" > "${sampleOutdir}/${sample_name}.informative.bed"
    echo "There are no reads with unmapped mates to analyze." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    # Not generating bam files for this case.

    ## Cleanup
    rm -r ${INTERMEDIATEDIR}
    date
    echo "Done with merge_bamintersect.sh"
    exit 0
fi

first="True"
for i in ${dir_names[@]}; do
    if [ ${first} = "True" ]; then
        cp "${i}dsgrep_out.csv" "${INTERMEDIATEDIR}/unsorted_dsgrep.bed"
        first="False"
    else
        cat "${i}dsgrep_out.csv" >> "${INTERMEDIATEDIR}/unsorted_dsgrep.bed"
    fi
done

date
echo Sorting...
sort-bed "${INTERMEDIATEDIR}/unsorted_dsgrep.bed" > "${sampleOutdir}/${sample_name}.bed"
date

debug_fa "Finished sorting dsgrep.bed"

##########################################################################################################
## Create the final output tables.

if [ ${make_csv} != "--make_csv" ]; then
    # Need the bed file to make the table.
    # Maybe make both bed and bam files, and delete bed file later?
    echo "Making 2 bam files, so exiting merge_bamintersect.sh prior to calling filter_tsv.bash"
    exit 0
fi

if [ ${make_table} != "True" ]; then
    echo "LP integration table not required."
    exit 0
fi

# Initialize the counts output file with a header.
echo -e "chromosome-bam2\tStart_Pos\tEnd_Pos\tWidth\tNearest_Gene\tPost-filter_Reads\tchromosome-bam1\tSample" > "${sampleOutdir}/${sample_name}.counts.txt"
echo -e "chromosome-bam2\tStart_Pos\tEnd_Pos\tWidth\tNearest_Gene\tPost-filter_Reads\tchromosome-bam1\tSample" > "${sampleOutdir}/${sample_name}.informative.counts.txt"

# Make sure this has not been left over from a previous run.
rm -f "${sampleOutdir}/${sample_name}.informative.bed"

debug_fa "Starting the main_chrom while loop"

# Merge the output related to each bam1 chromosome:
while read main_chrom ; do
    debug_fa "Top of the main_chrom while loop"

    bedops --chrom ${main_chrom} --everything "${sampleOutdir}/${sample_name}.bed" > "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.bed"

    echo "Working on:  ${main_chrom}" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

    echo "Running filter_tsv without filters." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    ${src}/filter_tsv.sh ${sampleOutdir} ${sample_name} ${bam2genome} "NA" ${main_chrom} ${INTERMEDIATEDIR} ${integrationsite} \
                         all_reads_counts ${verbose}

    debug_fa "Completed filter_tsv for all_reads_counts"

    echo "Running filter_tsv with HA filters (if any) and with the Exclude Regions file (which may be empty)." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    ${src}/filter_tsv.sh ${sampleOutdir} ${sample_name} ${bam2genome} ${genome2exclude} ${main_chrom} \
                         ${INTERMEDIATEDIR} ${integrationsite} informative_reads_counts ${verbose}

    debug_fa "Completed filter_tsv for informative_reads_counts"

    cat "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.informative.bed" >> "${sampleOutdir}/${sample_name}.informative.bed"

    ## Number of reads
    n1=$(wc -l < "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.bed")
    ## Number of informative reads
    n2=$(wc -l < "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.informative.bed")

    echo -e "Mapped reads with unmapped mates from ${main_chrom}:\t${n1}\tof which\t${n2}\tpassed through the Exclude Regions filters and over the HAs (if any), and are potentially informative.\n" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

done < "${sampleOutdir}/log/${sample_name}.chrom_list1"

debug_fa "Out of the main_chrom while loop"


# Subset bam files for uploading UCSC custom tracks:
cut -f4 "${sampleOutdir}/${sample_name}.informative.bed" > "${INTERMEDIATEDIR}/${sample_name}.readNames.txt"

if [ -s "${INTERMEDIATEDIR}/${sample_name}.readNames.txt" ]; then
    debug_fa "${sample_name}.readNames.txt was found"

    samtools view -H ${bamname1} > "${INTERMEDIATEDIR}/${sample_name}.bam1_informative_reads.sam"
    samtools view -F ${BAM_E1} ${bamname1} | "${src}/subsetBAM.py" "${INTERMEDIATEDIR}/${sample_name}.readNames.txt" >> "${INTERMEDIATEDIR}/${sample_name}.bam1_informative_reads.sam"
    samtools view -b "${INTERMEDIATEDIR}/${sample_name}.bam1_informative_reads.sam" > "${INTERMEDIATEDIR}/${sample_name}.bam1_informative_reads.unsorted.bam"
    samtools sort -o "${sampleOutdir}/${sample_name}.bam1_informative_reads.bam" "${INTERMEDIATEDIR}/${sample_name}.bam1_informative_reads.unsorted.bam"
    samtools index "${sampleOutdir}/${sample_name}.bam1_informative_reads.bam"
	rm "${INTERMEDIATEDIR}/${sample_name}.bam1_informative_reads.unsorted.bam"
	
    debug_fa "Completed samtools lines for bam1"

    samtools view -H ${bamname2} > "${INTERMEDIATEDIR}/${sample_name}.bam2_informative_reads.sam"
    samtools view -F ${BAM_E2} ${bamname2} | "${src}/subsetBAM.py" "${INTERMEDIATEDIR}/${sample_name}.readNames.txt" >> "${INTERMEDIATEDIR}/${sample_name}.bam2_informative_reads.sam"
    samtools view -b "${INTERMEDIATEDIR}/${sample_name}.bam2_informative_reads.sam" > "${INTERMEDIATEDIR}/${sample_name}.bam2_informative_reads.unsorted.bam"
    samtools sort -o "${sampleOutdir}/${sample_name}.bam2_informative_reads.bam" "${INTERMEDIATEDIR}/${sample_name}.bam2_informative_reads.unsorted.bam"
    samtools index "${sampleOutdir}/${sample_name}.bam2_informative_reads.bam"
	rm "${INTERMEDIATEDIR}/${sample_name}.bam2_informative_reads.unsorted.bam"
		
    debug_fa "Completed samtools lines for bam2"

    echo "bam files with informative reads are at:" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    echo "    ${sampleOutdir}/${sample_name}.bam1_informative_reads.bam" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    echo "    ${sampleOutdir}/${sample_name}.bam2_informative_reads.bam" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
else
    debug_fa "${sample_name}.readNames.txt was not found"

    # No informative reads
    echo "For ${sample_name}, there are no bam files with informative reads." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
fi

debug_fa "Done with the samtools section."

#############################################################################
## Cleanup
rm -r ${INTERMEDIATEDIR}
#############################################################################

date
echo "Leaving merge_bamintersect.sh"

