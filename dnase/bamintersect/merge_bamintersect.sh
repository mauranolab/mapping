#!/bin/bash
set -eu -o pipefail

## Variables are passed in via sbatch export.
##########################################################################################################
## Merge the bam_intersect output files:

dir_names=($(ls -d ${INTERMEDIATEDIR}/bamintersectPyOut/*/))      ## These have a trailing /

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
echo -e "chromosome-bam2\tStart_Pos\tEnd_Pos\tWidth\tNearest_Gene\tPost-filter_Reads\tchromosome-bam1" > "${sampleOutdir}/${sample_name}.counts.txt"
echo -e "chromosome-bam2\tStart_Pos\tEnd_Pos\tWidth\tNearest_Gene\tPost-filter_Reads\tchromosome-bam1" > "${sampleOutdir}/${sample_name}.informative.counts.txt"

# Make sure this has not been left over from a previous run.
rm -f "${sampleOutdir}/${sample_name}.informative.bed"

# Merge the output related to each bam1 chromosome:
while read main_chrom ; do
    bedops --chrom ${main_chrom} --everything "${sampleOutdir}/${sample_name}.bed" > "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.bed"

    echo "Running filter_tsv without filters." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    ${src}/filter_tsv.sh ${sampleOutdir} ${sample_name} ${bam2genome} "NA" ${main_chrom} ${INTERMEDIATEDIR} ${integrationsite} all_reads_counts

    echo "Running filter_tsv with HA filters (if any) and with the Exclude Regions file (which may be empty)." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    ${src}/filter_tsv.sh ${sampleOutdir} ${sample_name} ${bam2genome} ${genome2exclude} ${main_chrom} \
                         ${INTERMEDIATEDIR} ${integrationsite} informative_reads_counts

    cat "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.informative.bed" >> "${sampleOutdir}/${sample_name}.informative.bed"

    ## Number of reads
    n1=$(wc -l < "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.bed")
    ## Number of informative reads
    n2=$(wc -l < "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.informative.bed")

    echo -e "Mapped reads with unmapped mates from ${main_chrom}:\t${n1}\tof which\t${n2}\tpassed through the Exclude Regions filters and over the HAs (if any), and are potentially informative.\n" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"


    # Subset bam files for UCSC loading:
    cut -f4 "${sampleOutdir}/${sample_name}.informative.bed" > "${INTERMEDIATEDIR}/${sample_name}.readNames.txt"

    if [ -s "${INTERMEDIATEDIR}/${sample_name}.readNames.txt" ]; then
        samtools view -H ${bamname1} > "${INTERMEDIATEDIR}/${sample_name}.bam1_informative_reads.sam"
        samtools view ${bamname1} | grep --file="${INTERMEDIATEDIR}/${sample_name}.readNames.txt" >> "${INTERMEDIATEDIR}/${sample_name}.bam1_informative_reads.sam"
        samtools view -b "${INTERMEDIATEDIR}/${sample_name}.bam1_informative_reads.sam" > "${sampleOutdir}/${sample_name}.bam1_informative_reads.bam"
        samtools index "${sampleOutdir}/${sample_name}.bam1_informative_reads.bam"
    
        samtools view -H ${bamname2} > "${INTERMEDIATEDIR}/${sample_name}.bam2_informative_reads.sam"
        samtools view ${bamname2} | grep --file="${INTERMEDIATEDIR}/${sample_name}.readNames.txt" >> "${INTERMEDIATEDIR}/${sample_name}.bam2_informative_reads.sam"
        samtools view -b "${INTERMEDIATEDIR}/${sample_name}.bam2_informative_reads.sam" > "${sampleOutdir}/${sample_name}.bam2_informative_reads.bam"
        samtools index "${sampleOutdir}/${sample_name}.bam2_informative_reads.bam"

        echo "bam files with informative reads are at:" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
        echo "    ${sampleOutdir}/${sample_name}.bam1_informative_reads.bam" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
        echo "    ${sampleOutdir}/${sample_name}.bam2_informative_reads.bam" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
        echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    else
        # No informative reads
        echo "For ${sample_name}, there are no bam files with informative reads." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
        echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    fi

done < "${sampleOutdir}/log/${sample_name}.chrom_list1"

#############################################################################
## Cleanup
rm -r ${INTERMEDIATEDIR}
#############################################################################

date
echo "Done with merge_bamintersect.sh"

