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

if [ ${make_csv} != "True" ]; then
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
    bedops --chrom ${main_chrom} --everything "${sampleOutdir}/${sample_name}.bed" > "${sampleOutdir}/${sample_name}.${main_chrom}.bed"

    ${src}/filter_tsv.sh ${sampleOutdir} ${bam1_5p_HA} ${bam1_3p_HA} ${sample_name} ${cegsgenome} ${annotationgenome} "NA" ${main_chrom} ${deletion_gene} ${deletion_range}

    ${src}/filter_tsv.sh ${sampleOutdir} ${bam1_5p_HA} ${bam1_3p_HA} ${sample_name} ${cegsgenome} ${annotationgenome} \
                       ${exclude_regions_from_counts} ${main_chrom} ${deletion_gene} ${deletion_range}

    cat "${sampleOutdir}/${sample_name}.${main_chrom}.informative.bed" >> "${sampleOutdir}/${sample_name}.informative.bed"

    ## Number of reads
    n1=$(wc -l < "${sampleOutdir}/${sample_name}.${main_chrom}.bed")
    ## Number of informative reads
    n2=$(wc -l < "${sampleOutdir}/${sample_name}.${main_chrom}.informative.bed")

    # Don't need to keep these.
    rm "${sampleOutdir}/${sample_name}.${main_chrom}.informative.bed"
    rm "${sampleOutdir}/${sample_name}.${main_chrom}.bed"

    echo -e "Mapped reads with unmapped mates from ${main_chrom}:\t${n1}\tof which\t${n2}\tpassed through the filters and HAs, and are potentially informative.\n" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
done < "${sampleOutdir}/log/chrom_list1"

#############################################################################
## Cleanup
rm -r ${INTERMEDIATEDIR}
#############################################################################

date
echo "Done with merge_bamintersect.sh"

