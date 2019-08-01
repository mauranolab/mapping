#!/bin/bash
set -eu -o pipefail

#SBATCH --mail-type=END
#SBATCH --mail-user=cadley.nyulangone@gmail.com

## Variables are passed in via sbatch export.
##########################################################################################################
## Merge the bam_intersect output files:

dir_names=($(ls -d ${sampleOutdir}/TEMP_DIR_CL1/*/))      ## These have a trailing /

first="True"
for i in ${dir_names[@]}; do
    if [ ${first} = "True" ]; then
        cp "${i}dsgrep_out.csv" "${sampleOutdir}/f1.bed"
        first="False"
    else
        cat "${i}dsgrep_out.csv" >> "${sampleOutdir}/f1.bed"
    fi
done

date
echo Sorting...
sort-bed "${sampleOutdir}/f1.bed" > "${sampleOutdir}/${sample_name}.bed"

date
echo Deleting files and directories...

rm "${sampleOutdir}/f1.bed"
rm -r "${sampleOutdir}/TEMP_DIR_CL1"
rm -r "${sampleOutdir}/bams"

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

# Clear out some old output files, if any.
rm -f "${sampleOutdir}/${sample_name}.informative.bed"
rm -f "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

# Initialize the counts output file with a header.
echo -e "chromosome-bam2\tStart_Pos\tEnd_Pos\tWidth\tNearest_Gene\tInf Reads\tchromosome-bam1" > "${sampleOutdir}/${sample_name}.counts.txt"

# Merge the output related to each bam1 chromosome:
while read main_chrom ; do
    bedops --chrom ${main_chrom} --everything "${sampleOutdir}/${sample_name}.bed" > "${sampleOutdir}/${sample_name}.${main_chrom}.bed"

    ${src}/filter_tsv.sh ${sampleOutdir} ${bam1_5p_HA} ${bam1_3p_HA} ${sample_name} ${cegsgenome} ${annotationgenome} \
                       ${exclude_regions_from_counts} ${main_chrom}

    cat "${sampleOutdir}/${sample_name}.${main_chrom}.informative.bed" >> "${sampleOutdir}/${sample_name}.informative.bed"

    ## Number of reads
    n1=$(wc -l < "${sampleOutdir}/${sample_name}.${main_chrom}.bed")
    ## Number of informative reads
    n2=$(wc -l < "${sampleOutdir}/${sample_name}.${main_chrom}.informative.bed")

    echo -e "Reads from ${main_chrom}:\t${n1}\tof which\t${n2}\tare informative.\n" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

    # Don't need to keep these anymore.
    rm "${sampleOutdir}/${sample_name}.${main_chrom}.informative.bed"
    rm "${sampleOutdir}/${sample_name}.${main_chrom}.bed"
done < "${sampleOutdir}/log/chrom_list1"

date
echo "Done with merge_bamintersect.sh"

