#!/bin/bash
set -eu -o pipefail

#SBATCH --job-name=cleanup
#SBATCH --mail-type=END
#SBATCH --mail-user=cadley.nyulangone@gmail.com

## Variables are passed in via sbatch export.
##########################################################################################################
## Merge the bam_intersect output files:

dir_names=($(ls -d ${FINAL_OUTDIR}/TEMP_DIR_CL1/*/))      ## These have a trailing /

first="True"
for i in ${dir_names[@]}; do
    if [ ${first} = "True" ]; then
        cp "${i}dsgrep_out.csv" "${FINAL_OUTDIR}/f1.bed"
        first="False"
    else
        cat "${i}dsgrep_out.csv" >> "${FINAL_OUTDIR}/f1.bed"
    fi
done

date
echo Sorting...
sort-bed "${FINAL_OUTDIR}/f1.bed" > "${FINAL_OUTDIR}/${sample_name}.all.bed"

date
echo Deleting files and directories...

rm "${FINAL_OUTDIR}/f1.bed"
rm -r "${FINAL_OUTDIR}/TEMP_DIR_CL1"
rm -r "${FINAL_OUTDIR}/bams"

##########################################################################################################
## Create the final output tables.

if [ ${make_csv} != "True" ]; then
    # Need the bed file to make the table.
    # Maybe make both bed and bam files, and delete bed file later?
    echo "Making 2 bam files, so exiting cleanup.sbatch prior to calling filter_tsv.bash"
    exit 0
fi

if [ ${make_table} != "True" ]; then
    echo "LP integration table not required."
    exit 0
fi

echo -e "chromosome-bam2\tStart_Pos\tEnd_Pos\tWidth\tNearest_Gene\tInf Reads\tchromosome-bam1" >> "${FINAL_OUTDIR}/${sample_name}.counts.txt"

while read main_chrom ; do
    bedops --chrom ${main_chrom} --everything "${FINAL_OUTDIR}/${sample_name}.all.bed" > "${FINAL_OUTDIR}/${sample_name}.${main_chrom}.bed"

    ${src}/filter_tsv.bash ${FINAL_OUTDIR} ${bam1_5p_HA} ${bam1_3p_HA} ${sample_name} ${cegsgenome} ${annotationgenome} \
                       ${exclude_regions_from_counts} ${main_chrom}

    cat "${FINAL_OUTDIR}/${sample_name}.${main_chrom}.informative.bed" >> "${FINAL_OUTDIR}/${sample_name}.informative.bed"

    ## Number of reads
    n1=$(wc -l < "${FINAL_OUTDIR}/${sample_name}.${main_chrom}.bed")
    ## Number of informative reads
    n2=$(wc -l < "${FINAL_OUTDIR}/${sample_name}.${main_chrom}.informative.bed")

    echo -e "Reads from ${main_chrom}:\t${n1}\tof which\t${n2}\tare informative.\n" >> "${FINAL_OUTDIR}/${sample_name}.counts.anc_info.txt"

    # Don't need to keep these anymore.
    rm "${FINAL_OUTDIR}/${sample_name}.${main_chrom}.informative.bed"
    rm "${FINAL_OUTDIR}/${sample_name}.${main_chrom}.bed"
done < "${FINAL_OUTDIR}/log/chrom_list1"

date
echo "Done with cleanup.sbatch"

