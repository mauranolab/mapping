#!/bin/bash
set -eu -o pipefail

sampleOutdir=$1
bam1_5p_HA=$2
bam1_3p_HA=$3
sample_name=$4
cegsgenome=$5
annotationgenome=$6
exclude_regions_from_counts=$7
main_chrom=$8
deletion_gene=$9
deletion_range=${10}
INTERMEDIATEDIR=${11}

echo "main_chrom: ${main_chrom}       TMPDIR is: ${TMPDIR}"

#####################################################################################
# Parse homology arm boundaries.  These are the regions beyond the edge of the HAs.
outer_HAs5p="${TMPDIR}/outer_HAs5p.bed"
outer_HAs3p="${TMPDIR}/outer_HAs3p.bed"

IFS=':' read -r chr range <<< "${bam1_5p_HA}"
IFS='-' read -r r1 r2 <<< "${range}"
echo "${chr}"$'\t'"${r1}"$'\t'"${r2}" > ${outer_HAs5p}

IFS=':' read -r chr range <<< "${bam1_3p_HA}"
IFS='-' read -r r1 r2 <<< "${range}"
echo "${chr}"$'\t'"${r1}"$'\t'"${r2}" > ${outer_HAs3p}

bedops -u ${outer_HAs5p} ${outer_HAs3p} > "${TMPDIR}/outer_HAs.bed"

deletion_range_f="${TMPDIR}/deletion_range.bed"
IFS=':' read -r chr range <<< "${deletion_range}"
IFS='-' read -r r1 r2 <<< "${range}"
echo "${chr}"$'\t'"${r1}"$'\t'"${r2}" > ${deletion_range_f}
#####################################################################################

filter_csv_output="${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.informative.bed"

# Maybe we want to filter out some of the reads from our universe:
if [ ${exclude_regions_from_counts} = "NA" ];then
    # We're not deleting any reads in this scenario, so sort and move on.
    # The "NA" is assigned in merge_bamintersect.sh when the "exclude_regions_from_counts" input file field is blank.
    # "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.bed" = The universe of all reads from a bam1 chromosome, in bed12 format
                                                        # [chr start end readID flag +/-] x 2
    sort-bed "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.bed" > "${filter_csv_output}"
else
    # We want to delete "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.bed" reads that overlap the ranges in this bed3 file.
    exclude_regions_filename=${exclude_regions_from_counts##*/}

    # Remove the ".bed", and append "_sorted.bed"
    sorted_exclude_regions_filename="${exclude_regions_filename%.bed}_sorted.bed"

    # Get rid of comment lines prior to sorting.
    grep -v '^#' ${exclude_regions_from_counts} | sort-bed - > "${TMPDIR}/${sorted_exclude_regions_filename}"   # This is a sorted bed3 file of the ranges we want to delete.


    # Comments for the 5 piped lines below:
    # 1) Delete reads that overlap the ranges defined in sorted_exclude_regions_filename (which are defined with respect to bam1 coordinates).
    # 2) Now delete reads that are in the HAs.
    #     2a) The output of step 1 has 12 columns: [chr start end readID flag +/-] x 2, the bam1 data being in 1-6, and the bam2 data being in 7-12.
    #         The HA coordinates are with respect to bam2, so we need to switch the 6-column halves to use the bedops functions on this file.
    # 3) Then sort by the bam2 reads, then keep only the bam2 reads (and their bam1 mates) that do not overlap the HAs:
    # 4) Now put the surviving bam1 reads back on the left hand side,
    # 5) and sort.

    sort-bed "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.bed" | bedops --not-element-of 1 - "${TMPDIR}/${sorted_exclude_regions_filename}" | \
    awk 'BEGIN {OFS="\t"} ; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | \
    sort-bed - | bedops --not-element-of 1 - "${TMPDIR}/outer_HAs.bed" | \
    awk 'BEGIN {OFS="\t"} ; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | \
    sort-bed - > "${filter_csv_output}"
fi
# filter_csv.output contains the reads we want to keep, in bed12 format.

#####################################################################################
# Start building the output table:

# Get just the hg38/mm10 bed3 part, then sort it:
cut -f7-9 "${filter_csv_output}" | sort-bed - > "${TMPDIR}/sorted_speciesReadCoords.bed3"

# The below makes +/- 500 bps regions around each hg38/mm10 read, and merges them when possible. It will return bed3 output.
bedops --range 500 --merge "${TMPDIR}/sorted_speciesReadCoords.bed3" > "${TMPDIR}/all_regions.bed"

# Find the gene closest to each region:
if [ ${annotationgenome} = "mm10" ]; then
    closest-features --delim "\t" --closest "${TMPDIR}/all_regions.bed" /vol/isg/annotation/bed/mm10/gencodev17/GencodevM17.gene.bed > "${TMPDIR}/reads_and_geneNames.bed"
else
    closest-features --delim "\t" --closest "${TMPDIR}/all_regions.bed" /vol/isg/annotation/bed/hg38/gencodev25/Gencodev25.gene.bed > "${TMPDIR}/reads_and_geneNames.bed"
fi


# Cut out just the gene name column, then look for blank gene names.
# Sometimes the gene name is blank (like when the region is not one of the usual chr types). 
# Replace the blank with a dash.
# Then add the gene name column to the bed3 output of "bedops --range".
cut -d $'\t' -f7 "${TMPDIR}/reads_and_geneNames.bed" |
awk -f <(cat << "AWK_HEREDOC_01"
BEGIN{}
{
   if($0=="")
   {
      print "-"
   }else
   {
      print $0
   }
}
END{}
AWK_HEREDOC_01
) |
paste -d $'\t' "${TMPDIR}/all_regions.bed" - > "${TMPDIR}/all_regions_and_geneNames.bed"

# Count the number of reads in each region
# Then compute the size of each region, and adjust for the extra 500 bps "bedops --range" added to each end
# Then sort the table
# Then remove line items with only one read
bedmap --count "${TMPDIR}/all_regions.bed" "${TMPDIR}/sorted_speciesReadCoords.bed3" |
paste -d $'\t' "${TMPDIR}/all_regions_and_geneNames.bed" - |
awk -f <(cat << "AWK_HEREDOC_01"
BEGIN{FS="\t"; OFS="\t"}
{
   chromStart = $2 + 500
   chromEnd = $3 - 500
   width = $3 - $2 - 1000
   print $1, chromStart, chromEnd, width, $4, $5
}
END{}
AWK_HEREDOC_01
) |
sort -V -t $'\t' -k 6,6rn -k 1,1 -k 2,2n - |
grep -v $'\t'1$ - > "${TMPDIR}/short_sorted_table.bed" || true

# Append the bam1 chromosome name onto the last column.
sed "s/$/\t${main_chrom}/" "${TMPDIR}/short_sorted_table.bed" > "${TMPDIR}/short_sorted_table_chromName.bed"

# Dump data into main output file.
if [ ${exclude_regions_from_counts} = "NA" ];then
   cat "${TMPDIR}/short_sorted_table_chromName.bed" >> "${sampleOutdir}/${sample_name}.counts.txt"
else
   cat "${TMPDIR}/short_sorted_table_chromName.bed" >> "${sampleOutdir}/${sample_name}.informative.counts.txt"
fi

#####################################################################################
#####################################################################################
# Find the number of "deletion_gene" reads that cross over the edges of the HAs.

# Get the boundaries of the "deletion_gene" reads.
grep "${deletion_gene}" "${TMPDIR}/short_sorted_table_chromName.bed" |
awk -f <(cat << "AWK_HEREDOC_03"
BEGIN{FS="\t"; OFS="\t"}
{
   print $1,$2,$3
}
END{}
AWK_HEREDOC_03
) > "${TMPDIR}/deletion_gene_boundary_reads.bed" || true

echo -e "${main_chrom} (Exclude_Regions file: ${exclude_regions_from_counts}):" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

while read -r line_in; do
    # Extract the "deletion_gene" reads from our universe of retained reads, built above.
    # Then get the reads that are in the deletion range.
    echo ${line_in} | sed 's/ /\t/g' |
    bedops --element-of 1 "${TMPDIR}/sorted_speciesReadCoords.bed3" - |
    bedops --element-of 1 - "${deletion_range_f}" > "${TMPDIR}/del_range_reads_out"

    num_lines=$(wc -l < "${TMPDIR}/del_range_reads_out")
    echo -e "    Number of [${deletion_gene} ${line_in}] reads in the Deletion Range:\t${num_lines}" >> \
            "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
done < "${TMPDIR}/deletion_gene_boundary_reads.bed"

echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

#####################################################################################
#####################################################################################
# For testing
# touch "${TMPDIR}/${main_chrom}"
# cp -r "${TMPDIR}" "${sampleOutdir}"

