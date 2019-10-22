#!/bin/bash
set -eu -o pipefail

sampleOutdir=$1
sample_name=$2
bam2genome=$3
genome2exclude=$4
main_chrom=$5
INTERMEDIATEDIR=$6
integrationsite=$7
counts_type=$8
verbose=$9

#####################################################################################
## This function implements verbose output for debugging purposes.
debug_fa() {
    if [ ${verbose} = "True" ]; then
        echo "${1}"
    fi
}
#####################################################################################

echo "main_chrom: ${main_chrom}       TMPDIR is: ${TMPDIR}"

filter_csv_output="${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.informative.bed"

# Maybe we want to filter out some of the reads from our universe:
if [ ${genome2exclude} = "NA" ];then
    # We're not deleting any reads in this scenario, so copy to the output file and move on.
    # "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.bed" = The universe of all reads from a bam1 chromosome, in bed12 format
                                                              # [chr start end readID flag +/-] x 2

    # Don't use mv, because we need this file to survive for the second call to filter_tsv, in which execution will drop into the below "else" section.
    cp "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.bed" "${filter_csv_output}"
else
    # Comments for the 5 piped lines below:
    # 1) Delete reads that overlap the ranges defined in genome2exclude (which are defined with respect to bam1 coordinates).
    # 2) Now delete reads that are in the HAs.
    #     2a) The output of step 1 has 12 columns: [chr start end readID flag +/-] x 2, the bam1 data being in 1-6, and the bam2 data being in 7-12.
    #         The HA coordinates are with respect to bam2, so we need to switch the 6-column halves to use the bedops functions on this file.
    # 3) Then sort by the bam2 reads, then keep only the bam2 reads (and their bam1 mates) that do not overlap the genome1exclude's (HAs):
    # 4) Now put the surviving bam1 reads back on the left hand side,
    # 5) and sort.

    bedops --not-element-of 1 "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.bed" "${genome2exclude}" | \
    awk 'BEGIN {OFS="\t"} ; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | \
    sort-bed - | bedops --not-element-of 1 - "${INTERMEDIATEDIR}/genome1exclude.bed" | \
    awk 'BEGIN {OFS="\t"} ; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | \
    sort-bed - > "${filter_csv_output}"
fi
# filter_csv.output contains the reads we want to keep, in bed12 format.

echo "${main_chrom} (Exclude_Regions file is $(basename ${genome2exclude})):" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

#####################################################################################
# Start building the output table:

# Get just the hg38/mm10 bed3 part, then sort it:
cut -f7-9 "${filter_csv_output}" | sort-bed - > "${TMPDIR}/sorted_speciesReadCoords.bed3"

debug_fa "Making 500 bps ranges with bedops."

# The below makes +/- 500 bps regions around each hg38/mm10 read, and merges them when possible. It will return bed3 output.
bedops --range 500 --merge "${TMPDIR}/sorted_speciesReadCoords.bed3" > "${TMPDIR}/all_regions.bed"

debug_fa "Starting call to closest-features"

# Find the gene closest to each region:
if [ ${bam2genome} = "mm10" ]; then
    closest-features --delim "\t" --closest "${TMPDIR}/all_regions.bed" /vol/isg/annotation/bed/mm10/gencodev17/GencodevM17.gene.bed > "${TMPDIR}/reads_and_geneNames.bed"
elif [ ${bam2genome} = "rn6" ]; then
    closest-features --delim "\t" --closest "${TMPDIR}/all_regions.bed" /vol/isg/annotation/bed/rn6/ensembl96/Ensemblv96_Rnor.gene.bed > "${TMPDIR}/reads_and_geneNames.bed"
else
    closest-features --delim "\t" --closest "${TMPDIR}/all_regions.bed" /vol/isg/annotation/bed/hg38/gencodev25/Gencodev25.gene.bed > "${TMPDIR}/reads_and_geneNames.bed"
fi

debug_fa "closest-features complete."

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

# Append the bam1 chromosome name and the short sample name onto the last column.
IFS='.' read -r -a array <<< "${sample_name}"    # Needed because: sample_name="${sample_name}.${bam1genome}_vs_${bam2genome}"
short_sample_name="${array[0]}"
sed "s/$/\t${main_chrom}\t${short_sample_name}/" "${TMPDIR}/short_sorted_table.bed" > "${TMPDIR}/short_sorted_table_chromName.bed"

# Dump data into main output file.
if [ ${counts_type} = "all_reads_counts" ];then
   cat "${TMPDIR}/short_sorted_table_chromName.bed" >> "${sampleOutdir}/${sample_name}.counts.txt"
else
   # Here counts_type == "informative_reads_counts"
   if [ ${genome2exclude} = "NA" ]; then
       echo "No HAs are available for this scenario"  >> "${sampleOutdir}/${sample_name}.informative.counts.txt"
   else
       cat "${TMPDIR}/short_sorted_table_chromName.bed" >> "${sampleOutdir}/${sample_name}.informative.counts.txt"
   fi
fi

#####################################################################################
# Find the number of reads in the deletion range.

bedops --element-of 1 "${TMPDIR}/sorted_speciesReadCoords.bed3" "${INTERMEDIATEDIR}/deletion_range.bed" > \
                      "${TMPDIR}/del_range_reads_out.bed"

num_lines=$(wc -l < "${TMPDIR}/del_range_reads_out.bed")
echo -e "    Number of reads in the Deletion Range:\t${num_lines}" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

echo "Leaving filter_tsv"

#####################################################################################
#####################################################################################
# For testing
# touch "${TMPDIR}/${main_chrom}"
# cp -r "${TMPDIR}" "${sampleOutdir}"

