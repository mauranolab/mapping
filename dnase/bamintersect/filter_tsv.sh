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

TMPDIR_CSV=`mktemp -d`   # TMPDIR_CSV has no trailing slash
echo "main_chrom: ${main_chrom}       TMPDIR_CSV is: ${TMPDIR_CSV}"

#####################################################################################
# Parse homology arm boundaries.  These are the regions beyond the edge of the HAs.
outer_HAs5p="${TMPDIR_CSV}/outer_HAs5p"
outer_HAs3p="${TMPDIR_CSV}/outer_HAs3p"

IFS=':' read -r chr range <<< "${bam1_5p_HA}"
IFS='-' read -r r1 r2 <<< "${range}"
echo "${chr}"$'\t'"${r1}"$'\t'"${r2}" > ${outer_HAs5p}

IFS=':' read -r chr range <<< "${bam1_3p_HA}"
IFS='-' read -r r1 r2 <<< "${range}"
echo "${chr}"$'\t'"${r1}"$'\t'"${r2}" > ${outer_HAs3p}

cp "${TMPDIR_CSV}/outer_HAs5p" "${TMPDIR_CSV}/outer_HAs"
cat "${TMPDIR_CSV}/outer_HAs3p" >> "${TMPDIR_CSV}/outer_HAs"

deletion_range_f="${TMPDIR_CSV}/deletion_range"
IFS=':' read -r chr range <<< "${deletion_range}"
IFS='-' read -r r1 r2 <<< "${range}"
echo "${chr}"$'\t'"${r1}"$'\t'"${r2}" > ${deletion_range_f}
#####################################################################################

file_1="${sampleOutdir}/${sample_name}.${main_chrom}.bed"  # The universe of all reads from a bam1 chromosome, in bed12 format
                                                           # [chr start end readID flag +/-] x 2

# Maybe we want to filter out some of the reads from our universe:
if [ ${exclude_regions_from_counts} = "NA" ];then
    # We're not deleting any reads in this scenario, so sort and move on.
    # The "NA" is assigned in merge_bamintersect.sh when the "exclude_regions_from_counts" input file field is blank.
    sort-bed ${file_1} > "${TMPDIR_CSV}/filter_csv.output"
else
    # We want to delete file_1 reads that overlap the ranges in this bed3 file.
    file_2=${exclude_regions_from_counts##*/}

    # Remove the ".bed", and append "_sorted.bed"
    file_3="${file_2%.bed}_sorted.bed"

    # Get rid of comment lines prior to sorting.
    grep -v '^#' ${exclude_regions_from_counts} | sort-bed - > "${TMPDIR_CSV}/${file_3}"   # This is a sorted bed3 file of the ranges we want to delete.

    # Delete reads that overlap the ranges defined in file_3 (which are defined with respect to bam1 coordinates).
    sort-bed ${file_1} | bedops --not-element-of 1 - "${TMPDIR_CSV}/${file_3}" > "${TMPDIR_CSV}/filter_csv.output"

    # Now delete reads that are in the HAs.
    #     filter_csv.output has 12 columns: [chr start end readID flag +/-] x 2, the bam1 data being in 1-6, and the bam2 data being in 7-12.
    #     The HA coordinates are with respect to bam2, so we need to switch the 6-column halves to use the bedops functions on this file.
    cut -f1-6 "${TMPDIR_CSV}/filter_csv.output" > "${TMPDIR_CSV}/filter_csv_bam1.output"
    cut -f7-12 "${TMPDIR_CSV}/filter_csv.output" > "${TMPDIR_CSV}/filter_csv_bam2.output"
    paste "${TMPDIR_CSV}/filter_csv_bam2.output" "${TMPDIR_CSV}/filter_csv_bam1.output" >  "${TMPDIR_CSV}/filter_csv_bamBoth.output"

    #     Sort by the bam2 reads, then keep only the bam2 reads (and their bam1 mates) that do not overlap the HAs:
    sort-bed "${TMPDIR_CSV}/filter_csv_bamBoth.output" | bedops --not-element-of 1 - "${TMPDIR_CSV}/outer_HAs" > "${TMPDIR_CSV}/filter_csv.output"

    #     Now put the surviving bam1 reads back on the left hand side, and sort:
    cut -f1-6 "${TMPDIR_CSV}/filter_csv.output" > "${TMPDIR_CSV}/filter_csv_bam2.output"
    cut -f7-12 "${TMPDIR_CSV}/filter_csv.output" > "${TMPDIR_CSV}/filter_csv_bam1.output"
    paste "${TMPDIR_CSV}/filter_csv_bam1.output" "${TMPDIR_CSV}/filter_csv_bam2.output" > "${TMPDIR_CSV}/filter_csv_bamBoth.output"
    sort-bed "${TMPDIR_CSV}/filter_csv_bamBoth.output" > "${TMPDIR_CSV}/filter_csv.output"

    cp "${TMPDIR_CSV}/filter_csv.output" "${sampleOutdir}/${sample_name}.${main_chrom}.informative.bed"
fi
# filter_csv.output contains the reads we want to keep, in bed12 format.

#####################################################################################
# Start building the output table:

# Get just the hg38 bed3 part, then sort it:
cut -f7-9 "${TMPDIR_CSV}/filter_csv.output" | sort-bed - > "${TMPDIR_CSV}/filter_csv.bed3"

# The below makes +/- 500 bps regions around each hg38 read, and merges them when possible. It will return bed3 output.
bedops --range 500 --merge "${TMPDIR_CSV}/filter_csv.bed3" > "${TMPDIR_CSV}/all_regions.bed"

# Find the gene closest to each region:
if [ ${annotationgenome} = "mm10" ]; then
    closest-features --delim "\t" --closest "${TMPDIR_CSV}/all_regions.bed" /vol/isg/annotation/bed/mm10/gencodev17/GencodevM17.gene.bed > "${TMPDIR_CSV}/tmp1.out"
else
    closest-features --delim "\t" --closest "${TMPDIR_CSV}/all_regions.bed" /vol/isg/annotation/bed/hg38/gencodev25/Gencodev25.gene.bed > "${TMPDIR_CSV}/tmp1.out"
fi

# Cut out just the gene name column:
cut -d $'\t' -f7 "${TMPDIR_CSV}/tmp1.out" > "${TMPDIR_CSV}/tmp2.out"

# Sometimes the gene name is blank (like when the region is not one of the usual chr types). 
# Replace the blank with a dash.
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
) "${TMPDIR_CSV}/tmp2.out" > "${TMPDIR_CSV}/tmp3.out"


# Add the gene name column to the bed3 output of "bedops --range".
paste -d $'\t' "${TMPDIR_CSV}/all_regions.bed" "${TMPDIR_CSV}/tmp3.out" > "${TMPDIR_CSV}/tmp4.out"

# Count the number of reads in each region
bedmap --count "${TMPDIR_CSV}/all_regions.bed" "${TMPDIR_CSV}/filter_csv.bed3" > "${TMPDIR_CSV}/tmp5.out"
paste -d $'\t' "${TMPDIR_CSV}/tmp4.out" "${TMPDIR_CSV}/tmp5.out" > "${TMPDIR_CSV}/tmp6.out"

# Compute the size of each region, and adjust for the extra 500 bps "bedops --range" added to each end:
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
) "${TMPDIR_CSV}/tmp6.out" > "${TMPDIR_CSV}/tmp7.out"

# Sort the table:
sort -V -t $'\t' -k 6,6rn -k 1,1 -k 2,2n "${TMPDIR_CSV}/tmp7.out" > "${TMPDIR_CSV}/tmp8a.out"

# Remove line items with only one read.
# grep -P  -v '\t'1$ "${TMPDIR_CSV}/tmp8a.out" > "${TMPDIR_CSV}/tmp8b.out" || true
grep  -v $'\t'1$ "${TMPDIR_CSV}/tmp8a.out" > "${TMPDIR_CSV}/tmp8b.out" || true
# mv "${TMPDIR_CSV}/tmp8a.out" "${TMPDIR_CSV}/tmp8b.out"   # ... or not.

# Append the bam1 chromosome name onto the last column.
sed "s/$/\t${main_chrom}/" "${TMPDIR_CSV}/tmp8b.out" > "${TMPDIR_CSV}/tmp8.out"

# Dump data into main output file.
if [ ${exclude_regions_from_counts} = "NA" ];then
   cat "${TMPDIR_CSV}/tmp8.out" >> "${sampleOutdir}/${sample_name}.counts.txt"
else
   cat "${TMPDIR_CSV}/tmp8.out" >> "${sampleOutdir}/${sample_name}.informative.counts.txt"
fi

#####################################################################################
#####################################################################################
# Find the number of "deletion_gene" reads that cross over the edges of the HAs.

# Get the boundaries of the "deletion_gene" reads.
grep "${deletion_gene}" "${TMPDIR_CSV}/tmp8.out" |
awk -f <(cat << "AWK_HEREDOC_03"
BEGIN{FS="\t"; OFS="\t"}
{
   print $1,$2,$3
}
END{}
AWK_HEREDOC_03
) > "${TMPDIR_CSV}/del_gen1" || true

cat ${outer_HAs5p} > "${TMPDIR_CSV}/outer_HAs"
cat ${outer_HAs3p} >> "${TMPDIR_CSV}/outer_HAs"

echo -e "${main_chrom} (Exclude_Regions file: ${exclude_regions_from_counts}):" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

while read -r line_in; do
    # IFS=$'\t' read -r  x y z rest_of_line <<< $line_in
    echo ${line_in} | sed 's/ /\t/g' > "${TMPDIR_CSV}/del_gen2"

    # Extract the "deletion_gene" reads from our universe of retained reads, built above.
    bedops --element-of 1 "${TMPDIR_CSV}/filter_csv.bed3" "${TMPDIR_CSV}/del_gen2" > "${TMPDIR_CSV}/del_gen_reads"

    # Now get the reads that are in the deletion range.
    bedops --element-of 1 "${TMPDIR_CSV}/del_gen_reads" "${deletion_range_f}" > "${TMPDIR_CSV}/del_range_reads_out"

    num_lines=$(wc -l < "${TMPDIR_CSV}/del_range_reads_out")
    echo -e "    Number of [${deletion_gene} ${line_in}] reads in the Deletion Range:\t${num_lines}" >> \
            "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
done < "${TMPDIR_CSV}/del_gen1"

echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

#####################################################################################
#####################################################################################
# For testing
# touch "${TMPDIR_CSV}/${main_chrom}"
# cp -r "${TMPDIR_CSV}" "${sampleOutdir}"

