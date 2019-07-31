#!/bin/bash
set -eu -o pipefail

FINAL_OUTDIR=$1
bam1_5p_HA=$2
bam1_3p_HA=$3
sample_name=$4
cegsgenome=$5
annotationgenome=$6
exclude_regions_from_counts=$7
main_chrom=$8

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
#####################################################################################

file_1="${FINAL_OUTDIR}/${sample_name}.${main_chrom}.bed"  # The universe of all reads from a bam1 chromosome, in bed12 format
                                                           # [chr start end readID flag +/-] x 2

# Maybe we want to filter out some of the reads from our universe:
if [ ${exclude_regions_from_counts} = "NA" ];then
    # We're not deleting any reads in this scenario, so sort and move on.
    # The "NA" is assigned in cleanup.sbatch when the "exclude_regions_from_counts" input file field is blank.
    sort-bed ${file_1} > "${TMPDIR_CSV}/filter_csv.output"
else
    # We want to delete file_1 reads that overlap the ranges in this bed3 file.
    file_2=${exclude_regions_from_counts##*/}

    # Remove the ".bed", and append "_sorted.bed"
    file_3="${file_2%.bed}_sorted.bed"

    # Get rid of comment lines prior to sorting.
    grep -v '^#' ${exclude_regions_from_counts} | sort-bed - > "${TMPDIR_CSV}/${file_3}"   # This is a sorted bed3 file of the ranges we want to delete.

    # Delete reads that overlap the ranges defined in file_3.
    sort-bed ${file_1} | bedops --not-element-of 1 - "${TMPDIR_CSV}/${file_3}" > "${TMPDIR_CSV}/filter_csv.output"
fi
# filter_csv.output contains the reads we want to keep, in bed12 format.

cp "${TMPDIR_CSV}/filter_csv.output" "${FINAL_OUTDIR}/${sample_name}.${main_chrom}.informative.bed"

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

# Compute the size of each region:
awk -f <(cat << "AWK_HEREDOC_01"
BEGIN{FS="\t"; OFS="\t"}
{
   width = $3 - $2
   print $1, $2, $3, width, $4, $5
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
cat "${TMPDIR_CSV}/tmp8.out" >> "${FINAL_OUTDIR}/${sample_name}.counts.txt"

#####################################################################################
#####################################################################################
# Find the number of HPRT1 reads that cross over the edges of the HAs.

# Get the boundaries of the HPRT1 reads.
grep HPRT1 "${TMPDIR_CSV}/tmp8.out" |
awk -f <(cat << "AWK_HEREDOC_03"
BEGIN{FS="\t"; OFS="\t"}
{
   print $1,$2,$3
}
END{}
AWK_HEREDOC_03
) > "${TMPDIR_CSV}/hprt1" || true

# Need to sort them for bedops.
sort-bed "${TMPDIR_CSV}/hprt1" > "${TMPDIR_CSV}/hprt2"

# Extract the HPRT1 reads from our universe of retained reads, built above.
bedops --element-of 1 "${TMPDIR_CSV}/filter_csv.bed3" "${TMPDIR_CSV}/hprt2" > "${TMPDIR_CSV}/hprt_reads"

# Now get the reads that cross the HA edges.
bedops --element-of 1 "${TMPDIR_CSV}/hprt_reads" ${outer_HAs5p} > "${TMPDIR_CSV}/hprt_reads5p"
bedops --element-of 1 "${TMPDIR_CSV}/hprt_reads" ${outer_HAs3p} > "${TMPDIR_CSV}/hprt_reads3p"

num_lines_5p=$(wc -l < "${TMPDIR_CSV}/hprt_reads5p")
num_lines_3p=$(wc -l < "${TMPDIR_CSV}/hprt_reads3p")

echo -e "${main_chrom}:" >> "${FINAL_OUTDIR}/${sample_name}.counts.anc_info.txt"
echo -e "    Number of HPRT1 5p HA reads:\t${num_lines_5p}" >> "${FINAL_OUTDIR}/${sample_name}.counts.anc_info.txt"
echo -e "    Number of HPRT1 3p HA reads:\t${num_lines_3p}" >> "${FINAL_OUTDIR}/${sample_name}.counts.anc_info.txt"
#####################################################################################
#####################################################################################
# For testing
# cp -r "${TMPDIR_CSV}" "${FINAL_OUTDIR}"
