#!/bin/bash
set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'


## This function implements verbose output for debugging purposes.
debug_fa() {
    if [ ${verbose} = "True" ]; then
        echo "[DEBUG] ${1}"
    fi
}
#####################################################################################


sampleOutdir=$1
sample_name=$2
bam2genome=$3
genome2exclude=$4
main_chrom=$5
INTERMEDIATEDIR=$6
integrationsite=$7
counts_type=$8
verbose=$9


debug_fa "main_chrom: ${main_chrom}       TMPDIR is: ${TMPDIR}"

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
    
    bedops -n 1 "${INTERMEDIATEDIR}/${sample_name}.${main_chrom}.bed" "${genome2exclude}" | \
    awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | \
    sort-bed - | bedops -n 1 - "${INTERMEDIATEDIR}/genome1exclude.bed" | \
    awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | \
    sort-bed - > "${filter_csv_output}"
fi
# filter_csv.output contains the reads we want to keep, in bed12 format.

echo "${main_chrom} (Exclude_Regions file is $(basename ${genome2exclude})):" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

#####################################################################################
# Start building the output table:

# Get just the hg38/mm10 bed3 part, then sort it:
cut -f7-9 "${filter_csv_output}" | sort-bed - > ${TMPDIR}/sorted_speciesReadCoords.bed3


# Find the gene closest to each region:
case "${bam2genome}" in
mm10)
    geneAnnotationFile=/vol/isg/annotation/bed/mm10/gencodev17/GencodevM17.gene.bed
    ;;
rn6)
    geneAnnotationFile=/vol/isg/annotation/bed/rn6/ensembl96/Ensemblv96_Rnor.gene.bed
    ;;
hg38)
    geneAnnotationFile=/vol/isg/annotation/bed/hg38/gencodev31/Gencodev31.gene.bed
    ;;
*)
    echo "ERROR: Don't recognize genome ${bam2genome}";
    exit 1;;
esac


# The below makes +/- 500 bps regions around each hg38/mm10 read, and merges them when possible. It will return bed3 output.
bedops --range 500 -m ${TMPDIR}/sorted_speciesReadCoords.bed3 > ${TMPDIR}/all_regions.bed

awk -F "\t" 'BEGIN {OFS="\t"} $7=="protein_coding"' ${geneAnnotationFile} |
closest-features --delim "|" --closest ${TMPDIR}/all_regions.bed - |
# Sometimes the gene name is blank (like when the region is not one of the usual chr types); Replace the blank with a dash.
awk -F "|" 'BEGIN {OFS="\t"} {split($2, gene, "\t"); if(gene[4]=="") {gene[4]="-"} print $1, gene[4]}' |

# Count the number of reads in each region
# Then compute the size of each region, and adjust for the extra 500 bps "bedops --range" added to each end
# Then remove line items with only one read
bedmap --delim "\t" --echo --count - ${TMPDIR}/sorted_speciesReadCoords.bed3 |
awk -v minReadsCutoff=2 -F "\t" 'BEGIN {OFS="\t"} \
$5>=minReadsCutoff { \
   chromStart = $2 + 500; \
   chromEnd = $3 - 500; \
   width = $3 - $2 - 1000; \
   print $1, chromStart, chromEnd, width, $4, $5; \
}' |
#Sort by Post-filter_Reads, then chromosome-bam2, then Start_Pos
sort -V -t $'\t' -k 6,6rn -k 1,1 -k 2,2n - |
# Append the bam1 chromosome name and the short sample name onto the last column.
awk -v short_sample_name=`echo "${sample_name}" | cut -d "." -f1` -v main_chrom=${main_chrom} -F "\t" 'BEGIN {OFS="\t"} {print $0, main_chrom, short_sample_name}' > ${TMPDIR}/short_sorted_table.bed

if [ ${counts_type} = "all_reads_counts" ];then
    cat ${TMPDIR}/short_sorted_table.bed >> ${sampleOutdir}/${sample_name}.counts.txt
else
    cat ${TMPDIR}/short_sorted_table.bed >> ${sampleOutdir}/${sample_name}.informative.counts.txt
fi


#####################################################################################
# Find the number of reads in the deletion range.

bedops -e 1 "${TMPDIR}/sorted_speciesReadCoords.bed3" "${INTERMEDIATEDIR}/deletion_range.bed" > "${TMPDIR}/del_range_reads_out.bed"

num_lines=$(wc -l < "${TMPDIR}/del_range_reads_out.bed")
echo -e "    Number of reads in the Deletion Range:\t${num_lines}" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

echo "Finished filter_tsv.sh"
date

#####################################################################################
#####################################################################################
# For testing
# touch "${TMPDIR}/${main_chrom}"
# cp -r "${TMPDIR}" "${sampleOutdir}"

