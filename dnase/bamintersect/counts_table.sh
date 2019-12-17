#!/bin/bash

set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'

#####################################################################################
# counts_table.sh takes a bed12 file as input. This file is generated by bamintersect.py
# The 12 columns are:
#    [chr start end readID flag +/-] x 2, bam1 data being in 1-6, and bam2 data being in 7-12.
# The rows are individual reads.
#
# Lines in the bed12 file are
# aggregated to their nearest genes. Some additional values are calculated, such as quantity of
# reads, widths, and the range of related bam1 reads.
#
# The output file will be named ${OUTBASE}.counts.txt
#
# sampleName is in the form:       "<sample_name>.<bam1genome>_vs_<bam2genome>"
#
#####################################################################################

OUTBASE=$1
sample_name=$2
bam2genome=$3
bamintersectBED12=$4
INTERMEDIATEDIR=$5

echo "Making counts_table from ${bamintersectBED12}; will output to ${OUTBASE}.counts.txt"


echo -e -n "Number of total reads:\t"
cat ${bamintersectBED12} | wc -l



# Set the appropriate gene annotation file.
# Note that hg38_full, etc. was converted to just "hg38" prior to calling bamintersect.
case "${bam2genome}" in
mm10)
    geneAnnotationFile=/vol/isg/annotation/bed/mm10/gencodev23/GencodevM23.gene.bed
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


# Initialize the counts output files with a header.
echo -e "#chrom_bam2\tchromStart_bam2\tchromEnd_bam2\tWidth_bam2\tNearestGene_bam2\tPost-filter_Reads\tchrom_bam1\tchromStart_bam1\tchromEnd_bam1\tWidth_bam1\tSample" > ${OUTBASE}.counts.txt

# Loop over the bam1 chromosome names.
cut -f1 ${bamintersectBED12} | sort | uniq > "${TMPDIR}/bam1_chroms"

while read main_chrom ; do
    # Get just the hg38/mm10/rn6 bed3 part, then sort it:
    awk -v mainChrom="${main_chrom}" -F "\t" 'BEGIN {OFS="\t"} $1==mainChrom {print $0; }' ${bamintersectBED12} |
    cut -f7-9 | sort-bed - > ${TMPDIR}/sorted_speciesReadCoords.bed3
    
    # The below makes +/- 500 bps regions around each hg38/mm10 read, and merges them when possible. It will return bed3 output.
    bedops --range 500 -m ${TMPDIR}/sorted_speciesReadCoords.bed3 > ${TMPDIR}/all_regions.bed
    
    cat ${geneAnnotationFile} |
    awk -F "\t" 'BEGIN {OFS="\t"} $7=="protein_coding"' |
    #Restrict to level 1 or 2 genes if available (Sox2 is level 2 for example)
    awk -F "\t" 'BEGIN {OFS="\t"} $5<=2 || $5=="."' |
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
       #BUGBUG temporary fix. doesnt properly revert coords if chromStart < 500
       if(chromEnd < chromStart) {chromStart = chromEnd-1} \
       width = $3 - $2 - 1000; \
       print $1, chromStart, chromEnd, width, $4, $5; \
    }' |
    #Sort by Post-filter_Reads, then chromosome-bam2, then Start_Pos
    sort -V -t $'\t' -k 6,6rn -k 1,1 -k 2,2n - |
    # Append the bam1 chromosome name and the short sample name onto the last column.
    awk -v short_sample_name=`echo "${sample_name}" | cut -d "." -f1` -v main_chrom=${main_chrom} -F "\t" 'BEGIN {OFS="\t"} {print $0, main_chrom, short_sample_name}' > ${TMPDIR}/short_sorted_table.bed

    # Get cluster range of the bam1 chromosome.
    while read line_in ; do
        read chromosome_bam2 Start_Pos End_Pos Width2 Nearest_Gene Post_filter_Reads chromosome_bam1 Sample <<< "${line_in}"
    
        echo -e "${chromosome_bam2}\t${Start_Pos}\t${End_Pos}" > "${TMPDIR}/cluster_1" 
    
        awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' < ${bamintersectBED12} | \
        sort-bed - | bedops -e 1 - "${TMPDIR}/cluster_1" | cut -f8-9 > ${TMPDIR}/cols_8_9
    
        min=`awk 'NR==1{min = $1 + 0; next} {if ($1 < min) min = $1;} END{print min}' ${TMPDIR}/cols_8_9`
        max=`awk 'NR==1{max = $2 + 0; next} {if ($2 > max) max = $2;} END{print max}' ${TMPDIR}/cols_8_9`
        let "Width1 = ${max} - ${min}"
    
        echo -e "${chromosome_bam2}\t${Start_Pos}\t${End_Pos}\t${Width2}\t${Nearest_Gene}\t${Post_filter_Reads}\t${chromosome_bam1}\t${min}\t${max}\t${Width1}\t${Sample}" >> ${OUTBASE}.counts.txt
    
    done < "${TMPDIR}/short_sorted_table.bed"
done < "${TMPDIR}/bam1_chroms"


#####################################################################################

echo "Finished counts_table.sh"
date

#####################################################################################
#####################################################################################
# For testing
# cp -r "${TMPDIR}" "${sampleOutdir}"

