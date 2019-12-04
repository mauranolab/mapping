#!/bin/bash
set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'

#####################################################################################

sampleOutdir=$1
sample_name=$2
bam2genome=$3
HAtableOutputBED12=$4
counts_type=$5
INTERMEDIATEDIR=$6
genome2exclude=$7

#####################################################################################
# Get reads that are in the deletion zone.

awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' < "${HAtableOutputBED12}" | sort-bed - | tee "${TMPDIR}/sorted_bam2onLeft.bed12" | \
bedops -e 1 - "${INTERMEDIATEDIR}/deletion_range.bed" > "${TMPDIR}/del_range_reads_out.bed"

filter_csv_output="${TMPDIR}/noDZreads.bed12"

bedops -n 1 "${TMPDIR}/sorted_bam2onLeft.bed12" "${INTERMEDIATEDIR}/deletion_range.bed" | \
awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' > ${filter_csv_output}

echo -e "${sample_name}:" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

if [ -s "${INTERMEDIATEDIR}/deletion_range.bed" ]; then
    num_lines=$(wc -l < "${TMPDIR}/del_range_reads_out.bed")
    echo -e "    Number of reads in the Deletion Range:\t${num_lines}" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
else
    echo -e "    There are no HAs, so there is no Deletion Range." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
fi
echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

#####################################################################################
# Maybe we want to filter out some reads.
if [ "${genome2exclude}" != "null" ];then
    # Comments for the 5 piped lines in the "else" part of the below if statement:
    # 1) Delete reads that overlap the ranges defined in genome2exclude (which are defined with respect to bam1 coordinates).
    #    These were attached to the genome1exclude.bed file via bedops -u in submit_bamintrsect.sh
    # 2) Now delete reads that are in the HAs.
    #     2a) The output of step 1 has 12 columns: [chr start end readID flag +/-] x 2, the bam1 data being in 1-6, and the bam2 data being in 7-12.
    #         The HA coordinates are with respect to bam2, so we need to switch the 6-column halves to use the bedops functions on this file.
    # 3) Then sort by the bam2 reads, then keep only the bam2 reads (and their bam1 mates) that do not overlap the genome1exclude's (HAs):
    # 4) Now put the surviving bam1 reads back on the left hand side,
    # 5) and sort.

    if [ "${genome2exclude}" = "HA_only" ]; then
        awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' "${filter_csv_output}" | \
        sort-bed - | bedops -n 1 - "${INTERMEDIATEDIR}/genome1exclude.bed" | \
        awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | \
        sort-bed - > "${filter_csv_output}_tmp"
        mv "${filter_csv_output}_tmp" "${filter_csv_output}"
    else
        sort-bed "${filter_csv_output}" | bedops -n 1 - "${INTERMEDIATEDIR}/genome1exclude.bed" | \
        awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | \
        sort-bed - | bedops -n 1 - "${INTERMEDIATEDIR}/genome1exclude.bed" | \
        awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | \
        sort-bed - > "${filter_csv_output}_tmp"
        mv "${filter_csv_output}_tmp" "${filter_csv_output}"
    fi

    cp "${filter_csv_output}" "${sampleOutdir}/${sample_name}.informative.bed"
fi
#####################################################################################

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

# Loop over the bam1 chromosome names.
cut -f1 "${filter_csv_output}" | sort | uniq > "${TMPDIR}/bam1_chroms"

# Clear this, which may exist due to previous calls to filter_tsv
rm -f ${TMPDIR}/short_sorted_table2.bed

while read main_chrom ; do
    # Get just the hg38/mm10/rn6 bed3 part, then sort it:
    awk -v mainChrom="${main_chrom}" -F "\t" 'BEGIN {OFS="\t"} $1==mainChrom {print $0; }' "${filter_csv_output}" |
    cut -f7-9 | sort-bed - > ${TMPDIR}/sorted_speciesReadCoords.bed3
    
    # The below makes +/- 500 bps regions around each hg38/mm10 read, and merges them when possible. It will return bed3 output.
    bedops --range 500 -m ${TMPDIR}/sorted_speciesReadCoords.bed3 > ${TMPDIR}/all_regions.bed
    
    cat ${geneAnnotationFile} |
    awk -F "\t" 'BEGIN {OFS="\t"} $7=="protein_coding"' |
    #Restrict to level 1 genes if available
    awk -F "\t" 'BEGIN {OFS="\t"} $5==1 || $5=="."' |
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

    # Get cluster range of the bam1 chromosome.
    while read line_in ; do
        read chromosome_bam2 Start_Pos End_Pos Width2 Nearest_Gene Post_filter_Reads chromosome_bam1 Sample <<< "${line_in}"
     
        echo -e "${chromosome_bam2}\t${Start_Pos}\t${End_Pos}" > "${TMPDIR}/cluster_1" 
    
        awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' < "${filter_csv_output}" | \
        sort-bed - | bedops -e 1 - "${TMPDIR}/cluster_1" | cut -f8-9 > ${TMPDIR}/cols_8_9
    
        min=`awk 'NR==1{min = $1 + 0; next} {if ($1 < min) min = $1;} END{print min}' ${TMPDIR}/cols_8_9`
        max=`awk 'NR==1{max = $2 + 0; next} {if ($2 > max) max = $2;} END{print max}' ${TMPDIR}/cols_8_9`
        let "Width1 = ${max} - ${min}"
    
        echo -e "${chromosome_bam2}\t${Start_Pos}\t${End_Pos}\t${Width2}\t${Nearest_Gene}\t${Post_filter_Reads}\t${chromosome_bam1}\t${min}\t${max}\t${Width1}\t${Sample}" >> ${TMPDIR}/short_sorted_table2.bed
    
    done < "${TMPDIR}/short_sorted_table.bed"
done < "${TMPDIR}/bam1_chroms"


# Output the final counts files.
if [ -f "${TMPDIR}/short_sorted_table2.bed" ]; then
    if [ ${counts_type} = "all_reads_counts" ]; then
        echo -e "chromosome-bam2\tStart_Pos\tEnd_Pos\tWidth\tNearest_Gene\tPost-filter_Reads\tchromosome-bam1\tStart_Pos\tEnd_Pos\tWidth\tSample" > "${sampleOutdir}/${sample_name}.counts.txt"
        cat ${TMPDIR}/short_sorted_table2.bed >> ${sampleOutdir}/${sample_name}.counts.txt
    else
        echo -e "chromosome-bam2\tStart_Pos\tEnd_Pos\tWidth\tNearest_Gene\tPost-filter_Reads\tchromosome-bam1\tStart_Pos\tEnd_Pos\tWidth\tSample" > "${sampleOutdir}/${sample_name}.informative.counts.txt"
        cat ${TMPDIR}/short_sorted_table2.bed >> ${sampleOutdir}/${sample_name}.informative.counts.txt
    fi
else
    # Generate blank files when there are no reads.
    if [ ${counts_type} = "all_reads_counts" ]; then
        echo -e "chromosome-bam2\tStart_Pos\tEnd_Pos\tWidth\tNearest_Gene\tPost-filter_Reads\tchromosome-bam1\tStart_Pos\tEnd_Pos\tWidth\tSample" > "${sampleOutdir}/${sample_name}.counts.txt"
    else
        echo -e "chromosome-bam2\tStart_Pos\tEnd_Pos\tWidth\tNearest_Gene\tPost-filter_Reads\tchromosome-bam1\tStart_Pos\tEnd_Pos\tWidth\tSample" > "${sampleOutdir}/${sample_name}.informative.counts.txt"
    fi
fi

#####################################################################################

echo "Finished filter_tsv.sh"
date

#####################################################################################
#####################################################################################
# For testing
# cp -r "${TMPDIR}" "${sampleOutdir}"

