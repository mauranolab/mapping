#!/bin/bash
set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'

module load miller/5.4.0

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
# sampleName is in the form: "<sample_name>.<bam1genome>_vs_<bam2genome>"
#
#####################################################################################

OUTBASE=$1
sample_name=$2
bam1genome=$3
bam2genome=$4
bamintersectBED12=$5

echo "[counts_table] Making counts_table from ${bamintersectBED12}; will output to ${OUTBASE}.counts.txt"


echo -e -n "Number of total reads:\t"
cat ${bamintersectBED12} | wc -l


if [ ! -s ${bamintersectBED12} ]; then
    echo "[counts_table] There are no reads to analyze, quitting successfully."
    touch ${OUTBASE}.counts.txt
    exit 0
fi

#Annotate a file by adding a column with the nearest gene name. Inputfile can be stdin (as -); outputs to stdout
annotateNearestGeneName() {
    local inputfile=$1
    local genome=$2
    
    # Set the appropriate gene annotation file.
    # Note that hg38_full, etc. was converted to just "hg38" prior to calling bamintersect.
    local geneAnnotationFile=""
    case "${genome}" in
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
        geneAnnotationFile="CEGS"
        ;;
    esac
    
    if [[ "${geneAnnotationFile}" == "CEGS" ]]; then
        cat /vol/cegs/sequences/cegsvectors/cegsvectors/cegsvectors.bed | perl -pe 's/ /_/g;' > $TMPDIR/geneAnnotationFile.${genome}.bed
    else
        #Needs prefiltering
        #Add gene name
        cat ${geneAnnotationFile} |
        awk -F "\t" 'BEGIN {OFS="\t"} $7=="protein_coding"' |
        #Restrict to level 1 or 2 genes if available (Sox2 is level 2 for example)
        awk -F "\t" 'BEGIN {OFS="\t"} $5<=2 || $5=="."' > $TMPDIR/geneAnnotationFile.${genome}.bed
    fi
    geneAnnotationFile="$TMPDIR/geneAnnotationFile.${genome}.bed"
    
    cat ${inputfile} |
    #Replace main bed coordinates with the center of each element to refine annotation ID
    awk -F "\t" 'BEGIN {OFS="\t"} {width=$3-$2; center=int($2+width/2); print $1, center-1, center, $0}' | sort-bed - |
    closest-features --delim "|" --closest --dist - ${geneAnnotationFile} |
    #Replace empty gene names with a dash
    #Add distance to the gene name; note distance is relative to the genomic location (+ is to the right) not the transcriptional direction
    awk -F "|" 'BEGIN {OFS="\t"} { \
        split($2, gene, "\t"); if(gene[4]=="") {gene[4]="-"} \
        if($3=="NA" || $3==0) { \
            distance = ""; \
        } else if($3>0) { \
            distance = "-" $3 "bp"; \
        } else { \
            #Expression needs to be surrounded by parentheses to avoid getting interpreted as subtraction (thus swallowing the prior token)
            distance = "+" (-1*$3) "bp"; \
        } \
        print $1, gene[4] distance; \
    }' |
    #Drop the centered coordinates in the first three columns
    cut -f4- | sort-bed -
}


#Takes a bed file and outputs a list of minimal regions covering all elements grouped by fixed distance to stdout
merge_tight() {
    local file=$1
    local range=$2
    
    #Generate bam1 regions
    cat ${file} | bedops --range ${range} -m - | 
    bedmap --echo-map-range - ${file} |
    #Put the region coordinates into the ID column
    bedmap --delim "\t" --echo --echo-ref-name -
}

bamintersectBED12_byStrand="${TMPDIR}/${sample_name}.bamintersectBED12_byStrand.bed"
for strand_bam1 in "+" "-"; do
    for strand_bam2 in "+" "-"; do
        # Get reads from just one strand.
        awk -v strand1=${strand_bam1} -v strand2=${strand_bam2} 'BEGIN {FS="\t";OFS="\t"} {if($6==strand1 && $12==strand2){print $0}}' ${bamintersectBED12} > ${bamintersectBED12_byStrand}
        if [ ! -s ${bamintersectBED12_byStrand} ]; then
            # No reads for this strand combination are present.
            continue
        fi
        
        #Generate bam1 regions
        merge_tight ${bamintersectBED12_byStrand} 500 > ${TMPDIR}/${sample_name}.regions.bam1.bed
        
        #Split output into separate files by bam1 region
        mkdir $TMPDIR/bam2readsByBam1region
        bedmap --bp-ovr 1 --delim "|" --multidelim "|" --echo --echo-map ${TMPDIR}/${sample_name}.regions.bam1.bed ${bamintersectBED12_byStrand} |
        #Output format: bam2 chrom chromStart, chromEnd, bam1region
        awk -v outdir=$TMPDIR/bam2readsByBam1region -F "|" 'BEGIN {OFS="\t"} {split($1, bam1region, "\t"); for(i=2; i<=NF; i++) {split($i, bam2read, "\t"); print bam2read[7], bam2read[8], bam2read[9], bam1region[4] > outdir "/" NR ".unsorted.bam2.bed"}}'
        
        for file in $TMPDIR/bam2readsByBam1region/*.unsorted.bam2.bed; do
            cur=`basename ${file} .unsorted.bam2.bed`
            cat ${file} | sort-bed - > $TMPDIR/bam2readsByBam1region/${cur}.sorted.bam2.bed
            
            #Generate bam2 regions for the set of reads from each bam1 region
            merge_tight $TMPDIR/bam2readsByBam1region/${cur}.sorted.bam2.bed 500 | cut -f1-3 |
            #Add a bam1 region id for each read that maps to this bam2 region
            bedmap --delim "\t" --multidelim ";" --echo --echo-map-id - $TMPDIR/bam2readsByBam1region/${cur}.sorted.bam2.bed |
            awk -F "\t" 'BEGIN {OFS="\t"} {split($4, bam1region, ";"); for(i=1; i<=length(bam1region); i++) {print $1, $2, $3, bam1region[i]}}'
        done |
        #Each line is the pair of bam1/bam2 regions spanned by a single mate pair, so now count them
        sort -k1,1 -k2,2 -k4,4 | uniq -c |
        #Now filter for minimum reads and finalize output columns
        awk -v minReadsCutoff=2 'BEGIN {OFS="\t"} $1>=minReadsCutoff {split($5, bam1region, /[:-]/); print $2, $3, $4, $4-$3, $1, bam1region[1], bam1region[2], bam1region[3], bam1region[3]-bam1region[2]}' | sort-bed - |
        #Add gene name to bam2
        annotateNearestGeneName - ${bam2genome} |
        #switch to bam1 coordinates, leaving the bam2 gene name at the end where it belongs
        awk -F "\t" 'BEGIN {OFS="\t"} {print $6, $7, $8, $9, $5, $1, $2, $3, $4, $10}' | sort-bed - |
        #Add gene name to bam1
        annotateNearestGeneName - ${bam1genome} |
        #Sort by Reads, chrom_bam2, chromStart_bam2
        #tried mlr but it doesn't like the field name below
        #mlr --tsv sort -nr Reads -f ${#chrom_bam2} -n chromStart_bam2
        sort -k5,5nr -k1,1 -k2,2n |
        #Reorder columns to put bam1 gene where it belongs
        awk -v short_sample_name=`echo "${sample_name}" | cut -d "." -f1` -v strand1=${strand_bam1} -v strand2=${strand_bam2} -F "\t" 'BEGIN {OFS="\t"; print "#chrom_bam1", "chromStart_bam1", "chromEnd_bam1", "Width_bam1", "NearestGene_bam1", "Strand_bam1", "Reads", "chrom_bam2", "chromStart_bam2", "chromEnd_bam2", "Width_bam2", "NearestGene_bam2", "Strand_bam2", "Sample"} {print $1, $2, $3, $4, $11, strand1, $5, $6, $7, $8, $9, $10, strand2, short_sample_name}' >> ${OUTBASE}.counts.tmp.txt
        
        #Cleanup
        rm -rf $TMPDIR/bam2readsByBam1region
    done
done

# Delete extra column header lines, and sort.
awk 'BEGIN {FS="\t";OFS="\t"} {if(NR==1 || $1!="#chrom_bam1"){print $0}else{next}}' ${OUTBASE}.counts.tmp.txt |
mlr --tsv sort -f \#chrom_bam1 -n chromStart_bam1 > ${OUTBASE}.counts.txt
rm ${OUTBASE}.counts.tmp.txt

#####################################################################################

echo "[counts_table] Done"
date
