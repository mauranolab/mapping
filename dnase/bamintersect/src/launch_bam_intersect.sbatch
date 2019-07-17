#!/bin/bash
set -eu -o pipefail

#SBATCH --job-name=launch_bam_intersect

########################################################
date
module load samtools/1.9
module load picard/2.18.15
module load bedops/2.4.35
## https://bedops.readthedocs.io/en/latest/content/reference/set-operations/bedops.html
########################################################
## We are testing the status of SOURCE_FAIL manually, and don't need to enable "set -u" here.
set +u
size=${#SOURCE_FILE}
set -u   ## Turn it back on now.
if [ ${size} -lt 1 ]; then
    echo No source file provided via sbatch export.
    sample_name=sample
    SOURCE_FILE=$1
    FINAL_OUTDIR=$2
    src=$3
    ## These should be in the form of:
    ##     SOURCE_FILE=<some full path>/<the name of an existing file, containing input values>
    ##     FINAL_OUTDIR=<some full path>/<a directory name, which will be created via mkdir below>
else
    ## $1 is empty.
    echo We have a source file via sbatch export.
    sample_name=${SOURCE_FILE##*/}
    ## Passed in via sbatch export:
    ##    SOURCE_FILE
    ##    FINAL_OUTDIR
    ##    src
fi
########################################################
source $SOURCE_FILE

x=${bamname1%.bam}    # Cut off the trailing "bam"
y=${x##*.}            # Assign the base to y. This will be one of genotypes, like "hg38_full"
z=${x%.${y}}          # Cut off the trailing ${y}
sample_name=${z##*/}  # The remaining base is the sample name.
echo "sample_name is: ${sample_name}"

x=${bamname2%.bam}
annotationgenome=${x##*.}
echo "annotationgenome is: ${annotationgenome}"

## SOURCE_FILE is normally stored in ${src}/input_library
## It sets these variables (the example assignments are just for illustration):

## bamname1
## bam1_keep_flags
## bam1_exclude_flags
## bamname2
## bam2_keep_flags
## bam2_exclude_flags
## final_csv_delete

## Examples:
## H1_dPIGA_dHPRT_PL1_Capture-BS02660A_LP087_noEF1a

## cegsvectors bam file:
## bamname1=/vol/cegs/mapped/FCH7NKKBGXB/capture/H1_dPIGA_dHPRT1_PL1_Capture-BS02660A/H1_dPIGA_dHPRT1_PL1_Capture-BS02660A.cegsvectors_LP087.bam

## Require read is paired.
## Require mate is unmapped.
## bam1_keep_flags="9"

## Exclude read unmapped.
## Exclude PCR or optical duplicates.
## Exclude supplementary aligments.
## bam1_exclude_flags="3076"

## hg38 bam file:
## bamname2=/vol/cegs/mapped/FCH7NKKBGXB/capture/H1_dPIGA_dHPRT1_PL1_Capture-BS02660A/H1_dPIGA_dHPRT1_PL1_Capture-BS02660A.hg38_full.bam

## Require read is paired.
## bam2_keep_flags="1"

## Exclude read unmapped.
## Exclude PCR or optical duplicates.
## Exclude supplementary aligments.
## bam2_exclude_flags="3076"

## This parameter tells bam_intersect.py how to match the reads from each bam file.
## It can be:  'True', 'False', or 'None'
## Should the matches be read1/read1 (True), or read1/read2 (False)?  None designates unpaired reads.
## reads_match=False

## add "noEF1a" kinds of ranges here. Numbers are BED style, not UCSC positions.
## LP087   60  1239    EF-1-alpha promoter
## LP087   6130    6536    hEIF1pA
## final_csv_delete="LP087:60-1239|LP087:6130-6536"
##          or
## final_csv_delete=""

################################################################################################
## Get the chromosomes for this sample.

samtools idxstats ${bamname1} | awk '$1 != "*" {print $1}' > "${FINAL_OUTDIR}/log/chrom_list1"
num_lines=$(wc -l < "${FINAL_OUTDIR}/log/chrom_list1")
echo "chrom_list1 num_lines is: ${num_lines}"
cp "${FINAL_OUTDIR}/log/chrom_list1" "${FINAL_OUTDIR}/inputs.chrom_list1.txt"

if [ "${num_lines}" -ge "5" ]; then
    ## Get the simple chromosome name mask
    grep -E '^chr[0-9]*$|^chr[XYM]$' "${FINAL_OUTDIR}/log/chrom_list1" | sed 's/^/^/' | sed 's/$/$/' > "${FINAL_OUTDIR}/log/chrom_list1_simple_mask"
    
    ## Apply mask to chomosome names
    grep -v -f "${FINAL_OUTDIR}/log/chrom_list1_simple_mask" "${FINAL_OUTDIR}/log/chrom_list1" > "${FINAL_OUTDIR}/log/chrom_list1_long"
    
    ## We'll strip the pipes out in sort_chrom.sbatch
    chrom_list1_input_to_samtools=$(cat "${FINAL_OUTDIR}/log/chrom_list1_long" | paste -sd "|" -)
    
    ## Get simple chromosome names
    grep -E '^chr[0-9]*$|^chr[XYM]$' "${FINAL_OUTDIR}/log/chrom_list1" >  "${FINAL_OUTDIR}/log/chrom_list1_simple"
    num_lines=$(wc -l < "${FINAL_OUTDIR}/log/chrom_list1_long")
    if [ "${num_lines}" -ge "1" ]; then
        echo "all_other" >> "${FINAL_OUTDIR}/log/chrom_list1_simple"
    fi
else
    # A small number of chromosomes is an indicator that this is a LP/PL file.
    # There will be no "all_other" in this file.
    cp "${FINAL_OUTDIR}/log/chrom_list1" "${FINAL_OUTDIR}/log/chrom_list1_simple"
fi


samtools idxstats ${bamname2} | awk '$1 != "*" {print $1}' > "${FINAL_OUTDIR}/log/chrom_list2"
num_lines=$(wc -l < "${FINAL_OUTDIR}/log/chrom_list2")
echo "chrom_list2 num_lines is: ${num_lines}"
cp "${FINAL_OUTDIR}/log/chrom_list2" "${FINAL_OUTDIR}/inputs.chrom_list2.txt"

## Get the simple chromosome name mask
grep -E '^chr[0-9]*$|^chr[XYM]$' "${FINAL_OUTDIR}/log/chrom_list2" | sed 's/^/^/' | sed 's/$/$/' > "${FINAL_OUTDIR}/log/chrom_list2_simple_mask"

## Apply mask to chomosome names
grep -v -f "${FINAL_OUTDIR}/log/chrom_list2_simple_mask" "${FINAL_OUTDIR}/log/chrom_list2" > "${FINAL_OUTDIR}/log/chrom_list2_long"

## We'll strip the pipes out in sort_chrom.sbatch
chrom_list2_input_to_samtools=$(cat "${FINAL_OUTDIR}/log/chrom_list2_long" | paste -sd "|" -)

## Get simple chromosome names
grep -E '^chr[0-9]*$|^chr[XYM]$' "${FINAL_OUTDIR}/log/chrom_list2" >  "${FINAL_OUTDIR}/log/chrom_list2_simple"
num_lines=$(wc -l < "${FINAL_OUTDIR}/log/chrom_list2_long")
if [ "${num_lines}" -ge "1" ]; then
    echo "all_other" >> "${FINAL_OUTDIR}/log/chrom_list2_simple"
fi
################################################################################################
## Create the small chromosome bam files, and sort them:

export_vars="FINAL_OUTDIR=${FINAL_OUTDIR}"
export_vars="${export_vars},BAM=${bamname1}"
export_vars="${export_vars},BAM_K=${bam1_keep_flags}"      # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_E=${bam1_exclude_flags}"   # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_N=1"

JOB_ID0=""

while read -r chrom; do
    export_vars2="${export_vars},chrom=${chrom}"

    if [ ${chrom} = "all_other" ]; then
        export_vars2="${export_vars2},input_to_samtools=${chrom_list1_input_to_samtools}"
    else
        export_vars2="${export_vars2},input_to_samtools=${chrom}"
    fi

    JOB_ID1=$(sbatch --parsable --export=ALL,${export_vars2} \
            --output="${FINAL_OUTDIR}/log/sort_chrom.${sample_name}.${chrom}.1.o_%j" "${src}/sort_chrom.sbatch")

    JOB_ID0="${JOB_ID0}:${JOB_ID1}"
done < "${FINAL_OUTDIR}/log/chrom_list1_simple"; 

export_vars="FINAL_OUTDIR=${FINAL_OUTDIR}"
export_vars="${export_vars},BAM=${bamname2}"
export_vars="${export_vars},BAM_K=${bam2_keep_flags}"      # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_E=${bam2_exclude_flags}"   # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_N=2"


while read -r chrom; do
    export_vars2="${export_vars},chrom=${chrom}"

    if [ ${chrom} = "all_other" ]; then
        export_vars2="${export_vars2},input_to_samtools=${chrom_list2_input_to_samtools}"
    else
        export_vars2="${export_vars2},input_to_samtools=${chrom}"
    fi

    JOB_ID1=$(sbatch --parsable --export=ALL,${export_vars2} \
            --output="${FINAL_OUTDIR}/log/sort_chrom.${sample_name}.${chrom}.2.o_%j" "${src}/sort_chrom.sbatch")

    JOB_ID0="${JOB_ID0}:${JOB_ID1}"
done < "${FINAL_OUTDIR}/log/chrom_list2_simple"; 


export_vars="FINAL_OUTDIR=${FINAL_OUTDIR}"
export_vars="${export_vars},src=${src}"
export_vars="${export_vars},reads_match=${reads_match}"

## Need to send these into cleanup when we enable filter_tsv.bash.
## They also need to be added to SOURCE_FILE.
## export_vars="${export_vars},final_csv_delete=${final_csv_delete}"
## export_vars="${export_vars},bam1_5p_HA=${bam1_5p_HA}"
## export_vars="${export_vars},bam1_3p_HA=${bam1_3p_HA}"
## export_vars="${export_vars},sample_name=${sample_name}"
## export_vars="${export_vars},annotationgenome=${annotationgenome}"
## export_vars="${export_vars},make_csv=${make_csv}"

sbatch --parsable --dependency=afterok${JOB_ID0} --export=ALL,${export_vars} \
       --output="${FINAL_OUTDIR}/log/cleanup.${sample_name}.o_%j" "${src}/cleanup.sbatch"

date

