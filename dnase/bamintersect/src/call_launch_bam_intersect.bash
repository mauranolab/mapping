#!/bin/bash
set -eu -o pipefail

src="/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/src"

# Maybe this should be passed in?
outdirs="/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/outdirs"

module load samtools/1.9
module load picard/2.18.15
module load bedops/2.4.35

################################################
# Start of "getopt" section
#
# No colons = option does not have an argument.
# One colon = option has a required argument.
# Two colons = option has an optional argument.
# Note: Use of the colons is not the same as making the option itself required/optional !
# Also, optional (two colon) arguments MUST be passed with an equal sign.
#    This works:            option=argument
#    This does not work:    option argument
# The "optional" feature is a GNU extension to the basic getopt.

long_arg_list=(
    sample_name:
    bamname1:
    bam1_keep_flags:
    bam1_exclude_flags:
    bamname2:
    bam2_keep_flags:
    bam2_exclude_flags:
    reads_match
    do_not_make_csv
    do_not_make_table
    help
)

long_args=$(printf "%s," "${long_arg_list[@]}")   # Turn the arg array into a string with commas.
long_args=$(echo ${long_args} | sed 's/,$/ /')    # Get rid of final comma.

# There must be at least one option associated with "-o", or you need to deal with getopt's default behavior.
CMD_LINE=$(getopt -o h --long "${long_args}" -n "$(basename "$0")" -- "$@")
# Catch getopt errors here if not using "set -e"

eval set -- "$CMD_LINE"

# getopt default values:
    make_csv=True
    make_table=True
    reads_match=False
    
    # Require read is paired.
    # Require mate is unmapped.
    bam1_keep_flags="9"
    
    # Exclude read unmapped.
    # Exclude PCR or optical duplicates.
    # Exclude supplementary aligments.
    bam1_exclude_flags="3076"
    
    # Require read is paired.
    bam2_keep_flags="1"
    
    # Exclude read unmapped.
    # Exclude PCR or optical duplicates.
    # Exclude supplementary aligments.
    bam2_exclude_flags="3076"


# Extract options and their arguments into variables.
while true ; do
    case "$1" in
        --sample_name)
            sample_name=$2 ; shift 2 ;;
        --bamname1)
            bamname1=$2 ; shift 2 ;;
        --bam1_keep_flags)
            bam1_keep_flags=$2 ; shift 2 ;;
        --bam1_exclude_flags)
            bam1_exclude_flags=$2 ; shift 2 ;;
        --bamname2)
            bamname2=$2 ; shift 2 ;;
        --bam2_keep_flags)
            bam2_keep_flags=$2 ; shift 2 ;;
        --bam2_exclude_flags)
            bam2_exclude_flags=$2 ; shift 2 ;;
        --reads_match)
            reads_match=True ; shift 1 ;;
        --do_not_make_csv)
            make_csv=False ; shift 1 ;;
        --do_not_make_table)
            make_table=False ; shift 1 ;;
        -h|--help) echo "Some helpful text." ; shift ; exit 0 ;;
        --) shift ; break ;;
        *) echo "getopt internal error!" ; exit 1 ;;
    esac
done

# End of getopt section.
################################################
# We now have:  the standard sample_name, bamname1, bamname2
# Derive some other required values from what we have.

    x=${bamname1%.bam}
    cegsgenome=${x##*.}

    x=${bamname2%.bam}
    annotationgenome=${x##*.}

    sample_name="${sample_name}.${cegsgenome}_vs_${annotationgenome}"
    echo "final sample_name is: ${sample_name}"

    # Modify this for use in the flowcell pipeline ?
    FINAL_OUTDIR="${outdirs}/${sample_name}"

    # Get the correct chrX non-HA zones (0-based bed style numbers, not UCSC positions):
    # The 5p and 3p outer boundaries of the HAs are as of Mar 29, 2019.
    # The assumed HPRT1 deletion zone is:  chrX    134459946   134501642
    # See: /vol/cegs/sequences/hg38/HPRT1/HPRT1_assembly.bed
    # Also: \\research-cifs.nyumc.org\Research\CEGS\Sequences\Landing Pads Summary Table.xlsx (8 JUL 2019)
    if echo "${cegsgenome}" | grep -q "LP058" ; then
        # HA lengths: 1032/1014 (LP058)
        bam1_5p_HA="chrX:0-134458914"
        bam1_3p_HA="chrX:134502656-156040895"
    elif echo "${cegsgenome}" | grep -q "LP087" ; then
        # HA lengths: 251/250 (LP061)
        bam1_5p_HA="chrX:0-134459695"
        bam1_3p_HA="chrX:134501892-156040895"
    elif echo "${cegsgenome}" | grep -q "LP123" ; then
        # HA lengths: 251/250 (LP061)
        bam1_5p_HA="chrX:0-134459695"
        bam1_3p_HA="chrX:134501892-156040895"
    elif echo "${cegsgenome}" | grep -q "LP062" ; then
        # HA lengths: 100/100 (LP062)
        bam1_5p_HA="chrX:0-134459846"
        bam1_3p_HA="chrX:134501742-156040895"
    else
        bam1_5p_HA="NA"
        bam1_3p_HA="NA"
        echo "No HA match for ${sample_name}     genome is: ${cegsgenome}"
    fi


    # A bed file of ranges to delete at the end of the process.
    # exclude_regions_from_counts="/vol/mauranolab/cadlej01/projects/LP_Integration/LP_uninformative_regions/LP_vs_hg38.bed"
    #          or
    # exclude_regions_from_counts="NA"
    if echo "${cegsgenome}" | grep -q "LP058" ; then
        exclude_regions_from_counts="/vol/mauranolab/cadlej01/projects/LP_Integration/LP_uninformative_regions/LP_vs_hg38.bed"
        cegsLP=True
    elif echo "${cegsgenome}" | grep -q "LP087" ; then
        exclude_regions_from_counts="/vol/mauranolab/cadlej01/projects/LP_Integration/LP_uninformative_regions/LP_vs_hg38.bed"
        cegsLP=True
    elif echo "${cegsgenome}" | grep -q "LP123" ; then
        exclude_regions_from_counts="/vol/mauranolab/cadlej01/projects/LP_Integration/LP_uninformative_regions/LP_vs_hg38.bed"
        cegsLP=True
    elif echo "${cegsgenome}" | grep -q "LP062" ; then
        exclude_regions_from_counts="/vol/mauranolab/cadlej01/projects/LP_Integration/LP_uninformative_regions/LP_vs_hg38.bed"
        cegsLP=True
    elif echo "${cegsgenome}" | grep -q "pSpCas9" ; then
        exclude_regions_from_counts="/vol/mauranolab/cadlej01/projects/LP_Integration/LP_uninformative_regions/pSpCas9_vs_hg38.bed"
        cegsLP=False
    elif echo "${cegsgenome}" | grep -q "PL1" ; then
        exclude_regions_from_counts="/vol/mauranolab/cadlej01/projects/LP_Integration/LP_uninformative_regions/PL1_payload_vs_hg38.bed"
        cegsLP=False
    else
        echo "No uninformative regions match for ${sample_name}"
        echo "genome is: ${cegsgenome}"
        exclude_regions_from_counts="NA"
        cegsLP=False
    fi


    # Some special cases which modify the above:
    if [ "${cegsLP}" = "True" ]; then
        if echo "${annotationgenome}" | grep -q "pSpCas9" ; then
            exclude_regions_from_counts="/vol/mauranolab/cadlej01/projects/LP_Integration/LP_uninformative_regions/LP_vs_pSpCas9.bed"
        elif echo "${annotationgenome}" | grep -q "PL1" ; then
            exclude_regions_from_counts="/vol/mauranolab/cadlej01/projects/LP_Integration/LP_uninformative_regions/LP_vs_PL1_payload.bed"
        fi
    fi

################################################
# Make some output directories.

# Make the master output directory for this analysis
rm -rf ${FINAL_OUTDIR}  # Clear out the old one, if necessary.
mkdir ${FINAL_OUTDIR}

# Make a place for the related slurm reports.
mkdir "${FINAL_OUTDIR}/log"

# Make a place for all the small chromosome bam files.
mkdir "${FINAL_OUTDIR}/bams"

########################################################
# Make the chromosome lists which will eventually drive the array jobs.

samtools idxstats ${bamname1} | awk '$1 != "*" {print $1}' > "${FINAL_OUTDIR}/log/chrom_list1"
num_lines=$(wc -l < "${FINAL_OUTDIR}/log/chrom_list1")
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
    chrom_list1_input_to_samtools="NA"   # To avoid a pipefail.
fi


samtools idxstats ${bamname2} | awk '$1 != "*" {print $1}' > "${FINAL_OUTDIR}/log/chrom_list2"
num_lines=$(wc -l < "${FINAL_OUTDIR}/log/chrom_list2")
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
export_vars="${export_vars},input_to_samtools2=${chrom_list1_input_to_samtools}"

num_lines=$(wc -l < "${FINAL_OUTDIR}/log/chrom_list1_simple")

JOB_ID0=$(sbatch --parsable --export=ALL,${export_vars} --array="1-${num_lines}" \
        --output="${FINAL_OUTDIR}/log/sort_chrom_1.${sample_name}.o_%A_%a" "${src}/sort_chrom.sbatch")


export_vars="FINAL_OUTDIR=${FINAL_OUTDIR}"
export_vars="${export_vars},BAM=${bamname2}"
export_vars="${export_vars},BAM_K=${bam2_keep_flags}"      # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_E=${bam2_exclude_flags}"   # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_N=2"
export_vars="${export_vars},input_to_samtools2=${chrom_list2_input_to_samtools}"

num_lines=$(wc -l < "${FINAL_OUTDIR}/log/chrom_list2_simple")

JOB_ID1=$(sbatch --parsable --export=ALL,${export_vars} --array="1-${num_lines}" \
        --output="${FINAL_OUTDIR}/log/sort_chrom_2.${sample_name}.o_%A_%a" "${src}/sort_chrom.sbatch")

JOB_ID0=":${JOB_ID0}:${JOB_ID1}"

################################################################################################
## Prepare for the big array job:

export_vars="FINAL_OUTDIR=${FINAL_OUTDIR}"

JOB_ID1=$(sbatch --parsable --dependency=afterok${JOB_ID0} --export=ALL,${export_vars} \
          --output="${FINAL_OUTDIR}/log/prep_array_job.${sample_name}.o_%j" "${src}/prep_big_array_job.sbatch")

################################################################################################
## Launch the big array job:

export_vars="FINAL_OUTDIR=${FINAL_OUTDIR}"
export_vars="${export_vars},src=${src}"
export_vars="${export_vars},reads_match=${reads_match}"
export_vars="${export_vars},make_csv=${make_csv}"

## "big_array_job.sbatch" does the actual launching of the array job. The array job itself does not return
## till it is complete. So JOB_ID2 will not be marked complete till the array job finishes, yet the code does
## not need to pause at this point. "big_array_job.sbatch" is launched, and we move right on to "cleanup.sbatch", below.
JOB_ID2=$(sbatch --parsable --dependency=afterok:${JOB_ID1} --export=ALL,${export_vars} \
          --output="${FINAL_OUTDIR}/log/big_array_job.${sample_name}.o_%j" "${src}/big_array_job.sbatch")

################################################################################################
## Cleanup:

export_vars="FINAL_OUTDIR=${FINAL_OUTDIR}"
export_vars="${export_vars},src=${src}"
## export_vars="${export_vars},reads_match=${reads_match}"
export_vars="${export_vars},exclude_regions_from_counts=${exclude_regions_from_counts}"
export_vars="${export_vars},bam1_5p_HA=${bam1_5p_HA}"
export_vars="${export_vars},bam1_3p_HA=${bam1_3p_HA}"
export_vars="${export_vars},sample_name=${sample_name}"
export_vars="${export_vars},cegsgenome=${cegsgenome}"
export_vars="${export_vars},annotationgenome=${annotationgenome}"
export_vars="${export_vars},make_csv=${make_csv}"
export_vars="${export_vars},make_table=${make_table}"

sbatch --parsable --dependency=afterok:${JOB_ID2} --export=ALL,${export_vars} \
       --output="${FINAL_OUTDIR}/log/cleanup.${sample_name}.o_%j" "${src}/cleanup.sbatch"

