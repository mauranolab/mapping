#!/bin/bash
set -eu -o pipefail

src="/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6"

# Maybe this should be passed in?
outdirs="/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/outdirs"

module load samtools/1.9

module load picard/2.18.15 1>/dev/null
# Redirects this output from picard:
#     Loading picard/2.18.15
#       Loading requirement: jdk/jdk1.8.0_101

module load bedops/2.4.35

############################################################
# Start of "getopt" section
############################################################
# For use with the getopt "-h" flag:

function usage {
cat << USAGE

Usage: $(basename "$0") [Options]
  Required Options:
     --sample_name {sample name}
     --bam1 {full path to bam file #1}
     --bam1genome {something like LP123, mm10, or hg38. Others are valid as well.}
     --bam2 {full path to bam file #2}
     --bam2genome {something like LP123, mm10, or hg38. Others are valid as well.}

  Optional Options:
          --bam1_keep_flags {a sam format flag value required of all bam1 reads to be 
                             included in the analysis. Default = 9}
                            
          --bam1_exclude_flags {a sam format flag value required of all bam1 reads to be 
                                excluded from the analysis. Default = 3076}
                  Examples:

                      # Require read is paired.
                      # Require mate is unmapped.
                      bam1_keep_flags="9"

                      # Exclude read unmapped.
                      # Exclude PCR or optical duplicates.
                      # Exclude supplementary aligments.
                      bam1_exclude_flags="3076"

          --bam2_keep_flags {a sam format flag value required of all bam2 reads to be 
                             included in the analysis. Default = 1}
                            
          --bam2_exclude_flags {a sam format flag value required of all bam2 reads to be 
                                excluded from the analysis. Default = 3076}

          --reads_match {no argument}
                 # This flag tells bam_intersect.py how to match the reads from each bam file.
                 # If set, reads should be matched like [read1 read1], or [read2 read2].
                 # If unset, reads are paired like [read1 read2], or [read2 read1].
                 # Set this option for unpaired reads.
                 # The default is to have it unset.

           --do_not_make_csv {no argument}
                 # Main output is a single bed file (unset), or 2 bam files (set).
                 # The default is to have it unset.

           --do_not_make_table {no argument}
                 # Output a "counts.txt" table summarizing results (unset), or no table (set).
                 # The default is to have it unset.

           --clear_logs {no argument}
                 # If set, clears out the log file subdirectories, since they can accumulate logs from previous runs.
                 # The default is to have it unset.

           -h,--help    Print this help message.

USAGE
}
############################################################
# Defining the getopt arguments:
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
    bam1:
    bam1genome:
    bam1_keep_flags:
    bam1_exclude_flags:
    bam2:
    bam2genome:
    bam2_keep_flags:
    bam2_exclude_flags:
    reads_match
    do_not_make_csv
    do_not_make_table
    clear_logs
    help
)

############################################################
# Parsing the getopt arguments:

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
    clear_logs=False
    
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
        --bam1)
            bamname1=$2 ; shift 2 ;;
        --bam1genome)
            cegsgenome=$2 ; shift 2 ;;
        --bam1_keep_flags)
            bam1_keep_flags=$2 ; shift 2 ;;
        --bam1_exclude_flags)
            bam1_exclude_flags=$2 ; shift 2 ;;
        --bam2)
            bamname2=$2 ; shift 2 ;;
        --bam2genome)
            annotationgenome=$2 ; shift 2 ;;
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
        --clear_logs)
            clear_logs=True ; shift 1 ;;
        -h|--help) usage; shift ; exit 0 ;;
        --) shift ; break ;;
        *) echo "getopt internal error!" ; exit 1 ;;
    esac
done

# End of getopt section.
################################################
# We now have:  the standard sample_name, bamname1, bamname2, cegsgenome, annotationgenome
# Derive some other required values from what we have.

    sample_name="${sample_name}.${cegsgenome}_vs_${annotationgenome}"
    echo "final sample_name is: ${sample_name}"

    # Modify this for use in the flowcell pipeline ?
    sampleOutdir="${outdirs}/${sample_name}"

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

if [ "${clear_logs}" = "True" ]; then
    rm -rf "${sampleOutdir}/log"  # Clear out the log directory.
fi

if [ -d "${sampleOutdir}" ]; then
    if [ ! -d "${sampleOutdir}/log" ]; then
        mkdir "${sampleOutdir}/log"
    fi

    if [ ! -d "${sampleOutdir}/bams" ]; then
        mkdir "${sampleOutdir}/bams"
    fi
else
    # Make the master output directory for this analysis
    mkdir ${sampleOutdir}

    # Make a place for the related slurm reports.
    mkdir "${sampleOutdir}/log"

    # Make a place for all the small chromosome bam files.
    mkdir "${sampleOutdir}/bams"
fi

########################################################
# Make the chromosome lists which will eventually drive the array jobs.

samtools idxstats ${bamname1} | awk '$1 != "*" {print $1}' > "${sampleOutdir}/log/chrom_list1"
num_lines=$(wc -l < "${sampleOutdir}/log/chrom_list1")
cp "${sampleOutdir}/log/chrom_list1" "${sampleOutdir}/inputs.chrom_list1.txt"

if [ "${num_lines}" -ge "5" ]; then
    ## Get the simple chromosome name mask
    grep -E '^chr[0-9]*$|^chr[XYM]$' "${sampleOutdir}/log/chrom_list1" | sed 's/^/^/' | sed 's/$/$/' > "${sampleOutdir}/log/chrom_list1_simple_mask"

    ## Apply mask to chomosome names
    grep -v -f "${sampleOutdir}/log/chrom_list1_simple_mask" "${sampleOutdir}/log/chrom_list1" > "${sampleOutdir}/log/chrom_list1_long"

    ## We'll strip the pipes out in sort_bamintersect.sh
    chrom_list1_input_to_samtools=$(cat "${sampleOutdir}/log/chrom_list1_long" | paste -sd "|" -)

    ## Get simple chromosome names
    grep -E '^chr[0-9]*$|^chr[XYM]$' "${sampleOutdir}/log/chrom_list1" >  "${sampleOutdir}/log/chrom_list1_simple"
    num_lines=$(wc -l < "${sampleOutdir}/log/chrom_list1_long")
    if [ "${num_lines}" -ge "1" ]; then
        echo "all_other" >> "${sampleOutdir}/log/chrom_list1_simple"
    fi
else
    # A small number of chromosomes is an indicator that this is a LP/PL file.
    # There will be no "all_other" in this file.
    cp "${sampleOutdir}/log/chrom_list1" "${sampleOutdir}/log/chrom_list1_simple"
    chrom_list1_input_to_samtools="NA"   # To avoid a pipefail.
fi


samtools idxstats ${bamname2} | awk '$1 != "*" {print $1}' > "${sampleOutdir}/log/chrom_list2"
num_lines=$(wc -l < "${sampleOutdir}/log/chrom_list2")
cp "${sampleOutdir}/log/chrom_list2" "${sampleOutdir}/inputs.chrom_list2.txt"

## Get the simple chromosome name mask
grep -E '^chr[0-9]*$|^chr[XYM]$' "${sampleOutdir}/log/chrom_list2" | sed 's/^/^/' | sed 's/$/$/' > "${sampleOutdir}/log/chrom_list2_simple_mask"

## Apply mask to chomosome names
grep -v -f "${sampleOutdir}/log/chrom_list2_simple_mask" "${sampleOutdir}/log/chrom_list2" > "${sampleOutdir}/log/chrom_list2_long"

## We'll strip the pipes out in sort_bamintersect.sh
chrom_list2_input_to_samtools=$(cat "${sampleOutdir}/log/chrom_list2_long" | paste -sd "|" -)

## Get simple chromosome names
grep -E '^chr[0-9]*$|^chr[XYM]$' "${sampleOutdir}/log/chrom_list2" >  "${sampleOutdir}/log/chrom_list2_simple"
num_lines=$(wc -l < "${sampleOutdir}/log/chrom_list2_long")
if [ "${num_lines}" -ge "1" ]; then
    echo "all_other" >> "${sampleOutdir}/log/chrom_list2_simple"
fi
################################################################################################
# Setup directories and files for the array jobs.

for BAM_N in $(seq 1 2); do
    if [ "${BAM_N}" = "1" ]; then
        BAM=${bamname1}
    else
        BAM=${bamname2}
    fi

    while read chrom;  do
        x=${BAM%.bam}         # Cut off the trailing "bam"
        y=${x##*.}            # Assign the base to y. This will be one of genotypes, like "hg38_full"
        z=${x%.${y}}          # Cut off the trailing ${y}
        short_name=${z##*/}   # The remaining base is the short name.

        BAM_OUT="${sampleOutdir}/bams/${short_name}.${chrom}.${BAM_N}.bam"
        echo ${BAM_OUT} >> "${sampleOutdir}/bams/chr_list_${BAM_N}"
    done < "${sampleOutdir}/log/chrom_list${BAM_N}_simple"
done

TEMP_DIR_CL1="${sampleOutdir}/TEMP_DIR_CL1"
rm -rf ${TEMP_DIR_CL1}   ## In case it got left over from a previous failed run.
echo TEMP_DIR_CL1 is: ${TEMP_DIR_CL1}
mkdir ${TEMP_DIR_CL1}

## bam1s is a full path to the sorted, single chromosome bam file generated by sort_chrom.sbatch
## BASE1 looks like this (for each chromosome):             chr10.1.bam
## bam2s and BASE2 are similar.  BASE2 can look like this:  chr10_KL568015v1_random.2.bam
while read -r bam1s; do
    while read -r bam2s; do
        BASE1=${bam1s##*/}    ## Get rid of the path before the bam file name.
        BASE2=${bam2s##*/}

        BASE1=${BASE1#*\.}   ## Get rid of the sample name.
        BASE2=${BASE2#*\.}

        echo ${bam1s} ${bam2s} "${TEMP_DIR_CL1}/${BASE1}___${BASE2}"  >> "${sampleOutdir}/bams/array_list"
    done < "${sampleOutdir}/bams/chr_list_2"
done < "${sampleOutdir}/bams/chr_list_1"

################################################################################################
## Create the small chromosome bam files, and sort them:

export_vars="sampleOutdir=${sampleOutdir}"
export_vars="${export_vars},BAM=${bamname1}"
export_vars="${export_vars},BAM_K=${bam1_keep_flags}"      # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_E=${bam1_exclude_flags}"   # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_N=1"
export_vars="${export_vars},input_to_samtools2=${chrom_list1_input_to_samtools}"

num_lines=$(wc -l < "${sampleOutdir}/log/chrom_list1_simple")

JOB_ID0=$(sbatch --parsable --export=ALL,${export_vars} --array="1-${num_lines}" --job-name=sort_bamintersect \
        --output="${sampleOutdir}/log/sort_bamintersect_1.${sample_name}.o_%A_%a" "${src}/sort_bamintersect.sh")


export_vars="sampleOutdir=${sampleOutdir}"
export_vars="${export_vars},BAM=${bamname2}"
export_vars="${export_vars},BAM_K=${bam2_keep_flags}"      # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_E=${bam2_exclude_flags}"   # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_N=2"
export_vars="${export_vars},input_to_samtools2=${chrom_list2_input_to_samtools}"

num_lines=$(wc -l < "${sampleOutdir}/log/chrom_list2_simple")

JOB_ID1=$(sbatch --parsable --export=ALL,${export_vars} --array="1-${num_lines}" --job-name=sort_bamintersect \
        --output="${sampleOutdir}/log/sort_bamintersect_2.${sample_name}.o_%A_%a" "${src}/sort_bamintersect.sh")

JOB_ID0=":${JOB_ID0}:${JOB_ID1}"

################################################################################################
## The big array job:

n1=$(wc -l < "${sampleOutdir}/log/chrom_list1_simple")
n2=$(wc -l < "${sampleOutdir}/log/chrom_list2_simple")
let "array_size = ${n1} * ${n2}"

export_vars="sampleOutdir=${sampleOutdir}"
export_vars="${export_vars},src=${src}"
export_vars="${export_vars},reads_match=${reads_match}"
export_vars="${export_vars},make_csv=${make_csv}"

JOB_ID1=$(sbatch --parsable --dependency=afterok${JOB_ID0} --export="ALL,${export_vars}" --array="1-${array_size}" \
          --output="${sampleOutdir}/log/bamintersect.${sample_name}.o_%A_%a" --job-name=bamintersect "${src}/bamintersect.sh")

################################################################################################
## Cleanup:

export_vars="sampleOutdir=${sampleOutdir}"
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

sbatch --parsable --dependency=afterok:${JOB_ID1} --export=ALL,${export_vars} --job-name=merge_bamintersect \
       --output="${sampleOutdir}/log/merge_bamintersect.${sample_name}.o_%j" "${src}/merge_bamintersect.sh"

