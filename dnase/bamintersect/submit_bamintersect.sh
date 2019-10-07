#!/bin/bash
############################################################
# Setting up.

set -eu -o pipefail

############################################################
# Where is this file located? We use this info to find and set other resources.
#
# This version of the code for setting src does not necessarily return a full path to the source directory.
# It returns the path that was used to call submit_bamintersect.sh, which could be like this:
#    ./submit_bamintersect.sh
#              or
#    ../submit_bamintersect.sh
#              or 
#    some other non-full path
#
# In the above scenarios, src="." or src=".."
# The current versions of the bamintersect related scripts will still work properly, but
# if future revisions entail changing working directories or parsing src, then this should
# be kept in mind.
src=$( dirname "${BASH_SOURCE[0]}" )
echo "src is: ${src}"

# The below is for cases where this script is not being called as part of a larger pipeline.
# Avoids conflicts should "sample_name" be identical in simultaneous runs of this script.
if [ ${TMPDIR} = "/tmp" ]; then
    TMPDIR=`mktemp -d`   # TMPDIR has no trailing slash
fi
############################################################

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
     --integrationsite {sitename_HA#1_HA#2   sitename is like HPRT1. HA#s refer to numbers in the HomologyArms.bed file}
     --bam1 {full path to bam file #1}
     --bam1genome {something like LP123, mm10, or hg38. Others are valid as well.}
     --bam2 {full path to bam file #2}
     --bam2genome {something like LP123, mm10, or hg38. Others are valid as well.}
     --outdir {outdir} output directory

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
    integrationsite:
    outdir:
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
    make_csv="--make_csv"
    make_table=True
    reads_match=""
    integrationsite="NA"
    
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
        --integrationsite)
            integrationsite=$2 ; shift 2 ;;
        --outdir)
            sampleOutdir=$2 ; shift 2 ;;
        --bam1)
            bamname1=$2 ; shift 2 ;;
        --bam1genome)
            bam1genome=$2 ; shift 2 ;;
        --bam1_keep_flags)
            bam1_keep_flags=$2 ; shift 2 ;;
        --bam1_exclude_flags)
            bam1_exclude_flags=$2 ; shift 2 ;;
        --bam2)
            bamname2=$2 ; shift 2 ;;
        --bam2genome)
            bam2genome=$2 ; shift 2 ;;
        --bam2_keep_flags)
            bam2_keep_flags=$2 ; shift 2 ;;
        --bam2_exclude_flags)
            bam2_exclude_flags=$2 ; shift 2 ;;
        --reads_match)
            reads_match="--same" ; shift 1 ;;
        --do_not_make_csv)
            make_csv="" ; shift 1 ;;
        --do_not_make_table)
            make_table=False ; shift 1 ;;
        -h|--help) usage; shift ; exit 0 ;;
        --) shift ; break ;;
        *) echo "getopt internal error!" ; exit 1 ;;
    esac
done

# End of getopt section.
################################################
sample_name="${sample_name}.${bam1genome}_vs_${bam2genome}"
echo "final sample_name is: ${sample_name}"
################################################
# Derive the "excluded regions" from bam1genome and bam2genome: 

# A bed file of ranges to delete at the end of the process.  cegsvector reads that fall in these ranges may actually be
# a misleading mapping of genomic reads, so they are considered uninformative.

if echo "${bam1genome}" | grep -q "^LP[0-9]\+" ; then
    genome2exclude="/vol/mauranolab/cadlej01/projects/bamintersect/LP_uninformative_regions/LP_vs_${bam2genome}.bed"
else
    genome2exclude="/vol/mauranolab/cadlej01/projects/bamintersect/LP_uninformative_regions/${bam1genome}_vs_${bam2genome}.bed"
fi

########################################################
# Make some directories to hold the final output for this sample.

mkdir -p "${sampleOutdir}/log"

########################################################
# Make a directory structure to hold files passed between jobs.
# It is built within ${sampleOutdir} so it can be persistent.
# It gets removed at the end of the merge_bamintersect.sh script.
INTERMEDIATEDIR="${sampleOutdir}/tmp.${sample_name}.${bam1genome}"   # INTERMEDIATEDIR contains no trailing slash.
mkdir ${INTERMEDIATEDIR}

# Subdirectories will hold output from bamintersect.py
mkdir ${INTERMEDIATEDIR}/bamintersectPyOut

# Make a directory for all the small sorted, chromosome bam files generated by sort_bamintersect.sh
mkdir ${INTERMEDIATEDIR}/sorted_bams

########################################################
# Make the chromosome lists which will eventually drive the array jobs.

samtools idxstats ${bamname1} | awk '$1 != "*" {print $1}' > "${sampleOutdir}/log/${sample_name}.chrom_list1"
num_lines=$(wc -l < "${sampleOutdir}/log/${sample_name}.chrom_list1")

if [ "${num_lines}" -ge "5" ]; then
    ## Apply the simple chromosome name mask
    grep -E '^chr[0-9]*$|^chr[XYM]$' "${sampleOutdir}/log/${sample_name}.chrom_list1" | sed 's/^/^/' | sed 's/$/$/' |
    grep -v -f - "${sampleOutdir}/log/${sample_name}.chrom_list1" > "${TMPDIR}/${sample_name}.chrom_list1_long"
    
    ## We'll strip the pipes out in sort_bamintersect.sh
    chrom_list1_input_to_samtools=$(cat "${TMPDIR}/${sample_name}.chrom_list1_long" | paste -sd "|" -)
    
    ## Get simple chromosome names
    grep -E '^chr[0-9]*$|^chr[XYM]$' "${sampleOutdir}/log/${sample_name}.chrom_list1" >  "${sampleOutdir}/log/${sample_name}.chrom_list1_simple"
    num_lines=$(wc -l < "${TMPDIR}/${sample_name}.chrom_list1_long")
    if [ "${num_lines}" -ge "1" ]; then
        echo "all_other" >> "${sampleOutdir}/log/${sample_name}.chrom_list1_simple"
    fi
else
    # A small number of chromosomes is an indicator that this is a LP/PL file.
    # There will be no "all_other" in this file.
    cp "${sampleOutdir}/log/${sample_name}.chrom_list1" "${sampleOutdir}/log/${sample_name}.chrom_list1_simple"
    chrom_list1_input_to_samtools="NA"   # To avoid a pipefail.
fi

samtools idxstats ${bamname2} | awk '$1 != "*" {print $1}' > "${sampleOutdir}/log/${sample_name}.chrom_list2"
num_lines=$(wc -l < "${sampleOutdir}/log/${sample_name}.chrom_list2")

## Apply the simple chromosome name mask
grep -E '^chr[0-9]*$|^chr[XYM]$' "${sampleOutdir}/log/${sample_name}.chrom_list2" | sed 's/^/^/' | sed 's/$/$/' |
grep -v -f - "${sampleOutdir}/log/${sample_name}.chrom_list2" > "${TMPDIR}/${sample_name}.chrom_list2_long"

## We'll strip the pipes out in sort_bamintersect.sh
chrom_list2_input_to_samtools=$(cat "${TMPDIR}/${sample_name}.chrom_list2_long" | paste -sd "|" -)

## Get simple chromosome names
grep -E '^chr[0-9]*$|^chr[XYM]$' "${sampleOutdir}/log/${sample_name}.chrom_list2" >  "${sampleOutdir}/log/${sample_name}.chrom_list2_simple"
num_lines=$(wc -l < "${TMPDIR}/${sample_name}.chrom_list2_long")
if [ "${num_lines}" -ge "1" ]; then
    echo "all_other" >> "${sampleOutdir}/log/${sample_name}.chrom_list2_simple"
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
        BAM_OUT="${INTERMEDIATEDIR}/sorted_bams/${sample_name}.${chrom}.${BAM_N}.bam"
        echo ${BAM_OUT} >> "${INTERMEDIATEDIR}/sorted_bams/chr_list_${BAM_N}"
    done < "${sampleOutdir}/log/${sample_name}.chrom_list${BAM_N}_simple"
done


## bam1s is a full path to the sorted, single chromosome bam file generated by sort_chrom.sbatch
## BASE1 looks like this (for each chromosome):             chr10.1.bam
## bam2s and BASE2 are similar.  BASE2 can look like this:  chr10_KL568015v1_random.2.bam
while read -r bam1s; do
    while read -r bam2s; do
        BASE1=${bam1s##*/}    ## Get rid of the path before the bam file name.
        BASE2=${bam2s##*/}

        BASE1=${BASE1#*\.}   ## Get rid of the sample name.
        BASE2=${BASE2#*\.}

        # Not sure if this is a problem.  Got inconsistent results.  Once saw the
        # merge sbatch have "Dependency never satisfied".  But then it went away next 
        # time.  Maybe due to buffering here? stdbuf should flush all the non-sbatch output
        # from here and above, if it becomes a problem. The sbatch's will flush when they exit.
        # stdbuf --output=0 echo ${bam1s} ${bam2s} "${INTERMEDIATEDIR}/bamintersectPyOut/${BASE1}___${BASE2}"  >> "${TMPDIR}/TMPDIR_bams/inputs.array_list.txt"
        echo ${bam1s} ${bam2s} "${INTERMEDIATEDIR}/bamintersectPyOut/${BASE1}___${BASE2}"  >> "${INTERMEDIATEDIR}/sorted_bams/inputs.array_list.txt"

    done < "${INTERMEDIATEDIR}/sorted_bams/chr_list_2"
done < "${INTERMEDIATEDIR}/sorted_bams/chr_list_1"

################################################################################################
## Create the small chromosome bam files, and sort them:
echo "Submitting sort jobs"

export_vars="sampleOutdir=${sampleOutdir}"
export_vars="${export_vars},BAM=${bamname1}"
export_vars="${export_vars},BAM_K=${bam1_keep_flags}"      # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_E=${bam1_exclude_flags}"   # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_N=1"
export_vars="${export_vars},input_to_samtools2=${chrom_list1_input_to_samtools}"
export_vars="${export_vars},INTERMEDIATEDIR=${INTERMEDIATEDIR}"
export_vars="${export_vars},sample_name=${sample_name}"
export_vars="${export_vars},logdir=${sampleOutdir}/log"

num_lines=$(wc -l < "${sampleOutdir}/log/${sample_name}.chrom_list1_simple")

# /vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/test_java contains code and comments regarding java thread
# problems generated by multiple calls of sort_bamintersect. qos=low probably is good enough to handle the issues.
qsub -S /bin/bash -cwd -terse -j y --export=ALL,${export_vars} -N sort_bamintersect_1.${sample_name} -o ${sampleOutdir}/log --qos=low -t 1-${num_lines} ${src}/sort_bamintersect.sh > ${sampleOutdir}/sgeid.sort_bamintersect_1

echo done with sort_

export_vars="sampleOutdir=${sampleOutdir}"
export_vars="${export_vars},BAM=${bamname2}"
export_vars="${export_vars},BAM_K=${bam2_keep_flags}"      # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_E=${bam2_exclude_flags}"   # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_N=2"
export_vars="${export_vars},input_to_samtools2=${chrom_list2_input_to_samtools}"
export_vars="${export_vars},INTERMEDIATEDIR=${INTERMEDIATEDIR}"
export_vars="${export_vars},sample_name=${sample_name}"
export_vars="${export_vars},logdir=${sampleOutdir}/log"

num_lines=$(wc -l < "${sampleOutdir}/log/${sample_name}.chrom_list2_simple")

qsub -S /bin/bash -cwd -terse -j y --export=ALL,${export_vars} -N sort_bamintersect_2.${sample_name} -o ${sampleOutdir}/log --qos=low -t 1-${num_lines} ${src}/sort_bamintersect.sh > ${sampleOutdir}/sgeid.sort_bamintersect_2

echo done with 2nd sort_

################################################################################################
# Call bamintersect.py repeatedly via the bamintersect.sh array job:
echo "Submitting bamintersect job"

n1=$(wc -l < "${sampleOutdir}/log/${sample_name}.chrom_list1_simple")
n2=$(wc -l < "${sampleOutdir}/log/${sample_name}.chrom_list2_simple")
let "array_size = ${n1} * ${n2}"

export_vars="src=${src}"
export_vars="${export_vars},reads_match=${reads_match}"
export_vars="${export_vars},make_csv=${make_csv}"
export_vars="${export_vars},INTERMEDIATEDIR=${INTERMEDIATEDIR}"

qsub -S /bin/bash -cwd -terse -j y -hold_jid `cat ${sampleOutdir}/sgeid.sort_bamintersect_1 ${sampleOutdir}/sgeid.sort_bamintersect_2 | perl -pe 's/\n/,/g;'` --export=ALL,${export_vars} -N bamintersect.${sample_name} -o ${sampleOutdir}/log -t 1-${array_size} ${src}/bamintersect.sh > ${sampleOutdir}/sgeid.merge_bamintersect
rm -f ${sampleOutdir}/sgeid.sort_bamintersect_1 ${sampleOutdir}/sgeid.sort_bamintersect_2

################################################################################################
## Send some summary info to the anc file:
echo "bam files:" > "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "First:  ${bamname1}" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "Second: ${bamname2}" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

# Make the filter files for merge_bamintersect.sh and filter_tsv.sh.orig
if [ "${integrationsite}" != "null" ] && [ "${integrationsite}" != "NA" ]; then
    IFS='_' read integrationSiteName HA1 HA2 <<< "${integrationsite}"
    echo "${integrationSiteName} HAs:" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    HA_file="/vol/cegs/sequences/${bam2genome}/${integrationSiteName}/${integrationSiteName}_HomologyArms.bed"
    grep "${integrationSiteName}_${HA1}$" ${HA_file} >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    grep "${integrationSiteName}_${HA2}$" ${HA_file} >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

    grep "${integrationSiteName}_${HA1}$" ${HA_file} > "${INTERMEDIATEDIR}/genome1exclude.bed"
    grep "${integrationSiteName}_${HA2}$" ${HA_file} >> "${INTERMEDIATEDIR}/genome1exclude.bed"

    # Get the deletion range coordinates:
    IFS=$'\t' read chrom HA1_5p HA1_3p all_other <<< "$(head -n1 "${INTERMEDIATEDIR}/genome1exclude.bed")"
    IFS=$'\t' read chrom HA2_5p HA2_3p all_other <<< "$(tail -n1 "${INTERMEDIATEDIR}/genome1exclude.bed")"
    echo "${chrom}"$'\t'"${HA1_3p}"$'\t'"${HA2_5p}" > "${INTERMEDIATEDIR}/deletion_range.bed"
else
    echo "integrationsite is null/NA" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    echo "No HAs are available, and so there is no Deletion Range." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

    integrationsite="NA"
    touch "${INTERMEDIATEDIR}/genome1exclude.bed"
    touch "${INTERMEDIATEDIR}/deletion_range.bed"
fi

if [ ! -f ${genome2exclude} ]; then
    touch "${INTERMEDIATEDIR}/empty_file.bed"
    genome2exclude="${INTERMEDIATEDIR}/empty_file.bed"
fi
echo "The Exclude Regions file is: $(basename ${genome2exclude})  [${genome2exclude}]" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

################################################################################################
samtools idxstats ${bamname1} > "${TMPDIR}/${sample_name}.counts.anc_info.txt"
num_lines=$(wc -l < "${TMPDIR}/${sample_name}.counts.anc_info.txt")
let "num_lines = num_lines - 1"
echo "samtools idx output for first bam file: reference sequence name, sequence length, # mapped reads and # unmapped reads." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
head "-${num_lines}" "${TMPDIR}/${sample_name}.counts.anc_info.txt" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

if [ ! -f ${genome2exclude} ]; then
    touch "${INTERMEDIATEDIR}/empty_file.bed"
    genome2exclude="${INTERMEDIATEDIR}/empty_file.bed"
fi
echo "The Exclude Regions file is: $(basename ${genome2exclude})  [${genome2exclude}]" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

################################################################################################
samtools idxstats ${bamname1} > "${TMPDIR}/${sample_name}.counts.anc_info.txt"
num_lines=$(wc -l < "${TMPDIR}/${sample_name}.counts.anc_info.txt")
let "num_lines = num_lines - 1"
echo "samtools idx output for first bam file: reference sequence name, sequence length, # mapped reads and # unmapped reads." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
head "-${num_lines}" "${TMPDIR}/${sample_name}.counts.anc_info.txt" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

################################################################################################
## Merge output from the array jobs.
echo "Submitting merge job"

export_vars="sampleOutdir=${sampleOutdir}"
export_vars="${export_vars},src=${src}"
export_vars="${export_vars},genome2exclude=${genome2exclude}"
export_vars="${export_vars},sample_name=${sample_name}"
export_vars="${export_vars},bam1genome=${bam1genome}"
export_vars="${export_vars},bam2genome=${bam2genome}"
export_vars="${export_vars},make_csv=${make_csv}"
export_vars="${export_vars},make_table=${make_table}"
export_vars="${export_vars},INTERMEDIATEDIR=${INTERMEDIATEDIR}"
export_vars="${export_vars},integrationsite=${integrationsite}"
export_vars="${export_vars},logdir=${sampleOutdir}/log"

qsub -S /bin/bash -cwd -terse -j y -hold_jid `cat ${sampleOutdir}/sgeid.merge_bamintersect` --export=ALL,${export_vars} -N merge_bamintersect.${sample_name} -o ${sampleOutdir}/log ${src}/merge_bamintersect.sh
rm -f ${sampleOutdir}/sgeid.merge_bamintersect

