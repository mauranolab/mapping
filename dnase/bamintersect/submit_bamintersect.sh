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

TMPDIR=`mktemp -d`
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
    make_csv="--make_csv"
    make_table=True
    reads_match=""
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
        --outdir)
            sampleOutdir=$2 ; shift 2 ;;
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
            reads_match="--same" ; shift 1 ;;
        --do_not_make_csv)
            make_csv="" ; shift 1 ;;
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

    # Get the correct HA zones (0-based bed style numbers, not UCSC positions):
    #     The 5p and 3p outer boundaries of the HPRT1 HAs are as of Mar 29, 2019.
    #     See: /vol/cegs/sequences/hg38/HPRT1/HPRT1_assembly.bed
    #
    #     The 5p and 3p outer boundaries of the Sox2 HAs are as of Aug 9, 2019.
    #     See: /vol/cegs/sequences/mm10/Sox2/Sox2_assembly.bed
    #
    # Also see: \\research-cifs.nyumc.org\Research\CEGS\Sequences\Landing Pads Summary Table.xlsx (8 JUL 2019)
    #
    # Also get the name of the deletion gene, which is used with the bedops closest-features function to 
    # determine which non-cegsvectors reads are candidates for being "informative".
    if echo "${cegsgenome}" | grep -q "LP058" ; then
        # HA lengths: 1032/1014 (LP058)
        bam1_5p_HA="chrX:134458914-134459946"
        bam1_3p_HA="chrX:134501642-134502656"
        deletion_gene=HPRT1
    elif echo "${cegsgenome}" | grep -q "LP087" ; then
        # HA lengths: 251/250 (LP061)
        bam1_5p_HA="chrX:134459695-134459946"
        bam1_3p_HA="chrX:134501642-134501892"
        deletion_gene=HPRT1
    elif echo "${cegsgenome}" | grep -q "LP123" ; then
        # HA lengths: 251/250 (LP061)
        bam1_5p_HA="chrX:134459695-134459946"
        bam1_3p_HA="chrX:134501642-134501892"
        deletion_gene=HPRT1
    elif echo "${cegsgenome}" | grep -q "LP062" ; then
        # HA lengths: 100/100 (LP062)
        bam1_5p_HA="chrX:134459846-134459946"
        bam1_3p_HA="chrX:134501642-134501742"
        deletion_gene=HPRT1
    elif echo "${cegsgenome}" | grep -q "LP097a" ; then
        # HA lengths: 142/278 (LP097a)
        bam1_5p_HA="chr3:34631310-34631452"
        bam1_3p_HA="chr3:34768108-34768386"
        deletion_gene=Sox2ot
    elif echo "${cegsgenome}" | grep -q "LP131a" ; then
        # HA lengths: 142/278 (LP131a)
        bam1_5p_HA="chr3:34631310-34631452"
        bam1_3p_HA="chr3:34768108-34768386"
        deletion_gene=Sox2ot
    elif echo "${cegsgenome}" | grep -q "LP097b" ; then
        # HA lengths: 142/183 (LP097b)
        bam1_5p_HA="chr3:34631310-34631452"
        bam1_3p_HA="chr3:34774117-34774300"
        deletion_gene=Sox2ot
    elif echo "${cegsgenome}" | grep -q "LP131b" ; then
        # HA lengths: 142/183 (LP131b)
        bam1_5p_HA="chr3:34631310-34631452"
        bam1_3p_HA="chr3:34774117-34774300"
        deletion_gene=Sox2ot
    else
        bam1_5p_HA="NA"
        bam1_3p_HA="NA"
        deletion_gene="NA"
        deletion_range="NA"
        echo "No HA match for ${sample_name}     genome is: ${cegsgenome}"
    fi


    # Compute deletion range coordinates.
    # This is still used in filter_tsv.sh, but I'd like to discuss that.
    if [ "${bam1_5p_HA}" != "NA" ]; then
        IFS=':' read -ra ADDR <<< "${bam1_5p_HA}"
        del_chrom=${ADDR[0]}
        del_range=${ADDR[1]}
        IFS='-' read -ra ADDR <<< "${del_range}"
        del_5p=${ADDR[1]}

        IFS=':' read -ra ADDR <<< "${bam1_3p_HA}"
        del_range=${ADDR[1]}
        IFS='-' read -ra ADDR <<< "${del_range}"
        del_3p=${ADDR[0]}

        deletion_range="${del_chrom}:${del_5p}-${del_3p}"
    fi


    # A bed file of ranges to delete at the end of the process.  cegsvector reads that fall in these ranges may actually be
    # a misleading mapping of genomic reads, so they are considered uninformative.
    if echo "${cegsgenome}" | grep -q "LP058" ; then
        exclude_regions_from_counts="/vol/mauranolab/cadlej01/projects/LP_Integration/LP_uninformative_regions/LP_vs_hg38.bed"
        cegsLP=True
    elif echo "${cegsgenome}" | grep -q "LP097a" ; then
        exclude_regions_from_counts="/vol/mauranolab/cadlej01/projects/LP_Integration/LP_uninformative_regions/LP_vs_mm10.bed"
        cegsLP=True
    elif echo "${cegsgenome}" | grep -q "LP097b" ; then
        exclude_regions_from_counts="/vol/mauranolab/cadlej01/projects/LP_Integration/LP_uninformative_regions/LP_vs_mm10.bed"
        cegsLP=True
    elif echo "${cegsgenome}" | grep -q "LP131a" ; then
        exclude_regions_from_counts="/vol/mauranolab/cadlej01/projects/LP_Integration/LP_uninformative_regions/LP_vs_mm10.bed"
        cegsLP=True
    elif echo "${cegsgenome}" | grep -q "LP131b" ; then
        exclude_regions_from_counts="/vol/mauranolab/cadlej01/projects/LP_Integration/LP_uninformative_regions/LP_vs_mm10.bed"
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
        # No uninformative regions
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
    else
        if echo "${annotationgenome}" | grep -q "mm10" ; then
            # Arrive here when $cegsgenome == pSpCas9, and $annotationgenome == mm10
            exclude_regions_from_counts="NA"
        fi
    fi

# Done getting the derived data.
########################################################
# Make some directories to hold the final output for this sample.

if [ "${clear_logs}" = "True" ]; then
    rm -rf "${sampleOutdir}/log"  # Clear out the log directory, if it's there.  Avoids piling up old log files during development.
fi
mkdir -p ${sampleOutdir}/log


########################################################
# Make a directory structure to hold files passed between jobs.
# It is built within ${sampleOutdir} so it can be persistent.
# It gets removed at the end of the merge_bamintersect.sh script.
INTERMEDIATEDIR=`mktemp -d --tmpdir=${sampleOutdir}`   # INTERMEDIATEDIR contains no trailing slash.

# Subdirectories will hold output from bamintersect.py
mkdir ${INTERMEDIATEDIR}/bamintersectPyOut

# Make a directory for all the small sorted, chromosome bam files generated by sort_bamintersect.sh
mkdir ${INTERMEDIATEDIR}/sorted_bams

########################################################
# Make the chromosome lists which will eventually drive the array jobs.

samtools idxstats ${bamname1} | awk '$1 != "*" {print $1}' > "${sampleOutdir}/log/${sample_name}.chrom_list1"
num_lines=$(wc -l < "${sampleOutdir}/log/${sample_name}.chrom_list1")

if [ "${num_lines}" -ge "5" ]; then
    ## Get the simple chromosome name mask
    grep -E '^chr[0-9]*$|^chr[XYM]$' "${sampleOutdir}/log/${sample_name}.chrom_list1" | sed 's/^/^/' | sed 's/$/$/' > "${TMPDIR}/${sample_name}.chrom_list1_simple_mask"

    ## Apply mask to chomosome names
    grep -v -f "${TMPDIR}/${sample_name}.chrom_list1_simple_mask" "${sampleOutdir}/log/${sample_name}.chrom_list1" > "${TMPDIR}/${sample_name}.chrom_list1_long"

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

## Get the simple chromosome name mask
grep -E '^chr[0-9]*$|^chr[XYM]$' "${sampleOutdir}/log/${sample_name}.chrom_list2" | sed 's/^/^/' | sed 's/$/$/' > "${TMPDIR}/${sample_name}.chrom_list2_simple_mask"

## Apply mask to chomosome names
grep -v -f "${TMPDIR}/${sample_name}.chrom_list2_simple_mask" "${sampleOutdir}/log/${sample_name}.chrom_list2" > "${TMPDIR}/${sample_name}.chrom_list2_long"

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
        x=${BAM%.bam}         # Cut off the trailing "bam"
        y=${x##*.}            # Assign the base to y. This will be one of genotypes, like "hg38_full"
        z=${x%.${y}}          # Cut off the trailing ${y}
        short_name=${z##*/}   # The remaining base is the short name.

        BAM_OUT="${INTERMEDIATEDIR}/sorted_bams/${short_name}.${chrom}.${BAM_N}.bam"
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

num_lines=$(wc -l < "${sampleOutdir}/log/${sample_name}.chrom_list1_simple")

# /vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/test_java contains code and comments regarding java thread
# problems generated by multiple calls of sort_bamintersect. qos=low probably is good enough to handle the issues.
qsub -S /bin/bash -cwd -terse -j y --export=ALL,${export_vars} -N sort_bamintersect_1.${sample_name} -o ${sampleOutdir}/log --qos=low -t 1-${num_lines} ${src}/sort_bamintersect.sh > ${sampleOutdir}/sgeid.sort_bamintersect_1


export_vars="sampleOutdir=${sampleOutdir}"
export_vars="${export_vars},BAM=${bamname2}"
export_vars="${export_vars},BAM_K=${bam2_keep_flags}"      # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_E=${bam2_exclude_flags}"   # Required. Use "0" if necessary.
export_vars="${export_vars},BAM_N=2"
export_vars="${export_vars},input_to_samtools2=${chrom_list2_input_to_samtools}"
export_vars="${export_vars},INTERMEDIATEDIR=${INTERMEDIATEDIR}"
export_vars="${export_vars},sample_name=${sample_name}"

num_lines=$(wc -l < "${sampleOutdir}/log/${sample_name}.chrom_list2_simple")

qsub -S /bin/bash -cwd -terse -j y --export=ALL,${export_vars} -N sort_bamintersect_2.${sample_name} -o ${sampleOutdir}/log --qos=low -t 1-${num_lines} ${src}/sort_bamintersect.sh > ${sampleOutdir}/sgeid.sort_bamintersect_2


################################################################################################
## The big array job:
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
echo "${bamname1}" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "${bamname2}" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

echo "Deletion Gene: ${deletion_gene}" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

echo "${deletion_gene} Deletion Range is: ${deletion_range}" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

echo "${deletion_gene} HAs:" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "    ${bam1_5p_HA}" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "    ${bam1_3p_HA}" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

samtools idxstats ${bamname1} > "${INTERMEDIATEDIR}/${sample_name}.counts.anc_info.txt"
num_lines=$(wc -l < "${INTERMEDIATEDIR}/${sample_name}.counts.anc_info.txt")
let "num_lines = num_lines - 1"
echo "samtools idx output for first bam file: reference sequence name, sequence length, # mapped reads and # unmapped reads." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
head "-${num_lines}" "${INTERMEDIATEDIR}/${sample_name}.counts.anc_info.txt" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo "working file" > "${INTERMEDIATEDIR}/${sample_name}.counts.anc_info.txt"

################################################################################################
## Merge output from the array jobs.
echo "Submitting merge job"

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
export_vars="${export_vars},deletion_gene=${deletion_gene}"
export_vars="${export_vars},deletion_range=${deletion_range}"
export_vars="${export_vars},INTERMEDIATEDIR=${INTERMEDIATEDIR}"

qsub -S /bin/bash -cwd -terse -j y -hold_jid `cat ${sampleOutdir}/sgeid.merge_bamintersect` --export=ALL,${export_vars} -N merge_bamintersect.${sample_name} -o ${sampleOutdir}/log ${src}/merge_bamintersect.sh
rm -f ${sampleOutdir}/sgeid.merge_bamintersect

