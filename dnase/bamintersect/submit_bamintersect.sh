#!/bin/bash
############################################################
# Setting up.

set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'


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

# The below is for cases where this script is not being called as part of a larger pipeline.
# Avoids conflicts should "sample_name" be identical in simultaneous runs of this script.
if [ ${TMPDIR} = "/tmp" ]; then
    TMPDIR=`mktemp -d`   # TMPDIR has no trailing slash
fi
############################################################

module load samtools
module load bedops

############################################################
# Start of "getopt" section
############################################################
# For use with the getopt "-h" flag:

function usage {
cat << USAGE

Usage: $(basename "$0") [Options]
  Required Arguments:
    --sample_name {sample name}
    --integrationsite {sitename_HA#1_HA#2   Must be canonical name (sitename e.g. "HPRT1"; HA#s refer to numbers in the HomologyArms.bed file) or null.}
    --bam1 {full path to bam file #1}
    --bam1genome [ LP123, mm10, hg38, rn6, ... ]
    --bam2 {full path to bam file #2}
    --bam2genome [ LP123, mm10, hg38, rn6, ... ]
    --outdir {outdir} output directory
    
    Genome notes:
    1) The mammalian genome is customarily passed as bam1. This affects:
        Files for uninformative read masking differ
        HA table is run on bam1, backbone analysis on bam2
        Reads must have unmapped mates in bam2 (we can't require this for bam1 since typically the payload matches some/all of the reference so it would mask internal reads)
        Reads for counts tables are clustered first in bam1, then separate bam2region clusters are defined for each bam1region
    
    2) Mammalian genomes should use the annotationgenome name (hg38, not hg38_full)
    
  Optional Arguments:
    --max_mismatches {The max number of mismatches a read is allowed to have vs the reference. The value is in the read's NM tag. Default = 1}
    
    --noReqFullyAligned {If set, the default requirement that reads must be fully aligned is disabled. Default = Not set}
    
    --bam1_keep_flags {a sam format flag value required of all bam1 reads to be included in the analysis. Default = 1}
        # Require read is paired
    
    --bam1_exclude_flags {a sam format flag value for bam1 reads to be excluded from the analysis. Default = 3078}
        # Exclude mapped in proper pair
        # Exclude read unmapped.
        # Exclude PCR or optical duplicates.
        # Exclude supplementary aligments.
    
    --bam2_keep_flags {a sam format flag value required of all bam2 reads to be included in the analysis. Default = 9}
    
    --bam2_exclude_flags {a sam format flag value for bam2 reads to be excluded from the analysis. Default = 3078}
    
    --reads_match {no argument}
           # This flag tells bam_intersect.py how to match the reads from each bam file.
           # If set, reads should be matched like [read1 read1], or [read2 read2].
           # If unset, reads are paired like [read1 read2], or [read2 read1].
           # Set this option for unpaired reads.
           # The default is to have it unset.
    
     --no_bed {no argument}
           # Main output is a single bed file (unset), or 2 bam files (set).
           # The default is to have it unset.
    
     --no_countstable {no argument}
           # Output a "counts.txt" table summarizing results (unset), or no table (set).
           # The default is to have it unset.
    
     --verbose {no argument}
           # Prints out debugging statements
    
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
    max_mismatches:
    noReqFullyAligned
    reads_match
    no_csv
    no_countstable
    verbose
    help
)

# These arguments are required arguments, as opposed to the options with required arguments described above:
rqd_arg_names=(
    sample_name
    integrationsite
    outdir
    bam1
    bam1genome
    bam2
    bam2genome
)

declare -A rqd_arg_status
for i in "${rqd_arg_names[@]}"; do
    # We start by assuming the argument is not there ("0").
    rqd_arg_status["${i}"]="0"
done

############################################################
# Parsing the getopt arguments:

long_args=$(printf "%s," "${long_arg_list[@]}")   # Turn the arg array into a string with commas.
long_args=$(echo ${long_args} | sed 's/,$/ /')    # Get rid of final comma.

# There must be at least one option associated with "-o", or you need to deal with getopt's default behavior.
# Default behavior: The first parameter of getopt that does not start with a '-' (and is not an option argument) is used as the short options string.
# getopt returns a string where unmanaged input parameters are removed (which could be used elsewhere), and arguments for options are filled in
options=$(getopt -o h --long "${long_args}" -n "[submit_bamintersect] ERROR" -- "$@")
# Catch getopt errors here if not using "set -e"

echo "Parameters: ${options}"


#Reset the positional parameters based on the output of getopt
eval set -- "${options}"

# getopt default values:
    make_bed="--make_bed"
    make_table=True
    reads_match=""
    integrationsite="null"
    bam1_keep_flags="1"
    bam1_exclude_flags="3078"
    bam2_keep_flags="9"
    bam2_exclude_flags="3078"
    verbose=False
    max_mismatches=1
    ReqFullyAligned="--ReqFullyAligned"


# Extract options and their arguments into variables.
while true ; do
    case "$1" in
        --sample_name)
            sample_name=$2 ; shift 1 ; rqd_arg_status[sample_name]="1" ;;
        --integrationsite)
            integrationsite=$2 ; shift 1 ; rqd_arg_status[integrationsite]="1" ;;
        --outdir)
            sampleOutdir=$2 ; shift 1 ; rqd_arg_status[outdir]="1" ;;
        --bam1)
            bam1=$2 ; shift 1 ; rqd_arg_status[bam1]="1" ;;
        --bam1genome)
            bam1genome=$2 ; shift 1 ; rqd_arg_status[bam1genome]="1" ;;
        --bam1_keep_flags)
            bam1_keep_flags=$2 ; shift 1 ;;
        --bam1_exclude_flags)
            bam1_exclude_flags=$2 ; shift 1 ;;
        --bam2)
            bam2=$2 ; shift 1 ; rqd_arg_status[bam2]="1" ;;
        --bam2genome)
            bam2genome=$2 ; shift 1 ; rqd_arg_status[bam2genome]="1" ;;
        --bam2_keep_flags)
            bam2_keep_flags=$2 ; shift 1 ;;
        --bam2_exclude_flags)
            bam2_exclude_flags=$2 ; shift 1 ;;
        --max_mismatches)
            max_mismatches=$2 ; shift 1 ;;
        --noReqFullyAligned)
            ReqFullyAligned="" ;;
        --reads_match)
            reads_match="--reads_match" ;;
        --no_bed)
            make_bed="" ;;
        --no_countstable)
            make_table=False ;;
        --verbose)
            verbose=True ;;
        -h|--help) usage; exit 0 ;;
        --) shift ; break ;;
        *) echo "ERROR getopt internal error $1!" ; exit 1 ;;
    esac
    shift
done


# Check that all required arguments have been provided.
ierr=0
for i in "${rqd_arg_names[@]}"; do
    if [ ${rqd_arg_status["${i}"]} != "1" ]; then
        echo "ERROR Missing required arg to submit_bamintersect.sh :  ${i}"
        ierr=1
    fi
done

if [ ${ierr} = "1" ]; then
    # We are missing at least one required argument.
    echo "ERROR Missing required required argument."
    exit 1
fi

# End of getopt section.
################################################
sample_name="${sample_name}.${bam1genome}_vs_${bam2genome}"


# Make some directories to hold the final output for this sample.
mkdir -p "${sampleOutdir}/log"

# Make a directory structure to hold files passed between jobs.
# It is built within ${sampleOutdir} so it can be persistent.
# It gets removed at the end of the merge_bamintersect.sh script.
INTERMEDIATEDIR="${sampleOutdir}/tmp.${sample_name}"   # INTERMEDIATEDIR contains no trailing slash.

# Subdirectories will hold output from bamintersect.py
mkdir -p ${INTERMEDIATEDIR}/bamintersectPyOut

# Make a directory for all the small sorted, chromosome bam files generated by sort_bamintersect.sh
mkdir -p ${INTERMEDIATEDIR}/sorted_bams


echo "Running on $HOSTNAME. Using $TMPDIR as tmp"


if [ ! -s ${bam1} ]; then
    echo "ERROR: Could not find bam1 file ${bam1}"
    exit 1
fi

if [ ! -s ${bam2} ]; then
    echo "ERROR: Could not find bam2 file ${bam2}"
    exit 1
fi


################################################################################################
# Make a two-column chromosome list for the sort array jobs, where first column is short job name and second is the list of chromsomes for samtools view
echo
echo "sort_bamintersect"

#Outputs array job inputfile to stdout by collapsing chromosomes if too many
#Output is tab-delimited file with two columns: short name, list of chromosomes separated by spaces
function make_chrom_inputfile {
    local bamfile=$1
    local prefix=$2
    
    samtools idxstats ${bamfile} | awk -F "\t" 'BEGIN {OFS="\t"} $1!="*"' > ${TMPDIR}/idxstats.txt
    n=$(wc -l < ${TMPDIR}/idxstats.txt)
    
    #Chunk bam file based on >5M read chunks of chromosomes to reduce array job size for large runs
    cat ${TMPDIR}/idxstats.txt | awk -v maxchunksize=5000000 -v prefix=${prefix} -F "\t" 'BEGIN {OFS="\t"; chunknum=0; chroms=""; chunksize=0} {\
        chroms = chroms " " $1; \
        #$3+$4 adds mapped and unmapped reads \
        chunksize += $3+$4; \
        if(chunksize>maxchunksize) {\
            chunknum+=1;\
            print prefix chunknum, chroms;\
            chroms = ""
            chunksize=0;
        }
    }\
    END {if(chunksize>0) {chunknum+=1; print prefix chunknum, chroms}}'
}

make_chrom_inputfile ${bam1} "bam1_" > ${INTERMEDIATEDIR}/inputs.sort.bam1.txt
make_chrom_inputfile ${bam2} "bam2_" > ${INTERMEDIATEDIR}/inputs.sort.bam2.txt

if [ ! -s "${INTERMEDIATEDIR}/inputs.sort.bam1.txt" ]; then
    echo "No reads left in bam1; quitting successfully"
    exit 0
fi
if [ ! -s "${INTERMEDIATEDIR}/inputs.sort.bam2.txt" ]; then
    echo "No reads left in bam2; quitting successfully"
    exit 0
fi


num_lines=$(wc -l < ${INTERMEDIATEDIR}/inputs.sort.bam1.txt)
qsub -S /bin/bash -cwd -terse -j y -N sort_bamintersect_1.${sample_name} -o ${sampleOutdir}/log -t 1-${num_lines} --qos normal "${src}/sort_bamintersect.sh ${INTERMEDIATEDIR} ${sampleOutdir} ${sample_name} ${bam1} 1 ${bam1_keep_flags} ${bam1_exclude_flags}" > ${sampleOutdir}/sgeid.sort_bamintersect_1.${sample_name}

num_lines=$(wc -l < ${INTERMEDIATEDIR}/inputs.sort.bam2.txt)
qsub -S /bin/bash -cwd -terse -j y -N sort_bamintersect_2.${sample_name} -o ${sampleOutdir}/log -t 1-${num_lines} --qos normal "${src}/sort_bamintersect.sh ${INTERMEDIATEDIR} ${sampleOutdir} ${sample_name} ${bam2} 2 ${bam2_keep_flags} ${bam2_exclude_flags}" > ${sampleOutdir}/sgeid.sort_bamintersect_2.${sample_name}


################################################################################################
# Call bamintersect.py on each pairwise chrom combination via the bamintersect.sh array job:
echo
echo "bamintersect"
for bam1chromname in `cut -f1 ${INTERMEDIATEDIR}/inputs.sort.bam1.txt`; do
    for bam2chromname in `cut -f1 ${INTERMEDIATEDIR}/inputs.sort.bam2.txt`; do
        echo -e "${INTERMEDIATEDIR}/sorted_bams/${sample_name}.${bam1chromname}.bam\t${INTERMEDIATEDIR}/sorted_bams/${sample_name}.${bam2chromname}.bam\t${INTERMEDIATEDIR}/bamintersectPyOut/${bam1chromname}_vs_${bam2chromname}.bed" >> ${INTERMEDIATEDIR}/inputs.bamintersect.txt
    done
done


n=$(wc -l < ${INTERMEDIATEDIR}/inputs.bamintersect.txt)
qsub -S /bin/bash -cwd -terse -j y -hold_jid `cat ${sampleOutdir}/sgeid.sort_bamintersect_1.${sample_name} ${sampleOutdir}/sgeid.sort_bamintersect_2.${sample_name} | perl -pe 's/\n/,/g;'` -N bamintersect.${sample_name} -o ${sampleOutdir}/log -t 1-${n} "${src}/bamintersect.sh ${src} ${INTERMEDIATEDIR} ${ReqFullyAligned} ${max_mismatches} ${make_bed} ${reads_match}" > ${sampleOutdir}/sgeid.merge_bamintersect.${sample_name}
rm -f ${sampleOutdir}/sgeid.sort_bamintersect_1.${sample_name} ${sampleOutdir}/sgeid.sort_bamintersect_2.${sample_name}


################################################################################################
# Make the filter files for merge_bamintersect.sh and filter_tsv.sh.orig
if [ "${integrationsite}" != "null" ]; then
    # In the next line we rely on HA1 being the 5p side HA.
    IFS='_' read integrationSiteName HA1 HA2 <<< "${integrationsite}"
    echo "Looking up HA coordinates for ${integrationsite}:"
    HA_file="/vol/cegs/sequences/${bam1genome}/${integrationSiteName}/${integrationSiteName}_HomologyArms.bed"
    
    if [ -s "${HA_file}" ]; then
        if grep -q "${integrationSiteName}_${HA1}$" ${HA_file} && grep -q "${integrationSiteName}_${HA2}$" ${HA_file}; then
            # Get the HA coordinates:
            grep "${integrationSiteName}_${HA1}$" ${HA_file} > "${INTERMEDIATEDIR}/HA_coords.bed"
            grep "${integrationSiteName}_${HA2}$" ${HA_file} >> "${INTERMEDIATEDIR}/HA_coords.bed"
            
            echo -n "Found homology arm coordinates:"
            cat ${INTERMEDIATEDIR}/HA_coords.bed
            
            # Get the deletion range coordinates:
            cat ${INTERMEDIATEDIR}/HA_coords.bed | awk -F "\t" 'BEGIN {OFS="\t"} NR==1 {HA1_3p=$3} NR==2 {print $1, HA1_3p, $2}' > ${INTERMEDIATEDIR}/deletion_range.bed
            echo -n "Found coordinates for deletion spanning (bp): "
            awk -F "\t" 'BEGIN {OFS="\t"} {print $3-$2}' ${INTERMEDIATEDIR}/deletion_range.bed
        else
            echo "WARNING could not find homology arm coordinates for ${integrationSiteName}_${HA1} or ${integrationSiteName}_${HA2} in ${HA_file}"
        fi
    else
        echo "WARNING No HA file exists for integrationsite: ${integrationsite} so there is no Deletion Range"
        echo ""
        
        touch "${INTERMEDIATEDIR}/HA_coords.bed"
        touch "${INTERMEDIATEDIR}/deletion_range.bed"
    fi
else
    echo "No integration site was provided, and so there is no Deletion Range."
    echo ""
    
    touch "${INTERMEDIATEDIR}/HA_coords.bed"
    touch "${INTERMEDIATEDIR}/deletion_range.bed"
fi


################################################################################################
## Merge output from the array jobs.
echo
echo "merge_bamintersect"

export_vars="sampleOutdir=${sampleOutdir}"
export_vars="${export_vars},src=${src}"
export_vars="${export_vars},sample_name=${sample_name}"
export_vars="${export_vars},bam1genome=${bam1genome}"
export_vars="${export_vars},bam2genome=${bam2genome}"
export_vars="${export_vars},make_bed=${make_bed}"
export_vars="${export_vars},make_table=${make_table}"
export_vars="${export_vars},INTERMEDIATEDIR=${INTERMEDIATEDIR}"
export_vars="${export_vars},bam1=${bam1}"
export_vars="${export_vars},bam2=${bam2}"
export_vars="${export_vars},verbose=${verbose}"
export_vars="${export_vars},ReqFullyAligned=${ReqFullyAligned}"

qsub -S /bin/bash -cwd -terse -j y -hold_jid `cat ${sampleOutdir}/sgeid.merge_bamintersect.${sample_name}` --export=ALL,${export_vars} -N merge_bamintersect.${sample_name} -o ${sampleOutdir}/log "${src}/merge_bamintersect.sh" > /dev/null
rm -f ${sampleOutdir}/sgeid.merge_bamintersect.${sample_name}

echo
echo "Done!!!"
date

