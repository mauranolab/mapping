#!/bin/bash
set -eu -o pipefail

################################################
# Start of "getopt" section
#
# One colon=option has a required argument
# Two colons=option has an optional argument.
# This is not the same as making the option itself required/optional !
# Also, optional arguments MUST be passed with an equal sign.
#    This works:            option=argument
#    This does not work:    option argument
# The "optional" feature is a GNU extension to the basic getopt.

long_arg_list=(
    src:
    outdir:
    sample_name:
    bamname1:
    bam1_keep_flags:
    bam1_exclude_flags:
    bamname2:
    bam2_keep_flags:
    bam2_exclude_flags:
    reads_match:
    final_csv_delete::
    bam1_5p_HA::
    bam1_3p_HA::
    make_csv:
    make_table:
    help
)

long_args=$(printf "%s," "${long_arg_list[@]}")   # Turn the arg array into a string with commas.
long_args=$(echo ${long_args} | sed 's/,$/ /')    # Get rid of final comma.

# There must be at least one option associated with "-o", or you need to deal with getopt's default behavior.
CMD_LINE=$(getopt -o h --long "${long_args}" -n "$(basename "$0")" -- "$@")
# Catch getopt errors here if not using "set -e"

eval set -- "$CMD_LINE"

# Extract options and their arguments into variables.
while true ; do
    case "$1" in
        --src) 
            src=$2 ; shift 2 ;;
        --outdir)
            FINAL_OUTDIR=$2 ; shift 2 ;;
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
            reads_match=$2 ; shift 2 ;;
        --final_csv_delete)
            case "$2" in
                 "") final_csv_delete="NA" ; shift 2 ;;
                 *)  final_csv_delete=$2 ; shift 2 ;;
            esac ;;
        --bam1_5p_HA)
            case "$2" in
                 "") bam1_5p_HA="NA" ; shift 2 ;;
                 *)  bam1_5p_HA=$2 ; shift 2 ;;
            esac ;;
        --bam1_3p_HA)
            case "$2" in
                 "") bam1_3p_HA="NA" ; shift 2 ;;
                 *)  bam1_3p_HA=$2 ; shift 2 ;;
            esac ;;
        --make_csv)
            make_csv=$2 ; shift 2 ;;
        --make_table)
            make_table=$2 ; shift 2 ;;
        -h|--help) echo "Some helpful text." ; shift ; exit 0 ;;
        --) shift ; break ;;
        *) echo "getopt internal error!" ; exit 1 ;;
    esac
done

# End of getopt section.
################################################

# Make the master output directory for this analysis
rm -rf ${FINAL_OUTDIR}  # Clear out the old one, if necessary.
mkdir ${FINAL_OUTDIR}

# Make a place for the related slurm reports.
mkdir "${FINAL_OUTDIR}/log"

# Make a place for all the small chromosome bam files.
mkdir "${FINAL_OUTDIR}/bams"

export_vars="src=${src}"
export_vars="${export_vars},FINAL_OUTDIR=${FINAL_OUTDIR}"
export_vars="${export_vars},sample_name=${sample_name}"
export_vars="${export_vars},bamname1=${bamname1}"
export_vars="${export_vars},bam1_keep_flags=${bam1_keep_flags}"
export_vars="${export_vars},bam1_exclude_flags=${bam1_exclude_flags}"
export_vars="${export_vars},bamname2=${bamname2}"
export_vars="${export_vars},bam2_keep_flags=${bam2_keep_flags}"
export_vars="${export_vars},bam2_exclude_flags=${bam2_exclude_flags}"
export_vars="${export_vars},reads_match=${reads_match}"
export_vars="${export_vars},final_csv_delete=${final_csv_delete}"
export_vars="${export_vars},bam1_5p_HA=${bam1_5p_HA}"
export_vars="${export_vars},bam1_3p_HA=${bam1_3p_HA}"
export_vars="${export_vars},make_csv=${make_csv}"
export_vars="${export_vars},make_table=${make_table}"

sbatch --export=${export_vars} --output="${FINAL_OUTDIR}/log/launch_bam_intersect.${sample_name}.o_%j" \
       "${src}/launch_bam_intersect.sbatch"

