#!/bin/bash
set -eu -o pipefail

SOURCE_base=/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/input_library
FINAL_base=/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/outdirs
src=/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/src

while read -r input_file; do
    # Skip blank lines and comments in the input manifest.
    input_file_len=${#input_file}
    if [ ${input_file_len} -lt 1 ]; then
        # echo "Skipping blank line."
        continue
    else
        first_chr=${input_file:0:1}
        if [ ${first_chr} = "#" ] ; then
            # echo "Skipping: ${input_file}"
            continue
        fi
    fi

    SOURCE_FILE="${SOURCE_base}/${input_file}"
    FINAL_OUTDIR="${FINAL_base}/${input_file}"
    echo
    echo $SOURCE_FILE
    echo $FINAL_OUTDIR

    export_vars="SOURCE_FILE=${SOURCE_FILE}"
    export_vars="${export_vars},FINAL_OUTDIR=${FINAL_OUTDIR}"
    export_vars="${export_vars},src=${src}"

    # Make the master output directory for this analysis
    rm -rf ${FINAL_OUTDIR}  # Clear out the old one, if necessary.
    mkdir ${FINAL_OUTDIR}

    # Make a place for the related slurm reports.
    mkdir "${FINAL_OUTDIR}/log"

    # Make a place for all the small chromosome bam files.
    mkdir "${FINAL_OUTDIR}/bams"

    sbatch --export=${export_vars} --output="${FINAL_OUTDIR}/log/launch_bam_intersect.${input_file}.o_%j"  "${src}/launch_bam_intersect.sbatch"
done < launch_pipeline.manifest

