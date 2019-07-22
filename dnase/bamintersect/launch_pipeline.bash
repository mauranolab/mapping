#!/bin/bash
set -eu -o pipefail

src="/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/src"

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

    ######################################################################################
    # These values in this section need to get set via something.
    # However, it doesn't have to be via sourcing a file.
    source "/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/input_library/${input_file}"
    # Set in the above source file:
    #    bamname1
    #    bam1_keep_flags
    #    bam1_exclude_flags
    #    bamname2
    #    bam2_keep_flags
    #    bam2_exclude_flags
    #    reads_match
    #    final_csv_delete
    #    bam1_5p_HA
    #    bam1_3p_HA
    #    make_csv
    #    make_table

    # Could also derive sample_name from the bam files, for example.
    sample_name=${input_file}
    outdir="/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/outdirs/${sample_name}"

    ######################################################################################

    # Notice the equal signs, which are needed for getopt!
    "${src}/call_launch_bam_intersect.bash" \
        --src ${src} \
        --outdir ${outdir} \
        --sample_name ${sample_name} \
        --bamname1 ${bamname1} \
        --bam1_keep_flags ${bam1_keep_flags} \
        --bam1_exclude_flags ${bam1_exclude_flags} \
        --bamname2 ${bamname2} \
        --bam2_keep_flags ${bam2_keep_flags} \
        --bam2_exclude_flags ${bam2_exclude_flags} \
        --reads_match ${reads_match} \
        --final_csv_delete=${final_csv_delete} \
        --bam1_5p_HA=${bam1_5p_HA} \
        --bam1_3p_HA=${bam1_3p_HA} \
        --make_csv ${make_csv} \
        --make_table ${make_table}

done < launch_pipeline.manifest

