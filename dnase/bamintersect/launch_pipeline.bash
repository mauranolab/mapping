#!/bin/bash
set -eu -o pipefail

while read -r input_file; do
    echo

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

    # sam format flag settings. Doing it this way to show usage. They aren't needed if the defaults are desired.

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
    
    # Some flags. Doing it this way to show usage. They aren't needed if the defaults are desired.

        reads_match=""    # This is the default.
        # Could also be: reads_match="--reads_match"
        # This parameter tells bam_intersect.py how to match the reads from each bam file.
        # If set, reads should be matched like [read1 read1], or [read2 read2].
        # If unset, reads are paired like [read1 read2], or [read2 read1].
        # Set this option for unpaired reads.
        # The default is to have it unset.

        do_not_make_csv=""   # This is the default.
        # Could also be: do_not_make_csv="--do_not_make_csv"
        # Make a single csv file (the default), or 2 bam files.

        do_not_make_table=""  # This is the default.
        # Could also be: do_not_make_table="--do_not_make_table"
        # Make the LP Integration read count table (the default), or not.

    source "/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/input_library/${input_file}"
    # Set in the above source file:
    #    bamname1
    #    bamname2
    #    plus anything from the above that needs to be overwritten (none in my test case).
    # The bamnames are full path names.

    ##########################################################################################
    # Construct the standard sample name. Could also just get it from somewhere, if available.
    sample_name=${input_file}

    x=${sample_name##*.}   # Assign the base to x. This base will be some descriptor made to differentiate
                           # input files with the same standard sample name, but different cegsvector IDs.
                           # Example:
                           #     H1_dHPRT1_LP123_Capture-BS02831A.pSpCas9_hg38
                           #     H1_dHPRT1_LP123_Capture-BS02831A.LP123_hg38
    sample_name=${sample_name%.${x}}   # Cut off the trailing ${x}

    echo "Launching: ${sample_name}"

    ##########################################################################################

    # The last three flags aren't needed when using the defaults. Just showing usage here.
    "/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/src/call_launch_bam_intersect.bash" \
        --sample_name ${sample_name} \
        --bamname1 ${bamname1} \
        --bam1_keep_flags ${bam1_keep_flags} \
        --bam1_exclude_flags ${bam1_exclude_flags} \
        --bamname2 ${bamname2} \
        --bam2_keep_flags ${bam2_keep_flags} \
        --bam2_exclude_flags ${bam2_exclude_flags} \
        ${reads_match} \
        ${do_not_make_csv} \
        ${do_not_make_table}

done < launch_pipeline.manifest


# What a minimal version of the call looks like, if all the defaults are accepted:
#    "/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/src/call_launch_bam_intersect.bash" \
#        --sample_name ${sample_name} \
#        --bamname1 ${bamname1} \
#        --bamname2 ${bamname2} 

