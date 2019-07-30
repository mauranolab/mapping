#!/bin/bash

# An example of a single call to call_launch_bam_intersect.bash
# This might be embedded in a larger program.
# The input variables need to be set first, from some source.

bamname1="/vol/cegs/mapped/FCH2YVKBGXB/chipseq/mESC_RnHoxa_Clone2_2DEbs-H3K27ac-BS02585A/mESC_RnHoxa_Clone2_2DEbs-H3K27ac-BS02585A.mm10.bam"
bam1_keep_flags="0"   # Different from the default.

bamname2="/vol/cegs/mapped/FCH2YVKBGXB/chipseq/mESC_RnHoxa_Clone2_2DEbs-H3K27ac-BS02585A/mESC_RnHoxa_Clone2_2DEbs-H3K27ac-BS02585A.rn6.bam"
bam2_keep_flags="0"   # Different from the default.

# Construct the standard sample name:
sample_name="mESC_RnHoxa_Clone2_2DEbs-H3K27ac-BS02585A"

echo
echo "Launching: ${sample_name}"

"/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/src/call_launch_bam_intersect.bash" \
    --sample_name ${sample_name} \
    --bamname1 ${bamname1} \
    --bam1_keep_flags ${bam1_keep_flags} \
    --bamname2 ${bamname2} \
    --bam2_keep_flags ${bam2_keep_flags} \
    --reads_match \
    --do_not_make_table

