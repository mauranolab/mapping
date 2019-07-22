#!/bin/bash

# An example of a single call to call_launch_bam_intersect.bash
# This might be embedded in a larger program.
# The input variables need to be set first, from some source.

src="/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/src"

sample_name="mESC_RnHoxa_Clone2_2DEbs-H3K27ac-BS02585A"
outdir="/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/outdirs/${sample_name}"

bamname1="/vol/cegs/mapped/FCH2YVKBGXB/chipseq/mESC_RnHoxa_Clone2_2DEbs-H3K27ac-BS02585A/mESC_RnHoxa_Clone2_2DEbs-H3K27ac-BS02585A.mm10.bam"
bam1_keep_flags="0"
bam1_exclude_flags="3076"

bamname2="/vol/cegs/mapped/FCH2YVKBGXB/chipseq/mESC_RnHoxa_Clone2_2DEbs-H3K27ac-BS02585A/mESC_RnHoxa_Clone2_2DEbs-H3K27ac-BS02585A.rn6.bam"
bam2_keep_flags="0"
bam2_exclude_flags="3076"

reads_match=True
final_csv_delete=""

bam1_5p_HA=""
bam1_3p_HA=""

make_csv=True
make_table=False

# Notice the equal signs!
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

