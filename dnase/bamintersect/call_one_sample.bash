#!/bin/bash

sample_name=H1_dHPRT1_LP123_Capture-BS02831A_NEW

bam1=/vol/cegs/mapped/FCHC7T7BGXB/capture/H1_dHPRT1_LP123_Capture-BS02831A/H1_dHPRT1_LP123_Capture-BS02831A.cegsvectors_LP123.bam

x=${bam1%.bam}
y=${x##*.}
bam1genome=${y##*_} 


bam2=/vol/cegs/mapped/FCHC7T7BGXB/capture/H1_dHPRT1_LP123_Capture-BS02831A/H1_dHPRT1_LP123_Capture-BS02831A.hg38_full.bam
bam2genome=hg38


/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/submit_bamintersect.sh -h
exit 0

# The "--clear_logs" below is optional.
/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/submit_bamintersect.sh \
       --sample_name ${sample_name} \
       --bam1 ${bam1} \
       --bam1genome ${bam1genome} \
       --bam2 ${bam2} \
       --bam2genome ${bam2genome} \
       --clear_logs

