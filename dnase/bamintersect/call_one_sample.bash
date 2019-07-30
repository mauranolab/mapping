#!/bin/bash

sample_name=H1_dHPRT1_LP123_Capture-BS02831A_NEW

bamname1=/vol/cegs/mapped/FCHC7T7BGXB/capture/H1_dHPRT1_LP123_Capture-BS02831A/H1_dHPRT1_LP123_Capture-BS02831A.cegsvectors_LP123.bam

bamname2=/vol/cegs/mapped/FCHC7T7BGXB/capture/H1_dHPRT1_LP123_Capture-BS02831A/H1_dHPRT1_LP123_Capture-BS02831A.hg38_full.bam


/vol/mauranolab/cadlej01/projects/Sud_mm10_rn6/src/call_launch_bam_intersect.bash \
       --sample_name ${sample_name} \
       --bamname1 ${bamname1} \
       --bamname2 ${bamname2}

