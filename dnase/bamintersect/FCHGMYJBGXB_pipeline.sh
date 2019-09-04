#!/bin/bash

./submit_bamintersect.sh \
       --sample_name BI6_Cast_dSox2_LP131b_Clone1-BS03022A \
       --bam1 /vol/cegs/mapped/FCHGMYJBGXB/capture/BI6_Cast_dSox2_LP131b_Clone1-BS03022A/BI6_Cast_dSox2_LP131b_Clone1-BS03022A.cegsvectors_LP131.bam \
       --bam1genome LP131b \
       --bam2 /vol/cegs/mapped/FCHGMYJBGXB/capture/BI6_Cast_dSox2_LP131b_Clone1-BS03022A/BI6_Cast_dSox2_LP131b_Clone1-BS03022A.mm10.bam \
       --bam2genome mm10 \
       --clear_logs

exit 0

./submit_bamintersect.sh \
       --sample_name BI6_Cast_dSox2_LP131b_Clone1-BS03022A \
       --bam1 /vol/cegs/mapped/FCHGMYJBGXB/capture/BI6_Cast_dSox2_LP131b_Clone1-BS03022A/BI6_Cast_dSox2_LP131b_Clone1-BS03022A.cegsvectors_pSpCas9.bam \
       --bam1genome pSpCas9 \
       --bam2 /vol/cegs/mapped/FCHGMYJBGXB/capture/BI6_Cast_dSox2_LP131b_Clone1-BS03022A/BI6_Cast_dSox2_LP131b_Clone1-BS03022A.mm10.bam \
       --bam2genome mm10 \
       --clear_logs

./submit_bamintersect.sh \
       --sample_name BI6_Cast_dPiga-BS03020A \
       --bam1 /vol/cegs/mapped/FCHGMYJBGXB/capture/BI6_Cast_dPiga-BS03020A/BI6_Cast_dPiga-BS03020A.cegsvectors_pSpCas9.bam \
       --bam1genome pSpCas9 \
       --bam2 /vol/cegs/mapped/FCHGMYJBGXB/capture/BI6_Cast_dPiga-BS03020A/BI6_Cast_dPiga-BS03020A.mm10.bam \
       --bam2genome mm10 \
       --clear_logs

./submit_bamintersect.sh \
       --sample_name BI6_Cast_dSox2_LP97b_Clone1-BS03021A \
       --bam1 /vol/cegs/mapped/FCHGMYJBGXB/capture/BI6_Cast_dSox2_LP97b_Clone1-BS03021A/BI6_Cast_dSox2_LP97b_Clone1-BS03021A.cegsvectors_LP097.bam \
       --bam1genome LP097b \
       --bam2 /vol/cegs/mapped/FCHGMYJBGXB/capture/BI6_Cast_dSox2_LP97b_Clone1-BS03021A/BI6_Cast_dSox2_LP97b_Clone1-BS03021A.mm10.bam \
       --bam2genome mm10 \
       --clear_logs

./submit_bamintersect.sh \
       --sample_name BI6_Cast_dSox2_LP97b_Clone1-BS03021A \
       --bam1 /vol/cegs/mapped/FCHGMYJBGXB/capture/BI6_Cast_dSox2_LP97b_Clone1-BS03021A/BI6_Cast_dSox2_LP97b_Clone1-BS03021A.cegsvectors_pSpCas9.bam \
       --bam1genome pSpCas9 \
       --bam2 /vol/cegs/mapped/FCHGMYJBGXB/capture/BI6_Cast_dSox2_LP97b_Clone1-BS03021A/BI6_Cast_dSox2_LP97b_Clone1-BS03021A.mm10.bam \
       --bam2genome mm10 \
       --clear_logs

