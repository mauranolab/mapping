#!/bin/bash

#BUGBUG weblogo breaking pipefail
set -e # -o pipefail

#FCHWLJCBGXY_2016nov04
#iPCR

#iPCR
#full seq is CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGGGATCTTTGTCCAAACTCATCGAGCTCGGGA
sample=BS00067A-5xBGlo_K562d4_2hDpn_iPCR
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
qsub -S /bin/bash -j y -N submit.${ds}-${name} -b y "~/scratch/transposon/src/submitIntegrations.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGG 0 0 16 /vol/mauranolab/flowcells/fastq/FCHGMLVBGXY/Project_Lab/Sample_${ds}/"

sample=BS00068A-5xBGlo_K562d6_2hDpn_iPCR
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
qsub -S /bin/bash -j y -N submit.${ds}-${name} -b y "~/scratch/transposon/src/submitIntegrations.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGG 0 0 16 /vol/mauranolab/flowcells/fastq/FCHGMLVBGXY/Project_Lab/Sample_${ds}/"

sample=BS00069A-5xBGlo_K562d10_2hDpn_iPCR
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
qsub -S /bin/bash -j y -N submit.${ds}-${name} -b y "~/scratch/transposon/src/submitIntegrations.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGG 0 0 16 /vol/mauranolab/flowcells/fastq/FCHGMLVBGXY/Project_Lab/Sample_${ds}/"

sample=BS00070A-5xBGlo_K562d4_16hDpn_iPCR
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
qsub -S /bin/bash -j y -N submit.${ds}-${name} -b y "~/scratch/transposon/src/submitIntegrations.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGG 1 0 16 /vol/mauranolab/flowcells/fastq/FCHGMLVBGXY/Project_Lab/Sample_${ds}/"

sample=BS00071A-5xBGlo_K562d6_16hDpn_iPCR
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
qsub -S /bin/bash -j y -N submit.${ds}-${name} -b y "~/scratch/transposon/src/submitIntegrations.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGG 1 0 16 /vol/mauranolab/flowcells/fastq/FCHGMLVBGXY/Project_Lab/Sample_${ds}/"

sample=BS00072A-5xBGlo_K562d10_16hDpn_iPCR
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
qsub -S /bin/bash -j y -N submit.${ds}-${name} -b y "~/scratch/transposon/src/submitIntegrations.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGG 1 0 16 /vol/mauranolab/flowcells/fastq/FCHGMLVBGXY/Project_Lab/Sample_${ds}/"

sample=BS00075A-5xBGlo_K562d4_2x16hDpn_iPCR
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
qsub -S /bin/bash -j y -N submit.${ds}-${name} -b y "~/scratch/transposon/src/submitIntegrations.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGG 0 1 16 /vol/mauranolab/flowcells/fastq/FCHGMLVBGXY/Project_Lab/Sample_${ds}/"

sample=BS00076A-5xBGlo_K562d6_2x16hDpn_iPCR
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
qsub -S /bin/bash -j y -N submit.${ds}-${name} -b y "~/scratch/transposon/src/submitIntegrations.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGG 0 1 16 /vol/mauranolab/flowcells/fastq/FCHGMLVBGXY/Project_Lab/Sample_${ds}/"

sample=BS00077A-5xBGlo_K562d10_2x16hDpn_iPCR
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
qsub -S /bin/bash -j y -N submit.${ds}-${name} -b y "~/scratch/transposon/src/submitIntegrations.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGG 0 1 16 /vol/mauranolab/flowcells/fastq/FCHGMLVBGXY/Project_Lab/Sample_${ds}/"
