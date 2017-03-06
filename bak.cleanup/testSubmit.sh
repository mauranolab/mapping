#!/bin/bash
shopt -s expand_aliases

#BUGBUG weblogo breaking pipefail
set -e # -o pipefail


#sample=BS00105A-Bglo_K562d5_Max_DNA
#ds=`echo $sample | cut -d '-' -f1`
#name=`echo $sample | cut -d '-' -f2`
#echo "$ds $name"
#qsub -S /bin/bash -j y -N submit.${ds}-${name} -b y "~/scratch/transposon/src/submitBCcounts.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGG 4 0 16 R2 /vol/mauranolab/flowcells/fastq/FCHWG37BGXY/Project_Lab/Sample_${ds}/"


sample=BS00106A-Bglo_K562d5_JM109_DNA
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
qsub -S /bin/bash -j y -N submit.${ds}-${name} -b y "~/scratch/transposon/src/submitBCcounts.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGG 5 1 16 R2 /vol/mauranolab/flowcells/fastq/FCHWG37BGXY/Project_Lab/Sample_${ds}/"

sample=BS00104A-BGloGFPBC4_JM109_Plasmid
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
bqsub -j y -N submit.${ds}-${name} -b y "~/scratch/transposon/src/submitBCcounts.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGG 5 1 16 R2 /vol/mauranolab/flowcells/fastq/FCHWG37BGXY/Project_Lab/Sample_${ds}/"


sample=BS00073B-5xBGlo_K562d6_RNA
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
bqsub -j y -N submit.${ds}-${name} -b y "~/scratch/transposon/src/submitBCcounts.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGG 0 0 16 R1 /vol/mauranolab/flowcells/fastq/FCHWG37BGXY/Project_Lab/Sample_${ds}/"


sample=BS00103A-BGloGFPBC4_MAX_Plasmid
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
bqsub -j y -N submit.${ds}-${name} -b y "~/scratch/transposon/src/submitBCcounts.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGG 4 0 16 R2 /vol/mauranolab/flowcells/fastq/FCHWG37BGXY/Project_Lab/Sample_${ds}/"
