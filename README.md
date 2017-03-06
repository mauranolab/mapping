#Transposon scripts 
Submit as below for corresponding DNA, RNA, or iPCR library

DNA
```
sample=BS00NNN-SampleName_DNA
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
qsub -S /bin/bash -j y -N submit.${ds}-${name} -b y "./submitBCcounts.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGGGATC trimR1 trimR2 bclen R2 /vol/mauranolab/flowcells/fastq/FLOWCELL/Project_Lab/Sample_${ds}/"
```

iPCR
```
sample=BS00NNN-SampleName_iPCR
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
qsub -S /bin/bash -j y -N submit.${ds}-${name} -b y "./submitIntegrations.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGGGATC trimR1 trimR2 bclen /vol/mauranolab/flowcells/fastq/FLOWCELL/Project_Lab/Sample_${ds}/"
```

RNA
```
sample=BS00NNN-SampleName_RNA
ds=`echo $sample | cut -d '-' -f1`
name=`echo $sample | cut -d '-' -f2`
echo "$ds $name"
qsub -S /bin/bash -j y -N submit.${ds}-${name} -b y "bash ./submitBCcounts.sh ${ds}-${name} CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGGGATC trimR1 trimR2 bclen R1 /vol/mauranolab/flowcells/fastq/FLOWCELL/Project_Lab/Sample_${ds}/"
```
