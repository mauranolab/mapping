# Transposon scripts 
This repository contains all the scripts for the Transposon pipeline.

## Structure of Transposon folder
All new flowcells receive a separate folder name aligned.FLOWCELL  
A submit script is created for each flowcell submitFLOWCELL.sh

## Submit jobs from new flowcell
The createFlowcellSubmit.py script is used to submit samples from the newly sequenced flowcell.  
The script requires the flowcell name and writes the files to stdout.  
It reads in the info.txt in /vol/mauranolab/flowcells/data/FLOWCELL,  
which is a copy of the flowcell from the googlesheet

```
./createFlowcellSubmit.py --Flowcell FLOWCELL >../../submitFLOWCELL.sh
```

The resulting script has a separate line for each sample.  
Sample types allowed are Transposon RNA, Transposon DNA, Transposon iPCR, HiC, 3C, DNase-seq, and ChIP-seq.  
Below are examples of Transposon RNA, Transposon DNA, and Transposon iPCR samples.  
```
cd aligned.FLOWCELL
qsub -S /bin/bash -j y -N submit.BS00770B-RDL_20171106_K562_T0053_Hsp68_Day10_iPCR -b y "/home/maagj01/scratch/transposon/src/submitIntegrations.sh BS00770B-RDL_20171106_K562_T0053_Hsp68_Day10_iPCR CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGGGATCTTTGTCCAAACTCATCGAGCTCGGGA 0 0 16 CTTCCGACTTCAACTGTA /vol/mauranolab/flowcells/fastq/FCHM73JBGX5/Project_Maurano/Sample_BS00770B/"
qsub -S /bin/bash -j y -N submit.BS00932A-NEM_20180206_K562_Gglo_Hem_2_RNA -b y "/home/maagj01/scratch/transposon/src/submitBCcounts.sh BS00932A-NEM_20180206_K562_Gglo_Hem_2_RNA CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGGGATCTTTGTCCAAACTCATCGAGCTCGGGA 1 9 16 R1 AGCTGCACAGCAACACCGAGCTGGGCATCGTGGAGTACCAGCACGCCTTCAAGACCCCCATCGCCTTCGCCAGATC /vol/mauranolab/flowcells/fastq/FCHM73JBGX5/Project_Maurano/Sample_BS00932A/"
qsub -S /bin/bash -j y -N submit.BS00879A-NEM_20180206_K562_Gglo_NT_DNA -b y "/home/maagj01/scratch/transposon/src/submitBCcounts.sh BS00879A-NEM_20180206_K562_Gglo_NT_DNA CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGGGATCTTTGTCCAAACTCATCGAGCTCGGGA 8 2 16 R2 AGCTGCACAGCAACACCGAGCTGGGCATCGTGGAGTACCAGCACGCCTTCAAGACCCCCATCGCCTTCGCCAGATC /vol/mauranolab/flowcells/fastq/FCHHJVNBGX5/Project_Maurano/Sample_BS00879A/"
```


## Quality control of flowcell
To generate a flowcellInfo directory for the flowcell, go into the aligned.FLOWCELL folder.  
Run the FlowcellInfo.sh script and specify the output.  
```
bash ../src/Flowcell_Info.sh ~/public_html/blog/flowcellInfo/
cd ..
```

This will create a folder with the FLOWCELL name.  
Within this folder will be different folders with different quality metrics
```
FlowcellSummary -- Summary statistics for the flowcell
EditDist -- Hamming distance from expected plasmid/primer sequence
Weblogos_raw -- Weblogos without trimming
Weblogos_processed -- Weblogos with trimming
Weblogos_Barcode -- Weblogos for barcodes
Weblogos_UMI -- Weblogos for UMI
QC_RNA_DNA -- Quaility assessment of RNA and DNA samples
BarcodeFreq -- Barcode frequencies
BarcodeFreq_UMI -- Barcode frequencies adjusted for unique molecular identifiers.
UMI_distribution -- Distributions of UMIs per barcode
SaturationCurve -- Saturation curves of number of UMI adjusted Barcodes
```

## All vs all overlap between samples
The script below runs an all vs all comparison of specified samples.  
This generates a tsv file with all overlaps, as well as a heatmap of the jaccard index between each sample comparison.  
The analysis is good to check for contaminations between samples.  
First create a txt file with all samples to compare preceded my FLOWCELL folder (can so far only be relative path)  
<br>

```
cat FCHKWNLBGX5Samples.txt
aligned.FCHKWNLBGX5/BS00953A-RDL_20180207_K562_T0073_CMVNull_DNA
aligned.FCHKWNLBGX5/BS00989A-RDL_20180207_K562_T0074_Hsp68Low_RNA
aligned.FCHKWNLBGX5/BS00980A-MSH_20180129_K562_T0075_CMVNull_DNA
aligned.FCHKWNLBGX5/BS00986A-RDL_20180207_K562_T0073_CMVLow_RNA
aligned.FCHKWNLBGX5/BS00956A-RDL_20180207_K562_T0074_Hsp68Null_DNA
```

It can also use a comma separated list of all samples
```
aligned.FCHKWNLBGX5/BS00953A-RDL_20180207_K562_T0073_CMVNull_DNA,aligned.FCHKWNLBGX5/BS00989A-RDL_20180207_K562_T0074_Hsp68Low_RNA,aligned.FCHKWNLBGX5/BS00980A-MSH_20180129_K562_T0075_CMVNull_DNA,aligned.FCHKWNLBGX5/BS00986A-RDL_20180207_K562_T0073_CMVLow_RNA,aligned.FCHKWNLBGX5/BS00956A-RDL_20180207_K562_T0074_Hsp68Null_DNA
```

To run the script use the following command in the Transposon directory
```
./src/OverlapAnySampleWithRest.R -s FCHKWNLBGX5Samples.txt -f -t 1 -o ~/public_html/blog/20180423/FCHKWNLBGX5Overlap
-s -- Samplefile
-f -- Samples are in a file (Do not use if comma separated file)
-t -- Threshold over UMI corrected barcodes to compare overlap with
-o -- Output name (createas a OUTPUT.tsv and OUTPUT.pdf)
```

## Overlap between specific samples
This is more used when we want to compare overlap of DNA, RNA, and iPCR samples from the same transfection.  
Here we can select the threshold for UMI corrected barcodes for DNA and RNA samples, while iPCR remains unfiltered.

```
./src/OverlapSamples.R \
\
-a aligned.FCHFMK5BGX5/BS00726E-RDL_20171211_K562_LTR7_day1_DNA,\
aligned.FCHFMK5BGX5/BS00727D-RDL_20171211_K562_CMV_day1_DNA \
\
-b aligned.FCHFMK5BGX5/BS00743A-RDL_20171211_K562_LTR7_day1_RNA,\
aligned.FCHFMK5BGX5/BS00744A-RDL_20171211_K562_CMV_day1_RNA \
\
-c aligned.FCHFNJJBGX5/BS00772A-RDL_20171211_K562_T0055_LTR7_d1sort_iPCR,\
aligned.FCHFNJJBGX5/BS00773A-RDL_20171211_K562_T0056_CMV_d1sort_iPCR \
-t 10 \
-o ~/public_html/blog/20180430/Ovarlap_10.tsv 
```

This results in a output file
```
A	B	C	#A	#B	#C_mapped	#A_10	#B_10	#C_1	#Intersect_A10_B10	#Intersect_A10_C1	#Intersect_B10_C1	#Intersect_A10_B10_C1
BS00726E-RDL-20171211-K562-LTR7-day1-DNA	BS00743A-RDL-20171211-K562-LTR7-day1-RNA	BS00772A-RDL-20171211-K562-T0055-LTR7-d1sort-iPCR	54,104	5,655	4,079	43,864,658	4,079	2,160	2,913	249	231
BS00727D-RDL-20171211-K562-CMV-day1-DNA	BS00744A-RDL-20171211-K562-CMV-day1-RNA	BS00773A-RDL-20171211-K562-T0056-CMV-d1sort-iPCR	14,422	1,785	125	6,032	896	125	130	24	32	18
```

Options
```
Usage: src/OverlapSamples.R [options]


Options:
        -a CHARACTER, --samplesA=CHARACTER
                Comma seperated list of SamplesA

        -b CHARACTER, --samplesB=CHARACTER
                Comma seperated list of SamplesB

        -c CHARACTER, --samplesC=CHARACTER
                Comma seperated list of SamplesC (assumed to be iPCR)

        -d CHARACTER, --samplesC2=CHARACTER
                Comma seperated list of SamplesC (assumed to be iPCR)

        -t CHARACTER, --threshold=CHARACTER
                Comma seperated list of SamplesB

        -o CHARACTER, --output=CHARACTER
                output folder for all comparisons

        -h, --help
                Show this help message and exit
```



## Run transposon analysis
To analysis the transposon samples, we need to submit DNA, RNA, and iPCR samples as follows.

```
./src/AnalyseOverlappedBCs.sh  aligned.Merged/BS00829-K562_T0058_CMV_d5sort_High_DNA_Merged/BS00829-K562_T0058_CMV_d5sort_High_DNA_Merged.barcode.counts.UMI.corrected.txt \
aligned.Merged/BS00746-K562_CMV_day5_GFPHigh_RNA_Merged/BS00746-K562_CMV_day5_GFPHigh_RNA_Merged.barcode.counts.UMI.corrected.txt \
aligned.FCHGLHJBGX5/BS00846A-RDL_20171211_K562_T0058_CMV_d5sort_High_iPCR/BS00846A-RDL_20171211_K562_T0058_CMV_d5sort_High_iPCR.barcodes.coords.bed \
~/public_html/blog/GFP_sort/CMV_High_UMI

```

