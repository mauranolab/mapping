#!/bin/env python3.5
import re
import argparse
import pandas as pd

parser = argparse.ArgumentParser(prog = "createFlowcellSubmit", description = "Outputs submit commands for all samples from the info.txt located in the flowcell data folder", add_help=True)
parser.add_argument('--Flowcell', action='store', help='Flowcell name')

try:
       args = parser.parse_args()
except argparse.ArgumentError as exc:
       print(exc.message, '\n', exc.argument)
       sys.exit(2)



#Transposon pipeline
def transSamples(sampleID, sampleName, subLibrary, lab, sampleType, r1Trim, r2Trim):
       fullName = sampleID +  subLibrary + "-" + sampleName
       
       
       fileLocation = '/vol/mauranolab/flowcells/fastq/' + Flowcell + "/Project_" + lab + "/Sample_" + sampleID + subLibrary + '/"'
       
       qsub = 'qsub -S /bin/bash -j y -N submit.' + fullName + ' -b y "/home/maagj01/scratch/transposon/src/'
       Primer = ' CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGGGATCTTTGTCCAAACTCATCGAGCTCGGGA ' + str(int(round(float(r1Trim)))) + ' ' + str(int(round(float(r2Trim))))
       
       #Unique line to each library 
       if sampleType == "Transposon DNA":
              Script= 'submitBCcounts.sh '
              sampleUnique = ' 16 R2 AGCTGCACAGCAACACCGAGCTGGGCATCGTGGAGTACCAGCACGCCTTCAAGACCCCCATCGCCTTCGCCAGATC '
              
       if sampleType == "Transposon RNA":
              Script = 'submitBCcounts.sh '
              sampleUnique = ' 16 R1 AGCTGCACAGCAACACCGAGCTGGGCATCGTGGAGTACCAGCACGCCTTCAAGACCCCCATCGCCTTCGCCAGATC '
              
       if sampleType == "Transposon iPCR":
              Script = 'submitIntegrations.sh '
              sampleUnique = ' 16 CTTCCGACTTCAACTGTA '
              
       
       submitCommand = qsub + Script + fullName + Primer + sampleUnique + fileLocation
       return(submitCommand)


#Chromosome configuration capture
#TODO fix organism
def chromConfCapture(sampleID, sampleName, subLibrary, lab):
       organism = "hg38" 
       fullName = sampleID +  subLibrary + "-" + sampleName
       
       qsub = 'qsub -S /bin/bash -j y -N submit.' + fullName + \
       ' -b y "/home/maagj01/scratch/transposon/src/submitHiC.sh ' +fullName 
       
       sampleUnique = '/home/maagj01/scratch/transposon/captureC/config/config_Dpn.txt  /home/maagj01/scratch/transposon/captureC/genomeFrag/hg38_dpnii.bed  K562_Dpn ' + organism
       
       fileLocation = ' /vol/mauranolab/flowcells/fastq/' + Flowcell + "/Project_" + lab + "/Sample_" + sampleID + subLibrary + '/"'
       
       submitCommand = qsub + sampleUnique + fileLocation
       return(submitCommand)


#DNase/ ChIP-seq
#TODO fix organism for DNase and HiC
def DNase(sampleID, sampleName, subLibrary, lab):
       organism = "hg38" 
       fullName = sampleID +  subLibrary + "-" + sampleName
       submitCommand = '/vol/isg/encode/dnase/src/submit.sh ' + organism + " " + sampleName + " " + sampleID + subLibrary
        
       return(submitCommand)


#Submit right library
def annotateSample(sampleID, sampleName, sublibrary, lab, sampleType, r1Trim, r2Trim):
       if sampleType == "Transposon DNA" or sampleType == "Transposon RNA" or sampleType == "Transposon iPCR":
              qsubLine = transSamples(sampleID, sampleName, sublibrary, lab, sampleType, r1Trim, r2Trim)
       elif sampleType == "Hi-C" or sampleType == "3C-seq" or sampleType == "Capture-C":
              qsubLine = chromConfCapture(sampleID, sampleName, sublibrary, lab)
       elif sampleType == "DNase" or sampleType == "Nano-DNase" or sampleType == "ChIP-seq" or sampleType == "DNase-seq":
              qsubLine = DNase(sampleID, sampleName, sublibrary, lab)
       return(qsubLine)



####
#Read in flowcell and create a submit script for all samples with BS number
####
Flowcell = args.Flowcell

flow = open('/vol/mauranolab/flowcells/data/' + Flowcell + '/info.txt', 'r')
#flow = open('/home/maagj01/scratch/transposon/Analysis/Nick_Mamrak/RnaAsDna/march12_RnaAsDNA.tsv', 'r')
flow = flow.readlines()
flow = [line.rstrip().split('\t') for line in flow]
flowcellFile = pd.DataFrame(flow[20:], columns=flow[20]).iloc[1:]



BS = re.compile(r'BS[0-9]{5}')
for index, Sample in flowcellFile.iterrows():
       try:
              BSnum = BS.search(str(Sample["Sample #"]))
              if BSnum is not None:
                     qsubLine = annotateSample(Sample["Sample #"], Sample["#Sample Name"], Sample["Sub-library"], Sample["Lab"], Sample["Sample Type"], str(Sample["R1 Trim (P5)"]), str(Sample["R2 Trim (P7)"]))
                     print(qsubLine)
       except:
              print('ERROR in ', Sample["#Sample Name"])


