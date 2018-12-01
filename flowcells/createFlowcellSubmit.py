#!/bin/env python3.5
import sys
import re
import argparse
import pandas as pd

version="1.1"

doDNaseCleanup = False

parser = argparse.ArgumentParser(prog = "createFlowcellSubmit", description = "Outputs submit commands for all samples from the info.txt located in the flowcell data folder", allow_abbrev=False)
parser.add_argument('--flowcellID', action='store', type = str, help='Flowcell ID (will look for /vol/mauranolab/flowcells/FCxxx/info.txt)', required=True)
parser.add_argument('--project', action='store', type = str, help='Only process samples in this project', required=False)
parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s ' + version)

basedir="/vol/mauranolab/flowcells"

try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument)
    sys.exit(2)



###Handlers for individual data types
#Transposon pipeline
def transposonSamples(sampleName, sampleID, lab, sampleType, species, r1Trim, r2Trim):
    fullSampleName = sampleID + "-" + sampleName
    fileLocation = basedir + "/fastq/" + flowcellID + "/Project_" + lab + "/Sample_" + sampleID + '/"'
    
    qsub = "qsub -S /bin/bash -j y -N submit." + fullSampleName + ' -b y "/vol/mauranolab/transposon/src/'
    
    #default parameters that can be overriden below
    BCreadSeq = "CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGGGATCTTTGTCCAAACTCATCGAGCTCGGGA"
    bclen = "16"
    
    #Unique line to each library 
    if sampleType == "Transposon DNA":
        program= "submitBCcounts.sh"
        bcread = "R2"
        plasmidSeq="AGCTGCACAGCAACACCGAGCTGGGCATCGTGGAGTACCAGCACGCCTTCAAGACCCCCATCGCCTTCGCCAGATC"
        extractBCargs="--no-BCrevcomp"
    elif sampleType == "Transposon RNA":
        program = "submitBCcounts.sh"
        bcread = "R1"
        plasmidSeq="AGCTGCACAGCAACACCGAGCTGGGCATCGTGGAGTACCAGCACGCCTTCAAGACCCCCATCGCCTTCGCCAGATC"
        extractBCargs="--no-BCrevcomp"
    elif sampleType == "Transposon iPCR":
        program = "submitIntegrations.sh"
        bcread = "R1"
        plasmidSeq="CTTCCGACTTCAACTGTA"
        extractBCargs="--no-BCrevcomp"
    elif sampleType == "Transposon iPCR Capture":
        program = "submitIntegrations.sh"
        BCreadSeq="AAAGATCCCACAACTAGCBBBBBBBBBBBBBBBBCCTAATTGCCTAGAAAGGAGCAGACG"
        bcread = "R2"
        plasmidSeq="None"
        extractBCargs="--BCrevcomp"
    else:
        raise Exception("Can't handle sampleType" + sampleType)
    
    submitCommand = ' '.join([qsub + program, fullSampleName, BCreadSeq, str(r1Trim), str(r2Trim), bclen, bcread, plasmidSeq, extractBCargs, fileLocation])
    return(submitCommand)


#Chromosome configuration capture
#TODO fix organism
def chromConfCapture(sampleName, sampleID, lab, sampleType, species):
    organism = "hg38" 
    fullSampleName = sampleID + "-" + sampleName
    
    qsub = "qsub -S /bin/bash -j y -N submit." + fullSampleName + \
    ' -b y "/home/maagj01/scratch/transposon/src/submitHiC.sh ' +fullSampleName 
    
    sampleUnique = "/home/maagj01/scratch/transposon/captureC/config/config_Dpn.txt  /home/maagj01/scratch/transposon/captureC/genomeFrag/hg38_dpnii.bed K562_Dpn " + organism
    
    fileLocation = basedir + "/fastq/" + flowcellID + "/Project_" + lab + "/Sample_" + sampleID + '/"'
    
    submitCommand = qsub + sampleUnique + " " + fileLocation
    return(submitCommand)


#DNase/ ChIP-seq
def DNase(sampleName, sampleID, lab, sampleType, species):
    speciesToGenomeReference = {
        'Human': 'hg38_noalt',
        'Mouse': 'mm10'
    }
    reference = speciesToGenomeReference[species]
    
    fullSampleName = sampleID + "-" + sampleName
    submitCommand = "/vol/mauranolab/mapped/src/submit.sh " + reference + " mapBwaAln,dnase " + sampleName + " " + sampleID
    
    global doDNaseCleanup
    doDNaseCleanup = True
    
    return(submitCommand)


def DNA(sampleName, sampleID, lab, sampleType, species):
    speciesToGenomeReference = {
        'Human': 'hg38_full',
        'Mouse': 'mm10',
        'Rat': 'rn6',
        'Human+mouse': 'hg38_full,mm10',
        'Human+yeast': 'hg38_sacCer3',
        'Mouse+yeast': 'mm10_sacCer3',
        'Rat+yeast': 'rn6_sacCer3'
    }
    reference = speciesToGenomeReference[species]
    fullSampleName = sampleID + "-" + sampleName
    submitCommand = "/vol/mauranolab/mapped/src/submit.sh " + reference + " mapBwaMem,callsnps " + sampleName + " " + sampleID
    
    global doDNaseCleanup
    doDNaseCleanup = True
    
    return(submitCommand)


###Dispatch appropriate function handler per sample line
#Takes dict representing a single row from the FC file iterator and returns the processing command line
def annotateSample(Sample, project):
    try:
        sampleName=Sample["#Sample Name"]
        sampleID=Sample["Sample #"]
        lab=Sample["Lab"]
        sampleType=Sample["Sample Type"]
        species=Sample["Species"]
        r1Trim=str(Sample["R1 Trim (P5)"])
        r2Trim=str(Sample["R2 Trim (P7)"])
        
        if project is not None and lab != project:
            return
        
        BS = re.compile(r'BS[0-9]{5}')
        BSnum = BS.search(sampleID)
        if BSnum is None:
            raise Exception("Can't parse " + sampleID + "as BS number")
        
        if sampleType != "Pool":
            if sampleType == "Transposon DNA" or sampleType == "Transposon RNA" or sampleType == "Transposon iPCR" or sampleType == "Transposon iPCR Capture":
                qsubLine = transposonSamples(sampleName, sampleID, lab, sampleType, species, r1Trim, r2Trim)
            elif sampleType == "Hi-C" or sampleType == "3C-seq" or sampleType == "Capture-C":
                qsubLine =  "#" + chromConfCapture(sampleName, sampleID, lab, sampleType, species)
            elif sampleType == "Nano-DNase" or sampleType == "ChIP-seq" or sampleType == "DNase-seq":
                qsubLine = DNase(sampleName, sampleID, lab, sampleType, species)
            elif sampleType == "DNA" or sampleType == "DNA Capture":
                qsubLine = DNA(sampleName, sampleID, lab, sampleType, species)
            else:
                raise Exception("Don't know how to process " + sampleType)
            
            #We have successfully generated the command line
            print(qsubLine)
    except Exception as e:
        print("WARNING for sample ", sampleName, ":", e, file=sys.stderr)
        print("#Couldn't process " + sampleName + "-" + sampleID)



####
#Read in flowcell and create a submit script for all samples with BS number
####
flowcellID = args.flowcellID
project = args.project

flowcellInfoFile = open(basedir + "/data/" + flowcellID + "/info.txt", 'r')
#flowcellInfoFile = open("/home/maagj01/scratch/transposon/Analysis/Nick_Mamrak/RnaAsDna/march12_RnaAsDNA.tsv", 'r')
flowcellInfoFile = flowcellInfoFile.readlines()
flowcellInfoFile = [line.rstrip().split('\t') for line in flowcellInfoFile]

startRow=0
while(flowcellInfoFile[startRow][0] != "#Sample Name"):
    startRow+=1

flowcellFile = pd.DataFrame(flowcellInfoFile[startRow:], columns=flowcellInfoFile[startRow]).iloc[1:]

#TODO Iterate through to generate inputs.txt

for index, Sample in flowcellFile.iterrows():
    if args.verbose:
        print("\nParsing:\n" + str(Sample), file=sys.stderr)
    annotateSample(Sample, project)


if doDNaseCleanup:
    print()
    print("rqsub -b y -j y -N analyzeInserts -hold_jid `cat sgeid.analysis | perl -pe 's/\\n/,/g;'` \"/vol/mauranolab/mapped/src/analyzeInserts.R\"")
    print("bqsub -j y -N mapped_readcounts -hold_jid `cat sgeid.analysis | perl -pe 's/\\n/,/g;'` /vol/mauranolab/mapped/src/mapped_readcounts.sh")
