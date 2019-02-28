#!/bin/env python3.5
import sys
import re
import argparse
import pandas as pd


basedir="/vol/mauranolab/flowcells"


###Argument parsing
version="1.2"

parser = argparse.ArgumentParser(prog = "createFlowcellSubmit", description = "Outputs submit commands for all samples from the info.txt located in the flowcell data folder", allow_abbrev=False)

parser.add_argument('--flowcellIDs', action='store', type = str, help='Flowcell IDs (will look for /vol/mauranolab/flowcells/FCxxx/info.txt, multiple FCs separated by comma)', required=True)
parser.add_argument('--sampletypes', action='store', type = str, help='Only process samples of this type (multiple projects separated by comma)', required=False)
parser.add_argument('--projects', action='store', type = str, help='Only process samples in these projects (multiple projects separated by comma)', required=False)
parser.add_argument('--samples', action='store', type = str, help='Only process samples matching this BS number (multiple samples separated by comma)', required=False)

#TODO
#parser.add_argument("--aggregate", action='store_true', default=False, help = "Aggregate samples")
#parser.add_argument("--aggregate-sublibrary", action='store_true', default=False, help = "Aggregate samples")

parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s ' + version)

try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument)
    sys.exit(2)


if args.flowcellIDs is None:
    flowcellIDs = None
else:
    flowcellIDs = args.flowcellIDs.split(",")

if args.projects is None:
    projects = None
else:
    projects = args.projects.split(",")

if args.sampletypes is None:
    sampletypes = None
else:
    sampletypes = args.sampletypes.split(",")

if args.samples is None:
    samples = None
else:
    samples = args.samples.split(",")


###Handlers for individual data types
#Transposon pipeline
def transposonSamples(sampleName, sampleID, lab, sampleType, species, r1Trim, r2Trim):
    fullSampleName = sampleID + "-" + sampleName
    fileLocation = basedir + "/fastq/" + flowcellID + "/Project_" + lab + "/Sample_" + sampleID + '/"'
    
    #Would like to put submit logfile in sample directory bu would need to mkdir first, or change qsub to do mkdir on -o: ' -o ' + fullSampleName + 
    qsub = "qsub -S /bin/bash -j y -N submit." + fullSampleName + ' -b y "/vol/mauranolab/transposon/src/'
    
    #default parameters that can be overriden below
    BCreadSeq = "CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGGGATCTTTGTCCAAACTCATCGAGCTCGGGA"
    bclen = "16"
    
    #Unique line to each library 
    if sampleType == "Transposon DNA" or sampleType == "Transposon RNA":
        program= "submitBCcounts.sh"
        plasmidSeq="AGCTGCACAGCAACACCGAGCTGGGCATCGTGGAGTACCAGCACGCCTTCAAGACCCCCATCGCCTTCGCCAGATC"
        extractBCargs="--no-BCrevcomp"
        if sampleType == "Transposon DNA":
            bcread = "R2"
        elif sampleType == "Transposon RNA":
            bcread = "R1"
        else:
            raise Exception("IMPOSSIBLE")
    elif sampleType == "Transposon iPCR" or sampleType == "Transposon iPCR Capture":
        program = "submitIntegrations.sh"
        if sampleType == "Transposon iPCR":
            bcread = "R1"
            plasmidSeq="CTTCCGACTTCAACTGTA"
            extractBCargs="--no-BCrevcomp"
        elif sampleType == "Transposon iPCR Capture":
            BCreadSeq="AAAGATCCCACAACTAGCBBBBBBBBBBBBBBBBCCTAATTGCCTAGAAAGGAGCAGACG"
            bcread = "R2"
            plasmidSeq="None"
            extractBCargs="--BCrevcomp"
        else:
            raise Exception("IMPOSSIBLE")
    else:
        raise Exception("Can't handle sampleType" + sampleType)
    
    global doTransposonCleanup
    doTransposonCleanup = True
    
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
    submitCommand = "/vol/mauranolab/mapped/src/submit.sh " + reference + " mapBwaAln," + ("chipseq" if sampleType == "ChIP-seq" else "dnase") + " " + sampleName + " " + sampleID
    
    global doDNaseCleanup
    doDNaseCleanup = True
    
    global DNaseFragmentLengthPlot
    if fragmentLengths[sampleType] > DNaseFragmentLengthPlot:
        DNaseFragmentLengthPlot = fragmentLengths[sampleType]
    
    return(submitCommand)


def DNA(sampleName, sampleID, lab, sampleType, species):
    speciesToGenomeReference = {
        'Human': 'hg38_full',
        'Mouse': 'mm10',
        'Rat': 'rn6',
        'Human+mouse': 'hg38_full,mm10',
        'Human+yeast': 'hg38_sacCer3',
        'Mouse+yeast': 'mm10_sacCer3',
        'Rat+yeast': 'rn6_sacCer3',
        'Mouse+rat': 'mm10,rn6'
    }
    reference = speciesToGenomeReference[species]
    fullSampleName = sampleID + "-" + sampleName
    submitCommand = "/vol/mauranolab/mapped/src/submit.sh " + reference + " mapBwaMem,callsnps " + sampleName + " " + sampleID
    
    #No specific cleanup yet
    global doDNaseCleanup
    doDNaseCleanup = True
    
    global DNaseFragmentLengthPlot
    if fragmentLengths[sampleType] > DNaseFragmentLengthPlot:
        DNaseFragmentLengthPlot = fragmentLengths[sampleType]
    
    if sampleType=="DNA Capture":
        global doDNACaptureCleanup
        doDNACaptureCleanup = True
    
    return(submitCommand)


###Dispatch appropriate function handler per sample line
#Takes dict representing a single row from the FC file iterator and returns the processing command line
def annotateSample(Sampleline, projects, sampletypes, samples):
    try:
        sampleName=Sampleline["#Sample Name"]
        sampleID=Sampleline["Sample #"]
        lab=Sampleline["Lab"]
        sampleType=Sampleline["Sample Type"]
        species=Sampleline["Species"]
        r1Trim=str(Sampleline["R1 Trim (P5)"])
        r2Trim=str(Sampleline["R2 Trim (P7)"])
        
        if sampletypes is not None and sampleType not in sampletypes:
            return
        
        if projects is not None and lab not in projects:
            return
        
        if re.compile(r'BS[0-9]{5}[A-Z]').search(sampleID) is None:
            raise Exception("Can't parse " + sampleID + "as BS number")
        
        if samples is not None and sampleID not in samples:
            return
        
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
doDNaseCleanup = False
doDNACaptureCleanup = False
fragmentLengths = { 'DNase-seq': 300, 'ChIP-seq': 500, 'DNA': 750 , 'DNA Capture': 750 }
DNaseFragmentLengthPlot = 0
doTransposonCleanup = False


print("#", ' '.join(sys.argv), sep="")
#TODO Iterate through to generate inputs.txt


for flowcellID in flowcellIDs:
    print("\n#Samples from ", flowcellID, sep="")
    flowcellInfoFileName = basedir + "/data/" + flowcellID + "/info.txt"
    flowcellInfoFile = open(flowcellInfoFileName, 'r')
    #flowcellInfoFile = open("/vol/mauranolab/flowcells/data/FCHM2LKBGX9/info.txt", 'r')
    flowcellInfoFile = flowcellInfoFile.readlines()
    flowcellInfoFile = [line.rstrip().split('\t') for line in flowcellInfoFile]
    
    startRow=0
    while(flowcellInfoFile[startRow][0] != "#Sample Name"):
        startRow+=1
    
    if args.verbose:
        print("Reading", flowcellInfoFileName, "starting at line", startRow, file=sys.stderr)
    
    flowcellFile = pd.DataFrame(flowcellInfoFile[startRow:], columns=flowcellInfoFile[startRow]).iloc[1:]
    
    #Just handle this simple case for now
    if flowcellIDs is not None and len(flowcellIDs) is 1 and projects is not None and len(projects) is 1:
        if 'DNA' in set(flowcellFile['Sample Type']) or 'DNA Capture' in set(flowcellFile['Sample Type']) or 'DNase-seq' in set(flowcellFile['Sample Type']) or 'Nano-DNase' in set(flowcellFile['Sample Type']) or 'ChIP-seq' in set(flowcellFile['Sample Type']):
            print('find /vol/mauranolab/flowcells/fastq/' + str(flowcellIDs[0]) + '/' + str(projects[0]) + '/* -maxdepth 1 -name "*.fastq.gz" | sort > inputs.txt')
            print()
    
    for index, line in flowcellFile.iterrows():
        if args.verbose:
            print("\nParsing:\n" + str(line), file=sys.stderr)
        annotateSample(line, projects, sampletypes, samples)


if doDNaseCleanup:
    print()
    #Just use max fragment length for now
    #TODO run separate plots for each sample type?
    print("qsub -b y -j y -N analyzeInserts -hold_jid `cat sgeid.analysis | perl -pe 's/\\n/,/g;'` \"/vol/mauranolab/mapped/src/analyzeInserts.R " + str(DNaseFragmentLengthPlot) + "\"")
    print("qsub -S /bin/bash -j y -N mapped_readcounts -hold_jid `cat sgeid.analysis | perl -pe 's/\\n/,/g;'` /vol/mauranolab/mapped/src/mapped_readcounts.sh")


#TODO placeholder for now
if doDNACaptureCleanup:
    print()
    print("#qsub -S /bin/bash -j y -N analyzeCapture -hold_jid `cat sgeid.analysis | perl -pe 's/\\n/,/g;'` /vol/mauranolab/mapped/src/analyzeCapture.sh mappedgenome bait dirs")


#BUGBUG doesn't work since we don't have the pid of the final job at submission time
#if doTransposonCleanup:
#    print()
#    print("qsub -S /bin/bash -b y -j y -N transposon_info -hold_jid `cat sgeid.merge | perl -pe 's/\\n/,/g;'` \"/vol/mauranolab/transposon/src/Flowcell_Info.sh /home/mauram01/public_html/flowcellInfo/\"")
#    print("rm -f sgeid.merge")

print("rm -f sgeid.analysis")
