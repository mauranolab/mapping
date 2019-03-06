#!/bin/env python3.5
import sys
import re
import argparse
import pandas as pd


###Argument parsing
version="1.3"

parser = argparse.ArgumentParser(prog = "createFlowcellSubmit", description = "Outputs submit commands for all samples from the info.txt located in the flowcell data folder", allow_abbrev=False)

parser.add_argument('--flowcellIDs', action='store', type = str, help='Flowcell IDs (will look for /vol/mauranolab/flowcells/FCxxx/info.txt, multiple FCs separated by comma)', required=True)
parser.add_argument('--sampletypes', action='store', type = str, help='Only process samples of this type (multiple projects separated by comma)', required=False)
parser.add_argument('--projects', action='store', type = str, help='Only process samples in these projects (multiple projects separated by comma)', required=False)
parser.add_argument('--samples', action='store', type = str, help='Only process samples matching this BS number (multiple samples separated by comma)', required=False)
parser.add_argument("--aggregate", action='store_true', default=False, help = "Aggregate samples with the same BS number")
parser.add_argument("--aggregate-sublibraries", action='store_true', default=False, help = "Aggregate samples with the same BS number and sublibrary, and remark duplicates")
parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s ' + version)

try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument)
    sys.exit(2)

if args.verbose:
    print("\nargs:" + str(args), file=sys.stderr)

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


#File locations based on lab/sampleType
def getBasedir(lab, sampleType):
    if args.aggregate or args.aggregate_sublibraries:
        if sampleType in ['DNA', 'DNA Capture', 'DNase-seq', 'Nano-DNase', 'ChIP-seq']:
            if lab == "Maurano":
                basedir = "/vol/mauranolab/mapped/"
            elif lab == "CEGS":
                basedir = "/vol/cegs/mapped/"
            else:
                raise Exception("Don't know how to aggregate data from " + lab)
        elif sampleType in ['Transposon DNA', 'Transposon RNA', 'Transposon iPCR', 'Transposon iPCR Capture']:
           basedir = "/vol/mauranolab/transposon/"
        else:
            raise Exception("Don't know how to aggregate data for " + sampleType) 
    else:
        basedir="/vol/mauranolab/flowcells/fastq/"
    return basedir


###Handlers for individual data types
#Each function takes dict representing a single row from the FC file iterator and returns the processing command line
#Transposon pipeline
def aggregateTransposonSamples(lines):
    #Exact parameters don't matter
    basedir = getBasedir(None, "Transposon DNA")
    
    return '/vol/mauranolab/transposon/src/submitMerge.sh ' + lines.iloc[0]['Original Sample #'] + "-" + lines.iloc[0]['#Sample Name'] + " " + ' '.join([ basedir + fc + "/" + lines.iloc[0]['Original Sample #'] + "-" + lines.iloc[0]['#Sample Name'] + '/' for fc in lines['FC'].tolist() ])


def transposonSamples(line):
    fullSampleName = line["Sample #"] + "-" + line["#Sample Name"]
    sampleType = line["Sample Type"]
    
    #Might like to put submit logfile in sample directory but would need to mkdir first, or change qsub to do mkdir on -o: ' -o ' + fullSampleName + 
    qsub = "qsub -S /bin/bash -j y -N submit." + fullSampleName + ' -b y "/vol/mauranolab/transposon/src/'
    
    fileLocation = getBasedir(line['Lab'], sampleType) + line["FC"] + "/Project_" + line["Lab"] + "/Sample_" + line["Sample #"] + '/"'
    
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
    submitCommand = ' '.join([qsub + program, fullSampleName, BCreadSeq, line["R1 Trim (P5)"], line["R2 Trim (P7)"], bclen, bcread, plasmidSeq, extractBCargs, fileLocation])
    
    global doTransposonCleanup
    doTransposonCleanup = True
    
    return(submitCommand)


#Chromosome configuration capture
#TODO fix organism
def chromConfCapture(line):
    organism = "hg38" 
    fullSampleName = line["Sample #"] + "-" + line["#Sample Name"]
    
    qsub = "qsub -S /bin/bash -j y -N submit." + fullSampleName + \
    ' -b y "/home/maagj01/scratch/transposon/src/submitHiC.sh ' +fullSampleName 
    
    sampleUnique = "/home/maagj01/scratch/transposon/captureC/config/config_Dpn.txt  /home/maagj01/scratch/transposon/captureC/genomeFrag/hg38_dpnii.bed K562_Dpn " + organism
    
    fileLocation = getBasedir(line['Lab'], line["Sample Type"]) + flowcellID + "/Project_" + line["Lab"] + "/Sample_" + line["Sample #"] + '/"'
    
    submitCommand = qsub + sampleUnique + " " + fileLocation
    return(submitCommand)


def bwaPipeline(sampleName, sampleID, sampleType, mappedgenome, processingCommand, analysisCommand):
    submitCommand = "/vol/mauranolab/mapped/src/submit.sh " + mappedgenome + " " + processingCommand + "," + analysisCommand + " " + sampleName + " " + sampleID
    
    global doDNaseCleanup
    doDNaseCleanup = True
    
    global DNaseFragmentLengthPlot
    fragmentLengths = { 'DNase-seq': 300, 'ChIP-seq': 500, 'DNA': 750 , 'DNA Capture': 750 }
    if fragmentLengths[sampleType] > DNaseFragmentLengthPlot:
        DNaseFragmentLengthPlot = fragmentLengths[sampleType]
    
    return(submitCommand)


#DNase/ ChIP-seq
def DNase(line):
    speciesToGenomeReference = {
        'Human': 'hg38_noalt',
        'Mouse': 'mm10'
    }
    mappedgenome = speciesToGenomeReference[line["Species"]]
    
    if args.aggregate:
        processingCommand="aggregate"
    elif args.aggregate_sublibraries:
        processingCommand="aggregateRemarkDups"
    else:
        processingCommand="mapBwaAln"
    
    sampleType = line["Sample Type"]
    
    return(bwaPipeline(line["#Sample Name"], line["Sample #"], sampleType, mappedgenome, processingCommand, "chipseq" if sampleType == "ChIP-seq" else "dnase"))


def DNA(line):
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
    mappedgenome = speciesToGenomeReference[line["Species"]]
    
    if args.aggregate:
        processingCommand="aggregate"
    elif args.aggregate_sublibraries:
        processingCommand="aggregateRemarkDups"
    else:
        processingCommand="mapBwaMem"
    
    sampleType = line["Sample Type"]
    
    if sampleType=="DNA Capture":
        global doDNACaptureCleanup
        doDNACaptureCleanup = True
    
    return(bwaPipeline(line["#Sample Name"], line["Sample #"], sampleType, mappedgenome, processingCommand, "callsnps"))



####
#Main body - read in flowcell and create a submit script for all samples
####

###Parse flowcell info
flowcellFile = pd.DataFrame()
for flowcellID in flowcellIDs:
    flowcellInfoFileName = "/vol/mauranolab/flowcells/data/" + flowcellID + "/info.txt"
    
    flowcellInfoFile = open(flowcellInfoFileName, 'r')
    #flowcellInfoFile = open("/vol/mauranolab/flowcells/data/FCHM2LKBGX9/info.txt", 'r')
    flowcellInfoFile = flowcellInfoFile.readlines()
    flowcellInfoFile = [line.rstrip().split('\t') for line in flowcellInfoFile]
    
    startRow=0
    while(flowcellInfoFile[startRow][0] != "#Sample Name"):
        startRow+=1
    
    if args.verbose:
        print("Reading", flowcellInfoFileName, "starting at line", startRow, file=sys.stderr)
    
    curFlowcellFile = pd.DataFrame(flowcellInfoFile[startRow:], columns = flowcellInfoFile[startRow]).iloc[1:]
    curFlowcellFile['FC'] = flowcellID
    flowcellFile = pd.concat([flowcellFile, curFlowcellFile])

#Subset to samples we're interested in
flowcellFile = flowcellFile[flowcellFile['Sample Type' ] != "Pool"]
if sampletypes is not None:
    flowcellFile = flowcellFile[flowcellFile['Sample Type'].isin(sampletypes)]
if projects is not None:
    flowcellFile = flowcellFile[flowcellFile['Lab'].isin(projects)]
if samples is not None:
    flowcellFile = flowcellFile[flowcellFile['Sample #'].isin(samples)]

#Adjust sample IDs and drop duplicate rows to handle aggregations
if args.aggregate or args.aggregate_sublibraries:
    flowcellFile['Original Sample #'] = flowcellFile['Sample #']
    if args.aggregate:
        flowcellFile['Sample #'] = flowcellFile['Sample #'].apply(lambda bs: re.sub("(BS[0-9]{5})[A-Z]", "\\1", bs))
    #BWA pipeline just needs last entry per sample
    if len(set(flowcellFile['Sample Type']).intersection(set(['DNA', 'DNA Capture', 'DNase-seq', 'Nano-DNase', 'ChIP-seq']))) > 0:
        flowcellFile = flowcellFile[flowcellFile.duplicated('Sample #', keep='last')]


###Pre-processing
print("#", ' '.join(sys.argv), sep="")

#Initialize inputs.txt
if len(set(flowcellFile['Sample Type']).intersection(set(['DNA', 'DNA Capture', 'DNase-seq', 'Nano-DNase', 'ChIP-seq']))) > 0:
    if flowcellIDs is not None and projects is not None and len(projects) is 1:
        #Just handle the single-project case
        #getBasedir returns the same for all these types, so just hardcode DNA
        if args.aggregate or args.aggregate_sublibraries:
            print('find ' + ' '.join([(getBasedir(projects[0], 'DNA')  + '%/Project_' + str(projects[0]) + '/*').replace("%", s) for s in flowcellIDs]) + ' -maxdepth 1 -name "*.bam" | sort > inputs.txt')
        else:
            print('find ' + ' '.join([(getBasedir(projects[0], 'DNA')  + '%/Project_' + str(projects[0]) + '/*').replace("%", s) for s in flowcellIDs]) + ' -maxdepth 1 -name "*.fastq.gz" | sort > inputs.txt')
        print()


#Some aggregations must be processed in groups of samples
if args.aggregate or args.aggregate_sublibraries:
    if len(set(flowcellFile['Sample Type']).intersection(set(['Transposon DNA', 'Transposon RNA', 'Transposon iPCR', 'Transposon iPCR Capture']))) > 0:
        for curSample in set(flowcellFile[flowcellFile['Sample Type'].isin(['Transposon DNA', 'Transposon RNA', 'Transposon iPCR', 'Transposon iPCR Capture'])]['Sample #']):
            print(aggregateTransposonSamples(flowcellFile[flowcellFile['Sample #']==curSample]))


###Dispatch appropriate function handler per sample line
doDNaseCleanup = False
doDNACaptureCleanup = False
DNaseFragmentLengthPlot = 0
doTransposonCleanup = False

for index, line in flowcellFile.iterrows():
    if args.verbose:
        print("\nParsing:\n" + str(line), file=sys.stderr)
    sampleName = line["#Sample Name"]
    sampleID = line["Sample #"]
    sampleType = line["Sample Type"]
    
    try:
        if re.compile(r'^BS[0-9]{5}[A-Z]?$').search(sampleID) is None:
            raise Exception("Can't parse " + sampleID + " as BS number")
        
        if sampleType == "Transposon DNA" or sampleType == "Transposon RNA" or sampleType == "Transposon iPCR" or sampleType == "Transposon iPCR Capture":
            if  args.aggregate or args.aggregate_sublibraries:
                continue
            qsubLine = transposonSamples(line)
        elif sampleType == "Hi-C" or sampleType == "3C-seq" or sampleType == "Capture-C":
            qsubLine =  "#" + chromConfCapture(line)
        elif sampleType == "Nano-DNase" or sampleType == "ChIP-seq" or sampleType == "DNase-seq":
            qsubLine = DNase(line)
        elif sampleType == "DNA" or sampleType == "DNA Capture":
            qsubLine = DNA(line)
        else:
            raise Exception("Don't know how to process " + sampleType)
        
        #We have successfully generated the command line
        print(qsubLine)
    except Exception as e:
        print("WARNING for sample ", sampleName, ":", e, file=sys.stderr)
        print("#Couldn't process " + sampleName + "-" + sampleID)


###Cleanup calls
if doDNaseCleanup:
    print()
    #Just use max fragment length for now
    #TODO run separate plots for each sample type?
    print("qsub -b y -j y -N analyzeInserts -hold_jid `cat sgeid.analysis | perl -pe 's/\\n/,/g;'` \"/vol/mauranolab/mapped/src/analyzeInserts.R " + str(DNaseFragmentLengthPlot) + "\"")
    print("qsub -S /bin/bash -j y -N mapped_readcounts -hold_jid `cat sgeid.analysis | perl -pe 's/\\n/,/g;'` /vol/mauranolab/mapped/src/mapped_readcounts.sh")
    print("rm -f sgeid.analysis")


#TODO placeholder for now
if doDNACaptureCleanup:
    print()
    print("#qsub -b y -S /bin/bash -j y -N analyzeCapture -hold_jid `cat sgeid.analysis | perl -pe 's/\\n/,/g;'` \"/vol/mauranolab/mapped/src/analyzeCapture.sh mappedgenome bait dirs\"")


#BUGBUG doesn't work since we don't have the pid of the final job at submission time
#if doTransposonCleanup:
#    print()
#    print("qsub -S /bin/bash -b y -j y -N transposon_info -hold_jid `cat sgeid.merge | perl -pe 's/\\n/,/g;'` \"/vol/mauranolab/transposon/src/Flowcell_Info.sh /home/mauram01/public_html/flowcellInfo/\"")
#    print("rm -f sgeid.merge")

