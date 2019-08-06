#!/bin/env python3.5
import sys
import re
import argparse
import pandas as pd
import pygsheets
import os
import glob

#Get LIMS info
sys.path.append("/vol/mauranolab/mapped/src/flowcells")
from lims import getLIMSsheet, getValueFromLIMS


###Argument parsing
version="1.5"

parser = argparse.ArgumentParser(prog = "createFlowcellSubmit", description = "Outputs submit commands for all samples from the info.txt located in the flowcell data folder", allow_abbrev=False)

parser.add_argument('--flowcellIDs', action='store', type = str, help='Flowcell IDs (will look for /vol/mauranolab/flowcells/FCxxx/info.txt, multiple FCs separated by comma)', required=True)
parser.add_argument('--samples', action='store', type = str, help='Only process samples matching this BS number (multiple samples separated by comma, done as string matching)', required=False)
parser.add_argument('--projects', action='store', type = str, help='Only process samples in these projects (multiple projects separated by comma)', required=False)
parser.add_argument('--people', action='store', type = str, help='Only process samples made by these people (multiple names separated by comma)', required=False)
parser.add_argument('--sampletypes', action='store', type = str, help='Only process samples of this type (multiple projects separated by comma)', required=False)
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

if args.people is None:
    people = None
else:
    people = args.people.split(",")

if args.sampletypes is None:
    sampletypes = None
else:
    sampletypes = args.sampletypes.split(",")

if args.samples is None:
    samples = None
else:
    samples = args.samples.split(",")


###Utility functions
def quoteStringsWithSpaces(str):
    if " " in str:
        str = "\"" + str + "\""
    return str


#Get base location of fastq or mapped files based on lab/sampleType
def getBasedir(lab, sampleType, fc):
    if args.aggregate or args.aggregate_sublibraries:
        if sampleType in bwaPipelineAnalysisCommandMap.keys():
            if lab == "Maurano":
                basedir = "/vol/mauranolab/mapped/"
            elif lab == "CEGS":
                basedir = "/vol/cegs/mapped/"
            else:
                raise Exception("Don't know how to aggregate data from " + lab)
            basedir += fc + "/" + bwaPipelineAnalysisCommandMap[sampleType] + "/"
        elif sampleType in ['Transposon DNA', 'Transposon RNA', 'Transposon 10xRNA', 'Transposon iPCR', 'Transposon iPCR Capture']:
           basedir = "/vol/mauranolab/transposon/" + fc
        else:
            raise Exception("Don't know how to aggregate data for " + sampleType) 
    else:
        basedir="/vol/mauranolab/flowcells/fastq/" + fc + "/Project_" + lab + "/"
    return basedir


###Handlers for individual data types
#Each function takes dict representing a single row from the FC file iterator and returns the processing command line
#Transposon pipeline
def aggregateTransposonSamples(lines):
    return '/vol/mauranolab/mapped/src/transposon/submitMerge.sh ' + lines.iloc[0]['Sample #'] + "-" + lines.iloc[0]['#Sample Name'] + " " + ' '.join([ getBasedir(None, line['Sample Type'], line['FC']) + "/" + line['Original Sample #'] + "-" + line['#Sample Name'] + '/' for index, line in lines.iterrows() ])


def transposonSamples(line):
    fullSampleName = line["Sample #"] + "-" + line["#Sample Name"]
    sampleType = line["Sample Type"]
    sampleTypeShort= sampleType.split(' ')[1]
    
    #Might like to put submit logfile in sample directory but would need to mkdir first, or change qsub to do mkdir on -o: ' -o ' + fullSampleName + 
    qsub = "qsub -S /bin/bash -j y -b y -N submit." + fullSampleName + ' "/vol/mauranolab/mapped/src/transposon/submit.sh '
    
    fileLocation = getBasedir(line['Lab'], sampleType, line['FC']) + "Sample_" + line["Sample #"] + '/"'
    
    #default parameters that can be overriden below
    BCreadSeq = "CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGGGATCTTTGTCCAAACTCATCGAGCTCGGGA"
    bclen = "16"
    
    #Unique line to each library 
    if sampleType in ['Transposon DNA', 'Transposon RNA', 'Transposon 10xRNA']:
        plasmidSeq="AGCTGCACAGCAACACCGAGCTGGGCATCGTGGAGTACCAGCACGCCTTCAAGACCCCCATCGCCTTCGCCAGATC"
        extractBCargs="--no-BCrevcomp"
        if sampleType == "Transposon DNA":
            bcread = "R2"
        elif sampleType == "Transposon RNA":
            bcread = "R1"
        elif sampleType == "Transposon 10xRNA":
            BCreadSeq = "TTGGACAAAGATCCCACAACTAGCBBBBBBBBBBBBBBBBCCTAATTGCCTAGAAAGGAGCAGACGATATGGCGTCGCTCC"
            bcread = "R2"
            extractBCargs="--BCrevcomp"
            plasmidSeq="None" #Nothing to align to in the cell BC / UMI (just the pT tail after).
        else:
            raise Exception("IMPOSSIBLE")
    elif sampleType == "Transposon iPCR" or sampleType == "Transposon iPCR Capture":
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
    submitCommand = ' '.join([qsub + fullSampleName, sampleTypeShort, line["R1 Trim (P5)"], line["R2 Trim (P7)"], bcread, bclen, BCreadSeq, plasmidSeq, extractBCargs, fileLocation])
    
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
    
    fileLocation = getBasedir(line['Lab'], line["Sample Type"], line['FC']) + "Sample_" + line["Sample #"] + '/"'
    
    submitCommand = qsub + sampleUnique + " " + fileLocation
    return(submitCommand)


#BWA pipeline for DNA / DNase / ChIP-seq
bwaPipelineAnalysisCommandMap = { "DNase-seq": "dnase", "Nano-DNase": "dnase", "ChIP-seq": "chipseq", "DNA": "dna", "DNA Capture": "capture" }
bwaPipelineFragmentLengthsMap = { 'DNase-seq': 300, 'Nano-DNase': 300, 'ChIP-seq': 500, 'DNA': 750 , 'DNA Capture': 750 }


def getBwaPipelineOutdir(sampleType):
    if args.aggregate or args.aggregate_sublibraries:
        bwaPipelineOutdir = ""
    else:
        #Segregate FC directories by assay type
        bwaPipelineOutdir = bwaPipelineAnalysisCommandMap[sampleType] + "/"
    return(bwaPipelineOutdir)


def addCEGSgenomes(line):
    if line['Lab'] != "CEGS":
        return []
    else:
        sampleName = line["#Sample Name"]
        
#        Disable finding reference -- only support getting it through google sheets
#        if re.search(r'_(Amplicon|BAC|Yeast)$', sampleName) is not None:
#            return ",cegsvectors_" + "_".join(sampleName.split("_")[0:3])
        
        additionalGenomesToMap = []
        geneticModification = getValueFromLIMS(lims, line['Original Sample #'], 'Genetic modification')
        for curGeneticModification in geneticModification.split(","):
            for curCegsGenome in cegsGenomes:
                #Genome must match entire contents of a field (e.g., HPRT1 doesn't match dHPRT1), excepting the square brackets denoting integration site
                if curCegsGenome == re.sub(r'\[.+\]$', '', curGeneticModification):
                    #Add the genome including the cegsvectors_ prefix
                    additionalGenomesToMap.append("cegsvectors_" + curCegsGenome)
                    #Not possible to match any further genomes since wildcards are not allowed
                    break
        return additionalGenomesToMap

def bwaPipeline(line):
    sampleType = line["Sample Type"]
    
    if sampleType in ["DNA", "DNA Capture"]:
        speciesToGenomeReference = {
            'Human': 'hg38_full',
            'Mouse': 'mm10',
            'Rat': 'rn6',
            'Human+yeast': 'hg38_sacCer3',
            'Mouse+yeast': 'mm10_sacCer3',
            'Rat+yeast': 'rn6_sacCer3',
        }
    elif sampleType in ["Nano-DNase", "ChIP-seq", "DNase-seq"]:
        speciesToGenomeReference = {
            'Human': 'hg38_noalt',
            'Mouse': 'mm10',
            'Rat': 'rn6'
        }
    else:
        raise Exception("Can't parse " + sampleType)
    mappedgenomes = [ speciesToGenomeReference[curSpecies ] for curSpecies in getValueFromLIMS(lims, line['Original Sample #'], 'Species').split(",") ]
    mappedgenomes += addCEGSgenomes(line)
    
    if args.aggregate:
        processingCommand="aggregate"
    elif args.aggregate_sublibraries:
        processingCommand="aggregateRemarkDups"
    else:
        if sampleType in ["DNA", "DNA Capture"]:
            processingCommand="mapBwaMem"
        elif sampleType in ["Nano-DNase", "ChIP-seq", "DNase-seq"]:
            processingCommand="mapBwaAln"
        else:
            raise Exception("Can't parse " + sampleType)
    
    if line["Sample Type"]=="DNA Capture":
        global doDNACaptureCleanup
        doDNACaptureCleanup = True
    
    sampleAnnotation = [ ]
    sex = getValueFromLIMS(lims, line['Original Sample #'], 'Sex')
    if sex is not "":
        sampleAnnotation.append("Sex=" + sex)
    bait_set = getValueFromLIMS(lims, line['Original Sample #'], 'Bait set')
    if bait_set is not "":
        sampleAnnotation.append("Bait_set=" + bait_set)
    
    submitCommand = "/vol/mauranolab/mapped/src/dnase/submit.sh " + ",".join(sorted(mappedgenomes)) + " " + processingCommand + "," + bwaPipelineAnalysisCommandMap[sampleType] + " " + getBwaPipelineOutdir(sampleType) + line["#Sample Name"] + " " + line["Sample #"] + " " + ','.join(sampleAnnotation)
    
    global doBwaCleanup
    doBwaCleanup = True
    
    return(submitCommand)


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


###Pre-processing
print("#", ' '.join([ quoteStringsWithSpaces(arg) for arg in sys.argv ]), sep="")
print()

#Subset to samples we're interested in
flowcellFile = flowcellFile[flowcellFile['Sample Type' ] != "Pool"]
if sampletypes is not None:
    flowcellFile = flowcellFile[flowcellFile['Sample Type'].isin(sampletypes)]
if projects is not None:
    flowcellFile = flowcellFile[flowcellFile['Lab'].isin(projects)]
if people is not None:
    flowcellFile = flowcellFile[flowcellFile['Made By'].isin(people)]
if samples is not None:
    #Do partial string matching rather than exact equality to allow flexible subsetting of certain sublibraries or of all sublibraries for a given BS number
    flowcellFile = flowcellFile[flowcellFile['Sample #'].apply(lambda bs: any([s in bs for s in samples]))]


#Initialize inputs.txt - must be done before we drop duplicate sample rows
if len(set(flowcellFile['Sample Type']).intersection(set(bwaPipelineFragmentLengthsMap.keys()))) > 0:
    dirs = set([])
    for index, line in flowcellFile[flowcellFile['Sample Type'].isin(bwaPipelineFragmentLengthsMap.keys())].iterrows():
        dirs.add(getBasedir(line['Lab'], line['Sample Type'], line['FC']) + '*')
    findCmd =  'find ' + ' '.join(dirs) + ' -maxdepth 1 -regextype posix-awk -regex "^.+(' + '|'.join(sorted(flowcellFile['Sample #'].unique())) + ').+$"'
    if args.aggregate or args.aggregate_sublibraries:
        findCmd += ' -name "*.bam"'
    else:
        findCmd += ' -name "*.fastq.gz"'
    print(findCmd + ' | sort > inputs.txt')
    print()


#Adjust sample IDs and drop duplicate rows to handle aggregations
flowcellFile['Original Sample #'] = flowcellFile['Sample #']
if args.aggregate or args.aggregate_sublibraries:
    #Sort aggregations by BS number since sample sheet order doesn't matter
    flowcellFile = flowcellFile.sort_values(by='Sample #')
    if args.aggregate:
        flowcellFile['Sample #'] = flowcellFile['Sample #'].apply(lambda bs: re.sub("(BS[0-9]{5})[A-Z]", "\\1", bs))
    #BWA pipeline just needs last entry per sample
    if len(set(flowcellFile['Sample Type']).intersection(set(bwaPipelineFragmentLengthsMap.keys()))) > 0:
        flowcellFile = flowcellFile[~flowcellFile.duplicated('Sample #', keep='last')]


#Some aggregations must be processed in groups of samples
if args.aggregate or args.aggregate_sublibraries:
    if len(set(flowcellFile['Sample Type']).intersection(set(['Transposon DNA', 'Transposon RNA', 'Transposon 10xRNA', 'Transposon iPCR', 'Transposon iPCR Capture']))) > 0:
        for curSample in set(flowcellFile[flowcellFile['Sample Type'].isin(['Transposon DNA', 'Transposon RNA', 'Transposon 10xRNA', 'Transposon iPCR', 'Transposon iPCR Capture'])]['Sample #']):
            print(aggregateTransposonSamples(flowcellFile[flowcellFile['Sample #']==curSample]))


###Dispatch appropriate function handler per sample line
limsWks, lims, limsMask = getLIMSsheet("LIMS")

#Will map to these custom genomes when specified, stored as they appear in LIMS (without cegsvectors_ prefix)
cegsGenomes = [ re.sub(r'^cegsvectors_', '', os.path.basename(x)) for x in glob.glob("/vol/cegs/sequences/cegsvectors_*") ]

doBwaCleanup = False
doDNACaptureCleanup = False
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
        
        if sampleType in ['Transposon DNA', 'Transposon RNA', 'Transposon 10xRNA', 'Transposon iPCR', 'Transposon iPCR Capture']:
            if args.aggregate or args.aggregate_sublibraries:
                continue
            qsubLine = transposonSamples(line)
        elif sampleType in ["Hi-C", "3C-seq", "Capture-C"]:
            qsubLine =  "#" + chromConfCapture(line)
        elif sampleType in bwaPipelineAnalysisCommandMap:
            qsubLine = bwaPipeline(line)
        else:
            raise Exception("Don't know how to process " + sampleType)
        
        #We have successfully generated the command line
        print(qsubLine)
    except Exception as e:
        print("WARNING for sample ", sampleName, ":", e, sep="", file=sys.stderr)
        print("#Couldn't process " + sampleName + "-" + sampleID)


###Cleanup calls
if doBwaCleanup:
    for sampleType in flowcellFile['Sample Type'][flowcellFile['Sample Type'].isin(bwaPipelineAnalysisCommandMap.keys())].unique():
        basedir = getBwaPipelineOutdir(sampleType)
        if args.aggregate or args.aggregate_sublibraries:
            sgeoutput = ""
        else:
            sgeoutput = " -o " + basedir
        
        print()
        print("#" + sampleType)
        print("qsub -b y -j y -N analyzeInserts" + sgeoutput + " -hold_jid `cat " + basedir + "sgeid.analysis | perl -pe 's/\\n/,/g;'` \"/vol/mauranolab/mapped/src/dnase/analyzeInserts.R " + str(bwaPipelineFragmentLengthsMap[sampleType]) + " " + basedir + "\"")
        print("qsub -S /bin/bash -j y -N mapped_readcounts" + sgeoutput + " -hold_jid `cat " + basedir + "sgeid.analysis | perl -pe 's/\\n/,/g;'` \"/vol/mauranolab/mapped/src/dnase/mapped_readcounts.sh " + basedir + "\"")
    #Leave this around so the hold_jid args don't fail
    #    print("rm -f sgeid.analysis")


#TODO placeholder for now
if doDNACaptureCleanup:
    print()
    print("#qsub -b y -S /bin/bash -j y -N analyzeCapture -o " + bwaPipelineAnalysisCommandMap[sampleType] + " -hold_jid `cat " + getBwaPipelineOutdir('DNA Capture') + "sgeid.analysis | perl -pe 's/\\n/,/g;'` \"/vol/mauranolab/mapped/src/dnase/analyzeCapture.sh mappedgenome bait dirs\"")


#BUGBUG doesn't work since we don't have the pid of the final job at submission time
#if doTransposonCleanup:
#    print()
#    print("qsub -S /bin/bash -b y -j y -N transposon_info -hold_jid `cat sgeid.merge | perl -pe 's/\\n/,/g;'` \"/vol/mauranolab/mapped/src/transposon/Flowcell_Info.sh /home/mauram01/public_html/flowcellInfo/\"")
#    print("rm -f sgeid.merge")

