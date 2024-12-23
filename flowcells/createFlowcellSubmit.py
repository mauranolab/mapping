#!/bin/env python
import sys
import re
import argparse
import pandas as pd
import pygsheets
import os
import glob

#Get LIMS info
sys.path.append("/gpfs/data/mauranolab/mapped/src/flowcells")
from lims import getLIMSsheet, getValueFromLIMS


###Argument parsing
version="1.6"

parser = argparse.ArgumentParser(prog = "createFlowcellSubmit", description = "Outputs submit commands for all samples from the info.txt located in the flowcell data folder", allow_abbrev=False)

parser.add_argument("--flowcellIDs", action="store", type = str, help="Flowcell IDs (will look for /gpfs/data/isg_sequencing/FCxxx/info.txt, multiple FCs separated by comma) [%(default)s]", required=True)
parser.add_argument("--samplenames", action="store", type = str, help="Only process samples with matching names (multiple samples separated by comma, done as string matching) [%(default)s]", required=False)
parser.add_argument("--samples", action="store", type = str, help="Only process samples matching this BS number (multiple samples separated by comma, done as string matching) [%(default)s]", required=False)
parser.add_argument("--projects", action="store", type = str, help="Only process samples in these projects (multiple projects separated by comma) [%(default)s]", required=False)
parser.add_argument("--people", action="store", type = str, help="Only process samples made by these people (multiple names separated by comma) [%(default)s]", required=False)
parser.add_argument("--sampletypes", action="store", type = str, help="Only process samples of this type (multiple projects separated by comma) [%(default)s]", required=False)
parser.add_argument("--aggregate", action="store_true", default=False, help = "Aggregate samples with the same BS number [%(default)s]")
parser.add_argument("--aggregate-sublibraries", action="store_true", default=False, help = "Aggregate samples with the same BS number and sublibrary, and remark duplicates [%(default)s]")
parser.add_argument("--verbose", action="store_true", default=False, help = "Verbose mode")
parser.add_argument("--version", action="version", version="%(prog)s " + version)

try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, "\n", exc.argument)
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

if args.samplenames is None:
    samplenames = None
else:
    samplenames = args.samplenames.split(",")


###Utility functions
def quoteStringsWithSpaces(str):
    if " " in str or "&" in str:
        str = "\"" + str + "\""
    return str


#Get base location of fastq or mapped files based on lab/sampleType
def getBasedir(lab, sampleType, fc):
    if args.aggregate or args.aggregate_sublibraries:
        if sampleType in bwaPipelineAnalysisCommandMap.keys():
            if lab == "Maurano":
                basedir = "/gpfs/data/mauranolab/mapped/"
            elif lab == "SARS":
                basedir = "/gpfs/data/mauranolab/sars/mapped/"
            elif lab == "CEGS":
                basedir = "/gpfs/data/cegs/mapped/"
            else:
                raise Exception("Don't know how to aggregate data from " + lab)
            basedir += fc + "/" + bwaPipelineAnalysisCommandMap[sampleType] + "/"
        elif sampleType in ["Transposon DNA", "Transposon RNA", "Transposon 10xRNA", "Transposon iPCR", "Transposon iPCR Capture"]:
           basedir = "/gpfs/data/mauranolab/transposon/" + fc
        else:
            raise Exception("Don't know how to aggregate data for " + sampleType) 
    else:
        basedir="/gpfs/data/isg_sequencing/fastq/" + fc + "/Project_" + lab + "/"
    return basedir


###FC loading
def getFlowcellInfoFromFile(flowcellID):
    flowcellInfoFileName = "/gpfs/data/isg_sequencing/data/" + flowcellID + "/info.txt"
    
    flowcellInfoFile = open(flowcellInfoFileName, "r")
    #flowcellInfoFile = open("/gpfs/data/isg_sequencing/data/FCHM2LKBGX9/info.txt", "r")
    flowcellInfoFileContent = flowcellInfoFile.readlines()
    
    flowcellInfoFile.close()
    
    flowcellInfoFileContent = [line.rstrip().split("\t") for line in flowcellInfoFileContent]
    
    startRow=0
    while(flowcellInfoFileContent[startRow][0] != "#Sample Name"):
        startRow+=1
    
    if args.verbose:
        print("Reading", flowcellInfoFileName, "starting at line", startRow, file=sys.stderr)
    
    curFlowcellFile = pd.DataFrame(flowcellInfoFileContent[startRow:], columns = flowcellInfoFileContent[startRow]).iloc[1:]
    curFlowcellFile["FC"] = flowcellID
    
    return curFlowcellFile


###Handlers for individual data types
#Each function takes dict representing a single row from the FC file iterator and returns the processing command line
#Transposon pipeline
def aggregateTransposonSamples(lines):
    return "/gpfs/data/mauranolab/mapped/src/transposon/submitMerge.sh " + lines.iloc[0]["Sample #"] + "-" + lines.iloc[0]["#Sample Name"] + " " + " ".join([ getBasedir(None, line["Sample Type"], line["FC"]) + "/" + line["Original Sample #"] + "-" + line["#Sample Name"] + "/" for index, line in lines.iterrows() ])


def transposonSamples(line):
    fullSampleName = line["Sample #"] + "-" + line["#Sample Name"]
    sampleType = line["Sample Type"]
    sampleTypeShort= sampleType.split(" ")[1]
    
    #Might like to put submit logfile in sample directory but would need to mkdir first, or change qsub to do mkdir on -o: ' -o ' + fullSampleName + 
    qsub = "qsub --time 4:00:00 --mem-per-cpu 4G -j y -S /bin/bash -N submit." + fullSampleName + ' "/gpfs/data/mauranolab/mapped/src/transposon/submit.sh '
    
    fileLocation = getBasedir(line["Lab"], sampleType, line["FC"]) + "Sample_" + line["Sample #"] + '/"'
    
    #default parameters that can be overriden below
    bclen = "16"
    
    if line["R1 Trim (P5)"]=="" or line["R2 Trim (P7)"]=="":
        raise Exception("Missing trim values")
    
    #Unique line to each library 
    if sampleType in ["Transposon DNA", "Transposon RNA"]:
        BCreadSeq = "CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGGGATCTTTGTCCAAACTCATCGAGCTCGGGA"
        plasmidSeq="AGCTGCACAGCAACACCGAGCTGGGCATCGTGGAGTACCAGCACGCCTTCAAGACCCCCATCGCCTTCGCCAGATC"
        extractBCargs="--no-BCrevcomp"
        if sampleType == "Transposon DNA":
            bcread = "R2"
        elif sampleType == "Transposon RNA":
            bcread = "R1"
        else:
            raise Exception("IMPOSSIBLE")
    elif sampleType == "Transposon iPCR":
        bcread = "R1"
        BCreadSeq = "CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBGCTAGTTGTGGGATC" #Shorter sequence ends at the DpnII site
        plasmidSeq="CTTCCGACTTCAACTGTA"
        extractBCargs="--no-BCrevcomp"
    elif sampleType == "Transposon iPCR Capture":
        bcread = "R2"
        BCreadSeq = "AAAGATCCCACAACTAGCBBBBBBBBBBBBBBBBCCTAATTGCCTAGAAAGGAGCAGACG"
        plasmidSeq = "None"
        extractBCargs="--BCrevcomp"
    elif sampleType == "Transposon 10xRNA":
        bcread = "R2"
        BCreadSeq = "TTGGACAAAGATCCCACAACTAGCBBBBBBBBBBBBBBBBCCTAATTGCCTAGAAAGGAGCAGACGATATGGCGTCGCTCC"
        plasmidSeq="None" #Nothing to align to in the cell BC / UMI (just the pT tail after).
        #NB hardcoded to 10x GEM lane 1, needs to be changed if from a different lane
        extractBCargs="\'--BCrevcomp  --cellBCsuffix 1\'"
    else:
        raise Exception("Can't handle sampleType" + sampleType)
    submitCommand = " ".join([qsub + fullSampleName, sampleTypeShort, line["R1 Trim (P5)"], line["R2 Trim (P7)"], bcread, bclen, BCreadSeq, plasmidSeq, extractBCargs, fileLocation])
    
    global doTransposonCleanup
    doTransposonCleanup = True
    
    return(submitCommand)


#Chromosome configuration capture
#TODO fix organism
def chromConfCapture(line):
    organism = "hg38" 
    fullSampleName = line["Sample #"] + "-" + line["#Sample Name"]
    
    qsub = "qsub --time 4:00:00 --mem-per-cpu 4G -j y -N submit." + fullSampleName + ' -b y -S /bin/bash "/home/maagj01/scratch/transposon/src/submitHiC.sh ' +fullSampleName
    
    sampleUnique = "/home/maagj01/scratch/transposon/captureC/config/config_Dpn.txt  /home/maagj01/scratch/transposon/captureC/genomeFrag/hg38_dpnii.bed K562_Dpn " + organism
    
    fileLocation = getBasedir(line["Lab"], line["Sample Type"], line["FC"]) + "Sample_" + line["Sample #"] + '/"'
    
    submitCommand = qsub + sampleUnique + " " + fileLocation
    return(submitCommand)


#BWA pipeline for DNA / DNase / ChIP-seq
bwaPipelineAnalysisCommandMap = { "DNase-seq": "dnase", "Nano-DNase": "dnase", "ATAC-seq": "atac", "ChIP-seq": "chipseq", "CUT&RUN": "chipseq", "DNA": "dna", "DNA Capture": "capture", "Amplicon": "amplicon" }
bwaPipelineFragmentLengthsMap = { "DNase-seq": 300, "Nano-DNase": 300, "ATAC-seq": 300, "ChIP-seq": 500, "CUT&RUN": 500, "DNA": 750 , "DNA Capture": 750, "Amplicon": 500 }


def getBwaPipelineOutdir(sampleType):
    if args.aggregate or args.aggregate_sublibraries:
        bwaPipelineOutdir = ""
    else:
        #Segregate FC directories by assay type
        bwaPipelineOutdir = bwaPipelineAnalysisCommandMap[sampleType] + "/"
    return(bwaPipelineOutdir)


def checkForCEGSreferences(requestedReferencesForMapping, warn=False):
    genomesToMap = []
    for curRequestedReference in requestedReferencesForMapping:
        #Note if geneticModification/customReference are empty strings, then the list will contain an empty string
        if curRequestedReference != "":
            #Genome must match entire contents of a field (e.g., HPRT1 doesn't match dHPRT1), excepting the square brackets denoting integration site
            if curRequestedReference in cegsGenomes:
                #Add the genome including the cegsvectors_ prefix
                genomesToMap.append("cegsvectors_" + curRequestedReference)
            else:
                if re.search(r'^(dHprt|dIgf2|dSox2|loxP|dPiga|dHoxa|dH2d|dH2k|dB2m|dHba|dPigaExon2|dLsp1|dSyt8|dHIDAD|dRet|dANK1|dCacna1c|dApobec3|dTAF1~SVA2|dIfn|dTAF1~SINE1|dTAF1~VNTR1|dTAF1~Alu1|dCACNA1C)', curRequestedReference) is None:
                    print("WARNING could not find reference for '" + curRequestedReference, "'", sep="", file=sys.stderr)
    return(genomesToMap)


def addCEGSgenomes(line):
    if line["Lab"] not in ["CEGS", "Chakravarti"]:
        return []
    else:
        sampleName = line["#Sample Name"]
        
#        Disable finding reference from name -- only support getting it through LIMS annotation
#        if re.search(r'_(Amplicon|BAC|Yeast)$', sampleName) is not None:
#            return ",cegsvectors_" + "_".join(sampleName.split("_")[0:3])
        
        geneticModifications = getValueFromLIMS(lims, line["Original Sample #"], "Genetic Modification").split(",")
        geneticModifications = [ re.sub(r"\[.+\]$", "", cur) for cur in geneticModifications ]
        #Take only first three (MenDel) fields in case there are comments/constructs (no longer supported -- must now match exactly)
        #geneticModifications = [ "_".join(cur.split("_")[0:3]) for cur in geneticModifications ]
        
        customReferences = getValueFromLIMS(lims, line["Original Sample #"], "Custom Reference").split(",")
        
        if args.verbose:
            print("\nWill look for genetic modifications [ ", ",".join(geneticModifications), " ] and custom references [ ",  ",".join(customReferences), " ].", sep="", file=sys.stderr)
        
        additionalGenomesToMap = checkForCEGSreferences(geneticModifications)
        additionalGenomesToMap = additionalGenomesToMap + checkForCEGSreferences(customReferences)
        
        return additionalGenomesToMap


def getBwaReference(species):
    if species == "":
        print("WARNING Species was not provided", sep="", file=sys.stderr)
        return None
    
    if sampleType in ["DNA", "DNA Capture", "Amplicon"]:
        speciesToGenomeReference = {
            "Human": "hg38_full",
            "Mouse": "mm10",
            "Castaneus": "mm10all_CASTEiJ_female",
            "Rat": "rn6",
            "Human+yeast": "hg38_sacCer3",
            "Mouse+yeast": "mm10_sacCer3",
            "Rat+yeast": "rn6_sacCer3",
            "Yeast": "sacCer3",
            "SARS-Cov2": "wuhCor1",
            "Human+SARS-Cov2": "hg38_full_wuhCor1",
            "Mole": "talOcc4",
            "Mole+yeast": "talOcc4_sacCer3",
            "Pufferfish+yeast": "sacCer3", #Do not have Fugu assembly at the moment
            "Pufferfish": None, #Do not have Fugu assembly at the moment
            "Plasmid": None,
            "Fly": None, #Do not have Drosophila assembly at the moment
            "T2T": "t2t",
        }
    elif sampleType in ["Nano-DNase", "ChIP-seq", "CUT&RUN", "DNase-seq", "ATAC-seq"]:
        speciesToGenomeReference = {
            "Human": "hg38_noalt",
            "Mouse": "mm10",
            "Castaneus": "mm10all_CASTEiJ_female",
            "Rat": "rn6",
            "Yeast": "sacCer3",
            "Plasmid": None,
            }
    else:
        raise Exception("Can't parse " + sampleType)
    
    if species not in speciesToGenomeReference:
        raise Exception("Can't find reference for species " + species, ".")
    
    return(speciesToGenomeReference[species])


def bwaPipeline(line):
    sampleType = line["Sample Type"]
    
    mappedgenomes = [ getBwaReference(curSpecies) for curSpecies in getValueFromLIMS(lims, line["Original Sample #"], "Species").split(",") ]
    if None in mappedgenomes:
        #None is used for a special case of species we silently do not map to, e.g. Plasmid is not handled right now except as a custom reference
        mappedgenomes.remove(None)
    mappedgenomes += addCEGSgenomes(line)
    
    if len(mappedgenomes) != len(set(mappedgenomes)):
        print("WARNING Duplicate genomes specified: " + str(mappedgenomes), sep="", file=sys.stderr)
    
    if len(mappedgenomes) == 0:
        raise Exception("No genomes specified for mapping")
    
    if args.aggregate:
        processingCommand="aggregate"
    elif args.aggregate_sublibraries:
        processingCommand="aggregateRemarkDups"
    else:
        if sampleType in ["DNA", "DNA Capture", "Amplicon"]:
            processingCommand="mapBwaMem"
        elif sampleType in ["Nano-DNase", "ChIP-seq", "CUT&RUN", "DNase-seq", "ATAC-seq"]:
            processingCommand="mapBwaAln"
        else:
            raise Exception("Can't parse " + sampleType)
    
    if line["Sample Type"] in ["DNA Capture"]:
        global doDNACaptureCleanup
        doDNACaptureCleanup = True
    
    if line["Lab"] in ["SARS"]:
        global doSARScleanup
        doSARScleanup = True
    
    sampleAnnotation = [ ]
    sex = getValueFromLIMS(lims, line["Original Sample #"], "Sex")
    if sex != "":
        sampleAnnotation.append("Sex=" + sex)
    bait_set = getValueFromLIMS(lims, line["Original Sample #"], "Bait set")
    if bait_set != "":
        sampleAnnotation.append("Bait_set=" + bait_set)
    geneticModification = getValueFromLIMS(lims, line["Original Sample #"], "Genetic Modification")
    if geneticModification != "":
        sampleAnnotation.append("Genetic_Modification=" + geneticModification)
    
    submitCommand = "/gpfs/data/mauranolab/mapped/src/dnase/submit.sh " + ",".join(sorted(set(mappedgenomes))) + " " + processingCommand + "," + bwaPipelineAnalysisCommandMap[sampleType] + " " + getBwaPipelineOutdir(sampleType) + line["#Sample Name"] + " " + line["Sample #"]
    if len(sampleAnnotation) > 0:
        submitCommand += " \"" + ";".join(sampleAnnotation) + "\""
    
    global doBwaCleanup
    doBwaCleanup = True
    
    return(submitCommand)


####
#Main body - read in flowcell and create a submit script for all samples
####

###Parse flowcell info
limsWks, lims, limsMask = getLIMSsheet("LIMS")

flowcellFile = pd.DataFrame()
for flowcellID in flowcellIDs:
    flowcellFile = pd.concat([flowcellFile, getFlowcellInfoFromFile(flowcellID)])

###Pre-processing
print("#", " ".join([ quoteStringsWithSpaces(arg) for arg in sys.argv ]), sep="")
print()

#Subset to samples we"re interested in based on command line arguments
flowcellFile = flowcellFile[flowcellFile["Sample Type" ] != "Pool"]
if sampletypes is not None:
    flowcellFile = flowcellFile[flowcellFile["Sample Type"].isin(sampletypes)]
if projects is not None:
    flowcellFile = flowcellFile[flowcellFile["Lab"].isin(projects)]
if people is not None:
    flowcellFile = flowcellFile[flowcellFile["Made By"].isin(people)]
if samples is not None:
    #Do partial string matching rather than exact equality to allow flexible subsetting of certain sublibraries or of all sublibraries for a given BS number
    flowcellFile = flowcellFile[flowcellFile["Sample #"].apply(lambda bs: any([s in bs for s in samples]))]
if samplenames is not None:
    #Do partial string matching rather than exact equality to allow flexible subsetting of certain sublibraries or of all sublibraries for a given BS number
    flowcellFile = flowcellFile[flowcellFile["#Sample Name"].apply(lambda bs: any([s in bs for s in samplenames]))]


#Initialize inputs.txt - must be done before we drop duplicate sample rows
if len(set(flowcellFile["Sample Type"]).intersection(set(bwaPipelineFragmentLengthsMap.keys()))) > 0:
    dirs = set([])
    for index, line in flowcellFile[flowcellFile["Sample Type"].isin(bwaPipelineFragmentLengthsMap.keys())].iterrows():
        dirs.add(getBasedir(line["Lab"], line["Sample Type"], line["FC"]) + "*")
    findCmd =  "find " + " ".join(dirs) + ' -maxdepth 1 -regextype posix-awk -regex "^.+(' + '|'.join(sorted(flowcellFile['Sample #'].unique())) + ').+$"'
    if args.aggregate or args.aggregate_sublibraries:
        findCmd += ' -name "*.bam"'
    else:
        findCmd += ' -name "*.fastq.gz"'
    print(findCmd + " | sort > inputs.txt")
    print()


#Adjust sample IDs and drop duplicate rows to handle aggregations
flowcellFile["Original Sample #"] = flowcellFile["Sample #"]
if args.aggregate or args.aggregate_sublibraries:
    #Sort aggregations by BS number since sample sheet order doesn"t matter
    flowcellFile = flowcellFile.sort_values(by="Sample #")
    if args.aggregate:
        flowcellFile["Sample #"] = flowcellFile["Sample #"].apply(lambda bs: re.sub("(BS[0-9]{5})[A-Z]", "\\1", bs))
    #BWA pipeline just needs last entry per sample
    if len(set(flowcellFile["Sample Type"]).intersection(set(bwaPipelineFragmentLengthsMap.keys()))) > 0:
        flowcellFile = flowcellFile[~flowcellFile.duplicated("Sample #", keep="last")]


#Some aggregations must be processed in groups of samples
if args.aggregate or args.aggregate_sublibraries:
    if len(set(flowcellFile["Sample Type"]).intersection(set(["Transposon DNA", "Transposon RNA", "Transposon 10xRNA", "Transposon iPCR", "Transposon iPCR Capture"]))) > 0:
        for curSample in set(flowcellFile[flowcellFile["Sample Type"].isin(["Transposon DNA", "Transposon RNA", "Transposon 10xRNA", "Transposon iPCR", "Transposon iPCR Capture"])]["Sample #"]):
            print(aggregateTransposonSamples(flowcellFile[flowcellFile["Sample #"]==curSample]))


###Dispatch appropriate function handler per sample line
#Will map to these custom genomes when specified, stored as they appear in LIMS (without cegsvectors_ prefix)
cegsGenomes = [ re.sub(r"^cegsvectors_", "", os.path.basename(x)) for x in glob.glob("/gpfs/data/cegs/sequences/cegsvectors_*") ]
if args.verbose:
    print("\nFound custom references:" + ",".join(cegsGenomes), file=sys.stderr)


doBwaCleanup = False
doDNACaptureCleanup = False
doSARScleanup = False
doTransposonCleanup = False

for index, line in flowcellFile.iterrows():
    if args.verbose:
        print("\nParsing:\n" + str(line), file=sys.stderr)
    sampleName = line["#Sample Name"]
    sampleID = line["Sample #"]
    sampleType = line["Sample Type"]
    
    try:
        if re.compile(r"^BS[0-9]{5}[A-Z]?$").search(sampleID) is None:
            raise Exception("Can't parse " + sampleID + " as BS number")
        
        if sampleType in ["Transposon DNA", "Transposon RNA", "Transposon 10xRNA", "Transposon iPCR", "Transposon iPCR Capture"]:
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
        print("ERROR for sample ", sampleName, ":", e, sep="", file=sys.stderr)
        print("#Couldn't process " + sampleName + "-" + sampleID)


###Cleanup calls
if doBwaCleanup:
    for sampleType in flowcellFile["Sample Type"][flowcellFile["Sample Type"].isin(bwaPipelineAnalysisCommandMap.keys())].unique():
        basedir = getBwaPipelineOutdir(sampleType)
        if args.aggregate or args.aggregate_sublibraries:
            sgeoutput = ""
        else:
            sgeoutput = " -o " + basedir
        
        print()
        print("#" + sampleType)
        print("qsub --time 4:00:00 --mem-per-cpu 4G -j y -b y -N analyzeInserts" + sgeoutput + " -hold_jid `cat " + basedir + "sgeid.analysis | perl -pe 's/\\n/,/g;'` \"/gpfs/data/mauranolab/mapped/src/dnase/analyzeInserts.R " + str(bwaPipelineFragmentLengthsMap[sampleType]) + " " + basedir + "\"")
        print("qsub --time 4:00:00 --mem-per-cpu 4G -j y -b y -N analyzeGC" + sgeoutput + " -hold_jid `cat " + basedir + "sgeid.analysis | perl -pe 's/\\n/,/g;'` \"/gpfs/data/mauranolab/mapped/src/dnase/analyzeGC.R " + basedir + "\"")
        print("qsub --time 4:00:00 --mem-per-cpu 4G -S /bin/bash -j y -N mapped_readcounts" + sgeoutput + " -hold_jid `cat " + basedir + "sgeid.analysis | perl -pe 's/\\n/,/g;'` \"/gpfs/data/mauranolab/mapped/src/dnase/mapped_readcounts.sh " + basedir + "\"")
    #Leave this around so the hold_jid args don't fail
    #    print("rm -f sgeid.analysis")


if doSARScleanup:
    print()
    print("SARS")
    for sampleType in flowcellFile["Sample Type"][flowcellFile["Sample Type"].isin(bwaPipelineAnalysisCommandMap.keys())].unique():
        basedir = getBwaPipelineOutdir(sampleType)
        print("qsub --time 4:00:00 --mem-per-cpu 4G -b y -S /bin/bash -j y -N analyzeSARS -o " + basedir + " -hold_jid `cat " + basedir + "sgeid.analysis | perl -pe 's/\\n/,/g;'` \"/gpfs/data/mauranolab/mapped/src/dnase/analyzeSARS.sh " + basedir + "\"")


#TODO placeholder for now
if doDNACaptureCleanup:
    print()
    basedir=getBwaPipelineOutdir("DNA Capture")
    print("#qsub --time 4:00:00 --mem-per-cpu 4G -j y -b y -S /bin/bash -N analyzeCapture -o " + basedir + " -hold_jid `cat " + basedir + "sgeid.analysis | perl -pe 's/\\n/,/g;'` \"/gpfs/data/mauranolab/mapped/src/dnase/analyzeCapture.sh mappedgenome bait dirs\"")


#BUGBUG doesn't work since we don't have the pid of the final job at submission time
#if doTransposonCleanup:
#    print()
#    print("qsub --time 4:00:00 --mem-per-cpu 4G -j y -S /bin/bash -b y -N transposon_info -hold_jid `cat sgeid.merge | perl -pe 's/\\n/,/g;'` \"/gpfs/data/mauranolab/mapped/src/transposon/Flowcell_Info.sh /home/mauram01/public_html/flowcellInfo/\"")
#    print("rm -f sgeid.merge")
