#!/bin/env python
from trackhub import Track, CompositeTrack, ViewTrack, default_hub, SuperTrack
from trackhub.track import SubGroupDefinition
from sys import argv
import sys
import re
import argparse
import io
from collections import OrderedDict
import locale
import csv
import os
import urllib.parse


locale.setlocale(locale.LC_ALL, 'en_US')


#####################################################################################################
# MakeTrackhub.py requires the forked daler/trackhub from mauranolab. It can be installed locally:
#     pip install --upgrade --user .
#
# Or directly from git:
#     pip install --upgrade --user git+git://github.com/mauranolab/trackhub.git
#####################################################################################################

parser = argparse.ArgumentParser(prog = "MakeTrackhub.py", description = "Creates a trackhub with composite and subtracks for all tracks.", add_help = True)
parser.add_argument('Input', action = 'store', help = 'Input tsv file with sample metadata. Name, DS, Replicate, Color, Assay, analyzed_reads, SPOT, Num_hotspots, Exclude, Genomic_coverage, Group, Age, filebase')
#NB Exclude not used right now 
parser.add_argument('--genome', action = 'store', required = True, help = 'genome assembly name')
parser.add_argument('--URLbase', action = 'store', required = False, default="", help = 'URL base at which track files are hosted')
parser.add_argument('--includeSampleIDinSampleCol', action = 'store_true', default = False, help = 'Append the Sample ID in the Sample column')
parser.add_argument('--checksamples', action = 'store_true', default = False, help = 'Mark Density and Coverage tracks for display by turning on composite track without going to configuration page (NB ChIP-seq input samples are never displayed by default)')
parser.add_argument('--supertrack', action = 'store', required = False, help = 'Encompass all composite tracks generated within supertrack. Supertrack name specified as parameter.')
parser.add_argument('--supertrackPriority', type=float, action = 'store', required = False, help = 'Priority values set the order in which supertracks are displayed in the browser control panel.')
parser.add_argument('--generateHTMLdescription', action = 'store_true', default = False, help = 'Link to HTML descriptions for composite tracks, assumed to be present at [genome]/descriptions/[group name].html')
parser.add_argument('--tracknameprefix', action = 'store', required = False, default="", help = 'Add prefix within track names (e.g. to permit unique names).')
parser.add_argument('--subgroupnames', action = 'store', required = False, default="SampleID", help = 'Comma-separated list of columns in input to add as subgroups; will be sorted in order specified')
parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")


try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument, file=sys.stderr)
    sys.exit(2)

# Used to keep track of what track names have been assigned to samples, in case they appear in multiple flow cells.
sampleName_dict = dict()

#track name must contain only the following chars: [a-zA-Z0-9_]
def cleanTrackName(x):
    x = re.sub(r'\W+', '', x)
    x = re.sub(r'-', '_', x)
    x = re.sub(r'[^a-zA-Z0-9_]', '', x)
    return x


#Helper functions for managing subGroups
def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

def internalSubgroupName(label):
    return label.lower().replace(' ', '_').replace('=', '')

def cleanFactorForSubGroup(label):
    return label.replace(' ', '_').replace('~', '_').replace('=', '').replace('+', '').replace('[', '').replace(']', '').replace('(', '').replace(')', '')

def createSubGroup(subGroups, subGroupNames, keys, label):
    name = internalSubgroupName(label)
    
    #What does natural_key do?
    uniqkeys = set(sorted([cleanFactorForSubGroup(key) for key in keys], key=natural_key))
    if uniqkeys and not all(x == 'NA' for x in uniqkeys):
        if len(uniqkeys) > 1000:
            print("[MakeTrackhub.py] WARNING Subgroup ", label, " exceeds maximum 1000 members (has ", len(uniqkeys), ")", sep="", file=sys.stderr)
        
        subGroupDict = dict(zip(uniqkeys, uniqkeys))
        subGroupDef = SubGroupDefinition(
            # The "label" below needs to have the underscore in it, else only the part up to the first
            # space will be displayed. However, the underscore will not be displayed in the browser. 
            # Kent converts the underscores to spaces for display purposes.
            # Hyphens seem to be retained though.
            name=name,
            label=label,
            mapping=OrderedDict(sorted(subGroupDict.items(), key=lambda x: x[1]))
        )
        subGroups.append(subGroupDef)
        subGroupNames[label] = len(uniqkeys)


#For an array, return a dict mapping the elements of the array to their shortest unique prefix
#https://stackoverflow.com/questions/22468435/how-would-i-look-for-the-shortest-unique-subsequence-from-a-set-of-words-in-pyth
def shortest_unique_strings(array, minlength=1):
    ret = {}
    for k in set(array):
        #Initialize with the full element in case k is fully a prefix of another
        ret[k] = k
        for ix in range(min(len(k), minlength), len(k)):
            #When array has only 1 element, return identity map rather than truncating element
            if len([key for key in array if key.startswith(k[:ix])]) == 1 and len(array) > 1:
                ret[k] = k[:ix]
                break
    return ret


########
#import metadata
########
if args.Input=="-":
    inputfile = sys.stdin
else:
    inputfile = open(args.Input, 'r') 
#inputfile = open('/vol/isg/encode/dnase/SamplesForTrackhub.tsv', 'r') 

inputfile_reader = csv.DictReader(inputfile, delimiter='\t')
input_data_all = []
for line in inputfile_reader:
    input_data_all.append(line)


# If new Assay types are added in samplesforTrackhub.R, then changes must also be made
# in the following for loop, as well as in the various track construction blocks below.
assays = set([])
for assay_type in set([line.get('Assay') for line in input_data_all]):
    if assay_type in["DNase-seq", "DNA", "DNA Capture", "Amplicon", "None", "ATAC-seq"]:
        assays.add(assay_type)
    else:
        #ChIP-seq tracks have the epitope as the assay type
        #NB this will pick up any samples whose type doesn't match the above
        assays.add("ChIP-seq")

if args.checksamples:
    DensCovTracksDefaultDisplayMode="on"
else:
    DensCovTracksDefaultDisplayMode="off"

if args.subgroupnames is not None:
    customSubGroupNames = args.subgroupnames.split(',')
else:
    customSubGroupNames = []


########
#Create hub file
########
#The actual values here don't get used since we just print the trackDb file but are required
hub, genomes_file, genome, trackdb = default_hub(
    hub_name="ENCODE_DNase",
    genome=args.genome,
    short_label="ENCODE DNase-seq",
    long_label="Encode DNase-seq",
    email="maurano@nyu.edu")


# Initialize the supertrack
# "supertrackonoff" requires the forked daler/trackhub
if args.supertrack is not None:
    supertrack = SuperTrack(
        name=cleanTrackName(args.supertrack),
        short_label=args.supertrack,
        long_label=args.supertrack,
        supertrackonoff="show")
    if args.genome in ["cegsvectors", "t2t"]:
        supertrack.add_params(group=args.genome)
    if args.supertrackPriority is not None:
        #NB: trackhub requires that priority be a float or it will silently drop the line
        supertrack.add_params(priority=args.supertrackPriority)
    trackdb.add_tracks(supertrack)


# Now build the composites, and add them all to the supertrack.
for assay_type in assays:
    #Matches the options allowed in submit.sh and createFlowcellSubmit.py
    bwaPipelineAnalysisCommandMap = { "DNase-seq": "dnase", "Nano-DNase": "dnase", "ATAC-seq":"atac", "ChIP-seq": "chipseq", "DNA": "dna", "DNA Capture": "capture", "Amplicon": "amplicon" }
    assay_suffix = bwaPipelineAnalysisCommandMap[assay_type] if assay_type in bwaPipelineAnalysisCommandMap else "none"
    
    input_data = []
    for line in input_data_all:
        if line['Assay'] == assay_type or assay_type == "ChIP-seq" and line['Assay'] not in assays:
            input_data.append(line)
        else:
            continue
    
    if args.verbose:
        print("[MakeTrackhub.py] Processing assay", assay_type, file=sys.stderr)
    
    ########
    # Create composite track for major groups of samples (e.g. Duke, UW, etc.)
    # These are flowcell IDs for assay_type == "DNA"
    ########
    #Find all unique sampleNames and create a supertrack for each
    Group = sorted(set([line['Group'] for line in input_data]), reverse=True)
    
    for curGroup in Group:
        ########
        #Create subgroup definitions
        ########
        matchingSamples = [line for line in input_data if curGroup == line['Group']]
        
        if args.verbose:
            print("[MakeTrackhub.py] Processing", len(matchingSamples), "samples in group", curGroup, file=sys.stderr)
        
        subGroupDefs = []
        #Dict containing the active subgroup labels and the number of unique levels
        subGroupNames = dict()
        
        createSubGroup(subGroupDefs, subGroupNames, [re.sub(r'_L$|_R$', '', line['Name']) + ('-' + line['SampleID'] if args.includeSampleIDinSampleCol else '') for line in matchingSamples], "Sample")
        
        if not args.includeSampleIDinSampleCol:
            createSubGroup(subGroupDefs, subGroupNames, [line['SampleID'] for line in matchingSamples], 'SampleID')
        
        if assay_type == "ChIP-seq":
            createSubGroup(subGroupDefs, subGroupNames, [line['Assay'] for line in matchingSamples], 'Assay')
        
        for subGroupName in customSubGroupNames:
            createSubGroup(subGroupDefs, subGroupNames, [line[subGroupName] for line in matchingSamples], subGroupName)
        
        if args.genome == "cegsvectors":
            #Add cegsvectors to make sure a minimum label is printed if there is only a single Mapped_Genome
            cegsvectorsAbbreviatedNames = shortest_unique_strings(set(["cegsvectors"] + [line['Mapped_Genome'] for line in matchingSamples]), minlength=18)
            cegsvectorsAbbreviatedNames = {k: re.sub(r'^cegsvectors_', '', v) for k, v in cegsvectorsAbbreviatedNames.items()}
            createSubGroup(subGroupDefs, subGroupNames, [re.sub(r'^cegsvectors_', '', line['Mapped_Genome']) for line in matchingSamples], "Mapped_Genome")
        
        
        # SortOrder controls what is displayed, even if we retain all the subgroup definitions.
        SortOrder = ""
        for subGroupLabel in subGroupNames.keys():
            SortOrder += internalSubgroupName(subGroupLabel) + "=+ "
        
        
        curGroup_trackname = cleanTrackName(args.genome + args.tracknameprefix + "_" + curGroup + "_" + assay_suffix)
        # short labels are supposed to be 17 chars max. However, the browser does not seem to care.
        shrt_label = curGroup + '_' + assay_type
        # long labels are supposed to be 76 chars max. However, the browser does not seem to care.
        lng_label = curGroup + '_' +  assay_type
        
        composite = CompositeTrack(
            name=curGroup_trackname,
            short_label=shrt_label,
            long_label=lng_label,
            tracktype="bigBed",
            visibility="hide",
            parentonoff="off",
            sortOrder=SortOrder,
            autoScale = "off",
            bigDataUrl ='NULL',
            maxHeightPixels= "100:32:8")
        composite.add_subgroups(subGroupDefs)
        
        if args.generateHTMLdescription:
            if args.supertrack == "Flowcells":
                composite.add_params(html='descriptions/' + curGroup + '_' + assay_suffix + '.html')
            else:
                composite.add_params(html='descriptions/' + curGroup + '.html')
        
        params_dimensions = ""
        if assay_type == "ChIP-seq":
            if subGroupNames['Assay'] < subGroupNames['Sample']:
                params_dimensions="dimX=assay dimY=sample"
            else:
                params_dimensions="dimY=assay dimX=sample"
        else:
            params_dimensions="dimY=sample"
            if args.genome == "cegsvectors":
                params_dimensions = params_dimensions + " dimX=mapped_genome"
            elif 'Age' in subGroupNames:
                params_dimensions = params_dimensions + " dimX=age"
        
        #hardcoded for CEGS
        if 'Type' in subGroupNames:
            params_dimensions = params_dimensions + " dimA=type"
            
        if 'Replicate' in subGroupNames:
            params_dimensions = params_dimensions + " dimA=replicate"
            #NB requires mauranolab fork of daler/trackhub
            #Likely also compatible with pull request to revamp daler parameter support: https://github.com/daler/trackhub/pull/33
            composite.add_params(dimensionAchecked="rep1")
        composite.add_params(dimensions=params_dimensions)
        
        
        #Placeholders for view tracks, which will get initialized upon the first track below
        Reads_view = None
        Coverage_view = None
        Dens_view = None
        Hotspots_view = None
        Peaks_view = None
        Cuts_view = None
        Variants_view = None
        Delly_view = None
        Genotypes_view = None
        
        
        ########
        #Make tracks for all bigWigs in current dir
        ########
        for curSample in matchingSamples:
            sampleName = re.sub(r'_L$|_R$', '', curSample['Name'])
            
            
            ###Internal track name (never displayed) must be unique within the Genome Browser or dataHub and must begin with a letter
            if args.genome == "cegsvectors":
                sampleNameGenome = cegsvectorsAbbreviatedNames[curSample['Mapped_Genome']]
            else:
                if 'Annotation_Genome' in curSample and curSample['Annotation_Genome'] != "NA":
                    sampleNameGenome = curSample['Annotation_Genome']
                else:
                    sampleNameGenome = args.genome
            
            #There is a length limit of 128 characters, but in practice some are used for an internal prefix/suffix
            # We add 4 characters to sampleName_trackname later in the code ('_bam', '_cov', etc.)
            # The browser prepends 'hub_161_', and appends '.priority'
            # So sampleName_trackname is increased by 4+8+9=21 characters.
            # If sampleName_trackname starts with 107 characters, then 128 characters get sent to the server.  This causes an error.
            # So sampleName_trackname needs to be 106 characters or less. 
            sampleName_trackname = cleanTrackName(sampleNameGenome + "_" + curGroup + "_" + curSample['SampleID'])
            
            
            # Make sure there are no duplicate track names.
            if sampleName_trackname in sampleName_dict:
                if args.verbose:
                    print("[MakeTrackhub.py] ", sampleName_trackname, " duplicates an existing trackname", sep="", file=sys.stderr)
                sampleName_trackname += "_"
                i = 2
                while sampleName_trackname + str(i) in sampleName_dict:
                    i += 1
                sampleName_trackname += str(i)
            sampleName_dict[sampleName_trackname] = sampleName_trackname
            
            
            ###Create long label from analyzed reads, SPOTS and Hotspots. Format numbers into 1000's - longLabel must be <= 76 printable characters
            sampleDescription = sampleName + "-"
            if assay_type == "ChIP-seq":
                sampleDescription += curSample['Assay'] + "-"
            sampleDescription += curSample['SampleID']
            if curSample['analyzed_reads'] != "NA":
                sampleDescription += ' (' + locale.format_string("%d", int(curSample['analyzed_reads']), grouping=True) + ' analyzed reads, '
            if curSample['Genomic_coverage'] != "NA":
                sampleDescription += curSample['Genomic_coverage'] + 'x coverage)'
            else:
                if curSample['Num_hotspots'] == "NA":
                    sampleDescriptionNumHotspots = "no"
                else:
                    sampleDescriptionNumHotspots = locale.format_string("%d", int(curSample['Num_hotspots']), grouping=True)
                if curSample['SPOT'] == "NA":
                    sampleDescriptionSPOT = "NA"
                else:
                    sampleDescriptionSPOT = str(round(float(curSample['SPOT']), 2))
                sampleDescription += sampleDescriptionSPOT + ' SPOT, ' + sampleDescriptionNumHotspots + ' Hotspots)'
            
            if 'Age' in subGroupNames:
                sampleDescription += (', Age = ' + curSample['Age'] if curSample['Age'] != 'NA' else '')
            
            if 'Bait_set' in curSample and curSample['Bait_set'] != 'NA':
                #Remove "_bait" suffix formerly used in individual bait names to save space
                sampleDescription += ', Bait=' + ','.join([ re.sub('_bait$', '', bait) for bait in curSample['Bait_set'].split(",") ])
            
            if args.genome == "cegsvectors" and 'Mapped_Genome' in curSample and curSample['Mapped_Genome'] != 'NA':
                sampleDescription += ", Mapped to " + re.sub(r'^cegsvectors_', '', curSample['Mapped_Genome'])
            
            ####Set up subgroups
            #BUGBUG get error "Subgroup SampleID exceeds maximum 1000 members" for some of the ENCODE tracks
            sampleSubgroups = dict(sample=cleanFactorForSubGroup(sampleName + ('-'+curSample['SampleID'] if args.includeSampleIDinSampleCol else '') ))
            for subGroupLabel in ['SampleID', 'Assay'] + customSubGroupNames:
                if subGroupLabel in subGroupNames:
                    sampleSubgroups[internalSubgroupName(subGroupLabel)] = cleanFactorForSubGroup(curSample[subGroupLabel])
            
            #Keep these as short as possible to avoid running over track name maxlength
            if 'Mapped_Genome' in subGroupNames:
                sampleSubgroups['mapped_genome'] = cleanFactorForSubGroup(re.sub(r'^cegsvectors_', '', curSample['Mapped_Genome']))
            
            
            ###shortLabel must be <= 17 printable characters (sampleName + DS)
            sampleShortLabel = sampleName + "-" + curSample['SampleID']
            
            
            ###
            #adding viewUi="off" seems to expand the panel by default, opposite expectation
            if Reads_view is None:
                Reads_view = ViewTrack(
                    name="Reads_view_" + curGroup_trackname,
                    view="Reads",
                    visibility="hide",
                    parentonoff="off",
                    tracktype="bam",
                    short_label="Reads",
                    maxItems=50000, #Set high so that dense sequencing tracks can be displayed as this parameter can not be changed in the UI
        #            maxWindowToDraw=10000,
                    pairEndsByName=".",
                    minAliQual=20,
                    long_label="Reads")
                composite.add_view(Reads_view)
            track = Track(
                name=sampleName_trackname + '_bam',
                short_label=sampleShortLabel,
                long_label=assay_type + ' Reads ' + sampleDescription,
                url=args.URLbase + urllib.parse.quote(curSample['filebase']) + '.bam',
                subgroups=sampleSubgroups,
                tracktype='bam',
                parentonoff="off"
            )
            Reads_view.add_tracks(track)
            
            if curSample['analyzed_reads']!="NA":
                #Coverage_view
                if assay_type in ["DNA", "DNA Capture", "Amplicon"]:
                    if Coverage_view is None:
                        Coverage_view = ViewTrack(
                            name="Coverage_view_" + curGroup_trackname,
                            view="Coverage",
                            visibility="full",
                            parentonoff=DensCovTracksDefaultDisplayMode,
                            tracktype="bigWig",
                            viewLimits="0:500", #Keep high since it becomes a hard limit in the UI
                            autoScale='on',
                            alwaysZero='on',
                            maxHeightPixels="100:30:10",
                            short_label="Coverage",
                            long_label="Coverage",
                            windowingFunction="mean")
                        composite.add_view(Coverage_view)
                    track = Track(
                        name=sampleName_trackname + '_cov',
                        short_label=sampleShortLabel,
                        long_label=assay_type + ' Coverage ' + sampleDescription,
                        url=args.URLbase + urllib.parse.quote(curSample['filebase']) + '.coverage.bw',
                        subgroups=sampleSubgroups,
                        tracktype='bigWig',
                        color=curSample['Color'],
                        parentonoff=DensCovTracksDefaultDisplayMode
                    )
                    Coverage_view.add_tracks(track)
                
                #Dens_view
                if assay_type in ["DNase-seq", "ChIP-seq", "ATAC-seq"]:
                    if Dens_view is None:
                        Dens_view = ViewTrack(
                            name="Dens_view_" + curGroup_trackname,
                            view="Density",
                            visibility="full",
                            parentonoff=DensCovTracksDefaultDisplayMode,
                            tracktype="bigWig",
                            viewLimits="0:3",
                            autoScale="off",
                            maxHeightPixels="100:30:10",
                            short_label="Density",
                            long_label="Density",
                            windowingFunction="maximum")
                        composite.add_view(Dens_view)
                    track = Track(
                        name=sampleName_trackname + '_dens',
                        short_label=sampleShortLabel,
                        long_label=assay_type + ' Density ' + sampleDescription,
                        url=args.URLbase + urllib.parse.quote(curSample['filebase']) + '.density.bw',
                        subgroups=sampleSubgroups, 
                        tracktype='bigWig',
                        color=curSample['Color'],
                        parentonoff="off" if assay_type=='ChIP-seq' and curSample['Assay']=='input' else DensCovTracksDefaultDisplayMode
                    )
                    Dens_view.add_tracks(track)
                
                #Hotspots_view
                if assay_type in ["DNase-seq", "ChIP-seq", "ATAC-seq"]:
                    hotspot_base = os.path.basename(curSample['filebase'])
                    hotspot_dir = os.path.dirname(curSample['filebase'])
                    hotspot_path = hotspot_dir + "/hotspot2/" + hotspot_base
                    
                    if Hotspots_view is None:
                        Hotspots_view = ViewTrack(
                            name="Hotspots_view_" + curGroup_trackname,
                            view="Hotspots",
                            visibility="hide",
                            parentonoff="off",
                            tracktype="bigBed",
                            maxItems=250,
                            short_label="Hotspots",
                            long_label="Hotspots")
                        composite.add_view(Hotspots_view)
                    track = Track(
                        name=sampleName_trackname + '_hot',
                        short_label=sampleShortLabel,
                        long_label=assay_type + ' Hotspots (5% FDR) ' + sampleDescription,
                        url=args.URLbase + hotspot_path + '.hotspots.fdr0.05.bb',
                        subgroups=sampleSubgroups,
                        tracktype='bigBed 4',
                        color=curSample['Color'],
                        parentonoff="off"
                    )
                    Hotspots_view.add_tracks(track)
                    
                    #Peaks_View
                    if Peaks_view is None:
                        Peaks_view = ViewTrack(
                            name="Peaks_view_" + curGroup_trackname,
                            view="Peaks",
                            visibility="hide",
                            parentonoff="off",
                            tracktype="bigBed",
                            maxItems=250,
                            short_label="Peaks",
                            long_label="Peaks")
                        composite.add_view(Peaks_view)
                    track = Track(
                        name=sampleName_trackname + '_pks',
                        short_label=sampleShortLabel,
                        long_label=assay_type + ' Peaks ' + sampleDescription,
                        url=args.URLbase + hotspot_path  + '.peaks.bb',
                        subgroups=sampleSubgroups,
                        tracktype='bigBed 3',
                        color=curSample['Color'],
                        parentonoff="off"
                    )
                    Peaks_view.add_tracks(track)
                
                #PerBaseCutCount
                if assay_type in ["DNase-seq", "ATAC-seq"]:
                    if Cuts_view is None:
                        Cuts_view = ViewTrack(
                            name="Cuts_view_" + curGroup_trackname,
                            view="Cuts",
                            visibility="hide",
                            parentonoff="off",
                            tracktype="bigWig",
                            viewLimits="0:2",
                            autoScale="off",
                            maxHeightPixels="100:30:10",
                            short_label="Cut Counts",
                            long_label="Cut Counts",
                            windowingFunction="maximum")
                        composite.add_view(Cuts_view)
                    track = Track(
                        name=sampleName_trackname + '_cuts',
                        short_label=sampleShortLabel,
                        long_label=assay_type + ' Cut counts ' + sampleDescription,
                        url=args.URLbase + urllib.parse.quote(curSample['filebase']) + '.perBase.bw',
                        subgroups=sampleSubgroups,
                        tracktype='bigWig',
                        color=curSample['Color'],
                        parentonoff="off"
                    )
                    Cuts_view.add_tracks(track)
                
                #Variants_view
                if assay_type in ["DNA", "DNA Capture", "Amplicon"]:
                    if Variants_view is None:
                        Variants_view = ViewTrack(
                            name="Variants_view_" + curGroup_trackname,
                            view="Variants",
                            visibility="hide",
                            parentonoff="off",
                            tracktype='vcfTabix',
                            maxItems=250,
                            applyMinQual="true",
                            minQual=10,
                            short_label="Variants",
                            long_label="Variants")
                        composite.add_view(Variants_view)
                    track = Track(
                        name=sampleName_trackname + '_vcf',
                        short_label=sampleShortLabel,
                        long_label=assay_type + ' Variants ' + sampleDescription,
                        url=args.URLbase + urllib.parse.quote(curSample['filebase']) + '.filtered.vcf.gz',
                        subgroups=sampleSubgroups,
                        tracktype='vcfTabix',
                        parentonoff="off"
                    )
                    Variants_view.add_tracks(track)
                
                #Genotypes_view
                if assay_type in ["DNA", "DNA Capture", "Amplicon"]:
                    if Genotypes_view is None:
                        Genotypes_view = ViewTrack(
                            name="Genotypes_view_" + curGroup_trackname,
                            view="Genotypes",
                            visibility="dense",
                            parentonoff="off",
                            tracktype="bigBed",
                            maxItems=100000,
                            short_label="Genotypes",
                            long_label="Genotypes")
                        composite.add_view(Genotypes_view)
                    track = Track(
                        name=sampleName_trackname + '_gts',
                        short_label=sampleShortLabel,
                        long_label=assay_type + ' Genotypes ' + sampleDescription,
                        url=args.URLbase + urllib.parse.quote(curSample['filebase']) + '.genotypes.bb',
                        subgroups=sampleSubgroups,
                        tracktype='bigBed 9 .',
                        color=curSample['Color'],
                        itemRgb="on",
                        parentonoff="off"
                    )
                    Genotypes_view.add_tracks(track)
                
                #Delly_view
                if assay_type in ["DNA", "DNA Capture", "Amplicon"]:
                    if Delly_view is None:
                        Delly_view = ViewTrack(
                            name="Delly_view_" + curGroup_trackname,
                            view="DELLY",
                            visibility="pack",
                            parentonoff="off",
                            tracktype='vcfTabix',
                            maxItems=250,
                            applyMinQual="true",
                            minQual=0,
                            short_label="DELLY",
                            long_label="DELLY")
                        composite.add_view(Delly_view)
                    track = Track(
                        name=sampleName_trackname + '_delly',
                        short_label=sampleShortLabel,
                        long_label=assay_type + ' DELLY ' + sampleDescription,
                        url=args.URLbase + urllib.parse.quote(curSample['filebase']) + '.DELLY.filtered.vcf.gz',
                        subgroups=sampleSubgroups,
                        tracktype='vcfTabix',
                        parentonoff="off"
                    )
                    Delly_view.add_tracks(track)
        
        for view in composite.views:
            view.add_params(configurable="on")
        # Demonstrate some post-creation adjustments...here, just make control samples gray
        #for track in trackdb.tracks:
        #    if 'control' in track.name:
        #     track.add_params(color="100,100,100")
        
        if args.supertrack is not None:
            supertrack.add_tracks(composite)
        else:
            trackdb.add_tracks(composite)

print(trackdb)
print("\n")
