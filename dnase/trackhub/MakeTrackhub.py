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
parser.add_argument('--URLbase', action = 'store', required = True, help = 'URL base at which track files are hosted')
parser.add_argument('--includeSampleIDinSampleCol', action = 'store_true', default = False, help = 'Append the Sample ID in the Sample column')
parser.add_argument('--checksamples', action = 'store_true', default = False, help = 'Mark Density and Coverage tracks for display by turning on composite track without going to configuration page')
parser.add_argument('--supertrack', action = 'store', required = False, help = 'Encompass all composite tracks generated within supertrack. Supertrack name specified as parameter.')
parser.add_argument('--generateHTMLdescription', action = 'store_true', default = False, help = 'Link to HTML descriptions for composite tracks, assumed to be present at [genome]/descriptions/[group name].html')
parser.add_argument('--tracknameprefix', action = 'store', required = False, default="", help = 'Add prefix within track names (e.g. to permit unique names).')

try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument, file=sys.stderr)
    sys.exit(2)


#track name must contain only the following chars: [a-zA-Z0-9_]
def cleanTrackName(x):
    x = re.sub(r'\W+', '', x)
    x = re.sub(r'-', '_', x)
    x = re.sub(r'[^a-zA-Z0-9_]', '', x)
    return x


#Helper functions for managing subGroups
def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

def createSubGroup(subGroups, subGroupNames, keys, name, label):
    uniqkeys = set(keys)
    if uniqkeys and not all(x == 'NA' for x in uniqkeys):
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
        subGroupNames[name] = len(uniqkeys)


#For an array, return a dict mapping the elements of the array to their shortest unique prefix
#https://stackoverflow.com/questions/22468435/how-would-i-look-for-the-shortest-unique-subsequence-from-a-set-of-words-in-pyth
def shortest_unique_strings(array):
    ret = {}
    for k in set(array):
        for ix in range(len(k)):
            #Initialize with the full element in case k is fully a prefix of another
            ret[k] = k
            #When array has only 1 element, return identity map rather than truncating element
            if len([key for key in array if key.startswith(k[:ix])]) == 1 and len(array)>1:
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
    if assay_type in["DNase-seq", "DNA", "Capture", "None", "ATAC-seq"]:
        assays.add(assay_type)
    else:
        #ChIP-seq tracks have the epitope as the assay type
        #NB this will pick up any samples whose type doesn't match the above
        assays.add("ChIP-seq")

if args.checksamples:
    DensCovTracksDefaultDisplayMode="on"
else:
    DensCovTracksDefaultDisplayMode="off"

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
    if args.genome == "cegsvectors":
        supertrack.add_params(group="cegsvectors")
    trackdb.add_tracks(supertrack)


# Now build the composites, and add them all to the supertrack.
for assay_type in assays:
    #Matches the options allowed in submit.sh and createFlowcellSubmit.py
    bwaPipelineAnalysisCommandMap = { "DNase-seq": "dnase", "Nano-DNase": "dnase", "ATAC-seq":"atac", "ChIP-seq": "chipseq", "DNA": "dna", "Capture": "capture" }
    assay_suffix = bwaPipelineAnalysisCommandMap[assay_type] if assay_type in bwaPipelineAnalysisCommandMap else "none"
    
    input_data = []
    for line in input_data_all:
        if line['Assay'] == assay_type or assay_type == "ChIP-seq" and line['Assay'] not in assays:
            input_data.append(line)
        else:
            continue
    
    
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
        
        subGroupDefs = []
        #Dict containing the active subgroup internal names and the number of unique levels
        subGroupNames = dict()
        createSubGroup(subGroupDefs, subGroupNames, sorted([re.sub(r'_L$|_R$', '', line['Name']) + ('-' + line['DS'] if args.includeSampleIDinSampleCol else '') for line in matchingSamples]), "sampleName", "Sample")
        if not args.includeSampleIDinSampleCol:
            createSubGroup(subGroupDefs, subGroupNames, sorted([line['DS'] for line in matchingSamples]), "DS", "Sample_ID")
        createSubGroup(subGroupDefs, subGroupNames, sorted([line['Replicate'] for line in matchingSamples]), "replicate", "Replicate")
        createSubGroup(subGroupDefs, subGroupNames, sorted([line['Assay'] for line in matchingSamples]), "assay", "Assay")
        createSubGroup(subGroupDefs, subGroupNames, sorted([re.sub(' ', '_', line['Age']) for line in matchingSamples], key=natural_key), "age", "Age")
        if args.genome == "cegsvectors":
            #Add cegsvectors to make sure a minimum label is printed if there is only a single Mapped_Genome
            cegsvectorsAbbreviatedNames = shortest_unique_strings(set(["cegsvectors"] + [line['Mapped_Genome'] for line in matchingSamples]))
            #cegsvectorsAbbreviatedNames = {k: re.sub(r'^cegsvectors_', '', v) for k, v in cegsvectorsAbbreviatedNames.items()}
            createSubGroup(subGroupDefs, subGroupNames, sorted([re.sub(r'^cegsvectors_', '', line['Mapped_Genome']) for line in matchingSamples]), "mappedgenome", "Mapped_Genome")
        
        # SortOrder controls what is displayed, even if we retain all the subgroup definitions.
        SortOrder = "sampleName=+"
        if 'DS' in subGroupNames:
            SortOrder = SortOrder + " DS=+"
        if 'replicate' in subGroupNames:
            SortOrder = SortOrder + " replicate=+ "
        if 'age' in subGroupNames:
            SortOrder = SortOrder + " age=+ "
        SortOrder = SortOrder + " view=+"
        
        # Adding suffix
        curGroup_trackname = cleanTrackName(args.genome + args.tracknameprefix + "_" + curGroup + "_" + assay_suffix)
        
        # Do something like this later.
        # mydate = re.split(r'_', curGroup)[0]
        # priority=bignumber - mydate
        # Why is priority="2" not showing up in the output?
        # Should be a number anyway.
        
        # short labels are supposed to be 17 chars max. However, the browser does not seem to care.
        shrt_label = curGroup + '_' + assay_type
        # shrt_label = shrt_label[0:17]
        # long labels are supposed to be 76 chars max. However, the browser does not seem to care.
        lng_label = curGroup + '_' +  assay_type
        # lng_label = lng_label[0:76]
        
        composite = CompositeTrack(
            name=curGroup_trackname,
            short_label=shrt_label,
            long_label=lng_label,
            tracktype="bigBed",
            priority="2",
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
        if assay_type in ["DNase-seq", "DNA", "Capture"]:
            params_dimensions="dimY=sampleName"
            if args.genome == "cegsvectors":
                params_dimensions = params_dimensions + " dimX=mappedgenome"
            elif 'age' in subGroupNames:
                params_dimensions = params_dimensions + " dimX=age"
        else:
            # ChIP-seq
            if subGroupNames['assay'] < subGroupNames['sampleName']:
                params_dimensions="dimX=assay dimY=sampleName"
            else:
                params_dimensions="dimY=assay dimX=sampleName"

        if 'replicate' in subGroupNames:
            params_dimensions = params_dimensions + " dimA=replicate"
            #NB requires mauranolab fork of daler/trackhub
            #Likely also compatible with pull request to revamp daler parameter support: https://github.com/daler/trackhub/pull/33
            composite.add_params(dimensionAchecked="rep1")
        composite.add_params(dimensions=params_dimensions)
        
        
        ########
        #Create view track
        ########
        #Notes
        #adding viewUi="off" seems to expand the panel by default, opposite expectation
        
        Reads_view = ViewTrack(
            name="Reads_view_" + curGroup_trackname,
            view="Reads",
            visibility="hide",
            parentonoff="off",
            tracktype="bam",
            short_label="Reads",
            maxItems=10000, #Set high so that dense sequencing tracks can be displayed as this parameter can not be changed in the UI
#            maxWindowToDraw=10000,
            pairEndsByName="on",
            long_label="Reads")
        composite.add_view(Reads_view)
        
        if assay_type in ["DNA", "Capture"]:
            Coverage_view = ViewTrack(
                name="Coverage_view_" + curGroup_trackname,
                view="Coverage",
                visibility="full",
                parentonoff=DensCovTracksDefaultDisplayMode,
                tracktype="bigWig",
                viewLimits="0:500",
                autoScale='on',
                maxHeightPixels="100:30:10",
                short_label="Coverage",
                long_label="Coverage")
            composite.add_view(Coverage_view)
        
        if assay_type in ["DNase-seq", "ChIP-seq", "ATAC-seq"]:
            Dens_view = ViewTrack(
                name="Dens_view_" + curGroup_trackname,
                view="Density",
                visibility="full",
                parentonoff=DensCovTracksDefaultDisplayMode,
                tracktype="bigWig",
                viewLimits="0:3",
                autoScale='off',
                maxHeightPixels="100:30:10",
                short_label="Density",
                long_label="Density")
            composite.add_view(Dens_view)
            
        if assay_type in ["DNase-seq", "ChIP-seq", "ATAC-seq"]:
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
        
        if assay_type in ["DNase-seq", "ATAC-seq"]:
            Cuts_view = ViewTrack(
                name="Cuts_view_" + curGroup_trackname,
                view="Cuts",
                visibility="hide",
                parentonoff="off",
                tracktype="bigWig",
                viewLimits="0:2",
                autoScale='off',
                maxHeightPixels="100:30:10",
                short_label="Cut Counts",
                long_label="Cut Counts")
            composite.add_view(Cuts_view)
        
        if assay_type in ["DNA", "Capture"]:
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
        
        if assay_type in ["DNA", "Capture"]:
            Genotypes_view = ViewTrack(
                name="Genotypes_view_" + curGroup_trackname,
                view="Genotypes",
                visibility="hide",
                parentonoff="off",
                tracktype="bigBed",
                maxItems=250,
                short_label="Genotypes",
                long_label="Genotypes")
            composite.add_view(Genotypes_view)
        
        
        for view in composite.views:
            view.add_params(configurable="on")
        
        
        ########
        #Make tracks for all bigWigs in current dir
        ########
        #Create long label from analyzed reads, SPOTS and Hotspots. Format numbers into 1000's
        locale.setlocale(locale.LC_ALL, 'en_US')
        for curSample in matchingSamples:
            sampleName = re.sub(r'_L$|_R$', '', curSample['Name'])
            
            
            ###Internal track name (never displayed) must be unique within the Genome Browser or dataHub and must begin with a letter
            if args.genome == "cegsvectors":
                sampleNameGenome = cegsvectorsAbbreviatedNames[curSample['Mapped_Genome']]
            else:
                if 'Annotation_Genome' in curSample:
                    sampleNameGenome = curSample['Annotation_Genome']
                else:
                    sampleNameGenome = args.genome
            #There is a length limit of 128 characters, but in practice some are used for an internal prefix/suffix
            #Probably could omit sampleName here to save space
            sampleName_trackname = cleanTrackName(sampleNameGenome + "_" + curGroup + "_" + sampleName + "_" + curSample['DS'])
            
            
            ###Set up description - longLabel must be <= 76 printable characters
            sampleDescription = sampleName + "-"
            if assay_type == "ChIP-seq":
                sampleDescription += curSample['Assay']
            sampleDescription += curSample['DS'] + ' (' + locale.format("%d", int(curSample['analyzed_reads']), grouping=True) + ' analyzed reads, '
            if assay_type in ["DNA", "Capture"]:
                sampleDescription += curSample['Genomic_coverage'] + 'x genomic coverage)'
            else:
                if curSample['Num_hotspots'] == "NA":
                    sampleDescriptionNumHotspots = "no"
                else:
                    sampleDescriptionNumHotspots = locale.format("%d", int(curSample['Num_hotspots']), grouping=True) 
                sampleDescription += curSample['SPOT'] + ' SPOT, ' + sampleDescriptionNumHotspots + ' Hotspots)'
            
            if 'age' in subGroupNames:
                sampleDescription = sampleDescription + (', Age = ' + curSample['Age'] if curSample['Age'] != 'NA' else '')
            
            if 'Bait_set' in curSample and curSample['Bait_set'] != 'NA':
                sampleDescription = sampleDescription + ', Bait = ' + curSample['Bait_set']
            
            ####Set up subgroups
            sampleSubgroups = dict(sampleName=sampleName + ('-' + curSample['DS'] if args.includeSampleIDinSampleCol else ''), assay=curSample['Assay'], DS=curSample['DS'])
            if 'replicate' in subGroupNames:
                sampleSubgroups['replicate'] = curSample['Replicate']
            if 'age' in subGroupNames:
                sampleSubgroups['age'] = re.sub(' ', '_', curSample['Age'])
            if 'mappedgenome' in subGroupNames:
                sampleSubgroups['mappedgenome'] = re.sub(r'^cegsvectors_', '', curSample['Mapped_Genome'])
            
            
            ###shortLabel must be <= 17 printable characters (sampleName + DS)
            sampleShortLabel = sampleName + "-" + curSample['DS']
            
            
            ###
            track = Track(
                name=sampleName_trackname + '_bam',
                short_label=sampleShortLabel,
                long_label=assay_type + ' Reads ' + sampleDescription,
                url=args.URLbase + curSample['filebase'] + '.bam',
                subgroups=sampleSubgroups,
                tracktype='bam',
                parentonoff="off"
            )
            Reads_view.add_tracks(track)
            
            #Coverage_view
            if assay_type in ["DNA", "Capture"]:
                track = Track(
                    name=sampleName_trackname + '_cov',
                    short_label=sampleShortLabel,
                    long_label=assay_type + ' Coverage ' + sampleDescription,
                    url=args.URLbase + curSample['filebase'] + '.coverage.bw',
                    subgroups=sampleSubgroups,
                    tracktype='bigWig',
                    color=curSample['Color'],
                    parentonoff=DensCovTracksDefaultDisplayMode
                )
                Coverage_view.add_tracks(track)
            
            #Dens_view
            if assay_type in ["DNase-seq", "ChIP-seq", "ATAC-seq"]:
                track = Track(
                    name=sampleName_trackname + '_dens',
                    short_label=sampleShortLabel,
                    long_label=assay_type + ' Density ' + sampleDescription,
                    url=args.URLbase + curSample['filebase'] + '.density.bw',
                    subgroups=sampleSubgroups, 
                    tracktype='bigWig',
                    color=curSample['Color'],
                    parentonoff=DensCovTracksDefaultDisplayMode
                )
                Dens_view.add_tracks(track)
            
            #Hotspots_view
            if assay_type in ["DNase-seq", "ChIP-seq", "ATAC-seq"]:
                hotspot_base = os.path.basename(curSample['filebase'])
                hotspot_dir = os.path.dirname(curSample['filebase'])
                hotspot_path = hotspot_dir + "/hotspot2/" + hotspot_base
                
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
                track = Track(
                    name=sampleName_trackname + '_cuts',
                    short_label=sampleShortLabel,
                    long_label=assay_type + ' Cut counts ' + sampleDescription,
                    url=args.URLbase + curSample['filebase'] + '.perBase.bw',
                    subgroups=sampleSubgroups,
                    tracktype='bigWig',
                    color=curSample['Color'],
                    parentonoff="off"
                )
                Cuts_view.add_tracks(track)
            
            #Variants_view
            if assay_type in ["DNA", "Capture"]:
                track = Track(
                    name=sampleName_trackname + '_vcf',
                    short_label=sampleShortLabel,
                    long_label=assay_type + ' Variants ' + sampleDescription,
                    url=args.URLbase + curSample['filebase'] + '.filtered.vcf.gz',
                    subgroups=sampleSubgroups,
                    tracktype='vcfTabix',
                    parentonoff="off"
                )
                Variants_view.add_tracks(track)
            
            #Genotypes_view
            if assay_type in ["DNA", "Capture"]:
                track = Track(
                    name=sampleName_trackname + '_gts',
                    short_label=sampleShortLabel,
                    long_label=assay_type + ' Genotypes ' + sampleDescription,
                    url=args.URLbase + curSample['filebase'] + '.genotypes.bb',
                    subgroups=sampleSubgroups,
                    tracktype='bigBed 5',
                    color=curSample['Color'],
                    parentonoff="off"
                )
                Genotypes_view.add_tracks(track)
        
        
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
