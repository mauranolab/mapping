#!/bin/env python
from trackhub import Track, CompositeTrack, ViewTrack, default_hub
from sys import argv
import sys
import re
import argparse
import io
from collections import OrderedDict
import locale
import csv


parser = argparse.ArgumentParser(prog = "createTrackhub.py", description = "Creates a trackhub with composite and subtracks for all tracks.", add_help=True)
parser.add_argument('Input', action='store', help='Input tsv file with sample metadata. Name, DS, Replicate, Color, Assay, analyzed_reads, SPOT, Num.hotspots, Group, Age, filebase')
#NB Exclude not used right now
parser.add_argument('--genome', action='store', required=True, help='genome assembly name')
parser.add_argument('--assay', action='store', choices=['DNase-seq', 'ChIP-seq'], required=True, help='Name of assay')
parser.add_argument('--URLbase', action='store', required=True, help='URL base at which track files are hosted')


def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument, file=sys.stderr)
    sys.exit(2)


doDNase = args.assay == "DNase-seq"


########
#import metadata
########
if args.Input=="-":
    inputfile = sys.stdin
else:
    inputfile = open(args.Input, 'r') 
#inputfile = open('/vol/isg/encode/dnase/SamplesForTrackhub.tsv', 'r') 

inputfile_reader = csv.DictReader(inputfile, delimiter='\t')
input_data = []
for line in inputfile_reader:
    input_data.append(line)

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


########
#Create composite track for major groups of samples (e.g. Duke, UW, etc.)
########
#Find all unique cellTypes and create a supertrack for each
Group = sorted(set([line['Group'] for line in input_data]))

for curGroup in Group:
    composite = CompositeTrack(
        name=curGroup.replace(" ", "_"),
        short_label=curGroup,
        long_label=curGroup,
        tracktype="bigBed",
        priority="2",
        visibility="hide",
        sortOrder="cellType=+ Assay=+ view=+ replicate=+ DSnumber=+ Age=+",
        maxItems=10000,
        autoScale = "off",
        bigDataUrl ='NULL',
        maxHeightPixels= "6:32:128")
    
    
    ########
    #Create view track
    ########
    #Notes
    #adding viewUi="off" seems to expand the panel by default, opposite expectation
    
    Dens_view=ViewTrack(
        name="Dens_view_" + curGroup.replace(" ", "_"),
        view="Density",
        visibility="full",
        tracktype="bigWig",
        viewLimits = "0:3",
        autoScale='off',
        maxItems=10000,
        short_label="Density",
        long_label="Density")
    composite.add_view(Dens_view)
        
    if doDNase:
        Cuts_view=ViewTrack(
            name="Cuts_view_" + curGroup.replace(" ", "_"),
            view="Cuts",
            visibility="hide",
            tracktype="bigWig",
            viewLimits = "0:2",
            autoScale='off',
            maxItems=10000,
            short_label="Cut Counts",
            long_label="Cut Counts")
        composite.add_view(Cuts_view)
    
    Hotspots_view = ViewTrack(
        name="Hotspots_view_" + curGroup.replace(" ", "_"),
        view="Hotspots",
        visibility="hide",
        tracktype="bigBed 3",
        maxItems=10000,
        short_label="Hotspots",
        long_label="Hotspots")
    composite.add_view(Hotspots_view)
    
    Peaks_view = ViewTrack(
        name="Peaks_view_" + curGroup.replace(" ", "_"),
        view="Peaks",
        visibility="hide",
        tracktype="bigBed 3",
        maxItems=10000,
        short_label="Peaks",
        long_label="Peaks")
    composite.add_view(Peaks_view)
        
    Reads_view = ViewTrack(
        name="Reads_view_" + curGroup.replace(" ", "_"),
        view="Reads",
        visibility="hide",
        tracktype="bam",
        short_label="Reads",
        maxItems=10000,
        pairEndsByName="on",
        long_label="Reads")
    composite.add_view(Reads_view)
    
    for view in composite.views:
        view.add_params(configurable="on")
    
    
    ########
    #Create subgroup definitions
    ########
    matchingSamples = [line for line in input_data if curGroup == line['Group']]
    cellTypes = sorted(set([line['Name'] for line in matchingSamples]))
    cellTypes= [re.sub(r'_L$|_R$', '', i) for i in cellTypes]
    DSnumbers = sorted(set([line['DS'] for line in matchingSamples]))
    Replicates = sorted(set([line['Replicate'] for line in matchingSamples]))
    Assay = sorted(set([line['Assay'] for line in matchingSamples]))
    Age = sorted(set([re.sub(' ', '_', line['Age']) for line in matchingSamples]), key=natural_key)
    
    
    celltypeDict = dict(zip(cellTypes, cellTypes))
    dictAssay = dict(zip(Assay, Assay))
    DSdict = dict(zip(DSnumbers, DSnumbers))
    repDict = dict(zip(Replicates, Replicates))
    ageDict= dict(zip(Age, Age))
    
    from trackhub.track import SubGroupDefinition
    subgroups = [
        SubGroupDefinition(
            name="cellType",
            label="CellType",
            mapping=OrderedDict(sorted(celltypeDict.items(), key=lambda x: x[1]))),
        
        SubGroupDefinition(
            name="Assay",
            label="Assay",
            mapping=OrderedDict(sorted(dictAssay.items(), key=lambda x: x[1]))),
        
        SubGroupDefinition(
            name="replicate",
            label="Replicate",
            mapping=OrderedDict(sorted(repDict.items(), key=lambda x: x[1]))),
            #mapping=dict(rep1="rep1", rep2="rep2", repOther="repOther")),
            
        SubGroupDefinition(
            name="DSnumber",
            label="DSnumber",
            mapping=OrderedDict(sorted(DSdict.items(), key=lambda x: x[1]))),
            
        SubGroupDefinition(
            name="Age",
            label="Age",
            mapping=OrderedDict(sorted(ageDict.items(), key=lambda x: x[1]))),
    ]
    composite.add_subgroups(subgroups)
    
    
    ########
    #Make tracks for all bigWigs in current dir
    ########
    #Create long label from analyzed reads, SPOTS and Hotspots. Format numbers into 1000's
    locale.setlocale(locale.LC_ALL, 'en_US')
    for curSample in matchingSamples:
        cellTypeLR = re.sub(r'_L$|_R$', '', curSample['Name'])
        cellTypeLR_trackname = re.sub(r'\W+', '',cellTypeLR) + "-" + curSample['DS']
        
        #BUGBUG age not geting included for mouse ENCODE
        sampleDescription = cellTypeLR + "-" + curSample['DS'] + ' (' + locale.format("%d", int(curSample['analyzed_reads']), grouping=True) + ' analyzed reads, ' + curSample['SPOT'] + ' SPOT, ' + locale.format("%d", int(curSample['Num_hotspots']), grouping=True) + ' Hotspots)' + (', Age = ' + curSample['Age'] if curGroup == 'Fetal-REMC' else '')
        sampleSubgroups = dict(cellType=cellTypeLR, Assay=curSample['Assay'], replicate=curSample['Replicate'], DSnumber=curSample['DS'], Age=re.sub(' ', '_', curSample['Age']) if curGroup == 'Fetal-REMC' else 'NoAge')
        
        track = Track(
            name=cellTypeLR_trackname + '-reads',
            short_label=cellTypeLR + "-" + curSample['DS'],
            long_label=args.assay + ' Reads ' + sampleDescription,
            autoScale='off',
            url=args.URLbase + curSample['filebase'] + '.bam',
            subgroups=sampleSubgroups, 
            tracktype='bam',
            parentonoff="off",
            )
        Reads_view.add_tracks(track)
        
        #Dens_view
        track = Track(
            name=cellTypeLR_trackname + '-dens',
            short_label=cellTypeLR + "-" + curSample['DS'],
            long_label=args.assay + ' Density ' + sampleDescription,
            viewLimits='0:3',
            autoScale='off',
            url=args.URLbase + curSample['filebase'] + '.bw',
            subgroups=sampleSubgroups, 
            tracktype='bigWig',
            color=curSample['Color'],
            parentonoff="off",
            )
        Dens_view.add_tracks(track)
        
        if doDNase:
            #PerBaseCutCount
            track = Track(
                name=cellTypeLR_trackname + '-cuts',
                short_label=cellTypeLR + "-" + curSample['DS'],
                long_label=args.assay + ' Cut counts ' + sampleDescription,
                autoScale='off',
                url=args.URLbase + curSample['filebase'] + '.perBase.bw',
                subgroups=sampleSubgroups, 
                tracktype='bigWig',
                color=curSample['Color'],
                parentonoff="off",
                )
            Cuts_view.add_tracks(track)
        
        #Peaks_View
        track = Track(
            name=cellTypeLR_trackname + '-peaks',
            short_label=cellTypeLR + "-" + curSample['DS'],
            long_label=args.assay + ' Peaks ' + sampleDescription,
            autoScale='off',
            url=args.URLbase + curSample['filebase'].replace("/", "/hotspot2/") + '.peaks.bb',
            subgroups=sampleSubgroups, 
            tracktype='bigBed 3',
            color=curSample['Color'],
            parentonoff="off",
            )
        Peaks_view.add_tracks(track)
        
        #Hotspots_view
        track = Track(
            name=cellTypeLR_trackname + '-hotspots',
            short_label=cellTypeLR + "-" + curSample['DS'],
            long_label=args.assay + ' Hotspots (5% FDR) ' + sampleDescription,
            autoScale='off',
            url=args.URLbase + curSample['filebase'].replace("/", "/hotspot2/") + '.hotspots.fdr0.05.bb',
            subgroups=sampleSubgroups,
            tracktype='bigBed 3',
            color=curSample['Color'],
            parentonoff="off",
            )
        Hotspots_view.add_tracks(track)
        
    
    #TODO #daler github trackhub/trackhub/track.py doesn't seem to offer dimensionAchecked as option
    if doDNase:
        composite.add_params(dimensions="dimY=cellType dimA=replicate dimX=Age") #MTM 2018jul30 broke: , dimensionAchecked ="rep1"
    else:
        composite.add_params(dimensions="dimY=Assay dimA=replicate dimX=cellType")#, dimensionAchecked ="rep1",  dimensionBchecked="UW"
    
    
    # Demonstrate some post-creation adjustments...here, just make control samples gray
    #for track in trackdb.tracks:
    #    if 'control' in track.name:
    #     track.add_params(color="100,100,100")
    
    
    print(composite)
    print("\n")
    trackdb.add_tracks(composite)

