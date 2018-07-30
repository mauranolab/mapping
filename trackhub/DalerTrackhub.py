#!/bin/env python

#TODO convert to python3
#TODO missing newline between root track definitions

#http://genome.isg.med.nyu.edu/cgi-bin/hgTrackUi?hgsid=68_v7cxkrY8qgHBWt4NPi9xStLOu4XD&c=chr12&g=hub_48_composite_mouseencode_chipseq
#Try to follow outline of /vol/isg/encode/mouseencode/mapped/mm10/trackDb.txt
from trackhub import Track, default_hub
#from trackhub.upload import upload_hub, upload_track
from sys import argv
import sys
import re
import argparse
import io
from collections import OrderedDict
import locale

parser = argparse.ArgumentParser(prog = "createTrackhub", description = "Creates a trackhub with Composite and subtracks for all tracks. Columns should be: cellType, DSnumber, Replicate, Colour, Assay, nonredundant_Tags, SPOT, Hotspots", add_help=True)
parser.add_argument('Input', action='store', help='Input tsv file with all sample information')
parser.add_argument('-o','--output', action='store', dest='output',default='./',help='Output file for where to store trackhub.')


def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]
#parser.add_argument("--maxPolyG", action='store', type=int, default=75, help = "The minimum percentage of Gs to trim from the sequence [%(default)s].")

try:
       args = parser.parse_args()
except argparse.ArgumentError as exc:
       print(exc.message, '\n', exc.argument)
       sys.exit(2)

########
##open file handles
########
if args.Input=="-":
          inputfile = sys.stdin
else:
          inputfile = open(args.Input, 'r') 


print(args.output)
          
########
#import files
########
#inputfile = open('/vol/isg/encode/dnase201805/SamplesForTrackhub.tsv', 'r') 
#inputfile = open('/home/maagj01/projects/EncodeDNase/testOUT', 'r') 
input_data = inputfile.readlines()
input_data = [line.rstrip().split('\t') for line in input_data][1:]#Don't use header for now TODO use header for tracks

########
#Find all unique cellTypes and create a supertrack for each
########
Variable = sorted(set([i[9] for i in input_data]))


########
#Create hub file
########
hub, genomes_file, genome, trackdb = default_hub(
    hub_name="Encode_DNase",
    genome="hg38",
    short_label="Encode_DNase",
    long_label="Encode DNase",
    email="jesper.maag@med.nyu.edu")


########
#Create composite track for major groups of samples (e.g. Duke, UW, etc.)
########
for trackComp in Variable:
    from trackhub import CompositeTrack
    composite = CompositeTrack(
        name=trackComp.replace(" ", "_"),
        short_label=trackComp.replace(" ", "_"),
        long_label=trackComp.replace(" ", "_"),
        tracktype="bigBed",
        priority="2",
        visibility="hide",
        sortOrder="cellType=+ view=+ replicate=+ DSnumber=+ Age=+",
        maxItems=100000,
        autoScale = "off",
        bigDataUrl ='NULL',
        maxHeightPixels= "6:32:128")
    
    #composite.add_params()
    
    ########
    #Create view track
    ########
    from trackhub import ViewTrack
    NormalizedDensity_view=ViewTrack(
        name="NormalizedDensity_view_"+trackComp.replace(" ", "_"),
        view="NormalizedDensity",
        visibility="full",
        tracktype="bigWig",
        viewLimits = "0:3",
        autoScale='off',
        maxItems=int(100000),
        viewUi="off",
        short_label="NormalizedDensity",
        long_label="NormalizedDensity")
        
    perBaseCleavage_view=ViewTrack(
        name="perBaseCleavage_view_"+trackComp.replace(" ", "_"),
        view="perBaseCleavage",
        visibility="hide",
        tracktype="bigWig",
       viewLimits = "0:2",
        autoScale='off',
        maxItems=int(100000),
        viewUi="off",
        short_label="perBaseCleavage",
        long_label="perBaseCleavage")
    
    Hotspots_view = ViewTrack(
        name="Hotspots_view_"+trackComp.replace(" ", "_"),
        view="Hotspots",
        visibility="hide",
        tracktype="bigBed 3",
        maxItems=100000,
        viewUi="off",
        short_label="Hotspots",
        long_label="Hotspots")
    
    Peaks_view = ViewTrack(
        name="Peaks_view_"+trackComp.replace(" ", "_"),
        view="Peaks",
        visibility="hide",
        tracktype="bigBed 3",
        maxItems=100000,
        viewUi="off",
        short_label="Peaks",
        long_label="Peaks")
        
    Tags_view = ViewTrack(
        name="Tags_view_"+trackComp.replace(" ", "_"),
        view="Tags",
        visibility="hide",
        tracktype="bam",
        short_label="Tags",
        maxItems=100000,
        viewUi="off",
        #viewUI = "off",
        pairEndsByName="on",
        long_label="Tags")
    
    composite.add_view(perBaseCleavage_view)
    composite.add_view(NormalizedDensity_view)
    composite.add_view(Hotspots_view)
    composite.add_view(Peaks_view)
    composite.add_view(Tags_view)
    for view in composite.views:
        view.add_params(configurable="on")
    ########
    #Create subgroup definitions
    ########
    matching = [s for s in input_data if trackComp in s]
    cellTypes = sorted(set([i[0] for i in matching]))
    cellTypes= [re.sub(r'_L$|_R$','',i) for i in cellTypes]
    DSnumbers = sorted(set([i[1] for i in matching]))
    Replicates = sorted(set([i[2] for i in matching]))
    Age = sorted(set([re.sub(' ','_',i[10]) for i in matching]),key=natural_key)
    
    celltypeDict = dict(zip(cellTypes,cellTypes))
    
    DSdict = dict(zip(DSnumbers,DSnumbers))
    
    repDict = dict(zip(Replicates,Replicates))
    ageDict= dict(zip(Age,Age))
    from trackhub.track import SubGroupDefinition
    subgroups = [
    
        SubGroupDefinition(
            name="cellType",
            label="CellType",
            mapping=OrderedDict(sorted(celltypeDict.items(), key=lambda x: x[1]))),
    
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
    #Create long label from nonredundant tags, SPOTS and Hotspots.Format numbers into 1000's
    locale.setlocale(locale.LC_ALL, 'en_US')
    #input_data
    #f#n=input_data[0]
    for fn in matching:
        cellTypeLR = re.sub(r'_L$|_R$','',fn[0])
        track = Track(
            name=re.sub(r'\W+', '',cellTypeLR + '_' + fn[1] + '_tags'),
            short_label=cellTypeLR+'_'+fn[1],
            long_label="DNase Tags "+cellTypeLR+'_'+fn[1]+' ('+locale.format("%d", int(fn[5]), grouping=True)+' uniquely mapped tags, '+fn[6]+' SPOT, '+locale.format("%d", int(fn[7]), grouping=True)+' Hotspots)'+ ', Age = ' + fn[10],
            autoScale='off',
            url='https://cascade.isg.med.nyu.edu/mauranolab/encode/dnase201805/mapped/'+fn[0]+'-'+fn[1]+'.hg38/'+fn[0]+'-'+fn[1]+'.hg38.bam',
            subgroups=dict(cellType=cellTypeLR, replicate=fn[2], DSnumber=fn[1], Age=re.sub(' ','_',fn[10]) if trackComp == 'Fetal_Roadmap' else 'NoAge'),
            tracktype='bam',
            
            )
        Tags_view.add_tracks(track)
        #NormalizedDensity_view
        track = Track(
            name=re.sub(r'\W+', '',cellTypeLR + '_' + fn[1] + '_normalizedDensity'),
            short_label=cellTypeLR+'_'+fn[1],
            long_label='DNase Density ' + cellTypeLR+'_'+fn[1]+' ('+locale.format("%d", int(fn[5]), grouping=True)+' uniquely mapped tags, '+fn[6]+' SPOT, '+locale.format("%d", int(fn[7]), grouping=True)+' Hotspots)'+', Age = ' + fn[10],
            autoScale='off',
            url='https://cascade.isg.med.nyu.edu/mauranolab/encode/dnase201805/mapped/'+fn[0]+'-'+fn[1]+'.hg38/'+fn[0]+'-'+fn[1]+'.hg38.bw',
            subgroups=dict(cellType=cellTypeLR, replicate=fn[2], DSnumber=fn[1], Age=re.sub(' ','_',fn[10]) if trackComp == 'Fetal_Roadmap' else 'NoAge'),
            tracktype='bigWig',
            color=fn[3],
            )
        NormalizedDensity_view.add_tracks(track)
        #PerBaseCutCunt
        track = Track(
            name=re.sub(r'\W+', '',cellTypeLR + '_' + fn[1] + '_perBaseCleavage'),
            short_label=cellTypeLR+'_'+fn[1],
            long_label='DNase Cut counts ' + cellTypeLR + '_'+fn[1]+' ('+locale.format("%d", int(fn[5]), grouping=True)+' uniquely mapped tags, '+fn[6]+' SPOT, '+locale.format("%d", int(fn[7]), grouping=True)+' Hotspots)'+ ', Age = ' + fn[10],
            autoScale='off',
            url='https://cascade.isg.med.nyu.edu/mauranolab/encode/dnase201805/mapped/'+fn[0]+'-'+fn[1]+'.hg38/'+fn[0]+'-'+fn[1]+'.hg38.perBase.bw',
            subgroups=dict(cellType=cellTypeLR, replicate=fn[2], DSnumber=fn[1], Age=re.sub(' ','_',fn[10]) if trackComp == 'Fetal_Roadmap' else 'NoAge'),
            tracktype='bigWig',
            color=fn[3],
            )
        perBaseCleavage_view.add_tracks(track)
        #Peaks_View
        track = Track(
            name=re.sub(r'\W+', '',cellTypeLR + '_' + fn[1] + '_peaks'),
            short_label=cellTypeLR+'_'+fn[1],
            long_label='DNase Peaks ' + cellTypeLR+'_'+fn[1]+' ('+locale.format("%d", int(fn[5]), grouping=True)+' uniquely mapped tags, '+fn[6]+' SPOT, '+locale.format("%d", int(fn[7]), grouping=True)+' Hotspots)'+', Age = ' + fn[10],
            autoScale='off',
            url='https://cascade.isg.med.nyu.edu/mauranolab/encode/dnase201805/mapped/'+fn[0]+'-'+fn[1]+'.hg38'+'/hotspot2/'+fn[0]+'-'+fn[1]+'.hg38.peaks.bb',
            subgroups=dict(cellType=cellTypeLR, replicate=fn[2], DSnumber=fn[1], Age=re.sub(' ','_',fn[10]) if trackComp == 'Fetal_Roadmap' else 'NoAge'),
            tracktype='bigBed 3',
            color=fn[3],
            )
        Peaks_view.add_tracks(track)
        #Hotspots_view
        track = Track(
            name=re.sub(r'\W+', '',cellTypeLR + '_' + fn[1] + '_hotspots'),
            short_label=cellTypeLR + '_'+fn[1],
            long_label='DNase Hotspots ' + cellTypeLR +'_'+fn[1]+' ('+locale.format("%d", int(fn[5]), grouping=True)+' uniquely mapped tags, '+fn[6]+' SPOT, '+locale.format("%d", int(fn[7]), grouping=True)+' Hotspots)'+', Age = ' + fn[10],
            autoScale='off',
            url='https://cascade.isg.med.nyu.edu/mauranolab/encode/dnase201805/mapped/'+fn[0]+'-'+fn[1]+'.hg38'+'/hotspot2/'+fn[0]+'-'+fn[1]+'.hg38.hotspots.fdr0.01.bb',
            subgroups=dict(cellType=cellTypeLR, replicate=fn[2], DSnumber=fn[1], Age=re.sub(' ','_',fn[10]) if trackComp == 'Fetal_Roadmap' else 'NoAge'),
            tracktype='bigBed 3',
            color=fn[3],
            )
        Hotspots_view.add_tracks(track)
        
    for track in NormalizedDensity_view.subtracks:
        track.add_params(viewLimits='0:3', autoScale='off')
        
    composite.add_params( dimensions="dimY=cellType dimA=replicate dimX=Age") #MTM 2018jul30 broke: , dimensionAchecked ="rep1"
    print composite 
    trackdb.add_tracks(composite)

hub.local_fn = args.output
hub.render()

    #supertrack.add_track(track)
# Demonstrate some post-creation adjustments...here, just make control
# samples gray
for track in trackdb.tracks:
    if 'control' in track.name:
        track.add_params(color="100,100,100")


