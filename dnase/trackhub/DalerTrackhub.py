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

#TODO parameterize or infer
#mappedgenome = "mm10"
#doDNase = False
#assay="ChIP-seq"
#URLbase = 'https://cascade.isg.med.nyu.edu/mauranolab/encode/mouseencode_chipseq/mapped/'

mappedgenome = "hg38"
doDNase = True
assay="DNase-seq"
URLbase = 'https://cascade.isg.med.nyu.edu/mauranolab/encode/dnase/mapped/'


parser = argparse.ArgumentParser(prog = "createTrackhub", description = "Creates a trackhub with Composite and subtracks for all tracks.", add_help=True)
parser.add_argument('Input', action='store', help='Input tsv file with all sample information for grouping. Columns should be: cellType, DSnumber, Replicate, Colour, Assay, analyzed_reads, SPOT, Hotspots (TODO update)')
parser.add_argument('-o','--output', action='store', dest='output',default='./',help='Output file for where to store trackhub.')


def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

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
#inputfile = open('/vol/isg/encode/dnase/SamplesForTrackhub.tsv', 'r') 
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
    hub_name="Encode_DNase_201805",
    genome=mappedgenome,
    short_label="Encode_DNase",
    long_label="Encode DNase",
    email="maurano@nyu.edu")


########
#Create composite track for major groups of samples (e.g. Duke, UW, etc.)
########
for trackComp in Variable:
    from trackhub import CompositeTrack
    composite = CompositeTrack(
        name=trackComp.replace(" ", "_"),
        short_label=trackComp,
        long_label=trackComp,
        tracktype="bigBed",
        priority="2",
        visibility="hide",
        sortOrder="cellType=+ Assay=+ view=+ replicate=+ DSnumber=+ Age=+",
        maxItems=10000,
        autoScale = "off",
        bigDataUrl ='NULL',
        maxHeightPixels= "6:32:128")
    
    #composite.add_params()
    
    ########
    #Create view track
    ########
    #Notes
    #adding viewUi="off" seems to expand the panel by default, opposite expectations
    
    from trackhub import ViewTrack
    Dens_view=ViewTrack(
        name="Dens_view_" + trackComp.replace(" ", "_"),
        view="Dens",
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
            name="Cuts_view_" + trackComp.replace(" ", "_"),
            view="Cut counts",
            visibility="hide",
            tracktype="bigWig",
            viewLimits = "0:2",
            autoScale='off',
            maxItems=10000,
            short_label="Cut Counts",
            long_label="Cut Counts")
        composite.add_view(Cuts_view)
    
    Hotspots_view = ViewTrack(
        name="Hotspots_view_" + trackComp.replace(" ", "_"),
        view="Hotspots",
        visibility="hide",
        tracktype="bigBed 3",
        maxItems=10000,
        short_label="Hotspots",
        long_label="Hotspots")
    composite.add_view(Hotspots_view)
    
    Peaks_view = ViewTrack(
        name="Peaks_view_" + trackComp.replace(" ", "_"),
        view="Peaks",
        visibility="hide",
        tracktype="bigBed 3",
        maxItems=10000,
        short_label="Peaks",
        long_label="Peaks")
    composite.add_view(Peaks_view)
        
    Reads_view = ViewTrack(
        name="Reads_view_" + trackComp.replace(" ", "_"),
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
    matching = [s for s in input_data if trackComp in s]
    cellTypes = sorted(set([i[0] for i in matching]))
    cellTypes= [re.sub(r'_L$|_R$','',i) for i in cellTypes]
    DSnumbers = sorted(set([i[1] for i in matching]))
    Replicates = sorted(set([i[2] for i in matching]))
    Assay = sorted(set([i[4] for i in matching]))
    Age = sorted(set([re.sub(' ','_',i[10]) for i in matching]),key=natural_key)
#    Institute = sorted(set([re.sub(' ','_',i[11]) for i in matching]),key=natural_key)
    
    
    celltypeDict = dict(zip(cellTypes,cellTypes))
    dictAssay = dict(zip(Assay,Assay))
    DSdict = dict(zip(DSnumbers,DSnumbers))
    repDict = dict(zip(Replicates,Replicates))
    ageDict= dict(zip(Age,Age))
#    Instdict= dict(zip(Institute,Institute))
    
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
         
#    SubGroupDefinition(
#        name="Institute",
#        label="Institute",
#        mapping=OrderedDict(sorted(Instdict.items(), key=lambda x: x[1]))),
    
    ]
    composite.add_subgroups(subgroups)
    
    
    
    ########
    #Make tracks for all bigWigs in current dir
    ########
    #Create long label from analyzed reads, SPOTS and Hotspots.Format numbers into 1000's
    locale.setlocale(locale.LC_ALL, 'en_US')
    #input_data
    #f#n=input_data[0]
    for fn in matching:
     cellTypeLR = re.sub(r'_L$|_R$', '', fn[0])
     cellTypeLR_trackname = re.sub(r'\W+', '',cellTypeLR) + "-" + fn[1]
     track = Track(
         name=cellTypeLR_trackname + '-reads',
         short_label=cellTypeLR + "-" + fn[1],
         long_label=assay + ' Reads ' + cellTypeLR + "-" + fn[1] + ' (' + locale.format("%d", int(fn[5]), grouping=True)+' analyzed reads, ' + fn[6] + ' SPOT, ' + locale.format("%d", int(fn[7]), grouping=True)+' Hotspots)'+ ', Age = ' + fn[10],
         autoScale='off',
         url=URLbase+ fn[11] + '.bam',
         subgroups=dict(cellType=cellTypeLR, Assay=fn[4], replicate=fn[2], DSnumber=fn[1], Age=re.sub(' ','_',fn[10]) if trackComp == 'Fetal-REMC' else 'NoAge'),
         tracktype='bam',
         parentonoff="off",
         )
     Reads_view.add_tracks(track)
     
     #Dens_view
     track = Track(
         name=cellTypeLR_trackname + '-dens',
         short_label=cellTypeLR + "-" + fn[1],
         long_label=assay + ' Density ' + cellTypeLR + "-" + fn[1] + ' (' + locale.format("%d", int(fn[5]), grouping=True)+' analyzed reads, ' + fn[6] + ' SPOT, ' + locale.format("%d", int(fn[7]), grouping=True)+' Hotspots)'+', Age = ' + fn[10],
         viewLimits='0:3',
         autoScale='off',
         url=URLbase+ fn[11] + '.bw',
         subgroups=dict(cellType=cellTypeLR, Assay=fn[4], replicate=fn[2], DSnumber=fn[1], Age=re.sub(' ','_',fn[10]) if trackComp == 'Fetal-REMC' else 'NoAge'),
         tracktype='bigWig',
         color=fn[3],
        parentonoff="off",
         )
     Dens_view.add_tracks(track)
     
     if doDNase:
         #PerBaseCutCount
         track = Track(
             name=cellTypeLR_trackname + '-cuts',
             short_label=cellTypeLR + "-" + fn[1],
             long_label=assay + ' Cut counts ' + cellTypeLR + "-" + fn[1] + ' (' + locale.format("%d", int(fn[5]), grouping=True)+' analyzed reads, ' + fn[6] + ' SPOT, ' + locale.format("%d", int(fn[7]), grouping=True)+' Hotspots)'+ ', Age = ' + fn[10],
             autoScale='off',
             url=URLbase+ fn[11] + '.perBase.bw',
             subgroups=dict(cellType=cellTypeLR, Assay=fn[4], replicate=fn[2], DSnumber=fn[1], Age=re.sub(' ','_',fn[10]) if trackComp == 'Fetal-REMC' else 'NoAge'),
             tracktype='bigWig',
             color=fn[3],
            parentonoff="off",
             )
         Cuts_view.add_tracks(track)
     
     #Peaks_View
     track = Track(
         name=cellTypeLR_trackname + '-peaks',
         short_label=cellTypeLR + "-" + fn[1],
         long_label=assay + ' Peaks ' + cellTypeLR + "-" + fn[1] + ' (' + locale.format("%d", int(fn[5]), grouping=True)+' analyzed reads, ' + fn[6] + ' SPOT, ' + locale.format("%d", int(fn[7]), grouping=True)+' Hotspots)'+', Age = ' + fn[10],
         autoScale='off',
         url=URLbase + fn[11].replace("/", "/hotspot2/") + '.peaks.bb',
         subgroups=dict(cellType=cellTypeLR, Assay=fn[4], replicate=fn[2], DSnumber=fn[1], Age=re.sub(' ','_',fn[10]) if trackComp == 'Fetal-REMC' else 'NoAge'),
         tracktype='bigBed 3',
         color=fn[3],
        parentonoff="off",
         )
     Peaks_view.add_tracks(track)
     
     #Hotspots_view
     track = Track(
         name=cellTypeLR_trackname + '-hotspots',
         short_label=cellTypeLR + "-" + fn[1],
         long_label=assay + ' Hotspots (5% FDR) ' + cellTypeLR + "-" + fn[1] + ' (' + locale.format("%d", int(fn[5]), grouping=True)+' analyzed reads, ' + fn[6] + ' SPOT, ' + locale.format("%d", int(fn[7]), grouping=True)+' Hotspots)'+', Age = ' + fn[10],
         autoScale='off',
         url=URLbase + fn[11].replace("/", "/hotspot2/") + '.hotspots.fdr0.05.bb',
         subgroups=dict(cellType=cellTypeLR, Assay=fn[4], replicate=fn[2], DSnumber=fn[1], Age=re.sub(' ','_',fn[10]) if trackComp == 'Fetal-REMC' else 'NoAge'),
         tracktype='bigBed 3',
         color=fn[3],
        parentonoff="off",
         )
     Hotspots_view.add_tracks(track)
     
    
    #TODO #daler github trackhub/trackhub/track.py doesn't seem to offer dimensionAchecked as option
    if doDNase:
        composite.add_params(dimensions="dimY=cellType dimA=replicate dimX=Age") #MTM 2018jul30 broke: , dimensionAchecked ="rep1"
    else:
        composite.add_params(dimensions="dimY=Assay dimA=replicate dimX=cellType")#, dimensionAchecked ="rep1",  dimensionBchecked="UW"
    
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


