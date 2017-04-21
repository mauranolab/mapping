#http://genome.isg.med.nyu.edu/cgi-bin/hgTrackUi?hgsid=68_v7cxkrY8qgHBWt4NPi9xStLOu4XD&c=chr12&g=hub_48_composite_mouseencode_chipseq
#Try to follow outline of /vol/isg/encode/mouseencode/mapped/mm10/trackDb.txt
from trackhub import Track, default_hub
from trackhub.upload import upload_hub, upload_track
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
#inputfile = open('../SamplesForTrackhub.tsv', 'r') 
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
#Create composite track
########
for trackComp in Variable:
    from trackhub import CompositeTrack
    composite = CompositeTrack(
        name=trackComp.replace(" ", "_"),
        short_label=trackComp.replace(" ", "_"),
        long_label=trackComp.replace(" ", "_"),
        tracktype="bigBed",
        priority="2",
        filterComposite="on",
        dimensions="dimY=cellType dimA=replicate dimB=DSnumber",
        sortOrder="cellType=+ view=+ replicate=+ DSnumber=+",
        maxItems=100000,
        autoScale = "off",
        bigDataUrl ='NULL',
        maxHeightPixels= "6:32:128")
    
    composite.add_params(dragAndDrop='subtracks', visibility='full')
    
    ########
    #Create view track
    ########
    from trackhub import ViewTrack
    NormalizedDensity_view=ViewTrack(
        name="NormalizedDensity_view_"+trackComp.replace(" ", "_"),
        view="NormalizedDensity",
        visibility="full",
        tracktype="bigWig",
       # viewLimits = "0:5",
        autoScale='off',
        maxItems=int(100000),
        viewUi="off",
        short_label="NormalizedDensity",
        long_label="NormalizedDensity")
    
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
    DSnumbers = sorted(set([i[1] for i in matching]))
    
    dict1=dict(zip(cellTypes,cellTypes))
    
    DSdict=dict(zip(DSnumbers,DSnumbers))
    
    from trackhub.track import SubGroupDefinition
    subgroups = [
    
        SubGroupDefinition(
            name="cellType",
            label="CellType",
            mapping=OrderedDict(sorted(dict1.items(), key=lambda x: x[1]))),
    
        SubGroupDefinition(
            name="replicate",
            label="Replicate",
            mapping=dict(rep1="rep1", rep2="rep2", Other="Other")),
            
        SubGroupDefinition(
            name="DSnumber",
            label="DSnumber",
            mapping=OrderedDict(sorted(DSdict.items(), key=lambda x: x[1]))),
    
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
        track = Track(
            name=re.sub(r'\W+', '',fn[0]+'_'+fn[1]+'_tags'),
            short_label=fn[0]+'_'+fn[1],
            long_label=fn[0]+'_'+fn[1]+' ('+locale.format("%d", int(fn[5]), grouping=True)+' uniquly mapped tags, '+fn[6]+' SPOT, '+locale.format("%d", int(fn[7]), grouping=True)+' Hotspots)',
            autoScale='off',
            url='https://cascade.isg.med.nyu.edu/mauranolab/encode/dnase/mapped/'+fn[0]+'-'+fn[1]+'.hg38/'+fn[0]+'-'+fn[1]+'.hg38.bam',
            subgroups=dict(cellType=fn[0], replicate=fn[2], DSnumber=fn[1]),
            tracktype='bam',
            
            )
        Tags_view.add_tracks(track)
        #NormalizedDensity_view
        track = Track(
            name=re.sub(r'\W+', '',fn[0]+'_'+fn[1]+'_normalizedDensity'),
            short_label=fn[0]+'_'+fn[1],
            long_label=fn[0]+'_'+fn[1]+' ('+locale.format("%d", int(fn[5]), grouping=True)+' uniquly mapped tags, '+fn[6]+' SPOT, '+locale.format("%d", int(fn[7]), grouping=True)+' Hotspots)',
            autoScale='off',
            url='https://cascade.isg.med.nyu.edu/mauranolab/encode/dnase/mapped/'+fn[0]+'-'+fn[1]+'.hg38/'+fn[0]+'-'+fn[1]+'.hg38.bw',
            subgroups=dict(cellType=fn[0], replicate=fn[2], DSnumber=fn[1]),
            tracktype='bigWig',
            color=fn[3],
            )
        NormalizedDensity_view.add_tracks(track)
        #Peaks_View
        track = Track(
            name=re.sub(r'\W+', '',fn[0]+'_'+fn[1]+'_peaks'),
            short_label=fn[0]+'_'+fn[1],
            long_label=fn[0]+'_'+fn[1]+' ('+locale.format("%d", int(fn[5]), grouping=True)+' uniquly mapped tags, '+fn[6]+' SPOT, '+locale.format("%d", int(fn[7]), grouping=True)+' Hotspots)',
            autoScale='off',
            url='https://cascade.isg.med.nyu.edu/mauranolab/encode/dnase/mapped/'+fn[0]+'-'+fn[1]+'.hg38'+'/hotspots/'+fn[0]+'-'+fn[1]+'.hg38.fdr0.01.pks.bb ',
            subgroups=dict(cellType=fn[0], replicate=fn[2], DSnumber=fn[1]),
            tracktype='bigBed 3',
            color=fn[3],
            )
        Peaks_view.add_tracks(track)
        #Hotspots_view
        track = Track(
            name=re.sub(r'\W+', '',fn[0]+'_'+fn[1]+'_hotspots'),
            short_label=fn[0]+'_'+fn[1],
            long_label=fn[0]+'_'+fn[1]+' ('+locale.format("%d", int(fn[5]), grouping=True)+' uniquly mapped tags, '+fn[6]+' SPOT, '+locale.format("%d", int(fn[7]), grouping=True)+' Hotspots)',
            autoScale='off',
            url='https://cascade.isg.med.nyu.edu/mauranolab/encode/dnase/mapped/'+fn[0]+'-'+fn[1]+'.hg38'+'/hotspots/'+fn[0]+'-'+fn[1]+'.hg38.fdr0.01.hot.bb',
            subgroups=dict(cellType=fn[0], replicate=fn[2], DSnumber=fn[1]),
            tracktype='bigBed 3',
            color=fn[3],
            )
        Hotspots_view.add_tracks(track)
        
    for track in NormalizedDensity_view.subtracks:
        track.add_params(viewLimits='0:5', autoScale='off')
    
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


