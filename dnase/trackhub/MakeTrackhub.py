#!/bin/env python
from trackhub import Track, CompositeTrack, ViewTrack, default_hub, SuperTrack
from sys import argv
import sys
import re
import argparse
import io
from collections import OrderedDict
import locale
import csv
import os

parser = argparse.ArgumentParser(prog = "MakeTrackhub.py", description = "Creates a trackhub with composite and subtracks for all tracks.", add_help = True)
parser.add_argument('Input', action = 'store', help = 'Input tsv file with sample metadata. Name, DS, Replicate, Color, Assay, analyzed_reads, SPOT, Num_hotspots, Exclude, Genomic_coverage, Group, Age, filebase')
#NB Exclude not used right now 
parser.add_argument('--genome', action = 'store', required = True, help = 'genome assembly name')
parser.add_argument('--URLbase', action = 'store', required = True, help = 'URL base at which track files are hosted')
parser.add_argument('--includeSampleIDinSampleCol', action = 'store_true', default = False, help = 'Append the Sample ID in the Sample column')
parser.add_argument('--supertrack', action = 'store', required = False, help = 'Encompass all composite tracks generated within supertrack. Supertrack name specifiied as parameter.')
parser.add_argument('--generateHTMLdescription', action = 'store_true', default = False, help = 'Link to HTML descriptions for composite tracks, assumed to be present at [genome]/descriptions/[group name].html')

def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

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


assays = set([])
for assay_type in set([line.get('Assay') for line in input_data_all]):
    if(assay_type == "DNase-seq" or assay_type == "DNA"):
        assays.add(assay_type)
    else:
        #ChIP-seq tracks have the epitope as the assay type
        #NB this will pick up any samples whose type doesn't match the above
        assays.add("ChIP-seq")


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


# group="cegsvectors" keeps Flowcells, Aggregations, and cegsvectors out of the "Other" control group,
# and in the CEGS control group, when the selected genome is cegsvectors. When the selected geneome is
# one of the standard ones (hg38, mm10, rn6), then it is ignored. In that scenario, all our controls 
# are placed in a group titled with the short hub label.
if(args.supertrack is not None):
    if(args.genome == "cegsvectors"):
        supertrack = SuperTrack(
            name=args.supertrack,
            group="cegsvectors",
            short_label=args.supertrack,
            long_label=args.supertrack)
        trackdb.add_tracks(supertrack)
    else:
        supertrack = SuperTrack(
            name=args.supertrack,
            short_label=args.supertrack,
            long_label=args.supertrack)
        trackdb.add_tracks(supertrack)


# Now build the composites, and add them all to the supertrack.
for assay_type in assays:
    assay_suffix = assay_type.split("-")[0]
    
    input_data = []
    for line in input_data_all:
        if(line['Assay'] == assay_type or assay_type == "ChIP-seq" and line['Assay'] not in assays):
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
        # Suppress Replicate and/or Age in the Track Settings screen if they are all NAs.
        # SortOrder controls what is displayed, even if we retain all the subgroup definitions.
        SortOrder = "sampleName=+ DS=+ view=+ "
        
        uniqueReplicates = [line['Replicate'] for line in input_data if(line['Group'] == curGroup)]
        if(not all(x == 'NA' for x in uniqueReplicates)):
            SortOrder = SortOrder + "replicate=+ "
        
        uniqueAges = [line['Age'] for line in input_data if(line['Group'] == curGroup)]
        if(not all(x == 'NA' for x in uniqueAges)):
            SortOrder = SortOrder + "Age=+ "
        
        # Adding suffix
        curGroup_trackname = cleanTrackName(args.genome + "_" + curGroup + "_" + assay_suffix)
        
        # Do something like this later.
        # mydate = re.split(r'_', curGroup)[0]
        # priority=bignumber - mydate
        # Why is priority="2" not showing up in the output?
        # Should be a number anyway.
        
        # Change bigBed to 'bed 3' ?
        #    Typically, a single level composite has the same type as all of its subtracks 
        #    and offers user configuration options at the top level. But a multi-view composite
        #    is most often given the bare-bones "type bed 3", and offers user configuration
        #    options at the view level. Exceptions to this pattern do exist but they are rare.
        
        # short labels are supposed to be 17 chars max. However, the browser does not seem to care.
        shrt_label = curGroup + '_' + assay_type
        # shrt_label = shrt_label[0:17]
        # long labels are supposed to be 76 chars max. However, the browser does not seem to care.
        lng_label = curGroup + '_' +  assay_type
        # lng_label = lng_label[0:76]
        
        if(args.generateHTMLdescription):
            html_file = 'descriptions/' + curGroup + '.html'
        else:
            html_file = ' '
        
        composite = CompositeTrack(
            name=curGroup_trackname,
            short_label=shrt_label,
            long_label=lng_label,
            tracktype="bigBed",
            priority="2",
            visibility="hide",
            sortOrder=SortOrder,
            autoScale = "off",
            html=html_file,
            bigDataUrl ='NULL',
            maxHeightPixels= "100:32:8")
        
        ########
        #Create view track
        ########
        #Notes
        #adding viewUi="off" seems to expand the panel by default, opposite expectation
        
        if (assay_type == "DNase-seq") or (assay_type == "ChIP-seq"):
            Dens_view = ViewTrack(
                name="Dens_view_" + curGroup_trackname,
                view="Density",
                visibility="full",
                tracktype="bigWig",
                viewLimits = "0:3",
                autoScale='off',
                short_label="Density",
                long_label="Density")
            composite.add_view(Dens_view)
            
        if (assay_type == "DNase-seq"):
            Cuts_view = ViewTrack(
                name="Cuts_view_" + curGroup_trackname,
                view="Cuts",
                visibility="hide",
                tracktype="bigWig",
                viewLimits = "0:2",
                autoScale='off',
                short_label="Cut Counts",
                long_label="Cut Counts")
            composite.add_view(Cuts_view)
        
        if (assay_type == "DNase-seq") or (assay_type == "ChIP-seq"):
            Hotspots_view = ViewTrack(
                name="Hotspots_view_" + curGroup_trackname,
                view="Hotspots",
                visibility="hide",
                tracktype="bigBed 3",
                maxItems=10000,
                short_label="Hotspots",
                long_label="Hotspots")
            composite.add_view(Hotspots_view)
        
        if (assay_type == "DNase-seq") or (assay_type == "ChIP-seq"):
            Peaks_view = ViewTrack(
                name="Peaks_view_" + curGroup_trackname,
                view="Peaks",
                visibility="hide",
                tracktype="bigBed 3",
                maxItems=10000,
                short_label="Peaks",
                long_label="Peaks")
            composite.add_view(Peaks_view)
            
        if (assay_type == "DNA"):
            Coverage_view = ViewTrack(
                name="Coverage_view_" + curGroup_trackname,
                view="Coverage",
                visibility="full",
                tracktype="bigWig",
                viewLimits = "0:500",
                autoScale='on',
                short_label="Coverage",
                long_label="Coverage")
            composite.add_view(Coverage_view)
        
        if (assay_type == "DNA"):
            Variants_view = ViewTrack(
                name="Variants_view_" + curGroup_trackname,
                view="Variants",
                visibility="hide",
                tracktype="bigBed",
                autoScale='off',
                maxItems=10000,
                short_label="Variants",
                long_label="Variants")
            composite.add_view(Variants_view)
        
        Reads_view = ViewTrack(
            name="Reads_view_" + curGroup_trackname,
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
        sampleNames = sorted(set([re.sub(r'_L$|_R$', '', line['Name']) + ('-' + line['DS'] if args.includeSampleIDinSampleCol else '') for line in matchingSamples]))
        DSnumbers = sorted(set([line['DS'] for line in matchingSamples]))
        Replicates = sorted(set([line['Replicate'] for line in matchingSamples]))
        Assay = sorted(set([line['Assay'] for line in matchingSamples]))
        Age = sorted(set([re.sub(' ', '_', line['Age']) for line in matchingSamples]), key=natural_key)
        
        sampleNameDict = dict(zip(sampleNames, sampleNames))
        dictAssay = dict(zip(Assay, Assay))
        DSdict = dict(zip(DSnumbers, DSnumbers))
        repDict = dict(zip(Replicates, Replicates))
        ageDict= dict(zip(Age, Age))
        
        from trackhub.track import SubGroupDefinition
        subgroups = [
            SubGroupDefinition(
                name="sampleName",
                label="Sample",
                mapping=OrderedDict(sorted(sampleNameDict.items(), key=lambda x: x[1]))),
            
            SubGroupDefinition(
                name="Assay",
                label="Assay",
                mapping=OrderedDict(sorted(dictAssay.items(), key=lambda x: x[1]))),
            
            SubGroupDefinition(
                name="replicate",
                label="Replicate",
                mapping=OrderedDict(sorted(repDict.items(), key=lambda x: x[1]))),
                
            # The "label" below  needs to have the underscore in it, else only the part up to the first
            # space will be displayed. However, the underscore will not be displayed in the browser. 
            # Kent converts the underscores to spaces for display purposes.
            # Hyphens seem to be retained though.
            SubGroupDefinition(
                name="DS",
                label="Sample_ID",
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
            sampleName = re.sub(r'_L$|_R$', '', curSample['Name'])
            
            #shortLabel must be <= 17 printable characters (sampleName + DS)
            sampleShortLabel = sampleName + "-" + curSample['DS']
            
            #Internal track name (never displayed) must be unique within the Genome Browser or dataHub and must begin with a letter
            #There is a length limit of 128 characters, but in practice some are used for an internal prefix/suffix
            #Probably could omit sampleName here to save space
            sampleName_trackname = cleanTrackName(args.genome + "_" + curGroup + "_" + sampleName + "_" + curSample['DS'])
            
            #longLabel must be <= 76 printable characters
            if assay_type == "DNA":
                sampleDescription = sampleShortLabel + ' ' + curSample['Genomic_coverage'] + 'x Genomic Coverage (' + locale.format("%d", int(curSample['analyzed_reads']), grouping=True) + ' analyzed reads)' + (', Age = ' + curSample['Age'] if curSample['Age'] != 'NA' else '')
            else:
                if curSample['Num_hotspots'] == "NA":
                    sampleDescriptionNumHotspots = "no"
                else:
                    sampleDescriptionNumHotspots = locale.format("%d", int(curSample['Num_hotspots']), grouping=True) 
                sampleDescription = sampleShortLabel + ' (' + locale.format("%d", int(curSample['analyzed_reads']), grouping=True) + ' analyzed reads, ' + curSample['SPOT'] + ' SPOT, ' + sampleDescriptionNumHotspots + ' Hotspots)' + (', Age = ' + curSample['Age'] if curSample['Age'] != 'NA' else '')
            
            sampleSubgroups = dict(sampleName=sampleName + ('-' + curSample['DS'] if args.includeSampleIDinSampleCol else ''), Assay=curSample['Assay'], replicate=curSample['Replicate'], DS=curSample['DS'], Age=re.sub(' ', '_', curSample['Age']))
            
            
            track = Track(
                name=sampleName_trackname + '_bam',
                short_label=sampleShortLabel,
                long_label=assay_type + ' Reads ' + sampleDescription,
                url=args.URLbase + curSample['filebase'] + '.bam',
                subgroups=sampleSubgroups,
                tracktype='bam',
                parentonoff="off",
                )
            Reads_view.add_tracks(track)
            
            #Variants_view
            if (assay_type == "DNA"):
                track = Track(
                    name=sampleName_trackname + '_var',
                    short_label=sampleShortLabel,
                    long_label=assay_type + ' Variants ' + sampleDescription,
                    url=args.URLbase + curSample['filebase'] + '.variants.bb',
                    subgroups=sampleSubgroups,
                    tracktype='bigBed',
                    color=curSample['Color'],
                    parentonoff="off",
                    )
                Variants_view.add_tracks(track)
            
            #Coverage_view
            if (assay_type == "DNA"):
                track = Track(
                    name=sampleName_trackname + '_cov',
                    short_label=sampleShortLabel,
                    long_label=assay_type + ' Coverage ' + sampleDescription,
                    url=args.URLbase + curSample['filebase'] + '.coverage.bw',
                    subgroups=sampleSubgroups,
                    tracktype='bigWig',
                    color=curSample['Color'],
                    parentonoff="off",
                    )
                Coverage_view.add_tracks(track)
            
            #Dens_view
            if (assay_type == "DNase-seq") or (assay_type == "ChIP-seq"):
                track = Track(
                    name=sampleName_trackname + '_dens',
                    short_label=sampleShortLabel,
                    long_label=assay_type + ' Density ' + sampleDescription,
                    url=args.URLbase + curSample['filebase'] + '.density.bw',
                    subgroups=sampleSubgroups, 
                    tracktype='bigWig',
                    color=curSample['Color'],
                    parentonoff="off",
                    )
                Dens_view.add_tracks(track)
            
            #PerBaseCutCount
            if (assay_type == "DNase-seq"):
                track = Track(
                    name=sampleName_trackname + '_cuts',
                    short_label=sampleShortLabel,
                    long_label=assay_type + ' Cut counts ' + sampleDescription,
                    url=args.URLbase + curSample['filebase'] + '.perBase.bw',
                    subgroups=sampleSubgroups,
                    tracktype='bigWig',
                    color=curSample['Color'],
                    parentonoff="off",
                    )
                Cuts_view.add_tracks(track)
            
            #Peaks_View
            if (assay_type == "DNase-seq") or (assay_type == "ChIP-seq"):
                hotspot_base = os.path.basename(curSample['filebase'])
                hotspot_dir = os.path.dirname(curSample['filebase'])
                hotspot_path = hotspot_dir + "/hotspot2/" + hotspot_base
                
                track = Track(
                    name=sampleName_trackname + '_pks',
                    short_label=sampleShortLabel,
                    long_label=assay_type + ' Peaks ' + sampleDescription,
                    url=args.URLbase + hotspot_path  + '.peaks.bb',
                    subgroups=sampleSubgroups,
                    tracktype='bigBed 3',
                    color=curSample['Color'],
                    parentonoff="off",
                    )
                Peaks_view.add_tracks(track)
            
            #Hotspots_view
            if (assay_type == "DNase-seq") or (assay_type == "ChIP-seq"):
                hotspot_base = os.path.basename(curSample['filebase'])
                hotspot_dir = os.path.dirname(curSample['filebase'])
                hotspot_path = hotspot_dir + "/hotspot2/" + hotspot_base
                
                track = Track(
                    name=sampleName_trackname + '_hot',
                    short_label=sampleShortLabel,
                    long_label=assay_type + ' Hotspots (5% FDR) ' + sampleDescription,
                    url=args.URLbase + hotspot_path + '.hotspots.fdr0.05.bb',
                    subgroups=sampleSubgroups,
                    tracktype='bigBed 3',
                    color=curSample['Color'],
                    parentonoff="off",
                    )
                Hotspots_view.add_tracks(track)
        
        
        #TODO needs to pick the dimensions better
        #TODO #daler github trackhub/trackhub/track.py doesn't seem to offer dimensionAchecked as option
        if (assay_type == "DNase-seq"):
            composite.add_params(dimensions="dimY=sampleName dimA=replicate dimX=Age")
        elif (assay_type == "DNA"):
            composite.add_params(dimensions="dimY=sampleName")
        else:
            # ChIP-seq
            composite.add_params(dimensions="dimY=Assay dimA=replicate dimX=sampleName")
        
        
        # Demonstrate some post-creation adjustments...here, just make control samples gray
        #for track in trackdb.tracks:
        #    if 'control' in track.name:
        #     track.add_params(color="100,100,100")
        
        if(args.supertrack is not None):
            supertrack.add_tracks(composite)
        else:
            trackdb.add_tracks(composite)

#This would print the entire set of tracks if we took out the #'s.
#print("\n")
print(trackdb)
print("\n")
