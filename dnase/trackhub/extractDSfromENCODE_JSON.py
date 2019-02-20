#!/bin/env python3

import json
import collections
import re
import csv
import os
import sys
import argparse


parser = argparse.ArgumentParser(prog = "Decode ENCODE metadata", description = "Decodes ENCODE metadata and finds paired reads, biodata and matching libraries", add_help=True)
parser.add_argument('ENCODEmetadata', action='store', help='metadata.tsv from encodeproject.org')
parser.add_argument('-o', '--output', action='store', dest='output', help='Decoded metadata')
parser.add_argument('-j', '--JSONdirectory', action='store', dest='jsonD', help='JSON directory')

try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument, file=sys.stderr)
    sys.exit(2)

print("\nFile:", file=sys.stderr)
print(args, file=sys.stderr)


if args.ENCODEmetadata=="-":
    ENCODEmetadata = sys.stdin


def getDS(sample_line):
    sampleAccession = sample_line['File accession']
    #print(sample_line)
    
    data = open(args.jsonD+sampleAccession).read()
    sample_json = json.loads(data)
    
    regexp = re.compile(r'DS[0-9]{5}|DS[0-9]{4}')
    #I don't think this works for the mouse data (e.g. "C57BL/6 liver male adult (8 weeks)")
    before_keyword, keyword, after_keyword = sample_json['replicate']['experiment']['biosample_summary'].partition('human ')
    sampleCellType = after_keyword.split(' ')[0].replace(' ', '_').replace('-','_').replace('+','').replace(',','').replace('.','').replace('\\)','')
    if not sampleCellType and regexp.search(str(sample_json)) is not None:
        if 'biosample_synonyms' in sample_json['replicate']['experiment']:
            sampleCellType = sample_json['replicate']['experiment']['biosample_synonyms'][0].replace(' ', '_').replace('-','_').replace('+','').replace(',','').replace('.','').replace('\\)','')
        else:
            sampleCellType = sample_line['Biosample term name'].replace(' ', '_').replace('-','_').replace('+','').replace(',','').replace('.','').replace('\\)','')
    elif not sampleCellType and regexp.search(str(sample_json)) is None:
        sampleCellType = sample_line['Biosample term name'].replace(' ', '_').replace('-','_').replace('+','').replace(',','').replace('.','').replace('\\)','')
    #sample_json['replicate']['experiment']['description']
    #print(sample_json['replicate']['experiment']['description'])
    # submitted_file_name DS number is here
    #replicate/library/aliases OR in aliases
    if regexp.search(str(sample_json['replicate']['library']['aliases'])) is not None:
        m=re.search(regexp, str(sample_json['replicate']['library']['aliases']))
        DSid = DSid=m.group(0)
    elif regexp.search(str(sample_json['aliases'])) is not None:
        m=re.search(regexp, str(sample_json['aliases']))
        DSid = m.group(0)
    elif 'submitted_file_name' in sample_json and regexp.search(str(sample_json['submitted_file_name'])) is not None:
        m=re.search(regexp, str(sample_json['submitted_file_name']))
        DSid = m.group(0)
    elif regexp.search(str(sample_json)) is None:
        DSid = sample_json['replicate']['library']['@id'].split('/')[2]
    sampleCellType = re.sub("Entire_|entire_|Entire|entire", "", sampleCellType)
    if re.search('embryo', str(sample_json['replicate']['experiment']['biosample_summary'])) is not None:
        sampleCellType = 'f' + sampleCellType.capitalize()
    
    #BUGBUG I think the missing paired_with is obselete
    if sample_line['Run type'] == 'paired-ended':
        singleOrPaired = 'paired-ended'
        if sample_line['Paired end'] == '1' and sample_line['Audit ERROR'] != 'missing paired_with\n':
            R1sampleAccession = sampleAccession
            R2sampleAccession = sample_line['Paired with']
        elif sample_line['Paired end'] == '2' and sample_line['Audit ERROR'] != 'missing paired_with\n':
            R1sampleAccession = sample_line['Paired with']
            R2sampleAccession = sampleAccession
        #If paired_with is missing
        elif sample_line['Paired end'] == '1' and sample_line['Audit ERROR'] == 'missing paired_with\n':
            R1sampleAccession = sampleAccession
            R2sampleAccession = ''
    elif sample_line['Run type'] =='single-ended':
        singleOrPaired = 'single-ended'
        R1sampleAccession = sampleAccession
        R2sampleAccession = ''
    
    if sample_json['replicate']['experiment']['assay_title'] == "ChIP-seq":
        sampleAssay = sample_json['replicate']['experiment']['target']['label']
    else:
        sampleAssay = sample_json['replicate']['experiment']['assay_title']
    
    sampleInstitution = sample_json['lab']['institute_label']
    
    if sample_json['replicate']['experiment']['biosample_summary'] is not None:
        sampleAge = sample_json['replicate']['experiment']['biosample_summary']
        try: 
            sampleAge = sampleAge[sampleAge.index("(") + 1:sampleAge.rindex(")")]
        except:
            sampleAge = 'NA'
    
    libraryAccession = sample_json['replicate']['library']['accession']
    
    print(sampleAccession, DSid, sampleAge)
    wr.writerow([sampleAccession, libraryAccession, DSid, sampleCellType, sampleCellType, singleOrPaired, R1sampleAccession, R2sampleAccession, sampleAssay, sampleInstitution, sampleAge])

try:
    with open(args.ENCODEmetadata, 'r') as SampleFile:
        SampleFile_reader = csv.DictReader(SampleFile, delimiter='\t')
        
        outfile = open(args.output, 'w')
        # Pass file handle to csv.writer
        wr = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True)
        #N.B. Library provides ENCODE accession for samples with a DS number. paired-ended appears to be for legibility alone
        #Provide "Automatic Name" as information so that it's easy to recover names which have been manually changed from the automatic ones
        wr.writerow(['SampleAccession', 'LibraryAccession', 'DS', 'Automatic Name', 'Name', 'Single_or_Paired', 'R1', 'R2', 'Assay', 'Institution', 'Age'])
        for sample_line in SampleFile_reader:
            getDS(sample_line)
finally:
    outfile.close()
