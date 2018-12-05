#!/bin/env python3

import json
import collections
import re
import csv
import os
import sys
import argparse


parser = argparse.ArgumentParser(prog = "Decode ENCODE metadata", description = "Decodes ENCODE metadata and finds paired reads, biodata and matching libraries", add_help=True)
parser.add_argument('ENCODEmetadata', action='store', help='Barcode file')
parser.add_argument('-o','--output', action='store', dest='output',help='Decoded metadata')
parser.add_argument('-j','--JSONdirectory', action='store', dest='jsonD', help='JSON directory')

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
    sampleName=sample_line['File accession']
    print(sampleName)
    
    print(sample_line)
    try:
        #Hiercheal json python
        data=open(args.jsonD+sampleName).read()
        sample_json = json.loads(data)
        if sample_json['replicate']['experiment']['biosample_summary'] is not None:
            sampleAge = sample_json['replicate']['experiment']['biosample_summary']
            try: 
                sampleAge = sampleAge[sampleAge.index("(") + 1:sampleAge.rindex(")")]
            except:
                sampleAge = 'NA'
        print(sampleAge)
        #keyword = 'human '
        regexp = re.compile(r'DS[0-9]{5}|DS[0-9]{4}')
        before_keyword, keyword, after_keyword = sample_json['replicate']['experiment']['biosample_summary'].partition('human ')
        sampleCellType = after_keyword.split(' ')[0].replace(' ', '_').replace('-','_').replace('+','').replace(',','').replace('.','').replace('\\)','')
        if not sampleCellType and regexp.search(str(sample_json)) is not None:
            if len(sample_json['replicate']['experiment']['biosample_synonyms'])>=1:
                sampleCellType=sample_json['replicate']['experiment']['biosample_synonyms'][0].replace(' ', '_').replace('-','_').replace('+','').replace(',','').replace('.','').replace('\\)','')
            elif len(sample_json['replicate']['experiment']['biosample_synonyms'])==0:
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
            DSidOut = sampleName, DSid
            print(DSidOut)
        elif regexp.search(str(sample_json['aliases'])) is not None:
            m=re.search(regexp, str(sample_json['aliases']))
            DSid=m.group(0)
            DSidOut = sampleName, DSid 
            print(DSidOut)
        elif 'submitted_file_name' in sample_json and regexp.search(str(sample_json['submitted_file_name'])) is not None:
            m=re.search(regexp, str(sample_json['submitted_file_name']))
            DSid=m.group(0)
            DSidOut = sampleName, DSid
            print(DSidOut)

        elif regexp.search(str(sample_json)) is None:
            DSidOut = sampleName,  sample_json['replicate']['library']['@id'].split('/')[2]
            print(DSidOut)
            
            
        sampleCellType = re.sub("Entire_|entire_|Entire|entire", "", sampleCellType)
        if re.search('embryo', str(sample_json['replicate']['experiment']['biosample_summary'])) is not None:
            sampleCellType = 'f' + sampleCellType.capitalize()
        
        #BUGBUG I think the missing paired_with is obselete
        if sample_line['Run type'] =='paired-ended':
            if sample_line['Paired end'] =='1' and sample_line['Audit ERROR']!='missing paired_with\n':
                DSidOut1 = DSidOut + (sampleCellType, sample_line['Biosample term name'],  'paired-ended', sampleName, sample_line['Paired with'],sample_json['replicate']['experiment']['biosample_summary'],sample_json['lab']['institute_label'],sampleAge,  sample_json['replicate']['library']['accession'])
                
            elif sample_line['Paired end'] =='2' and sample_line['Audit ERROR']!='missing paired_with\n':
                DSidOut1 = DSidOut + (sampleCellType, sample_line['Biosample term name'], 'paired-ended', sample_line['Paired with'], sampleName,sample_json['replicate']['experiment']['biosample_summary'],sample_json['lab']['institute_label'],sampleAge, sample_json['replicate']['library']['accession'])
                
            #If paired_with is missing
            elif sample_line['Paired end'] =='1' and sample_line['Audit ERROR']=='missing paired_with\n':
                DSidOut1 = DSidOut + (sampleCellType, sample_line['Biosample term name'], 'paired-ended', sampleName, '',sample_json['replicate']['experiment']['biosample_summary'],sample_json['lab']['institute_label'],sampleAge, sample_json['replicate']['library']['accession'])
                
        elif sample_line['Run type'] =='single-ended':
            DSidOut1 = DSidOut + (sampleCellType, sample_line['Biosample term name'], 'single-ended', sampleName, '',sample_json['replicate']['experiment']['biosample_summary'],sample_json['lab']['institute_label'],sampleAge, sample_json['replicate']['library']['accession'])
   
  #print(SidOut1)
        #stripped_sample_line = [s.rstrip() for s in DSidOut]
        #print([stripped_sample_line])
        wr.writerows([DSidOut1])
    except Exception: 
        pass
        print('WARNING')
        if sample_line['Run type'] =='paired-ended':
            if sample_line['Paired end'] =='1' and sample_line['Audit ERROR']!='missing paired_with\n':
                DSidOut1 =  DSidOut + (sampleCellType, sample_line['Biosample term name'], 'paired-ended', sampleName, sample_line['Paired with'],sample_json['replicate']['experiment']['biosample_summary'],sample_json['lab']['institute_label'],sampleAge, sample_json['replicate']['library']['accession'])
                
            elif sample_line['Paired end'] =='2' and sample_line['Audit ERROR']!='missing paired_with\n':
                DSidOut1 = DSidOut + (sampleCellType, sample_line['Biosample term name'], 'paired-ended', sample_line['Paired with'], sampleName,sample_json['replicate']['experiment']['biosample_summary'],sample_json['lab']['institute_label'],sampleAge, sample_json['replicate']['library']['accession'])
                
            #If paired_with is missing
            elif sample_line['Paired end'] =='1' and sample_line['Audit ERROR']=='missing paired_with\n':
                DSidOut1 = DSidOut + (sampleCellType, sample_line['Biosample term name'], 'paired-ended', sampleName, '',sample_json['replicate']['experiment']['biosample_summary'],sample_json['lab']['institute_label'],sampleAge, sample_json['replicate']['library']['accession'])
                
        elif sample_line['Run type'] =='single-ended':
            DSidOut1 = DSidOut +  (sampleCellType, sample_line['Biosample term name'], 'single-ended', sampleName, '',sample_json['replicate']['experiment']['biosample_summary'],sample_json['lab']['institute_label'],sampleAge, sample_json['replicate']['library']['accession'])
            
        print(DSidOut1)
        wr.writerows([DSidOut1])

try:
    with open(args.ENCODEmetadata, 'r') as SampleFile:
        SampleFile_reader = csv.DictReader(SampleFile, delimiter='\t')
        
        outfile = open(args.output, 'w')
        # Pass file handle to csv.writer
        wr = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True)
        wr.writerow(['SampleID','GroupID','sampleCellType','Name','Single_or_Paired','R1','R2','Assay','Variable','Age', 'Library'])
        for sample_line in SampleFile_reader:
            getDS(sample_line)
finally:
    outfile.close()
