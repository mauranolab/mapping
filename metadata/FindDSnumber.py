import json,urllib.request,collections,re,csv,os,sys,argparse



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

try:
       SampleName = open(args.ENCODEmetadata, 'r')
       Sample_data = SampleName.readlines()
finally:
       SampleName.close()




hdr = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.11 (KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11',
       'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
       'Accept-Charset': 'ISO-8859-1,utf-8;q=0.7,*;q=0.3',
       'Accept-Encoding': 'none',
       'Accept-Language': 'en-US,en;q=0.8',
       'Connection': 'keep-alive'}

def getDS(inputmetadata):
    for line in Sample_data[1:]:
        SampleName = line.split('\t')[0]
        print(SampleName)
        ###IF JSON DATA IS NOT STORED LOCALLY, UNCOMMENT THE (1,2,3,6) LINES AND COMMENT THE 7TH LINE AFTER THIS 
        #url = "https://www.encodeproject.org/files/" + SampleName + "/?format=json"
        #print(SampleName)
        #req = urllib.request.Request(url, headers=hdr)
        try:
            #Hiercheal json python
            #data = urllib.request.urlopen(req).read()
            data=open(args.jsonD+SampleName).read()
            output = json.loads(data)
            #keyword = 'human '
            regexp = re.compile(r'DS[0-9]{5}|DS[0-9]{4}')
            befor_keyowrd, keyword, after_keyword = output['replicate']['experiment']['description'].partition('human ')
            cellType = after_keyword.split(' ')[0].replace(' ', '_').replace('-','_').replace('+','').replace(',','').replace('.','').replace('\\)','')
            if not cellType and regexp.search(str(output)) is not None:
            	if len(output['replicate']['experiment']['biosample_synonyms'])>=1:
            		cellType=output['replicate']['experiment']['biosample_synonyms'][0].replace(' ', '_').replace('-','_').replace('+','').replace(',','').replace('.','').replace('\\)','')
            	elif len(output['replicate']['experiment']['biosample_synonyms'])==0:
            		cellType = line.split('\t')[6].replace(' ', '_').replace('-','_').replace('+','').replace(',','').replace('.','').replace('\\)','')
            elif not cellType and regexp.search(str(output)) is None:
            	cellType = line.split('\t')[6].replace(' ', '_').replace('-','_').replace('+','').replace(',','').replace('.','').replace('\\)','')
            #output['replicate']['experiment']['description']
            #print(output['replicate']['experiment']['description'])
            # submitted_file_name DS number is here
             #replicate/library/aliases OR in aliases
            if regexp.search(str(output['replicate']['library']['aliases'])) is not None:
                m=re.search(regexp, str(output['replicate']['library']['aliases']))
                DSid = DSid=m.group(0)
                DSidOut = SampleName, DSid,output['replicate']['library']['@id'].split('/')[2]
                print(DSidOut)
            elif regexp.search(str(output['aliases'])) is not None:
                m=re.search(regexp, str(output['aliases']))
                DSid=m.group(0)
                DSidOut = SampleName, DSid ,output['replicate']['library']['@id'].split('/')[2]
                print(DSidOut)
            elif regexp.search(str(output['submitted_file_name'])) is not None:
                m=re.search(regexp, str(output['submitted_file_name']))
                DSid=m.group(0)
                DSidOut = SampleName, DSid,output['replicate']['library']['@id'].split('/')[2]
                print(DSidOut)

            elif regexp.search(str(output)) is None:
                DSidOut = SampleName,  output['replicate']['library']['@id'].split('/')[2]
                print(DSidOut)
                
                
            if line.split('\t')[32] =='paired-ended':
                if line.split('\t')[33] =='1' and line.split('\t')[47]!='missing paired_with\n':
                    DSidOut1 = DSidOut + (cellType, line.split('\t')[6],  'paired-ended', SampleName, line.split('\t')[34],output['replicate']['experiment']['description'],output['lab']['institute_label'])
                    
                elif line.split('\t')[33] =='2' and line.split('\t')[47]!='missing paired_with\n':
                    DSidOut1 = DSidOut + (cellType, line.split('\t')[6], 'paired-ended', line.split('\t')[34], SampleName,output['replicate']['experiment']['description'],output['lab']['institute_label'])
                    
                #If paired_with is missing
                elif line.split('\t')[33] =='1' and line.split('\t')[47]=='missing paired_with\n':
                    DSidOut1 = DSidOut + (cellType, line.split('\t')[6], 'paired-ended', SampleName, '',output['replicate']['experiment']['description'],output['lab']['institute_label'])
                    
            elif line.split('\t')[32] =='single-ended':
                DSidOut1 = DSidOut + (cellType, line.split('\t')[6], 'single-ended', SampleName, '',output['replicate']['experiment']['description'],output['lab']['institute_label'])
        
            print(DSidOut1)
            #stripped_line = [s.rstrip() for s in DSidOut]
            #print([stripped_line])
            wr.writerows([DSidOut1])
        except Exception: 
            pass
            print('wrong')
            if line.split('\t')[32] =='paired-ended':
                if line.split('\t')[33] =='1' and line.split('\t')[47]!='missing paired_with\n':
                    DSidOut1 =  DSidOut + (cellType, line.split('\t')[6], 'paired-ended', SampleName, line.split('\t')[34],output['replicate']['experiment']['description'],output['lab']['institute_label'])
                    
                elif line.split('\t')[33] =='2' and line.split('\t')[47]!='missing paired_with\n':
                    DSidOut1 = DSidOut + (cellType, line.split('\t')[6], 'paired-ended', line.split('\t')[34], SampleName,output['replicate']['experiment']['description'],output['lab']['institute_label'])
                    
                #If paired_with is missing
                elif line.split('\t')[33] =='1' and line.split('\t')[47]=='missing paired_with\n':
                    DSidOut1 = DSidOut + (cellType, line.split('\t')[6], 'paired-ended', SampleName, '',output['replicate']['experiment']['description'],output['lab']['institute_label'])
                    
            elif line.split('\t')[32] =='single-ended':
                DSidOut1 = DSidOut +  (cellType, line.split('\t')[6], 'single-ended', SampleName, '',output['replicate']['experiment']['description'],output['lab']['institute_label'])
                
            print(DSidOut1)
            wr.writerows([DSidOut1])


try:
        outfile = open(args.output, 'w')#use 
        # Pass file handle to csv.writer
        wr = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True)
        wr.writerow(['SampleID','GroupID','cellType','Name','Single_or_Paired','R1','R2','Assay','Variable'])
        getDS(SampleName)
finally:
        outfile.close()
