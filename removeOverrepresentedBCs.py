#!/bin/env python3
from sys import argv
import sys
#import re
import argparse
import csv
import collections
import os
from collections import defaultdict


def process_lines(input_data, wr):
    #Checks which barcodes are the most common before deduplication
    bcColNum = col
    myCounts = collections.Counter([line[bcColNum] for line in input_data])
    
    #Don't include failed BCs in further counts
    del myCounts['']
    totalBCs = sum(myCounts.values())
     
    if args.verbose:
        print("\nNumber of unique barcodes:", totalBCs, file=sys.stderr)
        print("The 10 most common BCs before removal:",myCounts.most_common(10), file=sys.stderr)
    
    failedBCs = [x for x in myCounts.keys() if myCounts[x] / totalBCs > freq]
    
    if args.verbose:
        print("\nNumber of barcodes removed:", len(failedBCs), file=sys.stderr)
    
    for line in input_data:
        oldBC = line[bcColNum]
        if oldBC in failedBCs:
            line[bcColNum] = ""
        wr.writerows([line])


###Main
parser = argparse.ArgumentParser(prog = "removeOverrepresentedBCs", description = "Removes overrepresented BCs based on frequency", allow_abbrev=False)
parser.add_argument('inputfilename', action='store', help='input filename. Format: tab-delimited with barcode sequences (other columns are passed through)')
parser.add_argument("--col", action='store', type=int, default=1, help = "Column with barcodes in it [%(default)s]")
parser.add_argument("--freq", action='store', type=float, default=0.01, help = "Mask barcodes present at this frequency or higher [%(default)s]")
parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-o','--output', action='store', dest='output',help='Deduplicated barcode file')


try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument, file=sys.stderr)
    sys.exit(2)

col = args.col-1
freq = args.freq


print("Removing overrepresented barcodes", file=sys.stderr)
print(args, file=sys.stderr)


## open file handles
if args.inputfilename=="-":
    inputfile = sys.stdin
else:
    inputfile = open(args.inputfilename, 'r') 

input_data = inputfile.readlines()
input_data = [line.rstrip().split('\t') for line in input_data]


if args.output=="-":
    outfile = sys.stdout
else:
    outfile = open(args.output, 'w')
wr = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True)

process_lines(input_data, wr)


print("\nAll barcodes have been processed", file=sys.stderr)


inputfile.close()
outfile.close()
