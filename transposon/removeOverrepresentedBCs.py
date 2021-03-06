#!/bin/env python3
from sys import argv
import sys
import argparse
import csv
import collections
import os
from collections import defaultdict


def process_lines(input_data, wr):
    #Checks which barcodes are the most common before deduplication
    myCounts = collections.Counter([line[col] for line in input_data])
    
    #Don't include failed BCs in further counts
    del myCounts['']
    totalBCs = sum(myCounts.values())
    
    if args.verbose:
        print("\nNumber of unique barcodes:", totalBCs, file=sys.stderr)
        print("The 10 most common BCs before removal:", myCounts.most_common(10), file=sys.stderr)
    
    failedBCs = [x for x in myCounts.keys() if myCounts[x] / totalBCs > freq and x != "NoBC"]
    
    print("[removeOverrepresentedBCs] Number of barcodes removed:", len(failedBCs), file=sys.stderr)
    
    for line in input_data:
        oldBC = line[col]
        if oldBC in failedBCs:
            line[col] = ""
            #check rather than throwing an exception in case a malformed file has non-uniform number of columns
            for maskcol in maskcols:
                if maskcol < len(line):
                    line[maskcol] = ""
        wr.writerows([line])


###Argument parsing
version="1.2"

parser = argparse.ArgumentParser(prog = "removeOverrepresentedBCs", description = "Removes overrepresented BCs based on frequency", allow_abbrev=False)
parser.add_argument('inputfilename', action='store', help='input filename. Format: tab-delimited with barcode sequences (other columns are passed through)')
parser.add_argument("--col", action='store', type=int, default=1, help = "Column with barcodes in it [%(default)s]")
parser.add_argument("--maskcols", action='store', type=str, help = "Comma-separated list of column numbers. For failed BCs, these columns to be masked with an empty screen but not analyzed [%(default)s]")
parser.add_argument("--freq", action='store', type=float, default=0.01, help = "Mask barcodes present at this frequency or higher [%(default)s]")
parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s ' + version)
parser.add_argument('-o','--output', action='store', dest='output',help='Deduplicated barcode file')


try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument, file=sys.stderr)
    sys.exit(2)

col = args.col-1
maskcols = [int(i)-1 for i in args.maskcols.split(",")]
maskcols.sort()
freq = args.freq


###Main
print("[removeOverrepresentedBCs] Removing overrepresented barcodes", file=sys.stderr)
print("[removeOverrepresentedBCs] " + str(args), file=sys.stderr)


## open file handles
if args.inputfilename=="-":
    inputfile = sys.stdin
else:
    inputfile = open(args.inputfilename, 'r') 

input_data = inputfile.readlines()
input_data = [line.rstrip("\n").split('\t') for line in input_data]


if args.output=="-":
    outfile = sys.stdout
else:
    outfile = open(args.output, 'w')
wr = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True)

process_lines(input_data, wr)


print("[removeOverrepresentedBCs] All barcodes have been processed", file=sys.stderr)


inputfile.close()
outfile.close()
