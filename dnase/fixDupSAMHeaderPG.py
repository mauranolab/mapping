#!/bin/env python3
import pysam
import argparse
import sys
import random


version="1.0"


parser = argparse.ArgumentParser(prog = "fixDupSAMHeaderPG", description = "Corrects duplicate PG entries in a SAM header", allow_abbrev=False)
parser.add_argument("raw_alignment", type = str, help = "Input raw alignment file (SAM)")
parser.add_argument("filtered_alignment", type = str, help = "Output filtered alignment file (uncompressed bam)")
#parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s ' + version)

#argparse does not set an exit code upon this error
#https://stackoverflow.com/questions/5943249/python-argparse-and-controlling-overriding-the-exit-status-code
try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument)
    sys.exit(2)


#print("[fixDupSAMHeaderPG.py] Parameters:", args, file=sys.stderr)


unfiltered_reads = pysam.AlignmentFile(args.raw_alignment, "r")
newheader = unfiltered_reads.header.to_dict()

#uniqueProgNames = set([ prog['ID'] for prog in newheader['PG'] ])

#import pprint
#pprint.pprint(unfiltered_reads.header.to_dict())

visitedProgIDs = []
for prog in newheader['PG']:
    if prog['ID'] in visitedProgIDs:
        while True:
            newID = prog['ID'] + "-dup" + str(random.getrandbits(16))
            if newID not in visitedProgIDs:
                break
        print("[fixDupSAMHeaderPG] Adjusted duplicate @PG header ID:", prog['ID'], " to ", newID, sep="", file=sys.stderr)
        prog['ID'] = newID
        #NB makes no attempt to fix PP tags referring to this tag
    visitedProgIDs.append(prog['ID'])


#Set PP tag to last @PG entry listed
#Reuses visitedProgIDs to find a unique ID for fixDupSAMHeaderPG record
newID = 'fixDupSAMHeaderPG'
while True:
    if newID not in visitedProgIDs:
        break
    newID = 'fixDupSAMHeaderPG'+ "-dup" + str(random.getrandbits(16))

#BUGBUG makes no attempt to follow PG chain, just assumes it's the last one
newheader['PG'].append({'ID':newID, 'PP':newheader['PG'][-1]['ID'], 'PN':'fixDupSAMHeaderPG.py', 'VN':version, 'CL':' '.join(sys.argv)})


try:
    filtered_reads = pysam.AlignmentFile(args.filtered_alignment, "wbu", header = newheader)
    
    allreads = unfiltered_reads.fetch(until_eof=True) #I believe .next would skip unmapped reads and isn't guaranteed to match file
    while(1):
        try:
            read = next(unfiltered_reads)
            ret = filtered_reads.write(read)
        except StopIteration:
            break
finally:
    unfiltered_reads.close()
    filtered_reads.close()

