#!/usr/bin/env python3

# subsetBAM.py - Select reads from a BAM file. Will filter all alignments for a given readname.

#Usage:
#subsetBAM.py in.bam out.bam

import sys
import argparse
import pysam

version="1.1"

parser = argparse.ArgumentParser(prog = "subsetBAM.py", description = "Returns a BAM file subsetted by a list of read names.", allow_abbrev=False, add_help = True)
parser.add_argument('raw_alignment', type=str, help='Input alignment filename (BAM/SAM)')
parser.add_argument('filtered_alignment', type=str, help='Output alignment filename (uncompressed BAM)')
parser.add_argument('--exclude_flags', type=int, default=0, help='SAM flags to filter out reads')
parser.add_argument('--include_readnames', type=str, required=False, help='Name of a file containing readnames to be retained')
parser.add_argument('--exclude_readnames', type=str, required=False, help='Name of a file containing readnames to be excluded')
parser.add_argument('--version', action='version', version='%(prog)s ' + version)

#argparse does not set an exit code upon this error
#https://stackoverflow.com/questions/5943249/python-argparse-and-controlling-overriding-the-exit-status-code
try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument)
    sys.exit(2)

print("[subsetBAM.py] Parameters:", args, file=sys.stderr)


# Open the read name file, strip off any trailing CRLF characters, then put all the read names into a set.
if args.include_readnames is not None:
    with open(args.include_readnames, 'r') as readNameFile:
        include_readnames = set(readName.rstrip('\r\n') for readName in readNameFile)

if args.exclude_readnames is not None:
    with open(args.exclude_readnames, 'r') as readNameFile:
        exclude_readnames = set(readName.rstrip('\r\n') for readName in readNameFile)

# Open the existing bam/sam file.
bamFile_in = pysam.AlignmentFile(args.raw_alignment, "r")


# Get the old header and make a new header and PG tag from it.
newheader = bamFile_in.header.to_dict()
pg={'ID':'subsetBAM.py', 'PP':newheader['PG'][-1]['ID'], 'PN':'subsetBAM.py', 'VN':version, 'CL':' '.join(sys.argv)}
newheader['PG'].append(pg)

# Open the new uncompressed bam file with the new header.
bamFile_out = pysam.AlignmentFile(args.filtered_alignment, "wbu", header=newheader)

# Read the bamfile line by line.
alignmentsRead = 0
alignmentsSkipped = 0
while True:
    try:
        bamFileRead = next(bamFile_in)
    except:
        # We ran out of lines, so exit.
         break
        
    # Exclude some reads.
    alignmentsRead += 1
    if bamFileRead.flag & args.exclude_flags > 0:
        alignmentsSkipped += 1
        continue
    
    readname = bamFileRead.query_name
    
    if args.exclude_readnames is not None and readname in exclude_readnames:
        alignmentsSkipped += 1
        continue
        
    
    if args.include_readnames is None or readname in include_readnames:
        bamFile_out.write(bamFileRead)
    else:
        alignmentsSkipped += 1

print("[subsetBAM.py] Finished reading ", alignmentsRead, " alignments, of which ", alignmentsSkipped, " were skipped.", file=sys.stderr, sep="")
