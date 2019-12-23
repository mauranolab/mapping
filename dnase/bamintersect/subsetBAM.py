#!/usr/bin/env python3
##################################################################################################################
# Usage:
#
# subsetBAM.py --exclude_flags <int> --readNames <text file containing a list of readnames> --bamFile_in <inputFile.bam> --bamFile_out <outputFile.bam>
##################################################################################################################

import sys
import argparse
import pysam

parser = argparse.ArgumentParser(prog = "subsetBAM.py",
                                 description = "Returns a BAM file subsetted by a list of read names.",
                                 add_help = True)

parser.add_argument('--exclude_flags', action='store', type=int, help='SAM flags to filter out certain reads.')
parser.add_argument('--readNames', action='store', type=str, help='Name of a file containing readnames to be extracted from bamFile_in.')
parser.add_argument('--bamFile_in', action='store', type=str, help='Name of the bam file from which bam lines will be retrieved.')
parser.add_argument('--bamFile_out', action='store', type=str, help='Name of the bam file to which subsetted bam lines will be written.')

args = parser.parse_args()

print("[subsetBAM.py] Parameters:", args, file=sys.stderr)


# Open the read name file, strip off any trailing CRLF characters, then put all the read names into a set.
with open(args.readNames, 'r') as readNameFile:
    readNameSet = set(readName.rstrip('\r\n') for readName in readNameFile)

# Open the existing bam file.
bamFile_in = pysam.AlignmentFile(args.bamFile_in, "rb")

# Get the old header and make a new header and PG tag from it.
header = pysam.AlignmentFile(args.bamFile_in,'rb').header
headerDict = header.to_dict()
version="1.0"
pg={'ID':'subsetBAM.py', 'PP':header['PG'][-1]['ID'], 'PN':'subsetBAM.py', 'VN':version, 'CL':' '.join(sys.argv)}
headerDict['PG'].append(pg)

# Open the new bam file with the new header.
bamFile_out = pysam.AlignmentFile(args.bamFile_out, "wb", header=headerDict)

# Read the bamfile, line by line.
while True:
    try:
        bamFileRead = next(bamFile_in)

        # Exclude some reads.
        if (bamFileRead.flag & args.exclude_flags) > 0:
            continue

        bamFileQname = bamFileRead.query_name
    except:
        # We ran out of lines, so exit.
         sys.exit(0)

    # Is the qname is in our readNameSet? 
    # Regarding speed, see: https://stackoverflow.com/questions/3949310/how-is-set-implemented
    if bamFileQname in readNameSet:
        bamFile_out.write(bamFileRead)

