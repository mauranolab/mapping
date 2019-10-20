#!/usr/bin/env python3
##################################################################################################################
# Usage:
#
#    ./subsetBAM.py [file containing a list of readnames] < [sam file] > subsetBAM.outputFile
#
#    samtools view [BAM file] | subsetBAM.py [file containing a list of readnames]  > subsetBAM.outputFile
##################################################################################################################

import sys

# Open the read name file, strip off any trailing CRLF characters, then put all the read names into a set.
with open(sys.argv[1], 'r') as indexfile:
    readNameSet = set(readName.rstrip('\r\n') for readName in indexfile)

# Read in the samfile, line by line.
for samLine in sys.stdin:
    # split the samfile line into the qname and "all other"
    qname, _ = samLine.split('\t', 1)
    # Is the qname is in our readNameSet? 
    # Regarding speed, see: https://stackoverflow.com/questions/3949310/how-is-set-implemented
    if qname in readNameSet:
        sys.stdout.write(samLine)

