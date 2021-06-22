#!/bin/env python

import sys
import pysam
import argparse

version="1.0"

def changeFlags(bam_input, bam_output):
    bamfile_in = pysam.AlignmentFile(bam_input, "rb")
    bamfile_out = pysam.AlignmentFile(bam_output, "wb", template=bamfile_in)
    
    for read in bamfile_in:
        if read.is_qcfail:
            read.flag = read.flag - args.clearflag   # Turn off qcfail flag.
        
        bamfile_out.write(read)
    
    bamfile_in.close()
    bamfile_out.close()


if (__name__ == '__main__'):
    parser = argparse.ArgumentParser(prog = "changeSAMflags.py", description = "Turns off specified flags for all reads in a bam file.", add_help = True)
    parser.add_argument('bam_input', type=str, help='A full path bam file name')
    parser.add_argument('bam_output', type=str, help='A full path bam file name')
    parser.add_argument("--clearflag", action = "store", type = int, default = 0, help = "Flag to clear in decimal [%(default)s]")
    parser.add_argument('--version', action='version', version='%(prog)s ' + version)
    
    args = parser.parse_args()
    print("[changeFlags.py] Parameters:", args, file=sys.stderr)
    
    changeFlags(args.bam_input, args.bam_output)

