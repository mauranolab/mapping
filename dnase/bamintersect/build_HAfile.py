#!/bin/env python3.5

import sys
import argparse
import pysam
import csv

#########################################################################################################################################

def build_HAfile_f(samfile_HA, samfile_mates, outputFile):

    # Open the input/output files.
    try:
        HA_samfile = pysam.AlignmentFile(samfile_HA, "r")
        Mate_samfile = pysam.AlignmentFile(samfile_mates, "r")

        bed12_output = open(outputFile, 'w', newline='')
        bed12_writer = csv.writer(bed12_output, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    except:
        print("[build_HAfile.py]", "Problem opening input and/or output files.", file=sys.stderr)
        return [101]
    
    
    # Read the samfile lines, then write the bed12 fields.
    while True:
        try:
            HA_read = next(HA_samfile)
            Mate_read = next(Mate_samfile)
        except:
            # Hit end of input files.
            bed12_output.close()
            return [0]
        
        try:
            HA_readID = HA_read.query_name
            Mate_readID = Mate_read.query_name
        
            if HA_readID != Mate_readID:
                # This should not happen. There must have been a problem in setting up the sam files.
                print("[build_HAfile.py]", "The problem that should not happen happened.", HA_readID, Mate_readID, file=sys.stderr)
                return [201]
            
            HA_chrom = HA_samfile.get_reference_name(HA_read.reference_id)
            Mate_chrom = Mate_samfile.get_reference_name(Mate_read.reference_id)
            
            # The pysam reference_start values are 0-based.
            # sam file format values are 1-based.
            HA_read_start = HA_read.reference_start       # This is 0-based !
            Mate_read_start = Mate_read.reference_start   # This is 0-based !
            
            # pysam's reference_end points to one past the last aligned residue. Returns None if
            # not available (read is unmapped or no cigar alignment present).
            #
            HA_read_end = HA_read.reference_end           # This is 0-based !
            Mate_read_end = Mate_read.reference_end       # This is 0-based !
            #
            # The above uses the sam CIGAR string. If start_pos=2 (0-based), and CIGAR=4S10M5S, then:
            #                                      reference_end=12 (0-based)
            #
            # A result of 'None' causes two tabs in a row to appear in the csv file, or just a trailing
            # tab if the field would have been the last in the row.
            
            HA_flag = HA_read.flag
            Mate_flag = Mate_read.flag
            
            if (HA_read.is_reverse):
                HA_strand = "-"
            else:
                HA_strand = "+"
            
            if (Mate_read.is_reverse):
                Mate_strand = "-"
            else:
                Mate_strand = "+"
        
            # Recall that the positions printed below are 0-based.
            bed12_writer.writerow([HA_chrom] + [HA_read_start] + [HA_read_end] +
                                  [HA_readID] + [HA_flag] + [HA_strand] +
                                  [Mate_chrom] + [Mate_read_start] + [Mate_read_end] +
                                  [Mate_readID] + [Mate_flag] + [Mate_strand] )
        except:
            print("[build_HAfile.py]", "Problem in main code body.", file=sys.stderr)
            bed12_output.close()
            return [202]

#########################################################################################################################################

if (__name__ == '__main__'):
    parser = argparse.ArgumentParser(prog = "build_HAfile.py",
                                     description = "Builds a bed12 file for HA/non-HA read pairs. Inputs are sam files containing the read mates.",
                                     add_help = True)

    parser.add_argument('--samfile_HA', action='store', type=str, help='sam file with reads in the HA.')
    parser.add_argument('--samfile_mates', action='store', type=str, help='sam file with reads not in the HA.')
    parser.add_argument('--outputFile', action='store', type=str, help='Name of the bed12 file which will be built.')

    args = parser.parse_args()

    print("[build_HAfile.py] Parameters:", args, file=sys.stderr)

    build_HAfile_out = build_HAfile_f(args.samfile_HA, args.samfile_mates, args.outputFile)
    sys.exit(build_HAfile_out[0])

