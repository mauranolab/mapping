#!/bin/env python3.5

import sys
import os
import argparse
import pysam
import csv
import re
from ctypes import CDLL


version="1.0"


def reads_match(fileA_read, fileB_read, same):
    if same:
        return fileA_read.is_read1 == fileB_read.is_read1
    else:
        return fileA_read.is_read1 != fileB_read.is_read1


def validateBothReads(max_mismatches, ReqFullyAligned, file1_read, file2_read):
    return validateSingleRead(max_mismatches, ReqFullyAligned, file1_read) and validateSingleRead(max_mismatches, ReqFullyAligned, file2_read)

def validateSingleRead(max_mismatches, ReqFullyAligned, read):
    # Make sure it is exact matches to the reference.
    if read.get_tag('NM') > max_mismatches:
        return False
    
    if ReqFullyAligned:
        # Before writing the read info to output, make sure they are both fully aligned.
        if re.search('[HSPDI]', read.cigarstring) is not None:
            return False
    
    return True


def write_output(bed_writer, file1_file, file2_file, file1_read, file2_read, bedout, bam1out_file, bam2out_file, max_mismatches, ReqFullyAligned):
    if not validateBothReads(max_mismatches, ReqFullyAligned, file1_read, file2_read):
        return
    
    if bam1out_file is not None:
        bam1out_file.write(file1_read)
    
    if bam2out_file is not None:
        bam2out_file.write(file2_read)
        return
    
    if bed_writer is not None:
        file1_readID = file1_read.query_name
        file2_readID = file2_read.query_name
        
        chrom_1 = file1_file.get_reference_name(file1_read.reference_id)
        chrom_2 = file2_file.get_reference_name(file2_read.reference_id)
        
        # The pysam reference_start values are 0-based.
        # sam file format values are 1-based.
        read_start_1 = file1_read.reference_start   # This is 0-based !
        read_start_2 = file2_read.reference_start   # This is 0-based !
        
        # pysam's reference_end points to one past the last aligned residue. Returns None if
        # not available (read is unmapped or no cigar alignment present).
        #
        read_end_1 = file1_read.reference_end      # This is 0-based !
        read_end_2 = file2_read.reference_end      # This is 0-based !
        #
        # The above uses the sam CIGAR string. If start_pos=2 (0-based), and CIGAR=4S10M5S, then:
        #                                      reference_end=12 (0-based)
        #
        # A result of 'None' causes two tabs in a row to appear in the csv file, or just a trailing
        # tab if the field would have been the last in the row.
        
        read_flag_1 = file1_read.flag
        read_flag_2 = file2_read.flag
        
        if file1_read.is_reverse:
            strand_1 = "-"
        else:
            strand_1 = "+"
        
        if file2_read.is_reverse:
            strand_2 = "-"
        else:
            strand_2 = "+"
        
        # Recall that the positions printed below are 0-based.
        bed_writer.writerow([chrom_1] + [read_start_1] + [read_end_1] + [file1_readID] + [read_flag_1] + [strand_1] + [chrom_2] + [read_start_2] + [read_end_2] + [file2_readID] + [read_flag_2] + [strand_2] )


def samtoolsCmpReadnames(cCode, read1, operator, read2):
    read1_b = read1.encode('utf-8')
    read2_b = read2.encode('utf-8')

    retValue = cCode.strnum_cmp(read1_b, read2_b)

    if operator == ">":
        return retValue > 0
    elif operator == "<":
        return retValue < 0
    elif operator == "==":
        return retValue == 0
    else:
        raise NameError('samtoolsCmpReadnames was called with an invalid operator.')


def openBamOutFile(bamin, bamout):
    if bamout is not None:
        newheader = bamin.header.to_dict()
        #Set PP tag to last @PG entry listed, also makes no attempt to follow PG chain in case it is run on a merged bam file.
        newheader['PG'].append({'ID':'bamintersect.py', 'PP':newheader['PG'][-1]['ID'], 'PN':'bamintersect.py', 'VN':version, 'CL':' '.join(sys.argv)})
        bamout_file = pysam.AlignmentFile(bamout, "wb", header=newheader)
    else:
        bamout_file = None

    return bamout_file

def bam_intersect_f(bam_name1, bam_name2, same, bedout, bam1out, bam2out, max_mismatches, ReqFullyAligned):
    # Get iterator handles for bam input files #1 and #2.
    file1_file = pysam.AlignmentFile(bam_name1, "rb")
    file2_file = pysam.AlignmentFile(bam_name2, "rb")
    
    if bedout is not None:
        bed_file = open(bedout, 'w', newline='')
        bed_writer = csv.writer(bed_file, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True, quotechar='"', quoting=csv.QUOTE_MINIMAL)
    else:
        bed_file = None
        bed_writer = None
    
    bam1out_file = openBamOutFile(file1_file, bam1out)
    bam2out_file = openBamOutFile(file2_file, bam2out)
    
    ##########################################################################
    # Inital read pulls. 'next' is an iterator method.
    try:
        # It has happened that the "LP backbone" has no mapped reads, for which returning an error here causes 
        file1_read = next(file1_file)  # Reads lines in file1_file in same order as they are in the file.
        file1_readID = file1_read.query_name
        
        # Returns an exception when the number of mapped reads for the chromosome is zero.
        # This is a common occurance for the non-standard chromosomes in all the genomes, and
        # can happen at random for the standard chromosomes as well.
        file2_read = next(file2_file)
        file2_readID = file2_read.query_name
    except:
        return
    
    # The logic in this loop handles scenarios in which file1 contains both a read1 and a read2 from a single readID.
    # So it handles both this scenario:
    #                               file1                     file2
    #                           readID_1<read1>           readID_1<read2>              <== a match
    #                           readID_1<read2>           readID_1<read1>              <== a match
    # And this scenario:
    #                               file1                     file2
    #                           readID_1<read1>           readID_1<read1>              <=== not a match
    #                           readID_1<read2>           readID_1<read2>              <== The file2 read matches
    #                                                                                      the previous file1 read, 
    #                                                                                      and the file1 read matches
    #                                                                                      the previous file2 read.
    #
    # It also handles scenarios in which one file contains both read1 and read2,
    # while the other file contains only a single read.
    
    try:
        while True:
            if samtoolsCmpReadnames(cCode, file1_readID, "<", file2_readID):
                # The file1 readID is less than the file2 readID. Get another file1 line, and check again.
                try: file1_read = next(file1_file)
                except: break  # All done
                file1_readID = file1_read.query_name
            elif file1_readID == file2_readID:
                if reads_match(file1_read, file2_read, same):
                    # Found a match. Print, then get new file1 and file2 reads.
                    write_output(bed_writer, file1_file, file2_file, file1_read, file2_read, bedout, bam1out_file, bam2out_file, max_mismatches, ReqFullyAligned)
                    
                    try: file1_read = next(file1_file)
                    except: break  # All done
                    file1_readID = file1_read.query_name
                    
                    try: file2_read = next(file2_file)
                    except: break  # All done
                    file2_readID = file2_read.query_name
                else:
                    # The readID's match, but the read1/read2 combination is wrong.
                    # Maybe the next reads will be the same readID, and have the desired read1/read2 values.
                    # So store the current reads, and look ahead at the next ones.
                    old_file1_read = file1_read
                    old_file2_read = file2_read
                    
                    try:
                        file1_read = next(file1_file)
                    except:
                        file1_eof = True    # We got to the end of file1. Later, we'll break out of this while loop.
                    else:
                        file1_eof = False   # It's OK to keep reading from file1.
                        file1_readID = file1_read.query_name
                    
                    try:
                        file2_read = next(file2_file)
                    except:
                        file2_eof = True     # We got to the end of file2. Later, we'll break out of this while loop.
                    else:
                        file2_eof = False    # It's OK to keep reading from file2.
                        file2_readID = file2_read.query_name
                    
                    if (not file1_eof) and (file1_readID == old_file1_read.query_name):
                        # We were able to get a new read from file1, and it matches the previous file1 readID.
                        # Let's see if it has the desired read1/read2 value.
                        if reads_match(file1_read, old_file2_read, same):
                            write_output(bed_writer, file1_file, file2_file, file1_read, old_file2_read, bedout, bam1out_file, bam2out_file, max_mismatches, ReqFullyAligned)
                            
                            # It did, and we wrote the reads to the csv file.
                            # Now get a new read for the next cycle of the while loop.
                            try:
                                file1_read = next(file1_file)
                            except: 
                                file1_eof = True
                            else:
                                file1_eof = False
                                file1_readID = file1_read.query_name
                    
                    if (not file2_eof) and (file2_readID == old_file2_read.query_name):
                        # We were able to get a new read from file2, and it matches the previous file2 readID.
                        # Let's see if it has the desired read1/read2 value.
                        if reads_match(old_file1_read, file2_read, same):
                            write_output(bed_writer, file1_file, file2_file, old_file1_read, file2_read, bedout, bam1out_file, bam2out_file, max_mismatches, ReqFullyAligned)
                            
                            # It did, and we wrote the reads to the csv file.
                            # Now get a new read for the next cycle of the while loop.
                            try:
                                file2_read = next(file2_file)
                            except:
                                file2_eof = True
                            else:
                                file2_eof = False
                                file2_readID = file2_read.query_name
                    
                    # Somewhere above, we ran out of reads. So it's time to exit the while loop.
                    if file1_eof or file2_eof: break    # All done
            else:
                # The file1 readID is bigger than the file2 readID. Get another file2 line, and check again.
                try:
                    file2_read = next(file2_file)
                except:
                    break  # All done
                file2_readID = file2_read.query_name
    finally:
        file1_file.close()
        file2_file.close()
        
        if bed_file is not None:
            bed_file.close()
        if bam1out is not None:
            bam1out_file.close()
        if bam2out is not None:
            bam2out_file.close()
    
    return


if (__name__ == '__main__'):
    parser = argparse.ArgumentParser(prog = "bam_intersect.py", description = "Given 2 bam files sorted by name (samtools sort, not picard or lexographical sort), produces a bed12 file which pairs read1 from bam file #1 with read2 from bam file #2 (with the same read IDs), and vice versa.", add_help = True)
    
    parser.add_argument('bam1', type=str, help='A full path bam file name, or stdin')
    parser.add_argument('bam2', type=str, help='A full path bam file name, or stdin')
    
    parser.add_argument('--bedout', type=str, default="bamintersect_out.bed", help='File name to output 12 column bed file')
    parser.add_argument('--bam1out', type=str, default=None, help='A full path bam file name')
    parser.add_argument('--bam2out', type=str, default=None, help='A full path bam file name')
    
    parser.add_argument('--reads_match', action='store_true', help='True if looking to match read1/read1, or if the reads are unpaired. False if looking to match read1/read2.')
    parser.add_argument('--max_mismatches', type=int, default=1, help='Maximum number of mismatches a read is allowed to have. The number of mismatches is the value of the read NM tag')
    parser.add_argument('--ReqFullyAligned', action='store_true', help='If set, require reads to be fully aligned.')
    
    parser.add_argument('--src', type=str, help='Full path to source code directory.', required = True) # src directory (needed to locate the C shared library).
    
    parser.add_argument('--version', action='version', version='%(prog)s ' + version)
    
    args = parser.parse_args()
    print("[bamintersect.py] Parameters:", args, file=sys.stderr)
    
    
    # Load our custom strnum_cmp c library:
    try:
        cCode = CDLL(args.src + r"/samtools_strnum_cmp.so")
        cCode.connect()   # Confirm connection was established.
    except:
        raise Exception("[bamintersect.py] Can't connect to samtools_strnum_cmp.so")
    
    
    bam_intersect_f(args.bam1, args.bam2, args.reads_match, args.bedout, args.bam1out, args.bam2out, args.max_mismatches, args.ReqFullyAligned)

