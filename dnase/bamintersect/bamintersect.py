#!/bin/env python

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


def write_output(bed_writer, bam1in_file, bam2in_file, file1_read, file2_read, bedout, bam1out_file, bam2out_file, max_mismatches, ReqFullyAligned):
    if not validateBothReads(max_mismatches, ReqFullyAligned, file1_read, file2_read):
        return
    
    global totalReadpairsOut
    totalReadpairsOut += 1
    
    if bam1out_file is not None:
        bam1out_file.write(file1_read)
    
    if bam2out_file is not None:
        bam2out_file.write(file2_read)
    
    if bed_writer is not None:
        file1_readID = file1_read.query_name
        file2_readID = file2_read.query_name
        
        chrom_1 = bam1in_file.get_reference_name(file1_read.reference_id)
        chrom_2 = bam2in_file.get_reference_name(file2_read.reference_id)
        
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


def openBamOutFile(headertemplate, bamout):
    if bamout is not None:
        newheader = headertemplate.header.to_dict()
        #Set PP tag to last @PG entry listed, also makes no attempt to follow PG chain in case it is run on a merged bam file.
        newheader['PG'].append({'ID':'bamintersect.py', 'PP':newheader['PG'][-1]['ID'], 'PN':'bamintersect.py', 'VN':version, 'CL':' '.join(sys.argv)})
        bamout_file = pysam.AlignmentFile(bamout, "wb", header=newheader)
    else:
        bamout_file = None

    return bamout_file


def getNextRead(bamfile):
    try:
        read = next(bamfile)
    except:
        return None
    return read


def bam_intersect_f(bam_name1, bam_name2, same, bedout, bam1out, bam2out, max_mismatches, ReqFullyAligned):
    # Get iterator handles for bam input files #1 and #2.
    bam1in_file = pysam.AlignmentFile(bam_name1, "rb")
    bam2in_file = pysam.AlignmentFile(bam_name2, "rb")
    
    if bedout is not None:
        bed_file = open(bedout, 'w', newline='')
        bed_writer = csv.writer(bed_file, delimiter='\t', lineterminator=os.linesep, skipinitialspace=True, quotechar='"', quoting=csv.QUOTE_MINIMAL)
    else:
        bed_file = None
        bed_writer = None
    
    bam1out_file = openBamOutFile(bam1in_file, bam1out)
    bam2out_file = openBamOutFile(bam2in_file, bam2out)
    
    ##########################################################################
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
    
    file1_read = getNextRead(bam1in_file)
    file2_read = getNextRead(bam2in_file)
    
    while file1_read is not None and file2_read is not None:
        if samtoolsCmpReadnames(cCode, file1_read.query_name, "<", file2_read.query_name):
            # The file1 readID is less than the file2 readID. Get another file1 read
            file1_read = getNextRead(bam1in_file)
        elif file1_read.query_name == file2_read.query_name:
            if reads_match(file1_read, file2_read, same):
                # Found a match. Print, then get new file1 and file2 reads.
                write_output(bed_writer, bam1in_file, bam2in_file, file1_read, file2_read, bedout, bam1out_file, bam2out_file, max_mismatches, ReqFullyAligned)
                
                file1_read = getNextRead(bam1in_file)
                file2_read = getNextRead(bam2in_file)
            else:
                # The readID's match, but the read1/read2 combination is wrong.
                # look ahead at  the next reads to see if they have the same readID, and have the desired read1/read2 values.
                old_file1_read = file1_read
                old_file2_read = file2_read
                
                file1_read = getNextRead(bam1in_file)
                file2_read = getNextRead(bam2in_file)
                
                if file1_read is not None and file1_read.query_name == old_file1_read.query_name:
                    # We were able to get a new read from file1, and it matches the previous file1 readID.
                    if reads_match(file1_read, old_file2_read, same):
                        write_output(bed_writer, bam1in_file, bam2in_file, file1_read, old_file2_read, bedout, bam1out_file, bam2out_file, max_mismatches, ReqFullyAligned)
                        file1_read = getNextRead(bam1in_file)
                
                if file2_read is not None and file2_read.query_name == old_file2_read.query_name:
                    # We were able to get a new read from file2, and it matches the previous file2 readID.
                    if reads_match(old_file1_read, file2_read, same):
                        write_output(bed_writer, bam1in_file, bam2in_file, old_file1_read, file2_read, bedout, bam1out_file, bam2out_file, max_mismatches, ReqFullyAligned)
                        file2_read = getNextRead(bam2in_file)
        else:
            # The file1 readID is bigger than the file2 readID. Get another file2 line, and check again.
            file2_read = getNextRead(bam2in_file)
    
    bam1in_file.close()
    bam2in_file.close()
    
    if bed_file is not None:
        bed_file.close()
    if bam1out is not None:
        bam1out_file.close()
    if bam2out is not None:
        bam2out_file.close()


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
    
    
    totalReadpairsOut = 0
    
    bam_intersect_f(args.bam1, args.bam2, args.reads_match, args.bedout, args.bam1out, args.bam2out, args.max_mismatches, args.ReqFullyAligned)
    
    print("[bamintersect.py] Wrote", totalReadpairsOut, "read pairs", file=sys.stderr)


