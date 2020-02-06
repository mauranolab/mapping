#!/bin/env python3.5
from sys import argv
import sys
import re
import argparse
import io
import gzip
from umi_tools._dedup_umi import edit_distance
import swalign


###Utility functions
#https://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python
#TODO switch to biopython
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def revcomp(seq):
    return "".join(complement.get(base, base) for base in reversed(seq))


def checkReadAgainstRefSequence(readseq, refseq, editDistTotals, mismatchByPos, maxEditDistRate):
    #Permit readthrough into barcode region by excising any BC from both readseq and refseq
#    print("DEBUG raw readseq=", readseq, "; refseq=", refseq, file=sys.stderr, sep="")
    bc_pos_match = re.search('B+', refseq)
    if bc_pos_match is not None:
        bc_start = bc_pos_match.start()
        bc_end = bc_pos_match.end()
        
        readseq = readseq[:bc_start] + readseq[bc_end:]
        refseq = refseq[:bc_start] + refseq[bc_end:]
    
    minLen = min(len(readseq), len(refseq))
    
#    print("DEBUG compare minLen=", minLen, "; readseq=", readseq, "; refseq=", refseq, file=sys.stderr, sep="")
    
    mismatchPosition = [i for i in range(minLen) if readseq[i] != refseq[i]]
    editDist = edit_distance(readseq[0:minLen].encode(), refseq[0:minLen].encode())
    
    readEditDistPassed = editDist <= (maxEditDistRate * minLen)
    
    editDistTotals[editDist] += 1
    
    #BUGBUG doesn't account for offset in --align mode. I guess this should this be relative to ref?
    for mm in list(map(int, mismatchPosition)):
        mismatchByPos[mm] += 1
    
    return editDist, readEditDistPassed


###Fixed parameters
###Command line arguments
version="1.1"

parser = argparse.ArgumentParser(prog = "extractBarcode.py", description = "Filter fastq.gz files to retain reads with barcodes matching expected sequence in both barcode read and in plasmid read", allow_abbrev=False)

bc_parser = parser.add_argument_group()
bc_parser.add_argument('--BCread', action='store', required=True, help='Read with barcode (fastq.gz file or - for stdin)')
bc_parser.add_argument('--bcRefSeq', action='store', required=True, help='Reference sequence for BC read including barcode (as B+ at BC position)')

plasmid_parser = parser.add_argument_group()
plasmid_parser.add_argument('--plasmidRead', action='store', default=None, help='Read containing the fixed plasmid (fastq.gz file)')
plasmid_parser.add_argument('--plasmidRefSeq', action='store', default=None, help='Reference sequence plasmid read')

parser.add_argument("--minBaseQ", action='store', type=int, default=30, help = "The minimum baseQ required over the BC region [%(default)s]. Assumes Phred 33. NB 2 bases are allowed to be below threshold")
parser.add_argument('--bclen', action='store', type=int, required=True, help='length of barcode (barcodes longer than this will be dropped) [%(default)s]')

parser.add_argument("--align", action='store_true', default=False, help = "Perform SW alignment of read to ref seq")
parser.add_argument("--maxOffset", action='store', type=int, default=3, help = "Max amount of bp shift allowed (only used with --align, must be non-negative) [%(default)s]")

BCrevcomp_parser = parser.add_mutually_exclusive_group(required=False)
parser.set_defaults(BCrevcomp=False)
BCrevcomp_parser.add_argument('--BCrevcomp', dest='BCrevcomp', action='store_true', help='BC should be reverse complemented from sequencing read')
BCrevcomp_parser.add_argument('--no-BCrevcomp', dest='BCrevcomp', action='store_false', help='BC should not be reverse complemented from sequencing read')

parser.add_argument("--maxEditDist", action='store', type=float, default=0.1, help = "Max edit (hamming) distance as a proportion of total sequence length. Applied to both BC and plasmid reads. [%(default)s]")

parser.add_argument("--cellBCsuffix", action='store', required=False, help = "Cell BC suffix which will be appended after an underscore, e.g. GEM well number. [%(default)s]")

parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s ' + version)

try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument)
    sys.exit(2)


if args.plasmidRead is None:
    if args.plasmidRefSeq is not None:
        parser.error("Can't specify plasmid sequence without giving plasmid read file")
    else:
        enforcePlasmidRead = False
else:
    if args.plasmidRefSeq is None:
        parser.error("Requires both plasmid read and plasmid sequence")
    else:
        enforcePlasmidRead = True


if args.BCread=="-":
    inputBCread = sys.stdin
else:
    inputBCread = gzip.open(args.BCread, 'rt') 


if enforcePlasmidRead:
    inputPlasmidRead = gzip.open(args.plasmidRead, 'rt')
    plasmidRefSeq = args.plasmidRefSeq
    plasmidRefSeq = plasmidRefSeq.upper()
    plasmidRefSeqLength = len(plasmidRefSeq)
    numWrongPlasmidSeq = 0

if args.cellBCsuffix is None:
    cellBCsuffix = ""
else:
    cellBCsuffix = "_" + args.cellBCsuffix

maxEditDistRate = args.maxEditDist

print("Extract barcodes\nParameters:", file=sys.stderr)
print(args, file=sys.stderr)


#Minimum baseq
minBaseQ = args.minBaseQ


#Reference sequence with Barcode
#bcRefSeq = 'CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBCTAGTTGTGGGATCTTTGTCCAAACTCATCGAGCTCGGGA'
#plasmidRefSeq='AGCTGCACAGCAACACCGAGCTGGGCATCGTGGAGTACCAGCACGCCTTCAAGACCCCCATCGCCTTCGCCAGATC'

bcRefSeq = args.bcRefSeq.upper()
bc_pos_match = re.search('B+', bcRefSeq)
bcRefSeq_bc_start = bc_pos_match.start()
bcRefSeq_bc_end = bc_pos_match.end()
bcRefSeq_bc_len = bcRefSeq_bc_end - bcRefSeq_bc_start
bcRefSeqLength = len(bcRefSeq)

#These will stay 0 unless we are in --align mode
BCreadOffset = 0
bcRefSeqOffset = 0

#If read doesn't get full BC, that's ok, we'll get as much as there is and check it against bclen at the end 
if args.verbose:
    print("Raw BC seq ", bcRefSeq_bc_start, "-", bcRefSeq_bc_end, ": ", bcRefSeq[bcRefSeq_bc_start : bcRefSeq_bc_end], "(", bcRefSeqLength, " nt)\n", file=sys.stderr)


if args.align:
       match = 2
       mismatch = -3
       scoring = swalign.NucleotideScoringMatrix(match, mismatch)
       gap = -100
       gapext= -100
       sw = swalign.LocalAlignment(scoring, gap, gapext, wildcard='B')  # you can also choose gap penalties, etc...
       
       barcodeOffsets = [0] * (2*args.maxOffset+1)


###BC extraction
if args.verbose:
    print("Starting extraction\n\n", file=sys.stderr)

numReads = 0
numReadsSkipped = 0
numReadsNotAligned = 0
numReadsWrongBCseq = 0
numReadsLowBCQual = 0
numReadsWrongBClength = 0
numReadsNsInBC = 0
numReadsNsInUMI = 0

#Initialize lists for worst-case of a read that is same length as the provided reference sequence
BCeditDistTotals = [0] * (bcRefSeqLength+1)
BCmismatchByPos = [0] * bcRefSeqLength

BCqualByPos = [0] * bcRefSeqLength
BCbasesByQual = [0] * 41 ## Assuming a PHRED score max of 40

if enforcePlasmidRead:
    PlasmidEditDistTotals = [0] * (plasmidRefSeqLength+1)
    PlasmidMismatchByPos = [0] * plasmidRefSeqLength

try:
    while(True):
        ###Get read (possibly with pair)
        BCread = [None] * 4
        BCread[0] = inputBCread.readline().rstrip('\n') #readname
        #Python has no way to detect EOF?!?
        if BCread[0] == "":
            break
        BCread[1] = inputBCread.readline().rstrip('\n').upper() #sequence
        BCread[2] = inputBCread.readline().rstrip('\n') #+
        BCread[3] = inputBCread.readline().rstrip('\n').upper() #baseQ
        
        if enforcePlasmidRead:
            PLread = [None] * 4
            PLread[0] = inputPlasmidRead.readline().rstrip('\n') #readname
            #Python has no way to detect EOF?!?
            if PLread[0] == "":
                break
            PLread[1] = inputPlasmidRead.readline().rstrip('\n').upper() #sequence
            PLread[2] = inputPlasmidRead.readline().rstrip('\n') #+
            PLread[3] = inputPlasmidRead.readline().rstrip('\n').upper() #baseQ
        
        #Increment counter only after successfully reading both BC and plasmid read
        numReads += 1
        
        
        ###Filters
        readAlignmentPassed = None
        readBCreadEditDistPassed = None
        readPlasmidReadEditDistPassed = None
        readBCMinBaseQpassed = None
        readBClengthPassed = None
        readNoNsInBCPassed = None
        readNoNsInUMIPassed = None
        
        
        #Optionally shift the BCsequence after alignment to reference
        #Don't think it would make sense to combine this with --plasmidRead / --plasmidRefSeq although it is not prohibited
        if args.align:
            alignment = sw.align(bcRefSeq, BCread[1])
            
            #Set this here, will be updated below in case we adjust BCread
            curReadLen = len(BCread[1])
            
            #Don't count barcode length in computing proportion of read that was matched to reference sequence.
            alignedLenLeftOfBC = max(bcRefSeq_bc_start - alignment.r_pos, 0)
            alignedLenRightOfBC = max(alignment.r_end - bcRefSeq_bc_end, 0)
            #alignedLen = alignment.r_end - alignment.r_pos - bcRefSeq_bc_len
            alignedLen = alignedLenLeftOfBC + alignedLenRightOfBC
            propAligned =  alignedLen / (curReadLen - bcRefSeq_bc_len)
            offset = alignment.q_pos - alignment.r_pos
            
            if alignedLen >= 20 and alignedLenLeftOfBC >= 5 and alignedLenRightOfBC >= 5 and propAligned > 0.8 and offset <= args.maxOffset and offset >= -args.maxOffset:
                if offset > 0:
                    BCreadOffset = offset
                    bcRefSeqOffset = 0
                else:
                    BCreadOffset = 0
                    bcRefSeqOffset = -1 * offset
                
#                print("DEBUG BCreadOffset=", BCreadOffset, "; bcRefSeqOffset=", bcRefSeqOffset, file=sys.stderr, sep="")
                
                barcodeOffsets[offset+args.maxOffset] += 1
                readAlignmentPassed = True
            else:
                readAlignmentPassed = False
            
            if args.verbose:
                print("Query: ", BCread[1], file=sys.stderr)
                print(alignment.dump(), file=sys.stderr)
                print("Alignment results: ", alignedLen, " total nt aligned (", propAligned *100, "%, read length ", curReadLen, "): ", alignedLenLeftOfBC, " nt left of BC, ", alignedLenRightOfBC, " nt right of BC; alignment offset=", offset, "; BCreadOffset=", BCreadOffset, ", bcRefSeqOffset=", bcRefSeqOffset, ".", file=sys.stderr, sep="")
        
        
        BCeditDist, readBCreadEditDistPassed = checkReadAgainstRefSequence(BCread[1][BCreadOffset:], bcRefSeq[bcRefSeqOffset:], BCeditDistTotals, BCmismatchByPos, maxEditDistRate)
        if enforcePlasmidRead:
            plasmidEditDist, readPlasmidReadEditDistPassed = checkReadAgainstRefSequence(PLread[1], plasmidRefSeq, PlasmidEditDistTotals, PlasmidMismatchByPos, maxEditDistRate)
        
        
        #Extract BC and check length
        bc_seq = BCread[1][(bcRefSeq_bc_start-bcRefSeqOffset) : (bcRefSeq_bc_end-bcRefSeqOffset)]
        #print("DEBUG found BC:", bc_seq, file=sys.stderr, sep="")
        if len(bc_seq) == args.bclen:
            readBClengthPassed = True
        else:
            readBClengthPassed = False
        
        if(args.BCrevcomp):
            bc_seq = revcomp(bc_seq)
        
        if re.search('[^ACGT]', bc_seq) is not None:
            readNoNsInBCPassed = False
        else:
            readNoNsInBCPassed = True
        
        #Check the baseQ of the barcode
        bc_baseQ = BCread[3][(bcRefSeq_bc_start-bcRefSeqOffset) : (bcRefSeq_bc_end-bcRefSeqOffset)]
        bc_baseQ = [ord(i)-33 for i in bc_baseQ]
        for i in range(len(bc_baseQ)):
            BCqualByPos[i] += bc_baseQ[i]
            BCbasesByQual[min(bc_baseQ[i], 40)] += 1
        
        #NB Permits 2 bases to be below threshold
        if sum([i >= minBaseQ for i in bc_baseQ]) >= args.bclen-2:
            readBCMinBaseQpassed = True
        else:
            readBCMinBaseQpassed = False
        
        #TODO 10x check cell BC baseq
        
        ###Extract UMI and cell BC if present
        #Read name format: [@...]_[cell_seq]?_[UMI_seq] [Illumina multiplexing BCs]
        readname = BCread[0].split(' ')[0][1:].split('_')
        UMI_seq = ""
        cell_seq = ""
        if(len(readname)>2):
            cell_seq = readname[1] + cellBCsuffix
            UMI_seq = readname[2]
        elif(len(readname)>1):
            UMI_seq = readname[1]
        
        if UMI_seq!="" and re.search('[^ACGT]', UMI_seq) is not None: #re will not find any thing in "" but just to be explicit
            readNoNsInUMIPassed = False
        else:
            readNoNsInUMIPassed = True
        
        
        if readBCreadEditDistPassed and readBCMinBaseQpassed and readBClengthPassed and readNoNsInBCPassed and (not enforcePlasmidRead or readPlasmidReadEditDistPassed) and (not args.align or readAlignmentPassed):
            #Read is good
            print(bc_seq, readname[0], UMI_seq, cell_seq, sep="\t")
        else:
            numReadsSkipped += 1
            print("", readname[0], "", "", sep="\t")
        
        
        ###Update statistics of failed reads
        if not readAlignmentPassed:
            numReadsNotAligned += 1
            if args.verbose:
                print("Bad alignment!", file=sys.stderr, sep="")
        if not readBCreadEditDistPassed:
            numReadsWrongBCseq += 1
            if args.verbose:
                print("Barcode read exceeds distance threshold; Hamming distance is ", BCeditDist, sep="", file=sys.stderr)
        if enforcePlasmidRead and not readPlasmidReadEditDistPassed:
            numWrongPlasmidSeq += 1
            if args.verbose:
                print("Plasmid sequence exceeds distance threshold;  Hamming distance is ", plasmidEditDist, sep="", file=sys.stderr)
        if not readBClengthPassed:
            numReadsWrongBClength += 1
            if args.verbose:
                print("BC not right length! ", file=sys.stderr)
        if not readNoNsInBCPassed:
            numReadsNsInBC += 1
            if args.verbose:
                print("Non-ACGT bases in BC", sep="", file=sys.stderr)
        if not readNoNsInUMIPassed:
            numReadsNsInUMI += 1
            if args.verbose:
                print("Non-ACGT bases in UMI", sep="", file=sys.stderr)
        if not readBCMinBaseQpassed:
            numReadsLowBCQual += 1
            if args.verbose:
                print("Low baseQ in BC", sep="", file=sys.stderr)


except KeyboardInterrupt:
    print("\n\nCaught CTRL+C. Quitting.", sep="", file=sys.stderr)

finally:
    inputBCread.close()
    if enforcePlasmidRead:
        inputPlasmidRead.close()
    
    print("\nextractBarcode.py statistics:", file=sys.stderr)
    print("Processed", numReads, "reads. Skipped", numReadsSkipped, "of these reads:", file=sys.stderr)
    if args.align:
        print("    Reads not aligned: ", numReadsNotAligned, " (", format(numReadsNotAligned/numReads*100, '.2f'), "%)", sep="", file=sys.stderr)
    print("    Reads with barcode read sequence exceeding distance threshold: ", numReadsWrongBCseq, " (", format(numReadsWrongBCseq/numReads*100, '.2f'), "%)", sep="", file=sys.stderr)
    if enforcePlasmidRead:
        print("    Reads with plasmid read sequence exceeding distance threshold: ", numWrongPlasmidSeq, " (", format(numWrongPlasmidSeq/numReads*100, '.2f'), "%)", sep="", file=sys.stderr)
    print("    Reads with >2 bp in BC with minBaseQ <", minBaseQ, ": ", numReadsLowBCQual, " (", format(numReadsLowBCQual/numReads*100, '.2f'), "%)", sep="", file=sys.stderr)
    print("    Reads with BC length not equal to ", args.bclen, ": ", numReadsWrongBClength, " (", format(numReadsWrongBClength/numReads*100, '.2f'), "%)", sep="", file=sys.stderr)
    print("    Reads with non-ACGT bases in BC: ", numReadsNsInBC, " (", format(numReadsNsInBC/numReads*100, '.2f'), "%)", sep="", file=sys.stderr)
    print("    Reads with non-ACGT bases in UMI: ", numReadsNsInUMI, " (", format(numReadsNsInUMI/numReads*100, '.2f'), "%)", sep="", file=sys.stderr)
    print("Percentage reads kept: ", format(((numReads-(numReadsSkipped))/numReads)*100, '.2f'), '%', sep="", file=sys.stderr)
    print("\nBC read edit distances: ", list(range(0, len(BCeditDistTotals))), BCeditDistTotals, sep="", file=sys.stderr)
    print("BC read mismatches by position: ", list(range(0, len(BCmismatchByPos))), BCmismatchByPos, sep="", file=sys.stderr)
    print("BC avg. baseQ by position: ", list(range(0, len(BCqualByPos))), BCqualByPos, sep="", file=sys.stderr)
    print("BC baseQ frequency: ", list(range(0, len(BCbasesByQual))), BCbasesByQual, sep="", file=sys.stderr)

    if enforcePlasmidRead:
        print("Plasmid read edit distances: ", list(range(0, len(PlasmidEditDistTotals))), PlasmidEditDistTotals, sep="", file=sys.stderr)
        print("Plasmid read mismatches by position: ", list(range(0, len(PlasmidMismatchByPos))), PlasmidMismatchByPos, sep="", file=sys.stderr)
    if args.align:
        print("BC read offsets: ", list(range(-args.maxOffset, args.maxOffset+1)), barcodeOffsets, sep="", file=sys.stderr)
