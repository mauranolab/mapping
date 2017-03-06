#!/bin/env python3.5
from sys import argv
import sys
import re
import swalign
import argparse


#https://stackoverflow.com/questions/25188968/reverse-complement-of-dna-strand-using-python
#TODO switch to biopython
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def revcomp(seq):
       return "".join(complement.get(base, base) for base in reversed(seq))


parser = argparse.ArgumentParser(prog = "extractBarcode", description = "extract barcode from reads based on ref seq", add_help=True)
parser.add_argument('inputfilename', action='store', help='fastq file')
parser.add_argument('--referenceSeq', action='store', default='', help='Ref seq, includes B+ at BC position, should match primer sequence')
parser.add_argument("--align", action='store_true', default=False, help = "Perform SW alignment of read to ref seq")
parser.add_argument("--maxOffset", action='store', type=int, default=3, help = "Max amount of bp shift allowed (only used with --align, must be non-negative) [%(default)s]")
parser.add_argument('--bclen', action='store', type=int, help='length of barcode (barcodes longer than this will be trimmed) [%(default)s]')
parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

#argparse does not set an exit code upon this error
#https://stackoverflow.com/questions/5943249/python-argparse-and-controlling-overriding-the-exit-status-code
try:
       args = parser.parse_args()
except argparse.ArgumentError as exc:
       print(exc.message, '\n', exc.argument)
       sys.exit(2)


print("Extracting barcodes.\nParameters:", file=sys.stderr)
print(args, file=sys.stderr)


if args.inputfilename=="-":
       inputfile = sys.stdin
else:
       inputfile = open(args.inputfilename, 'r') 


#revcomp of R2 sequence
referenceSeq = args.referenceSeq
#referenceSeq = 'CCTTTCTAGGCAATTAGGABBBBBBBBBBBBBBBBCTAGTTGTGGGATCTTTGTCCAAACTCATCGAGCTCGGGA'
referenceSeq = referenceSeq.upper()
referenceSeqLen = len(referenceSeq)
bc_pos = re.search('B+', referenceSeq)
bc_start = bc_pos.start()
bc_end = bc_pos.end()
bc_len = bc_end - bc_start
#If read doesn't get full BC, that's ok, we'll get as much as there is
if args.verbose:
       print("Raw BC seq ", bc_start, "-", bc_end, ": ", referenceSeq[bc_start : bc_end], file=sys.stderr)


if args.align:
       match = 2
       mismatch = -3
       scoring = swalign.NucleotideScoringMatrix(match, mismatch)
       gap = -100
       gapext= -100
       sw = swalign.LocalAlignment(scoring, gap, gapext, wildcard='B')  # you can also choose gap penalties, etc...

barcodeOffsets = [0] * (2*args.maxOffset+1)
numskipped=0
numread=0
try:
       while(True):
              reads = [None] * 4
              reads[0] = inputfile.readline().rstrip('\n') #readname
              #Python has no way to detect EOF?!?
              if reads[0] == "":
                     break
              
              reads[1] = inputfile.readline().rstrip('\n').upper() #sequence
              reads[2] = inputfile.readline().rstrip('\n') #+
              reads[3] = inputfile.readline().rstrip('\n').upper() #baseQ
              
              numread += 1
              
              curReadLen = len(reads[1])
              if curReadLen > referenceSeqLen:
                     raise Exception("Read length (" + str(curReadLen) + ") is longer than ref sequence length (" + str(referenceSeqLen) + ")")
              
              if args.align:
                     #Ref=referenceSeq, Query= reads[1]
                     alignment = sw.align(referenceSeq, reads[1])
                     
                     if args.verbose:
                            print("Query: ", reads[1], file=sys.stderr)
                            print(alignment.dump(), file=sys.stderr)
                     
                     #Don't count barcode length in computing proportion of read that was matched to reference sequence, since it is all Ns.
                     alignedLen = alignment.r_end - alignment.r_pos - bc_len
                     propAligned =  alignedLen / (curReadLen - bc_len)
                     #BUGBUG matches seems wrong wrt wildcard?
                     propMatched = 1#alignment.matches / alignedLen
                     
                     offset = alignment.q_pos - alignment.r_pos
                     
                     if args.verbose:
                            print(alignedLen, " bp (", propAligned *100, "%, read length ", curReadLen, ") aligned; ", propMatched * 100, "% matched; offset=", offset, ".", file=sys.stderr, sep="")
                     
                     if propAligned > 0.8 and propMatched <= 1 and offset <= args.maxOffset and offset >= -args.maxOffset:
                            bc_seq = reads[1][(bc_start + offset) : (bc_end + offset)]
                            barcodeOffsets[offset+args.maxOffset] += 1
                     else:
                            if args.verbose:
                                   print("bad alignment!", file=sys.stderr, sep="")
                            bc_seq = ""
                            numskipped += 1
              else:
                     offset = 0
                     bc_seq = reads[1][(bc_start) : (bc_end)]
                     if args.verbose:
                            print("raw bc seq " + bc_seq, file=sys.stderr)
              
              if args.bclen is not None and bc_seq is not "":
                     bc_seq = bc_seq[0:args.bclen]#Truncate the BC
                     
                     if len(bc_seq) < args.bclen:#If the lenght of the truncated read is less than we want
                            if args.verbose:
                                   print("bc too short! ", file=sys.stderr)#Remove the barcode. Not enough BC to be useful
                            bc_seq = ""
                            numskipped += 1
              
              #Print BC, readname, UMI
              #Throw away 2nd part of readname (after space) and get rid of initial @. Split name from UMI (after _)
              readname = reads[0].split(' ')[0][1:].split('_')
              print(bc_seq, readname[0], readname[1], sep="\t")
              if args.verbose:
                     print("", file=sys.stderr, sep="")

except KeyboardInterrupt:
       print("\n\nCaught CTRL+C. Quitting.", sep="", file=sys.stderr)


finally:
       inputfile.close()
       print("\n\nextractBarcode.py statistics:", file=sys.stderr)
       print("Processed ", numread, " reads. Skipped", numskipped, "of these reads.", file=sys.stderr)
       if args.align:
              print("barcodeOffsets: ", list(range(-args.maxOffset, args.maxOffset+1)), barcodeOffsets, sep="", file=sys.stderr)


print("\nDone!", file=sys.stderr)
