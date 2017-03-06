#!/bin/env python3.5
from sys import argv
import sys
import re
import swalign
import argparse
from leven import levenshtein


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
parser.add_argument("--minBaseQ", action='store', type=int, default=30, help = "The minimum baseQ required over the BC region [%(default)s]. Assumes Phred 33")
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

#minBaseQ
minBaseQ = args.minBaseQ
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
       print("Raw BC seq ", bc_start, "-", bc_end, ": ", referenceSeq[bc_start : bc_end], "\n", file=sys.stderr)
       print("Starting extraction\n\n", file=sys.stderr)


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
numlowQual=0
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
              readname = reads[0].split(' ')[0][1:].split('_')
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
                     propAligned =  alignedLen / (curReadLen - bc_len )
                     
                     offset = alignment.q_pos - alignment.r_pos
                     
                     if alignedLen > 0:
                            alignedQueryWithoutBC = alignment.query[alignment.q_pos:(bc_start + offset)]
                            alignedRefWithoutBC = alignment.ref[alignment.r_pos:bc_start]
                            propMatched = 1 - levenshtein(alignedRefWithoutBC, alignedQueryWithoutBC) / alignedLen
                     else:
                            propMatched = 0
                     
                     if args.verbose:
                            print("Ref: ", alignedRefWithoutBC, " Query: ", alignedQueryWithoutBC, ". ", alignedLen, " bp aligned (", format(propAligned*100, '.0f'), "% of ", curReadLen - bc_len, " bp); ", format(propMatched*100, '.0f'), "% matched; offset=", offset, ".", file=sys.stderr, sep="")
                     
                     if propAligned == 1 and propMatched >= 0.94 and offset <= args.maxOffset and offset >= -args.maxOffset:
                            bc_baseQ = reads[3][(bc_start + offset) : (bc_end + offset)]#Retrieve the baseQ score for BC 
                            if  sum([ord(i)-33 >=minBaseQ for i in bc_baseQ])>=(args.bclen-2):
                            #if all(ord(i)-33 >=minBaseQ for i in bc_baseQ): #Check if all bases have a quality over minBaseQ
                                bc_seq = reads[1][(bc_start + offset) : (bc_end + offset)]
                                barcodeOffsets[offset+args.maxOffset] += 1
                                UMI_Seq=readname[1]
                            else:
                                bc_seq = ""
                                UMI_Seq=""
                                numlowQual += 1
                                if args.verbose:
                                    print("Low baseQ in BC",[ord(i)-33 for i in bc_baseQ], file=sys.stderr, sep="")
                     else:
                            if args.verbose:
                                   print("bad alignment!", file=sys.stderr, sep="")
                            bc_seq = ""
                            bc_baseQ = ""
                            UMI_Seq=""
                            numskipped += 1
              else:
                     offset = 0
                     bc_baseQ = reads[3][(bc_start) : (bc_end)]#ord('A')-33 bc_baseQ = 'EEEEEEEEEEEEAEE'
                     if  sum([ord(i)-33 >=minBaseQ for i in bc_baseQ])>=(args.bclen-2):
                     #if all(ord(i)-33 >=minBaseQ for i in bc_baseQ):
                            bc_seq = reads[1][(bc_start) : (bc_end)]
                            UMI_Seq=readname[1]
                            if args.verbose:
                                print("raw bc seq " + bc_seq, file=sys.stderr)
                     else:
                         bc_seq = ""
                         UMI_Seq=""
                         numlowQual += 1
                         if args.verbose:
                             print("Low baseQ in BC", file=sys.stderr, sep="")


              
              if args.bclen is not None and bc_seq is not "":
                     #Truncate the BC and check if the length of BC is long enough
                     bc_seq = bc_seq[0:args.bclen]

                     if len(bc_seq) < args.bclen:
                            if args.verbose:
                                   print("bc too short! ", file=sys.stderr)
                            bc_seq = ""
                            UMI_Seq=""
                            numskipped += 1
              
              #Print BC, readname, UMI
              #Throw away 2nd part of readname (after space) and get rid of initial @. Split name from UMI (after _)
              
              print(bc_seq, readname[0], UMI_Seq, sep="\t")
              if args.verbose:
                     print("\n", file=sys.stderr, sep="")

except KeyboardInterrupt:
       print("\n\nCaught CTRL+C. Quitting.", sep="", file=sys.stderr)


finally:
       inputfile.close()
       print("\n\nextractBarcode.py statistics:", file=sys.stderr)
       print("Processed ", numread, " reads. Skipped", numskipped, "of these reads.","Number of processed reads with <=2 bp with quality of <",minBaseQ," ",numlowQual, file=sys.stderr)
       print("Missed reads in total ",numskipped+numlowQual, file=sys.stderr)
       print("Percentage kept reads ", format(((numread-(numskipped+numlowQual))/numread)*100, '.2f'), '%', file=sys.stderr)
       if args.align:
              print("barcodeOffsets: ", list(range(-args.maxOffset, args.maxOffset+1)), barcodeOffsets, sep="", file=sys.stderr)


print("\nDone!", file=sys.stderr)
