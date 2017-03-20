#!/bin/env python3.5
from sys import argv
import sys
import re
import swalign
import argparse
import io
import gzip
import leven 
from leven import levenshtein
from umi_tools._dedup_umi import edit_distance

parser = argparse.ArgumentParser(prog = "keepPlasmidSeqBCs", description = "Filter fastq.gz files to retain reads with barcodes matching expected sequence in both Barcode read and in plasmid read", add_help=True)
parser.add_argument('--BCread', action='store', help='Read with barcode (fastq.gz file)')
parser.add_argument('--referenceSeq', action='store', default='', help='Reference sequence for BC read including barcode (as B+ at BC position)')
parser.add_argument('--plasmidRead', action='store', default=None, help='Read containing the fixed plasmid (fastq.gz file)')
parser.add_argument('--plasmidSeq', action='store', default=None, help='Plasmid reference sequence')
parser.add_argument("--minBaseQ", action='store', type=int, default=30, help = "The minimum baseQ required over the BC region [%(default)s]. Assumes Phred 33")
parser.add_argument('--bclen', action='store', type=int, help='length of barcode (barcodes longer than this will be trimmed) [%(default)s]')
parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

try:
       args = parser.parse_args()
except argparse.ArgumentError as exc:
       print(exc.message, '\n', exc.argument)
       sys.exit(2)


#inputBCread = gzip.open(args.BCread, 'rt') 
if args.BCread=="-":
       inputBCread = sys.stdin
else:
       inputBCread = gzip.open(args.BCread, 'rt') 

if args.plasmidRead is not None and args.plasmidSeq is None or args.plasmidRead is None and args.plasmidSeq is not None:
       print("Requries both plasmid read and plasmid sequence",file=sys.stderr)
       sys.exit(0)


if args.plasmidSeq is not None and  args.plasmidRead is not None:
       inputplasmidRead = gzip.open(args.plasmidRead, 'rt')
       plasmidSeq = args.plasmidSeq
       plasmidSeq = plasmidSeq.upper()
       plasmidSeqLength = len(plasmidSeq)
       #plasmidSeq = plasmidSeq.encode()
       numWrongPlasmid=0
        
#plasmidSeq='AGCTGCACAGCAACACCGAGCTGGGCATCGTGGAGTACCAGCACGCCTTCAAGACCCCCATCGCCTTCGCCAGATC'

print("Extract barcodes\nParameters:", file=sys.stderr)
print(args, file=sys.stderr)

#Minimum baseq
minBaseQ = args.minBaseQ

#Levenstein Distance
maxEditDist=2

#Reference sequence with Barcode
#referenceSeq = 'CCTTTCTAGGCAATTAGGBBBBBBBBBBBBBBBBCTAGTTGTGGGATCTTTGTCCAAACTCATCGAGCTCGGGA'

referenceSeq = args.referenceSeq
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
       


numskipped=0
numread=0
numWrongBCseq=0
numlowQual=0
numWrongBclength=0

BClevenDist = [0] * (bc_start+1)
PlasmidlevenDist = [0] * (60+1) #TODO limit to read and plasmid length

PlasmidmissmatchLoc = [0] * (60) 
BClmissmatchLoc = [0] * (bc_start)
BCHammingVsLevensthein = 0
PlasmidHammingVsLevensthein = 0
try:
       while(True):
              
              BCread = [None] * 4
              numread += 1
              BCread[0] = inputBCread.readline().rstrip('\n') #readname
              #Python has no way to detect EOF?!?
              if BCread[0] == "":
                     break
              BCread[1] = inputBCread.readline().rstrip('\n').upper() #sequence
              BCread[2] = inputBCread.readline().rstrip('\n') #+
              BCread[3] = inputBCread.readline().rstrip('\n').upper() #baseQ
              
              readname = BCread[0].split(' ')[0][1:].split('_')
              curReadLen = len(BCread[1])
              if curReadLen > referenceSeqLen:
                     raise Exception("Read length (" + str(curReadLen) + ") is longer than ref sequence length (" + str(referenceSeqLen) + ")")
              
              if args.plasmidSeq is not None and  args.plasmidRead is not None:
                     PLread = [None] * 4
                     PLread[0] = inputplasmidRead.readline().rstrip('\n') #readname
                     #Python has no way to detect EOF?!?
                     if PLread[0] == "":
                            break
                     PLread[1] = inputplasmidRead.readline().rstrip('\n').upper() #sequence
                     PLread[2] = inputplasmidRead.readline().rstrip('\n') #+
                     PLread[3] = inputplasmidRead.readline().rstrip('\n').upper() #baseQ
       
              #Filters
              readBCpassed = None
              readPlasmidpassed = None
              readMinbaseQpassed = None
              readLengthpassed = None
              
              bc_baseQ = BCread[3][(bc_start) : (bc_end)]
              
              #Minimum length of plasmid and BC
              PlasmidminLen = min(len(PLread[1]), len(plasmidSeq))
              BCminLen = len(BCread[1][0:bc_start])
              #Find mismatch position
              PlasmidmismatchPosition = [i for i in range(PlasmidminLen) if PLread[1][0:PlasmidminLen][i] != plasmidSeq[0:PlasmidminLen][i]]
              
              BCmismatchPosition = [i for i in range(BCminLen) if BCread[1][0:bc_start][i] != referenceSeq[:bc_start][i]]
             

              BCeditDist = edit_distance(BCread[1][0:bc_start].encode(), referenceSeq[:bc_start].encode())
              plasmidEditDist = edit_distance(PLread[1][0:PlasmidminLen].encode(), plasmidSeq[0:PlasmidminLen].encode())

              BCeditDistHamming = edit_distance(BCread[1][0:bc_start].encode(), referenceSeq[:bc_start].encode())
              plasmidEditDistHamming = edit_distance(PLread[1][0:PlasmidminLen].encode(), plasmidSeq[0:PlasmidminLen].encode())
              
              BClevenDist[BCeditDist]  += 1
              PlasmidlevenDist[plasmidEditDist] += 1
              
              for plasMis in list(map(int, PlasmidmismatchPosition)):
                     PlasmidmissmatchLoc[plasMis] += 1
              
              for BCmis in list(map(int, BCmismatchPosition)):
                     BClmissmatchLoc[BCmis] += 1
#
              #Check if the start of the barcode read matches the plasmid
              if  BCeditDist <= maxEditDist:
                     readBCpassed = True
              else:
                     readBCpassed = False

              if BCeditDistHamming != BCeditDist:
                     BCHammingVsLevensthein +=1
              if plasmidEditDistHamming != plasmidEditDist:
                     PlasmidHammingVsLevensthein +=1
              
              #Check if the plasmid read matches the plasmid
              if args.plasmidSeq is not None and  args.plasmidRead is not None:
                     if plasmidEditDist <= maxEditDist:
                            readPlasmidpassed = True
                     else:
                            readPlasmidpassed = False

                            
              
              
              #Check the baseQ of the barcode
              if sum([ord(i)-33 >=minBaseQ for i in bc_baseQ])>=(args.bclen-2):
                     readMinbaseQpassed = True
              else:
                     readMinbaseQpassed = False
              
              
              if len(BCread[1][(bc_start) : (bc_end)]) == args.bclen:
                     readLengthpassed = True
              else:
                     readLengthpassed = False



              if readBCpassed and readPlasmidpassed and readMinbaseQpassed and readLengthpassed and (args.plasmidSeq is not None and args.plasmidRead is not None or readPlasmidpassed):
                     bc_seq = BCread[1][(bc_start) : (bc_end)]
                     UMI_Seq=readname[1]
                     print(bc_seq, readname[0], UMI_Seq, sep="\t")
              else:
                     numskipped += 1 
                     print("", readname[0], "", sep="\t")
              
              
              #Get statistics of failed reads
              if not readBCpassed:
                     numWrongBCseq += 1
                     if args.verbose:
                            print("Barcode read doesn't match plasmid start. Levenshtein distance is ",BCeditDist, file=sys.stderr, sep="")
              if args.plasmidSeq is not None and  args.plasmidRead is not None:
                     if not readPlasmidpassed:
                            numWrongPlasmid += 1
                            if args.verbose:
                                   print("Plasmid sequence doesn't match. Levenshtein distance is ", plasmidEditDist, file=sys.stderr, sep="")
              if not readMinbaseQpassed:
                     numlowQual += 1
                     if args.verbose:
                            print("Low baseQ in BC", file=sys.stderr, sep="")
              if not readLengthpassed:
                     numWrongBclength += 1 
                     if args.verbose:
                            print("bc not right length! ", file=sys.stderr)

                     

                     
     
except KeyboardInterrupt:
       print("\n\nCaught CTRL+C. Quitting.", sep="", file=sys.stderr)

finally:
       inputBCread.close()
       inputplasmidRead.close()
       print("\n\nextractBarcode.py statistics:", file=sys.stderr)
       print("Processed ", numread, " reads. Skipped", numskipped, "of these reads. ", file=sys.stderr)
       print("Reads with wrong Barcode seq:", numWrongBCseq," (",format(numWrongBCseq/numread*100, '.2f'), '%',")",file=sys.stderr)
       if args.plasmidSeq is not None and  args.plasmidRead is not None:
              print("Reads with wrong Plasmid seq:", numWrongPlasmid," (",format(numWrongPlasmid/numread*100, '.2f'), '%',")",file=sys.stderr)
       print("Reads with >2 bp with minBaseQ < ",minBaseQ,":",numlowQual," (",format(numlowQual/numread*100, '.2f'), '%',")",file=sys.stderr)
       print("Reads with BC length not equal to ",args.bclen,":",numWrongBclength," (",format(numWrongBclength/numread*100, '.2f'), '%',")",file=sys.stderr)
       print("Percentage kept reads ", format(((numread-(numskipped))/numread)*100, '.2f'), '%',file=sys.stderr)
       print("BC Hamming distance: ", list(range(0,bc_start+1)), BClevenDist, sep="", file=sys.stderr)
       print("BC mismatch position: ", list(range(0,bc_start+1)), BClmissmatchLoc, sep="", file=sys.stderr)
       print("Plasmid Hamming distance: ", list(range(1,40+1)), PlasmidlevenDist, sep="", file=sys.stderr)
       print("Plasmid mismatch position: ", list(range(1,40+1)), PlasmidmissmatchLoc, sep="", file=sys.stderr)

print("\nDone!", file=sys.stderr)
