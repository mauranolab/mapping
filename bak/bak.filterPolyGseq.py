#!/bin/env python3.5
from sys import argv
import sys
import re
import swalign
import argparse


parser = argparse.ArgumentParser(prog = "filterPolyGseq", description = "Filter out polyG sequences from the fastq before trimmomatic", add_help=True)
parser.add_argument('--inputR1', action='store', help='R1 fastq file')
parser.add_argument('--inputR2', action='store', help='R1 fastq file')
parser.add_argument("--maxPolyG", action='store', type=int, default=75, help = "The minimum percentage of Gs to trim from the sequence [%(default)s].")
parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-oR1','--outputR1', action='store', dest='outputR1',help='Output file for R1')
parser.add_argument('-oR2','--outputR2', action='store', dest='outputR2',help='Output file for R2')


try:
       args = parser.parse_args()
except argparse.ArgumentError as exc:
       print(exc.message, '\n', exc.argument)
       sys.exit(2)


print("Extracting polyG reads.\nParameters:", file=sys.stderr)
print(args, file=sys.stderr)



inputR1 = open(args.inputR1, 'r') 
inputR2 = open(args.inputR2, 'r') 

maxPolyG = args.maxPolyG

## Create two output files
outfileR1 = open(args.outputR1, 'w')

outfileR2 = open(args.outputR2, 'w')


numread=0
polyGreads=0
try:
       while(True):
                #Get polyG ids from R1
                readsR1 = [None] * 4
                numread+=1
                readsR1[0] = inputR1.readline().rstrip('\n') #readname
                #Python has no way to detect EOF?!?
                if readsR1[0] == "":
                       break
                readsR1[1] = inputR1.readline().rstrip('\n').upper() #sequence
                readsR1[2] = inputR1.readline().rstrip('\n') #+
                readsR1[3] = inputR1.readline().rstrip('\n').upper() #baseQ
                ##Reads R2
                readsR2 = [None] * 4
                readsR2[0] = inputR2.readline().rstrip('\n') #readname
                #Python has no way to detect EOF?!?
                if readsR2[0] == "":
                       break
                readsR2[1] = inputR2.readline().rstrip('\n').upper() #sequence
                readsR2[2] = inputR2.readline().rstrip('\n') #+
                readsR2[3] = inputR2.readline().rstrip('\n').upper() #baseQ
                
                
               #if 'G'*round(len(readsR1[1])*maxPolyG/100) in readsR1[1] or 'G'*round(len(readsR2[1])*maxPolyG/100) in readsR2[1]: #75% string length of consecutive G 
                if readsR1[1].count('G')/len(readsR1[1])*100 >= maxPolyG or readsR2[1].count('G')/len(readsR2[1])*100 >= maxPolyG: #75% string is G 
                  #print(readsR1[0].split()[0],readsR1[1])
                  polyGreads+=1
                  if args.verbose:
                    print('Removing polyG sequence', file=sys.stderr, sep="")
                    #print(readsR1[0].split()[0], file=sys.stderr, sep="")#
                else:
                    #print(readsR1[0])
                    outfileR1.write('\n'.join([readsR1[0],readsR1[1],readsR1[2],readsR1[3]])+'\n')
                    outfileR2.write('\n'.join([readsR2[0],readsR2[1],readsR2[2],readsR2[3]])+'\n')

                    
            
except KeyboardInterrupt:
       print("\n\nCaught CTRL+C. Quitting.", sep="", file=sys.stderr)


finally:
       inputR1.close()
       inputR2.close()
       outfileR1.close()
       outfileR2.close()
       print("\nfilterPolyGseq.py statistics:", file=sys.stderr)
       print('Processed reads ', numread,":", "Poly-G reads removed",polyGreads, file=sys.stderr)
       #if args.align:
       #       print("barcodeOffsets: ", list(range(-args.maxOffset, args.maxOffset+1)), barcodeOffsets, sep="", file=sys.stderr)


print("\nDone!", file=sys.stderr)