#!/bin/env python
from sys import argv
import sys
import argparse
import io
import gzip

version="1.0"


#TODO
#handle non-gzip
#handle SE
#rescue R1 where only R2 is bad


def log_readfiltering(msg):
    if args.logfile is not None:
        logfile.write(msg + '\n')


parser = argparse.ArgumentParser(prog = "filterNextSeqReadsForPolyG.py", description = "Filter out polyG sequences from the fastq files before trimmomatic as polyG are given high quality scores from NextSeq", allow_abbrev=False)
parser.add_argument('--inputfileR1', type = str, help='input.R1.fastq.gz', required=True)
parser.add_argument('--inputfileR2', type = str, help='input.R2.fastq.gz', required=True)
parser.add_argument('--outputfileR1', type = str, help='output.R1.fastq.gz', required=True)
parser.add_argument('--outputfileR2', type = str, help='output.R2.fastq.gz', required=True)
parser.add_argument('--logfile', type = str, help='logfile.txt.gz - format: R1, R2, R1 pass, R2 pass', required=False)
parser.add_argument("--maxPolyG", type=int, default=75, help = "The minimum percentage of Gs to trim from the sequence -- [0-100] [%(default)s].")
parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s ' + version)


try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument)
    sys.exit(1)


print("Extracting polyG reads.\nParameters:", file=sys.stderr)
print(args, file=sys.stderr)


if args.maxPolyG < 0 or args.maxPolyG>100:
    print("ERROR: --maxPolyG ", args.maxPolyG, " is invalid", file=sys.stderr)
    sys.exit(2)

maxPolyG = args.maxPolyG / 100



numreads=0
polyGreads=0
R1polyG=0
R2polyG=0
R1andR2polyG=0

try:
    #Buffered reader/writer doesn't speed things up at all despite http://aripollak.com/pythongzipbenchmarks/
    inputfileR1 = gzip.open(args.inputfileR1, 'rt') 
    inputfileR2 = gzip.open(args.inputfileR2, 'rt') 

    ## Create output file
    outputfileR1 = gzip.open(args.outputfileR1, 'wt')
    outputfileR2 = gzip.open(args.outputfileR2, 'wt')
    
    if args.logfile is not None:
        logfile = gzip.open(args.logfile, 'wt')

    while(True):
          #Get polyG ids from R1
          readsR1 = [None] * 4
          readsR1[0] = inputfileR1.readline().rstrip('\n') #readname
          #Python has no way to detect EOF?!?
          if readsR1[0] == "":
              break
          readsR1[1] = inputfileR1.readline().rstrip('\n').upper() #sequence
          readsR1[2] = inputfileR1.readline().rstrip('\n') #+
          readsR1[3] = inputfileR1.readline().rstrip('\n').upper() #baseQ
          numreads+=1
          
          
          #Get polyG ids from R2
          readsR2 = [None] * 4
          readsR2[0] = inputfileR2.readline().rstrip('\n') #readname
          #Python has no way to detect EOF?!?
          if readsR2[0] == "":
              break
          readsR2[1] = inputfileR2.readline().rstrip('\n').upper() #sequence
          readsR2[2] = inputfileR2.readline().rstrip('\n') #+
          readsR2[3] = inputfileR2.readline().rstrip('\n').upper() #baseQ
          
          
          R1fail = readsR1[1].count('G')/len(readsR1[1]) >= maxPolyG
          R2fail = readsR2[1].count('G')/len(readsR2[1]) >= maxPolyG
          
          log_readfiltering('\t'.join([readsR1[1], readsR2[1], str(not R1fail), str(not R2fail)]))
          
          if R1fail or R2fail: #75% string is G 
            #print(readsR1[0].split()[0],readsR1[1])
            polyGreads+=1
            if R1fail and not R2fail:
                R1polyG+=1
            if R2fail and not R1fail:
                R2polyG+=1
            if R1fail and R2fail:
                R1andR2polyG+=1
            if args.verbose:
                print('Removing polyG sequence. R1: ' + readsR1[1] + ", R2: " + readsR2[1], file=sys.stderr, sep="")
                #print(readsR1[0].split()[0], file=sys.stderr, sep="")#
          else:
                #print(readsR1[0])
                outputfileR1.write('\n'.join([readsR1[0],readsR1[1],readsR1[2],readsR1[3]])+'\n')
                outputfileR2.write('\n'.join([readsR2[0],readsR2[1],readsR2[2],readsR2[3]])+'\n')



except KeyboardInterrupt:
    print("\n\nCaught CTRL+C. Quitting.", sep="", file=sys.stderr)


finally:
    inputfileR1.close()
    inputfileR2.close()
    outputfileR1.close()
    outputfileR2.close()
    if args.logfile is not None:
        logfile.close()
        
    print("\nfilterNextSeqReadsForPolyG.py statistics:", file=sys.stderr)
    print('Processed read pairs ', numreads,":", "Poly-G read pairs removed in total",polyGreads,":", 'R1 ',R1polyG,': R2 ',R2polyG,': Both ',R1andR2polyG, file=sys.stderr)


print("\nDone!", file=sys.stderr)
