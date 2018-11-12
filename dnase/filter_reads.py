#!/bin/env python3

#filter_reads.py - Set SAM flag 0x200 (QC-fail) for reads failing various criteria.
#
#See also fork in Github StamLab/stampipes
#
#Criteria:
#* Read must be end-to-end mapped without soft clipping or indels and meet cutoffs for MAPQ and NM
#* PE reads must have both reads present and meeting all above criteria. Additionally, they must face each other on the same reference and have an insert length within an upper and lower range cutoff (note the lower cutoff is bounded by read length but can be increased for additional stringency).
#
#Requirements:
#pySAM 0.8.2 or higher
#
#
#Input: assumed to be BAM for now
#
#Usage:
#Must be sorted by read name
#
#filter_reads.py in.bam out.bam
#
#or 
#
#samtools sort -@ $NSLOTS -O bam -T $TMPDIR/${sample}.sortbyname -l 1 -n - |
#filter_reads.py - - |
#samtools sort -@ $NSLOTS -O bam -T $TMPDIR/${sample}.sort -l 1 - > $sample.bam


#BUGBUG not handling setting unmapped reads when there are supplementary alignments, eg bwa mem. Fundamental problem as processing reads in pairs right now, would need to go through multiple lines each for R1 and R2.

import sys
if sys.version_info[0] < 3:
    print("Package requires Python 3")
    sys.exit(1)

import re
#import getopt
import pysam
import argparse


version="1.1"


parser = argparse.ArgumentParser(prog = "filter_reads", description = "manually corrects the flags in a single- or pair-end BAM alignment file", allow_abbrev=False)
parser.add_argument("raw_alignment", type = str, help = "Input raw alignment file (bam; must be sorted by name")
parser.add_argument("filtered_alignment", type = str, help = "Output filtered alignment file (uncompressed bam; sorted by name)")
parser.add_argument("--dropUnmappedReads", action = "store_true", default = False, help = "Omit reads where all segments are unmapped [%(default)s]")
parser.add_argument("--failUnwantedRefs", action = "store_true", default = False, help = "Reads mapping to unwanted reference sequences will be failed [%(default)s]")
parser.add_argument("--unwanted_refs_list", action = "store", type = str, default = "hap|random|^chrUn_|_alt$|scaffold|^C\d+", help = "Regex defining unwanted reference sequences [%(default)s]")
parser.add_argument("--min_mapq", action = "store", type = int, default = 0, help = "Reads must have at least this MAPQ to pass filter [%(default)s]")
parser.add_argument("--reqFullyAligned", action = "store_true", default = False, help = "Require full length of read to align without insertions, deletions, or clipping")
parser.add_argument("--max_mismatches", action = "store", type = float, default = 2, help = "Maximum mismatches to pass filter [%(default)s]")
parser.add_argument("--min_insert_size", action = "store", type = int, default = 0, help = "Minimum insert size to pass filter [%(default)s]")
parser.add_argument("--max_insert_size", action = "store", type = int, default = 500, help = "Maximum insert size to pass filter [%(default)s]")
parser.add_argument("--maxPermittedTrailingOverrun", action = "store", type = int, default = 2, help = "Effectively the number of bp of adapter allowed to be present at end of read (though we don't verify it matches adapter sequence or even mismatches the reference) [%(default)s]")
parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s ' + version)

#argparse does not set an exit code upon this error
#https://stackoverflow.com/questions/5943249/python-argparse-and-controlling-overriding-the-exit-status-code
try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument)
    sys.exit(2)


if args.verbose and args.filtered_alignment=="-":
    print("Can't do verbose and pipe output to STDOUT", file=sys.stderr)
    sys.exit()


if args.max_mismatches < 0:
    print("Can't do max_mismatches less than zero", file=sys.stderr)
    sys.exit()


print("[filter_reads.py] Parameters:", args, file=sys.stderr)

#Default values
verbose = args.verbose
reqFullyAligned = args.reqFullyAligned
maxNumMismatches = args.max_mismatches
min_MAPQ = args.min_mapq
minTemplateLength = args.min_insert_size
maxTemplateLength = args.max_insert_size
maxPermittedTrailingOverrun = args.maxPermittedTrailingOverrun
failUnwantedRefs = args.failUnwantedRefs
unwanted_refs_list = args.unwanted_refs_list
dropUnmappedReads = args.dropUnmappedReads



class ReadException(Exception):
    def __str__(self):
        if self.value == None:
            return self.msg + "\t"
        else:
            return self.msg + "\t" + self.value
    def __init__(self, msg, value=None):
        super(ReadException, self).__init__(msg)
        self.msg = msg
        self.value = value


'''
General function to set the flag field
'''
def set_read_flag(read, flag, mark):
    if mark:
        read.flag |= (1<<flag)
    else:
        read.flag &= ~(1<<flag)
    return read


def set_qc_fail(read, mark = True, qc_msg = None):
    global readPairFailureCodes
    global totalReadsFailed
    
    if(mark):
        if qc_msg in readPairFailureCodes:
            readPairFailureCodes[qc_msg] += 1
        else:
            readPairFailureCodes[qc_msg] = 1
    
        totalReadsFailed += 1
    
    return set_read_flag(read, 9, mark)


def isUnwantedChromosomeName(seqname):
    if re.search(unwanted_refs_list, seqname) is not None:
        return True
    else:
        return False


def isUnwantedChromosomeID(reference_id):
    seqname = unfiltered_reads.getrname(reference_id)
    return isUnwantedChromosomeName(seqname)


def parseRead(read):
    #Fix bwa unmapped bwa reads where MAPQ is incorrectly nonzero
    #fix validation errors: MAPQ should be 0 for unmapped read. prob with chrM reads with one marked as unmapped?
    #eg B6CASTF1J-m.B.cell-DS35978A.C57B6
    #HISEQ-2500-1:79:C62NNANXX:5:2313:20386:22637    609     chrM    16238   29      36M     =    16285   69      CAAAATCATGTTCCGTGAACCAAAACTCTAATCATA     BBBBBFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFF    X0:i:1  X1:i:0  MD:Z:36 PG:Z:MarkDuplicates     XG:i:0  AM:i:29 NM:i:0  SM:i:29 XM:i:0      XO:i:0  XT:A:U
    #HISEQ-2500-1:79:C62NNANXX:5:2313:20386:22637    661     chrM    16285   29      22M13S  =    16238   -69     AATAAACATTAACAAGTTAATGTAGCTTAATAACA FFFFFFBFBFFFBFFFFFFFFFFFF<FFFFBBBBB     MD:Z:22 PG:Z:MarkDuplicates     XG:i:0  AM:i:29 NM:i:0  SM:i:29 XM:i:0  XO:i:0  XT:A:M
    if read.is_unmapped: read.mapping_quality = 0
    
    global totalReads
    totalReads += 1
    
    return parseUMI(read)

def parseUMI(read):
    # strip off the umi, and place it in a custom tag (if it exists)
    
    try:
        #I think # is the the stam lab convention
        #/vol/mauranolab/flowcells/fastq/FCH3YVFAFXY/umi/demuxReadsByContent.py and /home/mauram01/.local/bin/umi_tools dedup expect them to be separated from read name by _ (but is configurable)
        #I think bcl2fastq uses :
        umi_loc = read.query_name.index('#')
    
    except:
        pass
    
    else:
        #Picard seems to expect RX tag by default https://broadinstitute.github.io/picard/command-line-overview.html#UmiAwareMarkDuplicatesWithMateCigar
        #UM?
        read.setTag("XD", read.query_name[umi_loc+1:])
        read.query_name = read.query_name[:umi_loc]
    
    return read;


def validateReadPair(read1, read2):
    if read1.reference_id != read2.reference_id: raise ReadException("Each read must map to the same reference sequence", unfiltered_reads.getrname(read1.reference_id) + "/" + unfiltered_reads.getrname(read2.reference_id))
    
    if read1.is_reverse == read2.is_reverse: raise ReadException("Must be mapped to opposite strands (F-R conformation)", str(read1.is_reverse) + "/" + str(read2.is_reverse))
    
    #Could also check for positive template lengths on reverse reads
    if read1.is_reverse and read1.reference_start < read2.reference_start or read2.is_reverse and read2.reference_start < read1.reference_start: raise ReadException("Reads must face each other", str(read1.reference_start) + "-" + str(read2.reference_start))
    
    #Verify that the flags of one read accurately match the next segment
    if read1.reference_id != read2.next_reference_id: raise ReadException("IMPOSSIBLE Read 1: mate disagrees on reference sequence" + unfiltered_reads.getrname(read1.reference_id) + "/" + unfiltered_reads.getrname(read2.reference_id))
    if read2.reference_id != read1.next_reference_id: raise ReadException("IMPOSSIBLE Read 2: mate disagrees on different reference sequence" + unfiltered_reads.getrname(read1.reference_id) + "/" + unfiltered_reads.getrname(read2.reference_id))
    
    if read1.mate_is_unmapped != read2.is_unmapped: raise ReadException("IMPOSSIBLE Read 1: mate disagrees on map/unmapped status")
    if read2.mate_is_unmapped != read1.is_unmapped: raise ReadException("IMPOSSIBLE Read 2: mate disagrees on map/unmapped status")
    
    if read1.mate_is_reverse != read2.is_reverse: raise ReadException("IMPOSSIBLE Read 1: mate disagrees on strand")
    if read1.mate_is_reverse != read2.is_reverse: raise ReadException("IMPOSSIBLE Read 2: mate disagrees on strand")
    
    if not abs(read1.template_length) == abs(read2.template_length): raise ReadException("IMPOSSIBLE Template lengths don't match!", str(abs(read1.template_length)) + "/" + str(abs(read2.template_length)))
    
    return;

def validateSingleRead(read):
    if read.mapping_quality < min_MAPQ: raise ReadException("Read must have MAPQ greater than cutoff", str(read.mapping_quality))
    
    #Small number of reads with MAPQ>0 but still unmapped
    if read.is_unmapped: raise ReadException("Read must be mapped")
    
    if maxNumMismatches < 1:
        #needs pysam 0.8.2
        if read.get_tag("NM") > maxNumMismatches * read.reference_length: raise ReadException("Too many mismatches", str(read.get_tag("NM")))
    else:
        if read.get_tag("NM") > maxNumMismatches: raise ReadException("Too many mismatches", str(read.get_tag("NM")))
    
    if failUnwantedRefs and isUnwantedChromosomeID(read.reference_id):
        raise ReadException("Maps to unwanted reference", unfiltered_reads.getrname(read.reference_id))
    
    if reqFullyAligned and read.is_supplementary: raise ReadException("Read is supplementary alignment")
    
    #Allow N for introns
    if reqFullyAligned and re.search('[HSPDI]', read.cigarstring) is not None: raise ReadException("Read contains indels or clipping", read.cigarstring)
    
    #Do some PE template length checks in here for simplicity but they are meaningless if the references don't match
    if read.is_paired and read.reference_id == read.next_reference_id:
        #TODO check reference_length > minReadLength
        
        if abs(read.template_length) > maxTemplateLength: raise ReadException("Each read pair must have an insert length below " + str(maxTemplateLength), str(abs(read.template_length)))
        
        #BUGBUG not sure how this works w/ supp aligns
        if abs(read.template_length) < read.query_length - maxPermittedTrailingOverrun: raise ReadException("Each read pair must have an insert length above read length (within " + str(maxPermittedTrailingOverrun) + " bp)", str(abs(read.template_length)))
        
        if abs(read.template_length) < minTemplateLength: raise ReadException("Each read pair must have an insert length above " + str(minTemplateLength), str(abs(read.template_length)))
    
    return;



###Initialization and processing
totalReads = 0
totalReadsFailed = 0
totalReadsDropped = 0
readPairFailureCodes = dict()

#TODO newer pysam takes threads argument but doesn't seem to do much
unfiltered_reads = pysam.AlignmentFile(args.raw_alignment, "rb")

print("[filter_reads.py] Processing header", file=sys.stderr)
newheader = unfiltered_reads.header.to_dict() #to_dict() is for the changes from version 0.14
newheader['PG'].append({'ID':'filter_reads.py', 'PN':'filter_reads.py', 'VN':version, 'CL':args})
filtered_reads = pysam.AlignmentFile(args.filtered_alignment, "wbu", header = newheader)


print("[filter_reads.py] Processing reads", file=sys.stderr)
allreads = unfiltered_reads.fetch(until_eof=True) #I believe .next would skip unmapped reads and isn't guaranteed to match file
nextread = parseRead(next(allreads))
while(1):
    if nextread is None:
        break
    
    curreads = []
    while(len(curreads)==0 or curreads[0].query_name == nextread.query_name):
        curreads.append(nextread)
        nextread = None
        try:
            nextread = parseRead(next(allreads))
        except StopIteration:
            break
    
    # filtering code
    try:
        read1=None
        read2=None
        allReadsUnmapped=True
        for read in curreads:
            if not read.is_supplementary and read.is_paired:
                #Only validate the non-supplementary reads as a pair
                if read.is_read1:
                    read1=read
                elif read.is_read2:
                    read2=read
            if not read.is_unmapped:
                allReadsUnmapped=False
        
        if read1 and read2:
            if verbose: print("PE (" + str(len(curreads)) + ")\t" + curreads[0].query_name, end="")
        else:
            if verbose: print("SE\t" + curreads[0].query_name, end="")
        
        for read in curreads:
            validateSingleRead(read)
        
        if read1 and read2:
            validateReadPair(read1, read2)
        
    except ReadException as e:
        if verbose: print("\tFAIL:", e, end="")
        qc_fail = True
        qc_msg = e.msg
    else:
        if verbose: print("\tPASS", end="")
        qc_fail = False
        qc_msg = None
    
    finally:
        for read in curreads:
            set_qc_fail(read, qc_fail, qc_msg)
            
            if dropUnmappedReads and allReadsUnmapped:
                totalReadsDropped+=1
            else:
                # write to file
                filtered_reads.write(read)
        
        if dropUnmappedReads and allReadsUnmapped:
            if verbose: print(" (dropped)")
        else:
            if verbose: print("")


# clean-up and close files
unfiltered_reads.close()
filtered_reads.close()


#BUGBUG getting lost sys.stderr at end now

print("\n[filter_reads.py] Failure codes by read pair (", totalReads, " reads processed):", sep="", file=sys.stderr)
print("[filter_reads.py]", readPairFailureCodes, file=sys.stderr)


print("[filter_reads.py] Done!", totalReadsFailed, "failed reads, ", totalReadsDropped, "reads dropped", file=sys.stderr)
