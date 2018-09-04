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


import sys
if sys.version_info[0] < 3:
    print("Package requires Python 3")
    sys.exit(1)


import re
#import getopt

#for local python2 setup
#sigh
#sys.path.remove("/net/lebowski/vol1/sw/python/2.7.3/lib/python2.7/site-packages/pysam-0.6-py2.7-linux-x86_64.egg")
#sys.path.insert(0, "/home/maurano/.local/lib/python2.7/site-packages/pysam")
import pysam


def parseUMI(read):
    #Fix bwa unmapped bwa reads where MAPQ is incorrectly nonzero
    #fix validation errors: MAPQ should be 0 for unmapped read. prob with chrM reads with one marked as unmapped?
    #eg B6CASTF1J-m.B.cell-DS35978A.C57B6
    #HISEQ-2500-1:79:C62NNANXX:5:2313:20386:22637    609     chrM    16238   29      36M     =    16285   69      CAAAATCATGTTCCGTGAACCAAAACTCTAATCATA     BBBBBFFFFFFFFFFFFFFFFFFFBFFFFFFFFFFF    X0:i:1  X1:i:0  MD:Z:36 PG:Z:MarkDuplicates     XG:i:0  AM:i:29 NM:i:0  SM:i:29 XM:i:0      XO:i:0  XT:A:U
    #HISEQ-2500-1:79:C62NNANXX:5:2313:20386:22637    661     chrM    16285   29      22M13S  =    16238   -69     AATAAACATTAACAAGTTAATGTAGCTTAATAACA FFFFFFBFBFFFBFFFFFFFFFFFF<FFFFBBBBB     MD:Z:22 PG:Z:MarkDuplicates     XG:i:0  AM:i:29 NM:i:0  SM:i:29 XM:i:0  XO:i:0  XT:A:M
    if read.is_unmapped: read.mapping_quality = 0
    
    global totalReads
    totalReads += 1
    
    # strip off the umi, and place it in a custom tag (if it exists)
    
    try:
        umi_loc = read.query_name.index('#')
    
    except:
        pass
    
    else:
        read.setTag("XD", read.query_name[umi_loc+1:])
        read.query_name = read.query_name[:umi_loc]
    
    return read;


def isUnwantedChromosomeName(seqname):
    if re.search('hap', seqname) is not None or re.search('random', seqname) is not None or re.search('chrUn_', seqname) is not None or re.search('_alt', seqname) is not None or re.search('^scaffold', seqname) is not None or re.search('^C\d+', seqname) is not None:
        return True
    else:
        return False


def isUnwantedChromosomeID(reference_id):
    #check for missing mate chromosome name for safety -- sometime bwa forgets to set mate unmapped flag
    if reference_id == -1:
        return True
    
    seqname = unfiltered_reads.getrname(reference_id)
    
    return isUnwantedChromosomeName(seqname)


def filterUnwantedChromosomes(read):
    if not read.is_unmapped:
        if isUnwantedChromosomeID(read.reference_id):
            read.reference_id = -1
            read.reference_start = -1
            read.is_unmapped = True
            read.mapping_quality = 0
            #TODO clear CIGAR, other flags, TLEN
#        else:
#why BUGBUG not working??
#            read.reference_id = filtered_reads.gettid(unfiltered_reads.getrname(read.reference_id))
    if read.is_paired and not read.mate_is_unmapped:
        if isUnwantedChromosomeID(read.next_reference_id):
            read.next_reference_id = -1
            read.next_reference_start = -1
            read.mate_is_unmapped = True
            #TODO set mate_is_reverse, TLEN?
#        else:
#why BUGBUG not working??
#            read.next_reference_id = filtered_reads.gettid(unfiltered_reads.getrname(read.next_reference_id))
    return read;


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


def set_proper_pair(read, mark = True):
    return set_read_flag(read, 1, mark)


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


def validateSingleRead(read):
    if read.mapping_quality < min_MAPQ: raise ReadException("Read must have MAPQ greater than cutoff", str(read.mapping_quality))
    
    #Small number of reads with MAPQ>0 but still unmapped
    if read.is_unmapped: raise ReadException("Read must be mapped")
    
    #needs pysam 0.8.2
    if read.get_tag("NM") > maxNumMismatches: raise ReadException("Too many mismatches", str(read.get_tag("NM")))
    
    #Allow N for introns
    if reqFullyAligned and re.search('[HSPDI]', read.cigarstring) is not None: raise ReadException("Read contains indels or clipping", read.cigarstring)
    
    return;


import argparse


#Can add allow_abbrev=False in python 3.5+
parser = argparse.ArgumentParser(prog = "filter_reads", description = "manually corrects the flags in a single- or pair-end BAM alignment file")
parser.add_argument("raw_alignment", type = str, help = "Input raw alignment file (must be sorted by name")
parser.add_argument("filtered_alignment", type = str, help = "Output filtered alignment file (sorted by name)")
parser.add_argument("--min_mapq", action = "store", type = int, default = 0, help = "Reads must have at least this MAPQ to pass filter [%(default)s]")
parser.add_argument("--reqFullyAligned", action = "store_true", default = False, help = "Require full length of read to align without insertions, deletions, or clipping")
parser.add_argument("--max_mismatches", action = "store", type = int, default = 2, help = "Maximum mismatches to pass filter [%(default)s]")
parser.add_argument("--min_insert_size", action = "store", type = int, default = 0, help = "Minimum insert size to pass filter [%(default)s]")
parser.add_argument("--max_insert_size", action = "store", type = int, default = 500, help = "Maximum insert size to pass filter [%(default)s]")
parser.add_argument("--maxPermittedTrailingOverrun", action = "store", type = int, default = 2, help = "Effectively the number of bp of adapter allowed to be present at end of read (though we don't verify it matches adapter sequence or even mismatches the reference) [%(default)s]")
parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

#argparse does not set an exit code upon this error
#https://stackoverflow.com/questions/5943249/python-argparse-and-controlling-overriding-the-exit-status-code
try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument)
    sys.exit(2)

print("[filter_reads.py] Parameters:", args, file=sys.stderr)


#Default values
verbose = args.verbose
reqFullyAligned = args.reqFullyAligned
maxNumMismatches = args.max_mismatches
min_MAPQ = args.min_mapq
minTemplateLength = args.min_insert_size
maxTemplateLength = args.max_insert_size
maxPermittedTrailingOverrun = args.maxPermittedTrailingOverrun


if verbose and args.filtered_alignment=="-":
    print("Can't do verbose and pipe output to STDOUT", file=sys.stderr)
    sys.exit()


totalReads = 0
totalReadsFailed = 0
readPairFailureCodes = dict()

unfiltered_reads = pysam.AlignmentFile(args.raw_alignment, "rb")

print("[filter_reads.py] Processing header", file=sys.stderr)
newheader = unfiltered_reads.header

newheaderRefSequences = list()
for curSeq in newheader['SQ']:
    if not isUnwantedChromosomeName(curSeq['SN']):
         newheaderRefSequences.append(curSeq)
#why BUGBUG not working??
#newheader['SQ'] = newheaderRefSequences

#TODO add more detail
newheader['PG'].append({'ID':'filter_reads.py', 'PN':'filter_reads.py', 'CL':args})

filtered_reads = pysam.AlignmentFile(args.filtered_alignment, "wbu", header = newheader) #template = unfiltered_reads


print("[filter_reads.py] Processing reads", file=sys.stderr)
read1 = None
# pull the reads
allreads = unfiltered_reads.fetch(until_eof=True) #I believe .next would skip unmapped reads and isn't guaranteed to match file
while(1):
    try:
        if read1 == None:
            read1 = next(allreads)
    except StopIteration:
        break
    read1 = parseUMI(read1)
    read1 = filterUnwantedChromosomes(read1)
    
    try:
        read2 = parseUMI(next(allreads))
        read2 = filterUnwantedChromosomes(read2)
    except StopIteration:
        read2 = None
    
    
    if (read1 and read2) and read1.is_paired and read2.is_paired and read1.query_name == read2.query_name:
        (read1, read2) = (read1, read2) if read1.is_read1 else (read2, read1) 
        
        #Only need to print first name for PE reads
        if verbose: print("PE\t" + read1.query_name, end="")
        
        # filtering code
        try:
            validateSingleRead(read1)
            validateSingleRead(read2)
            
            #PE-specific validation
            if read1.reference_id != read2.reference_id: raise ReadException("Each read must map to the same reference sequence", unfiltered_reads.getrname(read1.reference_id) + "/" + unfiltered_reads.getrname(read2.reference_id))
            
            if read1.is_reverse == read2.is_reverse: raise ReadException("Must be mapped to opposite strands (F-R conformation)", str(read1.is_reverse) + "/" + str(read2.is_reverse))
            
            #Could also check for positive template lengths on reverse reads
            if read1.is_reverse and read1.reference_start < read2.reference_start or read2.is_reverse and read2.reference_start < read1.reference_start: raise ReadException("Reads must face each other", str(read1.reference_start) + "-" + str(read2.reference_start))
            
            if read1.reference_id != read2.next_reference_id: raise ReadException("Read 1: mate claims to map to different reference sequence" + unfiltered_reads.getrname(read1.reference_id) + "/" + unfiltered_reads.getrname(read2.reference_id))
            if read2.reference_id != read1.next_reference_id: raise ReadException("Read 2: mate claims to map to different reference sequence" + unfiltered_reads.getrname(read1.reference_id) + "/" + unfiltered_reads.getrname(read2.reference_id))
            
            #TODO
            #check reference_length > minReadLength
            
            if not abs(read1.template_length) == abs(read2.template_length): raise ReadException("IMPOSSIBLE Template lengths don't match!", str(abs(read1.template_length)) + "/" + str(abs(read2.template_length)))
            
            if abs(read1.template_length) > maxTemplateLength or abs(read2.template_length) > maxTemplateLength: raise ReadException("Each read pair must have an insert length below " + str(maxTemplateLength), str(abs(read1.template_length)) + "/" + str(abs(read2.template_length)))
            
            if abs(read1.template_length) < read1.query_length - maxPermittedTrailingOverrun or abs(read2.template_length) < read2.query_length - maxPermittedTrailingOverrun: raise ReadException("Each read pair must have an insert length above read length (within " + str(maxPermittedTrailingOverrun) + " bp)", str(abs(read1.template_length)) + "/" + str(abs(read2.template_length)))
            
            if abs(read1.template_length) < minTemplateLength or abs(read2.template_length) < minTemplateLength: raise ReadException("Each read pair must have an insert length above " + str(minTemplateLength), str(abs(read1.template_length)) + "/" + str(abs(read2.template_length)))
            
        except ReadException as e:
            if verbose: print("\tFAIL:", e)
            
            # failed a test above, not properly paired
            
            proper_pair = False
            qc_fail = True
            qc_msg = e.msg
        else:
            if verbose: print("\tPASS")
            proper_pair = True
            qc_fail = False
            qc_msg = None
        
        finally:
            set_qc_fail(read1, qc_fail, qc_msg)
            set_proper_pair(read1, proper_pair)
            
            set_qc_fail(read2, qc_fail, qc_msg)
            set_proper_pair(read2, proper_pair)
            
            # write to file
            filtered_reads.write(read1)
            read1 = None
            filtered_reads.write(read2)
            read2 = None
        
    else:
        #Could be single-ended read or PE read without an aligned mate
        if verbose: print("SE\t" + read1.query_name, end="")
        
        # filtering code
        try:
            validateSingleRead(read1)
            
            if read1.is_paired:
                if not read1.mate_is_unmapped:
                    raise ReadException("PE read without mapped mate (though flag says mate was mapped)")
                else:
                    raise ReadException("PE read without mapped mate")
        
        except ReadException as e:
            if verbose: print("\tFAIL:", e)
            qc_fail = True
            qc_msg = e.msg
        
        else:
            if verbose: print("\tPASS")
            qc_fail = False
            qc_msg = None
        
        finally:
            set_qc_fail(read1, qc_fail, qc_msg)
            #can't be properly paired (it's SE)
            set_proper_pair(read1, False)
        
            # write to file
            filtered_reads.write(read1)
        
            read1 = read2
            read2 = None

# clean-up and close files
unfiltered_reads.close()
filtered_reads.close()


#BUGBUG getting lost sys.stderr at end now

print("\n[filter_reads.py] Failure codes by read pair (", totalReads, " reads processed):", sep="", file=sys.stderr)
print("[filter_reads.py]", readPairFailureCodes, file=sys.stderr)


print("[filter_reads.py] Done!", totalReadsFailed, "failed reads", file=sys.stderr)
