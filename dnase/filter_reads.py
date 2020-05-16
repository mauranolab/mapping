#!/bin/env python3

#filter_reads.py - Set SAM flag 0x200 (QC-fail) for reads failing various criteria.
#QC fail is set with |, so does not overwrite existing fail flag.
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
import copy
import pysam
import argparse
import collections


version="1.5"


parser = argparse.ArgumentParser(prog = "filter_reads", description = "manually corrects the flags in a single- or pair-end BAM alignment file", allow_abbrev=False)
parser.add_argument("raw_alignment", type = str, help = "Input raw alignment file (bam; must be sorted by name")
parser.add_argument("filtered_alignment", type = str, help = "Output filtered alignment file (uncompressed bam; sorted by name)")
parser.add_argument("--reheader", action = "store", type = str, help = "Merge in header information from bam or text file (SAM format; alignment section ignored) [%(default)s]")
parser.add_argument("--dropUnmappedReads", action = "store_true", default = False, help = "Omit reads where all segments are unmapped [%(default)s]")
parser.add_argument("--copyCigar", action = "store_true", default = False, help = "Copy CIGAR string to OC tag [%(default)s]")

parser.add_argument("--failUnwantedRefs", action = "store_true", default = False, help = "Reads mapping to unwanted reference sequences will be failed [%(default)s]")
parser.add_argument("--unwanted_refs_list", action = "store", type = str, default = "hap|random|^chrUn_|_alt$|scaffold|^C\d+", help = "Regex defining unwanted reference sequences [%(default)s]")
parser.add_argument("--min_mapq", action = "store", type = int, default = 0, help = "Reads must have at least this MAPQ to pass filter [%(default)s]")
parser.add_argument("--max_mismatches", action = "store", type = float, help = "Maximum mismatches to pass filter, either an integer representing a fixed cutoff or a decimal [0-1] for a variable cutoff as a proportion of read length [%(default)s]")
parser.add_argument("--reqFullyAligned", action = "store_true", default = False, help = "Require full length of read to align without insertions, deletions, or clipping")
parser.add_argument("--min_reference_length", action = "store", type = int, default = 0, help = "Reads must have at least this many reference bp in their alignment [%(default)s]")
parser.add_argument("--min_insert_size", action = "store", type = int, default = 0, help = "Minimum insert size to pass filter [%(default)s]")
parser.add_argument("--max_insert_size", action = "store", type = int, help = "Maximum insert size to pass filter [%(default)s]")
parser.add_argument("--maxPermittedTrailingOverrun", action = "store", type = int, help = "Effectively the number of bp of adapter allowed to be present at end of read (though we don't verify it matches adapter sequence or even mismatches the reference) [%(default)s]")

parser.add_argument("--verbose", action='store_true', default=False, help = "Verbose mode")
parser.add_argument('--version', action='version', version='%(prog)s ' + version)

#argparse does not set an exit code upon this error
#https://stackoverflow.com/questions/5943249/python-argparse-and-controlling-overriding-the-exit-status-code
try:
    args = parser.parse_args()
except argparse.ArgumentError as exc:
    print(exc.message, '\n', exc.argument)
    sys.exit(2)


verbose = args.verbose
if verbose and args.filtered_alignment=="-":
    print("Can't do verbose and pipe output to STDOUT", file=sys.stderr)
    sys.exit()


if args.max_mismatches is not None and args.max_mismatches < 0:
    print("Can't do max_mismatches less than zero", file=sys.stderr)
    sys.exit()


print("[filter_reads.py] Parameters:", args, file=sys.stderr)



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
    if re.search(args.unwanted_refs_list, seqname) is not None:
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
    
    #Strip off the umi and cell barcode, and place them in a custom tag (if it exists).
    #This also loses any info in the read name after the space (In the case of bwa, these are already gone)
    #Too late to retain UMI/cellBC qualities
    #Read name format: [...]_[cell_seq]?_[UMI_seq] [Illumina multiplexing BCs]
    #Hash (#) was the the stam lab read name convention, and we used the XD tag circa 2015 but we haven't ever used either here
    #See SAM optional fields spec at https://samtools.github.io/hts-specs/SAMtags.pdf
    #cellranger uses UR for raw UMI
    readname = read.query_name.split(' ')[0].split('_')
    if(len(readname)>2):
        cell_seq = readname[1]
        UMI_seq = readname[2]
        
        read.setTag("RX", UMI_seq)
        read.setTag("CR", cell_seq) #CB would be corrected
    elif(len(readname)>1):
        UMI_seq = readname[1]
        read.setTag("RX", UMI_seq)
    
    read.query_name = readname[0]
    
    if args.copyCigar:
        read.setTag("OC", read.cigarstring)
    
    return read;


def validateReadPair(read1, read2):
    #Whether a read is unmapped gets checked in validateSingleRead; here we don't want to repeat that here, just check whether fully mapped reads have additional inconsistencies
    if read1.is_unmapped or read2.is_unmapped:
        return
    
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
    if args.reqFullyAligned and read.is_supplementary: raise ReadException("Read is supplementary alignment")
    
    #TODO should we have a separate command line flag to pass unmapped reads? Technically it is orthogonal to MAPQ==0. 
    #Small number of reads with MAPQ>0 but still unmapped
    if read.is_unmapped: raise ReadException("Read must be mapped")
    
    if read.mapping_quality < args.min_mapq: raise ReadException("Read must have MAPQ greater than cutoff", str(read.mapping_quality))
    
    if args.failUnwantedRefs and isUnwantedChromosomeID(read.reference_id): raise ReadException("Maps to unwanted reference", unfiltered_reads.getrname(read.reference_id))
    
    if args.max_mismatches is not None:
        if args.max_mismatches < 1:
            #needs pysam 0.8.2
            if read.get_tag("NM") > args.max_mismatches * read.reference_length: raise ReadException("Too many mismatches", str(read.get_tag("NM")))
        else:
            if read.get_tag("NM") > args.max_mismatches: raise ReadException("Too many mismatches", str(read.get_tag("NM")))
    
    #Allow N for introns
    if args.reqFullyAligned and re.search('[HSPDI]', read.cigarstring) is not None: raise ReadException("Read contains indels or clipping", read.cigarstring)
    
    if read.reference_length < args.min_reference_length: raise ReadException("Read has too few reference bp aligned", read.reference_length)
    
    #Do some PE template length checks in here for simplicity but they are meaningless if the references don't match
    if not read.is_unmapped and not read.mate_is_unmapped and read.is_paired and read.reference_id == read.next_reference_id:
        if args.max_insert_size is not None and abs(read.template_length) > args.max_insert_size: raise ReadException("Each read pair must have an insert length below " + str(args.max_insert_size), str(abs(read.template_length)))
        
        #BUGBUG not sure how this works w/ supp aligns
        if args.maxPermittedTrailingOverrun is not None and abs(read.template_length) < read.query_length - args.maxPermittedTrailingOverrun: raise ReadException("Each read pair must have an insert length above read length (within " + str(args.maxPermittedTrailingOverrun) + " bp)", str(abs(read.template_length)))
        
        if abs(read.template_length) < args.min_insert_size: raise ReadException("Each read pair must have an insert length above " + str(args.min_insert_size), str(abs(read.template_length)))
    
    return;


###Helper functions for merging SAM headers
#Converts a list of dicts to a dict of dicts where the key is the required SAM header fields
commentnum = 0
def samHeaderTagListToDict(header, tag):
    #Globally unique ID for comment lines
    global commentnum
    
    #Deep copy for help in debugging
    if tag in header:
        header = copy.deepcopy(header[tag])
    else:
        header = []
    
    #Don't use dict comprehension since dicts are unordered in python <3.6
    headerDict = collections.OrderedDict()
    for value in header:
        if tag == "SQ":
            key = value['SN'] + "_" + str(value['LN'])
        elif tag == "CO":
            #Use simple incremented to get unique ID so that comments are never merged
            key = value
            commentnum += 1
        else:
            key = value['ID']
        headerDict[key] = value
    return headerDict


def mergeSamHeadersAsDict(header, reheader):
    #Deep copy for help in debugging
    header = copy.deepcopy(header)
    reheader = copy.deepcopy(reheader)
    
    #Don't want to overwrite the main header tag
    if 'HD' in reheader:
        del reheader['HD']
    #Go through @SQ, @RG, @PG, @CO tags individually
    for tag in reheader:
        if tag in reheader:
            #Now we have a list of dicts
            header_tag_dict = samHeaderTagListToDict(header, tag)
            reheader_tag_dict = samHeaderTagListToDict(reheader, tag)
            
            for fieldID in reheader_tag_dict:
                if fieldID in header_tag_dict:
                    #Now we simply have two dicts that can be merged with update
                    header_tag_dict[fieldID].update(reheader_tag_dict[fieldID])
                else:
                    #Otherwise just add it
                    header_tag_dict[fieldID] = reheader_tag_dict[fieldID]
            #Now replace the entire tag list with the updated one from the dict
            header[tag] = list(header_tag_dict.values())
        else:
            #Our bam file does not have any of this tag at all, e.g. comments
            header[tag] = reheader[tag]
    return header



###Initialization and processing
totalReads = 0
totalReadsFailed = 0
totalReadsDropped = 0
readPairFailureCodes = dict()

#newer pysam takes threads argument but doesn't seem to do much
unfiltered_reads = pysam.AlignmentFile(args.raw_alignment, "rb")


print("[filter_reads.py] Processing header", file=sys.stderr)
newheader = unfiltered_reads.header.to_dict()

if args.reheader is not None:
    try:
        if re.search('\.bam$', args.reheader) is not None:
            reheaderFile = pysam.AlignmentFile(args.reheader, "rb")
            reheaderDict = reheaderFile.header.to_dict()
        else:
            #Read in header lines from a text file, in sam format but we ignore the non-header lines
            reheaderFile = open(args.reheader, 'rt')
            headerLineRE = re.compile(r'^@')
            reheaderDict = pysam.AlignmentHeader.from_text(''.join(filter(headerLineRE.search, reheaderFile.readlines()))).to_dict()
        newheader = mergeSamHeadersAsDict(newheader, reheaderDict)
    finally:
        reheaderFile.close()

#Set PP tag to last @PG entry listed (should be bwa, but technically order is not specified), also makes no attempt to follow PG chain in case it is run on a merged bam file.
newheader['PG'].append({'ID':'filter_reads.py', 'PP':newheader['PG'][-1]['ID'], 'PN':'filter_reads.py', 'VN':version, 'CL':' '.join(sys.argv)})

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
            if not read.is_supplementary and not read.is_secondary and read.is_paired:
                #Only validate the primary, non-supplementary reads as a pair
                #TODO should we enforce this for other reads? How can one identify pairs?
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
            
            if args.dropUnmappedReads and allReadsUnmapped:
                totalReadsDropped+=1
            else:
                # write to file
                filtered_reads.write(read)
        
        if args.dropUnmappedReads and allReadsUnmapped:
            if verbose: print(" (dropped)")
        else:
            if verbose: print("")


# clean-up and close files
unfiltered_reads.close()
filtered_reads.close()


print("\n[filter_reads.py] Failure codes by read pair (", totalReads, " reads processed):", sep="", file=sys.stderr)
print("[filter_reads.py]", readPairFailureCodes, file=sys.stderr)


print("[filter_reads.py] Done!", totalReadsFailed, "failed reads,", totalReadsDropped, "reads dropped", file=sys.stderr)
