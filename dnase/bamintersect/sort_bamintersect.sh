#!/bin/bash
set -eu -o pipefail

bam=${1}
inputsfile=${2}
bam_keep_flags=${3} # Required. Use "0" if necessary
bam_exclude_flags=${4} # Required. Use "0" if necessary
uninformativeRegionFile=${5}
OUTBASE=${6}
src=${7}


jobid=$SLURM_ARRAY_TASK_ID
jobname=`awk -F "\t" -v jobid=$jobid 'NR==jobid {print $1}' ${inputsfile}`
chroms=`awk -F "\t" -v jobid=$jobid 'NR==jobid {print $2}' ${inputsfile}`

bam_out="${OUTBASE}.${jobname}.bam"


echo "Running on $HOSTNAME. Using $TMPDIR as tmp"
echo "Sorting ${chroms} of ${bam} to ${bam_out}"
date
echo


#Likely more IO efficient to cache the subsetted BAM locally
samtools view -@ ${NSLOTS} -f ${bam_keep_flags} -F ${bam_exclude_flags} -O BAM ${bam} ${chroms} > $TMPDIR/sort.${jobname}.bam

#Identify read names overlapping the uninformativeRegionFile for removal
#We could implement bed input directly in subsetBAM.py to avoid a second pass, but I think we would need to require the file be sorted
samtools view -L ${uninformativeRegionFile} $TMPDIR/sort.${jobname}.bam -O SAM | cut -f1 | uniq > $TMPDIR/uninformativeReads.txt

## Sort bam file by read name.  Done by samtools via strnum_cmp.c
#     Query names are split into alternating subfields of pure nondigits and pure digits.
#     
#     To sort a pair of query names, the corresponding subfields of each query name are compared sequentially:
#         When a nondigit subfield is compared to a nondigit subfield, the sort is lexicographic.
#         When a all-digit subfield is compared to a all-digit subfield, the sort is numeric.
#         When a all-digit subfield is compared to a nondigit subfield, the sort is lexicographic.
#       
#     Lexicographic sorts by left justifying both subfields, and comparing them from left to right, character by character.
#     The character comparisons are made by ASCII value.
#     
#     Numeric sorts do not entail character by character comparisons.  The entire subfield is considered as one integer, and
#     the integer values of the two subfields are compared.  Leading zeroes are only used as tie breakers. The subfield with
#     more leading zeroes is placed before the subfield with fewer leading zeroes. 

#Since the next set of jobs will be reading this heavily, it's likely worthwhile to use high compression for big bam files; also note that the sort is likely limiting
${src}/subsetBAM.py --exclude_readnames $TMPDIR/uninformativeReads.txt $TMPDIR/sort.${jobname}.bam - |
samtools sort -@ $NSLOTS -O bam -m 4000M -T $TMPDIR/sortbyname -l 8 -n -o ${bam_out}

echo "Done!!!"
date
