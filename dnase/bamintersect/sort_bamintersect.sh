#!/bin/bash
set -eu -o pipefail

########################################################
## Variables passed in via sbatch export.
########################################################
src=$1
INTERMEDIATEDIR=$2
sampleOutdir=$3
sample_name=$4
BAM=$5
BAM_N=$6
BAM_K=$7 # Required. Use "0" if necessary.
BAM_E=$8 # Required. Use "0" if necessary.


jobid=$SLURM_ARRAY_TASK_ID
jobname=`awk -F "\t" -v jobid=$jobid 'NR==jobid {print $1}' ${INTERMEDIATEDIR}/inputs.sort.bam${BAM_N}.txt`
chroms=`awk -F "\t" -v jobid=$jobid 'NR==jobid {print $2}' ${INTERMEDIATEDIR}/inputs.sort.bam${BAM_N}.txt`

BAM_OUT="${INTERMEDIATEDIR}/sorted_bams/${sample_name}.${jobname}.bam"

echo "Sorting ${chroms} of ${BAM}"
echo "Running on $HOSTNAME. Using $TMPDIR as tmp"
date
echo


#Likely more IO efficient to cache the subsetted BAM locally
samtools view -@ ${NSLOTS} -f ${BAM_K} -F ${BAM_E} -O BAM ${BAM} ${chroms} > $TMPDIR/${sample_name}.${jobname}.bam

#Identify read names overlapping the uninformativeRegionFile for removal
#We could implement bed input directly in subsetBAM.py to avoid a second pass, but I think we would need to require the file be sorted
samtools view -L ${INTERMEDIATEDIR}/uninformativeRegionFile.bed $TMPDIR/${sample_name}.${jobname}.bam -O SAM | cut -f1 | uniq > $TMPDIR/uninformativeReads.txt

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
${src}/subsetBAM.py --exclude_readnames $TMPDIR/uninformativeReads.txt $TMPDIR/${sample_name}.${jobname}.bam - |
samtools sort -@ $NSLOTS -O bam -m 4000M -T $TMPDIR/sortbyname -l 8 -n -o ${BAM_OUT}

echo "Done!!!"
date

