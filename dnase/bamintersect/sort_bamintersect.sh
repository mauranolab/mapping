#!/bin/bash
set -eu -o pipefail

########################################################
## Variables passed in via sbatch export.
########################################################
INTERMEDIATEDIR=$1
sampleOutdir=$2
sample_name=$3
BAM=$4
BAM_N=$5
BAM_K=$6 # Required. Use "0" if necessary.
BAM_E=$7 # Required. Use "0" if necessary.


jobid=$SLURM_ARRAY_TASK_ID
jobname=`awk -F "\t" -v jobid=$jobid 'NR==jobid {print $1}' ${INTERMEDIATEDIR}/inputs.sort.bam${BAM_N}.txt`
chroms=`awk -F "\t" -v jobid=$jobid 'NR==jobid {print $2}' ${INTERMEDIATEDIR}/inputs.sort.bam${BAM_N}.txt`

BAM_OUT="${INTERMEDIATEDIR}/sorted_bams/${sample_name}.${jobname}.bam"

echo "Sorting ${chroms}"
echo "Running on $HOSTNAME. Using $TMPDIR as tmp"
date
echo

## Sort bam file by read name.  Done by samtools via strnum_cmp.c
#     Query names are split into alternating subfields of pure nondigits and pure digits.
#     
#     To sort a pair of query names, the corresponding subfields of each query name are compared sequentially:
#         When a nondigit subfield is compared to a nondigit subfield, the sort is lexigraphic.
#         When a all-digit subfield is compared to a all-digit subfield, the sort is numeric.
#         When a all-digit subfield is compared to a nondigit subfield, the sort is lexigraphic.
#       
#     Lexigraphic sorts by left justifying both subfields, and comparing them from left to right, character by character.
#     The character comparisons are made by ASCII value.
#     
#     Numeric sorts do not entail character by character comparisons.  The entire subfield is considered as one integer, and
#     the integer values of the two subfields are compared.  Leading zeroes are only used as tie breakers. The subfield with
#     more leading zeroes is placed before the subfield with fewer leading zeroes. 

#Since the next set of jobs will be reading this heavily, might be worthwhile to use high compression (-l 9) for big bam files
samtools view -h -u -f ${BAM_K} -F ${BAM_E} ${BAM} ${chroms} | samtools sort -@ $NSLOTS -O bam -m 4000M -T $TMPDIR/sortbyname -l 1 -n -o ${BAM_OUT}

echo "Done!!!"
date

