#!/bin/bash
set -eu -o pipefail

########################################################
## Variables passed in via sbatch export.
########################################################

## The "|| true" prevents the SIGPIPE signal problem. It's only needed when set -eo pipefail is enabled.
jobid=$SLURM_ARRAY_TASK_ID
chrom=`awk -v jobid=$jobid 'NR==jobid' ${sampleOutdir}/log/${sample_name}.chrom_list${BAM_N}_simple`

if [ ${chrom} = "all_other" ]; then
    input_to_samtools=${input_to_samtools2}
else
    input_to_samtools=${chrom}
fi

BAM_OUT="${INTERMEDIATEDIR}/sorted_bams/${sample_name}.${chrom}.${BAM_N}.bam"

## Replace pipes with spaces.
input_to_samtools3=${input_to_samtools//|/ }

echo "Beginning sort"
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
#     The charcter comparisons are made by ASCII value.
#     
#     Numeric sorts do not entail character by character comparisons.  The entire subfield is considered as one integer, and
#     the integer values of the two subfields are compared.  Leading zeroes are only used as tie breakers. The subfield with
#     more leading zeroes is placed before the subfield with fewer leading zeroes. 
samtools view -h -b -f ${BAM_K} -F ${BAM_E} ${BAM} ${input_to_samtools3} | samtools sort -n -o ${BAM_OUT}

echo "Done!!!"
date

