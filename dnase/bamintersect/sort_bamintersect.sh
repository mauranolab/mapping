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

## Sort bam file by read name.  samtools does it "naturally", via strnum_cmp.c
samtools view -h -b -f ${BAM_K} -F ${BAM_E} ${BAM} ${input_to_samtools3} | samtools sort -n -o ${BAM_OUT}

echo "Done!!!"
date

