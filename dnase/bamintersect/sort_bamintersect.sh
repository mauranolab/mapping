#!/bin/bash
set -eu -o pipefail

########################################################
## Variables passed in via sbatch export.
########################################################

## The "|| true" prevents the SIGPIPE signal problem. It's only needed when set -eo pipefail is enabled.
chrom=$(tail -n+${SLURM_ARRAY_TASK_ID} "${sampleOutdir}/log/chrom_list${BAM_N}_simple" | head -n 1) || true

if [ ${chrom} = "all_other" ]; then
    input_to_samtools=${input_to_samtools2}
else
    input_to_samtools=${chrom}
fi

echo Starting to sort chrom: ${chrom}

x=${BAM%.bam}         # Cut off the trailing "bam"
y=${x##*.}            # Assign the base to y. This will be one of genotypes, like "hg38_full"
z=${x%.${y}}          # Cut off the trailing ${y}
sample_name=${z##*/}  # The remaining base is the sample name.

BAM_OUT="${TMPDIR_bams}/${sample_name}.${chrom}.${BAM_N}.bam"

## Replace pipes with spaces.
input_to_samtools3=${input_to_samtools//|/ }

## Sort bam file by read name.  picard does it lexigraphically.  samtools does it "naturally".
## See: https://github.com/samtools/hts-specs/pull/361 for clarification and pointers to background about this.
samtools view -h -b -f ${BAM_K} -F ${BAM_E} ${BAM} ${input_to_samtools3} | \
java -XX:ParallelGCThreads=2 -Dpicard.useLegacyParser=false -jar $PICARDPATH/picard.jar SortSam \
     -I=/dev/stdin \
     -O=${BAM_OUT} \
     -VERBOSITY=WARNING \
     -QUIET=true \
     -SORT_ORDER=queryname

echo Done with ${chrom}

