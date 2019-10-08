#!/bin/bash
set -eu -o pipefail

########################################################
## Variables passed in via sbatch export.
########################################################

## The "|| true" prevents the SIGPIPE signal problem. It's only needed when set -eo pipefail is enabled.
chrom=$(tail -n+${SLURM_ARRAY_TASK_ID} "${sample_name}/log/${sample_name}.chrom_list${BAM_N}_simple" | head -n 1) || true

if [ ${chrom} = "all_other" ]; then
    input_to_samtools=${input_to_samtools2}
else
    input_to_samtools=${chrom}
fi

echo Starting to sort chrom: ${chrom}

BAM_OUT="${INTERMEDIATEDIR}/sorted_bams/${sample_name}.${chrom}.${BAM_N}.bam"

## Replace pipes with spaces.
input_to_samtools3=${input_to_samtools//|/ }

## Regarding this java error, which is sometimes generated by picard SortSam:
##         There is insufficient memory for the Java Runtime Environment to continue.
##         Cannot create GC thread. Out of system resources.
##         ... etc.
##
## See: https://dzone.com/articles/troubleshoot-outofmemoryerror-unable-to-create-new
##                                      and
##      https://javapapers.com/java/types-of-java-garbage-collectors
## 
## I tried various combinations of -XX:ParallelGCThreads=<#> -XX:ConcGCThreads=<#> -Xms<#>m -Xmx<#>m -Xss<#>k
##
## This:  -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 and also this: -XX:+UseSerialGC seem to work equally well.
##
## The various memory adjustments can allow for a small amount of extra SortSam jobs, but does not seem to be
## material in solving the thread problem, which occurs when other jobs are running on the cluster.

## Sort bam file by read name.  picard does it lexigraphically.  samtools does it "naturally".
## See: https://github.com/samtools/hts-specs/pull/361 for clarification and pointers to background about this.
samtools view -h -b -f ${BAM_K} -F ${BAM_E} ${BAM} ${input_to_samtools3} | \
java -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -Dpicard.useLegacyParser=false -jar $PICARDPATH/picard.jar SortSam \
     -I=/dev/stdin \
     -O=${BAM_OUT} \
     -VERBOSITY=WARNING \
     -QUIET=true \
     -SORT_ORDER=queryname

echo "Done with ${chrom}"

