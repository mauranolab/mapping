#!/bin/bash
set -eu -o pipefail

########################################################
## src is passed in via an sbatch export
## reads_match is passed in via an sbatch export
## make_bed is passed in via an sbatch export
## INTERMEDIATEDIR is passed in via an sbatch export ($INTERMEDIATEDIR/sorted_bams contains output of sort_bamintersect.sh)
########################################################

echo "Running on $HOSTNAME. Using $TMPDIR as tmp"

jobid=${SLURM_ARRAY_TASK_ID}
bam1=`awk -F "\t" -v jobid=${jobid} 'NR==jobid {print $1}' ${INTERMEDIATEDIR}/inputs.bamintersect.txt`
bam2=`awk -F "\t" -v jobid=${jobid} 'NR==jobid {print $2}' ${INTERMEDIATEDIR}/inputs.bamintersect.txt`
outfile=`awk -F "\t" -v jobid=${jobid} 'NR==jobid {print $3}' ${INTERMEDIATEDIR}/inputs.bamintersect.txt`


## reads_match=True means that we are looking to pair up read1's with read1's, and read2's with read2's.
##             True is also used for unpaired reads.
## reads_match=False means that we are looking to pair up read1's with read2's.
#BUGBUG doesn't respect ${make_bed}
${src}/bamintersect.py ${bam1} ${bam2} --src ${src} --bedout ${TMPDIR}/bamintersect_out.bed --max_mismatches ${max_mismatches} ${reads_match} ${ReqFullyAligned}
sort-bed ${TMPDIR}/bamintersect_out.bed > ${outfile}

echo
echo "Done!!!"
date

