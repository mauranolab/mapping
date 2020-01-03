#!/bin/bash
set -eu -o pipefail

src=$1
INTERMEDIATEDIR=$2
ReqFullyAligned=$3
max_mismatches=$4
make_bed=${5-} #optional
reads_match=${6-} #optional, BUGBUG if --make_bed not specified


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

