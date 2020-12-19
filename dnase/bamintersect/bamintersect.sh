#!/bin/bash
set -eu -o pipefail

max_mismatches=${1}
INTERMEDIATEDIR=${2}
src=${3}
reads_match=${4-} #optional


echo "Running on $HOSTNAME. Using $TMPDIR as tmp"

jobid=${SLURM_ARRAY_TASK_ID}
bam1=`awk -F "\t" -v jobid=${jobid} 'NR==jobid {print $1}' ${INTERMEDIATEDIR}/inputs.bamintersect.txt`
bam2=`awk -F "\t" -v jobid=${jobid} 'NR==jobid {print $2}' ${INTERMEDIATEDIR}/inputs.bamintersect.txt`
outfile=`awk -F "\t" -v jobid=${jobid} 'NR==jobid {print $3}' ${INTERMEDIATEDIR}/inputs.bamintersect.txt`

echo "bam1: ${bam1}"
echo "bam2: ${bam2}"
echo "outfile: ${outfile}"

## reads_match=True means that we are looking to pair up read1's with read1's, and read2's with read2's.
##             True is also used for unpaired reads.
## reads_match=False means that we are looking to pair up read1's with read2's.
${src}/bamintersect.py ${bam1} ${bam2} --src ${src} --bedout ${TMPDIR}/bamintersect_out.bed --max_mismatches ${max_mismatches} ${reads_match} --ReqFullyAligned --bam1out ${TMPDIR}/bamintersect_out.bam1.bam --bam2out ${TMPDIR}/bamintersect_out.bam2.bam

echo
echo "Sorting bed file"
sort-bed ${TMPDIR}/bamintersect_out.bed > ${outfile}

echo
echo "Sorting bam file"
samtools sort -@ $NSLOTS -O bam -m 4000M -T $TMPDIR/sam.sort -l 1 -o `echo ${outfile} | perl -pe 's/\.bed$/.bam1.bam/g;'` ${TMPDIR}/bamintersect_out.bam1.bam
samtools sort -@ $NSLOTS -O bam -m 4000M -T $TMPDIR/sam.sort -l 1 -o `echo ${outfile} | perl -pe 's/\.bed$/.bam2.bam/g;'` ${TMPDIR}/bamintersect_out.bam2.bam


echo
echo "Done!!!"
date

