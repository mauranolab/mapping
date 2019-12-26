#!/bin/bash
set -eu -o pipefail

########################################################
## src is passed in via an sbatch export
## reads_match is passed in via an sbatch export
## make_bed is passed in via an sbatch export
## INTERMEDIATEDIR is passed in via an sbatch export ($INTERMEDIATEDIR/sorted_bams contains output of sort_bamintersect.sh)
########################################################

echo "Running on $HOSTNAME. Using $TMPDIR as tmp"


## The "|| true" prevents the SIGPIPE signal problem. It's only needed when set -eo pipefail is enabled.
bam_intersect_data=$(tail -n+${SLURM_ARRAY_TASK_ID} "${INTERMEDIATEDIR}/inputs.bamintersect.txt" | head -n 1) || true
read -r bam1 bam2 outdir <<< ${bam_intersect_data}

####################################################################################################################
## reads_match=True means that we are looking to pair up read1's with read1's, and read2's with read2's.
##             True is also used for unpaired reads.
## reads_match=False means that we are looking to pair up read1's with read2's.
${src}/bamintersect.py --src ${src} --bam1 ${bam1} --bam2 ${bam2} --outdir ${TMPDIR} --max_mismatches ${max_mismatches} ${reads_match} ${make_bed} ${ReqFullyAligned}

## See if there is anything in dsgrep_out.csv, which is produced by bamintersect.py
## If not, then let TMPDIR die. This reduces the number of merge operations needed in merge_bamintersect.sh
if [ -s "${TMPDIR}/dsgrep_out.csv" ]; then
    mkdir -p "${outdir}"
    mv "${TMPDIR}/dsgrep_out.csv" ${outdir}
fi

echo
echo "Done!!!"
date

