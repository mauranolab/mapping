#!/bin/bash
set -eu -o pipefail

########################################################
## src is passed in via an sbatch export
## reads_match is passed in via an sbatch export
## make_csv is passed in via an sbatch export
## INTERMEDIATEDIR is passed in via an sbatch export ($INTERMEDIATEDIR/sorted_bams contains output of sort_bamintersect.sh)
########################################################

echo "SLURM_ARRAY_JOB_ID=$SLURM_ARRAY_JOB_ID"
echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
echo

## The "|| true" prevents the SIGPIPE signal problem. It's only needed when set -eo pipefail is enabled.
bam_intersect_data=$(tail -n+${SLURM_ARRAY_TASK_ID} "${INTERMEDIATEDIR}/sorted_bams/inputs.array_list.txt" | head -n 1) || true

read -r bam1 bam2 outdir <<< ${bam_intersect_data}
echo ${bam1}
echo ${bam2}
echo ${outdir}      # "${INTERMEDIATEDIR}/bamintersectPyOut/${BASE1}___${BASE2}"

####################################################################################################################
## reads_match=True means that we are looking to pair up read1's with read1's, and read2's with read2's.
##             True is also used for unpaired reads.
## reads_match=False means that we are looking to pair up read1's with read2's.

## Set to True in the input file.
if [ ${reads_match} = "False" ]; then
    same=""
else
    same="--same"
fi

if [ ${make_csv} = "False" ]; then
    make_csv_flag=""
else
    make_csv_flag="--make_csv"
fi

"${src}/bamintersect.py" --bam1 ${bam1} --bam2 ${bam2} --outdir ${TMPDIR} ${same} ${make_csv_flag}

## See if there is anything in dsgrep_out.csv, which is produced by bamintersect.py
num_chars=$(wc -c < "${TMPDIR}/dsgrep_out.csv")
echo "num_chars -c is:  ${num_chars}"

## Let TMPDIR die if dsgrep_out.csv contains no output.
if [ ${num_chars} -ge 1 ]; then
    mv "${TMPDIR}" "${outdir}"

    ## Sort to prepare for consolidation later.
    sort-bed "${outdir}/dsgrep_out.csv" > "${outdir}/f1"
    mv "${outdir}/f1" "${outdir}/dsgrep_out.csv"
fi

