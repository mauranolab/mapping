#!/bin/bash
set -eu -o pipefail

#SBATCH --job-name=bam_intersect_array

########################################################
## FINAL_OUTDIR is passed in via an sbatch export
## src is passed in via an sbatch export
## reads_match is passed in via an sbatch export
## make_csv is passed in via an sbatch export
########################################################

echo "SLURM_ARRAY_JOB_ID=$SLURM_ARRAY_JOB_ID"
echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
echo "FINAL_OUTDIR is: ${FINAL_OUTDIR}"   # Passed in via sbatch export
echo

## The "|| true" prevents the SIGPIPE signal problem. It's only needed when set -eo pipefail is enabled.
bam_intersect_data=$(tail -n+${SLURM_ARRAY_TASK_ID} "${FINAL_OUTDIR}/bams/array_list" | head -n 1) || true

read -r bam1 bam2 outdir <<< ${bam_intersect_data}
echo ${bam1}
echo ${bam2}
echo ${outdir}

mkdir ${outdir}    ## "${TEMP_DIR_CL1}/${BASE1}___${BASE2}" 

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

"${src}/bam_intersect.py" --bam1 ${bam1} --bam2 ${bam2} --outdir ${outdir} ${same} ${make_csv_flag}

## Delete directories with zero results from bam_intersect.py
num_chars=$(wc -c < "${outdir}/dsgrep_out.csv")
echo "num_chars -c is:  ${num_chars}"
if [ ${num_chars} -lt 1 ]; then
    rm -rf "${outdir}"
else
    ## Sort to prepare for consolidation later.
    sort-bed "${outdir}/dsgrep_out.csv" > "${outdir}/f1"
    mv "${outdir}/f1"  "${outdir}/dsgrep_out.csv"
fi

