#!/bin/bash
set -eu -o pipefail

#SBATCH --job-name=big_array_job

num_lines=$(wc -l < "${FINAL_OUTDIR}/bams/array_list")

## Regarding the size of the sbatch array job we are about to launch:
## N1 = (# of bam1 numbered chromosomes) + (4 [for X,Y,M, and "all_other"])
## N2 = (# of bam2 numbered chromosomes) + (4 [for X,Y,M, and "all_other"])
## Number of array jobs = N1 * N2.
## As an example, mm10 vs rn6 generates 552 jobs:  (19 + 3 + 1) * (20 + 3 + 1)
##
## As long as the number of array jobs is less than 5,000 we will pass through the
## below while loop just once. sbatch will launch all (N1 * N2) jobs in a single array call.

array_incr=5000
array_start=1
array_end=${array_incr}

while :
do
    if [ ${array_end} -gt ${num_lines} ]; then
        array_end=${num_lines}
    fi

    export_vars="FINAL_OUTDIR=${FINAL_OUTDIR}"
    export_vars="${export_vars},src=${src}"
    export_vars="${export_vars},reads_match=${reads_match}"
    export_vars="${export_vars},make_csv=${make_csv}"

    ## Max possible value for array end is 100,000
    ## Max possible number of jobs is 50,000
    sbatch --wait --export="ALL,${export_vars}" --array="${array_start}-${array_end}" \
           --output="${FINAL_OUTDIR}/log/bam_intersect_array.o_%A_%a" "${src}/bam_intersect_array.sbatch"

    if [ ${array_end} -eq ${num_lines} ]; then
        break
    else
        let "array_start = ${array_start} + ${array_incr}"
        let "array_end = ${array_end} + ${array_incr}"
    fi
done

