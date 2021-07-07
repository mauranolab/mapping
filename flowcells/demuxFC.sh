#!/bin/bash
#set -eu -o pipefail
set -u

#Used for output directory name
fc=$1
shift
userDemuxOptions=$@

#--use-bases-mask Y*,I*,n*,Y* demuxes on BC1 alone for all samples


module add bcl2fastq2/2.20


echo
echo "Demultiplexing ${fc}; userDemuxOptions=${userDemuxOptions}"
date


if [ -e "../../fastq/${fc}" ]; then
    echo "ERROR: ../../fastq/${fc} already exists!"
    exit 1
fi


#https://hpc.nih.gov/apps/bcl2fastq.html says optimal relative allocations but doesn't seem to change runtime much:- r 0.25 -w 0.25 -p 0.875
#bumps up processing-threads over NSLOTS up based on observed CPU usage
#Info on --use-bases-mask https://www.biostars.org/p/344768/
#TODO we don't split collaborator or undetermined fastqs; otherwise could compress there
bcl2fastq --fastq-compression-level 9 --no-lane-splitting --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0 --barcode-mismatches 0 --loading-threads 3 --writing-threads 3 --processing-threads $((${NSLOTS}/7*8)) --output-dir ../../fastq/${fc} ${userDemuxOptions} > bcl2fastq.log 2>&1


#BUGBUG seems to set nonzero exit code?


echo "Restricting permissions on non-mauranolab/CEGS/SARS projects"
find ../../fastq/${fc}  -maxdepth 1 -name "Project_*" -type d -not -name "Project_Maurano" -not -name "Project_CEGS" -not -name "Project_SARS" -print0 | xargs -0 --no-run-if-empty echo chgrp mauranolab-seq


#Running fastqc
#cd ../../fastq/${fc}
#mkdir -p fastqc
#fc="${fc}"
#for f in $files; do
#    sample1=`basename $f .fastq.gz`
#    echo "$sample1"
#    fastQcOutdir="fastqc/${fc}_${sample1}_qc"
#    if [ ! -d "$fastQcOutdir" ]; then
#        qsub -cwd -V -N ${sample1}.qc -o fastqc/${sample1}.qc -S /bin/bash -j y -b y -p -500 "mkdir -p $fastQcOutdir; fastqc --outdir $fastQcOutdir $f"
#    fi
#done


echo
echo "Done!!!"
date
