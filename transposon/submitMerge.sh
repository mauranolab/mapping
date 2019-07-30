#!/bin/bash
set -euo pipefail

#TODO merge with submit.sh a la dnase pipeline -- main issue is whether to remove files or not


#Limit thread usage by python processes using OPENBLAS (esp. scipy). Set here and will be inherited by spawned jobs
#https://stackoverflow.com/questions/51256738/multiple-instances-of-python-running-simultaneously-limited-to-35
export OPENBLAS_NUM_THREADS=1

module load pigz
module load samtools/1.9

src=/vol/mauranolab/mapped/src/transposon

#Hardcoded right now rather than as parameters like dnase pipeline
runMerge=1

###Parse command line args
sample=$1
shift
indivsamples=$@


OUTDIR=${sample}
mkdir -p ${OUTDIR}


#export NSLOTS=1


if [[ ${indivsamples} =~ _iPCR ]]; then
    sampleType="iPCR"
    minReadCutoff=2
    bamfiles=`echo ${indivsamples} | perl -pe 's/\/+$//g;' -e 's/\/+( +)/\1/g;' | awk 'BEGIN {ORS=" "} {for(i=1; i<=NF; i++) {split($i, path, "/"); print $i "/" path[length(path)] ".bam" }}'`
    echo -e "Will merge bam files: ${bamfiles}\n"
else
    #As long as it's not iPCR
    sampleType="DNA"
    minReadCutoff=10
    bamfiles=""
fi


if [ ${runMerge} -eq 1 ]; then
    bcfiles=`echo ${indivsamples} | perl -pe 's/\/+$//g;' -e 's/\/+( +)/\1/g;' | awk 'BEGIN {ORS=" "} {for(i=1; i<=NF; i++) {split($i, path, "/"); print $i "/" path[length(path)] ".barcodes.preFilter.txt.gz" }}'`
    echo -e "Will merge barcode files: ${bcfiles}\n"
    
    cat <<EOF | qsub -S /bin/bash -j y -b y -N ${sample} -o ${sample} -terse > sgeid.merge.${sample}
    set -eu -o pipefail
    echo "Merging barcodes"
    zcat -f ${bcfiles} | pigz -p ${NSLOTS} -c -9 > ${OUTDIR}/${sample}.barcodes.preFilter.txt.gz
    
    if [[ ${sampleType} == "iPCR" ]]; then
        echo "Merging bam files"
        if [[ `echo ${bamfiles} | wc | awk '{print $2}'` -gt 1 ]]; then 
            samtools merge -f -l 9 $OUTDIR/${sample}.bam ${bamfiles}
        else 
            cp ${bamfiles} $OUTDIR/${sample}.bam
        fi
        #TODO sort in mapIntegrations.sh?
        #samtools index $OUTDIR/${sample}.bam
    fi
EOF
fi

if [ ${runMerge} -eq 1 ]; then
    analysisHold="-hold_jid `cat sgeid.merge.${sample}`"
else
    analysisHold=""
fi


#-o ${sample} breaks Jesper's Flowcell_Info.sh
cat <<EOF | qsub -S /bin/bash -j y -b y -N ${sample} -terse ${analysisHold} | perl -pe 's/[^\d].+$//g;'  #> sgeid.analysis
set -eu -o pipefail
${src}/analyzeBCcounts.sh ${minReadCutoff} ${sample}

if [[ ${sampleType} == "iPCR" ]]; then
    ${src}/analyzeIntegrations.sh ${sample}
fi
EOF
rm -f sgeid.merge.${sample}

