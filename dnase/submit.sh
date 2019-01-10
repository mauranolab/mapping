#!/bin/bash
set -eu -o pipefail


#Load modules
module load picard/2.18.15
module load FastQC/0.11.4
module load bedtools/2.25
module load bedops/2.4.35
module load bwa/0.7.15
module load htslib/1.9
module load samtools/1.9
module load bcftools/1.9
module load trimmomatic/0.38
module load python/3.5.0
module load samblaster/0.1.22
module load ucsckentutils/12152017
module load hotspot/4.1
module load hotspot2/2.1.1
module load pigz


###Hardcoded configuration options
#Common
qsubargs=""
mapThreads=4
mergeThreads=3

#For big jobs:
#qsubargs="--qos normal -p -500"
#Seem to still have some memory problems with mapThreads <3?
#mapThreads=3
#mergeThreads=1


genomesToMap=$1
analysisType=$2
sample=$3
BS=$4
#Files for a given sample will be output to a folder named "${sample}-${BS}-$genome"


processingCommand=`echo "${analysisType}" | awk -F "," '{print $1}'`
analysisCommand=`echo "${analysisType}" | awk -F "," '{print $2}'`

if [[ "${processingCommand}" != "none" ]] && [[ "${processingCommand}" != "aggregate" ]] && [[ "${processingCommand}" != "aggregateRemarkDups" ]] && [[ "${processingCommand}" != "mapBwaAln" ]] && [[ "${processingCommand}" != "mapBwaMem" ]]; then
    echo "ERROR submit: unknown processing command ${processingCommand} in analysisType ${analysisType}"
    exit 1
fi

if [[ "${analysisCommand}" != "atac" ]] && [[ "${analysisCommand}" != "dnase" ]] && [[ "${analysisCommand}" != "chipseq" ]] && [[ "${analysisCommand}" != "callsnps" ]] && [[ "${analysisCommand}" != "none" ]]; then 
    echo "ERROR submit: unknown analysis command ${analysisCommand} in analysisType ${analysisType}"
    exit 2
fi

if ! grep -q ${BS} inputs.txt; then
    echo "submit: Can't find ${BS}"
    exit 3
fi


echo "Will process ${sample} (${BS}) using ${analysisType} pipeline for genomes ${genomesToMap}"


sampleOutdir="${sample}-${BS}"
mkdir -p ${sampleOutdir}

#Run the job from a local copy of the pipeline for archival and to prevent contention issues if the main tree is modified mid-job
src=`pwd`/${sampleOutdir}/.src
mkdir -p ${src}
#use cp -p to preserve timestamps
find /vol/mauranolab/mapped/src -type f | xargs -I {} cp -p {} ${src}


if [[ "${analysisType}" =~ ^map ]]; then
    grep ${BS} inputs.txt > ${sampleOutdir}/inputs.map.txt
    n=`cat ${sampleOutdir}/inputs.map.txt | wc -l`
    mapname="${sample}-${BS}.${genomesToMap}"
    #SGE doesn't accept a -t specification with gaps, so we'll start R2 jobs that will die instantly rather than prune here
    echo "Mapping ${n} jobs for ${mapname}"
    #NB running many jobs with threads < 4 can sometimes get memory errors, not quite sure of the culprit
    qsub -S /bin/bash -cwd -V $qsubargs -pe threads ${mapThreads} -terse -j y -b y -t 1-${n} -o ${sampleOutdir} -N map.${mapname} "${src}/map.sh ${genomesToMap} ${analysisType} ${sample}-${BS} ${BS} ${src}" | perl -pe 's/[^\d].+$//g;' > ${sampleOutdir}/sgeid.map.${mapname}
fi


if [[ "${processingCommand}" =~ ^map ]] || [[ "${processingCommand}" =~ ^aggregate ]]; then
    if [[ "${processingCommand}" =~ ^map ]]; then
        mergeHold="-hold_jid `cat ${sampleOutdir}/sgeid.map.${mapname}`"
    else
        mergeHold=""
    fi
    
    for curGenome in `echo ${genomesToMap} | perl -pe 's/,/ /g;'`; do
        mergename="${sample}-${BS}.${curGenome}"
        echo "${mergename} merge"
        #NB compression threads are multithreaded. However, samblaster and index are not; former is ~1/4 of time and latter is trivial
        qsub -S /bin/bash -cwd -V $qsubargs -pe threads ${mergeThreads} -terse -j y -b y ${mergeHold} -o ${sampleOutdir} -N merge.${mergename} "${src}/merge.sh ${analysisType} ${sample}-${BS} ${BS} ${curGenome}" | perl -pe 's/[^\d].+$//g;' > ${sampleOutdir}/sgeid.merge.${mergename}
    done
fi
#Clean up sgeid even if we don't run merge
if [[ "${processingCommand}" =~ ^map ]]; then
    rm -f ${sampleOutdir}/sgeid.map.${mapname}
fi



for curGenome in `echo ${genomesToMap} | perl -pe 's/,/ /g;'`; do
    analysisname="${sample}-${BS}.${curGenome}"
    
    #Submit job even if none so that the read count stats get run
    if [[ "1" == "1" ]]; then
        if [[ "${processingCommand}" =~ ^map ]] || [[ "${processingCommand}" =~ ^aggregate ]]; then
            analysisHold="-hold_jid `cat ${sampleOutdir}/sgeid.merge.${analysisname}`"
        else
            analysisHold=""
        fi
        
        echo "${analysisname} analysis"
        qsub -S /bin/bash -cwd -V $qsubargs -terse -j y -b y ${analysisHold} -o ${sampleOutdir} -N analysis.${analysisname} "${src}/analysis.sh ${curGenome} ${analysisType} ${sample}-${BS} ${BS} ${src}" | perl -pe 's/[^\d].+$//g;' > ${sampleOutdir}/sgeid.analysis.${analysisname}
        
        cat ${sampleOutdir}/sgeid.analysis.${analysisname} >> sgeid.analysis
        
        rm -f ${sampleOutdir}/sgeid.analysis.${analysisname}
    fi
    
    #Clean up sgeid even if we don't run analysis
    rm -f ${sampleOutdir}/sgeid.merge.${analysisname}
done


echo
