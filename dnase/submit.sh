#!/bin/bash
set -eu -o pipefail


#Load modules
module load picard/2.18.15
module load FastQC/0.11.4
module load bedops/2.4.35
module load bwa/0.7.15
module load htslib/1.9
module load samtools/1.9
module load bcftools/1.9
module load trimmomatic/0.38
module load python/3.5.0
module load samblaster/0.1.22
module load ucsckentutils/10132015
module load hotspot/4.1
module load hotspot2/2.1.1


###Hardcoded configuration options
#Common
qsubargs=""


genomesToMap=$1
analysisType=$2
sample=$3
DS=$4
#Files for a given sample will be output to a folder named "${sample}-${DS}-$genome"

src=/vol/mauranolab/mapped/src/


processingCommand=`echo "${analysisType}" | awk -F "," '{print $1}'`
analysisCommand=`echo "${analysisType}" | awk -F "," '{print $2}'`

if [[ "${processingCommand}" != "none" ]] && [[ "${processingCommand}" != "aggregate" ]] && [[ "${processingCommand}" != "aggregateRemarkDups" ]] && [[ "${processingCommand}" != "mapBwaAln" ]] && [[ "${processingCommand}" != "mapBwaMem" ]]; then
    echo "ERROR submit: unknown processing command ${processingCommand} in analysisType ${analysisType}"
    exit 1
fi

if [[ "${analysisCommand}" != "atac" ]] && [[ "${analysisCommand}" != "dnase" ]] && [[ "${analysisCommand}" != "callsnps" ]] && [[ "${analysisCommand}" != "none" ]]; then 
    echo "ERROR submit: unknown analysis command ${analysisCommand} in analysisType ${analysisType}"
    exit 2
fi

if ! grep -q ${DS} inputs.txt; then
    echo "Can't find ${DS}"
    exit 3
fi


echo "Will process ${sample} (${DS}) using ${analysisType} pipeline for genomes ${genomesToMap}"


sampleOutdir="${sample}-${DS}"
mkdir -p ${sampleOutdir}


if [[ "${analysisType}" =~ ^map ]]; then
    numlines=`grep ${DS} inputs.txt | wc -l`
    
    firstline=`grep -n ${DS} inputs.txt | head -1 | awk -F ":" '{print $1}'`
    lastline=`grep -n ${DS} inputs.txt | tail -1 | awk -F ":" '{print $1}'`
    expectednumlines=`echo "$lastline-$firstline+1" | bc -l -q`
    if [ "$expectednumlines" -ne "$numlines" ]; then
        echo "Wrong number of lines (expected $expectednumlines, found $numlines), are all fastq files contiguous in inputs.txt?"
        exit 4
    fi
    
    mapname="${sample}-${DS}.${genomesToMap}"
    #SGE doesn't accept a -t specification with gaps, so we'll start R2 jobs that will die instantly rather than prune here
    echo "Mapping ${mapname} input.txt lines $firstline-$lastline"
    qsub -S /bin/bash -cwd -V $qsubargs -pe threads 8 -terse -j y -b y -t $firstline-$lastline -o ${sampleOutdir} -N map.${mapname} "${src}/map.sh ${genomesToMap} ${analysisType} ${sample}-${DS} ${DS} ${src}" | perl -pe 's/[^\d].+$//g;' > ${sampleOutdir}/sgeid.map.${mapname}
fi


if [[ "${processingCommand}" =~ ^map ]] || [[ "${processingCommand}" =~ ^aggregate ]]; then
    if [[ "${processingCommand}" =~ ^map ]]; then
        mergeHold="-hold_jid `cat ${sampleOutdir}/sgeid.map.${mapname}`"
    else
        mergeHold=""
    fi
    
    for curGenome in `echo ${genomesToMap} | perl -pe 's/,/ /g;'`; do
        mergename="${sample}-${DS}.${curGenome}"
        echo "${mergename} merge"
        qsub -S /bin/bash -cwd -V $qsubargs -terse -j y -b y ${mergeHold} -o ${sampleOutdir} -N merge.${mergename} "${src}/merge.sh ${analysisType} ${sample}-${DS} ${DS} ${curGenome}" | perl -pe 's/[^\d].+$//g;' > ${sampleOutdir}/sgeid.merge.${mergename}
    done
fi
#Clean up sgeid even if we don't run merge
if [[ "${processingCommand}" =~ ^map ]]; then
    rm -f ${sampleOutdir}/sgeid.map.${mapname}
fi



for curGenome in `echo ${genomesToMap} | perl -pe 's/,/ /g;'`; do
    analysisname="${sample}-${DS}.${curGenome}"
    
    #Submit job even if none so that the read count stats get run
    if [[ "1" == "1" ]]; then
#    if [[ "${analysisCommand}" != "none" ]]; then
        if [[ "${processingCommand}" =~ ^map ]] || [[ "${processingCommand}" =~ ^aggregate ]]; then
            analysisHold="-hold_jid `cat ${sampleOutdir}/sgeid.merge.${analysisname}`"
        else
            analysisHold=""
        fi
        
        echo "${analysisname} analysis"
        qsub -S /bin/bash -cwd -V $qsubargs -terse -j y -b y ${analysisHold} -o ${sampleOutdir} -N analysis.${analysisname} "${src}/analysis.sh ${curGenome} ${analysisType} ${sample}-${DS} ${DS} ${src}" | perl -pe 's/[^\d].+$//g;' > ${sampleOutdir}/sgeid.analysis.${analysisname}
        
        cat ${sampleOutdir}/sgeid.analysis.${analysisname} >> sgeid.analysis
        
        rm -f ${sampleOutdir}/sgeid.analysis.${analysisname}
    fi
    
    #Clean up sgeid even if we don't run analysis
    rm -f ${sampleOutdir}/sgeid.merge.${analysisname}
done


echo
