#!/bin/bash
set -eu -o pipefail

#NB submit.sh must be called with a full path so we can infer srcbase

#Limit thread usage by python processes using OPENBLAS (esp. scipy). Set here and will be inherited by spawned jobs
#https://stackoverflow.com/questions/51256738/multiple-instances-of-python-running-simultaneously-limited-to-35
export OPENBLAS_NUM_THREADS=1


#Load modules
module load picard/2.22.1
module load FastQC/0.11.4
module load bedtools/2.29.2
module load bedops/2.4.39
module load bwa/0.7.17
module load minimap2/2.17
module load htslib/1.10.2
module load samtools/1.10
module load bcftools/1.10.2
module load trimmomatic/0.39
module load python/3.8.1
module load samblaster/0.1.24
module load ucsckentutils/379
module load hotspot/4.1
module load hotspot2/2.1.1
module load pigz
module load primerclip/0.3.8


###Hardcoded configuration options
#Common
qsubargs=""
mapThreads=4
mergeThreads=3
runBamIntersect=1

#For big jobs:
#qsubargs="--qos normal -p -500"
#Seem to still have some memory problems with mapThreads <3? Maybe analysis.sh has a transient memory peak, say in hotspot?
#mapThreads=3
#mergeThreads=1


genomesToMap=$1
analysisType=$2
samplePrefix=$3
BS=$4
if [ "$#" -ge 5 ]; then
    sampleAnnotation=$5
else
    #Dummy value to preserve argument position
    sampleAnnotation="NULL"
fi
#Files for a given sample will be output to a folder named "${samplePrefix}-${BS}" (samplePrefix may include subsequent subdirectories)


processingCommand=`echo "${analysisType}" | awk -F "," '{print $1}'`
sampleType=`echo "${analysisType}" | awk -F "," '{print $2}'`

if [[ "${processingCommand}" != "none" ]] && [[ "${processingCommand}" != "aggregate" ]] && [[ "${processingCommand}" != "aggregateRemarkDups" ]] && [[ "${processingCommand}" != "mapBwaAln" ]] && [[ "${processingCommand}" != "mapBwaMem" ]] && [[ "${processingCommand}" != "mapMinimap" ]] &&  [[ "${processingCommand}" != "bamintersect" ]]; then
    echo "ERROR submit: unknown processing command ${processingCommand} in analysisType ${analysisType}"
    exit 1
fi

if [[ "${sampleType}" != "atac" ]] && [[ "${sampleType}" != "dnase" ]] && [[ "${sampleType}" != "chipseq" ]] && [[ "${sampleType}" != "dna" ]] && [[ "${sampleType}" != "capture" ]] && [[ "${sampleType}" != "amplicon" ]] && [[ "${sampleType}" != "none" ]]; then
    echo "ERROR submit: unknown sample type ${sampleType} in analysisType ${analysisType}"
    exit 2
fi

if ! grep -q ${BS} inputs.txt; then
    echo "ERROR submit: Can't find ${BS}"
    exit 3
fi


echo "Submitting jobs for ${samplePrefix} (${BS})"
echo "Pipeline: ${analysisType}"
echo "Genomes ${genomesToMap}"

sampleOutdir="${samplePrefix}-${BS}"
name=`basename ${sampleOutdir}`
mkdir -p ${sampleOutdir}

#Run the job from a local copy of the pipeline for archival and to prevent contention issues if the main tree is modified mid-job
srcbase=$( dirname "${BASH_SOURCE[0]}" )
src=`pwd`/${sampleOutdir}/.src
mkdir -p ${src}
#Avoid slowing things down by copying just the files in the base directory (i.e. not the trackhub directory)
#use cp -p to preserve timestamps
find ${srcbase} -maxdepth 1 -type f | xargs -I {} cp -p {} ${src}
#Only subdir we need is bamintersect
cp -rp ${srcbase}/bamintersect ${src}


###Map
if [[ "${processingCommand}" =~ ^map ]]; then
    grep ${BS} inputs.txt > ${sampleOutdir}/inputs.map.txt
    n=`cat ${sampleOutdir}/inputs.map.txt | wc -l`
    mapname="${name}.${genomesToMap}"
    #SGE doesn't accept a -t specification with gaps, so we'll start R2 jobs that will die instantly rather than prune here
    echo
    echo "Mapping ${n} jobs"
    echo "+ ${mapname}"
    #NB running many jobs with threads < 4 can sometimes get memory errors, not quite sure of the culprit
    qsub -S /bin/bash -cwd -V ${qsubargs} -pe threads ${mapThreads} -terse -j y -b y -t 1-${n} -o ${sampleOutdir} -N map.${mapname} "${src}/map.sh ${genomesToMap} ${analysisType} ${sampleOutdir} ${BS} ${src}" | perl -pe 's/[^\d].+$//g;' > ${sampleOutdir}/sgeid.map.${mapname}
fi


###Merge
if [[ "${processingCommand}" =~ ^map ]] || [[ "${processingCommand}" =~ ^aggregate ]]; then
    if [[ "${processingCommand}" =~ ^map ]]; then
        mergeHold="-hold_jid `cat ${sampleOutdir}/sgeid.map.${mapname}`"
    else
        mergeHold=""
    fi
    
    echo
    echo "Merge"
    for curGenome in `echo ${genomesToMap} | perl -pe 's/,/ /g;'`; do
        mergename="${name}.${curGenome}"
        echo "+ ${mergename}"
        #NB compression threads are multithreaded. However, samblaster and index are not; former is ~1/4 of time and latter is trivial
        qsub -S /bin/bash -cwd -V ${qsubargs} -pe threads ${mergeThreads} -terse -j y -b y ${mergeHold} -o ${sampleOutdir} -N merge.${mergename} "${src}/merge.sh ${analysisType} ${sampleOutdir} ${BS} ${curGenome} ${src}" | perl -pe 's/[^\d].+$//g;' > ${sampleOutdir}/sgeid.merge.${mergename}
    done
fi
#Clean up sgeid even if we don't run merge
if [[ "${processingCommand}" =~ ^map ]]; then
    rm -f ${sampleOutdir}/sgeid.map.${mapname}
fi


###bamintersect
#Parses Genetic Modification entry from LIMS to return the integration site (in square brackets) for the given entry
function getIntegrationSite {
    local sampleAnnotationGeneticModification=$1
    local entry=$2
    
    #Parse out integration site. If we can't find it for whatever reason, use null
    #Parse each genetic modification on separate line, and move part in [] to $2
    echo "${sampleAnnotationGeneticModification}" | perl -pe 's/,/\n/g;' | perl -pe 's/\[/\t/g;' -e 's/\]/\t/g;' | awk -v entry=${entry} -F "\t" 'BEGIN {integrationsite="null"} $1==entry {integrationsite=$2} END {print integrationsite}'
}


if [ "${runBamIntersect}" -eq 1 ] && ([[ "${processingCommand}" == "bamintersect" ]] || [[ "${sampleType}" == "dna" ]] || [[ "${sampleType}" == "capture" ]]); then
    echo
    echo "bamintersect"
    
    sampleAnnotationGeneticModification=`echo "${sampleAnnotation}" | awk -v key="Genetic_Modification" -F ";" '{for(i=1; i<=NF; i++) { split($i, cur, "="); if(cur[1]==key) {print cur[2]; exit}}}'`
    #TODO should this also include Custom Reference? That is missing an integration site
    
    #Put LPICE into both so we get LPICE vs. payload.
    #Manually skip rn6 since we haven't integrated anything there
    mammalianGenomes=`echo "${genomesToMap}" | perl -pe 's/,/\n/g;' | awk '$0!~/^cegsvectors_/ && $1!="rn6" || $0=="cegsvectors_LPICE"'`
    cegsGenomes=`echo "${genomesToMap}" | perl -pe 's/,/\n/g;' | awk '$0~/^cegsvectors_/'`
    for mammalianGenome in ${mammalianGenomes}; do
        #Same as genomeinfo.sh
        mammalianAnnotationGenome=`echo ${mammalianGenome} | perl -pe's/^cegsvectors_LPICE$/LPICE/g;' -e 's/_.+$//g;' -e 's/all$//g;'`
        for cegsGenome in ${cegsGenomes}; do
            cegsGenomeShort="${cegsGenome/cegsvectors_/}"
            #Limit this to references listed as genetic modifications
            #If we are doing ICE, make sure we don't do LPICE vs LPICE
            #Need to do LPICE vs. assemblon but manually skip LPICE vs. pSpCas9 or rtTA
            if [[ "${sampleAnnotationGeneticModification}" =~ ${cegsGenomeShort} ]] && [[ "${mammalianGenome}" != "${cegsGenome}" ]] && ( [[ "${mammalianGenome}" != "cegsvectors_LPICE" ]] || ([[ ! "${cegsGenome}" =~ "cegsvectors_pSpCas9" ]] && [[ "${cegsGenome}" != "cegsvectors_rtTA" ]]) ); then
                integrationsite=`getIntegrationSite ${sampleAnnotationGeneticModification} ${cegsGenomeShort}`
                #For assemblons this yields the LP name rather looking up the LP integrationsite, so look up the site for the LP itself
                if [[ "${integrationsite}" =~ ^LP ]]; then
                    #echo "found ${integrationsite}, looking up again"
                    integrationsite=`getIntegrationSite ${sampleAnnotationGeneticModification} ${integrationsite}`
                    if [[ "${integrationsite}" == "null" ]]; then
                        echo "WARNING could not find LP integration site for ${cegsGenomeShort}"
                    fi
                fi
                echo "+ ${mammalianAnnotationGenome} vs. ${cegsGenomeShort}[${integrationsite}]"
                
                if [[ "${processingCommand}" =~ ^map ]] || [[ "${processingCommand}" =~ ^aggregate ]]; then
                    bamIntersectHold="-hold_jid `cat ${sampleOutdir}/sgeid.merge.${name}.${mammalianGenome} ${sampleOutdir}/sgeid.merge.${name}.${cegsGenome} | perl -pe 's/\n/,/g;'`"
                else
                    bamIntersectHold=""
                fi
                
                mkdir -p ${sampleOutdir}/bamintersect/log
                qsub -S /bin/bash -cwd -V ${qsubargs} -terse -j y -b y ${bamIntersectHold} -o ${sampleOutdir}/bamintersect/log -N submit_bamintersect.${name}.${mammalianAnnotationGenome}_vs_${cegsGenomeShort} "${src}/bamintersect/submit_bamintersect.sh --sample_name ${name} --outdir ${sampleOutdir}/bamintersect --bam1 ${sampleOutdir}/${name}.${mammalianGenome}.bam --bam1genome ${mammalianAnnotationGenome} --bam2 ${sampleOutdir}/${name}.${cegsGenome}.bam --bam2genome ${cegsGenomeShort} --integrationsite ${integrationsite}" > /dev/null
            fi
        done
    done
fi


###Analysis
#Always run analysis job, but specifying none will have analysis.sh skip most analyses
if [[ "${processingCommand}" != bamintersect ]]; then
    echo
    echo "Analysis"
fi

for curGenome in `echo ${genomesToMap} | perl -pe 's/,/ /g;'`; do
    analysisname="${name}.${curGenome}"
    
    #Always submit analysis job, unless we are just running bamintersect. If processingCommand is none, analysis.sh will just run the read count stats and nothing else
    if [[ "${processingCommand}" != bamintersect ]]; then
        if [[ "${processingCommand}" =~ ^map ]] || [[ "${processingCommand}" =~ ^aggregate ]]; then
            analysisHold="-hold_jid `cat ${sampleOutdir}/sgeid.merge.${analysisname}`"
        else
            analysisHold=""
        fi
        
        echo "+ ${analysisname}"
        qsub -S /bin/bash -cwd -V ${qsubargs} -terse -j y -b y ${analysisHold} -o ${sampleOutdir} -N analysis.${analysisname} "${src}/analysis.sh ${curGenome} ${analysisType} ${sampleOutdir} ${BS} \"${sampleAnnotation}\" ${src}" | perl -pe 's/[^\d].+$//g;' > ${sampleOutdir}/sgeid.analysis.${analysisname}
        
        cat ${sampleOutdir}/sgeid.analysis.${analysisname} >> `dirname ${sampleOutdir}`/sgeid.analysis
        
        rm -f ${sampleOutdir}/sgeid.analysis.${analysisname}
    fi
    
    #Clean up sgeid even if we don't run analysis
    rm -f ${sampleOutdir}/sgeid.merge.${analysisname}
done


echo
echo
