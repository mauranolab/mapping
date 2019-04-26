#!/bin/bash
set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'

analysisType=$1
sampleOutdir=$2
BS=$3
mappedgenome=$4

processingCommand=`echo "${analysisType}" | awk -F "," '{print $1}'`
#analysisCommand=`echo "${analysisType}" | awk -F "," '{print $2}'`

name=`basename ${sampleOutdir}`

echo "Running on $HOSTNAME. Using $TMPDIR as tmp"
echo "Merging bams (${processingCommand} analysis); output name ${name} (${BS}) against genome ${mappedgenome}"
date


files=""
numfiles=0
if [[ "${processingCommand}" =~ ^map ]]; then
    #Here we generate a list of expected filenames
    #TODO switch to wildcard matching based on sample, DS and mappedgenome?
    #NB originally checked grep R1 but some old solexa lanes didn't include it
    for f1 in `grep ${BS} inputs.txt | grep -v _R2`; do
        #First get FC
        #NB This must match exactly what is in map.sh
        #BUGBUG a bit fragile
        fc=`readlink -f ${f1} | xargs dirname | xargs dirname | xargs dirname | xargs basename`
        if [[ ! "${fc}" =~ ^FC ]] ; then
            fc=""
        else
            fc="${fc}."
        fi
        
        
        f2=`echo ${f1} | perl -pe 's/_R1(_\d+)?/_R2$1/g;'`
        if echo ${f2} | grep -q R2 && grep -q ${f2} inputs.txt ; then
            f1=`echo ${f1} | perl -pe 's/_R1(_\d+)?/_R1R2\1/g;'`
        fi
        
        curOutputFile=`basename ${f1} .fastq.gz`
        curOutputFile="${sampleOutdir}/${fc}${curOutputFile}.${mappedgenome}.bam"
        
        if [[ -f "${curOutputFile}" ]]; then
            files="${files} ${curOutputFile}"
            numfiles=$((numfiles+1))
        else
            echo "WARNING: ${curOutputFile} doesn't exist"
            #Don't die to work around weirdness in pipeline FQ files
            #exit 1
        fi
    done
elif [[ "${processingCommand}" =~ ^aggregate ]]; then
    for curOutputFile in `grep ${BS} inputs.txt | grep ${mappedgenome}`; do
        if [[ -f "${curOutputFile}" ]]; then
            if [[ "${processingCommand}" == "aggregateRemarkDups" ]]; then
                curOutputFileBase=`basename ${curOutputFile} .bam`
                #Include random component to avoid name collisions; the usual method using tr returns nonzero exit code
                curOutputFileBase="${curOutputFileBase}."`head -200 /dev/urandom | cksum | cut -f1 -d " "`
                echo "sorting ${curOutputFile} by read name"
                samtools sort -@ $NSLOTS -O bam -m 5000M -T $TMPDIR/${curOutputFileBase}.sortbyname -l 1 -n ${curOutputFile} > ${TMPDIR}/${curOutputFileBase}.sorted.bam
                #OK to list temporary files since bams don't get deleted at the end of pipeline
                curOutputFile="${TMPDIR}/${curOutputFileBase}.sorted.bam"
            fi
            
            files="${files} ${curOutputFile}"
            numfiles=$((numfiles+1))
        else
            echo "ERROR: ${curOutputFile} doesn't exist!"
            exit 1
        fi
    done
else
    echo "ERROR: Unknown command processingCommand!"
    exit 1
fi


if [[ "${numfiles}" -eq 0 ]]; then
    echo "ERROR: No files found to merge!"
    exit 2
fi

echo "Will merge ${files}"


if [[ "${numfiles}" -eq 1 ]]; then
    echo "copying file"
    cp ${files} ${sampleOutdir}/${name}.${mappedgenome}.bam
else
    echo "merging files"
    
    if [[ "${processingCommand}" =~ ^map ]] || [[ "${processingCommand}" == "aggregateRemarkDups" ]]; then
        #mapping pipeline produces bams sorted by read name, or we sorted the individual bams to merge by name above. These will get sorted by coord after marking duplicates
        mergeOpts="-n -l 1"
    else
        mergeOpts="-l 9"
    fi
    
    #TODO this has massive memory usage for large runs (e.g. a full mouse genome uses 90% of 256 GB)
    #would -p be worth it for map runs to clean up the superfluous @PG tags in header?
    samtools merge -c -@ $NSLOTS ${mergeOpts} ${sampleOutdir}/${name}.${mappedgenome}.bam ${files}
fi

if [[ "${processingCommand}" =~ ^map ]] || [[ "${processingCommand}" == "aggregateRemarkDups" ]]; then
    echo
    echo "Marking duplicates"
    date
    ###Obsolete picard version
    #Used to need VALIDATION_STRINGENCY=LENIENT to avoid SAM validation error: ERROR...MAPQ should be 0 for unmapped read or CIGAR should have zero elements for unmapped read
    #http://seqanswers.com/forums/showthread.php?t=4246
    #BUGBUG Can make huge log files despite these options
    #http://sourceforge.net/p/samtools/mailman/message/32910359/
    #java -Xmx6g -jar /home/maurano/bin/picard-tools/MarkDuplicates.jar INPUT=${name}.${mappedgenome}.bam OUTPUT=${name}.markedDups.bam METRICS_FILE=$TMPDIR/${name}.picardDups.txt QUIET=TRUE VERBOSITY=ERROR COMPRESSION_LEVEL=9 ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT && mv ${name}.markedDups.bam ${name}.${mappedgenome}.bam

    ###Samblaster is faster
    #samblaster used an average of 1GB memory mapping ENCODE DNase data to hg38. 10/889 jobs used >5GB.
    samtools view -h ${sampleOutdir}/${name}.${mappedgenome}.bam |
    samblaster --addMateTags |
    #BUGBUG samblaster does not add a PP field in its @PG tag; when we are merging after map PP should just be bwa, but would be the prior tag SAMBLASTER when doing aggregateRemarkDups
    samtools view -Sb - > $TMPDIR/${name}.${mappedgenome}.markedDups.bam
    samtools sort -@ $NSLOTS -O bam -m 5000M -T $TMPDIR/${name}.sort -l 9 $TMPDIR/${name}.${mappedgenome}.markedDups.bam > ${sampleOutdir}/${name}.${mappedgenome}.markedDups.bam && mv ${sampleOutdir}/${name}.${mappedgenome}.markedDups.bam ${sampleOutdir}/${name}.${mappedgenome}.bam
fi


echo
echo "Indexing"
date
samtools index ${sampleOutdir}/${name}.${mappedgenome}.bam


if [[ "${processingCommand}" =~ ^map ]]; then
    echo "Removing source files"
    rm -f ${files}
fi


echo
echo -e "\nDone!"
date
