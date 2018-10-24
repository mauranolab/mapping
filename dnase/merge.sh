#!/bin/bash
set -eu -o pipefail

analysisType=$1
name=$2
DS=$3
mappedgenome=$4

processingCommand=`echo "${analysisType}" | awk -F "," '{print $1}'`
#analysisCommand=`echo "${analysisType}" | awk -F "," '{print $2}'`

echo "Merging bams (${processingCommand} analysis); output name ${name} (${DS}) against genome ${mappedgenome}"

date

sampleOutdir=${name}

files=""
numfiles=0
if [[ "${processingCommand}" =~ ^map ]]; then
    #Here we generate a list of expected filenames
    #TODO switch to wildcard matching based on sample, DS and mappedgenome?
    #NB originally checked grep R1 but some old solexa lanes didn't include it
    for f1 in `grep ${DS} inputs.txt | grep -v _R2`; do
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
    for curOutputFile in `grep ${DS} inputs.txt | grep ${mappedgenome}`; do
        if [[ -f "${curOutputFile}" ]]; then
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
    samtools merge -l 1 -@ $NSLOTS ${sampleOutdir}/${name}.${mappedgenome}.bam ${files}
fi

echo "Processing ${name}.${mappedgenome}.bam"
#TODO AddOrReplaceReadGroups or do bwa -r ""


if [[ "${processingCommand}" == "aggregateRemarkDups" ]] || [[ "${processingCommand}" =~ ^map ]]; then
    echo
    echo "mark dups"
    date
    ###Obsolete picard version
    #Used to need VALIDATION_STRINGENCY=LENIENT to avoid SAM validation error: ERROR...MAPQ should be 0 for unmapped read or CIGAR should have zero elements for unmapped read
    #http://seqanswers.com/forums/showthread.php?t=4246
    #BUGBUG Can make huge log files despite these options
    #Now getting new error: "Ignoring SAM validation error: ERROR: Record 11792763, Read name HISEQ-2500-1:94:C74YHANXX:2:1109:6562:81875, bin field of BAM record does not equal value computed based on alignment start and end, and length of sequence to which read is aligned"
    #http://sourceforge.net/p/samtools/mailman/message/32910359/
    #java -Xmx6g -jar /home/maurano/bin/picard-tools/MarkDuplicates.jar INPUT=${name}.${mappedgenome}.bam OUTPUT=${name}.markedDups.bam METRICS_FILE=$TMPDIR/${name}.picardDups.txt QUIET=TRUE VERBOSITY=ERROR COMPRESSION_LEVEL=9 ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT && mv ${name}.markedDups.bam ${name}.${mappedgenome}.bam

    ###Samblaster is faster
    #samblaster used an average of 1GB memory mapping ENCODE DNase data to hg38. 10/889 jobs used >5GB.
    samtools sort -@ $NSLOTS -O sam -m 1750M -T $TMPDIR/${name}.sortbyname -l 1 -n ${sampleOutdir}/${name}.${mappedgenome}.bam |
    samblaster |
    samtools view -Sb - > $TMPDIR/${name}.${mappedgenome}.markedDups.bam
    samtools sort -@ $NSLOTS -O bam -m 2500M -T $TMPDIR/${name}.sort -l 9 $TMPDIR/${name}.${mappedgenome}.markedDups.bam > ${sampleOutdir}/${name}.${mappedgenome}.markedDups.bam && mv ${sampleOutdir}/${name}.${mappedgenome}.markedDups.bam ${sampleOutdir}/${name}.${mappedgenome}.bam
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
