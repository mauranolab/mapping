#!/bin/bash
set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'

bam1=${1}
bam1genome=${2}
bam2=${3}
bam2genome=${4}
normbam=${5}
sampleOutdir=${6}
sample_name=${7}
verbose=${8}
INTERMEDIATEDIR=${9}
src=${10}

## This function implements verbose output for debugging purposes.
debug_fa() {
    if [ ${verbose} = "True" ]; then
        echo "[DEBUG] ${1}"
    fi
}


# Exclude reads with flags: read unmapped, mate unmapped, read is duplicate, or read is sup alignment
assemblyBackboneExcludeFlags=3084
# Exclude reads with flags: read unmapped, mate unmapped, read is duplicate, or read is sup alignment
HAExcludeFlags=3084


###################

echo "Running on $HOSTNAME. Using $TMPDIR as tmp"


## Variables are passed in via sbatch export.

echo "Merging bamintersect output files"
#xargs and bedops -u don't work well for some reason, so pipe everything through cat first; performance shouldn't be a huge deal
cat ${INTERMEDIATEDIR}/inputs.bamintersect.txt | cut -f3 | xargs cat | sort-bed - > ${TMPDIR}/${sample_name}.bed


if [ ! -s ${TMPDIR}/${sample_name}.bed ]; then
    # Don't generate an error. Create some appropriately empty output files.
    echo "" > ${sampleOutdir}/${sample_name}.bed
    echo "There are no reads left to analyze."
    # Not generating bam files for this case.
    
    ## Cleanup
    rm -r ${INTERMEDIATEDIR}
    echo "Done!!!"
    date
    exit 0
fi

debug_fa "Finished sorting dsgrep.bed"

###################
## Create the final output tables.
echo "Normalizing to 10M reads using read depth from ${normbam}"
num_bam1_reads=$(samtools view -c -F 512 ${normbam})

echo "Applying counts_table_mask"
#Starts out with bam1 coordinates
cat ${TMPDIR}/${sample_name}.bed |
#Remove bam1 reads that overlap the ranges defined in counts_table_mask.bed file via submit_bamintersect.sh
#Require 20 bp of read to be mapped outside masked region
bedmap --delim "|" --echo --bases-uniq - ${INTERMEDIATEDIR}/counts_table_mask.bed | awk -v minUniqBp=20 -F "|" 'BEGIN {OFS="\t"} {split($1, read, "\t"); readlength=read[3]-read[2]; if(readlength-$2>minUniqBp) {print $1}}' |
#Switch to bam2
awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | sort-bed - |
#Sort by the bam2 reads, then keep only the bam2 reads (and their bam1 mates) that do not overlap the counts_table_mask:
#Require 20 bp of read to be mapped outside masked region
bedmap --delim "|" --echo --bases-uniq - ${INTERMEDIATEDIR}/counts_table_mask.bed | awk -v minUniqBp=20 -F "|" 'BEGIN {OFS="\t"} {split($1, read, "\t"); readlength=read[3]-read[2]; if(readlength-$2>minUniqBp) {print $1}}' |
#Switch back to bam1 coordinates
awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | sort-bed - > ${sampleOutdir}/${sample_name}.bed


# Number of reads
n1=$(wc -l < ${TMPDIR}/${sample_name}.bed)
# Number of informative reads
n2=$(wc -l < ${sampleOutdir}/${sample_name}.bed)

echo
echo -e "Mapped reads with unmapped mates:\t${n1}\tof which\t${n2}\tare potentially informative."


echo
${src}/counts_table.sh ${sampleOutdir}/${sample_name}.bed ${bam1genome} ${bam2genome} ${num_bam1_reads} ${sampleOutdir}/${sample_name}.informative ${sample_name}


if [ -s "${INTERMEDIATEDIR}/HA_coords.bed" ]; then
    echo
    echo -e "Starting HA analysis"
    
    #Do each HA separately so we can detect spurious junctions between them (i.e. head-to-tail or other integrants)
    ${src}/bedOverlapPE.sh ${bam1} ${INTERMEDIATEDIR}/HA1_coords.bed ${HAExcludeFlags} ${TMPDIR}/${sample_name}.HA1.bed ${src}
    ${src}/bedOverlapPE.sh ${bam1} ${INTERMEDIATEDIR}/HA2_coords.bed ${HAExcludeFlags} ${TMPDIR}/${sample_name}.HA2.bed ${src}
    bedops -u ${TMPDIR}/${sample_name}.HA1.bed ${TMPDIR}/${sample_name}.HA2.bed |
    #uniq is in case of reciprocal hits that appear in both tables
    uniq > ${sampleOutdir}/${sample_name}.HA.bed
    
    ${src}/counts_table.sh ${sampleOutdir}/${sample_name}.HA.bed ${bam1genome} ${bam1genome} ${num_bam1_reads} ${sampleOutdir}/${sample_name}.HA ${sample_name}.HA
else
    echo
    echo "No HAs available, so there will be no HA analysis."
    touch ${sampleOutdir}/${sample_name}.HA.bed
fi


echo
echo "Starting backbone analysis"
#Parse first line of idxstats, which is required to be the payload
read -r payload payload_length all_other <<< $(samtools idxstats ${bam2})
# Require exactly two chromosomes, one of which ends in "_backbone"
numChroms=$(samtools idxstats ${bam2} | awk -F "\t" 'BEGIN {OFS="\t"} $1!="*"' | wc -l)
if [[ "${numChroms}" != 2 ]]; then
    echo "${payload} has wrong number of chromosomes, skipping table."
    touch ${sampleOutdir}/${sample_name}.assemblyBackbone.bed
else
    backbone=$(samtools idxstats ${bam2} | awk -F "\t" '$1~/_backbone$/ {print $1}')
    if [ `samtools idxstats ${bam2} | awk -F "\t" '$1~/_backbone$/ {print $1}' | wc -l` != 1 ]; then
        echo "WARNING: Assembly does not have exactly one backbone chromosome, quitting successfully: ${payload} ${backbone}"
        echo ""
        echo "samtools idxstats ${bam2} :"
        samtools idxstats ${bam2}
        touch ${sampleOutdir}/${sample_name}.assemblyBackbone.bed
    else
        # Success: the naming conventions are being followed
        
        echo "payload is: ${payload}"
        echo "backbone is: ${backbone}"
        
        # Define region of interest as payload chromosome.
        echo -e "${payload}\t0\t${payload_length}" > ${INTERMEDIATEDIR}/payload.bed
        
        ${src}/bedOverlapPE.sh ${bam2} ${INTERMEDIATEDIR}/payload.bed ${assemblyBackboneExcludeFlags} ${sampleOutdir}/${sample_name}.assemblyBackbone.bed ${src}
        ${src}/counts_table.sh ${sampleOutdir}/${sample_name}.assemblyBackbone.bed ${bam2genome} ${bam2genome} ${num_bam1_reads} ${sampleOutdir}/${sample_name}.assemblyBackbone ${sample_name}.assemblyBackbone
    fi
fi
echo


####################################################
echo
echo "Merging bam files"
bamfiles=`cat ${INTERMEDIATEDIR}/inputs.bamintersect.txt | cut -f3 | perl -pe 's/\.bed$/.bam1.bam/g;'; cat ${INTERMEDIATEDIR}/inputs.bamintersect.txt | cut -f3 | perl -pe 's/\.bed$/.bam2.bam/g;'`
samtools merge -@ $NSLOTS -c -p -O bam -l 1 $TMPDIR/${sample_name}.bam ${bamfiles}
cut -f4 ${sampleOutdir}/${sample_name}.bed > ${TMPDIR}/${sample_name}.readNames.txt
if [ -s "${TMPDIR}/${sample_name}.readNames.txt" ]; then
    debug_fa "${sample_name}.readNames.txt was found"
    
    ${src}/subsetBAM.py --include_readnames ${TMPDIR}/${sample_name}.readNames.txt $TMPDIR/${sample_name}.bam - |
    samtools sort -@ $NSLOTS -O bam -m 4000M -T $TMPDIR/${sample_name}.sort -l 9 -o ${sampleOutdir}/${sample_name}.bam
    
    samtools index ${sampleOutdir}/${sample_name}.bam
fi

#Prep for UCSC track links
#Remove "new" from the end of path so that we can reprocess data without affecting live data
projectdir=`pwd | perl -pe 's/^\/vol\/(cegs|mauranolab|isg\/encode)\///g;' | perl -pe 's/\/new$//g;'`
if [[ `pwd` =~ ^\/vol\/cegs\/ ]]; then
    UCSCbase="bigDataUrl=https://cascade.isg.med.nyu.edu/cegs/trackhub/${projectdir}/${sampleOutdir}"
else
    UCSCbase="bigDataUrl=https://cascade.isg.med.nyu.edu/~mauram01/${projectdir}/${sampleOutdir}"
fi

echo ""
echo "Paths for custom tracks to bam files with informative reads:"
echo "track name=\"${sample_name}\" description=\"${sample_name} Reads\" visibility=pack pairEndsByName=F maxItems=10000 type=bam ${UCSCbase}/${sample_name}.bam"


echo
echo "Subset bam files for HA"
cut -f4 ${sampleOutdir}/${sample_name}.HA.bed > ${TMPDIR}/${sample_name}.HA.readNames.txt

if [ -s "${TMPDIR}/${sample_name}.HA.readNames.txt" ]; then
    ${src}/subsetBAM.py --exclude_flags ${HAExcludeFlags} --include_readnames ${TMPDIR}/${sample_name}.HA.readNames.txt ${bam1} - |
    samtools sort -@ $NSLOTS -O bam -m 4000M -T $TMPDIR/${sample_name}.sort -l 9 -o ${sampleOutdir}/${sample_name}.HA.bam
    samtools index ${sampleOutdir}/${sample_name}.HA.bam
fi
echo


echo
echo "Subset bam files for assemblyBackbone"
cut -f4 ${sampleOutdir}/${sample_name}.assemblyBackbone.bed > ${TMPDIR}/${sample_name}.assemblyBackbone.readNames.txt

if [ -s "${TMPDIR}/${sample_name}.assemblyBackbone.readNames.txt" ]; then
    ${src}/subsetBAM.py --exclude_flags ${assemblyBackboneExcludeFlags} --include_readnames ${TMPDIR}/${sample_name}.assemblyBackbone.readNames.txt ${bam1} - |
    samtools sort -@ $NSLOTS -O bam -m 4000M -T $TMPDIR/${sample_name}.sort -l 9 -o ${sampleOutdir}/${sample_name}.assemblyBackbone.bam
    samtools index ${sampleOutdir}/${sample_name}.assemblyBackbone.bam
fi
echo

debug_fa "Done with the samtools section."



## Cleanup
rm -r ${INTERMEDIATEDIR}

echo
echo "Done!!!"
date

