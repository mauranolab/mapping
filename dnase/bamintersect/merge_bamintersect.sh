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
num_bam1_reads=${5}
sampleOutdir=${6}
sample_name=${7}
verbose=${8}
INTERMEDIATEDIR=${9}
src=${10}


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
    echo "There are no reads left to analyze; quitting successfully"
    # Don't generate an error. Create some appropriately empty output files.
    # Not generating bam files for this case.
    touch ${sampleOutdir}/${sample_name}.bed
    touch ${sampleOutdir}/${sample_name}.HA.bed
    touch ${sampleOutdir}/${sample_name}.assemblyBackbone.bed
    
    echo -e "#chrom_bam1\tchromStart_bam1\tchromEnd_bam1\tWidth_bam1\tNearestGene_bam1\tStrand_bam1\tReadsPer10M\tchrom_bam2\tchromStart_bam2\tchromEnd_bam2\tWidth_bam2\tNearestGene_bam2\tStrand_bam2\tSample\tGenomes" > ${sampleOutdir}/${sample_name}.informative.counts.txt
    #TODO ideally wouldn't have bam1genome vs bam2genome suffix for consistency
    echo -e "#NA\tNA\tNA\t0\t-\t\t0\tNA\tNA\tNA\t0\t-\t\t${sample_name}\t${bam1genome}_vs_${bam2genome}" >> ${sampleOutdir}/${sample_name}.informative.counts.txt
    cp ${sampleOutdir}/${sample_name}.informative.counts.txt ${sampleOutdir}/${sample_name}.HA.counts.txt
    cp ${sampleOutdir}/${sample_name}.informative.counts.txt ${sampleOutdir}/${sample_name}.assemblyBackbone.counts.txt
    
    ## Cleanup
    rm -r ${INTERMEDIATEDIR}
    echo "Done!!!"
    date
    exit 0
fi


###################
## Create the final output tables.

applyUninformativeMaskFile()
{
    #Reads from stdin
    countsTableMaskFile=$1
    
    #Starts out with bam1 coordinates
    #Remove bam1 reads that overlap the ranges defined in countsTableMaskFile.bed file via submit_bamintersect.sh
    #Require 20 bp of read to be mapped outside masked region
    bedmap --delim "|" --echo --bases-uniq - ${countsTableMaskFile} | awk -v minUniqBp=20 -F "|" 'BEGIN {OFS="\t"} {split($1, read, "\t"); readlength=read[3]-read[2]; if(readlength-$2>minUniqBp) {print $1}}' |
    #Switch to bam2
    awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | sort-bed - |
    #Sort by the bam2 reads, then keep only the bam2 reads (and their bam1 mates) that do not overlap the countsTableMaskFile:
    #Require 20 bp of read to be mapped outside masked region
    bedmap --delim "|" --echo --bases-uniq - ${countsTableMaskFile} | awk -v minUniqBp=20 -F "|" 'BEGIN {OFS="\t"} {split($1, read, "\t"); readlength=read[3]-read[2]; if(readlength-$2>minUniqBp) {print $1}}' |
    #Switch back to bam1 coordinates
    awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | sort-bed -
}

echo "Normalizing to 10M reads using read depth of ${num_bam1_reads} reads."

#TODO is it really necessary to repeat the mask here?
cat ${TMPDIR}/${sample_name}.bed | applyUninformativeMaskFile ${sampleOutdir}/log/countsTableMaskFile.${bam1genome}_vs_${bam2genome}.bed > ${sampleOutdir}/${sample_name}.bed


# Number of reads
n1=$(wc -l < ${TMPDIR}/${sample_name}.bed)
# Number of informative reads
n2=$(wc -l < ${sampleOutdir}/${sample_name}.bed)

echo
echo -e "Mapped reads with unmapped mates: ${n1} of which ${n2} are potentially informative."


echo
${src}/counts_table.sh ${sampleOutdir}/${sample_name}.bed ${bam1genome} ${bam2genome} ${num_bam1_reads} ${sampleOutdir}/${sample_name}.informative ${sample_name}


#NB could be run in submit, but keeping it here for clarity
if [ -s "${sampleOutdir}/log/HA_coords.${bam1genome}_vs_${bam2genome}.bed" ]; then
    echo
    echo -e "Starting HA analysis"
    
    #Do each HA separately so we can detect spurious junctions between them (i.e. head-to-tail or other integrants)
    ${src}/bedOverlapPE.sh ${bam1} ${INTERMEDIATEDIR}/HA1_coords.bed ${HAExcludeFlags} ${TMPDIR}/${sample_name}.HA1.bed ${src}
    ${src}/bedOverlapPE.sh ${bam1} ${INTERMEDIATEDIR}/HA2_coords.bed ${HAExcludeFlags} ${TMPDIR}/${sample_name}.HA2.bed ${src}
    bedops -u ${TMPDIR}/${sample_name}.HA1.bed ${TMPDIR}/${sample_name}.HA2.bed |
    #uniq is in case of reciprocal hits that appear in both tables
    uniq |
    applyUninformativeMaskFile ${sampleOutdir}/log/HAmaskFile.${bam1genome}_vs_${bam2genome}.bed > ${sampleOutdir}/${sample_name}.HA.bed
    
    ${src}/counts_table.sh ${sampleOutdir}/${sample_name}.HA.bed ${bam1genome} ${bam1genome} ${num_bam1_reads} ${sampleOutdir}/${sample_name}.HA ${sample_name}.HA
else
    echo
    echo "No HAs available, so there will be no HA analysis."
    touch ${sampleOutdir}/${sample_name}.HA.bed
fi


#NB could be run in submit, but keeping it here for clarity
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
        echo -e "${payload}\t0\t${payload_length}" > ${TMPDIR}/payload.bed
        
        ${src}/bedOverlapPE.sh ${bam2} ${TMPDIR}/payload.bed ${assemblyBackboneExcludeFlags} ${TMPDIR}/${sample_name}.assemblyBackbone.bed ${src}
        cat ${TMPDIR}/${sample_name}.assemblyBackbone.bed |
        applyUninformativeMaskFile ${sampleOutdir}/log/uninformativeRegionFile.${bam1genome}_vs_${bam2genome}.bed > ${sampleOutdir}/${sample_name}.assemblyBackbone.bed
        
        ${src}/counts_table.sh ${sampleOutdir}/${sample_name}.assemblyBackbone.bed ${bam2genome} ${bam2genome} ${num_bam1_reads} ${sampleOutdir}/${sample_name}.assemblyBackbone ${sample_name}.assemblyBackbone
    fi
fi
echo


####################################################
subsetBamFileToInformativeReads()
{
    inbed=$1
    inbam=$2
    outbam=$3
    base=`basename $1 .bed`
    
    cut -f4 ${inbed} > ${TMPDIR}/${base}.readNames.txt
    if [ -s "${TMPDIR}/${base}.readNames.txt" ]; then
        ${src}/subsetBAM.py --exclude_flags ${HAExcludeFlags} --include_readnames ${TMPDIR}/${base}.readNames.txt ${inbam} - |
        samtools sort -@ $NSLOTS -O bam -m 4000M -T $TMPDIR/${base}.sort -l 9 -o ${outbam}
        samtools index ${outbam}
    fi
}


echo
echo "Merging bam files"
bamfiles=`cat ${INTERMEDIATEDIR}/inputs.bamintersect.txt | cut -f3 | perl -pe 's/\.bed$/.bam1.bam/g;'; cat ${INTERMEDIATEDIR}/inputs.bamintersect.txt | cut -f3 | perl -pe 's/\.bed$/.bam2.bam/g;'`
samtools merge -@ $NSLOTS -c -p -O bam -l 1 $TMPDIR/${sample_name}.bam ${bamfiles}
subsetBamFileToInformativeReads ${sampleOutdir}/${sample_name}.bed $TMPDIR/${sample_name}.bam ${sampleOutdir}/${sample_name}.bam

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


echo "Subset bam files for HA"
subsetBamFileToInformativeReads ${sampleOutdir}/${sample_name}.HA.bed ${bam1} ${sampleOutdir}/${sample_name}.HA.bam
echo


echo "Subset bam files for assemblyBackbone"
subsetBamFileToInformativeReads ${sampleOutdir}/${sample_name}.assemblyBackbone.bed ${bam2} ${sampleOutdir}/${sample_name}.assemblyBackbone.bam
echo


## Cleanup
rm -r ${INTERMEDIATEDIR}

echo
echo "Done!!!"
date

