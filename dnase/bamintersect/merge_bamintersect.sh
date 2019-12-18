#!/bin/bash
set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'

## This function implements verbose output for debugging purposes.
debug_fa() {
    if [ ${verbose} = "True" ]; then
        echo "[DEBUG] ${1}"
    fi
}
##########################################################################################################


## Variables are passed in via sbatch export.

## Merge the bam_intersect output files:

debug_fa "Starting merge_bamintersect.sh"

# First check for the "no reads to merge" problem. Also deal with the usual pipefail issue via "true".
dir_names=($(ls -d ${INTERMEDIATEDIR}/bamintersectPyOut/*/ 2> /dev/null || true))      ## These have a trailing /
numElements="${#dir_names[@]}"
if [ "${numElements}" = "0" ]; then
    # Don't generate an error. Create some appropriately empty output files.
    echo "" > "${sampleOutdir}/${sample_name}.bed"
    echo "" > "${sampleOutdir}/${sample_name}.informative.bed"
    echo "There are no reads with unmapped mates to analyze."
    # Not generating bam files for this case.
    
    ## Cleanup
    rm -r ${INTERMEDIATEDIR}
    echo "Done!!!"
    date
    exit 0
fi

first="True"
for i in ${dir_names[@]}; do
    if [ ${first} = "True" ]; then
        cp "${i}dsgrep_out.csv" "${INTERMEDIATEDIR}/unsorted_dsgrep.bed"
        first="False"
    else
        cat "${i}dsgrep_out.csv" >> "${INTERMEDIATEDIR}/unsorted_dsgrep.bed"
    fi
done

sort-bed ${INTERMEDIATEDIR}/unsorted_dsgrep.bed > ${sampleOutdir}/${sample_name}.bed

debug_fa "Finished sorting dsgrep.bed"

##########################################################################################################
## Create the final output tables.

if [ ${make_bed} != "--make_bed" ]; then
    # Need the bed file to make the table.
    # Maybe make both bed and bam files, and delete bed file later?
    echo "Made 2 bam files, so exiting merge_bamintersect.sh"
    exit 0
fi

if [ ${make_table} != "True" ]; then
    echo "LP integration table not requested; quitting successfully."
    exit 0
fi


#Unfiltered
${src}/counts_table.sh ${sampleOutdir}/${sample_name} ${sample_name} ${bam2genome} ${sampleOutdir}/${sample_name}.bed
echo


#With filters
#Check for a curated uninformative regions file for this combination of genomes
if echo "${bam1genome}" | grep -q "^LP[0-9]\+" ; then
    genome2exclude="${src}/LP_uninformative_regions/LP_vs_${bam2genome}.bed"
#TODO hardcoded for now
elif echo "${bam1genome}" | grep -q "^Sox2_" ; then
    genome2exclude="${src}/LP_uninformative_regions/PL_vs_LP.bed"
else
    genome2exclude="${src}/LP_uninformative_regions/${bam1genome}_vs_${bam2genome}.bed"
fi
if [ ! -f "${genome2exclude}" ]; then
    echo "WARNING: Can't find ${genome2exclude}; will run without region mask"
    genome2exclude=""
else
    echo "Found uninformative regions file: $(basename ${genome2exclude}) [${genome2exclude}]" >> ${sampleOutdir}/${sample_name}.info.txt
fi
#Sort exclusion files and strip comments just in case
#BUGBUG hardcoded repeatmask for bam2genome
cat ${genome2exclude} ${INTERMEDIATEDIR}/HA_coords.bed | awk '$0 !~ /^#/' | sort-bed - | bedops -u - /vol/isg/annotation/bed/${bam2genome}/repeat_masker/Satellite.bed > ${TMPDIR}/uninformativeRegionFile.bed

#Starts out with bam1 coordinates
cat ${sampleOutdir}/${sample_name}.bed | 
#Remove bam1 reads that overlap the ranges defined in uninformativeRegionFile.bed file via submit_bamintrsect.sh
bedops -n 1 - ${TMPDIR}/uninformativeRegionFile.bed |
#Switch to bam2
awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | sort-bed - |
#Sort by the bam2 reads, then keep only the bam2 reads (and their bam1 mates) that do not overlap the uninformativeRegionFile:
bedops -n 1 - ${TMPDIR}/uninformativeRegionFile.bed |
#Only apply deletion_range to bam2
bedops -n 1 - ${INTERMEDIATEDIR}/deletion_range.bed | 
#Switch back to bam1 coordinates
awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | sort-bed - > ${sampleOutdir}/${sample_name}.informative.bed

${src}/counts_table.sh ${sampleOutdir}/${sample_name}.informative ${sample_name} ${bam2genome} ${sampleOutdir}/${sample_name}.informative.bed
echo


# HA analysis:
echo
if [ -s "${INTERMEDIATEDIR}/HA_coords.bed" ]; then
    echo -e "Starting HA analysis."
    ${src}/HA_table.sh ${bam2} 3076 ${sampleOutdir} ${sample_name} ${INTERMEDIATEDIR} ${src}
    ${src}/counts_table.sh ${sampleOutdir}/${sample_name}_HA "${sample_name}_HA" ${bam2genome} ${sampleOutdir}/${sample_name}_HA.bed
else
    echo -e "No HAs available, so there will be no HA analysis."
    
    # Create empty files.
    touch ${sampleOutdir}/${sample_name}_HA.counts.txt
    touch ${sampleOutdir}/${sample_name}_HA.bed
fi
echo


echo
echo "Look for reads spanning the assembly and the backbone."
${src}/assemblyBackbone_table.sh ${bam1} 3076 ${sampleOutdir} ${sample_name} $INTERMEDIATEDIR{} ${src}
${src}/counts_table.sh ${sampleOutdir}/${sample_name}.assemblyBackbone "${sample_name}.assemblyBackbone" ${bam2genome} ${sampleOutdir}/${sample_name}.assemblyBackbone.bed
echo

# Number of reads
n1=$(wc -l < "${sampleOutdir}/${sample_name}.bed")
# Number of informative reads
n2=$(wc -l < "${sampleOutdir}/${sample_name}.informative.bed")

echo -e "Summary:" >> ${sampleOutdir}/${sample_name}.info.txt
echo -e "Mapped reads with unmapped mates:\t${n1}\tof which\t${n2}\tare potentially informative." >> ${sampleOutdir}/${sample_name}.info.txt

###########################################################################################################################################
echo
echo "Subset bam files for uploading UCSC custom tracks"
cut -f4 ${sampleOutdir}/${sample_name}.assemblyBackbone.bed > ${TMPDIR}/${sample_name}.assemblyBackbone.readNames.txt

if [ -s "${TMPDIR}/${sample_name}.assemblyBackbone.readNames.txt" ]; then
    ${src}/subsetBAM.py --flags ${BAM_E1} --readNames ${TMPDIR}/${sample_name}.assemblyBackbone.readNames.txt --bamFile_in ${bam1} --bamFile_out - |
    samtools sort -@ $NSLOTS -O bam -m 5000M -T $TMPDIR/${sample_name}.sort -l 9 -o ${sampleOutdir}/${sample_name}.assemblyBackbone.bam
    samtools index ${sampleOutdir}/${sample_name}.assemblyBackbone.bam
fi

cut -f4 ${sampleOutdir}/${sample_name}.informative.bed > ${TMPDIR}/${sample_name}.readNames.txt
if [ -s "${TMPDIR}/${sample_name}.readNames.txt" ]; then
    debug_fa "${sample_name}.readNames.txt was found"
    
    ${src}/subsetBAM.py --flags ${BAM_E1} --readNames ${TMPDIR}/${sample_name}.readNames.txt --bamFile_in ${bam1} --bamFile_out - |
    samtools sort -@ $NSLOTS -O bam -m 5000M -T $TMPDIR/${sample_name}.sort -l 1 -o ${TMPDIR}/${sample_name}.informative.${bam1genome}.bam
    debug_fa "Completed samtools lines for bam1"
    
    ${src}/subsetBAM.py --flags ${BAM_E2} --readNames ${TMPDIR}/${sample_name}.readNames.txt --bamFile_in ${bam2} --bamFile_out - |
    samtools sort -@ $NSLOTS -O bam -m 5000M -T $TMPDIR/${sample_name}.sort -l 1 -o ${TMPDIR}/${sample_name}.informative.${bam2genome}.bam
    debug_fa "Completed samtools lines for bam2"
    
    samtools merge -@ $NSLOTS -O bam -l 9 ${sampleOutdir}/${sample_name}.informative.bam ${TMPDIR}/${sample_name}.informative.${bam1genome}.bam ${TMPDIR}/${sample_name}.informative.${bam2genome}.bam
    samtools index ${sampleOutdir}/${sample_name}.informative.bam
    
    #Prep for UCSC track links
    #Remove "new" from the end of path so that we can reprocess data without affecting live data
    projectdir=`pwd | perl -pe 's/^\/vol\/(cegs|mauranolab|isg\/encode)\///g;' | perl -pe 's/\/new$//g;'`
    if [[ `pwd` =~ ^\/vol\/cegs\/ ]]; then
        UCSCbase="bigDataUrl=https://cegs@cascade.isg.med.nyu.edu/cegs/${projectdir}/${sampleOutdir}"
    elif [[ `pwd` =~ ^\/vol\/isg\/encode\/ ]]; then
        UCSCbase="bigDataUrl=https://cascade.isg.med.nyu.edu/mauranolab/encode/${projectdir}/${sampleOutdir}"
    else
        UCSCbase="bigDataUrl=https://mauranolab@cascade.isg.med.nyu.edu/~mauram01/${projectdir}/${sampleOutdir}"
    fi
    
    echo "" >> ${sampleOutdir}/${sample_name}.info.txt
    echo "Paths for custom tracks to bam files with informative reads:" >> ${sampleOutdir}/${sample_name}.info.txt
    echo "track name=\"${sample_name}\" description=\"${sample_name} Reads\" visibility=pack pairEndsByName=F maxWindowToDraw=10000 maxItems=250 type=bam ${UCSCbase}/${sample_name}.informative.bam" >> ${sampleOutdir}/${sample_name}.info.txt
    
else
    # No informative reads
    echo "For ${sample_name}, there are no bam files with informative reads."
fi
echo

debug_fa "Done with the samtools section."

#############################################################################
## Cleanup
rm -r ${INTERMEDIATEDIR}
#############################################################################

echo
echo "Done!!!"
date

