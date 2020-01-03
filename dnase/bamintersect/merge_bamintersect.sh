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

echo "Running on $HOSTNAME. Using $TMPDIR as tmp"


## Variables are passed in via sbatch export.

echo "Merging the bam_intersect output files"
#xargs and bedops -u don't work well for some reason, so pipe everything through cat first; performance shouldn't be a huge deal
cat ${INTERMEDIATEDIR}/inputs.bamintersect.txt | cut -f3 | xargs cat | sort-bed - > ${sampleOutdir}/${sample_name}.bed


if [ ! -s ${sampleOutdir}/${sample_name}.bed ]; then
    # Don't generate an error. Create some appropriately empty output files.
    echo "" > "${sampleOutdir}/${sample_name}.informative.bed"
    echo "There are no reads left to analyze."
    # Not generating bam files for this case.
    
    ## Cleanup
    rm -r ${INTERMEDIATEDIR}
    echo "Done!!!"
    date
    exit 0
fi


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
    echo "Summary counts tables not requested; quitting successfully."
    exit 0
fi


#Unfiltered
echo
${src}/counts_table.sh ${sampleOutdir}/${sample_name} ${sample_name} ${bam1genome} ${bam2genome} ${sampleOutdir}/${sample_name}.bed
echo


echo "Generating uninformativeRegionFile"
#Check for a curated uninformative regions file for this combination of genomes
#For LP integrations, exclude HAs
#For assemblon (not LP) integrations, filter out reads in deletion (since we have the mapping to the custom assembly)
if echo "${bam2genome}" | grep -q "^LP[0-9]\+$"; then
    uninformativeRegionFiles="${src}/LP_uninformative_regions/LP_vs_${bam1genome}.bed ${INTERMEDIATEDIR}/HA_coords.bed"
#TODO payload masks hardcoded by name for now
elif echo "${bam2genome}" | grep -q "^Sox2_"; then
    uninformativeRegionFiles="${src}/LP_uninformative_regions/PL_vs_LP.bed ${INTERMEDIATEDIR}/deletion_range.bed"
elif echo "${bam2genome}" | egrep -q "^(Hoxa_|HPRT1)"; then
    uninformativeRegionFiles="${src}/LP_uninformative_regions/LPICE_vs_mm10.bed ${INTERMEDIATEDIR}/deletion_range.bed"
else
    uninformativeRegionFiles="${src}/LP_uninformative_regions/${bam2genome}_vs_${bam1genome}.bed"
fi

for curfile in ${uninformativeRegionFiles}; do
    if [ ! -f "${curfile}" ]; then
        #NB one missing file will kill the whole list
        echo "WARNING: Can't find ${curfile}; will run without any region mask"
        uninformativeRegionFiles=""
    else
        echo "Found uninformative regions file: $(basename ${curfile}) [${curfile}]"
    fi
done


#Sort exclusion files and strip comments just in case
cat ${uninformativeRegionFiles} | awk '$0 !~ /^#/' | sort-bed - |
#BUGBUG hardcoded repeatmask for bam1genome
bedops -u - /vol/isg/annotation/bed/${bam1genome}/repeat_masker/Satellite.bed > ${TMPDIR}/uninformativeRegionFile.bed


echo "Applying uninformativeRegionFile"
#Starts out with bam1 coordinates
cat ${sampleOutdir}/${sample_name}.bed | 
#Remove bam1 reads that overlap the ranges defined in uninformativeRegionFile.bed file via submit_bamintersect.sh
bedops -n 1 - ${TMPDIR}/uninformativeRegionFile.bed |
#Switch to bam2
awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | sort-bed - |
#Sort by the bam2 reads, then keep only the bam2 reads (and their bam1 mates) that do not overlap the uninformativeRegionFile:
bedops -n 1 - ${TMPDIR}/uninformativeRegionFile.bed |
#Switch back to bam1 coordinates
awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | sort-bed - > ${sampleOutdir}/${sample_name}.informative.bed

echo
${src}/counts_table.sh ${sampleOutdir}/${sample_name}.informative ${sample_name} ${bam1genome} ${bam2genome} ${sampleOutdir}/${sample_name}.informative.bed


# HA analysis:
echo
#Al
if [ -s "${INTERMEDIATEDIR}/HA_coords.bed" ]; then
    echo -e "Starting HA analysis."

    # Exclude reads with flags: read unmapped, mate unmapped, failed QC flag, read is duplicate, or read is sup alignment (3596).
    ${src}/HA_table.sh HA ${bam1} ${sampleOutdir} ${sample_name} ${INTERMEDIATEDIR} ${src} 3596 ${ReqFullyAligned}
    ${src}/counts_table.sh ${sampleOutdir}/${sample_name}.HA "${sample_name}.HA" ${bam1genome} ${bam1genome} ${sampleOutdir}/${sample_name}.HA.bed
else
    echo -e "No HAs available, so there will be no HA analysis."
fi


echo
echo "Look for reads spanning the assembly and the backbone."
# Exclude reads with flags: read unmapped, mate unmapped, read is duplicate, or read is sup alignment (3084).
${src}/HA_table.sh AB ${bam2} ${sampleOutdir} ${sample_name} $INTERMEDIATEDIR{} ${src} 3084 ${ReqFullyAligned}
${src}/counts_table.sh ${sampleOutdir}/${sample_name}.assemblyBackbone "${sample_name}.assemblyBackbone" ${bam2genome} ${bam2genome} ${sampleOutdir}/${sample_name}.assemblyBackbone.bed
echo

# Number of reads
n1=$(wc -l < "${sampleOutdir}/${sample_name}.bed")
# Number of informative reads
n2=$(wc -l < "${sampleOutdir}/${sample_name}.informative.bed")

echo
echo -e "Mapped reads with unmapped mates:\t${n1}\tof which\t${n2}\tare potentially informative."

###########################################################################################################################################
echo
echo "Subset bam files for uploading UCSC custom tracks"
cut -f4 ${sampleOutdir}/${sample_name}.assemblyBackbone.bed > ${TMPDIR}/${sample_name}.assemblyBackbone.readNames.txt

if [ -s "${TMPDIR}/${sample_name}.assemblyBackbone.readNames.txt" ]; then
    ${src}/subsetBAM.py --exclude_flags ${BAM_E1} --readNames ${TMPDIR}/${sample_name}.assemblyBackbone.readNames.txt --bamFile_in ${bam1} --bamFile_out - |
    samtools sort -@ $NSLOTS -O bam -m 4000M -T $TMPDIR/${sample_name}.sort -l 9 -o ${sampleOutdir}/${sample_name}.assemblyBackbone.bam
    samtools index ${sampleOutdir}/${sample_name}.assemblyBackbone.bam
fi

cut -f4 ${sampleOutdir}/${sample_name}.informative.bed > ${TMPDIR}/${sample_name}.readNames.txt
if [ -s "${TMPDIR}/${sample_name}.readNames.txt" ]; then
    debug_fa "${sample_name}.readNames.txt was found"
    
    ${src}/subsetBAM.py --exclude_flags ${BAM_E1} --readNames ${TMPDIR}/${sample_name}.readNames.txt --bamFile_in ${bam1} --bamFile_out - |
    samtools sort -@ $NSLOTS -O bam -m 4000M -T $TMPDIR/${sample_name}.sort -l 1 -o ${TMPDIR}/${sample_name}.informative.${bam1genome}.bam
    debug_fa "Completed samtools lines for bam1"
    
    ${src}/subsetBAM.py --exclude_flags ${BAM_E2} --readNames ${TMPDIR}/${sample_name}.readNames.txt --bamFile_in ${bam2} --bamFile_out - |
    samtools sort -@ $NSLOTS -O bam -m 4000M -T $TMPDIR/${sample_name}.sort -l 1 -o ${TMPDIR}/${sample_name}.informative.${bam2genome}.bam
    debug_fa "Completed samtools lines for bam2"
    
    samtools merge -@ $NSLOTS -O bam -l 9 ${sampleOutdir}/${sample_name}.informative.bam ${TMPDIR}/${sample_name}.informative.${bam1genome}.bam ${TMPDIR}/${sample_name}.informative.${bam2genome}.bam
    samtools index ${sampleOutdir}/${sample_name}.informative.bam
    
    #Prep for UCSC track links
    #Remove "new" from the end of path so that we can reprocess data without affecting live data
    projectdir=`pwd | perl -pe 's/^\/vol\/(cegs|mauranolab|isg\/encode)\///g;' | perl -pe 's/\/new$//g;'`
    if [[ `pwd` =~ ^\/vol\/cegs\/ ]]; then
        UCSCbase="bigDataUrl=https://cegs@cascade.isg.med.nyu.edu/cegs/${projectdir}/${sampleOutdir}"
    else
        UCSCbase="bigDataUrl=https://mauranolab@cascade.isg.med.nyu.edu/~mauram01/${projectdir}/${sampleOutdir}"
    fi
    
    echo ""
    echo "Paths for custom tracks to bam files with informative reads:"
    echo "track name=\"${sample_name}\" description=\"${sample_name} Reads\" visibility=pack pairEndsByName=F maxWindowToDraw=10000 maxItems=250 type=bam ${UCSCbase}/${sample_name}.informative.bam"
    
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

