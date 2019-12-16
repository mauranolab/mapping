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

if [ ${make_csv} != "--make_csv" ]; then
    # Need the bed file to make the table.
    # Maybe make both bed and bam files, and delete bed file later?
    echo "Making 2 bam files, so exiting merge_bamintersect.sh prior to calling counts_table.sh"
    exit 0
fi

if [ ${make_table} != "True" ]; then
    echo "LP integration table not requested; quitting successfully."
    exit 0
fi

#Unfiltered
${src}/counts_table.sh ${sampleOutdir}/${sample_name} ${sample_name} ${bam2genome} ${sampleOutdir}/${sample_name}.bed ${INTERMEDIATEDIR}
echo


#With filters
if [ "${uninformativeRegionFile}" != "null" ];then
    # Get reads that are in the deletion zone.
    
    #Starts out with bam1 coordinates
    cat ${sampleOutdir}/${sample_name}.bed | 
    #Remove bam1 reads that overlap the ranges defined in uninformativeRegionFile.bed file via submit_bamintrsect.sh
    bedops -n 1 - ${uninformativeRegionFile} |
    #Switch to bam2
    awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | sort-bed - |
    #Sort by the bam2 reads, then keep only the bam2 reads (and their bam1 mates) that do not overlap the uninformativeRegionFile:
    bedops -n 1 - ${uninformativeRegionFile} |
    #Only apply deletion_range to bam2
    bedops -n 1 - ${INTERMEDIATEDIR}/deletion_range.bed | 
    #Switch back to bam1 coordinates
    awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | sort-bed - > ${sampleOutdir}/${sample_name}.informative.bed
else
    cp ${sampleOutdir}/${sample_name}.bed ${sampleOutdir}/${sample_name}.informative.bed
fi
${src}/counts_table.sh ${sampleOutdir}/${sample_name}.informative ${sample_name} ${bam2genome} ${sampleOutdir}/${sample_name}.informative.bed ${INTERMEDIATEDIR}
echo


# HA analysis:
echo
if [ -s "${INTERMEDIATEDIR}/HA_coords.bed" ]; then
    echo -e "Starting HA analysis."
    ${src}/HA_table.sh ${bamname2} 3076 ${sampleOutdir} ${sample_name} ${INTERMEDIATEDIR} ${src}
    ${src}/counts_table.sh ${sampleOutdir}/${sample_name}_HA "${sample_name}_HA" ${bam2genome} ${sampleOutdir}/${sample_name}_HA.bed ${INTERMEDIATEDIR}
else
    echo -e "No HAs available, so there will be no HA analysis."
    
    # Create empty files.
    touch ${sampleOutdir}/${sample_name}_HA.counts.txt
    touch ${sampleOutdir}/${sample_name}_HA.bed
fi
echo


echo
echo "Look for reads spanning the assembly and the backbone."
${src}/assmblyBackbone_table.sh ${bamname1} 3076 ${sampleOutdir} ${sample_name} ${INTERMEDIATEDIR} ${src}
${src}/counts_table.sh ${sampleOutdir}/${sample_name}_assmblyBackbone "${sample_name}_assmblyBackbone" ${bam2genome} "${sampleOutdir}/${sample_name}_assmblyBackbone.bed" ${INTERMEDIATEDIR}
echo

# Number of reads
n1=$(wc -l < "${sampleOutdir}/${sample_name}.bed")
# Number of informative reads
n2=$(wc -l < "${sampleOutdir}/${sample_name}.informative.bed")

echo -e "Summary:" >> ${sampleOutdir}/${sample_name}.counts.anc_info.txt
echo -e "Mapped reads with unmapped mates:\t${n1}\tof which\t${n2}\tpassed through the Exclude Regions filters and over the HAs (if any), and are potentially informative." >> ${sampleOutdir}/${sample_name}.counts.anc_info.txt
echo -e "May include backbone reads, and reads not on the integration site chromosome. Some or all may not be in the \"informative.counts\" table, if they are over 500 bps from the nearest other informative read.\n" >> ${sampleOutdir}/${sample_name}.counts.anc_info.txt

###########################################################################################################################################
echo
echo "Subset bam files for uploading UCSC custom tracks"
cut -f4 "${sampleOutdir}/${sample_name}_assmblyBackbone.bed" > "${INTERMEDIATEDIR}/${sample_name}.readNames.txt"

if [ -s "${INTERMEDIATEDIR}/${sample_name}.readNames.txt" ]; then
    ${src}/subsetBAM.py --flags ${BAM_E1} --readNames ${INTERMEDIATEDIR}/${sample_name}.readNames.txt --bamFile_in ${bamname1} --bamFile_out ${TMPDIR}/${sample_name}.bam1_assmblyBackbone_reads.unsorted.bam
    samtools sort -@ $NSLOTS -O bam -m 5000M -T $TMPDIR/${sample_name}.sort -l 9 -o ${sampleOutdir}/${sample_name}.assmblyBackbone.${bam1genome}.bam ${TMPDIR}/${sample_name}.bam1_assmblyBackbone_reads.unsorted.bam
    samtools index ${sampleOutdir}/${sample_name}.assmblyBackbone.${bam1genome}.bam
fi

cut -f4 "${sampleOutdir}/${sample_name}.informative.bed" > "${INTERMEDIATEDIR}/${sample_name}.readNames.txt"
if [ -s "${INTERMEDIATEDIR}/${sample_name}.readNames.txt" ]; then
    debug_fa "${sample_name}.readNames.txt was found"
    
    ${src}/subsetBAM.py --flags ${BAM_E1} --readNames ${INTERMEDIATEDIR}/${sample_name}.readNames.txt --bamFile_in ${bamname1} --bamFile_out ${TMPDIR}/${sample_name}.bam1_informative_reads.unsorted.bam
    samtools sort -@ $NSLOTS -O bam -m 5000M -T $TMPDIR/${sample_name}.sort -l 9 -o ${sampleOutdir}/${sample_name}.informative.${bam1genome}.bam ${TMPDIR}/${sample_name}.bam1_informative_reads.unsorted.bam
    samtools index ${sampleOutdir}/${sample_name}.informative.${bam1genome}.bam
    debug_fa "Completed samtools lines for bam1"
    
    ${src}/subsetBAM.py --flags ${BAM_E2} --readNames ${INTERMEDIATEDIR}/${sample_name}.readNames.txt --bamFile_in ${bamname2} --bamFile_out ${TMPDIR}/${sample_name}.bam2_informative_reads.unsorted.bam
    samtools sort -@ $NSLOTS -O bam -m 5000M -T $TMPDIR/${sample_name}.sort -l 9 -o ${sampleOutdir}/${sample_name}.informative.${bam2genome}.bam ${TMPDIR}/${sample_name}.bam2_informative_reads.unsorted.bam
    samtools index ${sampleOutdir}/${sample_name}.informative.${bam2genome}.bam
    debug_fa "Completed samtools lines for bam2"
    
    
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
    
    echo "" >> ${sampleOutdir}/${sample_name}.counts.anc_info.txt
    echo "Paths for custom tracks to bam files with informative reads:" >> ${sampleOutdir}/${sample_name}.counts.anc_info.txt
    echo "track name=\"${sample_name}.${bam1genome}-reads\" description=\"${sample_name}.${bam1genome} Reads\" visibility=pack pairEndsByName=F maxWindowToDraw=10000 maxItems=250 type=bam ${UCSCbase}/${sample_name}.informative.${bam1genome}.bam" >> ${sampleOutdir}/${sample_name}.counts.anc_info.txt
    echo "track name=\"${sample_name}.${bam2genome}-reads\" description=\"${sample_name}.${bam2genome} Reads\" visibility=pack pairEndsByName=F maxWindowToDraw=10000 maxItems=250 type=bam ${UCSCbase}/${sample_name}.informative.${bam2genome}.bam" >> ${sampleOutdir}/${sample_name}.counts.anc_info.txt
    
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

