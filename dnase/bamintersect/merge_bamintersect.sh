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
    echo "There are no reads with unmapped mates to analyze." >> ${sampleOutdir}/${sample_name}.counts.anc_info.txt
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

sort-bed "${INTERMEDIATEDIR}/unsorted_dsgrep.bed" > "${sampleOutdir}/${sample_name}.bed"

debug_fa "Finished sorting dsgrep.bed"

##########################################################################################################
## Create the final output tables.

if [ ${make_csv} != "--make_csv" ]; then
    # Need the bed file to make the table.
    # Maybe make both bed and bam files, and delete bed file later?
    echo "Making 2 bam files, so exiting merge_bamintersect.sh prior to calling filter_tsv.bash"
    exit 0
fi

if [ ${make_table} != "True" ]; then
    echo "LP integration table not required."
    exit 0
fi

# First call to filter_tsv. Here we use all the reads.
echo "Running filter_tsv without filters." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
${src}/filter_tsv.sh ${sampleOutdir} ${sample_name} ${bam2genome} "${sampleOutdir}/${sample_name}.bed" all_reads_counts ${INTERMEDIATEDIR} null

# Second call to filter_tsv. Here we just use informative reads (not in the HAs, and not in the genome2exclude ranges.
echo "Running filter_tsv with HA filters (if any) and with the Exclude Regions file (which may be empty)." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
${src}/filter_tsv.sh ${sampleOutdir} ${sample_name} ${bam2genome} "${sampleOutdir}/${sample_name}.bed" informative_reads_counts ${INTERMEDIATEDIR} ${genome2exclude}

# HA analysis:
if [ -s "${INTERMEDIATEDIR}/genome1exclude.bed" ]; then
    echo -e "Starting HA analysis.\n" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    ${src}/HA_table.sh ${bamname2} 3076 ${sampleOutdir} ${sample_name} ${INTERMEDIATEDIR} ${src}
    ${src}/filter_tsv.sh ${sampleOutdir} "${sample_name}_HA" ${bam2genome} "${sampleOutdir}/${sample_name}_HA.bed" all_reads_counts ${INTERMEDIATEDIR} null
else
    echo -e "No HAs available, so there will be no HA analysis.\n" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    
    # Create empty files.
    touch "${sampleOutdir}/${sample_name}_HA.counts.txt"
    touch "${sampleOutdir}/${sample_name}_HA.bed"
    touch "${sampleOutdir}/${sample_name}_HA.counts.anc_info.txt"
fi

# Look for reads spanning the assembly and the backbone.
${src}/assmblyBackbone_table.sh ${bamname1} 3076 ${sampleOutdir} ${sample_name} ${INTERMEDIATEDIR} ${src}
${src}/filter_tsv.sh ${sampleOutdir} "${sample_name}_assmblyBackbone" ${bam2genome} "${sampleOutdir}/${sample_name}_assmblyBackbone.bed" all_reads_counts ${INTERMEDIATEDIR} null

# Number of reads
n1=$(wc -l < "${sampleOutdir}/${sample_name}.bed")
# Number of informative reads
n2=$(wc -l < "${sampleOutdir}/${sample_name}.informative.bed")

echo -e "Summary:" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo -e "Mapped reads with unmapped mates:\t${n1}\tof which\t${n2}\tpassed through the Exclude Regions filters and over the HAs (if any), and are potentially informative." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
echo -e "May include backbone reads, and reads not on the integration site chromosome. Some or all may not be in the \"informative.counts\" table, if they are over 500 bps from the nearest other informative read.\n" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"

###########################################################################################################################################
# Subset bam files for uploading UCSC custom tracks:
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
    
    echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    
    echo "Paths for custom tracks to bam files with informative reads:" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    
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
    
    echo "track name=\"${sample_name}.${bam1genome}-reads\" description=\"${sample_name}.${bam1genome} Reads\" visibility=pack pairEndsByName=F maxWindowToDraw=10000 maxItems=250 type=bam ${UCSCbase}/${sample_name}.informative.${bam1genome}.bam" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    
    echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    
    echo "track name=\"${sample_name}.${bam2genome}-reads\" description=\"${sample_name}.${bam2genome} Reads\" visibility=pack pairEndsByName=F maxWindowToDraw=10000 maxItems=250 type=bam ${UCSCbase}/${sample_name}.informative.${bam2genome}.bam" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    
    echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
else
    debug_fa "${sample_name}.readNames.txt was not found"
    
    # No informative reads
    echo "For ${sample_name}, there are no bam files with informative reads." >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
    echo "" >> "${sampleOutdir}/${sample_name}.counts.anc_info.txt"
fi

debug_fa "Done with the samtools section."

#############################################################################
## Cleanup
rm -r ${INTERMEDIATEDIR}
#############################################################################

echo
echo "Done!!!"
date

