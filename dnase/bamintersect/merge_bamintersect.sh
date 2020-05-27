#!/bin/bash
set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'

sampleOutdir=$1
src=$2
sample_name=$3
bam1genome=$4
bam2genome=$5
make_table=$6
INTERMEDIATEDIR=$7
bam1=$8
bam2=$9
verbose=$10
make_bed=${11-}          # BUGBUG Is blank if --make_bed is unset via --no_bed in submit_bamintersect. Then ReqFullyAligned may be shifted.
ReqFullyAligned=${12-}   # Is blank if --noReqFullyAligned in submit_bamintersect

## This function implements verbose output for debugging purposes.
debug_fa() {
    if [ ${verbose} = "True" ]; then
        echo "[DEBUG] ${1}"
    fi
}


# Exclude reads with flags: read unmapped, mate unmapped, read is duplicate, or read is sup alignment
assemblyBackboneExcludeFlags=3084
# Exclude reads with flags: read unmapped, mate unmapped, failed QC flag, read is duplicate, or read is sup alignment
homologyArmExcludeFlags=3596


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
if echo "${bam2genome}" | egrep -q "^LP[0-9]+$"; then
    uninformativeRegionFiles="${src}/LP_uninformative_regions/LP_vs_${bam1genome}.bed ${INTERMEDIATEDIR}/HA_coords.bed"
#TODO masks hardcoded by payload name for now
elif echo "${bam2genome}" | egrep -q "^(Sox2_|PL1)"; then
    uninformativeRegionFiles="${src}/LP_uninformative_regions/PL_vs_LP.bed ${INTERMEDIATEDIR}/deletion_range.bed"
elif echo "${bam2genome}" | egrep -q "^rtTA$"; then
    uninformativeRegionFiles="${src}/LP_uninformative_regions/rtTA_vs_mm10.bed ${INTERMEDIATEDIR}/deletion_range.bed"
elif echo "${bam2genome}" | egrep -q "^pSpCas9"; then
    #These files mask just mm10/hg38 right now, so ok to use the same for all the pSpCas9 derivative constructs
    uninformativeRegionFiles="${src}/LP_uninformative_regions/pSpCas9_vs_${bam1genome}.bed"
elif echo "${bam2genome}" | egrep -q "^(Hoxa_|HPRT1)"; then
    uninformativeRegionFiles="${src}/LP_uninformative_regions/PL_vs_LPICE.bed ${INTERMEDIATEDIR}/deletion_range.bed"
    if [[ "${bam1genome}" == "LPICE" ]]; then
        #Need also to mask the SV40pA on the assembly
        if [ ! -f "/vol/cegs/sequences/cegsvectors_${bam2genome}/cegsvectors_${bam2genome}.bed" ]; then
            echo "WARNING could not find bed file to mask SV40 poly(A) signal in ${bam2genome} genome"
        else
            awk -F "\t" '$4=="SV40 poly(A) signal"' /vol/cegs/sequences/cegsvectors_${bam2genome}/cegsvectors_${bam2genome}.bed > $TMPDIR/ICE_SV40pA.bed
            uninformativeRegionFiles="${uninformativeRegionFiles} $TMPDIR/ICE_SV40pA.bed"
        fi
    fi
else
    uninformativeRegionFiles="${src}/LP_uninformative_regions/${bam2genome}_vs_${bam1genome}.bed"
fi

#Genomic repeat annotation
if [ -f "/vol/isg/annotation/bed/${bam1genome}/repeat_masker/Satellite.bed" ]; then
    uninformativeRegionFiles="${uninformativeRegionFiles} /vol/isg/annotation/bed/${bam1genome}/repeat_masker/Satellite.bed"
fi
if [ -f "/vol/isg/annotation/bed/${bam2genome}/repeat_masker/Satellite.bed" ]; then
    uninformativeRegionFiles="${uninformativeRegionFiles} /vol/isg/annotation/bed/${bam2genome}/repeat_masker/Satellite.bed"
fi

for curfile in ${uninformativeRegionFiles}; do
    if [ ! -f "${curfile}" ]; then
        echo "ERROR: Can't find ${curfile}; will continue without any region mask"
        exit 1
    else
        echo "Found uninformative regions file: $(basename ${curfile}) [${curfile}]"
    fi
done

#Sort exclusion files and strip comments just in case
cat ${uninformativeRegionFiles} | awk '$0 !~ /^#/' | sort-bed - > ${TMPDIR}/uninformativeRegionFile.bed


echo "Applying uninformativeRegionFile"
#Starts out with bam1 coordinates
cat ${sampleOutdir}/${sample_name}.bed | 
#Remove bam1 reads that overlap the ranges defined in uninformativeRegionFile.bed file via submit_bamintersect.sh
#Require 20 bp mapped outside uninformative region
bedmap --delim "|" --echo --bases-uniq - ${TMPDIR}/uninformativeRegionFile.bed | awk -v minUniqBp=20 -F "|" 'BEGIN {OFS="\t"} {split($1, read, "\t"); readlength=read[3]-read[2]; if(readlength-$2>minUniqBp) {print $1}}' |
#Switch to bam2
awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | sort-bed - |
#Sort by the bam2 reads, then keep only the bam2 reads (and their bam1 mates) that do not overlap the uninformativeRegionFile:
#Require 20 bp mapped outside uninformative region
bedmap --delim "|" --echo --bases-uniq - ${TMPDIR}/uninformativeRegionFile.bed | awk -v minUniqBp=20 -F "|" 'BEGIN {OFS="\t"} {split($1, read, "\t"); readlength=read[3]-read[2]; if(readlength-$2>minUniqBp) {print $1}}' |
#Switch back to bam1 coordinates
awk -F "\t" 'BEGIN {OFS="\t"}; {print $7, $8, $9, $10, $11, $12, $1, $2, $3, $4, $5, $6}' | sort-bed - > ${sampleOutdir}/${sample_name}.informative.bed

echo
${src}/counts_table.sh ${sampleOutdir}/${sample_name}.informative ${sample_name} ${bam1genome} ${bam2genome} ${sampleOutdir}/${sample_name}.informative.bed


# HA analysis:
echo
#Al
if [ -s "${INTERMEDIATEDIR}/HA_coords.bed" ]; then
    echo -e "Starting HA analysis."
    
    ${src}/HA_table.sh HA ${bam1} ${sampleOutdir} ${sample_name} ${INTERMEDIATEDIR} ${src} ${homologyArmExcludeFlags} ${ReqFullyAligned}
    ${src}/counts_table.sh ${sampleOutdir}/${sample_name}.HA "${sample_name}.HA" ${bam1genome} ${bam1genome} ${sampleOutdir}/${sample_name}.HA.bed
else
    echo -e "No HAs available, so there will be no HA analysis."
fi


echo
echo "Look for reads spanning the assembly and the backbone."
${src}/HA_table.sh AB ${bam2} ${sampleOutdir} ${sample_name} $INTERMEDIATEDIR{} ${src} ${assemblyBackboneExcludeFlags} ${ReqFullyAligned}
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
echo "Merging bam files"
bamfiles=`cat ${INTERMEDIATEDIR}/inputs.bamintersect.txt | cut -f3 | perl -pe 's/\.bed$/.bam1.bam/g;'; cat ${INTERMEDIATEDIR}/inputs.bamintersect.txt | cut -f3 | perl -pe 's/\.bed$/.bam2.bam/g;'`
samtools merge -@ $NSLOTS -c -p -O bam -l 1 $TMPDIR/${sample_name}.bam ${bamfiles}
cut -f4 ${sampleOutdir}/${sample_name}.informative.bed > ${TMPDIR}/${sample_name}.readNames.txt
if [ -s "${TMPDIR}/${sample_name}.readNames.txt" ]; then
    debug_fa "${sample_name}.readNames.txt was found"
    
    ${src}/subsetBAM.py --readNames ${TMPDIR}/${sample_name}.readNames.txt --bamFile_in $TMPDIR/${sample_name}.bam --bamFile_out ${sampleOutdir}/${sample_name}.informative.bam
    samtools index ${sampleOutdir}/${sample_name}.informative.bam
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
echo "track name=\"${sample_name}\" description=\"${sample_name} Reads\" visibility=pack pairEndsByName=F maxItems=10000 type=bam ${UCSCbase}/${sample_name}.informative.bam"


echo
echo "Subset bam files for assemblyBackbone"
cut -f4 ${sampleOutdir}/${sample_name}.assemblyBackbone.bed > ${TMPDIR}/${sample_name}.assemblyBackbone.readNames.txt

if [ -s "${TMPDIR}/${sample_name}.assemblyBackbone.readNames.txt" ]; then
    ${src}/subsetBAM.py --exclude_flags ${assemblyBackboneExcludeFlags} --readNames ${TMPDIR}/${sample_name}.assemblyBackbone.readNames.txt --bamFile_in ${bam1} --bamFile_out - |
    samtools sort -@ $NSLOTS -O bam -m 4000M -T $TMPDIR/${sample_name}.sort -l 9 -o ${sampleOutdir}/${sample_name}.assemblyBackbone.bam
    samtools index ${sampleOutdir}/${sample_name}.assemblyBackbone.bam
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

