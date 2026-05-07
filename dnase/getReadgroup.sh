#!/bin/bash
set -u

BS=${1}
fc=${2}
src=${3}

#Will set a variable ${readgroup} for input into mapper

echo "[getReadgroup] looking up BS \"${BS}\" on FC \"${fc}\""

BS_nosuffix=`echo "${BS}" | perl -pe 's/[A-Z]$//g;'`


readgroup="@RG\\tID:${fc}${BS}\\tLB:${BS}\\tSM:${BS_nosuffix}"


FlowcellInfoFile="/gpfs/data/isg_sequencing/data/${fc}/info.txt"
#Only use LIMS-derived info sheet if lims.py succeeds; otherwise fall back on hardcoded info.txt
${src}/lims.py --getFCinfo "${fc}" > $TMPDIR/info.txt && FlowcellInfoFile="$TMPDIR/info.txt"

if [ -s "${FlowcellInfoFile}" ]; then
    instrument=`awk -F "\t" 'BEGIN {OFS="\t"} $1=="#Instrument" {print $2}' ${FlowcellInfoFile}`
    
    readgroup_date=`awk -F "\t" 'BEGIN {OFS="\t"; loaddate="NA"} $1=="#Load date" && $2!="" {loaddate=$2} END {print "DT:" loaddate}' ${FlowcellInfoFile}`
    #BUGBUG hardcoded column numbers
    sample_bcs=`awk -v ds=${BS} -F "\t" 'BEGIN {OFS="\t"} $0!~/^#/ && 0!="" && $2==ds {split($6, bc1, "_"); split($7, bc2, "_"); if(bc2[2]=="") {print bc1[2]} else {print bc1[2] "-" bc2[2]}}' ${FlowcellInfoFile}`
    
    case "${instrument}" in
    Balin)
        readgroup_instrument="PL:ILLUMINA\\tCN:Maurano_Lab\\tPM:NextSeq_500"
        ;;
    Gandalf)
        readgroup_instrument="PL:ILLUMINA\\tCN:Maurano_Lab\\tPM:NextSeq_2000"
        ;;
    ISG_GridION)
        readgroup_instrument="PL:ONT\\tCN:NYUMC_ISG\\tPM:GridION"
        ;;
    MSSM_Revio)
        readgroup_instrument="PL:PACBIO\\tCN:MSSM\\tPM:REVIO"
        ;;
    GTC_Revio)
        readgroup_instrument="PL:PACBIO\\tCN:NYUMC_GTC\\tPM:REVIO"
        ;;
    GTC_NovaSeq)
        readgroup_instrument="PL:ILLUMINA\\tCN:NYUMC_GTC\\tPM:NovaSeq_6000"
        ;;
    GTC_NextSeq)
        readgroup_instrument="PL:ILLUMINA\\tCN:NYUMC_GTC\\tPM:NextSeq_500"
        ;;
    GTC_MiSeq)
        readgroup_instrument="PL:ILLUMINA\\tCN:NYUMC_GTC\\tPM:MiSeq"
        ;;
    esac
    
    #BUGBUG BC: shows up in bwa command line but at some point disappears from the bam header
    readgroup="${readgroup}\\t${readgroup_instrument}\\t${readgroup_date}\\tBC:${sample_bcs}\\tPU:${fc/./}-${sample_bcs}"
fi

echo "[getReadgroup] ${readgroup}"
