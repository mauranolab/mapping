#!/bin/bash
set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'


###Parameters
mappedgenome=${1}
analysisType=${2}
sampleOutdir=${3}
sampleAnnotation=${4}
src=${5}

sampleType=`echo "${analysisType}" | awk -F "," '{print $2}'`

source ${src}/genomeinfo.sh ${mappedgenome}

jobid=$SGE_TASK_ID
#jobid=10

name=`basename ${sampleOutdir}`

jobname=`cat ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt | awk -F "\t" -v jobid=$jobid 'NR==jobid {print $1}'`
chroms=`cat ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt | awk -F "\t" -v jobid=$jobid 'NR==jobid {print $2}'`
echo "Running ${analysisType} analysis for ${chroms} of sample ${name} against genome ${mappedgenome}"
echo -e "SampleAnnotation\t${sampleAnnotation}"
date


#TMPDIR=`pwd`/tmp.callsnpsByChrom.${name}
#mkdir -p $TMPDIR
echo "Running on $HOSTNAME. Using $TMPDIR as tmp"


echo
echo "Extracting reads for coverage track"
date
#Coordinates are the positions covered by the read
#NB Doesn't count PCR duplicates or unmapped segments
samtools view -F 1028 ${sampleOutdir}/${name}.${mappedgenome}.bam ${chroms} |
awk -F "\t" 'BEGIN {OFS="\t"} $3!="chrEBV"' |
#sam2bed handles CIGAR string appropriately
sam2bed --do-not-sort |
#reformat as chrom, chromStart, chromEnd, flag
awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $7}' |
#Filter out reads that have no aligned reference bases (e.g. 4S2I23S)
awk -F "\t" 'BEGIN {OFS="\t"} $2<$3' |
sort-bed --max-mem 5G - | tee ${TMPDIR}/${name}.${mappedgenome}.${jobname}.reads.withflag.bed |
awk -F "\t" 'BEGIN {OFS="\t"} !and($4, 512) {print $1, $2, $3}' > ${TMPDIR}/${name}.${mappedgenome}.${jobname}.reads.passed.bed


echo
echo "Making coverage track"
date
#This track excludes PCR duplicates, unmapped segments, and QC fail reads
#BUGBUG does this really exclude QC fail reads? I only see -F 1028 above
#N.B. bedops --chop is only run at locations covered by 1+ read, therefore there will be no entries for bases with 0 coverage
bedops --chop 1 ${TMPDIR}/${name}.${mappedgenome}.${jobname}.reads.passed.bed |
bedmap --delim "\t" --bp-ovr 1 --echo --bases - ${TMPDIR}/${name}.${mappedgenome}.${jobname}.reads.passed.bed |
awk -F "\t" 'BEGIN {OFS="\t"} \
    lastChrom!=$1 || $2!=lastEnd || $4!=lastScore { \
        if(NR>1) { \
            curOutputLine++; \
            print lastChrom, firstStart, lastEnd, ".", lastScore; \
            lastPrinted=NR; \
        } \
        lastChrom=$1; firstStart=$2; lastScore=$4; \
    } \
    { \
        lastStart=$2; lastEnd=$3; \
    } \
END {curOutputLine++; if(NR!=lastPrinted) {print lastChrom, firstStart, lastEnd, ".", lastScore}}' | starch - > ${sampleOutdir}/${name}.${mappedgenome}.${jobname}.coverage.starch


echo
echo "Making allreads coverage track"
date
#This track excludes PCR duplicates and unmapped segments
bedops --chop 1 ${TMPDIR}/${name}.${mappedgenome}.${jobname}.reads.withflag.bed |
bedmap --delim "\t" --bp-ovr 1 --echo --bases - ${TMPDIR}/${name}.${mappedgenome}.${jobname}.reads.withflag.bed |
awk -F "\t" 'BEGIN {OFS="\t"} \
    lastChrom!=$1 || $2!=lastEnd || $4!=lastScore { \
        if(NR>1) { \
            curOutputLine++; \
            print lastChrom, firstStart, lastEnd, ".", lastScore; \
            lastPrinted=NR; \
        } \
        lastChrom=$1; firstStart=$2; lastScore=$4; \
    } \
    { \
        lastStart=$2; lastEnd=$3; \
    } \
END {curOutputLine++; if(NR!=lastPrinted) {print lastChrom, firstStart, lastEnd, ".", lastScore}}' | starch - > ${sampleOutdir}/${name}.${mappedgenome}.${jobname}.coverage.allreads.starch


echo
echo "Making windowed coverage track"
date
#Make windowed coverage track of number of overlapping reads per fixed 100-bp window
#This track excludes PCR duplicates, unmapped segments, and QC fail reads
#BUGBUG does this really exclude QC fail reads? I only see -F 1028 above
for chrom in ${chroms}; do
    awk -v chrom=${chrom} -F "\t" 'BEGIN {OFS="\t"} $1==chrom' ${chromsizes}
done |
awk -F "\t" '{OFS="\t"; print $1, 0, $2}' | sort-bed - | cut -f1,3 | awk -v step=100 -v binwidth=100 'BEGIN {OFS="\t"} {for(i=0; i<=$2-binwidth; i+=step) {print $1, i, i+binwidth, "."} }' |
bedmap --delim "\t" --bp-ovr 1 --echo --bases - ${TMPDIR}/${name}.${mappedgenome}.${jobname}.reads.passed.bed |
#Transform score column into average coverage of base pairs within bin
awk -v binwidth=100 -F "\t" 'BEGIN {OFS="\t"} {$NF=$NF/binwidth; print}' |
starch - > ${sampleOutdir}/${name}.${mappedgenome}.${jobname}.coverage.binned.starch


echo
echo "Calling variants for genome ${mappedgenome} using ploidy ${ploidy} and reference ${referencefasta}. Will annotate dbSNP IDs from ${dbsnpvcf}"
date
#Current documentation at https://samtools.github.io/bcftools/howtos/index.html

minSNPQ=10
minGQ=99
minDP=10

#Set up bcftools call --ploidy argument. ${ploidy} is initialized in genomeinfo.sh based on mappedgenome. For most mamalian genomes, --ploidy is set to depend on sample sex. For cegsvectors_*, it is fixed at 1.
sampleAnnotationSex=`echo "${sampleAnnotation}" | awk -v key="Sex" -F ";" '{for(i=1; i<=NF; i++) { split($i, cur, "="); if(cur[1]==key) {print cur[2]; exit}}}'`
case "${sampleAnnotationSex}" in
Male)
    sex="M";;
Female)
    sex="F";;
*)
    sex="M"
    echo "WARNING assuming sample sex=${sex}"
    ;;
esac

samtools view -H ${sampleOutdir}/${name}.${mappedgenome}.bam | awk -v sex=${sex} -F "\t" 'BEGIN {OFS="\t"} $1=="@RG" {for(i=2; i<=NF; i++) {split($i, tag, ":"); if (tag[1]=="SM") {print tag[2], sex}}}' > $TMPDIR/samplesfile.txt
ploidy="${ploidy} --samples-file $TMPDIR/samplesfile.txt"

#TODO --max-depth 10000 was carried over from 2015 nat genet paper -- still useful? Handling of the latter changed in samtools 1.9
#from Iyer et al PLoS Genet 2018: -C50 -pm2 -F0.05 -d10000
#(explanation of selected options)
#  -C, --adjust-MQ INT     adjust mapping quality; recommended:50, disable:0 [0]
#  -F, --gap-frac FLOAT    minimum fraction of gapped reads [0.002]
#2020mar21 raised max-idepth to permit calling indels from samples with high coverage (i.e. capture)
if [[ "${sampleType}" == "amplicon" ]]; then
    #Not sure why but bcftools throws away regions completely with higher threshold for swift amplicon data
    #https://www.biostars.org/p/384808/
    pileupParams=""
else
    pileupParams="--adjust-MQ 50"
fi

#manually set --ns to avoid dropping reads with QC-fail set
bcftools mpileup -r `echo ${chroms} | perl -pe 's/ /,/g;'` -f ${referencefasta} --ns UNMAP,SECONDARY,DUP --redo-BAQ ${pileupParams} --gap-frac 0.05 --max-depth 10000 --max-idepth 200000 -a DP,AD --output-type u ${sampleOutdir}/${name}.${mappedgenome}.bam |
#NB for some reason if the intermediate file is saved instead of piped, bcftools call outputs a GQ of . for everything
#Iyer et al PLoS Genet 2018 uses --multiallelic-caller
#https://sourceforge.net/p/samtools/mailman/message/32931405/
#https://samtools.github.io/bcftools/call-m.pdf
#TODO would --gvcf merge hzref calls into blocks?
bcftools call --threads $NSLOTS --keep-alts ${ploidy} --multiallelic-caller -f GQ --output-type u |
#Apply only the DP filter now to keep bcf file size down
bcftools filter -i "INFO/DP>=${minDP}" --output-type b > ${sampleOutdir}/${name}.${mappedgenome}.${jobname}.bcf
bcftools index ${sampleOutdir}/${name}.${mappedgenome}.${jobname}.bcf


echo "Filter and normalize variants"
date

#Normalize and split multiallelics
bcftools norm --threads $NSLOTS --check-ref w -m - --fasta-ref ${referencefasta} --output-type u ${sampleOutdir}/${name}.${mappedgenome}.${jobname}.bcf |
#Single ampersand requires all filters to be met in the same sample, FORMAT/DP checks per-sample depth, and --set-GTs masks genotypes failing filters. It doesn't usually matter here since there's only one sample, but it does if this code gets applied to multisample VCF file
bcftools filter -i "QUAL>=${minSNPQ} & GQ>=${minGQ} & FORMAT/DP>=${minDP}" --set-GTs . --output-type u |
#Keep only SNPs with a nonref genotype. --trim-alt-alleles cleans up after --keep-alts above
bcftools view -i 'GT="alt"' --trim-alt-alleles --output-type u - |
#Apply SnpGap/IndelGap after all other filters (particularly --trim-alt-alleles which, in conjunction with bcftools norm -m, can lead to suppression of multiallelic sites otherwise).
#TODO not clear these filters are relevant to certain engineered variants; the main argument for keeping them is that illumina may have poor sensitivity regardless.
bcftools filter --SnpGap 3 --IndelGap 10 --set-GTs . --output-type z > $TMPDIR/${name}.${mappedgenome}.${jobname}.filtered.vcf.gz
bcftools index --tbi $TMPDIR/${name}.${mappedgenome}.${jobname}.filtered.vcf.gz


#Annotate vcf with rsIDs
#Can't be piped since it wants an index
if [[ -f "${dbsnpvcf}" && "${dbsnpvcf}" != "/dev/null" ]]; then
    echo "Adding dbSNP IDs and compressing"
    date
    bcftools annotate -r `echo ${chroms} | perl -pe 's/ /,/g;'` --annotations ${dbsnpvcf} --columns ID --output-type z $TMPDIR/${name}.${mappedgenome}.${jobname}.filtered.vcf.gz -o ${sampleOutdir}/${name}.${mappedgenome}.${jobname}.filtered.vcf.gz
else
    echo "No dbSNP IDs to add -- just compressing"
    date
    bcftools view --output-type z $TMPDIR/${name}.${mappedgenome}.${jobname}.filtered.vcf.gz -o ${sampleOutdir}/${name}.${mappedgenome}.${jobname}.filtered.vcf.gz
fi

echo
echo -e "\nDone!"
date
