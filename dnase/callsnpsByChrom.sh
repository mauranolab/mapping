#!/bin/bash
set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header'
#alias starch='starch --header'
alias closest-features='closest-features --header'


###Parameters
mappedgenome=$1
analysisType=$2
name=$3
src=$4

source ${src}/genomeinfo.sh ${mappedgenome}

jobid=$SGE_TASK_ID
#jobid=10

sampleOutdir=${name}
chrom=`cat ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt | awk -v jobid=$jobid 'NR==jobid'`
echo "Running ${analysisType} analysis for ${chrom} of sample ${name} against genome ${mappedgenome}"
date


#TMPDIR=`pwd`/tmp.makeTracks.${name}
#mkdir -p $TMPDIR
echo "using $TMPDIR as TMPDIR"


echo
echo "Making coverage track"
date
#Coordinates are the positions covered by the read
#NB includes reads marked PCR duplicates
unstarch ${chrom} ${sampleOutdir}/${name}.${mappedgenome}.reads.starch |
awk -F "\t" 'BEGIN {OFS="\t"} NR>1 {readlength=$5; if($6=="+") {chromStart=$2} else {chromStart=$2-readlength+1} print $1, chromStart, chromStart+readlength}' |
sort-bed --max-mem 5G - | 
starch - > ${TMPDIR}/${name}.${mappedgenome}.${chrom}.reads.withlength.starch

unstarch ${TMPDIR}/${name}.${mappedgenome}.${chrom}.reads.withlength.starch | bedops --chop 1 - |
bedmap --delim "\t" --bp-ovr 1 --echo --count - ${TMPDIR}/${name}.${mappedgenome}.${chrom}.reads.withlength.starch |
awk -F "\t" 'BEGIN {OFS="\t"} \
    lastChrom!=$1 || $2!=lastEnd || $4!=lastScore { \
        if(NR>1) { \
            curOutputLine++; \
            print lastChrom, firstStart, lastEnd, lastScore; \
            lastPrinted=NR; \
        } \
        lastChrom=$1; firstStart=$2; lastScore=$4; \
    } \
    { \
        lastStart=$2; lastEnd=$3; \
    } \
END {curOutputLine++; if(NR!=lastPrinted) {print lastChrom, firstStart, lastEnd, lastScore}}' > $TMPDIR/${name}.${mappedgenome}.${chrom}.coverage.bedGraph

cat $TMPDIR/${name}.${mappedgenome}.${chrom}.coverage.bedGraph | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' | starch - > ${sampleOutdir}/${name}.${mappedgenome}.${chrom}.coverage.starch


echo
echo "Calling variants for genome ${mappedgenome} using ploidy ${ploidy} and reference ${referencefasta}. Will annotate dbSNP IDs from $dbsnpvcf"
date

#Current documentation at https://samtools.github.io/bcftools/howtos/index.html

minSNPQ=200
minGQ=99
#hack to get VCF at minDP8, but starch at 12
minDP=20


samtools mpileup -r ${chrom} -u --min-BQ 20 --max-depth 10000 --redo-BAQ -f ${referencefasta} --BCF -t DP,AD ${sampleOutdir}/${name}.${mappedgenome}.bam |
#NB for some reason if the intermediate file is saved instead of piped, bcftools call outputs a GQ of . for everything
#TODO should use --multiallelic-caller ?
bcftools call ${ploidy} --keep-alts --consensus-caller --variants-only -f GQ --output-type v | bgzip -c -@ $NSLOTS > ${sampleOutdir}/${name}.${mappedgenome}.${chrom}.vcf.gz
bcftools index ${sampleOutdir}/${name}.${mappedgenome}.${chrom}.vcf.gz
tabix -p vcf ${sampleOutdir}/${name}.${mappedgenome}.${chrom}.vcf.gz


echo "Filter and normalize variants"
date
bcftools filter -O u -i "INFO/DP>=${minDP} && QUAL>=${minSNPQ} && GQ>=${minGQ}" ${sampleOutdir}/${name}.${mappedgenome}.${chrom}.vcf.gz |
bcftools norm --threads $NSLOTS --check-ref w --fasta-ref ${referencefasta} --output-type z - > $TMPDIR/${name}.${mappedgenome}.${chrom}.filtered.vcf.gz
bcftools index $TMPDIR/${name}.${mappedgenome}.${chrom}.filtered.vcf.gz

#Annotate vcf with rsIDs
#Can't be piped since it wants an index
if [[ -f "${dbsnpvcf}" && "${dbsnpvcf}" != "/dev/null" ]]; then
    echo "Adding dbSNP IDs and compressing"
    date
    bcftools annotate -r ${chrom} --annotations ${dbsnpvcf} --columns ID --output-type v $TMPDIR/${name}.${mappedgenome}.${chrom}.filtered.vcf.gz | 
    bgzip -c -@ $NSLOTS > ${sampleOutdir}/${name}.${mappedgenome}.${chrom}.filtered.vcf.gz
else
    echo "No dbSNP IDs to add -- just compressing"
    date
    bcftools view --output-type v $TMPDIR/${name}.${mappedgenome}.${chrom}.filtered.vcf.gz | 
    bgzip -c -@ $NSLOTS > ${sampleOutdir}/${name}.${mappedgenome}.${chrom}.filtered.vcf.gz
fi

rm -f ${sampleOutdir}/${name}.${mappedgenome}.${chrom}.vcf.gz.csi ${sampleOutdir}/${name}.${mappedgenome}.${chrom}.vcf.gz.tbi

echo
echo -e "\nDone!"
date
