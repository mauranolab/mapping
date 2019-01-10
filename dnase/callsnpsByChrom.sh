#!/bin/bash
set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
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
#NB Doesn't count PCR duplicates
unstarch ${chrom} ${sampleOutdir}/${name}.${mappedgenome}.reads.starch |
awk -F "\t" 'BEGIN {OFS="\t"} !and($7, 1024)  {readlength=$5; if($6=="+") {chromStart=$2} else {chromStart=$2-readlength+1} print $1, chromStart, chromStart+readlength}' |
sort-bed --max-mem 5G - | 
starch - > ${TMPDIR}/${name}.${mappedgenome}.${chrom}.reads.withlength.starch

unstarch ${chrom} ${TMPDIR}/${name}.${mappedgenome}.${chrom}.reads.withlength.starch | bedops --chop 1 - |
bedmap --chrom ${chrom} --delim "\t" --bp-ovr 1 --echo --bases - ${TMPDIR}/${name}.${mappedgenome}.${chrom}.reads.withlength.starch |
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
END {curOutputLine++; if(NR!=lastPrinted) {print lastChrom, firstStart, lastEnd, ".", lastScore}}' | starch - > ${sampleOutdir}/${name}.${mappedgenome}.${chrom}.coverage.starch


echo
echo "Making windowed coverage track"
date
#Make windowed coverage track of number of overlapping reads per fixed 100-bp window
awk -v chrom=${chrom} -F "\t" 'BEGIN {OFS="\t"} $1==chrom' ${chromsizes} | 
awk '{OFS="\t"; $3=$2; $2=0; print}' | sort-bed - | cut -f1,3 | awk -v step=100 -v binwidth=100 'BEGIN {OFS="\t"} {for(i=0; i<=$2-binwidth; i+=step) {print $1, i, i+binwidth, "."} }' | 
bedmap --chrom ${chrom} --delim "\t" --bp-ovr 1 --echo --bases -  ${TMPDIR}/${name}.${mappedgenome}.${chrom}.reads.withlength.starch |
#Transform score column into average coverage of base pairs within bin
awk -v binwidth=100 -F "\t" 'BEGIN {OFS="\t"} {$NF=$NF/binwidth; print}' |
starch - > ${sampleOutdir}/${name}.${mappedgenome}.${chrom}.coverage.binned.starch


echo
echo "Calling variants for genome ${mappedgenome} using ploidy ${ploidy} and reference ${referencefasta}. Will annotate dbSNP IDs from ${dbsnpvcf}"
date

#Current documentation at https://samtools.github.io/bcftools/howtos/index.html
echo "TODO WARNING assuming male sample"
samtools view -H ${sampleOutdir}/${name}.${mappedgenome}.bam | awk -F "\t" 'BEGIN {OFS="\t"} $1=="@RG" {for(i=2; i<=NF; i++) {split($i, tag, ":"); if (tag[1]=="SM") {print tag[2], "M"}}}' > $TMPDIR/samplesfile.txt
ploidy="${ploidy} --samples-file $TMPDIR/samplesfile.txt"


#TODO --min-BQ 20 --max-depth 10000 were carried over from 2015 nat genet paper -- still useful? Handling of the latter changed in samtools 1.9
#from Iyer et al PLoS Genet 2018
#-C50 -pm2 -F0.05 -d10000
#  -C, --adjust-MQ INT     adjust mapping quality; recommended:50, disable:0 [0]
#  -F, --gap-frac FLOAT    minimum fraction of gapped reads [0.002]
bcftools mpileup -r ${chrom} --redo-BAQ -f ${referencefasta} --adjust-MQ 50 --gap-frac 0.05 --max-depth 10000 -a DP,AD -O u ${sampleOutdir}/${name}.${mappedgenome}.bam |
#NB for some reason if the intermediate file is saved instead of piped, bcftools call outputs a GQ of . for everything
#Iyer et al PLoS Genet 2018 uses --multiallelic-caller
#https://sourceforge.net/p/samtools/mailman/message/32931405/
#https://samtools.github.io/bcftools/call-m.pdf
bcftools call ${ploidy} --keep-alts --consensus-caller --variants-only -f GQ --output-type v | bgzip -c -@ $NSLOTS > ${sampleOutdir}/${name}.${mappedgenome}.${chrom}.vcf.gz
bcftools index ${sampleOutdir}/${name}.${mappedgenome}.${chrom}.vcf.gz
tabix -p vcf ${sampleOutdir}/${name}.${mappedgenome}.${chrom}.vcf.gz


echo "Filter and normalize variants"
date

minSNPQ=10
minGQ=99
minDP=10

bcftools filter -i "INFO/DP>=${minDP} && QUAL>=${minSNPQ} && GQ>=${minGQ}" --SnpGap 3 --IndelGap 10 -O u ${sampleOutdir}/${name}.${mappedgenome}.${chrom}.vcf.gz |
#Iyer et al PLoS Genet 2018 uses -m -any (split multiallelics)
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
