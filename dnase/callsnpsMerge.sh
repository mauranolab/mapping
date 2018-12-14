#!/bin/bash
set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'


###To re-run just the merge
#analysisType=callsnps
#for f in `ls *BS*/*.bam`; do
#    sampleOutdir=`dirname ${f} | xargs basename`
#    DS=`echo $sampleOutdir | cut -d "-" -f2`
#    name=`basename ${f} | cut -d "." -f1`
#    mappedgenome=`basename ${f} | cut -d "." -f2`
#    echo "${sampleOutdir}"
#    src="${sampleOutdir}/.src"
#    cp -p /vol/mauranolab/mapped/src/* ${src}
#    
#    qsub -S /bin/bash -cwd -V -terse -j y -b y -o ${sampleOutdir} -N merge.callsnps.${name}.${mappedgenome} "${src}/callsnpsMerge.sh ${mappedgenome} ${analysisType} ${name} ${BS} ${src}" #| perl -pe 's/[^\d].+$//g;'
#done


###Parameters
mappedgenome=$1
analysisType=$2
name=$3
src=$4


source ${src}/genomeinfo.sh ${mappedgenome}


sampleOutdir=${name}
echo "Merging ${analysisType} analysis for of sample ${name} against genome ${mappedgenome}"
date



#TMPDIR=`pwd`/tmp.makeTracks.${name}
#mkdir -p $TMPDIR
echo "using $TMPDIR as TMPDIR"
date


getFilesToMerge()
{
    local nexpected=`echo $* | perl -pe 's/ /\n/g;' | wc -l`
    local files=""
    for f in $*; do
        if [ -f "${f}" ]; then
            files="${files} ${f}"
        fi
    done
    local nfound=`echo ${files} | perl -pe 's/ /\n/g;' | wc -l`
    echo "Found ${nfound} of ${nexpected} files" > /dev/stderr
    echo "${files}"
}


echo "Merge coverage tracks"
#NB: for hg38_full, many jobs will not generate coverage tracks
covfiles=`cat ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt | xargs -I {} echo "${sampleOutdir}/${name}.${mappedgenome}.{}.coverage.starch"`
covfiles=$(getFilesToMerge $covfiles)
starchcat ${covfiles} > ${sampleOutdir}/${name}.${mappedgenome}.coverage.starch
rm -f ${covfiles}

unstarch ${sampleOutdir}/${name}.${mappedgenome}.coverage.starch | cut -f1-3,5 > $TMPDIR/${name}.${mappedgenome}.coverage.bedGraph

#Fix problems with reads running off end of chromosome
bedClip $TMPDIR/${name}.${mappedgenome}.coverage.bedGraph ${chromsizes} $TMPDIR/${name}.${mappedgenome}.coverage.clipped.bedGraph

#Kent tools can't use STDIN
#gets error if run on empty file: needLargeMem: trying to allocate 0 bytes (limit: 100000000000)
bedGraphToBigWig $TMPDIR/${name}.${mappedgenome}.coverage.clipped.bedGraph ${chromsizes} ${sampleOutdir}/${name}.${mappedgenome}.coverage.bw


echo "Merge SNPs"
date
fullvcffiles=`cat ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt | xargs -I {} echo "${sampleOutdir}/${name}.${mappedgenome}.{}.vcf.gz"`
fullvcffiles=$(getFilesToMerge $fullvcffiles)
bcftools concat  --output-type v ${fullvcffiles} | bgzip -c -@ $NSLOTS > ${sampleOutdir}/${name}.${mappedgenome}.vcf.gz
bcftools index ${sampleOutdir}/${name}.${mappedgenome}.vcf.gz
tabix -p vcf ${sampleOutdir}/${name}.${mappedgenome}.vcf.gz
rm -f ${fullvcffiles}


fltvcffiles=`cat ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt | xargs -I {} echo "${sampleOutdir}/${name}.${mappedgenome}.{}.filtered.vcf.gz"`
fltvcffiles=$(getFilesToMerge $fltvcffiles)
bcftools concat  --output-type v ${fltvcffiles} | bgzip -c -@ $NSLOTS > ${sampleOutdir}/${name}.${mappedgenome}.filtered.vcf.gz
bcftools index ${sampleOutdir}/${name}.${mappedgenome}.filtered.vcf.gz
tabix -p vcf ${sampleOutdir}/${name}.${mappedgenome}.filtered.vcf.gz
rm -f ${fltvcffiles}


echo
echo "Parsing VCF track"
date
#No additional filtering, just extract genotypes to starch
#NB repeated from perChrom analysis
minSNPQ=0
#
minTotalDP=0
minAlleleDP=0

${src}/parseSamtoolsGenotypesToBedFiles.pl ${sampleOutdir}/${name}.${mappedgenome}.filtered.vcf.gz $TMPDIR/variants ${minSNPQ} ${minTotalDP} ${minAlleleDP}

nvcfsamples=`bcftools query -l ${sampleOutdir}/${name}.${mappedgenome}.filtered.vcf.gz | wc -l`
#rsids
for sampleid in `bcftools query -l ${sampleOutdir}/${name}.${mappedgenome}.filtered.vcf.gz`; do
    if [[ "${nvcfsamples}" = 1 ]]; then
        vcfsamplename=""
    else
        #for multisample calling, include the sample name in the variants file name
        vcfsamplename=".${sampleid}"
    fi
    #BUGBUG labels indels with snv prefix
    cat $TMPDIR/variants.${sampleid}.txt | awk -F "\t" 'BEGIN {OFS="\t"; num=0} $4=="." {num++; $4="snv" num} {print}' | sort-bed - | starch - > ${sampleOutdir}/${name}${vcfsamplename}.${mappedgenome}.variants.starch
    
    #NB UCSC link from analysis.sh will be wrong for multisample calling
    #TODO no ... if truncated
    unstarch ${sampleOutdir}/${name}${vcfsamplename}.${mappedgenome}.variants.starch | awk -v maxlen=115 -F "\t" 'BEGIN {OFS="\t"} function truncgt(x) {if(length(x)>maxlen) {return substr(x, 1, maxlen) "..." } else {return x} } {split($8, gt, "/"); print $1, $2, $3, $4 "_" truncgt(gt[1]) "_" truncgt(gt[2]), $5}' > $TMPDIR/${name}${vcfsamplename}.variants.ucsc.bed

    bedToBigBed -type=bed5 $TMPDIR/${name}${vcfsamplename}.variants.ucsc.bed ${chromsizes} ${sampleOutdir}/${name}${vcfsamplename}.${mappedgenome}.variants.bb
done


echo
echo -e "\nDone!"
date
