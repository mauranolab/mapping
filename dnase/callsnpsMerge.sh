#!/bin/bash
set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'


###To re-run just the merge
#analysisType=dna
#for f in `ls *BS*/*.bam`; do
#    sampleOutdir=`dirname ${f} | xargs basename`
#    DS=`echo $sampleOutdir | cut -d "-" -f2`
#    name=`basename ${f} | cut -d "." -f1`
#    mappedgenome=`basename ${f} | cut -d "." -f2`
#    echo "${sampleOutdir}"
#    src="${sampleOutdir}/.src"
#    cp -p /vol/mauranolab/mapped/src/dnase/* ${src}
#    
#    qsub -S /bin/bash -cwd -V -terse -j y -b y -o ${sampleOutdir} -N merge.dna.${name}.${mappedgenome} "${src}/dnaMerge.sh ${mappedgenome} ${analysisType} ${name} ${BS} ${src}" #| perl -pe 's/[^\d].+$//g;'
#done


###Parameters
mappedgenome=$1
analysisType=$2
sampleOutdir=$3
sampleAnnotation=$4
src=$5


source ${src}/genomeinfo.sh ${mappedgenome}


name=`basename ${sampleOutdir}`
echo "Merging ${analysisType} analysis for of sample ${name} against genome ${mappedgenome}"
echo -e "SampleAnnotation\t${sampleAnnotation}"
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
            if [ ! -z "${files}" ]; then
                files="${files} "
            fi
            files="${files}${f}"
        fi
    done
    local nfound=`echo ${files} | perl -pe 's/ /\n/g;' | wc -l`
    #NB following line causes funkiness with SLURM output log, overwrites previously printed output
    #echo "Found ${nfound} of ${nexpected} files" > /dev/stderr
    echo "${files}"
}


echo
echo "Merge coverage tracks"
date
#NB: for hg38_full, many jobs will not generate coverage tracks
covfiles=`cat ${sampleOutdir}/inputs.dna.${mappedgenome}.txt | xargs -I {} echo "${sampleOutdir}/${name}.${mappedgenome}.{}.coverage.starch"`
covfiles=$(getFilesToMerge ${covfiles})
starchcat ${covfiles} > ${sampleOutdir}/${name}.${mappedgenome}.coverage.starch
rm -f ${covfiles}

unstarch ${sampleOutdir}/${name}.${mappedgenome}.coverage.starch | cut -f1-3,5 > $TMPDIR/${name}.${mappedgenome}.coverage.bedGraph

#Fix problems with reads running off end of chromosome
bedClip $TMPDIR/${name}.${mappedgenome}.coverage.bedGraph ${chromsizes} $TMPDIR/${name}.${mappedgenome}.coverage.clipped.bedGraph

#Kent tools can't use STDIN
#gets error if run on empty file: needLargeMem: trying to allocate 0 bytes (limit: 100000000000)
bedGraphToBigWig $TMPDIR/${name}.${mappedgenome}.coverage.clipped.bedGraph ${chromsizes} ${sampleOutdir}/${name}.${mappedgenome}.coverage.bw


echo
echo "Merge allreads coverage tracks"
date
#NB: for hg38_full, many jobs will not generate coverage tracks
allreadscovfiles=`cat ${sampleOutdir}/inputs.dna.${mappedgenome}.txt | xargs -I {} echo "${sampleOutdir}/${name}.${mappedgenome}.{}.coverage.allreads.starch"`
allreadscovfiles=$(getFilesToMerge ${allreadscovfiles})
starchcat ${allreadscovfiles} > ${sampleOutdir}/${name}.${mappedgenome}.coverage.allreads.starch
rm -f ${allreadscovfiles}

unstarch ${sampleOutdir}/${name}.${mappedgenome}.coverage.allreads.starch | cut -f1-3,5 > $TMPDIR/${name}.${mappedgenome}.coverage.allreads.bedGraph

#Fix problems with reads running off end of chromosome
bedClip $TMPDIR/${name}.${mappedgenome}.coverage.allreads.bedGraph ${chromsizes} $TMPDIR/${name}.${mappedgenome}.coverage.allreads.clipped.bedGraph

#Kent tools can't use STDIN
#gets error if run on empty file: needLargeMem: trying to allocate 0 bytes (limit: 100000000000)
bedGraphToBigWig $TMPDIR/${name}.${mappedgenome}.coverage.allreads.clipped.bedGraph ${chromsizes} ${sampleOutdir}/${name}.${mappedgenome}.coverage.allreads.bw


#TODO would going through wig instead of bedGraph save space here?
echo
echo "Merge windowed coverage tracks"
date
#NB: for hg38_full, many jobs will not generate coverage tracks
windowedcovfiles=`cat ${sampleOutdir}/inputs.dna.${mappedgenome}.txt | xargs -I {} echo "${sampleOutdir}/${name}.${mappedgenome}.{}.coverage.binned.starch"`
windowedcovfiles=$(getFilesToMerge ${windowedcovfiles})
starchcat ${windowedcovfiles} > ${sampleOutdir}/${name}.${mappedgenome}.coverage.binned.starch
rm -f ${windowedcovfiles}


echo
echo "Merge SNPs"
date
fullvcffiles=`cat ${sampleOutdir}/inputs.dna.${mappedgenome}.txt | xargs -I {} echo "${sampleOutdir}/${name}.${mappedgenome}.{}.vcf.gz"`
fullvcffiles=$(getFilesToMerge ${fullvcffiles})
numfullvcffiles=`echo "${fullvcffiles}" | perl -pe 's/ /\n/g;' | wc -l`
case "${numfullvcffiles}" in
0)
    ;;
1)
    #bcftools concat fails when there is only one file
    cp ${fullvcffiles} ${sampleOutdir}/${name}.${mappedgenome}.vcf.gz;;
*)
    bcftools concat --output-type v ${fullvcffiles} | bgzip -c -@ $NSLOTS > ${sampleOutdir}/${name}.${mappedgenome}.vcf.gz;;
esac
bcftools index ${sampleOutdir}/${name}.${mappedgenome}.vcf.gz
tabix -p vcf ${sampleOutdir}/${name}.${mappedgenome}.vcf.gz
rm -f ${fullvcffiles}


fltvcffiles=`cat ${sampleOutdir}/inputs.dna.${mappedgenome}.txt | xargs -I {} echo "${sampleOutdir}/${name}.${mappedgenome}.{}.filtered.vcf.gz"`
fltvcffiles=$(getFilesToMerge ${fltvcffiles})
numfltvcffiles=`echo "${fullvcffiles}" | perl -pe 's/ /\n/g;' | wc -l`
case "${numfltvcffiles}" in
0)
    ;;
1)
    #bcftools concat fails when there is only one file
    cp ${fltvcffiles} ${sampleOutdir}/${name}.${mappedgenome}.filtered.vcf.gz;;
*)
    bcftools concat --output-type v ${fltvcffiles} | bgzip -c -@ $NSLOTS > ${sampleOutdir}/${name}.${mappedgenome}.filtered.vcf.gz;;
esac

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
    cat $TMPDIR/variants.${sampleid}.txt | awk -F "\t" 'BEGIN {OFS="\t"; num=0} $4=="." {num++; $4="var" num} {print}' | sort-bed - | starch - > ${sampleOutdir}/${name}${vcfsamplename}.${mappedgenome}.genotypes.starch
    
    #NB UCSC link from analysis.sh will be wrong for multisample calling
    unstarch ${sampleOutdir}/${name}${vcfsamplename}.${mappedgenome}.genotypes.starch | awk -v maxlen=115 -F "\t" 'BEGIN {OFS="\t"} function truncgt(x) {if(length(x)>maxlen) {return substr(x, 1, maxlen) "..." } else {return x} } {split($8, gt, "/"); gtout=truncgt(gt[1]); if(length(gt)>1) {gtout=gtout "/" truncgt(gt[2])} print $1, $2, $3, $4 "_" gtout, $5}' > $TMPDIR/${name}${vcfsamplename}.genotypes.ucsc.bed
    bedToBigBed -type=bed5 $TMPDIR/${name}${vcfsamplename}.genotypes.ucsc.bed ${chromsizes} ${sampleOutdir}/${name}${vcfsamplename}.${mappedgenome}.genotypes.bb
    
    #Personal Genome SNP format displays two alleles in vertical fashion and provides amino acid changes, but doesn't permit IDs
    #awk -F "\t" 'BEGIN {OFS="\t"} {split($8, gt, "/"); print $1, $2, $3, $8, length(gt), "0,0", "0,0"}'
    #https://genome.ucsc.edu/FAQ/FAQformat.html#format10
done


echo
echo -e "\nDone!"
date
