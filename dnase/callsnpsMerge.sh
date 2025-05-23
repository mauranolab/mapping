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
#    cp -p /gpfs/data/mauranolab/mapped/src/dnase/* ${src}
#    
#    qsub -S /bin/bash -cwd -V -terse -j y -b y -o ${sampleOutdir} -N merge.dna.${name}.${mappedgenome} "${src}/dnaMerge.sh ${mappedgenome} ${analysisType} ${name} ${BS} ${src}" #| perl -pe 's/[^\d].+$//g;'
#done


###Parameters
mappedgenome=${1}
analysisType=${2}
sampleOutdir=${3}
sampleAnnotation=${4}
src=${5}


source ${src}/genomeinfo.sh ${mappedgenome}


echo "Running on $HOSTNAME. Using $TMPDIR as tmp"

name=`basename ${sampleOutdir}`
echo "Merging ${analysisType} analysis for of sample ${name} against genome ${mappedgenome}"
echo -e "SampleAnnotation\t${sampleAnnotation}"
date



#TMPDIR=`pwd`/tmp.callsnpsMerge.${name}
#mkdir -p $TMPDIR
echo "Running on $HOSTNAME. Using $TMPDIR as tmp"
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
covfiles=`cut -f1 ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt | xargs -I {} echo "${sampleOutdir}/${name}.${mappedgenome}.{}.coverage.starch"`
covfiles=$(getFilesToMerge ${covfiles})
starchcat ${covfiles} > ${sampleOutdir}/${name}.${mappedgenome}.coverage.starch
rm -f ${covfiles}

unstarch ${sampleOutdir}/${name}.${mappedgenome}.coverage.starch | cut -f1-3,5 > $TMPDIR/${name}.${mappedgenome}.coverage.bedGraph

if [ ! -s "$TMPDIR/${name}.${mappedgenome}.coverage.bedGraph" ]; then
    # The coverage.bedGraph file is empty
    bedgraph_chrom=$(head -n 1 ${chromsizes} | cut -f1)
    echo -e "${bedgraph_chrom}\t0\t0\t0" > $TMPDIR/${name}.${mappedgenome}.coverage.clipped.bedGraph
else
    #Fix problems with reads running off end of chromosome
    bedClip $TMPDIR/${name}.${mappedgenome}.coverage.bedGraph ${chromsizes} $TMPDIR/${name}.${mappedgenome}.coverage.clipped.bedGraph
fi

#Kent tools can't use STDIN
#gets error if run on empty file: needLargeMem: trying to allocate 0 bytes (limit: 100000000000)
bedGraphToBigWig $TMPDIR/${name}.${mappedgenome}.coverage.clipped.bedGraph ${chromsizes} ${sampleOutdir}/${name}.${mappedgenome}.coverage.bw


echo
echo "Merge allreads coverage tracks"
date
#NB: for hg38_full, many jobs will not generate coverage tracks
allreadscovfiles=`cut -f1 ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt | xargs -I {} echo "${sampleOutdir}/${name}.${mappedgenome}.{}.coverage.allreads.starch"`
allreadscovfiles=$(getFilesToMerge ${allreadscovfiles})
starchcat ${allreadscovfiles} > ${sampleOutdir}/${name}.${mappedgenome}.coverage.allreads.starch
rm -f ${allreadscovfiles}

unstarch ${sampleOutdir}/${name}.${mappedgenome}.coverage.allreads.starch | cut -f1-3,5 > $TMPDIR/${name}.${mappedgenome}.coverage.allreads.bedGraph

if [ ! -s "$TMPDIR/${name}.${mappedgenome}.coverage.allreads.bedGraph" ]; then
    # The coverage.allreads.bedGraph file is empty
    bedgraph_chrom=$(head -n 1 ${chromsizes} | cut -f1)
    echo -e "${bedgraph_chrom}\t0\t0\t0" > $TMPDIR/${name}.${mappedgenome}.coverage.allreads.clipped.bedGraph
else
    #Fix problems with reads running off end of chromosome
    bedClip $TMPDIR/${name}.${mappedgenome}.coverage.allreads.bedGraph ${chromsizes} $TMPDIR/${name}.${mappedgenome}.coverage.allreads.clipped.bedGraph
fi

#Kent tools can't use STDIN
#gets error if run on empty file: needLargeMem: trying to allocate 0 bytes (limit: 100000000000)
bedGraphToBigWig $TMPDIR/${name}.${mappedgenome}.coverage.allreads.clipped.bedGraph ${chromsizes} ${sampleOutdir}/${name}.${mappedgenome}.coverage.allreads.bw


#TODO would going through wig instead of bedGraph save space here?
echo
echo "Merge windowed coverage tracks"
date
#NB: for hg38_full, many jobs will not generate coverage tracks
windowedcovfiles=`cut -f1 ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt | xargs -I {} echo "${sampleOutdir}/${name}.${mappedgenome}.{}.coverage.binned.starch"`
windowedcovfiles=$(getFilesToMerge ${windowedcovfiles})
starchcat ${windowedcovfiles} > ${sampleOutdir}/${name}.${mappedgenome}.coverage.binned.starch
rm -f ${windowedcovfiles}


echo
echo "Merge SNPs"
date
fullbcffiles=`cut -f1 ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt | xargs -I {} echo "${sampleOutdir}/${name}.${mappedgenome}.{}.bcf"`
fullbcffiles=$(getFilesToMerge ${fullbcffiles})
numfullbcffiles=`echo "${fullbcffiles}" | perl -pe 's/ /\n/g;' | wc -l`
case "${numfullbcffiles}" in
0)
    ;;
1)
    #bcftools concat fails when there is only one file
    cp ${fullbcffiles} ${sampleOutdir}/${name}.${mappedgenome}.bcf;;
*)
    bcftools concat --threads $NSLOTS --output-type b ${fullbcffiles} > ${sampleOutdir}/${name}.${mappedgenome}.bcf;;
esac
bcftools index ${sampleOutdir}/${name}.${mappedgenome}.bcf
rm -f ${fullbcffiles}

fullbcfcsifiles=`cut -f1 ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt | xargs -I {} echo "${sampleOutdir}/${name}.${mappedgenome}.{}.bcf.csi"`
fullbcfcsifiles=$(getFilesToMerge ${fullbcfcsifiles})
rm -f ${fullbcfcsifiles}


fltvcffiles=`cut -f1 ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt | xargs -I {} echo "${sampleOutdir}/${name}.${mappedgenome}.{}.filtered.vcf.gz"`
fltvcffiles=$(getFilesToMerge ${fltvcffiles})
numfltvcffiles=`echo "${fltvcffiles}" | perl -pe 's/ /\n/g;' | wc -l`
case "${numfltvcffiles}" in
0)
    ;;
1)
    #bcftools concat fails when there is only one file
    cp ${fltvcffiles} ${TMPDIR}/${name}.${mappedgenome}.filtered.vcf.gz;;
*)
    bcftools concat --output-type z ${fltvcffiles} -o ${TMPDIR}/${name}.${mappedgenome}.filtered.vcf.gz;;
esac
bcftools index --tbi ${TMPDIR}/${name}.${mappedgenome}.filtered.vcf.gz
rm -f ${fltvcffiles}


# DELLY does not have a --regions flag to allow for specific regions to be matched.
# this makes it difficult to parallelize by chromosome like the bcftools call command
echo
echo "Running DELLY"
date

## DELLY is hardcoded to skip reads with qcfail flag set, so generate new bam file with 512 cleared
${src}/changeSAMflags.py --clearflag 512 ${sampleOutdir}/${name}.${mappedgenome}.bam ${TMPDIR}/${name}.${mappedgenome}.QC_OK_bamfile.bam
samtools index ${TMPDIR}/${name}.${mappedgenome}.QC_OK_bamfile.bam


# If there are not enough valid reads, delly returns a non-zero exit code
#${dellyexclude} is set by genomeinfo.sh to exclude telomeric/centromeric regions (presumably more important for WGS)
set +e
delly call -t ALL -o ${sampleOutdir}/${name}.${mappedgenome}.delly.bcf ${dellyexclude} -g ${referencefasta} ${TMPDIR}/${name}.${mappedgenome}.QC_OK_bamfile.bam
dellyExitCode=$?
set -e

# If DELLY exited with a 0 code then create the vcf
if [ "$dellyExitCode" -eq 0 ]; then
    echo "Merging DELLY variants in to .filtered.vcf file"
    # Only accept variants that pass the FILTER test.
    # Don't bother masking excluded genotypes and just straight hard drop them
    bcftools filter --output-type z -i 'FILTER="PASS" & PE>10 & FORMAT/DV / (FORMAT/DV + FORMAT/DR) >= 0.20' ${sampleOutdir}/${name}.${mappedgenome}.delly.bcf -o ${sampleOutdir}/${name}.${mappedgenome}.DELLY.filtered.vcf.gz
    bcftools index --tbi ${sampleOutdir}/${name}.${mappedgenome}.DELLY.filtered.vcf.gz
    
    echo "Merging DELLY VCF"
    # Merge the DELLY calls in with the combined bcftools calls
    bcftools concat -a -O z -o ${sampleOutdir}/${name}.${mappedgenome}.filtered.vcf.gz ${sampleOutdir}/${name}.${mappedgenome}.DELLY.filtered.vcf.gz ${TMPDIR}/${name}.${mappedgenome}.filtered.vcf.gz
else
    #Use bcftools vcf directly
    cp --preserve=timestamps ${TMPDIR}/${name}.${mappedgenome}.filtered.vcf.gz ${sampleOutdir}/${name}.${mappedgenome}.filtered.vcf.gz
    echo "WARNING: DELLY Failed"
fi
bcftools index --tbi ${sampleOutdir}/${name}.${mappedgenome}.filtered.vcf.gz
#BUGBUG naming of "DELLY.filtered.vcf.gz" and "DELLY.filtered.vcf.gz.tbi" does not match ".delly.bcf" and ".delly.bcf.csi"


echo
echo "Parsing VCF track"
date
#NB does not include SVs from DELLY to ease visualization, especially for hybrid mouse samples
#No additional filtering, just extract genotypes to starch
#NB repeated from perChrom analysis
minSNPQ=0
#
minTotalDP=0
minAlleleDP=0

${src}/parseSamtoolsGenotypesToBedFiles.pl ${TMPDIR}/${name}.${mappedgenome}.filtered.vcf.gz  $TMPDIR/variants ${minSNPQ} ${minTotalDP} ${minAlleleDP}

nvcfsamples=`bcftools query -l ${sampleOutdir}/${name}.${mappedgenome}.filtered.vcf.gz | wc -l`
#rsids
for sampleid in `bcftools query -l ${sampleOutdir}/${name}.${mappedgenome}.filtered.vcf.gz`; do
    if [[ "${nvcfsamples}" = 1 ]]; then
        vcfsamplename=""
    else
        #for multisample calling, include the sample name in the variants file name
        vcfsamplename=".${sampleid}"
    fi
    #If there is no database ID (e.g. rsid), set up bed ID as chrom : pos _ alt (multiple alt alleles separated by "/").
    cat $TMPDIR/variants.${sampleid}.txt | awk -F "\t" 'BEGIN {OFS="\t"} $4=="." {split($8, gt, "/"); gtout=gt[1]; if(length(gt)>2) {gtout=gtout "/" gt[2]}; chromnum=$1; gsub(/^chr/, "", chromnum); $4=chromnum ":" $2+1 "_" gtout } {print}' | sort-bed - | starch - > ${sampleOutdir}/${name}${vcfsamplename}.${mappedgenome}.genotypes.starch
    
    #NB UCSC link from analysis.sh will be wrong for multisample calling
    #Start from .txt file to simplify logic even though we have to sort again
    #If there is no database ID (e.g. rsid), set up bed ID as chrom : pos _ alt (multiple alt alleles separated by "/"); each alt allele is truncated to max of 115 chars.
    #truncgt() is a fancy wrapper around substr for syntactic sugar and to add a "..."
    cat $TMPDIR/variants.${sampleid}.txt | awk -v maxlen=115 -F "\t" 'BEGIN {OFS="\t"} function truncgt(x) {if(length(x)>maxlen) {return substr(x, 1, maxlen) "..." } else {return x} } $4=="." {split($8, gt, "/"); gtout=truncgt(gt[1]); if(length(gt)>2) {gtout=gtout "/" truncgt(gt[2])}; chromnum=$1; gsub(/^chr/, "", chromnum); $4=chromnum ":" $2+1 "_" gtout } {print}' |
    awk -F "\t" 'BEGIN {OFS="\t"} {
        ref=$6;
        numAlleles=split($8, gt, "/");
        # Assign colors to the genotypes.bb lines.
        if (numAlleles>1 && gt[1]!=gt[2]) {
            color="19,165,220";        # het - blue
        } else {
            if (gt[1] == ref) {
                color="105,105,105";   # hz_ref - grey
            } else {
                color="253,199,0";     # hz_nonref - yellow
            }
        }
        # Enforce UCSC score max of 1000 (DELLY in particular has scores set to 10_000)
        score=$5;
        if (score>1000) {
            score=1000;
        }
        print $1, $2, $3, $4, score, ".", $2, $3, color;
    }' |
    sort-bed - > $TMPDIR/${name}${vcfsamplename}.genotypes.ucsc.bed
    
    bedToBigBed -type=bed9 $TMPDIR/${name}${vcfsamplename}.genotypes.ucsc.bed ${chromsizes} ${sampleOutdir}/${name}${vcfsamplename}.${mappedgenome}.genotypes.bb
    
    #Personal Genome SNP format displays two alleles in vertical fashion and provides amino acid changes, but doesn't permit IDs
    #awk -F "\t" 'BEGIN {OFS="\t"} {split($8, gt, "/"); print $1, $2, $3, $8, length(gt), "0,0", "0,0"}'
    #https://genome.ucsc.edu/FAQ/FAQformat.html#format10
done


echo
echo -e "\nDone!"
date
