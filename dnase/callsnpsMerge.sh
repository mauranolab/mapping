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
mappedgenome=${1}
analysisType=${2}
sampleOutdir=${3}
sampleAnnotation=${4}
src=${5}


source ${src}/genomeinfo.sh ${mappedgenome}


name=`basename ${sampleOutdir}`
echo "Merging ${analysisType} analysis for of sample ${name} against genome ${mappedgenome}"
echo -e "SampleAnnotation\t${sampleAnnotation}"
date



#TMPDIR=`pwd`/tmp.makeTracks.${name}
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
# FIXME: Should the .csi files be removed after merging the byChrom bcf's?
bcftools index ${sampleOutdir}/${name}.${mappedgenome}.bcf
rm -f ${fullbcffiles}


fltvcffiles=`cut -f1 ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt | xargs -I {} echo "${sampleOutdir}/${name}.${mappedgenome}.{}.filtered.vcf.gz"`
fltvcffiles=$(getFilesToMerge ${fltvcffiles})
numfltvcffiles=`echo "${fltvcffiles}" | perl -pe 's/ /\n/g;' | wc -l`
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
echo "Run Delly"
date

## Make a place for the temp bam files with qcfail turned off.
bam_output_qcOK=${TMPDIR}/${name}.${mappedgenome}.QC_OK_bamfile.bam

## Turn off all qcfail flags.
python ${src}/changeFlags.py ${sampleOutdir}/${name}.${mappedgenome}.bam ${bam_output_qcOK}

samtools index ${bam_output_qcOK}

set +e
## Look for variants (DEL INS DUP INV TRA).
delly call -t ALL \
    -o ${sampleOutdir}/${name}.${mappedgenome}.delly.bcf \
    -g ${referencefasta} \
    ${bam_output_qcOK}

# If delly exited with a 0 code then create the vcf
if [ $? -eq 0 ]; then
  bcftools index ${sampleOutdir}/${name}.${mappedgenome}.delly.bcf
  ## Create a vcf file, and set it up for viewing in vcfTabix format.
  ## Only accept variants that pass the FILTER test.
  bcftools filter -i 'FILTER="PASS" & PE>10' ${sampleOutdir}/${name}.${mappedgenome}.delly.bcf > ${sampleOutdir}/${name}.${mappedgenome}.delly.filtered.vcf

  bgzip -f ${sampleOutdir}/${name}.${mappedgenome}.delly.filtered.vcf

  bcftools index ${sampleOutdir}/${name}.${mappedgenome}.delly.filtered.vcf.gz
  tabix -p vcf ${sampleOutdir}/${name}.${mappedgenome}.delly.filtered.vcf.gz
else
  echo "ERROR: Delly Failed"
fi

set -e

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
    #If there is no database ID (e.g. rsid), set up bed ID as chrom : pos _ alt
    cat $TMPDIR/variants.${sampleid}.txt | awk -F "\t" 'BEGIN {OFS="\t"} $4=="." {split($8, gt, "/"); gtout=gt[1]; if(length(gt)>2) {gtout=gtout "/" gt[2]}; chromnum=$1; gsub(/^chr/, "", chromnum); $4=chromnum ":" $2+1 "_" gtout } {print}' | sort-bed - | starch - > ${sampleOutdir}/${name}${vcfsamplename}.${mappedgenome}.genotypes.starch
    
    #NB UCSC link from analysis.sh will be wrong for multisample calling
    #Start from .txt file to simplify logic even though we have to sort again
    #If there is no database ID (e.g. rsid), set up bed ID as chrom : pos _ alt; for database IDs, append alt allele
    #truncgt() is a fancy wrapper around substr for syntactic sugar and to add a "..."
    cat $TMPDIR/variants.${sampleid}.txt | awk -v maxlen=115 -F "\t" 'BEGIN {OFS="\t"} function truncgt(x) {if(length(x)>maxlen) {return substr(x, 1, maxlen) "..." } else {return x} } $4=="." {chromnum=$1; gsub(/^chr/, "", chromnum); $4=chromnum ":" $2+1 } {split($8, gt, "/"); gtout=truncgt(gt[1]); if(length(gt)>2) {gtout=gtout "/" truncgt(gt[2])}; $4=$4 "_" gtout } {print}' |
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
        print $1, $2, $3, $4, $5, ".", $2, $3, color;
    }' |
    sort-bed - > $TMPDIR/${name}${vcfsamplename}.genotypes.ucsc.bed
    
    bedToBigBed -type=bed9 $TMPDIR/${name}${vcfsamplename}.genotypes.ucsc.bed ${chromsizes} ${sampleOutdir}/${name}${vcfsamplename}.${mappedgenome}.genotypes.bb
    
    #Personal Genome SNP format displays two alleles in vertical fashion and provides amino acid changes, but doesn't permit IDs
    #awk -F "\t" 'BEGIN {OFS="\t"} {split($8, gt, "/"); print $1, $2, $3, $8, length(gt), "0,0", "0,0"}'
    #https://genome.ucsc.edu/FAQ/FAQformat.html#format10
done

echo
echo "Parsing Delly VCF track"
date
#No additional filtering, just extract genotypes to starch
#NB repeated from perChrom analysis
minSNPQ=0
#
minTotalDP=0
minAlleleDP=0

${src}/parseSamtoolsGenotypesToBedFiles.pl ${sampleOutdir}/${name}.${mappedgenome}.delly.filtered.vcf.gz $TMPDIR/delly_variants ${minSNPQ} ${minTotalDP} ${minAlleleDP}

nvcfsamples=`bcftools query -l ${sampleOutdir}/${name}.${mappedgenome}.delly.filtered.vcf.gz | wc -l`
#rsids
for sampleid in `bcftools query -l ${sampleOutdir}/${name}.${mappedgenome}.delly.filtered.vcf.gz`; do
    if [[ "${nvcfsamples}" = 1 ]]; then
        vcfsamplename=""
    else
        #for multisample calling, include the sample name in the variants file name
        vcfsamplename=".${sampleid}"
    fi
    #If there is no database ID (e.g. rsid), set up bed ID as chrom : pos _ alt
    cat $TMPDIR/delly_variants.${sampleid}.txt | awk -F "\t" 'BEGIN {OFS="\t"} $4=="." {split($8, gt, "/"); gtout=gt[1]; if(length(gt)>2) {gtout=gtout "/" gt[2]}; chromnum=$1; gsub(/^chr/, "", chromnum); $4=chromnum ":" $2+1 "_" gtout } {print}' | sort-bed - | starch - > ${sampleOutdir}/${name}${vcfsamplename}.${mappedgenome}.delly.genotypes.starch
    
    #NB UCSC link from analysis.sh will be wrong for multisample calling
    #Start from .txt file to simplify logic even though we have to sort again
    #If there is no database ID (e.g. rsid), set up bed ID as chrom : pos _ alt; for database IDs, append alt allele
    #truncgt() is a fancy wrapper around substr for syntactic sugar and to add a "..."
    cat $TMPDIR/delly_variants.${sampleid}.txt | awk -v maxlen=115 -F "\t" 'BEGIN {OFS="\t"} function truncgt(x) {if(length(x)>maxlen) {return substr(x, 1, maxlen) "..." } else {return x} } $4=="." {chromnum=$1; gsub(/^chr/, "", chromnum); $4=chromnum ":" $2+1 } {split($8, gt, "/"); gtout=truncgt(gt[1]); if(length(gt)>2) {gtout=gtout "/" truncgt(gt[2])}; $4=$4 "_" gtout } {print}' |
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
        # Delly has scores set to 10_000, this is outside of the 0->1_000 range
        # and these values need to be clamped
        score=$5;
        if (score>1000) {
          score=1000;
        }
        print $1, $2, $3, $4, score, ".", $2, $3, color;
    }' |
    sort-bed - > $TMPDIR/${name}${vcfsamplename}.delly.genotypes.ucsc.bed
    
    bedToBigBed -type=bed9 $TMPDIR/${name}${vcfsamplename}.delly.genotypes.ucsc.bed ${chromsizes} ${sampleOutdir}/${name}${vcfsamplename}.${mappedgenome}.delly.genotypes.bb
    
    #Personal Genome SNP format displays two alleles in vertical fashion and provides amino acid changes, but doesn't permit IDs
    #awk -F "\t" 'BEGIN {OFS="\t"} {split($8, gt, "/"); print $1, $2, $3, $8, length(gt), "0,0", "0,0"}'
    #https://genome.ucsc.edu/FAQ/FAQformat.html#format10
done


echo
echo -e "\nDone!"
date
