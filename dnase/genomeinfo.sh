#!/bin/bash
mappedgenome=$1

case "${mappedgenome}" in
#Reference sequence can be .fa.gz as long as bgzip was used to compress
#more info: bcftools call --ploidy list ?
hg19)
    #NB not the analysis set
    bwaIndex=/vol/isg/annotation/bwaIndex/hg19all/hg19all
    ploidy="--ploidy GRCh37"
    referencefasta=/vol/isg/annotation/fasta/hg19/hg19.fa
    dbsnpvcf=/dev/null
    ;;
hg38_noalt)
    bwaIndex=/vol/isg/annotation/bwaIndex/hg38_noalt/hg38_noalt
    ploidy="--ploidy GRCh38"
    referencefasta=/vol/isg/annotation/fasta/hg38_analysisSet_noalt/hg38_analysisSet_noalt.fa.gz
    dbsnpvcf=/vol/isg/annotation/bed/hg38/snp151/All_20180418.vcf.gz
    ;;
hg38_full)
    bwaIndex=/vol/isg/annotation/bwaIndex/hg38_full/hg38_full
    ploidy="--ploidy GRCh38"
    referencefasta=/vol/isg/annotation/fasta/hg38_analysisSet_full/hg38_analysisSet_full.fa.gz
    dbsnpvcf=/vol/isg/annotation/bed/hg38/snp151/All_20180418.vcf.gz
    ;;
mm10)
    bwaIndex=/vol/isg/annotation/bwaIndex/mm10_no_alt_analysis_set/mm10_no_alt_analysis_set
    ploidy=""
    referencefasta=/vol/isg/annotation/fasta/mm10_no_alt_analysis_set/mm10_no_alt_analysis_set.fa.gz
    dbsnpvcf=/dev/null
    ;;
rn6)
    bwaIndex=/vol/isg/annotation/bwaIndex/rn6/rn6
    ploidy=""
    referencefasta=/vol/isg/annotation/fasta/rn6/rn6.fa.gz
    dbsnpvcf=/dev/null
    ;;
hg38_sacCer3)
    bwaIndex=/vol/isg/annotation/bwaIndex/hg38_sacCer3/hg38_sacCer3
    ploidy="--ploidy GRCh38"
    referencefasta=/vol/isg/annotation/fasta/hg38_sacCer3/hg38_sacCer3.fa.gz
    dbsnpvcf=/vol/isg/annotation/bed/hg38/snp151/All_20180418.vcf.gz
    ;;
mm10_sacCer3)
    bwaIndex=/vol/isg/annotation/bwaIndex/mm10_sacCer3/mm10_sacCer3
    ploidy=""
    referencefasta=/vol/isg/annotation/fasta/mm10_sacCer3/mm10_sacCer3.fa.gz
    dbsnpvcf=/dev/null
    ;;
rn6_sacCer3)
    bwaIndex=/vol/isg/annotation/bwaIndex/rn6_sacCer3/rn6_sacCer3
    ploidy=""
    referencefasta=/vol/isg/annotation/fasta/rn6_sacCer3/rn6_sacCer3.fa.gz
    dbsnpvcf=/dev/null
    ;;
cegsvectors)
    bwaIndex=/vol/isg/annotation/bwaIndex/cegsvectors/cegsvectors
    ploidy="--ploidy 1"
    referencefasta=/vol/mauranolab/cegs/sequences/vectors/vectors.incells.fa
    dbsnpvcf=/dev/null
    ;;
*)
    echo "Don't recognize genome ${mappedgenome}";
    exit 1;;
esac


if [[ "${mappedgenome}" == "cegsvectors" ]]; then
    chromsizes="/vol/mauranolab/cegs/sequences/vectors/vectors.incells.chrom.sizes"
else
    chromsizes="/vol/isg/annotation/fasta/${mappedgenome}/${mappedgenome}.chrom.sizes"
fi


echo "genomeinfo for ${mappedgenome}: bwaIndex=${bwaIndex}, ploidy=${ploidy}, referencefasta=${referencefasta}, dbsnpvcf=${dbsnpvcf}, chromsizes=${chromsizes}"

