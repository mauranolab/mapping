#!/bin/bash
set -u

#On analysis sets
#https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
#https://gatkforums.broadinstitute.org/gatk/discussion/7857/reference-genome-components

#- Neither ENCODE hg38_no_alt_analysis_set nor hg38 have PAR masked
#PAR in hg38 from https://www.ncbi.nlm.nih.gov/grc/human:
#chrY	10001	2781479	PAR1
#chrY	56887903	57217415	PAR2
#chrX	10001	2781479	PAR1
#chrX	155701383	156030895	PAR2
#
#mm10 https://www.ncbi.nlm.nih.gov/grc/mouse uncliear if 1-indexed?
#chrY	90745845 	91644698	PAR
#chrX	169969759	170931299	PAR
#
#Rat is unclear which end of chrX it's on
#
#- There are no non-ACGTN in mm10, hg19/38 (from UCSC goldenpath, ENCODE analysis set, or the UCSC-provided NCBI analysis set) or rn6 (UCSC), but NCBI hg38 analysis set has 94
#zcat /vol/isg/annotation/fasta/mm10/mm10.fa.gz | grep -v "^>" | perl -pe 's/[ACGTN\n]+//ig;' | wc -c
#
#UCSC source: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/
#NCBI source (has non-N ambiguous bases and decoys):ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
#ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCA_000001405.27_GRCh38.p12
#
#- Decided not to use decoy hs38d1 -- couldn't find info on what exactly is in here but I think has EBV, HHV, etc. and other unmapping bits empirically found in HuRef and GM12878 by Heng Li and/or Broad

#Notes
#- No non-encode analysis set for mm10 seems available
#- None of the above seem to use latest hg38 (p12 as of 10/15/18), can't find good inventory of changes


mappedgenome=${1}


#Deal with some of the more complex reference index names
#NB this will call hotspots, etc. only on the first mammalian genome for the *_sacCer3 hybrid indices
annotationgenome=`echo ${mappedgenome} | perl -pe 's/_.+$//g;' -e 's/all$//g;'`

chromsizes="/vol/isg/annotation/fasta/${mappedgenome}/${mappedgenome}.chrom.sizes"

#requires delly module be loaded
dellypath=`which delly | xargs dirname | xargs dirname`
#default to no exclusion
dellyexclude=""

case "${mappedgenome}" in
#Reference sequence can be .fa.gz as long as bgzip was used to compress
hg19)
    #NB not the analysis set
    bwaIndex=/vol/isg/annotation/bwaIndex/hg19all/hg19all
    ploidy="--ploidy GRCh37"
    referencefasta=/vol/isg/annotation/fasta/hg19/hg19.fa
    dbsnpvcf=/vol/isg/annotation/bed/hg19/snp151/src/All_20180418.vcf.gz
    dellyexclude="-x ${dellypath}/excludeTemplates/human.hg19.excl.tsv"
    ;;
hg38_noalt|hg38_full|hg38_sacCer3)
    bwaIndex=/vol/isg/annotation/bwaIndex/${mappedgenome}/${mappedgenome}
    ploidy="--ploidy GRCh38"
    referencefasta=/vol/isg/annotation/fasta/${mappedgenome}/${mappedgenome}.fa.gz
    dbsnpvcf=/vol/isg/annotation/bed/hg38/snp151/src/All_20180418.vcf.gz
    dellyexclude="-x ${dellypath}/excludeTemplates/human.hg38.excl.tsv"
    ;;
mm10)
    bwaIndex=/vol/isg/annotation/bwaIndex/mm10_no_alt_analysis_set/mm10_no_alt_analysis_set
    ploidy="--ploidy-file /vol/isg/annotation/fasta/mm10_no_alt_analysis_set/mm10_no_alt_analysis_set.ploidy.txt"
    referencefasta=/vol/isg/annotation/fasta/mm10_no_alt_analysis_set/mm10_no_alt_analysis_set.fa.gz
    #non-human accession repository has switched from NCBI/dbSNP to EVA as of 9/2017 but I can't find a VCF export, so am using the Sanger mouse genome VCF to annotate variation that segregates among their strains rather the frozen dbSNP export (which I think is here: https://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/)
    #the v6 genotypes have changed dramatically so stick to v5 for now
    dbsnpvcf=/vol/mauranolab/mauram01/hybridmice/genotyping/v5/mgp.v5.merged.snps.indels.dbSNP142.vcf.gz
    dellyexclude="-x ${dellypath}/excludeTemplates/mouse.mm10.excl.tsv"
    ;;
mm10_sacCer3)
    bwaIndex=/vol/isg/annotation/bwaIndex/${mappedgenome}/${mappedgenome}
    ploidy="--ploidy-file /vol/isg/annotation/fasta/mm10_no_alt_analysis_set/mm10_no_alt_analysis_set.ploidy.txt"
    referencefasta=/vol/isg/annotation/fasta/${mappedgenome}/${mappedgenome}.fa.gz
    dbsnpvcf=/vol/mauranolab/mauram01/hybridmice/genotyping/v5/mgp.v5.merged.snps.indels.dbSNP142.vcf.gz
    dellyexclude="-x ${dellypath}/excludeTemplates/mouse.mm10.excl.tsv"
    ;;
mm10all_CASTEiJ_female)
    bwaIndex=/vol/isg/annotation/bwaIndex/${mappedgenome}/${mappedgenome}
    ploidy="--ploidy-file /vol/isg/annotation/fasta/mm10_no_alt_analysis_set/mm10_no_alt_analysis_set.ploidy.txt"
    referencefasta=/vol/isg/annotation/bwaIndex/${mappedgenome}/${mappedgenome}.fa.gz
    dbsnpvcf=/dev/null
    #BUGBUG no dellyexclude for strain-specific references
    chromsizes="/vol/isg/annotation/bwaIndex/${mappedgenome}/${mappedgenome}.chrom.sizes"
    ;;
rn6|rn6_sacCer3)
    bwaIndex=/vol/isg/annotation/bwaIndex/${mappedgenome}/${mappedgenome}
    ploidy="--ploidy-file /vol/isg/annotation/fasta/rn6/rn6.ploidy.txt"
    referencefasta=/vol/isg/annotation/fasta/${mappedgenome}/${mappedgenome}.fa.gz
    dbsnpvcf=/dev/null
    ;;
sacCer3)
    bwaIndex=/vol/isg/annotation/bwaIndex/${mappedgenome}/${mappedgenome}
    #Maybe better to leave as diploid?
    ploidy="--ploidy 1"
    referencefasta=/vol/isg/annotation/fasta/${mappedgenome}/${mappedgenome}.fa.gz
    dbsnpvcf=/dev/null
    dellyexclude="-x ${dellypath}/excludeTemplates/yeast.sacCer3.excl.tsv"
    ;;
wuhCor1|hg38_full_wuhCor1)
    bwaIndex=/vol/isg/annotation/bwaIndex/${mappedgenome}/${mappedgenome}
    ploidy="--ploidy 1"
    referencefasta=/vol/isg/annotation/fasta/${mappedgenome}/${mappedgenome}.fa.gz
    dbsnpvcf=/dev/null
    ;;
talOcc4|talOcc4_sacCer3)
    bwaIndex=/vol/isg/annotation/bwaIndex/${mappedgenome}/${mappedgenome}
    ploidy="--ploidy-file /vol/isg/annotation/fasta/talOcc4/talOcc4.ploidy.txt"
    referencefasta=/vol/isg/annotation/fasta/${mappedgenome}/${mappedgenome}.fa.gz
    dbsnpvcf=/dev/null
    ;;
cegsvectors*)
    #If this genome is an symlink, then substitute it with its target for the purposes of looking up bwa index and fasta file
    bwagenome=`readlink -f /vol/cegs/sequences/${mappedgenome} | xargs basename`
    bwaIndex=/vol/cegs/sequences/${bwagenome}/${bwagenome}
    #In case this genome is an symlink,  use its target name to get the right chrom.sizes name
    chromsizes="/vol/cegs/sequences/${bwagenome}/${bwagenome}.chrom.sizes"
    ploidy="--ploidy 1"
    #Tried using the full cegsvectors.fa.gz but picard (esp CollectMultipleMetrics) has trouble with it)
    referencefasta=/vol/cegs/sequences/${bwagenome}/${bwagenome}.fa.gz
    dbsnpvcf=/dev/null
    ;;
t2t)
    bwaIndex=/vol/isg/annotation/bwaIndex/${mappedgenome}/${mappedgenome}
    ploidy="--ploidy GRCh38"
    referencefasta=/vol/isg/annotation/fasta/${mappedgenome}/${mappedgenome}.fa.gz
    dbsnpvcf=/dev/null
    ;;
*)
    echo "ERROR: Don't recognize genome ${mappedgenome}";
    exit 1;;
esac


echo "genomeinfo for ${mappedgenome}: bwaIndex=${bwaIndex}, ploidy=${ploidy}, referencefasta=${referencefasta}, dbsnpvcf=${dbsnpvcf}, annotationgenome=${annotationgenome}, chromsizes=${chromsizes}, dellyexclude=${dellyexclude}"

if [[ ! -d `dirname "${bwaIndex}"` ]]; then
    echo "ERROR could not find required files for ${bwaIndex}!"
    exit 1
fi

if [[ ! -s "${referencefasta}" ]]; then
    echo "ERROR could not find required file ${referencefasta}!"
    exit 2
fi

if [[ "${dbsnpvcf}" != "/dev/null" ]] && [[ ! -s "${dbsnpvcf}" ]]; then
    echo "ERROR could not find required file ${dbsnpvcf}!"
    exit 3
fi

if [[ ! -s "${chromsizes}" ]]; then
    echo "ERROR could not find required file ${chromsizes}!"
    exit 4
fi
