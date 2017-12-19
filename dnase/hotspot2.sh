#!/bin/bash
set -e -o pipefail
#Load modules
module load picard/1.140
module load FastQC/0.11.4
module load bedops/2.4.26
module load bwa/0.7.15
module load htslib/1.2.1
module load samtools/1.3.1
module load trimmomatic/0.36
module load python/3.5.0
module load hotspot/4.1
module load samblaster/0.1.22
module load ucsckentutils/10132015
module load hotspot2/2.0 
#sample should be in this format cellType-DS.hg38
sample=$1 
echo $sample
sampleOutdir=$sample

mkdir -p ${sampleOutdir}/hotspot2
#mkdir -p ${sampleOutdir}/hotspots2_NEW_noQual
#Submit hotspots2

hotspotBAM=$TMPDIR/${sample}.bam

readsLength20=$(samtools view ${sampleOutdir}/${sample}.bam|awk 'length($10) ==20 {print $10}'| wc -l)
readsTotal=$(samtools view ${sampleOutdir}/${sample}.bam| wc -l)

propReads=$(echo $readsLength20/$readsTotal|bc -l)
if [ `echo "$propReads"|awk '{if ($1>0.25) print 1; else print 0}'` -ge 1 ];
then
       echo "More than 25% of reads are 20bp - q10"
       samflags="-q 10 -F 524"
else
       echo "Less than 25% of reads are 20bp - q20"
       samflags="-q 20 -F 524"
fi;



samtools view $samflags -b -1   ${sampleOutdir}/${sample}.bam > $hotspotBAM

hotspot2.sh  -c /vol/isg/annotation/bed/hg38/hotspots2/hg38.chrom.sizes -C /vol/isg/annotation/bed/hg38/hotspots2/hg38.CenterSites.starch  -F 0.20 -f 0.20 -M /vol/isg/annotation/bed/hg38/mappability/hg38.K20.mappable_only.starch  $hotspotBAM ${sampleOutdir}/hotspot2
#~/apps/hotspot2/scripts/hotspot2.sh  -c /vol/isg/annotation/bed/hg38/hotspots2/hg38.chrom.sizes -C /vol/isg/annotation/bed/hg38/hotspots2/hg38.CenterSites.starch  -F 0.20 -f 0.20 -M /vol/isg/annotation/bed/hg38/mappability/hg38.K36.mappable_only.starch  ${sampleOutdir}/${sample}.bam ${sampleOutdir}/hotspots2_NEW_noQual
#hotspot2.sh  -c /vol/isg/annotation/bed/hg38/hotspots2/hg38.chrom.sizes -C /vol/isg/annotation/bed/hg38/hotspots2/hg38.CenterSites.starch  -F 0.20 -f 0.20 -M /vol/isg/annotation/bed/hg38/mappability/hg38.K36.mappable_only.starch  ${sampleOutdir}/${sample}.bam ~/scratch/transposon/Analysis/hotspots2_NEW_noQual


if [ `echo "$propReads"|awk '{if ($1>0.25) print 1; else print 0}'` -ge 1 ];then 
    hotspot2.sh  -c /vol/isg/annotation/bed/hg38/hotspots2/hg38.chrom.sizes -C /vol/isg/annotation/bed/hg38/hotspots2/hg38.CenterSites.starch  -F 0.20 -f 0.20 -M /vol/isg/annotation/bed/hg38/mappability/hg38.K20.mappable_only.starch  $hotspotBAM ${sampleOutdir}/hotspot2

else
    hotspot2.sh  -c /vol/isg/annotation/bed/hg38/hotspots2/hg38.chrom.sizes -C /vol/isg/annotation/bed/hg38/hotspots2/hg38.CenterSites.starch  -F 0.20 -f 0.20 -M /vol/isg/annotation/bed/hg38/mappability/hg38.K36.mappable_only.starch  $hotspotBAM ${sampleOutdir}/hotspot2

fi;

for FDR in {0.05,0.01,0.005,0.001} 
do
       hsmerge.sh -f ${FDR} ${sampleOutdir}/hotspot2/${sample}.allcalls.starch ${sampleOutdir}/hotspot2/${sample}.hotspots.fdr${FDR}.starch
done


