#!/bin/bash
set -e -o pipefail

#sample should be in this format cellType-DS.hg38
sample=$1 
echo $sample
sampleOutdir=$sample

mkdir -p ${sampleOutdir}/hotspots2
#Submit hotspots2
~/hotspot2/scripts/hotspot2.sh  -c /home/maagj01/projects/Hotspot1vs2/hg38Sort.bed -C /home/maagj01/projects/Hotspot1vs2/hg38Sort.starch  -F 0.20 -f 0.20 -M /vol/isg/annotation/bed/hg38/mappability/hg38.K36.mappable_only.starch ${sampleOutdir}/${sample}.bam ${sampleOutdir}/hotspots2/

for FDR in {0.05,0.01,0.005,0.001} 
do
       ~/hotspot2/scripts/hsmerge.sh -f ${FDR} ${sampleOutdir}/hotspots2/${sample}.allcalls.starch ${sampleOutdir}/hotspots2/${sample}.hotspots.fdr${FDR}.starch
done


