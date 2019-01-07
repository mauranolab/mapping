#!/bin/bash
set -e # -o pipefail

module load trimmomatic/0.36 weblogo/3.5.0 ImageMagick picard/1.140 FastQC/0.11.4 samtools/1.3.1 bowtie2/2.2.9 

# the function "round()" was taken from 
# https://stempell.com/2009/08/rechnen-in-bash/
floor()
{
       echo $(printf %.$2f $(echo "scale=0;$1/1" | bc))
};

src=/home/maagj01/scratch/transposon/src

mkdir -p captureC

sample=$1
REconfig=$2
#R1/2trim clip before anything (i.e. for NNNN used to create diversity)
REfragment=$3
REjobName=$4
organism=$5

OUTDIR=${sample}
mkdir -p $OUTDIR

shift 5
basedir=$@
echo
echo "Looking for files in $basedir"
f1=`find $basedir/ -name "*_R1_*.fastq.gz"`
f2=`find $basedir/ -name "*_R2_*.fastq.gz"`






echo "Weblogo of raw reads"
zcat -f $f1 |   awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2'| shuf -n 1000000| awk '{print ">id-" NR; print}' |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} R1 raw sequence" --stacks-per-line 100 > $TMPDIR/${sample}.R1.raw.eps
convert $TMPDIR/${sample}.R1.raw.eps $OUTDIR/${sample}.R1.raw.png

zcat -f $f2 |   awk -F "\t" 'BEGIN {OFS="\t"} NR % 4 == 2'| shuf -n 1000000| awk '{print ">id-" NR; print}' |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} R2 raw sequence" --stacks-per-line 100 > $TMPDIR/${sample}.R2.raw.eps
convert $TMPDIR/${sample}.R2.raw.eps $OUTDIR/${sample}.R2.raw.png


echo 'Filtering poly-G reads'
/home/maagj01/scratch/transposon/src/filterNextSeqReadsForPolyG.py --inputR1 $f1 --inputR2 $f2 --maxPolyG 75 --outputR1 $TMPDIR/${sample}.R1.fastq.gz --outputR2 $TMPDIR/${sample}.R2.fastq.gz


#Regular illumina dsDNA protocol
illuminaAdapters=/home/maagj01/scratch/transposon/captureC/AdapterLinker/TrueSeq3-PE-2.fa 
#ATAC-seq
#illuminaAdapters=/cm/shared/apps/trimmomatic/0.36/adapters/NexteraPE-PE.fa
seedmis=0
#Pretty much anything below 10 works
PEthresh=0
SEthresh=0
mintrim=0
keepReverseReads=true
trimmomaticBaseOpts="-threads $NSLOTS"
echo "Trimming adapters"
trimmomaticSteps="ILLUMINACLIP:${illuminaAdapters}:0:10:3 MINLEN:20"
#trimmomaticSteps="ILLUMINACLIP:$illuminaAdapters:$seedmis:$PEthresh:$SEthresh:$mintrim:$keepReverseReads"

#MAXINFO:27:0.95 TRAILING:20

R1out=$OUTDIR/${sample}_R1.fastq.gz 
R1outUnpaired=$OUTDIR/${sample}_R1_unpaired.fastq.gz 

R2out=$OUTDIR/${sample}_R2.fastq.gz 
R2outUnpaired=$OUTDIR/${sample}_R2_unpaired.fastq.gz 
java org.usadellab.trimmomatic.TrimmomaticPE $TMPDIR/${sample}.R1.fastq.gz $TMPDIR/${sample}.R2.fastq.gz $R1out $R1outUnpaired $R2out $R2outUnpaired $trimmomaticSteps


mkdir -p captureC/${sample}/${sample}
module unload python/3.5.0 
module load python/2.7.10 
echo
echo 'Splitting files into 5M read chunks'
python ~/src/HiC-Pro/bin/utils/split_reads.py -r captureC/${sample}/${sample} --nreads 5000000 $R1out
python ~/src/HiC-Pro/bin/utils/split_reads.py -r captureC/${sample}/${sample} --nreads 5000000 $R2out

#if [[ `ls captureC/${sample}/${sample}/*.fastq| wc -l` -eq 0 ]]
for i in `ls captureC/${sample}/${sample}`; do  mv captureC/${sample}/${sample}/${i} captureC/${sample}/${sample}/${i}.fastq;done
#;done
#else 
#       rm captureC/${sample}/${sample}/*.fastq
#for i in `ls captureC/${sample}/${sample}`; do  mv captureC/${sample}/${sample}/${i} captureC/${sample}/${sample}/${i}.fastq ;done
#fi

CurDir=`pwd`
#sample=$(echo $sample|awk -F'-' '{print $1}')
echo
echo 'Creating submit config file'

~/src/HiC-Pro/bin/HiC-Pro  -c ${REconfig} -i ${CurDir}/captureC/${sample} -o ${CurDir}/${sample} -p <<-EOF
yes
EOF
module unload python/2.7.10 
module load python/3.5.0 

Step1=$(echo "HiCPro_step1_"${REjobName}".sh")
Step2=$(echo "HiCPro_step2_"${REjobName}".sh")


cd ${sample} 
numjobs=`ls -l  ${CurDir}/captureC/${sample}/${sample}|grep _R1_| wc -l`
echo
echo $numjobs
echo "Submitting bowtie2 alignment in ${numjobs} chunks"
qsub -S /bin/bash -t 1-${numjobs} -terse -j y --qos=full -N HiCstep1.${sample} -b y "module unload python/3.5.0; module load python/2.7.10; bash ${Step1}" 

module unload python/2.7.10 
module load python/3.5.0 


echo
echo "Submitting HiC analysis"
qsub -S /bin/bash  -j y --qos=full -hold_jid HiCstep1.${sample}   -N HiCstep2.${sample} -b y "module unload python/3.5.0; module load python/2.7.10; bash ${Step2}" | awk '{print $NF}' >../sgid.${sample}

cd ..

bamfilesR1=$(find captureC/${sample}/${sample} -name *_R1_*fastq|awk -F'/' '{print $NF}'|sed "s/.fastq/_${organism}all.bwt2merged.bam/g"| sed "s/^/${sample}\/bowtie_results\/bwt2\/${sample}\//g"|tr '\n' ' ')
bamfilesR2=$(find captureC/${sample}/${sample} -name *_R2_*fastq|awk -F'/' '{print $NF}'|sed "s/.fastq/_${organism}all.bwt2merged.bam/g"| sed "s/^/${sample}\/bowtie_results\/bwt2\/${sample}\//g"|tr '\n' ' ')

samfiles=$(find captureC/${sample}/${sample} -name *_R1_*fastq|awk -F'/' '{print $NF}'|sed 's/_R1_/_/g'| sed "s/.fastq/_${organism}all.bwt2pairs_interaction.sam/g"|sed "s/^/${sample}\/hic_results\/data\/${sample}\//g"|tr '\n' ' ')

BamfilesFromSam=$(echo $samfiles|sed 's/.sam/.bam/g')

cat <<EOF | qsub -S /bin/bash -terse -hold_jid `cat sgid.${sample}` -j y --qos=full -N Merging_${sample} -b y | perl -pe 's/[^\d].+$//g;'
#cat <<EOF | qsub -S /bin/bash -terse  -j y --qos=full -N ${sample} -b y | perl -pe 's/[^\d].+$//g;'
set -e -o pipefail
echo "Merging bam files to find location of mapped reads"
echo "Merging bam files"
echo $bamfilesR1


if [[ `echo $bamfilesR1|wc |awk '{print $2}'` -gt 1 ]]
then 
       samtools merge -f -l 9 - $bamfilesR1|samtools sort - $OUTDIR/${sample}_R1.sorted
       samtools merge -f -l 9 - $bamfilesR2|samtools sort - $OUTDIR/${sample}_R2.sorted
       

else 
       samtools sort $bamfilesR1 $OUTDIR/${sample}_R1.sorted
       samtools sort $bamfilesR2 $OUTDIR/${sample}_R2.sorted
fi
samtools index $OUTDIR/${sample}_R1.sorted.bam
samtools index $OUTDIR/${sample}_R2.sorted.bam

#echo 'Sorting interacting samfiles'
#for samFile in `echo $samfiles`; do bamFile=$(echo $samFile|sed 's/.sam//g'); samtools view -b $samFile| samtools sort - $bamFile ;done
#
#echo 'Merging samfiles'
#if [[ `echo $samfiles|wc |awk '{print $2}'` -gt 1 ]]
#then 
#       samtools merge -f -l 9 ${sample}/hic_results/data/${sample}/${sample}.interactions.bam $BamfilesFromSam
#       
#else 
#       cp $BamfilesFromSam ${sample}/hic_results/data/${sample}/${sample}.interactions.bam
#fi
#rm $BamfilesFromSam
EOF

#rm -f sgeid.Step2.${sample}
#rm -f Merging_${sample}
rm -f sgid.${sample}

echo
echo "Submitting HiC analysis"
qsub -S /bin/bash  -j y --qos=full -hold_jid HiCstep2.${sample} -N HiCSTAT.${sample} -b y -V "bash ~/scratch/transposon/captureC/src/Capture_Read1read2_mapping.sh ${sample}"


echo 'Done!!!'

#rm -r ${CurDir}/captureC/${sample}

#module unload python/2.7.10 
#module load python/3.5.0 
#qsub -S /bin/bash -j y --qos=full  -N viewPoint.${sample} -b y "awk '$2==$5' hic_results/data/${sample}/${sample}_allValidPairs  > hic_results/data/${sample}/${sample}_allValidPairs_sameChrom; module unload python/3.5.0; module load python/2.7.10; python ~/src/HiC-Pro/bin/utils/make_viewpoints.py -i hic_results/data/${sample}/${sample}_allValidPairs_sameChrom -f ~/scratch/transposon/captureC/genomeFrag/hg38_ncoI.bed -t ~/scratch/transposon/captureC/capture/NcoI_capture.bed  >ViewCapture.bedgraph"

