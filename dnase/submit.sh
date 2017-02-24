#!/bin/bash
set -e

#Load modules
module load picard/1.140
module load FastQC/0.11.4
module load bedops/2.4.19
module load bwa/0.7.15
module load htslib/1.2.1
module load samtools/1.2
module load trimmomatic/0.33
module load python/3.5.0
module load hotspot/4.1
module load samblaster/0.1.22


genome=$1
celltype=$2
DS=$3
name=${celltype}-${DS}

base=/vol/isg/encode/chipseq

mkdir -p tmp

if ! grep -q $DS inputs.txt; then
       echo "Can't find $DS"
       exit 1
fi

numlines=`grep $DS inputs.txt | wc -l`

firstline=`grep -n $DS inputs.txt | head -1 | awk -F ":" '{print $1}'`
lastline=`grep -n $DS inputs.txt | tail -1 | awk -F ":" '{print $1}'`
expectednumlines=`echo "$lastline-$firstline+1" | bc -l -q`
if [ "$expectednumlines" -ne "$numlines" ]; then
       echo "Wrong number of lines (expected $expectednumlines, found $numlines), are all fastq files contiguous in inputs.txt?"
       exit 2
fi


#SGE doesn't accept a complicated -t array, so we'll start R2 jobs that will die instantly rather than prune here
echo "Processing $name (input.txt lines $firstline-$lastline) for genome $genome"
qsub -p -450 -S /bin/bash -cwd -V  -pe threads 4 -terse -j y -b y -t $firstline-$lastline -N map.$name "$base/src/map.sh $celltype $DS $genome" | perl -pe 's/[^\d].+$//g;' > sgeid

echo -n "Your job "
cat sgeid | perl -pe 's/\n/ /g;'
echo "has been submitted"


name=$celltype-${DS}.${genome}

echo "$genome merge"
qsub -p -400 -S /bin/bash -cwd -V  -terse -j y -b y -hold_jid `cat sgeid` -N merge.$name "$base/src/merge.sh $name $DS $genome" | perl -pe 's/[^\d].+$//g;' > sgeid.merge.$genome

echo "$genome makeTracks"
qsub -p -400 -S /bin/bash -cwd -V  -terse -j y -b y -hold_jid `cat sgeid.merge.$genome` -N makeTracks.$name "$base/src/makeTracks.sh $name $DS $genome" | perl -pe 's/[^\d].+$//g;' > sgeid.makeTracks.$genome


rm -f sgeid sgeid.merge.$genome sgeid.makeTracks.$genome

echo
