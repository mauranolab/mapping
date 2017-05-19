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


genome=$1
celltype=$2
DS=$3
name=${celltype}-${DS}.${genome}

###Files will be outputted to a folder called name (celltype-DS-genome)


#scripts from https://github.com/mauranolab/pipelines/blob/master/dnase/
src=/vol/isg/encode/dnase/src/


sampleOutdir=$name
mkdir -p $sampleOutdir

#mkdir -p tmp

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
qsub -p -450 -S /bin/bash -cwd -V --qos=normal -pe threads 2 -terse -j y -b y -t $firstline-$lastline -o $sampleOutdir/ -N map.$name "$src/map.sh $celltype $DS $genome $src" | perl -pe 's/[^\d].+$//g;' > $sampleOutdir/sgeid

echo -n "Your job "
cat $sampleOutdir/sgeid | perl -pe 's/\n/ /g;'
echo "has been submitted"


echo "$name merge"
qsub -p -400 -S /bin/bash -cwd -V --qos=normal -terse -j y -b y -hold_jid `cat $sampleOutdir/sgeid` -o $sampleOutdir/ -N merge.$name "$src/merge.sh $name $DS $genome" | perl -pe 's/[^\d].+$//g;' > $sampleOutdir/sgeid.merge.$name

echo "$name makeTracks"
qsub -p -400 -S /bin/bash -cwd -V --qos=normal -terse -j y -b y -hold_jid `cat $sampleOutdir/sgeid.merge.$name` -o $sampleOutdir/ -N makeTracks.$name "$src/makeTracks.sh $name $DS $genome $src" | perl -pe 's/[^\d].+$//g;' > $sampleOutdir/sgeid.makeTracks.$name

echo "$name hotspots2"
qsub -p -450 -S /bin/bash -cwd -V --qos=normal -pe threads 1 -terse -j y -b y -hold_jid `cat $sampleOutdir/sgeid.makeTracks.$name` -o $sampleOutdir/ -N hotspot2.$name "$src/hotspot2.sh $name" |perl -pe 's/[^\d].+$//g;'  > $sampleOutdir/sgeid.hot2.$name

rm -f $sampleOutdir/sgeid $sampleOutdir/sgeid.merge.$name $sampleOutdir/sgeid.makeTracks.$name $sampleOutdir/sgeid.hot2.$name

echo
