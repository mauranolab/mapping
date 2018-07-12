#!/bin/bash
set -e -o pipefail

#Load modules
module load picard/1.140
module load FastQC/0.11.4
module load bedops/2.4.30
module load bwa/0.7.15
module load htslib/1.2.1
module load samtools/1.3.1
module load trimmomatic/0.36
module load python/3.5.0
module load samblaster/0.1.22
module load ucsckentutils/10132015
module load hotspot/4.1
module load hotspot2/2.1.1


runMap=1
qsubargs="-p 400 --qos=normal"


genome=$1
celltype=$2
DS=$3
name=${celltype}-${DS}.${genome}
#Files for a given sample will be output to a folder named "$celltype-$DS-$genome"

src=/vol/isg/encode/dnase/src/

if ! grep -q $DS inputs.txt; then
       echo "Can't find $DS"
       exit 1
fi


sampleOutdir=$name
mkdir -p $sampleOutdir


numlines=`grep $DS inputs.txt | wc -l`

firstline=`grep -n $DS inputs.txt | head -1 | awk -F ":" '{print $1}'`
lastline=`grep -n $DS inputs.txt | tail -1 | awk -F ":" '{print $1}'`
expectednumlines=`echo "$lastline-$firstline+1" | bc -l -q`
if [ "$expectednumlines" -ne "$numlines" ]; then
       echo "Wrong number of lines (expected $expectednumlines, found $numlines), are all fastq files contiguous in inputs.txt?"
       exit 2
fi


if [ "$runMap" -eq 1 ]; then
       #SGE doesn't accept a complicated -t array, so we'll start R2 jobs that will die instantly rather than prune here
       echo "Processing $name (input.txt lines $firstline-$lastline) for genome $genome"
       qsub -S /bin/bash -cwd -V $qsubargs -pe threads 2 -terse -j y -b y -t $firstline-$lastline -o $sampleOutdir/ -N map.$name "$src/map.sh $celltype $DS $genome $src" | perl -pe 's/[^\d].+$//g;' > $sampleOutdir/sgeid.map
       
       echo "$name merge"
       qsub -S /bin/bash -cwd -V $qsubargs -terse -j y -b y -hold_jid `cat $sampleOutdir/sgeid.map` -o $sampleOutdir -N merge.$name "$src/merge.sh $name $DS $genome" | perl -pe 's/[^\d].+$//g;' > $sampleOutdir/sgeid.merge.$name
       
       makeTracksHold="-hold_jid `cat $sampleOutdir/sgeid.merge.$name`"
else
       makeTracksHold=""
fi


echo "$name makeTracks"
qsub -S /bin/bash -cwd -V $qsubargs -terse -j y -b y $makeTracksHold -o $sampleOutdir -N makeTracks.$name "$src/makeTracks.sh $name $DS $genome $src" | perl -pe 's/[^\d].+$//g;' > $sampleOutdir/sgeid.makeTracks.$name


rm -f $sampleOutdir/sgeid.map $sampleOutdir/sgeid.merge.$name $sampleOutdir/sgeid.makeTracks.$name

echo
