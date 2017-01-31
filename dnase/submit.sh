#!/bin/bash
set -e

strain=$1
celltype=$2
DS=$3
name=${celltype}-${DS}

base=/vol/isg/encode/dnase

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
echo "Processing $name (input.txt lines $firstline-$lastline) for strain $strain"
qsub -p -450 -S /bin/bash -cwd -V -pe threads 4 -terse -j y -b y -t $firstline-$lastline -N map.$name "bash $base/src/map.sh $celltype $DS $strain" | perl -pe 's/[^\d].+$//g;' > sgeid

echo -n "Your job "
cat sgeid | perl -pe 's/\n/ /g;'
echo "has been submitted"


echo "$strain merge"
qsub -p -400 -S /bin/bash -cwd -V -terse -j y -b y -hold_jid `cat sgeid` -N merge.$name.$strain "bash $base/src/merge.sh $strain $celltype $DS $strain" | perl -pe 's/[^\d].+$//g;' > sgeid.$strain

rm -f sgeid sgeid.$strain

echo
