#!/bin/bash
set -eu -o pipefail

#TODO
#better command line arg processing


#Load modules
module load picard/1.140
module load FastQC/0.11.4
module load bedops/2.4.35
module load bwa/0.7.15
module load htslib/1.9
module load samtools/1.9
module load bcftools/1.9
module load trimmomatic/0.38
module load python/3.5.0
module load samblaster/0.1.22
module load ucsckentutils/10132015
module load hotspot/4.1
module load hotspot2/2.1.1


###Hardcoded configuration options
#Common
qsubargs="-p -400 --qos=normal"

#Mapping pipelines
runMap=1

#Aggregation pipelines
runMerge=1
markdups="markdups=TRUE"


genome=$1
analysisType=$2
celltype=$3
DS=$4
name=${celltype}-${DS}.${genome}
#Files for a given sample will be output to a folder named "$celltype-$DS-$genome"

src=/vol/isg/encode/dnase/src/

if [[ "$analysisType" != "aggregate_dnase" ]] && [[ "$analysisType" != "map_dnase" ]] && [[ "$analysisType" != "map_atac" ]] && [[ "$analysisType" != "map_callsnps" ]]; then 
    echo "ERROR submit: unknown analysisType $analysisType"
    exit 1
fi

if ! grep -q $DS inputs.txt; then
    echo "Can't find $DS"
    exit 2
fi


echo "Will process $name ($DS) using $analysisType pipeline for genome $genome"


sampleOutdir=$name
mkdir -p $sampleOutdir


if [[ "$analysisType" =~ "map" ]]; then
    if [ "$runMap" -eq 1 ]; then
        numlines=`grep $DS inputs.txt | wc -l`

        firstline=`grep -n $DS inputs.txt | head -1 | awk -F ":" '{print $1}'`
        lastline=`grep -n $DS inputs.txt | tail -1 | awk -F ":" '{print $1}'`
        expectednumlines=`echo "$lastline-$firstline+1" | bc -l -q`
        if [ "$expectednumlines" -ne "$numlines" ]; then
            echo "Wrong number of lines (expected $expectednumlines, found $numlines), are all fastq files contiguous in inputs.txt?"
            exit 3
        fi
    
        #SGE doesn't accept a -t specification with gaps, so we'll start R2 jobs that will die instantly rather than prune here
        echo "Will map input.txt lines $firstline-$lastline"
        qsub -S /bin/bash -cwd -V $qsubargs -pe threads 2 -terse -j y -b y -t $firstline-$lastline -o $sampleOutdir -N map.$name "$src/map.sh $genome $analysisType $celltype $DS $src" | perl -pe 's/[^\d].+$//g;' > $sampleOutdir/sgeid.map
    
        echo "$name merge"
        qsub -S /bin/bash -cwd -V $qsubargs -terse -j y -b y -hold_jid `cat $sampleOutdir/sgeid.map` -o $sampleOutdir -N merge.$name "$src/merge.sh $name $DS $genome markdups=TRUE" | perl -pe 's/[^\d].+$//g;' > $sampleOutdir/sgeid.merge.$name
    
        makeTracksHold="-hold_jid `cat $sampleOutdir/sgeid.merge.$name`"
    else
        makeTracksHold=""
    fi
elif [[ "$analysisType" =~ "aggregate" ]]; then
    if [ "$runMerge" -eq 1 ]; then
        echo "$name aggregate"
        qsub -S /bin/bash -cwd -V $qsubargs -terse -j y -b y -o $sampleOutdir -N agg.$name "$src/aggregate.sh $name $DS $genome $markdups" | perl -pe 's/[^\d].+$//g;' > $sampleOutdir/sgeid.agg.$name
    
        makeTracksHold="-hold_jid `cat $sampleOutdir/sgeid.agg.$name`"
    else
        makeTracksHold=""
    fi
else
    echo "ERROR: impossible!"
    exit 4
fi

echo "$name makeTracks"
qsub -S /bin/bash -cwd -V $qsubargs -terse -j y -b y $makeTracksHold -o $sampleOutdir -N makeTracks.$name "$src/makeTracks.sh $genome $analysisType $name $DS $src" | perl -pe 's/[^\d].+$//g;' > $sampleOutdir/sgeid.makeTracks.$name


rm -f $sampleOutdir/sgeid.map $sampleOutdir/sgeid.merge.$name $sampleOutdir/sgeid.makeTracks.$name

echo
