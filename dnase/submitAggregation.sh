#!/bin/bash
set -e -o pipefail

#Load modules
module load picard/1.140
module load bedops/2.4.35
module load htslib/1.9
module load samtools/1.9
module load python/3.5.0
module load samblaster/0.1.22
module load ucsckentutils/10132015
module load hotspot/4.1
module load hotspot2/2.1.1


runMerge=1
markdups="markdups=TRUE"
qsubargs=""


genome=$1
celltype=$2
DS=$3
name=${celltype}-${DS}.${genome}
#Files for a given sample will be output to a folder named "$celltype-$DS-$genome"

src=/vol/mauranolab/mapped/src/

if ! grep -q $DS inputs.txt; then
    echo "Can't find $DS"
    exit 1
fi


sampleOutdir=$name
mkdir -p $sampleOutdir


if [ "$runMerge" -eq 1 ]; then
    echo "$name aggregate"
    qsub -S /bin/bash -cwd -V $qsubargs -terse -j y -b y -o $sampleOutdir -N agg.$name "$src/aggregate.sh $name $DS $genome $markdups" | perl -pe 's/[^\d].+$//g;' > $sampleOutdir/sgeid.agg.$name
    
    makeTracksHold="-hold_jid `cat $sampleOutdir/sgeid.agg.$name`"
else
    makeTracksHold=""
fi


echo "$name makeTracks"
qsub -S /bin/bash -cwd -V $qsubargs -terse -j y -b y $makeTracksHold -o $sampleOutdir/ -N makeTracks.$name "$src/makeTracks.sh $name $DS $genome $src" | perl -pe 's/[^\d].+$//g;' > $sampleOutdir/sgeid.makeTracks.$name


rm -f $sampleOutdir/sgeid.agg.$name $sampleOutdir/sgeid.makeTracks.$name

echo
