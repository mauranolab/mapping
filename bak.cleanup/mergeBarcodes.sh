#qsub -S /bin/bash -hold_jid `cat sgeid.${sample}` -j y -N ${sample} -b y bash ~/scratch/transposon/src/mergeBarcodes.sh ${sample}"


#!/bin/bash
set -e -o pipefail
#Sort and merge barcodes before and after duplication

NSLOTS=1

sample=$1

#Create outdir
OUTDIR=$sample/barcodes
#Check how many chunks exists
NumberChunks=`ls $OUTDIR/$sample*barcodes.raw.txt| wc -l`
echo $NumberChunks

#Merge ordered chunks of raw barcodes
for i in `seq 1 $NumberChunks`
do
	echo $sample.$i.barcodes.raw.txt
    cat $OUTDIR/$sample.$i.barcodes.raw.txt >> $OUTDIR/$sample.MERGED.barcodes.raw.txt
done

#Merge ordered DEDUPED chunks of barcods
for i in `seq 1 $NumberChunks`
do
	echo $sample.$i.barcodes.txt
    cat $OUTDIR/$sample.$i.barcodes.txt >> $OUTDIR/$sample.MERGED.DEDUPED.barcodes.txt
done

