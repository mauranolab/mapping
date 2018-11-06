#!/bin/bash
set -eu -o pipefail


#NSLOTS=1

src=/vol/mauranolab/transposon/src


###Parse command line args
if [ "$#" -ne 6 ]; then
    echo "Wrong number of arguments"
    exit 1
fi

sample=$1
BCreadSeq=$2
bclen=$3
chunksize=$4
plasmidSeq=$5
extractBCargs=$6
#TODO extractBarcode.py has the $TMPDIR/${sample}.trimmed.plasmid.fastq.gz as input


#Save the name before we append jobid to $sample for simplicity
BCreadFile=$sample.trimmed.BC.fastq.gz
PlasmidreadFile=$sample.trimmed.plasmid.fastq.gz
jobid=${SGE_TASK_ID}
#jobid=1
OUTDIR=$sample
sample="${sample}.$jobid"


firstline=`echo "$chunksize * ($jobid - 1) + 1" | bc -l -q`
lastline=`echo "$chunksize * $jobid" | bc -l -q`
echo "Running on $HOSTNAME. Output to $OUTDIR. Jobid=$jobid (lines ${firstline}-${lastline})"
date


echo
echo "Extracting barcodes from $BCreadFile"
date
zcat -f $OUTDIR/$PlasmidreadFile | awk -v firstline=$firstline -v lastline=$lastline 'NR>=firstline && NR<=lastline' | gzip -1 -c - > $TMPDIR/$sample.plasmid.fastq.gz


if [[ "$plasmidSeq" == "None" ]]; then
    echo "No plasmid sequence provided, will extract barcodes from BC read only"
    plasmidcmd=""
else
    plasmidcmd="--plasmidSeq $plasmidSeq --plasmidRead $TMPDIR/$sample.plasmid.fastq.gz"
fi

zcat -f $OUTDIR/$BCreadFile | awk -v firstline=$firstline -v lastline=$lastline 'NR>=firstline && NR<=lastline' |
$src/extractBarcode.py --BCread - --referenceSeq $BCreadSeq --bclen $bclen $plasmidcmd --minBaseQ 30 $extractBCargs > $OUTDIR/$sample.barcodes.raw.txt
date


echo
echo "Merging similar barcodes"
date
cat $OUTDIR/$sample.barcodes.raw.txt | python $src/AdjacencyDeDup.py --col 1 -o $OUTDIR/$sample.barcodes.deduped.txt -  
date

#Only merge UMIs if the length is over 4
UMIlength=$(head -1 $OUTDIR/$sample.barcodes.deduped.txt | cut -f3)
if [[ ${#UMIlength} > "4" ]]; then
    echo "Deduping UMIs"
    echo "Merging similar UMIs per barcode"
    date
    cat $OUTDIR/$sample.barcodes.deduped.txt | 
    ##Here we stop being in order of original fastq
    sort -k1,1 |
    $src/AdjacencyDeDup.py --col 3 --groupcol 1 -o $OUTDIR/$sample.barcodes.txt - 
else
    #Skip UMI deduplication
    echo 'UMI too short'
    echo 'Skip UMI dedup'
    cp $OUTDIR/$sample.barcodes.deduped.txt $OUTDIR/$sample.barcodes.txt
fi


#Skip both dedups
#cp $OUTDIR/$sample.barcodes.raw.txt $OUTDIR/$sample.barcodes.txt


echo "Done!!!"
date
