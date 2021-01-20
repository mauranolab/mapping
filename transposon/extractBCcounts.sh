#!/bin/bash
set -eu -o pipefail


#NSLOTS=1

src=/vol/mauranolab/mapped/src/transposon


###Parse command line args
if [ "$#" -lt 6 ]; then
    echo "Not enough arguments"
    exit 1
fi

sample=$1
BCreadSeq=$2
bclen=$3
chunksize=$4
plasmidSeq=$5


#Just take the rest to simplify passing multiple arguments for extractBarcode.py
shift 5
extractBCargs=$@
#TODO extractBarcode.py has the $TMPDIR/${sample}.plasmid.fastq.gz as input


#Save the name before we append jobid to ${sample} for simplicity
BCreadFile=${sample}.BC.fastq.gz
PlasmidreadFile=${sample}.plasmid.fastq.gz
jobid=${SGE_TASK_ID}
#jobid=1
OUTDIR=${sample}
sample="${sample}.$jobid"


firstline=`echo "$chunksize * ($jobid - 1) + 1" | bc -l -q`
lastline=`echo "$chunksize * $jobid" | bc -l -q`
echo "Running on $HOSTNAME. Output to $OUTDIR. Jobid=$jobid (lines ${firstline}-${lastline})"
date
zcat -f $OUTDIR/${PlasmidreadFile} | awk -v firstline=$firstline -v lastline=$lastline 'NR>=firstline && NR<=lastline' | pigz -p ${NSLOTS} -c -1 - > $TMPDIR/${sample}.plasmid.fastq.gz


echo
echo "Extracting barcodes from ${BCreadFile}"
date
if [[ "${plasmidSeq}" == "None" ]]; then
    echo "No plasmid sequence provided, will extract barcodes using BC read only"
    plasmidcmd=""
else
    plasmidcmd="--plasmidRefSeq ${plasmidSeq} --plasmidRead $TMPDIR/${sample}.plasmid.fastq.gz"
fi

#TODO since extraction is fairly fast, it might be a bit more efficient to do this in the single thread submit script so that adjacency dedup can be run on a constant number of non-empty BCs
zcat -f $OUTDIR/${BCreadFile} | awk -v firstline=$firstline -v lastline=$lastline 'NR>=firstline && NR<=lastline' |
${src}/extractBarcode.py --BCread - --bcRefSeq ${BCreadSeq} --bclen ${bclen} ${plasmidcmd} --minBaseQ 30 ${extractBCargs} | pigz -p ${NSLOTS} -c -1 > $TMPDIR/${sample}.barcodes.raw.txt.gz
date


echo
echo "Correcting barcodes"
date
echo "Deduping BCs"
zcat $TMPDIR/${sample}.barcodes.raw.txt.gz | ${src}/AdjacencyDeDup.py --col 1 -o - -  | 
#Mask UMI and cell BC for reads where BC was ambiguous
awk -F "\t" 'BEGIN {OFS="\t"} $1=="" {$3=""; $4=""} {print}' |
pigz -p ${NSLOTS} -c -1 > $TMPDIR/${sample}.barcodes.deduped.txt.gz
date


echo
#Go through entire file with awk despite only looking at first line so zcat terminates properly
minCellBCLength=`zcat -f $TMPDIR/${sample}.barcodes.deduped.txt.gz | awk -F "\t" 'BEGIN {OFS="\t"} $1!="" {print length($4)}' | uniq | sort -n | awk 'NR==1'`
if [[ "${minCellBCLength}" -gt "4" ]]; then
    echo "Deduping cellBCs"
    date
    zcat -f $TMPDIR/${sample}.barcodes.deduped.txt.gz |
    ${src}/AdjacencyDeDup.py --col 4 -o - - | 
    #Mask BC and UMI for reads where cellBC was ambiguous
    awk -F "\t" 'BEGIN {OFS="\t"} $4=="" {$1=""; $3=""} {print}' |
    pigz -p ${NSLOTS} -c -1 > $TMPDIR/${sample}.barcodes.deduped.cellBCdeduped.txt.gz
    mv $TMPDIR/${sample}.barcodes.deduped.cellBCdeduped.txt.gz $TMPDIR/${sample}.barcodes.deduped.txt.gz
else
    #Skip cellBC deduplication
    echo "cellBC too short (${minCellBCLength}) -- skip cellBC deduplication"
fi


echo
#Only dedup UMIs if the length is over 4
#Go through entire file with awk despite only looking at first line so zcat terminates properly
minUMILength=`zcat -f $TMPDIR/${sample}.barcodes.deduped.txt.gz | awk -F "\t" 'BEGIN {OFS="\t"} $1!="" {print length($3)}' | uniq | sort -n | awk 'NR==1'`
if [[ "${minUMILength}" -gt "4" ]]; then
    echo "Deduping UMIs per barcode"
    #BUGBUG seems to take extremely long time in some cases -- e.g. BS01470A
    date
    zcat -f $TMPDIR/${sample}.barcodes.deduped.txt.gz |
    ##Here we stop being in order of original fastq
    sort -k1,1 |
    ${src}/AdjacencyDeDup.py --col 3 --groupcols 1 -o - - | 
    #Mask BC and cellBC for reads where UMI was ambiguous -- I don't think this affects counts since both alternate UMIs must be present and this would just be a duplicate, but I haven't checked. Spot checking shows the failed UMIs have lots of G
    awk -F "\t" 'BEGIN {OFS="\t"} $3=="" {$1=""; $4=""} {print}' |
    pigz -p ${NSLOTS} -c -9 > $OUTDIR/${sample}.barcodes.txt.gz
else
    #Skip UMI deduplication
    echo "UMI too short (${minUMILength}) -- skip UMI deduplication"
    zcat -f $TMPDIR/${sample}.barcodes.deduped.txt.gz | pigz -p ${NSLOTS} -c -9 > $OUTDIR/${sample}.barcodes.txt.gz
fi


echo
echo "Done!!!"
date
