#!/bin/bash
set -eu -o pipefail

src=/vol/mauranolab/mapped/src/transposon

#TODO inconsistently applied
minReadCutoff=$1
sample=$2


OUTDIR=${sample}

#Format: BC, readname, UMI, cellBC
if [ ! -s "${OUTDIR}/${sample}.barcodes.preFilter.txt.gz" ]; then
    echo "analyzeBCcounts.sh ERROR: barcode input file ${OUTDIR}/${sample}.barcodes.preFilter.txt.gz does not exist!"
    exit 1
fi


echo "Analyzing data for ${sample} (minReadCutoff=${minReadCutoff})"
date

echo
zcat ${OUTDIR}/${sample}.barcodes.preFilter.txt.gz | ${src}/removeOverrepresentedBCs.py --col 1 --maskcols 3,4 --freq 0.01 -o ${TMPDIR}/${sample}.barcodes.txt -


#Doesn't count length if the BC has been masked to ""
minUMILength=`zcat -f ${TMPDIR}/${sample}.barcodes.txt | awk -F "\t" 'BEGIN {OFS="\t"} $1!="" {print length($3)}' | uniq | sort -n | awk 'BEGIN {ret=0} NR==1 {ret=$0} END {print ret}'`
minCellBCLength=`zcat -f ${TMPDIR}/${sample}.barcodes.txt | awk -F "\t" 'BEGIN {OFS="\t"} $1!="" {print length($4)}' | uniq | sort -n | awk 'BEGIN {ret=0} NR==1 {ret=$0} END {print ret}'`


echo
echo -e "${sample}\tMinimum UMI length\t${minUMILength}"
echo -e "${sample}\tMinimum cellBC length\t${minCellBCLength}"

if [ "${minCellBCLength}" -gt 0 ]; then
    #The BC correction doesn't seem to do much, but the UMI corrections can have an effect
    #TODO apply whitelist ahead of time or correct cellBC/UMI to cellranger data here?
    echo
    echo "Secondary deduping BC by cellBC and UMI by BC+cellBC"
    date
    cat ${TMPDIR}/${sample}.barcodes.txt | sort -k4,4 -k1,1 |
    #Dedup BC by cellBC
    ${src}/AdjacencyDeDup.py --col 1 --groupcols 4 -o - - |
    awk -F "\t" 'BEGIN {OFS="\t"} $1=="" {$3=""; $4=""} {print}' |
    #Need to re-sort since some BCs have changed
    sort -k4,4 -k1,1 |
    #Dedup UMI by BC+cellBC
    ${src}/AdjacencyDeDup.py --col 3 --groupcols 4,1 -o - - |
    awk -F "\t" 'BEGIN {OFS="\t"} $3=="" {$1=""; $4=""} {print}' > ${TMPDIR}/${sample}.barcodes.txt.new
    mv ${TMPDIR}/${sample}.barcodes.txt.new ${TMPDIR}/${sample}.barcodes.txt
fi


echo
echo "Sorting and compressing"
date
sort -k1,1 -k3,3 -k4,4 ${TMPDIR}/${sample}.barcodes.txt | pigz -p ${NSLOTS} -c -9 - > ${OUTDIR}/${sample}.barcodes.txt.gz


if [ `cat ${TMPDIR}/${sample}.barcodes.txt | cut -f1 | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | wc -l` -eq 0 ]; then
    echo "analyzeBCcounts.sh WARNING: no barcodes left, so exiting!"
    echo "Done!!!"
    exit 0
fi


#Empty BCs haven't been filtered yet for non-UMI data
#Already sorted by BC, UMI, CellBC
zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | awk -F "\t" '$1!=""' | cut -f1 | uniq -c | awk 'BEGIN {OFS="\t"} {print $2, $1}' > ${OUTDIR}/${sample}.barcode.counts.txt


###Overall summary statistics
echo
numTotalReads=`zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | wc -l`
echo -e "${sample}\tNumber of total reads\t${numTotalReads}"
echo -n -e "${sample}\tNumber of total read barcodes\t"
zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | awk -F "\t" '$1!="" {count+=1} END {print count}'

#BC weblogo
zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | cut -f1 | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | shuf -n 1000000 | awk -F "\t" 'BEGIN {OFS="\t"} {print ">id-" NR; print}' |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} barcodes" --stacks-per-line 100 > $TMPDIR/${sample}.barcodes.eps
convert $TMPDIR/${sample}.barcodes.eps ${OUTDIR}/${sample}.barcodes.png

##UMIs weblogo
if [ "${minUMILength}" -gt 0 ]; then
    echo
    echo "Generating UMI weblogo"
    #Need to trim in case not all the same length
    zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | cut -f3 | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | shuf -n 1000000 | awk -v trim=${minUMILength} -F "\t" 'BEGIN {OFS="\t"} {print substr($0, 0, trim)}' | awk -F "\t" 'BEGIN {OFS="\t"} {print ">id-" NR; print}' |
    weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} UMI" --stacks-per-line 100 > $TMPDIR/${sample}.UMIs.eps
    convert $TMPDIR/${sample}.UMIs.eps ${OUTDIR}/${sample}.UMIs.png
fi


if [ "${minUMILength}" -ge 5 ]; then
    echo
    echo "Analyzing UMIs"
    date
    
    echo
    echo "Barcode counts (including UMI)"
    zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz |
    awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | cut -f1,3 |
    #Sort on BC, UMI and count reads per BC+UMI
    sort -k1,1 -k2,2 | uniq -c |
    awk 'BEGIN {OFS="\t"} {print $2, $3, $1}' |
    #Format: BC, UMI, n
    sort -k1,1 |
    pigz -p ${NSLOTS} -c -9 > ${OUTDIR}/${sample}.barcode.counts.withUMI.txt.gz
    
    echo "Barcode counts UMI-corrected"
    zcat -f ${OUTDIR}/${sample}.barcode.counts.withUMI.txt.gz | cut -f1 | sort | uniq -c | awk 'BEGIN {OFS="\t"} {print $2, $1}' > ${OUTDIR}/${sample}.barcode.counts.UMI_corrected.txt
    countfile="${OUTDIR}/${sample}.barcode.counts.UMI_corrected.txt"
    
    echo -n -e "${sample}\tNumber of unique barcodes+UMI\t"
    #TODO move to extract.py
    zcat -f ${OUTDIR}/${sample}.barcode.counts.withUMI.txt.gz | wc -l
    
    echo
    echo "Histogram of number of reads per barcode+UMI"
    zcat -f ${OUTDIR}/${sample}.barcode.counts.withUMI.txt.gz | cut -f3 | awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g
    
    echo
    echo "UMI lengths"
    #No more empty BCs are present so no need to check
    zcat -f ${OUTDIR}/${sample}.barcode.counts.withUMI.txt.gz | cut -f2 |
    awk '{print length($0)}' | sort -g | uniq -c | sort -k2,2 | awk '$2!=0'
else
    echo
    echo "No UMIs found"
    
    countfile="${OUTDIR}/${sample}.barcode.counts.txt"
fi


echo
echo "Barcode counts"

echo -n -e "${sample}\tNumber of unique barcodes\t"
#TODO move to extract.py
cat ${countfile} | wc -l

echo -n -e "${sample}\tNumber of unique barcodes passing minimum read cutoff\t"
#TODO move to extract.py
cat ${countfile} | awk -v minReadCutoff=${minReadCutoff} -F "\t" 'BEGIN {OFS="\t"} $2>=minReadCutoff' | wc -l

echo -n -e "${sample}\tNumber of analyzed reads\t"
cat ${countfile} | awk -v minReadCutoff=${minReadCutoff} -F "\t" 'BEGIN {OFS="\t"; sum=0; total=0} $2>=minReadCutoff {sum+=$2} {total+=$2} END {print sum}'

echo
echo "Histogram of number of reads per barcode"
cat ${countfile} | cut -f2 | awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g

echo
echo "Barcode lengths"
cut -f1 ${countfile} |
awk '{print length($0)}' | sort -g | uniq -c | sort -k2,2 | awk '$2!=0'


if [ "${minCellBCLength}" -gt 0 ]; then
    #Hardcoded for now
    #Aug 2019 Pilot
    #scRNAseqbase="/vol/mauranolab/transposon/scrnaseq/merged"
    #Dec 2019 scale
    scRNAseqbase="/vol/mauranolab/transposon/scrnaseq/FCH7Y3NBGXC"
    
    echo
    echo "Generating cellBC weblogo"
    #Need to trim in case not all the same length
    zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | cut -f4 | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | shuf -n 1000000 | awk -v trim=${minCellBCLength} -F "\t" 'BEGIN {OFS="\t"} {print substr($0, 0, trim)}' | awk -F "\t" 'BEGIN {OFS="\t"} {print ">id-" NR; print}' |
    weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} cellBC" --stacks-per-line 100 > $TMPDIR/${sample}.cellBCs.eps
    convert $TMPDIR/${sample}.cellBCs.eps ${OUTDIR}/${sample}.cellBCs.png
    
    echo "Barcode counts including cellBC"
    zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz |
    awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | cut -f1,3,4 |
    #Already sorted by BC, UMI, CellBC
    sort -k1,1 -k2,2 -k4,4 | 
    #Count reads per BC+UMI+CellBC
    uniq -c |
    awk 'BEGIN {OFS="\t"} {print $2, $3, $4, $1}' |
    #Format: BC, UMI, cellBC, n
    
    #Filter on minPropReadsForBC
    #Join in pre-UMI counts in column 5
    join -j 1 -a 1 - ${OUTDIR}/${sample}.barcode.counts.txt |
    #join -t "\t" doesn't work so get back to tabs
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' > $TMPDIR/${sample}.barcodes.minPropReadsForBC.txt
    #Format: BC, UMI, cellBC, n, total counts for this BC
    
    echo -n -e "${sample}\tNumber of reads passing minPropReadsForBC filter\t"
    cat $TMPDIR/${sample}.barcodes.minPropReadsForBC.txt | awk -v minPropReadsForBC=0.05 -F "\t" 'BEGIN {OFS="\t"} $4/$5>minPropReadsForBC {filtered+=$4} {total+=$4} END {print filtered, total}'
    
    cat $TMPDIR/${sample}.barcodes.minPropReadsForBC.txt |
    #Do the actual filter
    awk -v minPropReadsForBC=0 -F "\t" 'BEGIN {OFS="\t"} $4/$5>minPropReadsForBC' | cut -f1-4 |
    #Sort on BC and count reads perBC+CellBC
    sort -k1,1 -k3,3 |
    cut -f1,3 | sort | uniq -c | awk 'BEGIN {OFS="\t"} {print $2, $3, $1}' |
    #Format: BC, cellBC, n
    #TODO parameterize and do more precise filtering of 10x whitelist BCs. we do miss some by not allowing mismatches but it's pretty few reads
    #TODO report high-count cellBCs not on whitelist, mainly in case umi_tools dedup is not working well
    fgrep -w -f ${scRNAseqbase}/cells.whitelist.txt |
    fgrep -v -w -f ${scRNAseqbase}/cells.readcounts.excluded.txt |
    fgrep -v -w -f ${scRNAseqbase}/cells.pSB.excluded.txt |
    awk 'BEGIN {OFS="\t"; print "BC", "cellBC", "count"} {print}' > ${OUTDIR}/${sample}.barcode.counts.byCell.txt
    
    
    echo -n -e "${sample}\tNumber of cells\t"
    cat ${OUTDIR}/${sample}.barcode.counts.byCell.txt | mlr --headerless-csv-output --tsv cut -f cellBC | sort | uniq | wc -l
    
    echo -e "${sample}\tHistogram of number of cells per BC"
    cat ${OUTDIR}/${sample}.barcode.counts.byCell.txt | mlr --headerless-csv-output --tsv cut -f BC | sort | uniq -c | mlr -p --ofs "\t" --ofmt %.1lf histogram -f 1 --nbins 4 --auto -o nCells
    echo
    
    echo -e "${sample}\tHistogram of number of BCs per cell"
    cat ${OUTDIR}/${sample}.barcode.counts.byCell.txt | mlr --headerless-csv-output --tsv cut -f cellBC | sort | uniq -c | mlr -p --ofs "\t" --ofmt %.1lf histogram -f 1 --nbins 4 --auto -o nBCs
    echo
    
    ${src}/genotypeClones.py --inputfilename ${OUTDIR}/${sample}.barcode.counts.byCell.txt --outputwide ${OUTDIR}/${sample}.clones.txt --outputlong ${OUTDIR}/${sample}.clones.counts.filtered.txt --output - --printGraph ${OUTDIR}/clones | mlr --tsv sort -f clone,BC -nr count > ${OUTDIR}/${sample}.barcode.counts.byCell.filtered.txt
    echo
    
    echo -e "${sample}\tHistogram of number of cells per clone"
    cat ${OUTDIR}/${sample}.clones.txt | mlr --tsv --ofmt %.1lf histogram -f ncells --nbins 10 --auto
    echo
    
    echo -e "${sample}\tHistogram of number of BCs per clone"
    cat ${OUTDIR}/${sample}.clones.txt | mlr --tsv --ofmt %.1lf histogram -f nBCs --nbins 5 --auto
    echo
    
    echo -e "${sample}\tHistogram of number of cells per BCs"
    cat ${OUTDIR}/${sample}.barcode.counts.byCell.filtered.txt | mlr --headerless-csv-output --tsv cut -f BC | uniq -c | mlr --fs space --repifs --ocsv --ofs "\t" --ofmt %.1lf histogram -f 1 --nbins 5 --auto -o nCells
    echo
fi


echo
echo "Saturation curves"
#BUGBUG does not consider UMIs for thresholding
#NB hardcoded limit of 120M reads
for i in `echo {10000,50000,100000,250000,500000,1000000,2000000,3000000,4000000,5000000,10000000,15000000,20000000,25000000,30000000,35000000,40000000,45000000,50000000,55000000,60000000,70000000,80000000,90000000,100000000,110000000,120000000} | tr ' ' '\n' | gawk -v subset=${numTotalReads} '$1<=subset'`; do 
    if [ "${minUMILength}" -ge 5 ]; then
        saturation=$(zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | shuf -n $i | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $3}' | uniq | sort | uniq | awk -F "\t" '{print $1}' | sort | uniq -c | awk -v minReadCutoff=${minReadCutoff} '{if ($1>=minReadCutoff) print $2}' | wc -l)
    else
        saturation=$(zcat -f ${OUTDIR}/${sample}.barcodes.txt.gz | shuf -n $i | awk -F "\t" '{print $1}' | sort | uniq -c | awk -v minReadCutoff=${minReadCutoff} '{if ($1>=minReadCutoff) print $2}' | wc -l)
    fi
    echo "$i $saturation minreads${minReadCutoff}" 
done > ${OUTDIR}/${sample}.Saturation.minReads_${minReadCutoff}.txt


echo
echo "Done!!!"
date
