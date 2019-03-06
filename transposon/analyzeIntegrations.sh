#!/bin/bash
set -eu -o pipefail

src=/vol/mauranolab/transposon/src

sample=$1


#TODO inconsistently applied
minReadCutoff=2

OUTDIR=${sample}

hotspotfile=/vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38-final/K562-DS9764.hg38.fdr0.01.pks.bed



if [ ! -s "$OUTDIR/${sample}.barcodes.txt" ]; then
    echo "analyzeIntegrations.sh ERROR: barcode input file $OUTDIR/${sample}.barcodes.txt does not exist!"
    exit 1
fi

if [ ! -s "$OUTDIR/${sample}.bam" ]; then
    echo "analyzeIntegrations.sh ERROR: mapped reads input file $OUTDIR/${sample}.bam does not exist!"
    exit 2
fi


echo "Analyzing data for ${sample}"
echo "Analyzing read mapping (minReadCutoff=${minReadCutoff})"


echo
echo "SAMtools statistics for sample ${sample}"
samtools flagstat $OUTDIR/${sample}.bam | tee $TMPDIR/${sample}.flagstat.txt

echo
echo -n -e "${sample}\tTotal PF tags\t"
cat $TMPDIR/${sample}.flagstat.txt | grep "in total" | awk '{print $1}'


#Get the reads from the bam since we don't save the trimmed fastq
samtools view $OUTDIR/${sample}.bam | cut -f10 | cut -f1 | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | shuf -n 1000000 | awk -F "\t" 'BEGIN {OFS="\t"} {print ">id-" NR; print}' |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} genomic" --stacks-per-line 100 > $TMPDIR/${sample}.genomic.eps
convert $TMPDIR/${sample}.genomic.eps ${OUTDIR}/${sample}.genomic.png


echo
echo "Merge mapping and barcodes"
samflags="-F 512"
#Subtract additional 1 bp from reads on + strand so that the coordinates represent 1 nt to the left of the insertion site (i.e. for A^T, the coordinates point to T)
samtools view ${samflags} $OUTDIR/${sample}.bam | awk -F "\t" 'BEGIN {OFS="\t"} { \
    readlength = length($10); \
    if (and($2, 16)) { \
        strand="-"; \
        chromStart=$4+readlength-1; \
    } else { \
        strand="+"; \
        chromStart=$4-1; \
    } \
    chromEnd=chromStart+1; \
    print $3, chromStart, chromEnd, $1, "id-" NR, strand; \
}' | sort-bed - > $TMPDIR/${sample}.coords.bed

echo -e -n "Number of tags passing all filters\t"
cat $TMPDIR/${sample}.coords.bed | wc -l


echo
echo "read lengths:"
samtools view ${samflags} $OUTDIR/${sample}.bam | cut -f10 | awk 'BEGIN {ORS=", "} {lengths[length($0)]++} END {for (l in lengths) {print l " (" lengths[l] ")" }}' | perl -pe 's/,$//g;'


cat $OUTDIR/${sample}.barcodes.txt | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | 
sort -k2,2 > $TMPDIR/${sample}.barcodes.txt
cat $TMPDIR/${sample}.coords.bed | sort -k4,4 | join -1 4 -2 2 - $TMPDIR/${sample}.barcodes.txt | awk 'BEGIN {OFS="\t"} {print $2, $3, $4, $1, $6, $7, $8}' | sort-bed - > $TMPDIR/${sample}.barcodes.readnames.coords.raw.bed
#columns: chrom, start, end, readID, strand, BC seq, UMI
#NB strand in $5


echo -e -n "${sample}\tNumber of tags passing all filters and having barcodes assigned\t"
cat $TMPDIR/${sample}.barcodes.readnames.coords.raw.bed | wc -l


echo
echo "Histogram of barcode reads before coordinate-based deduping"
cat $TMPDIR/${sample}.barcodes.readnames.coords.raw.bed | cut -f6 | sort | uniq -c | awk '{print $1}' | awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g


echo
echo "Coordinate-based deduping"
cat $TMPDIR/${sample}.barcodes.readnames.coords.raw.bed | cut -f1-7 |
awk -v maxstep=5 -F "\t" 'BEGIN {OFS="\t"; groupnum=1} { \
    if( NR>1 && (last[1] != $1 || $2 - last[2] > maxstep || $3 - last[3] > maxstep || last[5] != $5 )) { \
        #5bp between the insertions group together\
        groupnum = groupnum+1; \
    } \
    $(NF+1) = groupnum; \
    print; \
    split($0, last); \
}' |
#NB no need to sort
${src}/AdjacencyDeDup.py --col 6 --groupcol 8 -o - - |
#BUGBUG occasionally spits out null barcodes , e.g. from FCHHLLWBGX7/BS01481A-RDL_20180722_K562_pMH022_T0098_GFPpos_iPCR
#chr12   109074738       109074739       NB501831:111:HHLLWBGX7:4:13407:1923:8064        -
awk -F "\t" 'BEGIN {OFS="\t"} $6!=""' |
#I think column 7 is the UMI, if any
cut -f1-7 > $OUTDIR/${sample}.barcodes.readnames.coords.bed
#NB Barcodes are now out of sync with those from analyzeBCcounts.sh


echo
echo "Histogram of barcode reads after coordinate-based deduping"
cat $OUTDIR/${sample}.barcodes.readnames.coords.bed | cut -f6 | sort | uniq -c | awk '{print $1}' | awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g


minUMIlength=5
if [[ `awk -v minUMIlength=${minUMIlength} -F "\t" 'BEGIN {OFS="\t"} length($3) >= minUMIlength {found=1} END {print found}' $OUTDIR/${sample}.barcodes.txt` == 1 ]]; then
    echo
    echo "Coordinate-based UMI deduping"
    
    cat $OUTDIR/${sample}.barcodes.readnames.coords.bed | cut -f1-2,5-7 | sort -k1,1 -k2,2g -k4,4 -k5,5 | uniq -c |
    awk 'BEGIN {OFS="\t"} {print $2, $3, $5, $1, $4}' > $TMPDIR/${sample}.coords.collapsedUMI.bed
    #(NB chromEnd is omitted)
    
    echo
    echo -n -e "${sample}\tNumber of unique coordinates+barcodes+UMI\t"
    #TODO move to extract.py
    cat $TMPDIR/${sample}.coords.collapsedUMI.bed | wc -l
else
    echo
    echo "No UMIs found"
    cat $OUTDIR/${sample}.barcodes.readnames.coords.bed | awk 'BEGIN {OFS="\t"} {print $1, $2, $6, 1, $5}' | sort -k1,1 -k2,2g -k3,3 > $TMPDIR/${sample}.coords.collapsedUMI.bed
fi
#columns: chrom, start, BC seq, count, strand


echo "Collapsing nearby insertions"
echo -n -e "${sample}\tNumber of insertions before collapsing nearby ones\t"
cat $TMPDIR/${sample}.coords.collapsedUMI.bed | wc -l

cut -f1-3,5 $TMPDIR/${sample}.coords.collapsedUMI.bed |
uniq -c |
awk 'BEGIN {OFS="\t"} {print $2, $3, $3+1, $4, $1, $5}' |
#columns: chrom, start, end, BC seq, count, strand
#Combine subsequent lines having the same barcode at slightly different sites
sort -k4,4 -k1,1 -k2,2g |
awk -v maxstep=5 -F "\t" 'BEGIN {OFS="\t"} \
NR==1 {split($0, last)} \
NR>1 { \
    if(last[1] != $1 || $2 - last[2] > maxstep || $3 - last[3] > maxstep || last[4] != $4 || last[6] != $6) { 
        #not the same site or BC, so print \
        print last[1], last[2], last[3], last[4], last[5], last[6]; \
        split($0, last); \
    } else { #the same BC \
        if($5 > last[5]) { \
            #Use the coords from this site since it has more reads \
            last[2]=$2; last[3]=$3; \
        } \
        last[5] += $5; #add the counts \
    } \
} \
END {print last[1], last[2], last[3], last[4], last[5], last[6]}' |
sort-bed - |
awk -v minReads=1 -F "\t" 'BEGIN {OFS="\t"} $5>=minReads' > $OUTDIR/${sample}.barcodes.coords.bed
#columns: chrom, start, end, readname, strand, BC seq, count (coords are of integration site)

echo -n -e "${sample}\tNumber of insertions after collapsing nearby ones\t"
cat $OUTDIR/${sample}.barcodes.coords.bed | wc -l


#TODO make second dedup pass after collapsing nearby sites? Or perhaps change group column above to include BCs at sites within maxstep?


echo "Generating UCSC track"
awk -v sample=${sample} -F "\t" 'BEGIN {OFS="\t"; print "track name=" sample " description=" sample "-integrations"} {print}' $OUTDIR/${sample}.barcodes.coords.bed > $OUTDIR/${sample}.barcodes.coords.ucsc.bed
projectdir=`pwd | perl -pe 's/^\/vol\/mauranolab\/transposon\///g;'`
UCSCbaseURL="https://mauranolab@cascade.isg.med.nyu.edu/~mauram01/transposon/${projectdir}/${OUTDIR}"
echo "${UCSCbaseURL}/${sample}.barcodes.coords.ucsc.bed"


cat $OUTDIR/${sample}.barcodes.coords.bed | awk -v minReadCutoff=${minReadCutoff} -F "\t" 'BEGIN {OFS="\t"} $5>minReadCutoff' | awk -F "\t" 'BEGIN {OFS="\t"} {$4="."; $5=0; print}' | uniq > $OUTDIR/${sample}.uniqcoords.bed


echo
echo "Integration sites by chrom"
cat $OUTDIR/${sample}.uniqcoords.bed | cut -f1 | sort -g | uniq -c | sort -k1,1g

echo
echo -e -n "${sample}\tTotal uniq sites\t"
cat $OUTDIR/${sample}.uniqcoords.bed | wc -l

echo -e -n "${sample}\tTotal uniq sites (within 5 bp, ignoring strand)\t"
cat $OUTDIR/${sample}.uniqcoords.bed | bedops --range 5 -m - | uniq | wc -l


echo
echo "Histogram of number of reads per barcode"
cat $OUTDIR/${sample}.barcodes.coords.bed | cut -f5 | awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g


echo
echo -n -e "${sample}\tProportion of unique insertion sites at TA\t"
cat $OUTDIR/${sample}.uniqcoords.bed | awk -F "\t" 'BEGIN {OFS="\t"} {$2-=2; $3-=1; print}' |  /home/mauram01/bin/bed2fasta.pl - /vol/isg/annotation/fasta/hg38 2>/dev/null | grep -v -e "^>" | tr '[a-z]' '[A-Z]' | awk -F "\t" 'BEGIN {OFS="\t"; count=0} $0=="TA" {count+=1} END {print count/NR}'


echo
echo "Histogram of number of insertion sites per barcode"
cat $OUTDIR/${sample}.barcodes.coords.bed | awk -v minReadCutoff=${minReadCutoff} -F "\t" 'BEGIN {OFS="\t"} $5>minReadCutoff' | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $6, $4}' | sort -k5,5 | cut -f5 | sort -g | uniq -c | sort -k1,1g | awk '{print $1}' | 
awk -v cutoff=2 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g


echo
echo "Histogram of number of insertions between two neighboring DNase sites"
uniqueIntervals=$(tail -n +2 ${hotspotfile} | paste ${hotspotfile} - | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $4, $5}' | sed \$d | awk -F "\t" 'BEGIN {OFS="\t"} {if ($1==$3) print $1, $2, $4}' | bedtools intersect -wa -a - -b $OUTDIR/${sample}.barcodes.coords.bed | sort | uniq | wc -l | awk '{print $1}')
allIntervals=$(wc -l ${hotspotfile} | awk '{print $1}')
zeroInsertions=`echo $allIntervals-$uniqueIntervals | bc -l`
echo " "$zeroInsertions 0

tail -n +2 ${hotspotfile} |
paste ${hotspotfile} - |
awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $4, $5}' |
sed \$d |
awk -F "\t" 'BEGIN {OFS="\t"} {if ($1==$3) print $1, $2, $4}' |
bedtools intersect -wa -a - -b $OUTDIR/${sample}.barcodes.coords.bed |
sort | uniq -c | awk '{print $1}' |
awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' |
sort | uniq -c | sort -k2,2g


echo
echo "doing DistToTSS"
closest-features --closest --dist --no-ref $OUTDIR/${sample}.uniqcoords.bed /vol/isg/annotation/bed/hg38/refseq_gene/refGene.CombinedTxStarts.bed | awk -F "|" 'BEGIN {OFS="\t"} function abs(value) {return (value<0?-value:value);} {print $2}' > $OUTDIR/DistToTSS.txt

R --quiet --no-save << EOF
#Work around "unable to start device PNG" on ISG cluster
options(bitmapType="cairo") 

filename <- "$OUTDIR/DistToTSS.txt"
cat("Reading", filename, "\n")
data <- read(filename)
if(is.null(data)) {
    quit(0, save = "no", status = 0)
}

colnames(data) <- c("DistToTSS")

data\$DistToTSS.bin <- cut(data\$DistToTSS, breaks=c(-10^7, -10^6, -10^5, -10^4, -2500, 0, 5000, 10^4, 10^5, 10^6, 10^7), right=F, include.lowest=T, labels=c("-10 Mb", "-1 Mb", "-100 kb", "-10 kb", "-2.5 kb", "+5 kb", "+10 kb", "+100 kb", "+1 Mb", "+10 Mb"))

library(ggplot2)
old <- theme_set(theme_classic(base_size=10)) #png
old <- theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

png(file="$OUTDIR/DistToTSS.png", width=400, height=400)
ggplot(data, aes(x=DistToTSS.bin)) +
geom_bar(color="black", size=0.25) +
theme_classic()+
theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(size=16),axis.title.y=element_text(size=20),axis.text.y=element_text(size=16))
dev.off()
EOF


echo
echo "doing DisttoDpn"
cat $OUTDIR/${sample}.uniqcoords.bed | closest-features --dist  --delim '\t' - /vol/isg/annotation/bed/hg38/REsites/Dpn/Dpn.bed | awk -F "\t" 'BEGIN {OFS="\t"} function abs(value) {return (value<0?-value:value);} {if ($6=="-") print $1, $2, $3, $10, abs($13); else if ($6=="+") print $1, $2, $3, $17, $20}' > $OUTDIR/DistDpn.bed

R --quiet --no-save << EOF
#Work around "unable to start device PNG" on ISG cluster
options(bitmapType="cairo") 

filename <- "$OUTDIR/DistDpn.bed"
cat("Reading", filename, "\n")
data <- read(filename)
if(is.null(data)) {
    quit(0, save = "no", status = 0)
}
colnames(data)[5] <- c("DinsToDpn")

library(ggplot2)
png(file="$OUTDIR/DistToDpnSite.png", width=400, height=400)
ggplot(data, aes(x=abs(DinsToDpn))) +
geom_histogram(color="black", size=0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(size=16),axis.title.y=element_text(size=20),axis.text.y=element_text(size=16))+
xlim(c(0,1500))+
ggtitle("${sample}")+
xlab('Distance to Dpn site')
dev.off()
EOF


echo
echo "doing DisttoMspI"
cat $OUTDIR/${sample}.uniqcoords.bed | closest-features --dist  --delim '\t' - /vol/isg/annotation/bed/hg38/REsites/MspI/MspI.bed | awk -F "\t" 'BEGIN {OFS="\t"} function abs(value) {return (value<0?-value:value);} {if ($6=="-") print $1, $2, $3, $10, abs($13); else if ($6=="+") print $1, $2, $3, $17, $20}' > $OUTDIR/DistDpn.bed

R --quiet --no-save << EOF
#Work around "unable to start device PNG" on ISG cluster
options(bitmapType="cairo") 

filename <- "$OUTDIR/DistDpn.bed"
cat("Reading", filename, "\n")
data <- read(filename)
if(is.null(data)) {
    quit(0, save = "no", status = 0)
}
colnames(data)[5] <- c("DinsToDpn")

library(ggplot2)
png(file="$OUTDIR/DistToDpnSite.png", width=400, height=400)
ggplot(data, aes(x=abs(DinsToDpn))) +
geom_histogram(color="black", size=0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(size=16),axis.title.y=element_text(size=20),axis.text.y=element_text(size=16))+
xlim(c(0,1500))+
ggtitle("${sample}")+
xlab('Distance to Dpn site')
dev.off()
EOF

echo
echo "doing DistToDNase"
cat $OUTDIR/${sample}.uniqcoords.bed | closest-features --dist  --delim '\t' --closest -  ${hotspotfile} | awk -F "\t" 'BEGIN {OFS="\t"} function abs(value) {return (value<0?-value:value);} {print $1, $2, $3, abs($10)}' > $OUTDIR/DistToDNase.bed

R --quiet --no-save << EOF
#Work around "unable to start device PNG" on ISG cluster
options(bitmapType="cairo") 

filename <- "$OUTDIR/DistToDNase.bed"
cat("Reading", filename, "\n")
data <- read(filename)
if(is.null(data)) {
    quit(0, save = "no", status = 0)
}
colnames(data)[4] <- c("DinsToDpn")

#data\$DistToTSS.bin <- cut(data\$DistToTSS, breaks=c(-10^7, -10^6, -10^5, -10^4, -2500, 0, 5000, 10^4, 10^5, 10^6, 10^7), right=F, include.lowest=T, labels=c("-10 Mb", "-1 Mb", "-100 kb", "-10 kb", "-2.5 kb", "+5 kb", "+10 kb", "+100 kb", "+1 Mb", "+10 Mb"))

library(ggplot2)
png(file="$OUTDIR/DistToDNase.png", width=400, height=400)
ggplot(data, aes(x=abs(DinsToDpn+1))) +
geom_histogram(color="black", size=0.25)+ 
theme_classic()+
theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(size=16),axis.title.y=element_text(size=20),axis.text.y=element_text(size=16))+
#xlim(c(0,1500))+
ggtitle("${sample}")+
xlab('Distance to DNase site')+
scale_x_log10(  breaks = scales::trans_breaks("log10", function(x) 10^x),  labels = scales::trans_format("log10", scales::math_format(10^.x)))+
annotation_logticks(sides = "b",long = unit(0.2, "cm"))    
dev.off()
EOF


echo
echo "doing GenicLocation"
gcat /vol/isg/annotation/bed/hg38/refseq_paint/refGene_hg38.bed6.starch | grep -v promoter > $TMPDIR/refGene_hg38.bed6

bedmap --delim '\t' --echo-map-id $OUTDIR/${sample}.uniqcoords.bed $TMPDIR/refGene_hg38.bed6 | awk -F "\t" 'BEGIN {OFS="\t"} {col=NF} $col~/coding/ {$col="coding"} $col~/TxS-1000/ {$col="TxS-1000"} $col~/TxStart-10kb/ {$col="TxStart-10kb"} $col~/3.UTR/ {$col="3UTR"} $col~/5.UTR/ {$col="5UTR"} $col~/intron/ {$col="intron"} $col=="" || $col~/3.proximal/ {$col="intergenic"} {print;}' | sort -g | uniq -c | sort -k1,1g


echo
echo "Done!!!"
date
