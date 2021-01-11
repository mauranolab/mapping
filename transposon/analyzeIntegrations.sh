#!/bin/bash
set -eu -o pipefail

src=/vol/mauranolab/mapped/src/transposon


getcolor () {
    local trackcolor
    
    trackcolor="51,128,195" #blue
    if [[ "$1" =~ A1fw ]] || [[ "$1" =~ C1fw ]] || [[ "$1" =~ A1rv ]] || [[ "$1" =~ C1rv ]]; then
        trackcolor="120,88,165" #purple
    elif [[ "$1" =~ HS2 ]]; then
        trackcolor="238,54,36" #red
    fi
    
    echo ${trackcolor}
}


sample=$1

trackcolor=$(getcolor $sample)

minReadCutoff=2

OUTDIR=${sample}


if [ ! -s "$OUTDIR/${sample}.barcodes.txt.gz" ]; then
    echo "analyzeIntegrations.sh ERROR: barcode input file $OUTDIR/${sample}.barcodes.txt does not exist!"
    exit 1
fi

if [ ! -s "$OUTDIR/${sample}.bam" ]; then
    echo "analyzeIntegrations.sh ERROR: mapped reads input file $OUTDIR/${sample}.bam does not exist!"
    exit 2
fi

# BUG: hardcoded and dependent of genome
curGenome="mm10"
case "${curGenome}" in
hg38_noalt)
    hotspotfile=/vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch;
    chromsizes=/vol/isg/annotation/fasta/hg38_noalt/hg38_noalt.chrom.sizes;;
mm10)
    # FIX hard-coded references to humam genome reference
    hotspotfile=/vol/isg/encode/dnase/mapped/K562-DS9764/hotspots/K562-DS9764.hg38_noalt-final/K562-DS9764.hg38_noalt.fdr0.01.pks.starch;
    chromsizes=/vol/isg/annotation/fasta/mm10/mm10.chrom.sizes;;
*)
    echo "Don't recognize genome ${curGenome}";
    exit 3;;
esac

## Ignore R1 in paired without excluding single-end reads
samflags="-F 64 -F 512"

echo "Analyzing data for ${sample}"
echo "Analyzing read mapping (minReadCutoff=${minReadCutoff})"
date


echo
echo "SAMtools statistics for sample ${sample}"
samtools flagstat $OUTDIR/${sample}.bam | tee $TMPDIR/${sample}.flagstat.txt

#Reproduce some of the statistics that filter_reads.py prints, but do it here for the full dataset. Note we use the 512 flag to hopefully speed these up before. I suppose we could just run filter_reads.py here without much performance penalty
echo -n -e "${sample}\tTotal PF reads\t"
cat $TMPDIR/${sample}.flagstat.txt | grep "in total" | awk '{print $1}'

echo -n -e "${sample}\tTotal reads mapping to pSB or pTR\t"
samtools view -c -f 512 $OUTDIR/${sample}.bam pTR pSB

echo -n -e "${sample}\tTotal mapped reads MAPQ<10\t"
samtools view -F 4 -f 512 $OUTDIR/${sample}.bam | awk -F "\t" 'BEGIN {OFS="\t"} $5<10' | wc -l

echo -n -e "${sample}\tTotal unmapped reads\t"
samtools view -f 516 -c $OUTDIR/${sample}.bam

echo -n -e "${sample}\tTotal reads mapped to unscaffolded contigs\t"
samtools view -f 512 $OUTDIR/${sample}.bam | cut -f3 | awk '$0 ~ /hap|random|^chrUn_|_alt$|scaffold|^C\d+/' | wc -l

echo
readlengths=`samtools view ${samflags} $OUTDIR/${sample}.bam | cut -f10 | awk 'BEGIN {ORS=", "} {lengths[length($0)]++} END {for (l in lengths) {print l " (" lengths[l] ")" }}' | perl -pe 's/, $//g;'`
echo -e "${sample}\tRead lengths (number of reads)\t${readlengths}"
minReadLength=`echo "${readlengths}" | perl -pe 's/ \([0-9]+\)//g;' -e 's/, /\n/g;' | sort -n | awk 'NR==1'`
echo -e "${sample}\tMinimum read length\t${minReadLength}"

hasRead=$(samtools view -c ${samflags} $OUTDIR/${sample}.bam)
if [ "$hasRead" -eq 0 ]; then
    echo "No reads passing all filters. Exiting earlier"
    echo
    echo "Done!!!"
    date
    exit 0
fi

#Get the reads from the bam since we don't save the trimmed fastq
#Need to trim in case not all the same length
samtools view ${samflags} $OUTDIR/${sample}.bam | cut -f10 | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | shuf -n 1000000 | awk -F "\t" -v trim=${minReadLength} 'BEGIN {OFS="\t"} {print substr($0, 0, trim)}' | awk -F "\t" 'BEGIN {OFS="\t"} {print ">id-" NR; print}' |
weblogo --datatype fasta --color-scheme 'classic' --size large --sequence-type dna --units probability --title "${sample} genomic" --stacks-per-line 100 > $TMPDIR/${sample}.genomic.eps
convert $TMPDIR/${sample}.genomic.eps ${OUTDIR}/${sample}.genomic.png

echo
echo "Extracting read coordinates"
date
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

echo -e -n "${sample}\tNumber of reads passing all filters\t"
cat $TMPDIR/${sample}.coords.bed | wc -l

echo -n -e "${sample}\tNumber of pairedreads mapped passing all filters\t"
samtools view -c ${samflags} -f 1 $OUTDIR/${sample}.bam

zcat -f $OUTDIR/${sample}.barcodes.txt.gz | awk -F "\t" 'BEGIN {OFS="\t"} $1!=""' | sort -k2,2 > $TMPDIR/${sample}.barcodes.txt
cat $TMPDIR/${sample}.coords.bed | sort -k4,4 | join -1 4 -2 2 - $TMPDIR/${sample}.barcodes.txt | awk 'BEGIN {OFS="\t"} {print $2, $3, $4, $1, $6, $7, $8}' | sort-bed - > $TMPDIR/${sample}.barcodes.readnames.coords.raw.bed
#columns: chrom, start, end, readID, strand, BC seq, UMI
#NB strand in $5
echo -e -n "${sample}\tNumber of reads passing all filters and having barcodes assigned\t"
cat $TMPDIR/${sample}.barcodes.readnames.coords.raw.bed | wc -l


echo
echo "Histogram of barcode reads before coordinate-based deduping"
cat $TMPDIR/${sample}.barcodes.readnames.coords.raw.bed | cut -f6 | sort | uniq -c | awk '{print $1}' | awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g


echo
echo "Coordinate-based deduping"
date
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
${src}/AdjacencyDeDup.py --col 6 --groupcols 8 -o - - |
#BUGBUG occasionally spits out null barcodes , e.g. from FCHHLLWBGX7/BS01481A-RDL_20180722_K562_pMH022_T0098_GFPpos_iPCR
#chr12   109074738       109074739       NB501831:111:HHLLWBGX7:4:13407:1923:8064        -
awk -F "\t" 'BEGIN {OFS="\t"} $6!=""' |
#I think column 7 is the UMI, if any
cut -f1-7 > $OUTDIR/${sample}.barcodes.readnames.coords.bed
#NB Barcodes are now out of sync with those from analyzeBCcounts.sh


echo
echo "Histogram of barcode reads after coordinate-based deduping"
cat $OUTDIR/${sample}.barcodes.readnames.coords.bed | cut -f6 | sort | uniq -c | awk '{print $1}' | awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g


#Only dedup UMIs if the length is over 5
#Go through entire file with awk despite only looking at first line so zcat terminates properly
minUMILength=`zcat -f $OUTDIR/${sample}.barcodes.txt.gz | awk -F "\t" 'BEGIN {OFS="\t"} $1!="" {print length($3)}' | uniq | sort -n | awk 'NR==1'`
if [[ "${minUMILength}" -ge "5" ]]; then
    #never implemented deduping
    echo
    echo "Analyzing UMIs"
    
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
#columns: chrom, start, BC seq, UMI count, strand


echo "Collapsing nearby insertions"
cut -f1-3,5 $TMPDIR/${sample}.coords.collapsedUMI.bed |
uniq -c |
awk 'BEGIN {OFS="\t"} {print $2, $3, $3+1, $4, $1, $5}' > $TMPDIR/${sample}.barcodes.coords.bed

echo -n -e "${sample}\tNumber of BC+insertions before collapsing nearby ones\t"
cat $TMPDIR/${sample}.barcodes.coords.bed | wc -l

echo -n -e "${sample}\tNumber of unique insertion sites before collapsing nearby ones\t"
cat $TMPDIR/${sample}.barcodes.coords.bed | awk -F "\t" 'BEGIN {OFS="\t"} {$4="."; $5=0; print}' | uniq | wc -l


#columns: chrom, start, end, BC seq, count, strand
#Combine subsequent lines having the same barcode at slightly different sites
if [ "$(cat $TMPDIR/${sample}.barcodes.coords.bed | wc -l)" -gt 0 ];
then
    cat $TMPDIR/${sample}.barcodes.coords.bed | sort -k4,4 -k1,1 -k2,2g |
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
    sort-bed - > $TMPDIR/${sample}.barcodes.coords.collapsed.bed
else
    cp $TMPDIR/${sample}.barcodes.coords.bed $TMPDIR/${sample}.barcodes.coords.collapsed.bed
fi

echo -n -e "${sample}\tNumber of BC+insertions after collapsing nearby ones\t"
cat $TMPDIR/${sample}.barcodes.coords.collapsed.bed | wc -l

echo
echo "Histogram of number of reads per barcode"
cat $TMPDIR/${sample}.barcodes.coords.collapsed.bed | cut -f5 | awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g


#Drop BCs which have less than 10% of the reads at each site
cat $TMPDIR/${sample}.barcodes.coords.collapsed.bed | bedmap --exact --delim '\t' --prec 0 --echo --sum - | awk -v minPropReadsAtSite=0.1 -F "\t" 'BEGIN {OFS="\t"} $5/$7>minPropReadsAtSite' | cut -f1-6 > $TMPDIR/${sample}.barcodes.coords.minPropReadsAtSite.bed
#columns: chrom, start, end, readname, strand, BC seq, count (coords are of integration site)

echo
echo -n -e "${sample}\tNumber of BC+insertions passing minPropReadsAtSite cutoff\t"
cat $TMPDIR/${sample}.barcodes.coords.minPropReadsAtSite.bed | wc -l

echo
echo "Histogram of number of reads per barcode"
cat $TMPDIR/${sample}.barcodes.coords.minPropReadsAtSite.bed | cut -f5 | awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g


#TODO make second dedup pass after collapsing nearby sites? Or perhaps change group column above to include BCs at sites within maxstep?


#Enforce minimum read cutoff
cat $TMPDIR/${sample}.barcodes.coords.minPropReadsAtSite.bed | awk -v minReadCutoff=${minReadCutoff} -F "\t" 'BEGIN {OFS="\t"} $5>minReadCutoff' > $TMPDIR/${sample}.barcodes.coords.minReadCutoff.bed


echo
echo "Histogram of number of insertion sites per barcode"
cat $TMPDIR/${sample}.barcodes.coords.minReadCutoff.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $6, $4}' | sort -k5,5 | cut -f5 | sort -g | uniq -c | sort -k1,1g | awk '{print $1}' |
awk -v cutoff=2 '{if($0>=cutoff) {print cutoff "+"} else {print}}' | sort -g | uniq -c | sort -k2,2g


#Identify iPCR insertions with more than one location
cat $TMPDIR/${sample}.barcodes.coords.minReadCutoff.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $4}' | sort | uniq -c | sort -nk1 | awk '$1==1 {print $2}' > $TMPDIR/${sample}.singleIns.txt

# save BC list of insertions with more than one location
cat $TMPDIR/${sample}.barcodes.coords.minReadCutoff.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $4}' | sort | uniq -c | sort -nk1 | awk '$1>1 {print $2}' > $OUTDIR/${sample}.multiIns.txt

#Identify insertion sites with more than one BC
cat $TMPDIR/${sample}.barcodes.coords.minReadCutoff.bed | awk -F "\t" 'BEGIN {OFS="\t"} {$4="."; $5=0; print}' | uniq -c | awk 'BEGIN {OFS="\t"} $1==1 {print $2, $3, $4, $5, $6, $7}' | sort-bed - > $TMPDIR/${sample}.singleBC.bed

if [ $(wc -l $TMPDIR/${sample}.singleIns.txt) -gt 0 ]; then
    cat $TMPDIR/${sample}.barcodes.coords.minReadCutoff.bed |
    #Remove BCs with more than one location
    awk -F "\t" 'BEGIN {OFS="\t"} NR==FNR{a[$1];next} ($4) in a' $TMPDIR/${sample}.singleIns.txt - |
    #Remove insertion sites with more than one BC
    bedmap --delim "|" --multidelim "|" --bp-ovr 1 --skip-unmapped --echo --echo-map - $TMPDIR/${sample}.singleBC.bed | awk -F "|" 'BEGIN {OFS="\t"} {split($0, main, "\t"); for(i=2; i<=NF; i++) {split($0, singleBC, "\t"); if(main[6]==singleBC[6]) {print $1; next}}}' > $OUTDIR/${sample}.barcodes.coords.bed
else
    cp $TMPDIR/${sample}.barcodes.coords.minReadCutoff.bed $OUTDIR/${sample}.barcodes.coords.bed
fi

#NB retains strand so a few sites are represented twice
cat $OUTDIR/${sample}.barcodes.coords.bed | awk -F "\t" 'BEGIN {OFS="\t"} {$4="."; $5=0; print}' | uniq | sort-bed - > $OUTDIR/${sample}.uniqcoords.bed


echo
echo "Generating UCSC track"
projectdir=`pwd | perl -pe 's/^\/vol\/mauranolab\/transposon\///g;'`
UCSCbaseURL="https://cascade.isg.med.nyu.edu/~mauram01/transposon/${projectdir}/${OUTDIR}"
echo "track name=${sample} description=${sample}-integrations db=hg38 color=$trackcolor"
echo "${UCSCbaseURL}/${sample}.barcodes.coords.bed"


###Analysis
echo
echo -n -e "${sample}\tNumber of BC+insertions after minimum read cutoff\t"
cat $OUTDIR/${sample}.barcodes.coords.bed | wc -l


echo
echo "Unique integration sites by chrom"
cat $OUTDIR/${sample}.uniqcoords.bed | cut -f1 | sort -g | uniq -c | sort -k1,1g

echo
echo -e -n "${sample}\tNumber of unique insertion sites\t"
cat $OUTDIR/${sample}.uniqcoords.bed | wc -l

echo -e -n "${sample}\tNumber of unique insertion sites (within 5 bp, ignoring strand)\t"
cat $OUTDIR/${sample}.uniqcoords.bed | bedops --range 5 -m - | uniq | wc -l


echo -n -e "${sample}\tProportion of unique insertion sites at TA\t"
if [ $(cat $OUTDIR/${sample}.uniqcoords.bed | wc -l) -gt 0 ];
then
    cat $OUTDIR/${sample}.uniqcoords.bed | awk -F "\t" 'BEGIN {OFS="\t"} {$2-=2; $3-=1; print}' |  /home/mauram01/bin/bed2fasta.pl - /vol/isg/annotation/fasta/hg38 2>/dev/null | grep -v -e "^>" | tr '[a-z]' '[A-Z]' | awk -F "\t" 'BEGIN {OFS="\t"; count=0} $0=="TA" {count+=1} END {print count/NR}'
else
    echo 0
fi

echo
echo "Histogram of number of insertions between two neighboring DNase sites"
unstarch ${hotspotfile} > $TMPDIR/hotspotfile.bed
uniqueIntervals=$(tail -n +2 $TMPDIR/hotspotfile.bed | paste $TMPDIR/hotspotfile.bed - | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $4, $5}' | sed \$d | awk -F "\t" 'BEGIN {OFS="\t"} {if ($1==$3) print $1, $2, $4}' | bedtools intersect -wa -a - -b $OUTDIR/${sample}.barcodes.coords.bed | sort | uniq | wc -l | awk '{print $1}')
allIntervals=$(wc -l $TMPDIR/hotspotfile.bed | awk '{print $1}')
zeroInsertions=`echo $allIntervals-$uniqueIntervals | bc -l`
echo " "$zeroInsertions 0

tail -n +2 $TMPDIR/hotspotfile.bed |
paste $TMPDIR/hotspotfile.bed - |
awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $4, $5}' |
sed \$d |
awk -F "\t" 'BEGIN {OFS="\t"} {if ($1==$3) print $1, $2, $4}' |
bedtools intersect -wa -a - -b $OUTDIR/${sample}.barcodes.coords.bed |
sort | uniq -c | awk '{print $1}' |
awk -v cutoff=10 '{if($0>=cutoff) {print cutoff "+"} else {print}}' |
sort | uniq -c | sort -k2,2g


echo
echo "Density track"

##Prep CNV track, set regions without data to baseline of 0 (which is really 2N)
#cat ${chromsizes} | 
#awk -F "\t" 'BEGIN {OFS="\t"} $1!="chrM" && $1!="chrEBV"' |
#awk -F "\t" '{OFS="\t"; print $1, 0, $2}' | sort-bed - |
#bedops -d - /vol/mauranolab/publicdata/K562_copynumber/CCLE.CNV.hg38.bed |
#awk -F "\t" '{OFS="\t"; print $1, $2, $3, "..", -1}' |
#bedops -u - /vol/mauranolab/publicdata/K562_copynumber/CCLE.CNV.hg38.bed |
#awk -F "\t" 'BEGIN {OFS="\t"} {$NF=2^($NF+1); print}' > $TMPDIR/K562.CCLE.CNV.hg38.bed

#score is insertion count per Mb bin. Disabled code normalizes to N=2 by CCLE CNV map
#TODO normalize to count of unique insertions?
cat ${chromsizes} | 
awk -F "\t" 'BEGIN {OFS="\t"} $1!="chrM" && $1!="chrEBV"' |
awk -F "\t" '{OFS="\t"; print $1, 0, $2}' | sort-bed - | cut -f1,3 | awk -v step=100000 -v binwidth=1000000 'BEGIN {OFS="\t"} {for(i=0; i<=$2-binwidth; i+=step) {print $1, i, i+binwidth, "."} }' | 
#--faster is ok since we are dealing with bins and read starts
bedmap --faster --delim "\t" --bp-ovr 1 --echo --count - $OUTDIR/${sample}.uniqcoords.bed |
#resize intervals down from full bin width to step size
awk -v step=100000 -v binwidth=1000000 -F "\t" 'BEGIN {OFS="\t"} {offset=(binwidth-step)/2 ; $2+=offset; $3-=offset; print}' |

##Take a distance-weighted average of the copy number over the interval and normalize to N=2
#bedmap --delim "|" --multidelim "|" --echo --echo-map - $TMPDIR/K562.CCLE.CNV.hg38.bed |
#awk -F "|" 'BEGIN {OFS="\t"} { \
#    split($1, dens, "\t"); avg=0; \
#    for(i=2; i<=NF; i++) { \
#        split($i, curCN, "\t");\
#        if(curCN[2] < dens[2]) {curCN[2]=dens[2]} \
#        if(curCN[3] > dens[3]) {curCN[3]=dens[3]} \
#        avg += curCN[5]*(curCN[3]-curCN[2])/(dens[3]-dens[2]); \
#    } \
#    print dens[1], dens[2], dens[3], dens[4], dens[5] / avg * 2}' |

tee $TMPDIR/${sample}.insertions.density.bed |
#Remember wig is 1-indexed
#NB assumes span == step
#bedgraph would be simpler but produces 5x larger file. Would variablestep result in smaller .bw file?
awk -F "\t" 'lastChrom!=$1 || $2 != lastChromEnd || $3-$2 != curspan {curstep=$3-$2; curspan=curstep; print "fixedStep chrom=" $1 " start=" $2+1 " step=" curstep " span=" curspan} {lastChrom=$1; lastChromStart=$2; lastChromEnd=$3; print $5}' > $TMPDIR/${sample}.insertions.density.wig

awk -F "\t" 'BEGIN {OFS="\t"} $5!=0' $TMPDIR/${sample}.insertions.density.bed | starch - > ${OUTDIR}/${sample}.insertions.density.starch

#Kent tools can't use STDIN
wigToBigWig $TMPDIR/${sample}.insertions.density.wig ${chromsizes} ${OUTDIR}/${sample}.insertions.density.bw

echo "track name=${sample}-dens description=\"${sample}-density\" maxHeightPixels=30 color=$trackcolor viewLimits=0:100 autoScale=on visibility=full db=hg38 type=bigWig bigDataUrl=${UCSCbaseURL}/${sample}.insertions.density.bw"


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
    quit(save="no", status=0)
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
cat $OUTDIR/${sample}.uniqcoords.bed | closest-features --dist  --delim '\t' - /vol/isg/annotation/bed/hg38/REsites/Dpn/Dpn.starch | awk -F "\t" 'BEGIN {OFS="\t"} function abs(value) {return (value<0?-value:value);} {if ($6=="-") print $1, $2, $3, $10, abs($13); else if ($6=="+") print $1, $2, $3, $17, $20}' > $OUTDIR/DistDpn.bed

R --quiet --no-save << EOF
#Work around "unable to start device PNG" on ISG cluster
options(bitmapType="cairo") 

filename <- "$OUTDIR/DistDpn.bed"
cat("Reading", filename, "\n")
data <- read(filename)
if(is.null(data)) {
    quit(save="no", status=0)
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
cat $OUTDIR/${sample}.uniqcoords.bed | closest-features --dist  --delim '\t' - /vol/isg/annotation/bed/hg38/REsites/MspI/MspI.starch | awk -F "\t" 'BEGIN {OFS="\t"} function abs(value) {return (value<0?-value:value);} {if ($6=="-") print $1, $2, $3, $10, abs($13); else if ($6=="+") print $1, $2, $3, $17, $20}' > $OUTDIR/DistDpn.bed

R --quiet --no-save << EOF
#Work around "unable to start device PNG" on ISG cluster
options(bitmapType="cairo") 

filename <- "$OUTDIR/DistDpn.bed"
cat("Reading", filename, "\n")
data <- read(filename)
if(is.null(data)) {
    quit(save="no", status=0)
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
    quit(save="no", status=0)
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
unstarch /vol/isg/annotation/bed/hg38/refseq_paint/refGene_hg38.bed6.starch | grep -v promoter > $TMPDIR/refGene_hg38.bed6

bedmap --delim '\t' --echo-map-id $OUTDIR/${sample}.uniqcoords.bed $TMPDIR/refGene_hg38.bed6 | awk -F "\t" 'BEGIN {OFS="\t"} {col=NF} $col~/coding/ {$col="coding"} $col~/TxS-1000/ {$col="TxS-1000"} $col~/TxStart-10kb/ {$col="TxStart-10kb"} $col~/3.UTR/ {$col="3UTR"} $col~/5.UTR/ {$col="5UTR"} $col~/intron/ {$col="intron"} $col=="" || $col~/3.proximal/ {$col="intergenic"} {print;}' | sort -g | uniq -c | sort -k1,1g


echo
echo "Done!!!"
date
