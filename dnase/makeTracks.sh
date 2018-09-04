#!/bin/bash
set -e -o pipefail

alias bedmap='bedmap --ec --sweep-all'

###Parameters
sample=$1
DS=$2
mappedgenome=$3
src=$4

chromsizes="/vol/isg/annotation/fasta/${mappedgenome}/${mappedgenome}.chrom.sizes"


###Setup
if [ -z "$NSLOTS" ]; then
    export NSLOTS=1
fi


# the function "round()" was taken from 
# http://stempell.com/2009/08/rechnen-in-bash/

# the round function:
round()
{
echo $(printf %.$2f $(echo "scale=$2;(((10^$2)*$1)+0.5)/(10^$2)" | bc))
};

bam2instrument()
{
    cut -f1 | cut -d ":" -f1 | 
    #Hack to clean up some encodeproject.org data that has underscore in place of colon after sequencer name
    perl -pe 's/_\d+_\d+_\d+_\d+$//g;' | 
    #SRR data
    perl -pe 's/^(SRR\d+)\.\d+$/$1/g;'
}


#Cheap hack
printfNA()
{
    msg=$1
    shift
    if [ "$1" == "NA" ]; then
        #https://stackoverflow.com/questions/10216340/authoritative-regular-expression-for-sprintf-format-string
        #https://en.wikipedia.org/wiki/Printf_format_string
        #NB don't parse % in type field (we don't want to do substitution on %%)
        msg=`echo -n "$msg" | perl -pe 's/%(\d+\$)?([-+ 0#])*(\d)*(\.?\d+|\.\*)?(ll|[ljqL])?[aAcdeEfFgGinopsuxX]/NA/g;'`
        printf "$msg"
    else
        printf "$msg" "$@"
    fi
}


getcolor () {
    case "$1" in
        B6129SF1J*)
                   trackcolor="238,54,36";; #red
        B6C3F1J*)
                   trackcolor="51,128,195";; #blue
        B6SPRETF1J*)
                   trackcolor="13,85,0";; #lt green
        B6CASTF1J*)
                   trackcolor="0,204,0";; #dk green
        B6PWKF1J)
                   trackcolor="120,88,165";; #purple
       *)
                   trackcolor="0,0,0";; #black
    esac
    echo $trackcolor
}


###Analysis
date

sampleOutdir=$sample
echo "Making tracks for sample $sample ($DS) against genome $mappedgenome"

#TMPDIR=`pwd`/tmp.makeTracks.$sample
#mkdir -p $TMPDIR
echo "using $TMPDIR as TMPDIR"


#BUGBUG define FC
fc="FC"

if grep $DS inputs.txt | grep -q .fastq ; then
    sequencedTags=`cat inputs.txt | grep $DS | sort | uniq | xargs zcat -f | awk 'END {print NR/4}'`
else
    echo "Not merging .fastq files; unable to count total number of sequenced reads"
    sequencedTags="NA"
fi


samflags="-q 20 -F 524"
readlength="36"

#For shorter old Duke data
#readsLength20bp was never defined here
#if [ "$sequencedTags" != "NA" ] && [ `echo "$readsLength20bp/$sequencedTags >= 0.25" | bc -l` == 1 ]; then 
#    echo "More than 25% of reads are 20bp - using q10"
#    samflags="-q 10 -F 524"
#    readlength="20"
#fi


echo "Making bed file"
samtools view $samflags $sampleOutdir/$sample.bam | 
awk -F "\t" 'BEGIN {OFS="\t"} $3!="chrEBV"' |
awk -F "\t" 'BEGIN {OFS="\t"} { \
      readlength = length($10); \
      insertlength = $9; \
#      tagSequence = $10; \
#      color=255; \
#      for(i=12; i<=NF; i++) { \
#         if (match($i, /NM:i:[0-9]+/)) { \
#            editdistance = substr($i, RSTART+5, RLENGTH-5); \
#         } \
#      } \
      if (and($2, 16)) { \
        strand="-"; \
#         colorString = color ",0,0" ;  \
        chromStart=$4-2+readlength; \
      } else { \
        strand="+"; \
#         colorString = "0,0," color; \
        chromStart=$4-1; \
      } \
#      chromStart=$4-1; \
#      chromEnd=chromStart+readlength; \
      chromEnd=chromStart+1; \
#      print $3, chromStart, chromEnd, tagSequence, editdistance, strand, 0, 0, colorString ; \
      print $3, chromStart, chromEnd, insertlength ; \
}' \
| sort-bed --max-mem 5G - | 
#tee $TMPDIR/$sample.tags.bed |
starch - > $sampleOutdir/$sample.tags.starch


echo
echo "Calculating tag counts"
echo "SAMtools statistics"
samtools flagstat $sampleOutdir/$sample.bam > $TMPDIR/$sample.flagstat.txt
cat $TMPDIR/$sample.flagstat.txt

#BUGBUG breaks for DSall or encode reps
#BUGBUG making multiple passes through bam now
PFalignments=`cat $TMPDIR/$sample.flagstat.txt | grep "in total" | awk '{print $1+$3}'`
numMappedReadsMitochondria=`samtools view -c -F 512 $sampleOutdir/$sample.bam chrM`
pctMappedReadsMitochondria=`echo $numMappedReadsMitochondria/$PFalignments*100 | bc -l -q`

#Exclude chrM reads in remaining counts
dupReads=`samtools view -F 512 -f 1024 $sampleOutdir/$sample.bam | awk -F "\t" 'BEGIN {count=0} $3!="chrM" {count+=1} END {print count}'`
analyzedReads=`unstarch $sampleOutdir/$sample.tags.starch | awk -F "\t" '$1!="chrM" {count+=1} END {print count}'`
pctdupReads=`echo "$dupReads/$analyzedReads*100" | bc -l -q`

if [ "$sequencedTags" != "NA" ]; then
    pctPFalignments=`echo "$PFalignments/$sequencedTags*100" | bc -l -q`
    #BUGBUG denominator wrong here
    pctanalyzedReads=`echo "$analyzedReads/$sequencedTags*100" | bc -l -q`
else
    pctPFalignments="NA"
    pctanalyzedReads="NA"
fi

#Tally how many reads were recovered from unpaired/SE reads (NB many of these may not even be PF, so are unrepresented)
PFalignmentsSE=`samtools view -F 1 $sampleOutdir/$sample.bam | wc -l`
analyzedReadsSE=`samtools view -F 513 $sampleOutdir/$sample.bam | awk -F "\t" 'BEGIN {count=0} $3!="chrM" {count+=1} END {print count}'`
if [ "$PFalignmentsSE" == "0" ]; then
    pctanalyzedReadsSE="NA"
else
    pctanalyzedReadsSE=`echo "$analyzedReadsSE/$PFalignmentsSE*100" | bc -l -q`
fi

echo
echo "Making density track"
#BUGBUG double counts fragments where both reads are in window

#Make density track of number of overlapping tags per 150-bp window
#Normalizes the density to 1M tags, ignores enrichment for now
#Note the last awk statement makes the exact intervals conform to Richard's convention that the counts are reported in 20bp windows including tags +/-75 from the center of that window
#Remember wig is 1-indexed (groan)
cat $chromsizes | 
grep -v random | grep -v _hap | grep -v _alt | grep -v chrUn | grep -v chrM | 
awk '{OFS="\t"; $3=$2; $2=0; print}' | sort-bed - | cut -f1,3 | awk 'BEGIN {OFS="\t"} {for(i=0; i<=$2-150; i+=20) {print $1, i, i+150} }' | 
bedmap --bp-ovr 1 --echo --count - $sampleOutdir/$sample.tags.starch | perl -pe 's/\|/\t\t/g;' | awk -F "\t" 'BEGIN {OFS="\t"} {$4="id-" NR; print}' |
awk -F "\t" 'BEGIN {OFS="\t"} {$2+=65; $3-=65; print}' |
awk -v analyzedReads=$analyzedReads -F "\t" 'BEGIN {OFS="\t"} {$5=$5/analyzedReads*1000000; print}' |
tee $TMPDIR/$sample.density.bed |
awk 'lastChr!=$1 {print "fixedStep chrom=" $1 " start=" $2+1 " step=" $3-$2 " span=" $3-$2; lastChr=$1} {print $5}' > $TMPDIR/$sample.wig

starch $TMPDIR/$sample.density.bed > $sampleOutdir/$sample.density.starch

#Kent tools can't use STDIN
wigToBigWig $TMPDIR/$sample.wig $chromsizes $sampleOutdir/$sample.bw

trackcolor=$(getcolor $sample)

analyzedReadsM=`echo "$analyzedReads/1000000" | bc -l -q` 
analyzedReadsM=$(round $analyzedReadsM 1)

#can't find a way to force autoscale=off. http://genome.ucsc.edu/goldenPath/help/bigWig.html implies it's not a parameter in this context
echo "track name=$sample description=\"$sample DNase Density (${analyzedReadsM}M analyzed tags; normalized to 1M)- BWA alignment\" maxHeightPixels=30 color=$trackcolor viewLimits=0:3 autoScale=off visibility=full type=bigWig bigDataUrl=https://cascade.isg.med.nyu.edu/~mauram01/mapped/$fc/$sampleOutdir/$sample.bw"

echo "Making cut count track"
samtools view $samflags $sampleOutdir/$sample.bam | sam2bed --do-not-sort | 
awk -F "\t" 'BEGIN {OFS="\t"} $1!="chrEBV"' |
#BUGBUG should just take this from $sampleOutdir/$sample.tags.starch, but would need to add strand to output, changing format a bit
awk '{if($6=="+"){s=$2; e=$2+1}else{s=$3; e=$3+1} print $1 "\t"s"\t"e"\tid\t1\t"$6 }' | sort-bed - | tee $TMPDIR/$sample.cuts.bed | 
bedops --chop 1 - | awk -F "\t" 'BEGIN {OFS="\t"} {$4="id-" NR; print}' > $TMPDIR/$sample.cuts.loc.bed
bedmap --delim '\t' --echo --count $TMPDIR/$sample.cuts.loc.bed $TMPDIR/$sample.cuts.bed | 
awk -v analyzedReads=$analyzedReads -F "\t" 'BEGIN {OFS="\t"} {$5=$5/analyzedReads*100000000; print}' |
starch - > $sampleOutdir/$sample.perBase.starch

#Skip chrM since UCSC doesn't like the cut count to the right of the last bp in a chromosome
unstarch $sampleOutdir/$sample.perBase.starch | cut -f1-3,5 | awk -F "\t" 'BEGIN {OFS="\t"} $1!="chrM"' > $TMPDIR/$sample.perBase.bedGraph

#Kent tools can't use STDIN
bedGraphToBigWig $TMPDIR/$sample.perBase.bedGraph $chromsizes $sampleOutdir/$sample.perBase.bw

echo "track name=$sample description=\"$sample cut counts (${analyzedReadsM}M nonredundant tags- BWA alignment\" maxHeightPixels=30 color=$trackcolor viewLimits=0:3 autoScale=off visibility=full type=bigWig bigDataUrl=https://cascade.isg.med.nyu.edu/~mauram01/mapped/$fc/$sampleOutdir/$sample.perBase.bw"


#echo "Making fragment coverage track"
#samtools view $samflags $sample.bam | sam2bed --do-not-sort | 
#awk -F "\t" 'BEGIN {OFS="\t"} $11>0 {print $1, $2, $2+$11}' | sort-bed - |
#bedmap --delim '\t' --echo --count $TMPDIR/$sample.cuts.loc.bed - | starch - > $sample.fragCoverage.starch
#right metric?


callHotspots1=0
callHotspots2=1 
if [ "$callHotspots1" == 1 ] || [ "$callHotspots2" == 1 ]; then
    echo
    echo "Will call hotspots"
    date
fi


mappableFile="/vol/isg/annotation/bed/${mappedgenome}/mappability/${mappedgenome}.K36.mappable_only.starch"
#For shorter old Duke data
#if (( $(bc -l <<<"$readsLength20bp/$sequencedTags >=0.25") )); then 
#      mappableFile="/vol/isg/annotation/bed/${mappedgenome}/mappability/${mappedgenome}.K20.mappable_only.starch"
#fi


if [ "$callHotspots1" == 1 ]; then
    echo "Preparing hotspot V1"
    date
    outbase=`pwd`
    mkdir -p $sampleOutdir/hotspots

    #Force creation of new density file (ours is normalized)
    hotspotBAM=$TMPDIR/${sample}.bam
    #250M is too much for hotspot1, 150M takes 12-24 hrs
    if [ $analyzedReads -gt 100000000 ]; then
        echo "$analyzedReads analyzed reads. Generating hotspots on subsample of 100M reads"
        #include -s directly in the variable as recent versions of samtools don't accept -s 1
        sampleAsProportionOfAnalyzedReads="-s "`echo "100000000/$analyzedReads" | bc -l -q`
    else
        echo "$analyzedReads analyzed reads. Generating hotspots on all reads"
    
        sampleAsProportionOfAnalyzedReads=""
    fi


    #BUGBUG hardcoded chromosome names
    samtools view $samflags -b -1 -@ $NSLOTS $sampleAsProportionOfAnalyzedReads $sampleOutdir/$sample.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > $hotspotBAM
    
    
    echo "Calling hotspots"
    #BUGBUG I think hotspot1 can use >40GB memory for some large datasets
    hotspotDens=$outbase/$sampleOutdir/hotspots/$sample.density.starch
    cd $sampleOutdir/hotspots
    $src/callHotspots.sh $hotspotBAM $hotspotDens $outbase/$sampleOutdir/hotspots $mappedgenome $mappableFile > $outbase/$sampleOutdir/hotspots/$sample.log 2>&1
    
    cd ../..
    
    
    hotspotfile=${sample}/hotspots/${sample}-final/${sample}.fdr0.01.hot.bed
    echo "Hotspots for UCSC browser"
    if [ -f "$hotspotfile" ]; then
        cut -f1-3 $hotspotfile > $TMPDIR/${sample}.fdr0.01.hot.bed
        bedToBigBed -type=bed3 $TMPDIR/${sample}.fdr0.01.hot.bed $chromsizes ${sample}/hotspots/${sample}.fdr0.01.hot.bb
    else
        echo "Can't find $hotspotfile to make bigBed"
    fi
    
    peakfile=${sample}/hotspots/${sample}-final/${sample}.fdr0.01.pks.bed
    if [ -f "$peakfile" ]; then
        cut -f1-3 $peakfile > $TMPDIR/${sample}.fdr0.01.pks.bed
        bedToBigBed -type=bed3 $TMPDIR/${sample}.fdr0.01.pks.bed $chromsizes ${sample}/hotspots/${sample}.fdr0.01.pks.bb
    else
        echo "Can't find $peakfile to make bigBed"
    fi
    
    spotout=${sample}/hotspots/$sample.spot.out
    
    
    #subsample for spot 
    if [ $analyzedReads -gt 10000000 ]; then
        echo
        echo "$analyzedReads analyzed reads. Calculating SPOT score on subsample of 10M reads"
        sampleAsProportionOfAnalyzedReads=`echo "10000000/$analyzedReads" | bc -l -q`
        samtools view $samflags -b -1 -@ $NSLOTS -s $sampleAsProportionOfAnalyzedReads ${sample}/$sample.bam > $TMPDIR/${sample}.10Mtags.bam
        
        mkdir -p $TMPDIR/$sample.hotspots.10Mtags
        cd $TMPDIR/$sample.hotspots.10Mtags
        
        #NB dens track doesn't exist
        $src/callHotspots.sh $TMPDIR/${sample}.10Mtags.bam $TMPDIR/${sample}.10Mtags.density.starch $TMPDIR/$sample.hotspots.10Mtags $mappedgenome $mappableFile > $outbase/$sampleOutdir/hotspots/$sample.10Mtags.log 2>&1
        
        spotout=$TMPDIR/${sample}.hotspots.10Mtags/${sample}.10Mtags.spot.out
        
        #NB otherwise prints pwd
        cd - > /dev/null
    fi
    
    echo "Done calling hotspots"
    date
    echo
fi

if [ "$callHotspots2" == 1 ]; then
    FDRhot2="0.05"
    if (( $(bc -l <<<"$FDRhot2 >=0.05") )); then
        #COMMENT Hotspot2 calls all hotspots for one FDR threshold. Further filtering is done through the for loop below
        echo "Calling hotspot2"
        date
        mkdir -p ${sampleOutdir}/hotspot2
        
        
        hotspot2.sh -c /vol/isg/annotation/bed/${mappedgenome}/hotspots2/${mappedgenome}.chrom.sizes -C /vol/isg/annotation/bed/${mappedgenome}/hotspots2/${mappedgenome}.CenterSites.starch -F $FDRhot2 -f $FDRhot2 -M ${mappableFile} $sampleOutdir/$sample.bam ${sampleOutdir}/hotspot2 > $sampleOutdir/hotspot2/$sample.log 2>&1
        for FDR in {0.05,0.01,0.005,0.001}; do
            echo "Thresholding hotspots at FDR $FDR"
            hsmerge.sh -f ${FDR} ${sampleOutdir}/hotspot2/${sample}.allcalls.starch ${sampleOutdir}/hotspot2/${sample}.hotspots.fdr${FDR}.starch
        done
        
        
        echo "Hotspots for UCSC browser"
        #Hotspots
        hotspot2file=${sampleOutdir}/hotspot2/${sample}.hotspots.fdr0.01.starch
        hotspot2fileFDR05=${sampleOutdir}/hotspot2/${sample}.hotspots.fdr0.05.starch
        if [ -f "$hotspot2fileFDR05" ] && [ `unstarch --elements $hotspot2fileFDR05` -gt 0 ]; then
            unstarch $hotspot2fileFDR05 | cut -f1-4 > $TMPDIR/${sample}.hotspots.fdr0.05.bed
            bedToBigBed -type=bed4 $TMPDIR/${sample}.hotspots.fdr0.05.bed $chromsizes ${sampleOutdir}/hotspot2/${sample}.hotspots.fdr0.05.bb
        fi
        
        #Peaks
        hotspot2peakfile=${sampleOutdir}/hotspot2/${sample}.peaks.starch
        if [ -f "$hotspot2peakfile" ] && [ `unstarch --elements $hotspot2peakfile` -gt 0 ]; then
            unstarch $hotspot2peakfile | cut -f1-4 > $TMPDIR/${sample}.peaks.bed
            bedToBigBed -type=bed4 $TMPDIR/${sample}.peaks.bed $chromsizes ${sampleOutdir}/hotspot2/${sample}.peaks.bb
        fi
        
        echo "Done calling hotspot2"
        date
        echo
    else 
        echo "WARNING: Hotspot2 should be run with FDR >= 0.05"
        echo "Skipping hotspot2"
        echo
    fi
fi


#Stats
echo
echo "*** Overall Stats ***"
echo
echo -e "Num_sequenced_reads\t$sequencedTags\t\t$sample"
printfNA "Num_pass_filter_alignments\t$PFalignments\t%.1f%%\t$sample\n" "$pctPFalignments"
printfNA "Num_mitochondria_reads\t$numMappedReadsMitochondria\t%.1f%%\t$sample\n" "$pctMappedReadsMitochondria"
printfNA "Num_analyzed_reads\t$analyzedReads\t%.1f%%\t$sample\n" "$pctanalyzedReads"
printfNA "Num_duplicate_reads\t$dupReads\t%.1f%%\t$sample\n" "$pctdupReads"

#Don't have denominator of unpaired reads we tried to map, so don't compute % for first
echo -e "Num_SE_pass_filter_alignments\t$PFalignmentsSE\t\t$sample"
#NB denominator is PF reads, not pctMappedReadsMitochondrias sequenced
printfNA "Num_SE_analyzed_reads\t$analyzedReadsSE\t%.1f%%\t$sample\n" "$pctanalyzedReadsSE"


if [ -f "$hotspotfile" ]; then
    numhotspots=`cat $hotspotfile | wc -l`
else
    numhotspots="NA"
fi
echo -e "Num_hotspots\t$numhotspots\t\t$sample"


if [ -f "$spotout" ]; then
    spot=`cat $spotout | awk 'BEGIN {OFS="\t"} NR==2 {print $3}'`
else
    spot="NA"
fi
echo -e "SPOT\t$spot\t\t$sample"


if [ -f "$hotspot2file" ]; then
    numhotspots2=`unstarch --elements $hotspot2file`
else
    numhotspots2="NA"
fi
echo -e "Num_hotspots2\t$numhotspots2\t\t$sample"


if [ -f "$hotspot2fileFDR05" ] && [ "$numhotspots2" -gt 0 ]; then
    numTotalReadsHotspot2=$(unstarch ${sampleOutdir}/hotspot2/${sample}.cutcounts.starch | awk '{print $5}'|  awk '{sum+=$1} END {print sum}')
    numReadsinDHSHotspot2=$(unstarch ${sampleOutdir}/hotspot2/${sample}.cutcounts.starch | bedops --header -e -1 - $hotspot2fileFDR05 | awk '{print $5}' | awk '{sum+=$1} END {print sum}')
    spot2=$(echo "scale=4; $numReadsinDHSHotspot2/$numTotalReadsHotspot2" | bc)
else
    spot2="NA"
fi
echo -e "SPOT2\t$spot2\t\t$sample"


echo
echo "Histogram of mismatches"
samtools view $samflags $sampleOutdir/$sample.bam | awk -F "\t" 'BEGIN {OFS="\t"} { \
    for(i=12; i<=NF; i++) { \
         if(match($i, /NM:i:/)) { \
            print substr($i, 6); \
         } \
    } \
}' | sort -g | uniq -c | sort -k2,2 -g


#NB XA tags are computed for the unpaired tags, while MAPQ reflects the final PE location
echo
echo "Histogram of number of best alignments"
samtools view  $samflags $sampleOutdir/$sample.bam | awk -F "\t" 'BEGIN {OFS="\t"} { \
    tag="NA"; \
    for(i=12; i<=NF; i++) { \
         if(match($i, /X0:i:/)) { \
            tag=substr($i, 6); \
#            if(tag>1) {print} \
         } \
#         if(match($i, /XA:Z:/)) { \
#            tag=length(split(substr($i, 6), xa, ";")); \
#         } \
    } \
    print tag; \
}' | awk -F "\t" 'BEGIN {OFS="\t"} $0>=10 {$0="10+"} {print}' | sort -g | uniq -c | sort -k2,2 -g


PEtags=`cat $TMPDIR/$sample.flagstat.txt | grep "paired in sequencing" | awk '{print $1+$3}'`
if [ $PEtags -gt 0 ]; then
    echo
    echo "Template lengths"
    samtools view $samflags $sampleOutdir/$sample.bam | awk -F "\t" 'BEGIN {OFS="\t"} $3!="chrM"' | awk -v sample=$sample -F "\t" 'BEGIN {OFS="\t"} $9>0 {print sample, $9}' | sort -k2,2n | tee $sampleOutdir/$sample.insertlengths.txt | cut -f2 |
    awk -F "\t" 'BEGIN {OFS="\t"} NR==1 {print "Minimum: " $0} {cum+=$0; lengths[NR]=$0; if($0<125) {lastLineUnder125=NR}} END {print "Number of tags: " NR; print "25th percentile: " lengths[int(NR*0.25)]; print "50th percentile: " lengths[int(NR*0.5)]; print "75th percentile: " lengths[int(NR*0.75)]; print "95th percentile: " lengths[int(NR*0.95)]; print "99th percentile: " lengths[int(NR*0.99)]; print "Maximum: " $0; print "Mean: " cum/NR; print "Prop. tags under 125 bp: " lastLineUnder125/NR}'
    gzip -f $sampleOutdir/$sample.insertlengths.txt
fi

echo
echo "Tag lengths:"
samtools view $samflags $sampleOutdir/$sample.bam | awk -F "\t" 'BEGIN {OFS="\t"} {print length($10)}' | sort -g | uniq -c | sort -k1,1g 

echo
echo "Tag count by sequencing instrument"
samtools view $samflags $sampleOutdir/$sample.bam | bam2instrument | sort -g | uniq -c | sort -k1,1g


echo
echo -e "\nDone!"
date
