#!/bin/bash
set -e # -o pipefail


###Parameters
sample=$1
DS=$2
mappedgenome=$3
src=$4

chromsizes="/vol/isg/annotation/fasta/${mappedgenome}/${mappedgenome}.chrom.sizes"

sampleOutdir=$sample
echo "Making tracks for sample $sample ($DS) against genome $mappedgenome"

#TMPDIR=`pwd`/tmp.makeTracks.$sample
#mkdir -p $TMPDIR
echo "using $TMPDIR as TMPDIR"

samflags="-q 20 -F 524"
#NB chrM being considered in most downstream analyses


###Setup
if [ -z "$NSLOTS" ]; then
       export NSLOTS=1
fi

alias bedmap='bedmap --ec --sweep-all'

# the function "round()" was taken from 
# http://stempell.com/2009/08/rechnen-in-bash/

# the round function:
round()
{
echo $(printf %.$2f $(echo "scale=$2;(((10^$2)*$1)+0.5)/(10^$2)" | bc))
};

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
echo "Making bed file"
samtools view $samflags $sampleOutdir/$sample.bam | awk -F "\t" 'BEGIN {OFS="\t"} { \
      readlength = length($10); \
      insertlength = $9; \
#      tagSequence = $10; \
#      color=255; \
#      for(i=12; i<=NF; i++) { \
#            if (match($i, /NM:i:[0-9]+/)) { \
#                  editdistance = substr($i, RSTART+5, RLENGTH-5); \
#            } \
#      } \
      if (and($2, 16)) { \
              strand="-"; \
#            colorString = color ",0,0" ;  \
              chromStart=$4-2+readlength; \
      } else { \
              strand="+"; \
#            colorString = "0,0," color; \
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
sequencedTags=`cat inputs.txt | grep $DS | sort | uniq | xargs zcat | awk 'END {print NR/4}'`
PFalignments=`cat $TMPDIR/$sample.flagstat.txt | grep "in total" | awk '{print $1+$3}'`
#NB used to count both columns of flagstat output ($1 + $3) for remaining metrics but it gives no guarantee of 1 line per tag
uniqMappedTags=`cat $TMPDIR/$sample.flagstat.txt | grep "mapped (" | awk '{print $1}'`
numMappedTagsMitochondria=`samtools view -c -F 512 $sampleOutdir/$sample.bam chrM`
propMappedTagsMitochondria=`echo $numMappedTagsMitochondria/$uniqMappedTags*100 | bc -l -q`
#NB now excludes chrM reads in duplicate counts
dupTags=`samtools view -c -F 512 -f 1024 $sampleOutdir/$sample.bam | awk -F "\t" '$3!="chrM" {count+=1} END {print count}'`
#dupTags=`cat $TMPDIR/$sample.flagstat.txt | grep "duplicates" | awk '{print $1}'`
analyzedTags=`unstarch --elements $sampleOutdir/$sample.tags.starch`
pctPFalignments=`echo $PFalignments/$sequencedTags*100 | bc -l -q`
pctUniqMappedTags=`echo $uniqMappedTags/$sequencedTags*100 | bc -l -q`
pctDupTags=`echo $dupTags/$uniqMappedTags*100 | bc -l -q`
pctAnalyzedTags=`echo $analyzedTags/$sequencedTags*100 | bc -l -q`

#Tally how many reads were recovered from unpaired/SE reads (NB many of these may not even be PF, so are unrepresented)
PFalignmentsSE=`samtools view -F 1 $sampleOutdir/$sample.bam | wc -l`
uniqMappedTagsSE=`samtools view -q 20 -F 525 $sampleOutdir/$sample.bam | wc -l`
pctUniqMappedTagsSE=`echo $uniqMappedTagsSE/$PFalignmentsSE*100 | bc -l -q`


echo
echo "Making density track"
#Make density track of number of overlapping tags per 150-bp window
#Normalizes the density to 1M tags, ignores enrichment for now
#Note the last awk statement makes the exact intervals conform to Richard's convention that the counts are reported in 20bp windows including tags +/-75 from the center of that window
#Remember wig is 1-indexed (groan)
cat $chromsizes | 
grep -v random | grep -v _hap | grep -v chrM | grep -v _alt|grep -v Un|
awk '{OFS="\t"; $3=$2; $2=0; print}' | sort-bed - | cut -f1,3 | awk 'BEGIN {OFS="\t"} {for(i=0; i<=$2-150; i+=20) {print $1, i, i+150} }' | 
bedmap --bp-ovr 1 --echo --count - $sampleOutdir/$sample.tags.starch | perl -pe 's/\|/\t\t/g;' | awk -F "\t" 'BEGIN {OFS="\t"} {$4="id-" NR; print}' |
awk -F "\t" 'BEGIN {OFS="\t"} {$2+=65; $3-=65; print}' |
awk -v analyzedTags=$analyzedTags -F "\t" 'BEGIN {OFS="\t"} {$5=$5/analyzedTags*1000000; print}' |
tee $TMPDIR/$sample.density.bed |
awk 'lastChr!=$1 {print "fixedStep chrom=" $1 " start=" $2+1 " step=" $3-$2 " span=" $3-$2; lastChr=$1} {print $5}' > $TMPDIR/$sample.wig

starch $TMPDIR/$sample.density.bed > $sampleOutdir/$sample.density.starch

wigToBigWig $TMPDIR/$sample.wig $chromsizes $sampleOutdir/$sample.bw

trackcolor=$(getcolor $sample)

analyzedTagsM=`echo $analyzedTags/1000000 | bc -l -q` 
analyzedTagsM=$(round $analyzedTagsM 1)

#can't find a way to force autoscale=off. http://genome.ucsc.edu/goldenPath/help/bigWig.html implies it's not a parameter in this context
echo "track name=$sample description=\"$sample DNase Density (${analyzedTagsM}M analyzed tags; normalized to 1M)- BWA alignment\" maxHeightPixels=30 color=$trackcolor viewLimits=0:3 autoScale=off visibility=full type=bigWig bigDataUrl=https://cascade.isg.med.nyu.edu/mauranolab/encode/mapped/$sampleOutdir/$sample.bw"

echo "Making cut count track"
samtools view $samflags $sampleOutdir/$sample.bam | sam2bed --do-not-sort | 
awk '{if($6=="+"){s=$2; e=$2+1}else{s=$3; e=$3+1} print $1 "\t"s"\t"e"\tid\t1\t"$6 }' | sort-bed - | tee $TMPDIR/$sample.cuts.bed | 
bedops --chop 1 - | awk -F "\t" 'BEGIN {OFS="\t"} {$4="id-" NR; print}' > $TMPDIR/$sample.cuts.loc.bed
bedmap --delim '\t' --echo --count $TMPDIR/$sample.cuts.loc.bed $TMPDIR/$sample.cuts.bed | 
awk -v analyzedTags=$analyzedTags -F "\t" 'BEGIN {OFS="\t"} {$5=$5/analyzedTags*100000000; print}' |
starch - > $sampleOutdir/$sample.perBase.starch

#Skip chrM since UCSC doesn't like the cut count to the right of the last bp in a chromosome
gcat $sampleOutdir/$sample.perBase.starch | cut -f1-3,5 | awk -F "\t" 'BEGIN {OFS="\t"} $1!="chrM"' > $TMPDIR/$sample.perBase.bedGraph

bedGraphToBigWig $TMPDIR/$sample.perBase.bedGraph $chromsizes $sampleOutdir/$sample.perBase.bw

echo "track name=$sample description=\"$sample cut counts (${analyzedTagsM}M analyzed tags- BWA alignment\" maxHeightPixels=30 color=$trackcolor viewLimits=0:3 autoScale=off visibility=full type=bigWig bigDataUrl=https://cascade.isg.med.nyu.edu/mauranolab/encode/mapped/$sampleOutdir/$sample.perBase.bw"


#echo "Making fragment coverage track"
#samtools view $samflags $sample.bam | sam2bed --do-not-sort | 
#awk -F "\t" 'BEGIN {OFS="\t"} $11>0 {print $1, $2, $2+$11}' | sort-bed - |
#bedmap --delim '\t' --echo --count $TMPDIR/$sample.cuts.loc.bed - | starch - > $sample.fragCoverage.starch
#right metric?


echo
echo "Calling hotspots"
outbase=`pwd`
mkdir -p $sampleOutdir/hotspots


#Force creation of new density file (ours is normalized)
hotspotDens=$outbase/$sampleOutdir/hotspots/$sample.density.starch
hotspotBAM=$TMPDIR/${sample}.bam
#250M is too much, 150M takes 12-24 hrs
if [ $uniqMappedTags -gt 100000000 ]; then
       echo "$uniqMappedTags uniquely mapped tags. Generating hotspots on subsample of 100M reads"
       
       sampleAsProportionOfUniqMappedTags=`echo 100000000/$uniqMappedTags | bc -l -q`
else
       echo "$uniqMappedTags uniquely mapped tags. Generating hotspots on all reads"
       
       sampleAsProportionOfUniqMappedTags=1
fi

samtools view $samflags -b -1 -@ $NSLOTS -s $sampleAsProportionOfUniqMappedTags $sampleOutdir/$sample.bam > $hotspotBAM

cd $sampleOutdir/hotspots

#BUGBUG I think hotspot1 can use >40GB memory for some large datasets
$src/callHotspots.sh $hotspotBAM $hotspotDens $outbase/$sampleOutdir/hotspots > $outbase/$sampleOutdir/hotspots/$sample.log 2>&1

cd ../..


hotspotfile=hotspots/${sample}/${sample}-final/${sample}.fdr0.01.hot.bed
echo "Hotspots for UCSC browser"
if [ -f "$hotspotfile" ]; then
       cut -f1-3 $hotspotfile > $TMPDIR/${sample}.fdr0.01.hot.bed
       bedToBigBed -type=bed3 $TMPDIR/${sample}.fdr0.01.hot.bed $chromsizes hotspots/${sample}/${sample}.fdr0.01.hot.bb
else
       echo "Can't find $hotspotfile to make bigBed"
fi

peakfile=hotspots/${sample}/${sample}-final/${sample}.fdr0.01.pks.bed
if [ -f "$peakfile" ]; then
       cut -f1-3 $peakfile > $TMPDIR/${sample}.fdr0.01.pks.bed
       bedToBigBed -type=bed3 $TMPDIR/${sample}.fdr0.01.pks.bed $chromsizes hotspots/${sample}/${sample}.fdr0.01.pks.bb
else
       echo "Can't find $peakfile to make bigBed"
fi


spotout=hotspots/${sample}/$sample.spot.out

#subsample for spot 
if [ $uniqMappedTags -gt 10000000 ]; then
       echo
       echo "$uniqMappedTags uniquely mapped tags. Calculating SPOT score on subsample of 10M reads"
       sampleAsProportionOfUniqMappedTags=`echo 10000000/$uniqMappedTags | bc -l -q`
       samtools view $samflags -b -1 -@ $NSLOTS -s $sampleAsProportionOfUniqMappedTags $sample.bam > $TMPDIR/${sample}.10Mtags.bam
       
       mkdir -p $TMPDIR/$sample.hotspots.10Mtags
       cd $TMPDIR/$sample.hotspots.10Mtags
       
       #NB dens track doesn't exist
       $src/callHotspots.sh $TMPDIR/${sample}.10Mtags.bam $TMPDIR/${sample}.10Mtags.density.starch $TMPDIR/$sample.hotspots.10Mtags > $outbase/$sampleOutdir/hotspots/$sample.10Mtags.log 2>&1
       
       #NB otherwise prints pwd
       cd - > /dev/null
       
       spotout=$TMPDIR/${sample}.hotspots.10Mtags/${sample}.10Mtags.spot.out
fi


#Stats
echo
echo "*** Overall Stats ***"
echo
printf "Num_sequenced_tags\t$sequencedTags\t\t%s\n" "$sample"
printf "Num_pass_filter_alignments\t$PFalignments\t%.1f%%\t%s\n" "$pctPFalignments" "$sample"
printf "Num_uniquely_mapped_tags\t$uniqMappedTags\t%.1f%%\t%s\n" "$pctUniqMappedTags" "$sample"
printf "Num_mitochondria_tags\t$numMappedTagsMitochondria\t%.1f%%\t%s\n" "$propMappedTagsMitochondria" "$sample"
printf "Num_duplicate_tags\t$dupTags\t%.1f%%\t%s\n" "$pctDupTags" "$sample"
printf "Num_analyzed_tags\t$analyzedTags\t%.1f%%\t%s\n" "$pctAnalyzedTags" "$sample"

#Don't have denominator of unpaired tags we tried to map, so don't compute % for first
printf "Num_SE_pass_filter_alignments\t$PFalignmentsSE\t\t%s\n" "$sample"
#NB denominator is PF tags, not tags sequenced
printf "Num_SE_uniquely_mapped_tags\t$uniqMappedTagsSE\t%.1f%%\t%s\n"  "$pctUniqMappedTagsSE" "$sample"


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
#                  if(tag>1) {print} \
            } \
#            if(match($i, /XA:Z:/)) { \
#                  tag=length(split(substr($i, 6), xa, ";")); \
#            } \
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


#Hack to deal with read names from SRA
if samtools view $sampleOutdir/$sample.bam | cut -f1 | head -10 | grep -v -q -e "^SRR"; then
       echo
       echo "Tag count by sequencing instrument"
       samtools view $samflags $sampleOutdir/$sample.bam | cut -f1 | cut -d ":" -f1 | sort | uniq -c | sort -k1,1g
fi

echo
echo -e "\nDone!"
date
