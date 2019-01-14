#!/bin/bash
set -eu -o pipefail

shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'

###Parameters
mappedgenome=$1
analysisType=$2
name=$3
BS=$4
src=$5


#processingCommand=`echo "${analysisType}" | awk -F "," '{print $1}'`
analysisCommand=`echo "${analysisType}" | awk -F "," '{print $2}'`


if [[ "${analysisCommand}" != "none" ]] && [[ "${analysisCommand}" != "atac" ]] && [[ "${analysisCommand}" != "dnase" ]] && [[ "${analysisCommand}" != "chipseq" ]] && [[ "${analysisCommand}" != "callsnps" ]]; then 
    echo "ERROR: unknown analysis command ${analysisCommand} in analysisType ${analysisType}"
    exit 1
fi

#TODO label capture vs. WGS
case "${analysisCommand}" in
    callsnps)
        ucscTrackDescriptionDataType="";;
    dnase)
        ucscTrackDescriptionDataType="DNase-seq";;
    atac)
        ucscTrackDescriptionDataType="ATAC-seq";;
    chipseq)
        ucscTrackDescriptionDataType="ChIP-seq";;
    *)
        ucscTrackDescriptionDataType="${analysisCommand}";;
esac


source ${src}/genomeinfo.sh ${mappedgenome}

#Deal with some of the more complex reference index names
#NB this will call hotspots only on the first mammalian genome for the *_sacCer3 hybrid indices
annotationgenome=`echo ${mappedgenome} | perl -pe 's/_.+$//g;' -e 's/all$//g;'`


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
    samtools view -F 2304 $* |
    cut -f1 | cut -d ":" -f1 |
    #Hack to clean up some encodeproject.org data that has a hyphen separating FC from read name
    cut -d "-" -f1 | 
    #Hack to clean up some encodeproject.org data that has underscore in place of colon after sequencer name
    perl -pe 's/_\d+_\d+_\d+_\d+$//g;' | 
    #SRR data
    perl -pe 's/^(SRR\d+)(\.\d+)+$/$1/g;'
}


#Cheap hack
printfNA()
{
    local msg=$1
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
    local trackcolor
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
#TMPDIR=`pwd`/tmp.analysis.${name}
#mkdir -p $TMPDIR
echo "using $TMPDIR as TMPDIR"

sampleOutdir=${name}
echo "Running ${analysisType} analysis for sample ${name} (${BS}) against genome ${mappedgenome}"


#Required files
if [ ! -s "${sampleOutdir}/${name}.${mappedgenome}.bam" ]; then
    echo "ERROR: can not find file ${sampleOutdir}/${name}.${mappedgenome}.bam"
    exit 2
fi

if [ ! -s "${sampleOutdir}/${name}.${mappedgenome}.bam.bai" ]; then
    echo "ERROR: can not find file ${sampleOutdir}/${name}.${mappedgenome}.bam.bai"
    exit 2
fi


if grep ${BS} inputs.txt | grep -q .fastq; then
    sequencedReads=`cat inputs.txt | grep ${BS} | grep .fastq | sort | uniq | xargs pigz -dc -f | awk 'END {print NR/4}'`
elif grep ${BS} inputs.txt | grep ${mappedgenome} | grep -q .bam; then
    #Look for inputs.txt one directory up from the bam files and get fastq files from there
    #Re-counting fastq files really too slow for large jobs so try to grab the counts from the individual analysis logfiles
    sequencedReads=`cat inputs.txt | grep ${BS} | grep .bam | grep ${mappedgenome} | xargs -I {} dirname {} | xargs -I {} find {} -name "analysis*.${mappedgenome}.o*" | xargs awk -F "\t" 'BEGIN {OFS="\t"; sum=0} $1== "Num_sequenced_reads" {sum+=$2} END {print sum}'`
else
    echo "Couldn't identify source of .fastq files; unable to count total number of sequenced reads"
    sequencedReads="NA"
fi


#Prep for UCSC track links
projectdir=`pwd | perl -pe 's/^\/vol\/(cegs|mauranolab|isg\/encode)\///g;'`
if [[ `pwd` =~ ^\/vol\/cegs\/ ]]; then
    UCSCbaseURL="https://cegs@cascade.isg.med.nyu.edu/cegs/${projectdir}/${sampleOutdir}"
elif [[ `pwd` =~ ^\/vol\/isg\/encode\/ ]]; then
    UCSCbaseURL="https://cascade.isg.med.nyu.edu/mauranolab/encode/${projectdir}/${sampleOutdir}"
else
    UCSCbaseURL="https://mauranolab@cascade.isg.med.nyu.edu/~mauram01/${projectdir}/${sampleOutdir}"
fi

trackcolor=$(getcolor ${name})


#NB using -F 512 delegates the thresholding on MAPQ/unpapped reads to filter_reads.py
samflags="-F 512"


#For shorter old Duke data
#readsLength20bp was never defined here
#if [ "${sequencedReads}" != "NA" ] && [ `echo "$readsLength20bp/${sequencedReads} >= 0.25" | bc -l` == 1 ]; then 
#    echo "More than 25% of reads are 20bp - using q10"
#    samflags="-q 10 -F 524"
#fi


echo
echo "Making bed file"
date
#Coordinates are the 5' end of the read
#Do per-chromosome to reduce size of sorts
for chrom in `samtools idxstats ${sampleOutdir}/${name}.${mappedgenome}.bam | cut -f1 | awk -F "\t" 'BEGIN {OFS="\t"} $1!="*" {print $1, 0, 1}' | sort-bed - | cut -f1`; do
    samtools view ${samflags} ${sampleOutdir}/${name}.${mappedgenome}.bam ${chrom} | 
    awk -F "\t" 'BEGIN {OFS="\t"} $3!="chrEBV"' |
    awk -F "\t" 'BEGIN {OFS="\t"} { \
        flag = $2; \
        softclipping = 0; \
        if(match($6, /^[0-9]+S/)) {softclipping+=substr($6, 1, RLENGTH-1)} \
        if(match($6, /[0-9]+S$/)) {softclipping+=substr($6, RSTART)} \
        readlength = length($10) - softclipping; \
        insertlength = $9; \
    #    readSequence = $10; \
    #    color=255; \
    #    for(i=12; i<=NF; i++) { \
    #        if (match($i, /NM:i:[0-9]+/)) { \
    #            editdistance = substr($i, RSTART+5, RLENGTH-5); \
    #        } \
    #    } \
        chromStart=$4-1; \
        if (and($2, 16)) { \
            strand="-"; \
            #colorString = color ",0,0" ;  \
            chromStart=chromStart+readlength-1; \
        } else { \
            strand="+"; \
            #colorString = "0,0," color; \
        } \
        #chromEnd=chromStart+readlength; \
        chromEnd=chromStart+1; \
        #print $3, chromStart, chromEnd, readSequence, editdistance, strand, 0, 0, colorString ; \
        print $3, chromStart, chromEnd, insertlength, readlength, strand, flag; \
    }' |
    sort-bed --max-mem 5G - | 
    #BUGBUG bedmap not working properly with header
    #awk -F "\t" 'BEGIN {OFS="\t"; print "#chrom", "chromStart", "chromEnd", "insertlength", "readlength", "strand"} {print}' |
    starch - > $TMPDIR/${name}.${mappedgenome}.reads.${chrom}.starch
done
starchcat $TMPDIR/${name}.${mappedgenome}.reads.*.starch > ${sampleOutdir}/${name}.${mappedgenome}.reads.starch


echo
echo "Calculating read counts"
date
echo "SAMtools statistics"
samtools flagstat ${sampleOutdir}/${name}.${mappedgenome}.bam > $TMPDIR/${name}.${mappedgenome}.flagstat.txt
cat $TMPDIR/${name}.${mappedgenome}.flagstat.txt

#BUGBUG breaks for DSall or encode reps
PFalignments=`cat $TMPDIR/${name}.${mappedgenome}.flagstat.txt | grep "in total" | awk '{print $1+$3}'`
numMappedReadsMitochondria=`samtools view -c -F 512 ${sampleOutdir}/${name}.${mappedgenome}.bam chrM`
pctMappedReadsMitochondria=`echo ${numMappedReadsMitochondria}/${PFalignments}*100 | bc -l -q`

#Exclude chrM reads in remaining counts
#BUGBUG making 4 passes through full bam below
#samtools view -c is faster than counting wc -l/awk; have a cleanup awk sum the lines in case xargs breaks it into multiple commands
dupReads=`samtools idxstats ${sampleOutdir}/${name}.${mappedgenome}.bam | awk -F "\t" 'BEGIN {OFS="\t"} $1!="chrM" && $1!="*" {print $1}' | xargs samtools view -F 512 -f 1024 -c ${sampleOutdir}/${name}.${mappedgenome}.bam | awk 'BEGIN {sum=0} {sum+=$0} END {print sum}'`
analyzedReads=`unstarch --list ${sampleOutdir}/${name}.${mappedgenome}.reads.starch | perl -pe 's/ +\|/\t/g;' | awk -F "\t" 'BEGIN {OFS="\t"; sum=0} NR==1 {for(i=1; i<=NF; i++) {if($i=="chr") {chromCol=i} else if($i=="uncompressedLineCount") {uncompressedLineCountCol=i} } } NR>1 && $chromCol!="chrM" {sum+=$uncompressedLineCountCol} END {print sum}'`


if [[ "${analyzedReads}" == 0 ]]; then
    pctdupReads="NA"
else
    pctdupReads=`echo "${dupReads}/${analyzedReads}*100" | bc -l -q`
fi

analyzedReadsM=`echo "${analyzedReads}/1000000" | bc -l -q` 
analyzedReadsM=$(round ${analyzedReadsM} 1)

if [ "${sequencedReads}" != "NA" ]; then
    pctPFalignments=`echo "${PFalignments}/${sequencedReads}*100" | bc -l -q`
    #BUGBUG denominator wrong here
    pctanalyzedReads=`echo "${analyzedReads}/${sequencedReads}*100" | bc -l -q`
else
    pctPFalignments="NA"
    pctanalyzedReads="NA"
fi

#Tally how many reads were recovered from unpaired/SE reads (NB many of these may not even be PF, so are unrepresented)
PFalignmentsSE=`samtools view -F 1 -c ${sampleOutdir}/${name}.${mappedgenome}.bam`
analyzedReadsSE=`samtools idxstats ${sampleOutdir}/${name}.${mappedgenome}.bam | awk -F "\t" 'BEGIN {OFS="\t"} $1!="chrM" && $1!="*" {print $1}' | xargs samtools view -F 513 -c ${sampleOutdir}/${name}.${mappedgenome}.bam | awk 'BEGIN {sum=0} {sum+=$0} END {print sum}'`
if [ "${PFalignmentsSE}" == "0" ]; then
    pctanalyzedReadsSE="NA"
else
    pctanalyzedReadsSE=`echo "${analyzedReadsSE}/${PFalignmentsSE}*100" | bc -l -q`
fi


#Now that we have analyzedReadsM we can print this track line, which is universal for all analysisCommand
echo
echo "Making BAM track"
echo "track type=bam name=${name}-reads description=\"${name} reads (${analyzedReadsM}M nonredundant reads)\" visibility=pack bigDataUrl=${UCSCbaseURL}/${name}.${mappedgenome}.bam pairEndsByName=T visibility=dense maxWindowToDraw=10000"


if [[ "${analysisCommand}" == "callsnps" ]]; then
    if [ "${PFalignments}" -lt 50000000 ] && [[ "${analyzedReads}" > 0 ]]; then
        #Call hotspot 2 for debugging weird issues in coverage tracks, but not on large datasets
        callHotspots1=0
        callHotspots2=1
    else
        echo "Will skip hotspot -- too many reads"
        callHotspots1=0
        callHotspots2=0
    fi
    if [[ "${analyzedReads}" == 0 ]]; then
        genomecov=0
    else
        echo
        echo "Calling SNPs"
        date
        
        #TODO switch to breaking up by coordinate rather than chromosome to split up larger jobs
        samtools idxstats ${sampleOutdir}/${name}.${mappedgenome}.bam | awk -F "\t" 'BEGIN {OFS="\t"} $1!="*" {print $1}' > ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt
        #unstarch --list-chromosomes ${sampleOutdir}/${name}.${mappedgenome}.reads.starch > ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt
        n=`cat ${sampleOutdir}/inputs.callsnps.${mappedgenome}.txt | wc -l`
        qsub -S /bin/bash -cwd -V -terse -j y -b y -t 1-${n} -o ${sampleOutdir} -N callsnps.${name}.${mappedgenome} "${src}/callsnpsByChrom.sh ${mappedgenome} ${analysisType} ${name} ${src}" | perl -pe 's/[^\d].+$//g;' > ${sampleOutdir}/sgeid.callsnps.${mappedgenome}
        qsub -S /bin/bash -cwd -V -terse -j y -b y -hold_jid `cat ${sampleOutdir}/sgeid.callsnps.${mappedgenome}` -o ${sampleOutdir} -N callsnpsMerge.${name}.${mappedgenome} "${src}/callsnpsMerge.sh ${mappedgenome} ${analysisType} ${name} ${src}" | perl -pe 's/[^\d].+$//g;'
        rm -f ${sampleOutdir}/sgeid.callsnps.${mappedgenome}
        
        
        echo
        echo "Collecting picard metrics"
        mkdir -p ${sampleOutdir}/picardmetrics
        java -Dpicard.useLegacyParser=false -jar ${PICARDPATH}/picard.jar CollectMultipleMetrics -INPUT ${sampleOutdir}/${name}.${mappedgenome}.bam -REFERENCE_SEQUENCE ${referencefasta} -OUTPUT ${sampleOutdir}/picardmetrics/${name}.${mappedgenome} -PROGRAM CollectGcBiasMetrics -VERBOSITY WARNING
        #Can use CollectHsMetrics but needs bait coordinates
        
        echo
        date
        #Need to make second pass as CollectWgsMetrics can't be run by CollectMultipleMetrics
        #Exclude flag 512 rather than it's default MAPQ/BQ filters
        samtools view -h -u ${samflags} ${sampleOutdir}/${name}.${mappedgenome}.bam | java -Dpicard.useLegacyParser=false -jar ${PICARDPATH}/picard.jar CollectWgsMetrics -INPUT /dev/stdin -REFERENCE_SEQUENCE ${referencefasta} -OUTPUT ${sampleOutdir}/picardmetrics/${name}.${mappedgenome}.wgsmetrics -VERBOSITY WARNING -COUNT_UNPAIRED=true
        
        sequencedbases=`awk -v col="PF_HQ_ALIGNED_Q20_BASES" -F "\t" 'BEGIN {OFS="\t"} $1=="CATEGORY" {for(i=1; i<=NF; i++) {if($i==col) {colnum=i; next}}} $1=="PAIR" || $1=="UNPAIRED" {sum+=$colnum} END {print sum}' ${sampleOutdir}/picardmetrics/${name}.${mappedgenome}.alignment_summary_metrics`
        #BUGBUG includes non-mappable regions, alt contigs etc.; mappable hg38 is about 82% of chrom sizes length
        genomelen=`awk -F "\t" 'BEGIN {OFS="\t"} {sum+=$2} END {print sum}' ${chromsizes}`
        genomecov=`echo "${sequencedbases}/${genomelen}" | bc -l -q`
        genomecov=$(round ${genomecov} 2)
        echo
        echo "*** Coverage Stats ***"
        echo -e "Sequenced_bases\t${sequencedbases}\t${name}\t${mappedgenome}"
        echo -e "Ref_genome_length\t${genomelen}\t${name}\t${mappedgenome}"
        echo -e "Genomic_coverage\t${genomecov}\t${name}\t${mappedgenome}"
        
        echo
        echo "UCSC track links (actual tracks will be generated in parallel)"
        date
        #Print track links here for convenience even if the files are not created yet
        echo
        echo "Making coverage track"
        echo "track name=${name}-cov description=\"${name} ${genomecov}x genomic coverage (${analyzedReadsM}M analyzed reads) )\" maxHeightPixels=30 color=$trackcolor viewLimits=0:500 on=off visibility=full type=bigWig bigDataUrl=${UCSCbaseURL}/${name}.${mappedgenome}.coverage.bw"
        
        echo
        echo "Making VCF track"
        echo "track type=vcfTabix name=${name}-vcf description=\"${name} VCF (${analyzedReadsM}M nonredundant reads)\" visibility=pack bigDataUrl=${UCSCbaseURL}/${name}.${mappedgenome}.filtered.vcf.gz"
        
        echo "Making variant track"
        echo "track type=bigBed name=${name}-variants description=\"${name} variants (${analyzedReadsM}M nonredundant reads)\" visibility=pack bigDataUrl=${UCSCbaseURL}/${name}.${mappedgenome}.variants.bb"
    fi
elif [[ "${analysisCommand}" == "dnase" ]] || [[ "${analysisCommand}" == "atac" ]] || [[ "${analysisCommand}" == "chipseq" ]]; then
    echo
    echo "Analyzing ${analysisCommand} data"
    date
    
    echo
    echo "Making density track"
    
    #Make density track of number of overlapping reads per 150-bp window
    #BUGBUG double counts fragments where both reads are in window
    cat ${chromsizes} | 
    awk -F "\t" 'BEGIN {OFS="\t"} $1!="chrM" && $1!="chrEBV"' |
    awk '{OFS="\t"; $3=$2; $2=0; print}' | sort-bed - | cut -f1,3 | awk -v step=20 -v binwidth=150 'BEGIN {OFS="\t"} {for(i=0; i<=$2-binwidth; i+=step) {print $1, i, i+binwidth, "."} }' | 
    #--faster is ok since we are dealing with bins and read starts
    bedmap --faster --skip-unmapped --delim "\t" --bp-ovr 1 --echo --count - ${sampleOutdir}/${name}.${mappedgenome}.reads.starch |
    #resize intervals down from full bin width to step size
    #Intervals then conform to Richard's convention that the counts are reported in 20bp windows including reads +/-75 from the center of that window
    awk -v step=20 -v binwidth=150 -F "\t" 'BEGIN {OFS="\t"} {offset=(binwidth-step)/2 ; $2+=offset; $3-=offset; print}' |
    #Normalizes the density to 1M reads
    awk -v analyzedReads=${analyzedReads} -F "\t" 'BEGIN {OFS="\t"} {$5=$5/analyzedReads*1000000; print}' |
    tee $TMPDIR/${name}.${mappedgenome}.density.bed |
    #Remember wig is 1-indexed
    awk 'lastChr!=$1 {print "fixedStep chrom=" $1 " start=" $2+1 " step=" $3-$2 " span=" $3-$2; lastChr=$1} {print $5}' > $TMPDIR/${name}.${mappedgenome}.wig
    
    starch $TMPDIR/${name}.${mappedgenome}.density.bed > ${sampleOutdir}/${name}.${mappedgenome}.density.starch
    
    #Kent tools can't use STDIN
    wigToBigWig $TMPDIR/${name}.${mappedgenome}.wig ${chromsizes} ${sampleOutdir}/${name}.${mappedgenome}.bw
    
    echo "track name=${name}-dens description=\"${name} ${ucscTrackDescriptionDataType} Density (${analyzedReadsM}M analyzed reads; normalized to 1M)\" maxHeightPixels=30 color=$trackcolor viewLimits=0:3 autoScale=off visibility=full type=bigWig bigDataUrl=${UCSCbaseURL}/${name}.${mappedgenome}.bw"
    
    
    if [[ "${analysisCommand}" != "chipseq" ]]; then
        echo
        echo "Making cut count track"
        unstarch ${sampleOutdir}/${name}.${mappedgenome}.reads.starch |
        awk -F "\t" 'BEGIN {OFS="\t"} $1!="chrEBV"' |
        awk -F "\t" 'BEGIN {OFS="\t"} {if($6=="+") {chromStart=$2; chromEnd=$2+1} else {chromStart=$3; chromEnd=$3+1} print $1, chromStart, chromEnd }' | sort-bed - | tee $TMPDIR/${name}.${mappedgenome}.cuts.bed | 
        bedops --chop 1 - | awk -F "\t" 'BEGIN {OFS="\t"} {$4="id-" NR; print}' > $TMPDIR/${name}.${mappedgenome}.cuts.loc.bed
        bedmap --faster --delim '\t' --bp-ovr 1 --echo --count $TMPDIR/${name}.${mappedgenome}.cuts.loc.bed $TMPDIR/${name}.${mappedgenome}.cuts.bed | 
        awk -v analyzedReads=${analyzedReads} -F "\t" 'BEGIN {OFS="\t"} {$5=$5/analyzedReads*100000000; print}' |
        starch - > ${sampleOutdir}/${name}.${mappedgenome}.perBase.starch
        
        #Skip chrM since UCSC doesn't like the cut count to the right of the last bp in a chromosome
        unstarch ${sampleOutdir}/${name}.${mappedgenome}.perBase.starch | cut -f1-3,5 | awk -F "\t" 'BEGIN {OFS="\t"} $1!="chrM"' > $TMPDIR/${name}.${mappedgenome}.perBase.bedGraph
        #Kent tools can't use STDIN
        bedGraphToBigWig $TMPDIR/${name}.${mappedgenome}.perBase.bedGraph ${chromsizes} ${sampleOutdir}/${name}.${mappedgenome}.perBase.bw
        
        echo "track name=${name}-cuts description=\"${name} ${ucscTrackDescriptionDataType} cut counts (${analyzedReadsM}M analyzed reads)\" maxHeightPixels=30 color=$trackcolor viewLimits=0:3 autoScale=off visibility=full type=bigWig bigDataUrl=${UCSCbaseURL}/${name}.${mappedgenome}.perBase.bw"
        
        
        #echo "Making fragment coverage track"
        #samtools view ${samflags} ${sampleOutdir}/${name}.${mappedgenome}.bam | sam2bed --do-not-sort | 
        #awk -F "\t" 'BEGIN {OFS="\t"} $11>0 {print $1, $2, $2+$11}' | sort-bed - |
        #bedmap --delim '\t' --echo --count $TMPDIR/${name}.${mappedgenome}.cuts.loc.bed - | starch - > ${sampleOutdir}/${name}.fragCoverage.starch
        #right metric?
        fi
    
    
    callHotspots1=1
    callHotspots2=1 
else
    callHotspots1=0
    callHotspots2=0 
fi


if ([ "${callHotspots1}" == 1 ] || [ "${callHotspots2}" == 1 ]) && [[ "${analyzedReads}" > 0 ]]; then
    echo
    echo "Will call hotspots"
    date
    
    mappableFile="/vol/isg/annotation/bed/${annotationgenome}/mappability/${annotationgenome}.K36.mappable_only.starch"
    #NB this will call hotspots only on the mammalian genome for the *_sacCer3 hybrid indices
    #BUGBUG sacCer3 genomes?
    #For shorter old Duke data
    #      mappableFile="/vol/isg/annotation/bed/${mappedgenome}/mappability/${mappedgenome}.K20.mappable_only.starch"
    
    #Neither hotspot filters out reads for -F 512
    hotspotBAM=$TMPDIR/${name}.${mappedgenome}.bam
    samtools view ${samflags} -b -1 -@ $NSLOTS ${sampleOutdir}/${name}.${mappedgenome}.bam `unstarch --list-chromosomes ${sampleOutdir}/${name}.${mappedgenome}.reads.starch | grep -E -v "hap|random|^chrUn_|_alt$|scaffold|^C\d+" | grep -v chrM | grep -v chrEBV` > ${hotspotBAM}
    
    if [ "${callHotspots1}" == 1 ]; then
        echo
        echo "Preparing hotspot V1"
        date
        outbase=`pwd`
        mkdir -p ${sampleOutdir}/hotspots
        
        #250M is too much for hotspot1, 150M takes 12-24 hrs
        if [ ${analyzedReads} -gt 100000000 ]; then
            echo "${analyzedReads} analyzed reads. Generating hotspots on subsample of 100M reads"
            date
            #hotspot uses name of bam file as output prefix
            mkdir -p $TMPDIR/hotspot1
            hotspot1BAM=$TMPDIR/hotspot1/${name}.${mappedgenome}.bam
            samtools view -b -1 -@ $NSLOTS -s `echo "100000000/${analyzedReads}" | bc -l -q` ${hotspotBAM} > ${hotspot1BAM}
        else
            echo "${analyzedReads} analyzed reads. Generating hotspots on all reads"
            hotspot1BAM=${hotspotBAM}
        fi
        
        echo
        echo "Calling hotspots"
        date
        #NB hotspot1 can use >40GB memory for some large datasets
        #Force creation of new density file (ours is normalized)
        hotspot1Dens=$TMPDIR/${name}.${mappedgenome}.hotspotDensity.starch
        cd ${sampleOutdir}/hotspots
        ${src}/callHotspots.sh ${hotspot1BAM} ${hotspot1Dens} ${outbase}/${sampleOutdir}/hotspots ${annotationgenome} ${mappableFile} ${chromsizes} > ${outbase}/${sampleOutdir}/hotspots/${name}.${mappedgenome}.log 2>&1
        
        #Check for results first to avoid nonzero exit code
        if grep -q -i -E "(error)" ${outbase}/${sampleOutdir}/hotspots/${name}.${mappedgenome}.log; then
            cat ${outbase}/${sampleOutdir}/hotspots/${name}.${mappedgenome}.log | grep -i -E "(error)"
        fi
        
        cd ../..
        
        
        hotspotfile=${sampleOutdir}/hotspots/${name}.${mappedgenome}-final/${name}.${mappedgenome}.fdr0.01.hot.bed
        echo "Hotspots for UCSC browser"
        if [ -s "${hotspotfile}" ]; then
            cut -f1-3 ${hotspotfile} > $TMPDIR/${name}.${mappedgenome}.fdr0.01.hot.bed
            bedToBigBed -type=bed3 $TMPDIR/${name}.${mappedgenome}.fdr0.01.hot.bed ${chromsizes} ${sampleOutdir}/hotspots/${name}.${mappedgenome}.fdr0.01.hot.bb
        else
            echo "WARNING could not find ${hotspotfile} to make bigBed"
        fi
    
        peakfile=${sampleOutdir}/hotspots/${name}.${mappedgenome}-final/${name}.${mappedgenome}.fdr0.01.pks.bed
        if [ -s "${peakfile}" ]; then
            cut -f1-3 ${peakfile} > $TMPDIR/${name}.${mappedgenome}.fdr0.01.pks.bed
            bedToBigBed -type=bed3 $TMPDIR/${name}.${mappedgenome}.fdr0.01.pks.bed ${chromsizes} ${sampleOutdir}/hotspots/${name}.${mappedgenome}.fdr0.01.pks.bb
        else
            echo "WARNING could not find ${peakfile} to make bigBed"
        fi
        
        spotout=${sampleOutdir}/hotspots/${name}.${mappedgenome}.spot.out
        
        
        #subsample for spot 
        if [ ${analyzedReads} -gt 10000000 ]; then
            echo
            echo "${analyzedReads} analyzed reads. Calculating SPOT score on subsample of 10M reads"
            date
            #Could use ${hotspot1bam} (which is capped at 100M reads but would need to change calculation
            samtools view -b -1 -@ $NSLOTS -s `echo "10000000/${analyzedReads}" | bc -l -q` ${hotspotBAM} > $TMPDIR/${name}.${mappedgenome}.10Mreads.bam
            
            mkdir -p $TMPDIR/${name}.${mappedgenome}.hotspots.10Mreads
            cd $TMPDIR/${name}.${mappedgenome}.hotspots.10Mreads
            
            #NB dens track doesn't exist yet
            ${src}/callHotspots.sh $TMPDIR/${name}.${mappedgenome}.10Mreads.bam $TMPDIR/${name}.${mappedgenome}.10Mreads.density.starch $TMPDIR/${name}.${mappedgenome}.hotspots.10Mreads ${annotationgenome} ${mappableFile} ${chromsizes} > ${outbase}/${sampleOutdir}/hotspots/${name}.${mappedgenome}.10Mreads.log 2>&1
            
            #Check for results first to avoid nonzero exit code
            if grep -q -i -E "(error)" ${outbase}/${sampleOutdir}/hotspots/${name}.${mappedgenome}.10Mreads.log; then
                cat ${outbase}/${sampleOutdir}/hotspots/${name}.${mappedgenome}.10Mreads.log | grep -i -E "(error)"
            fi
            
            spotout=$TMPDIR/${name}.${mappedgenome}.hotspots.10Mreads/${name}.${mappedgenome}.10Mreads.spot.out
            
            #NB otherwise prints pwd
            cd - > /dev/null
        fi
        
        echo "Done calling hotspots"
        date
        echo
    fi
    
    
    
    if [ "${callHotspots2}" == 1 ]; then
        #TODO does not exist for cegsvectors
        hotspot2centersites="/vol/isg/annotation/bed/${annotationgenome}/hotspots2/${annotationgenome}.CenterSites.starch"
        #hotspot2.sh seems to require a 0 in column 2, unlike the UCSC standard
        hotspot2chromsizes="/vol/isg/annotation/bed/${annotationgenome}/hotspots2/${annotationgenome}.chrom.sizes"
        
        if [[ ! -s "${hotspot2centersites}" ]]; then
            echo "WARNING: could not find CenterSites file for genome ${mappedgenome}. Skipping hotspot2"
            hotspot2file=""
            hotspot2fileFDR05=""
        else
            FDRhot2="0.05"
            if (( $(bc -l <<<"${FDRhot2} >=0.05") )); then
                #NB hotspot2 calls all hotspots for one FDR threshold. Further filtering is done through the for loop below
                echo
                echo "Calling hotspot2"
                date
                mkdir -p ${sampleOutdir}/hotspot2
                
                set +e
                
                
                if [ -s "${mappableFile}" ];then
                    hotspot2mappableFileArg="-M ${mappableFile}"
                else
                    echo "WARNING: could not find ${mappableFile}, calling hotspots without mappability info"
                    hotspot2mappableFileArg=""
                fi
                
                hotspot2.sh -c /vol/isg/annotation/bed/${annotationgenome}/hotspots2/${annotationgenome}.chrom.sizes -C ${hotspot2centersites} -F ${FDRhot2} -f ${FDRhot2} ${hotspot2mappableFileArg} ${hotspotBAM} ${sampleOutdir}/hotspot2 > ${sampleOutdir}/hotspot2/${name}.${mappedgenome}.log 2>&1
                
                #Check for results first to avoid nonzero exit code
                if grep -q -i -E "(error|warning)" ${sampleOutdir}/hotspot2/${name}.${mappedgenome}.log; then
                    cat ${sampleOutdir}/hotspot2/${name}.${mappedgenome}.log | grep -i -E "(error|warning)"
                fi
                
                for FDR in {0.05,0.01,0.005,0.001}; do
                    echo "Thresholding hotspots at FDR $FDR"
                    hsmerge.sh -f ${FDR} ${sampleOutdir}/hotspot2/${name}.${mappedgenome}.allcalls.starch ${sampleOutdir}/hotspot2/${name}.${mappedgenome}.hotspots.fdr${FDR}.starch
                done
                
                echo "Hotspots for UCSC browser"
                #Hotspots
                hotspot2file=${sampleOutdir}/hotspot2/${name}.${mappedgenome}.hotspots.fdr0.01.starch
                hotspot2fileFDR05=${sampleOutdir}/hotspot2/${name}.${mappedgenome}.hotspots.fdr0.05.starch
                if [ -s "${hotspot2fileFDR05}" ] && [ `unstarch --elements ${hotspot2fileFDR05}` -gt 0 ]; then
                    unstarch ${hotspot2fileFDR05} | cut -f1-4 > $TMPDIR/${name}.${mappedgenome}.hotspots.fdr0.05.bed
                    bedToBigBed -type=bed4 $TMPDIR/${name}.${mappedgenome}.hotspots.fdr0.05.bed ${chromsizes} ${sampleOutdir}/hotspot2/${name}.${mappedgenome}.hotspots.fdr0.05.bb
                else
                    echo "WARNING could not find ${hotspot2fileFDR05} to make bigBed"
                fi
                
                #Peaks
                hotspot2peakfile=${sampleOutdir}/hotspot2/${name}.${mappedgenome}.peaks.starch
                if [ -s "${hotspot2peakfile}" ] && [ `unstarch --elements ${hotspot2peakfile}` -gt 0 ]; then
                    unstarch ${hotspot2peakfile} | cut -f1-4 > $TMPDIR/${name}.${mappedgenome}.peaks.bed
                    bedToBigBed -type=bed4 $TMPDIR/${name}.${mappedgenome}.peaks.bed ${chromsizes} ${sampleOutdir}/hotspot2/${name}.${mappedgenome}.peaks.bb
                else
                    echo "WARNING could not find ${hotspot2peakfile} to make bigBed"
                fi
                
                set -e
                
                #hotspot2 generates a series of files redundant with our own tracks:
                #   density - bins start at 0 instead of 65, counts raw tags, includes 0 bins
                #   cutcounts - not 1M-read normalized
                rm -f ${sampleOutdir}/hotspot2/${name}.${mappedgenome}.density.starch ${sampleOutdir}/hotspot2/${name}.${mappedgenome}.density.bw
                #could also delete this but it's needed for SPOT calculation
                #rm -f ${sampleOutdir}/hotspot2/${name}.${mappedgenome}.cutcounts.starch
                
                echo "Done calling hotspot2"
                date
                echo
            else 
                echo "WARNING: hotspot2 should be run with FDR >= 0.05"
                echo "Skipping hotspot2"
                echo
                hotspot2file=""
                hotspot2fileFDR05=""
            fi
        fi
    fi
fi


#Stats
echo
echo "*** Overall Stats ***"
echo -e "Num_sequenced_reads\t${sequencedReads}\t\t${name}\t${mappedgenome}"
printfNA "Num_pass_filter_alignments\t${PFalignments}\t%.1f%%\t${name}\t${mappedgenome}\n" "${pctPFalignments}"
printfNA "Num_mitochondria_reads\t${numMappedReadsMitochondria}\t%.1f%%\t${name}\t${mappedgenome}\n" "${pctMappedReadsMitochondria}"
printfNA "Num_analyzed_reads\t${analyzedReads}\t%.1f%%\t${name}\t${mappedgenome}\n" "${pctanalyzedReads}"
printfNA "Num_duplicate_reads\t${dupReads}\t%.1f%%\t${name}\t${mappedgenome}\n" "$pctdupReads"

#Don't have denominator of unpaired reads we tried to map, so don't compute % for first
echo -e "Num_SE_pass_filter_alignments\t${PFalignmentsSE}\t\t${name}\t${mappedgenome}"
#NB denominator is PF reads, not pctMappedReadsMitochondrias sequenced
printfNA "Num_SE_analyzed_reads\t${analyzedReadsSE}\t%.1f%%\t${name}\t${mappedgenome}\n" "${pctanalyzedReadsSE}"


if [[ "${analysisCommand}" == "callsnps" ]]; then
    #Repeat it here for ease of parsing
    echo -e "Genomic_coverage\t\t${genomecov}\t${name}\t${mappedgenome}"
fi

if [[ "${callHotspots1}" == 1 ]]; then
    if [ -s "${hotspotfile}" ]; then
        numhotspots=`cat ${hotspotfile} | wc -l`
    else
        numhotspots="NA"
    fi
    echo -e "Num_hotspots\t${numhotspots}\t\t${name}\t${mappedgenome}"
    
    
    if [ -s "${spotout}" ]; then
        spot=`cat ${spotout} | awk 'BEGIN {OFS="\t"} NR==2 {print $3}'`
    else
        spot="NA"
    fi
    echo -e "SPOT\t${spot}\t\t${name}\t${mappedgenome}"
fi

if [[ "${callHotspots2}" == 1 ]]; then
    if [ -s "${hotspot2file}" ]; then
        numhotspots2=`unstarch --elements ${hotspot2file}`
    else
        numhotspots2="NA"
    fi
    echo -e "Num_hotspots2\t${numhotspots2}\t\t${name}\t${mappedgenome}"
    
    
    if [ -s "${hotspot2fileFDR05}" ] && [ "${numhotspots2}" -gt 0 ]; then
        #I believe you can't directly use .SPOT.txt since it is set to the maximum FDR chosen
        numTotalReadsHotspot2=$(unstarch ${sampleOutdir}/hotspot2/${name}.${mappedgenome}.cutcounts.starch | awk -F "\t" 'BEGIN {OFS="\t"; sum=0} {sum+=$5} END {print sum}')
        if [ "${numTotalReadsHotspot2}" -gt 0 ]; then
            numReadsinDHSHotspot2=$(unstarch ${sampleOutdir}/hotspot2/${name}.${mappedgenome}.cutcounts.starch | bedops -e -1 - ${hotspot2fileFDR05} | awk -F "\t" 'BEGIN {OFS="\t"; sum=0} {sum+=$5} END {print sum}')
            spot2=$(echo "scale=4; ${numReadsinDHSHotspot2}/${numTotalReadsHotspot2}" | bc)
        else
            spot2="NA"
        fi
    else
        spot2="NA"
    fi
    echo -e "SPOT2\t${spot2}\t\t${name}\t${mappedgenome}"
fi


echo
echo "Histogram of mismatches"
samtools view ${samflags} ${sampleOutdir}/${name}.${mappedgenome}.bam | awk -F "\t" 'BEGIN {OFS="\t"} { \
    for(i=12; i<=NF; i++) { \
         if(match($i, /NM:i:/)) { \
            print substr($i, 6); \
         } \
    } \
}' | awk '{names[$0]++} END {for (cur in names) {print names[cur], cur}}' | sort -k2,2g

#NB XA reads are computed for the unpaired reads, while MAPQ reflects the final PE location
echo
echo "Histogram of number of best alignments"
samtools view ${samflags} ${sampleOutdir}/${name}.${mappedgenome}.bam | awk -F "\t" 'BEGIN {OFS="\t"} { \
    read="NA"; \
    for(i=12; i<=NF; i++) { \
         if(match($i, /X0:i:/)) { \
            read=substr($i, 6); \
#            if(read>1) {print} \
         } \
#         if(match($i, /XA:Z:/)) { \
#            read=length(split(substr($i, 6), xa, ";")); \
#         } \
    } \
    print read; \
}' | awk -F "\t" 'BEGIN {OFS="\t"} $0>=10 {$0="10+"} {print}' | awk '{names[$0]++} END {for (cur in names) {print names[cur], cur}}' | sort -k2,2g


PEreads=`cat $TMPDIR/${name}.${mappedgenome}.flagstat.txt | grep "paired in sequencing" | awk '{print $1+$3}'`
if [[ "${PEreads}" > 0 ]] && [[ "${mappedgenome}" != "cegsvectors" ]]; then
    echo
    echo "Template lengths"
    samtools view ${samflags} ${sampleOutdir}/${name}.${mappedgenome}.bam | awk -F "\t" 'BEGIN {OFS="\t"} $3!="chrM"' | awk -v sample=${name}.${mappedgenome} -F "\t" 'BEGIN {OFS="\t"} $9>0 {print sample, $9}' | sort -k2,2n | tee ${sampleOutdir}/${name}.${mappedgenome}.insertlengths.txt | cut -f2 |
    awk -F "\t" 'BEGIN {OFS="\t"} NR==1 {print "Minimum: " $0} {cum+=$0; lengths[NR]=$0; if($0<125) {lastLineUnder125=NR}} END {print "Number of reads: " NR; print "25th percentile: " lengths[int(NR*0.25)]; print "50th percentile: " lengths[int(NR*0.5)]; print "75th percentile: " lengths[int(NR*0.75)]; print "95th percentile: " lengths[int(NR*0.95)]; print "99th percentile: " lengths[int(NR*0.99)]; print "Maximum: " $0; print "Mean: " cum/NR; print "Prop. reads under 125 bp: " lastLineUnder125/NR}'
    gzip -f ${sampleOutdir}/${name}.${mappedgenome}.insertlengths.txt
fi

echo
echo "read lengths:"
samtools view ${samflags} ${sampleOutdir}/${name}.${mappedgenome}.bam | cut -f10 | awk '{lengths[length($0)]++} END {for (l in lengths) {print lengths[l], l}}' | sort -k2,2g


echo
echo "read count by sequencing instrument"
#The sort -k1,1nr | head -100 caps the number of unique entries (taking the most common), helpful if we have not parsed the read name accurately
bam2instrument ${samflags} ${sampleOutdir}/${name}.${mappedgenome}.bam | awk '{names[$0]++} END {for (cur in names) {print names[cur], cur}}' | sort -k1,1nr | head -100 | sort -k2,2


echo
echo -e "\nDone!"
date
