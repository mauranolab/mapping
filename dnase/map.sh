#!/bin/bash
set -e -o pipefail


#Above not catching segfaults
#https://unix.stackexchange.com/questions/24307/how-can-i-trap-a-program-that-returns-139-segmentation-fault-in-bash
#not working
#https://unix.stackexchange.com/questions/24307/how-can-i-trap-a-program-that-returns-139-segmentation-fault-in-bash
#neither set -bm or set -o monitor help
#trap 'if [[ $? -eq 139 ]]; then echo "segfault !”; exit 1; fi' CHLD

#BTW can parse core dump with 'objdump -s core' or 'gdb prog core' (if you know prog)

date

permittedMismatches=3


#celltype is called sample in BWAmerge and submitBWA
celltype=$1
DS=$2
genome=$3
src=$4
echo "celltype: $celltype, DS: $DS, Genome: $genome"
shift
shift
shift
shift
userAlnOptions=$@


bam2instrument()
{
    cut -f1 | cut -d ":" -f1 | 
    #Hack to clean up some encodeproject.org data that has underscore in place of colon after sequencer name
    perl -pe 's/_\d+_\d+_\d+_\d+$//g;' | 
    #SRR data
    perl -pe 's/^(SRR\d+)\.\d+$/$1/g;'
}


jobid=$SGE_TASK_ID
readsFq=`cat inputs.txt | awk -v jobid=$jobid 'NR==jobid'`
if [ ! -f "$readsFq" ]; then
    echo "ERROR: Can not find file $readsFq"
    exit 1
fi
echo "Will process reads file $readsFq"


sample1=`echo $readsFq | awk '{print $2}'`
if [[ "$sample1" == "" ]] ; then
    sample1=`basename $readsFq | perl -pe 's/.fa(stq)?(.gz)?$//g;'`
fi


if [[ "$sample1" =~ "_R2" ]]; then
    echo "Won't process R2 file -- $sample1" > /dev/stderr 
    exit 0
fi


#NB needs to be on NFS if you want to run fastqc in separate job
#note sample1 is not unique (doesn't contain FC)
#TMPDIR=tmp/$sample1.$jobid
#mkdir -p $TMPDIR
echo "using $TMPDIR as TMPDIR"


echo
echo "Configuring trimming parameters"
#/home/jvierstra/proj/code/bio-tools/apps/trim-adapters-illumina/trim-adapters-illumina DS32747A_CTTGTA_L002_R1_001.fastq.gz DS32747A_CTTGTA_L002_R2_001.fastq.gz R1.jtrim.fastq.gz R2.jtrim.fastq.gz


#Trimmomatic options
#TODO Probably need different sequences per barcode. Note this fa file has 2 ident copies of left adapter and none of right adapter (with barcode).
#Regular illumina dsDNA protocol


#analysisType can either be DNase or ATAC
analysisType="DNase"
echo "Using $analysisType adapters for trimming"
if [ "$analysisType" == "DNase" ]; then 
    illuminaAdapters="/cm/shared/apps/trimmomatic/0.38/adapters/TruSeq3-PE-2.fa"
elif [ $analysisType == "ATAC" ]; then 
    illuminaAdapters="/cm/shared/apps/trimmomatic/0.38/adapters/NexteraPE-PE.fa"
else 
    echo "ERROR specify adapters"
    exit 2
fi

seedmis=2
#Pretty much anything below 10 works
PEthresh=5
SEthresh=5
mintrim=1
keepReverseReads=true
trimmomaticBaseOpts="-threads $NSLOTS -trimlog $TMPDIR/${sample1}.trim.log.txt"
trimmomaticSteps="TOPHRED33 ILLUMINACLIP:$illuminaAdapters:$seedmis:$PEthresh:$SEthresh:$mintrim:$keepReverseReads"
#MAXINFO:27:0.95 TRAILING:20

#For shorter old Duke data
#Check if samples contain DUKE adapter (TCGTATGCCGTCTTC) and trim to 20bp if more than 25% of reads do
#sequencedTags=$(zcat $readsFq | awk 'NR%4==2' | wc -l)
#BUGBUG - I think Jesper just took what looks like a readthrough sequence (TCGTATGCCGTCTTC). Not sure why the trimmer isn't dealing with this properly
#if [ `zcat $readsFq | awk -v thresh=0.25 -v sequencedTags=$sequencedTags 'NR%4==2 && $1~/TCGTATGCCGTCTTC/ {readsWithDukeSequence+=1} END {if (readsWithDukeSequence/sequencedTags>thresh) {print 1} else {print 0}}'` ]; then
#    echo "More than 25% of reads have DUKE sequence (TCGTATGCCGTCTTC) - Hard clip to 20bp reads"
#    trimmomaticSteps="CROP:20 $trimmomaticSteps"
#    minMAPQ=10
#else
#    echo "No DUKE sequence present"
    trimmomaticSteps=" $trimmomaticSteps MINLEN:27 CROP:36"
    minMAPQ=20
#fi


#BUGBUG For aggregate submission e.g. submitting multiple flowcells at once, there will be a collision if we have sequenced same sample on multiple flowcells. 
fc=""

#BUGBUG any older data w/o R1/R2 convention?
sample2=`echo $sample1 | perl -pe 's/_R1(_\d+)?$/_R2$1/g;'`
#BUGBUG note grep -q $fc dies on two PATSKI samples but they are SE
if echo "$sample1" | grep -q R1 && echo "$sample2" | grep -q R2 && grep "$sample2" inputs.txt | grep -q "$fc" ; then
    echo "Found R2 $sample2"
    if [ `grep "$sample2" inputs.txt | grep "$fc" | wc -l` -gt 1 ]; then
        echo "ERROR: Multiple R2 files found -- are there duplicate entries in inputs.txt?"
        exit 3
    fi
    
    reads2fq=`grep "$sample2" inputs.txt | grep "$fc"`
    if [ ! -f "$reads2fq" ]; then
        echo "ERROR: Can not find R2 file $reads2fq"
        exit 4
    fi
    
    echo "Will process reads file $reads2fq"
    
    PErun="TRUE"
    #NB provisional -- sample will be updated later with FC name
    sample=`echo $sample1 | perl -pe 's/_R1(_\d+)?/_R1R2\1/g;'`
    sample="${fc}$sample"
    
    echo "Filtering out reads with >75% G content"
    #TODO could potentially save R1 where only R2 is dark once it can handle single read. Could patch trimmomatic instead?
    $src/filterNextSeqReadsForPolyG.py --inputfileR1 $readsFq --inputfileR2 $reads2fq --outputfileR1 $TMPDIR/$sample1.pretrim.fastq.gz --outputfileR2 $TMPDIR/$sample2.pretrim.fastq.gz --maxPolyG 75
    
    java org.usadellab.trimmomatic.TrimmomaticPE $trimmomaticBaseOpts $TMPDIR/$sample1.pretrim.fastq.gz $TMPDIR/$sample2.pretrim.fastq.gz $TMPDIR/$sample1.fastq $TMPDIR/$sample1.unpaired.fastq $TMPDIR/$sample2.fastq $TMPDIR/$sample2.unpaired.fastq $trimmomaticSteps
    #TODO why does this output uncompressed fastq? I think it's just so one can test with -s below
    
    echo -n "Unpaired reads:"
    #Merge anything unpaired from either R1 or R2
    cat $TMPDIR/$sample1.unpaired.fastq $TMPDIR/$sample2.unpaired.fastq | tee $TMPDIR/$sample.unpaired.fastq | wc -l
    
    
    if [ ! -s "$TMPDIR/$sample1.fastq" ] && [ ! -s "$TMPDIR/$sample2.fastq" ]; then
        echo "No tags passed filtering, quitting successfully"
        exit 0
    fi
    
    
    mkdir -p fastqc
    fastQcOutdir="fastqc/${fc}${sample1}_qc"
    if [ ! -d "$fastQcOutdir" ]; then
        #qsub -cwd -V -N ${sample1}.qc -o fastqc/${sample1}.qc -S /bin/bash -j y -b y -p -500 "mkdir -p $fastQcOutdir; fastqc --outdir $fastQcOutdir $TMPDIR/$sample1.fastq"
        mkdir -p $fastQcOutdir; fastqc --outdir $fastQcOutdir $TMPDIR/$sample1.fastq
    fi
    
    fastQcOutdir="fastqc/${fc}${sample2}_qc"
    if [ ! -d "$fastQcOutdir" ]; then
        #qsub -cwd -V -N ${sample2}.qc -o fastqc/${sample1}.qc -S /bin/bash -j y -b y -p -500 "mkdir -p $fastQcOutdir; fastqc --outdir $fastQcOutdir $TMPDIR/$sample2.fastq"
        mkdir -p $fastQcOutdir; fastqc --outdir $fastQcOutdir $TMPDIR/$sample2.fastq
    fi
    
    echo
    echo "Histogram of read lengths for R1"
    cat $TMPDIR/$sample1.fastq | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' | sort -k1,1n
    
    echo
    echo "Histogram of read lengths for R1"
    cat $TMPDIR/$sample2.fastq | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' | sort -k1,1n
    
    echo
    echo "Histogram of read lengths for unpaired"
    cat $TMPDIR/$sample.unpaired.fastq | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' | sort -k1,1n
else
    PErun="FALSE"
    sample="${fc}$sample1"
    
    #TODO missing filterNextSeqReadsForPolyG.py
    
    #BUGBUG wrong adapter files
    java org.usadellab.trimmomatic.TrimmomaticSE $trimmomaticBaseOpts $readsFq $TMPDIR/$sample1.fastq $trimmomaticSteps
    
    
    if [ ! -s "$TMPDIR/$sample1.fastq" ]; then
        echo "No tags passed filtering, quitting successfully"
        exit 0
    fi
    
    
    mkdir -p fastqc
    fastQcOutdir="fastqc/${fc}${sample1}_qc"
    if [ ! -d "$fastQcOutdir" ]; then
        #qsub -cwd -V -N ${sample1}.qc -o fastqc/${sample1}.qc -S /bin/bash -j y -b y -p -500 "mkdir -p $fastQcOutdir; fastqc --outdir $fastQcOutdir $TMPDIR/$sample1.fastq"
        mkdir -p $fastQcOutdir; fastqc --outdir $fastQcOutdir $TMPDIR/$sample1.fastq
    fi
    
    echo
    echo "Histogram of read lengths for R1"
    cat $TMPDIR/$sample1.fastq | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' | sort -k1,1n
fi


date



echo
echo "Mapping $readsFq for $sample1"
echo "userAlnOptions=$userAlnOptions"

sampleOutdir=${celltype}-${DS}.${genome}
mkdir -p $sampleOutdir


#not sure what -R is, making it lower than samse/pe -n reduces mapped PE tags but not SE tags
#-Y filters sequences with \d+:Y:... after the space in the read name
#Originally -n 0.04 seemed to allow upto two mismatches at 36 bp (must have rounded up)
bwaAlnOpts="-n $permittedMismatches -l 32 $userAlnOptions -t $NSLOTS -Y"

#Other options notes:
#-q 0.20 does soft clip quality-based trimming of 3' end of reads, but only down to 35 bp
#http://seqanswers.com/forums/showthread.php?t=5628
#http://seqanswers.com/forums/showthread.php?t=6251


genomesToMap="$genome"
echo "Will map to genomes $genomesToMap"


#NB am losing about 15" to load index when submit multiple jobs
bwaIndexBase=/vol/isg/annotation/bwaIndex
for curGenome in $genomesToMap; do
    echo
    echo "Mapping to genome for genome $curGenome"
    case "$curGenome" in
    hg19)
        bwaIndex=$bwaIndexBase/hg19all/hg19all;;
    hg38)
#        bwaIndex=$bwaIndexBase/hg38all/hg38all;;
        bwaIndex=$bwaIndexBase/hg38_no_alt_analysis_set/hg38_no_alt_analysis_set;;
    mm10)
#        bwaIndex=$bwaIndexBase/mm10all/mm10all;;
        bwaIndex=$bwaIndexBase/mm10_no_alt_analysis_set/mm10_no_alt_analysis_set;;
    *)
        echo "Don't recognize genome $curGenome";
        exit 5;;
    esac
    
    echo "bwa aln $bwaAlnOpts $bwaIndex ..."
    bwa aln $bwaAlnOpts $bwaIndex $TMPDIR/$sample1.fastq > $TMPDIR/$sample1.$curGenome.sai

    date


    #Previously used -n 10 but never really used XA tag and maybe was causing sampe to occasionally truncate last line of output (dropping the tags)
    #TODO perhaps I should add -N to get complete set of hits?
    #Read group should properly be FC name, lane and barcode but we don't have the latter info
    DS_nosuffix=`echo $DS | perl -pe 's/[A-Z]$//g;'`
    bwaExtractOpts="-n 3 -r @RG\\tID:${sample}\\tLB:$DS\\tSM:${DS_nosuffix}"
    if [[ "$PErun" == "TRUE" ]] ; then
        echo -e "\nMapping R2 $reads2fq for $sample2"
        bwa aln $bwaAlnOpts $bwaIndex $TMPDIR/$sample2.fastq > $TMPDIR/$sample2.$curGenome.sai
        date
        
        #-P didn't have a major effect, but some jobs were ~10-40% faster but takes ~16GB RAM instead of 4GB
        extractcmd="sampe $bwaExtractOpts -a 500 $bwaIndex $TMPDIR/$sample1.$curGenome.sai $TMPDIR/$sample2.$curGenome.sai $TMPDIR/$sample1.fastq $TMPDIR/$sample2.fastq"
        
        #Only map unpaired reads if the file nonzero
        if [ -s "$TMPDIR/$sample.unpaired.fastq" ]; then
            echo -e "\nMapping unpaired $sample.unpaired.fastq for $sample1"
            bwa aln $bwaAlnOpts $bwaIndex $TMPDIR/$sample.unpaired.fastq > $TMPDIR/$sample.unpaired.$curGenome.sai
            date
            
            echo
            echo "Extracting unpaired reads"
            unpairedReadsSam="$TMPDIR/$sample.$curGenome.unpaired.sam"
            unpairedExtractcmd="samse $bwaExtractOpts $bwaIndex $TMPDIR/$sample.unpaired.$curGenome.sai $TMPDIR/$sample.unpaired.fastq"
            echo -e "unpairedExtractcmd=bwa $unpairedExtractcmd | (...)"
            bwa $unpairedExtractcmd | grep -v "^@" > $unpairedReadsSam
        else
            unpairedReadsSam=""
        fi
        #TODO merge headers instead of dropping
    else
        extractcmd="samse $bwaExtractOpts $bwaIndex $TMPDIR/$sample1.$curGenome.sai $TMPDIR/$sample1.fastq"
        unpairedReadsSam=""
    fi
    
    
    date
    echo
    echo "Extracting"
    echo -e "extractcmd=bwa $extractcmd | (...)"
    
    #BUGBUG wasting cores?
    #TODO use fork/wait
    bwa $extractcmd |
    cat - $unpairedReadsSam |
    \
    #calmd
    #Necessary?
    #BUGBUG I think this is the slow step
    #Fix NM/MD. Appears not to affect mpileup, but still worth it for ease of use (bwa still seems to set NM/MD wrong sporadically)
    #used to precalculate BAQ with samtools calmd -S -Abr -E, see http://www.biostars.org/p/1268
    #NB redirecting stderr since calmd is noisy, but you won't see real errors
    #samtools calmd -S -u - $fa 2> $TMPDIR/$sample.$curGenome.calmd.log |
    \
    #Use view instead of calmd
    #TODO prob better to do NSLOTS/2 or so
    #samtools view -@ $NSLOTS -S -u - | 
    samtools sort -@ $NSLOTS -m 1500M -O bam -T $TMPDIR/${sample}.sortbyname -l 1 -n - |
    $src/filter_reads.py --max_mismatches $permittedMismatches --min_mapq $minMAPQ - - |
    samtools sort -@ $NSLOTS -m 1500M -O bam -T $TMPDIR/${sample}.sort -l 1 - |
    #Add MC tag containing mate CIGAR for duplicate calling
    java -Xmx2g -jar /cm/shared/apps/picard/1.140/picard.jar FixMateInformation INPUT=/dev/stdin OUTPUT= $sampleOutdir/$sample.$curGenome.bam VERBOSITY=ERROR QUIET=TRUE COMPRESSION_LEVEL=1
    
#    echo
#    echo "Cleanup"
#    date
#    cp $sample.$curGenome.bam $TMPDIR/$sample.$curGenome.unclean.bam
#    #BUGBUG Should fix the ERROR... read errors, but doesn't do anything to first 100000 lines of test case except increment version to 1.4 "@HD   VN:1.4"
#    java -Xmx2g -jar /cm/shared/apps/picard/1.140/picard.jar/ CleanSam INPUT=$sample.$curGenome.bam OUTPUT=$sample.$curGenome.clean.bam COMPRESSION_LEVEL=1 && mv $sample.$curGenome.clean.bam $sample.$curGenome.bam
#
#rearranges flag order but doesn't fix problem either
#    java -Xmx2g -jar /cm/shared/apps/picard/1.140/picard.jar FixMateInformation INPUT=$sample.$curGenome.bam OUTPUT=$sample.$curGenome.clean.bam COMPRESSION_LEVEL=1 VALIDATION_STRINGENCY=LENIENT && mv $sample.$curGenome.clean.bam $sample.$curGenome.bam
    
    echo
    echo "SAMtools statistics for genome $genome"
    samtools flagstat $sampleOutdir/$sample.$curGenome.bam
done


echo "QC metrics to be done on un-merged data"
date
#TODO break out as separate job?


echo
echo "Mean quality by cycle"
#BUGBUG performs badly for SRR jobs -- some assumption not met?
java -Xmx3g -jar /cm/shared/apps/picard/1.140/picard.jar MeanQualityByCycle INPUT= $sampleOutdir/$sample.$curGenome.bam OUTPUT=$TMPDIR/$sample.baseQ.txt CHART_OUTPUT=$TMPDIR/$sample.baseQ.pdf VALIDATION_STRINGENCY=LENIENT

#BUGBUG slow
instrument=`samtools view $sampleOutdir/$sample.$curGenome.bam | bam2instrument | uniq | sort | uniq | perl -pe 's/\n$//g;' | perl -pe 's/\n/;/g;'`
awk -v instrument=$instrument -v fc=$fc -v sample=$sample -v ds=$DS -F "\t" 'BEGIN {OFS="\t"} $0!~/^#/ && $0!="" {if($1=="CYCLE") {$0=tolower($0); $(NF+1)="instrument\tfc\tsample\tDS"} else {$(NF+1)=instrument "\t" fc "\t" sample "\t" ds;} print}' $TMPDIR/$sample.baseQ.txt > $sampleOutdir/$sample.baseQ.txt


echo
#Prints positions
echo -n -e "Histogram of mismatches to reference in 36 bases:\t$instrument\t$fc\t$sample\t$DS\t"
samflags="-q 30 -F 1548"

#BUGBUG flag 2 not working?
samtools view $samflags $sampleOutdir/$sample.$curGenome.bam | awk -F "\t" 'BEGIN {OFS="\t"} { \
      readlength = length($10); \
      if (and($2, 16)) { \
         strand="-"; \
      } else { \
         strand="+"; \
      } \
      for(i=12; i<=NF; i++) { \
         if (match($i, /MD:Z:/)) { \
            #Hardcoded based on the length of the attribute name \
            curstart = 6; \
            curOffset = 0; \
            j = 0; \
            #A MM in the first position is preceded by 0, so we can use a single loop \
            while(match(substr($i, curstart), /^[0-9]+([ACGT]|\^[ACGT]+)/)) { \
               curBlockLength = RLENGTH; \
               #Find where the number ends \
               match(substr($i, curstart), /([ACGT]|\^[ACGT]+)/); \
               curOffset = curOffset + substr($i, curstart, RSTART-1); \
               #print NR ":cur parse offset: ", substr($i, curstart, RSTART-1), "( rstart=" RSTART ", rlength=" RLENGTH ")"; \
               \
               #Ignore indels for now \
               if ( RLENGTH == 1 ) { \
                  curOffset++; #Need to increment curOffset to count the polymorphic base (indels dont take up space) \
                  if(strand=="-") { \
                     mmOffsets[j] = readlength - (curOffset - 1); \
                  } else { \
                     mmOffsets[j] = curOffset; \
                  } \
                  # Print Line number, cycle number of MM, called base, quality\
                  print NR, mmOffsets[j], substr($10, curOffset, 1), substr($11, curOffset, 1); \
                  j++; \
               } \
               curstart = curstart + curBlockLength; \
            } \
         } \
      } \
}' | tee $TMPDIR/$sample.mm.txt |
awk -v minBaseQ=20 -F "\t" 'BEGIN { for (i=0; i<256; i++) { codeFor[sprintf("%c", i)] = i } } codeFor[$4]-33 > minBaseQ {print} ' | 
cut -f2 | sort -g | uniq -c | sort -k2,2 -g | awk 'BEGIN {ORS="\t"} {print $1}'
echo


echo
echo
echo "Gerald's call for mismatched positions"
awk 'BEGIN {OFS="\t"; print "cycle", "A", "C", "G", "T", "N"} {errors[$2,$3]++} END {for(i=1; i<=36; i++) {print i, errors[i, "A"], errors[i, "C"], errors[i, "G"], errors[i, "T"], errors[i, "N"]}}' $TMPDIR/$sample.mm.txt


echo
echo "Done!"
date
