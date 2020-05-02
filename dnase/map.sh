#!/bin/bash
set -eu -o pipefail

#Above not catching segfaults
#https://unix.stackexchange.com/questions/24307/how-can-i-trap-a-program-that-returns-139-segmentation-fault-in-bash
#not working
#https://unix.stackexchange.com/questions/24307/how-can-i-trap-a-program-that-returns-139-segmentation-fault-in-bash
#neither set -bm or set -o monitor help
#trap 'if [[ $? -eq 139 ]]; then echo "segfault !”; exit 1; fi' CHLD

#BTW can parse core dump with 'objdump -s core' or 'gdb prog core' (if you know prog)


shopt -s expand_aliases
alias bedops='bedops --ec --header'
alias bedmap='bedmap --ec --header --sweep-all'
#alias starch='starch --header'
alias closest-features='closest-features --header'


genomesToMap=$1
analysisType=$2
sampleOutdir=$3
BS=$4
src=$5

shift
shift
shift
shift
shift
userAlnOptions=$@


echo "Output directory: ${sampleOutdir}, BS: ${BS}, Genomes: ${genomesToMap}, analysisType:${analysisType}"
date


processingCommand=`echo "${analysisType}" | awk -F "," '{print $1}'`
sampleType=`echo "${analysisType}" | awk -F "," '{print $2}'`


if [[ "${sampleType}" == "dnase" ]] || [[ "${sampleType}" == "atac" ]] || [[ "${sampleType}" == "chipseq" ]]; then
    maxInsertSize=500
    permittedMismatches=3
else
    maxInsertSize=1000
    permittedMismatches="0.08"
fi


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


getReadgroup()
{
    local BS=$1
    
    local BS_nosuffix=`echo ${BS} | perl -pe 's/[A-Z]$//g;'`
    local readgroup="@RG\\tID:${fc}${BS}\\tLB:${BS}\\tSM:${BS_nosuffix}\\tPL:ILLUMINA"
    if [ -s "/vol/mauranolab/flowcells/data/${fc/./}/info.txt" ]; then
        local readgroup_instrument=`awk -F "\t" 'BEGIN {OFS="\t"} $1=="#Instrument" {print $2}' /vol/mauranolab/flowcells/data/${fc/./}/info.txt`
        
        local readgroup_date=`awk -F "\t" 'BEGIN {OFS="\t"; loaddate="NA"} $1=="#Load date" {loaddate=$2} END {print loaddate}' /vol/mauranolab/flowcells/data/${fc/./}/info.txt`
        #BUGBUG hardcoded column numbers
        local readgroup_bcs=`awk -v ds=${BS} -F "\t" 'BEGIN {OFS="\t"} $0!~/^#/ && 0!="" && $2==ds {split($6, bc1, "_"); split($7, bc2, "_"); print bc1[2] "-" bc2[2]}' /vol/mauranolab/flowcells/data/${fc/./}/info.txt`
        #BUGBUG BC: shows up in bwa command line but at some point disappears from the bam header
        readgroup="${readgroup}\\tDT:${readgroup_date}\\tBC:${readgroup_bcs}\\tPU:${fc/./}-${readgroup_bcs}"
        
        case "${readgroup_instrument}" in
        Balin)
            readgroup="${readgroup}\\tCN:Maurano_Lab\\tPM:NextSeq_500"
            ;;
        GTC_NovaSeq)
            readgroup="${readgroup}\\tCN:NYUMC_GTC\\tPM:NovaSeq_6000"
            ;;
        GTC_NextSeq)
            readgroup="${readgroup}\\tCN:NYUMC_GTC\\tPM:NextSeq_500"
            ;;
        esac
    fi
    
    echo "${readgroup}"
}


jobid=$SGE_TASK_ID
readsFq=`awk -v jobid=$jobid 'NR==jobid' ${sampleOutdir}/inputs.map.txt`
if [ ! -f "${readsFq}" ]; then
    echo "ERROR: Can not find file ${readsFq}"
    exit 1
fi
echo "Will process reads file ${readsFq}"


sample1=`echo ${readsFq} | awk '{print $2}'`
if [[ "${sample1}" == "" ]] ; then
    sample1=`basename ${readsFq} | perl -pe 's/.fa(stq)?(.gz)?$//g;'`
fi


if [[ "${sample1}" =~ "_R2" ]]; then
    echo "Won't process R2 file -- ${sample1}" > /dev/stderr 
    exit 0
fi


#NB needs to be on NFS if you want to run fastqc in separate job
#note sample1 is not unique (doesn't contain FC)
#TMPDIR=tmp/${sampleOutdir}.$jobid
#mkdir -p $TMPDIR
echo "Running on $HOSTNAME. Using $TMPDIR as tmp"

mkdir -p ${sampleOutdir}


echo
echo "Configuring trimming parameters"
#Trimmomatic options
if [[ "${sampleType}" == "atac" ]]; then 
    illuminaAdapters="/cm/shared/apps/trimmomatic/0.39/adapters/NexteraPE-PE.fa"
else 
    #TODO Probably need different sequences per barcode. Note this fa file has 2 ident copies of left adapter and none of right adapter (with barcode).
    illuminaAdapters="/cm/shared/apps/trimmomatic/0.39/adapters/TruSeq3-PE-2.fa"
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


#For old Duke data but too slow to enable by default
readsLongEnough=1
#Check if samples contain DUKE adapter (TCGTATGCCGTCTTC) and trim to 20bp if more than 25% of reads do
if [[ 0 == 1 ]]; then
    sequencedReads=$(zcat ${readsFq} | awk 'NR%4==2' | wc -l)
    #The proper adapter files are in /cm/shared/apps/trimmomatic/0.39/adapters/TruSeq2-SE.fa, but the Duke reads also have 8-9 As at the end after it reads through the adapter so it's probably better just to CROP:20
    if [ `zcat ${readsFq} | awk -v thresh=0.25 -v sequencedReads=$sequencedReads 'NR%4==2 && $1~/TCGTATGCCGTCTTC/ {readsWithDukeSequence+=1} END {if (readsWithDukeSequence/sequencedReads>thresh) {print 1} else {print 0}}'` ]; then
        echo "More than 25% of reads have DUKE sequence (TCGTATGCCGTCTTC) - Hard clip to 20bp reads"
        trimmomaticSteps="CROP:20 ${trimmomaticSteps}"
        readsLongEnough=0
    else
        echo "No DUKE sequence present"
        readsLongEnough=1
    fi
fi

if [[ "${sampleType}" == "dnase" ]] || [[ "${sampleType}" == "atac" ]] || [[ "${sampleType}" == "chipseq" ]]; then
    trimmomaticSteps="${trimmomaticSteps} CROP:36"
fi

if [ "${readsLongEnough}" -eq 1 ]; then
    trimmomaticSteps="${trimmomaticSteps} MINLEN:27"
fi

#BUGBUG a bit fragile
fc=`readlink -f ${readsFq} | xargs dirname | xargs dirname | xargs dirname | xargs basename`
if [[ ! "${fc}" =~ ^FC ]] ; then
    fc=""
else
    echo "Flowcell ${fc}"
    fc="${fc}."
fi


readgroup=$(getReadgroup ${BS})

sample2=`echo ${sample1} | perl -pe 's/_R1(_\d+)?$/_R2$1/g;'`
if echo "${sample1}" | grep -q _R1 && echo "${sample2}" | grep -q _R2 && grep "${sample2}" ${sampleOutdir}/inputs.map.txt | grep -q "${fc}" ; then
    echo "Found R2 ${sample2}"
    if [ `grep "${sample2}" ${sampleOutdir}/inputs.map.txt | grep "${fc}" | wc -l` -gt 1 ]; then
        echo "ERROR: Multiple R2 files found -- are there duplicate entries in ${sampleOutdir}/inputs.map.txt?"
        exit 3
    fi
    
    reads2fq=`grep "${sample2}" ${sampleOutdir}/inputs.map.txt | grep "${fc}"`
    if [ ! -f "${reads2fq}" ]; then
        echo "ERROR: Can not find R2 file ${reads2fq}"
        exit 4
    fi
    
    echo "Will process R2 reads file ${reads2fq}"
    
    PErun="TRUE"
    curfile=`echo ${sample1} | perl -pe 's/_R1(_\d+)?/_R1R2\1/g;'`
    curfile="${fc}${curfile}"
    
    
    echo "Filtering out reads with >75% G content"
    #TODO could potentially save R1 where only R2 is dark once it can handle single read. Could patch trimmomatic instead?
    #TODO super slow - replace with cutadapt --nextseq-trim?
    #BUGBUG only run on on PE data
    ${src}/filterNextSeqReadsForPolyG.py --inputfileR1 ${readsFq} --inputfileR2 ${reads2fq} --outputfileR1 $TMPDIR/${sample1}.pretrim.fastq.gz --outputfileR2 $TMPDIR/${sample2}.pretrim.fastq.gz --maxPolyG 75
    
    #seems to have fairly heavy memory requirements
    #java.lang.OutOfMemoryError: unable to create new native thread if run with just 2 threads
    java -XX:ParallelGCThreads=2 org.usadellab.trimmomatic.TrimmomaticPE ${trimmomaticBaseOpts} $TMPDIR/${sample1}.pretrim.fastq.gz $TMPDIR/${sample2}.pretrim.fastq.gz $TMPDIR/${sample1}.fastq.gz $TMPDIR/${sample1}.unpaired.fastq.gz $TMPDIR/${sample2}.fastq.gz $TMPDIR/${sample2}.unpaired.fastq.gz ${trimmomaticSteps}
    #BUGBUG java doesn't set nonzero exit code on trimmomatic exception
    
    rm -f $TMPDIR/${sample1}.pretrim.fastq.gz $TMPDIR/${sample2}.pretrim.fastq.gz
    
    #Merge anything unpaired from either R1 or R2
    zcat $TMPDIR/${sample1}.unpaired.fastq.gz $TMPDIR/${sample2}.unpaired.fastq.gz | gzip -c --fast > $TMPDIR/${curfile}.unpaired.fastq.gz
    rm -f $TMPDIR/${sample1}.unpaired.fastq.gz $TMPDIR/${sample2}.unpaired.fastq.gz
    unpairedReadLines=`zcat $TMPDIR/${curfile}.unpaired.fastq.gz | wc -l`
    echo "Unpaired read liness: ${unpairedReadLines}"
    
    
    if gzip -l $TMPDIR/${sample1}.fastq.gz | awk 'NR==2 {exit($2!=0)}' && gzip -l $TMPDIR/${sample1}.fastq.gz | awk 'NR==2 {exit($2!=0)}' && [[ "${unpairedReadLines}" == 0 ]]; then
        echo "WARNING: No tags passed filtering, quitting successfully"
        exit 0
    fi
    
    
    mkdir -p ${sampleOutdir}/fastqc
    fastQcOutdir="${sampleOutdir}/fastqc/${fc}${sample1}_qc"
    if [ ! -d "${fastQcOutdir}" ]; then
        #qsub -cwd -V -N ${sample1}.qc -o ${sampleOutdir}/fastqc/${sample1}.qc -S /bin/bash -j y -b y -p -500 "mkdir -p ${fastQcOutdir}; fastqc --outdir ${fastQcOutdir} $TMPDIR/${sample1}.fastq.gz"
        mkdir -p ${fastQcOutdir}; fastqc -t ${NSLOTS} --outdir ${fastQcOutdir} $TMPDIR/${sample1}.fastq.gz
    fi
    
    fastQcOutdir="${sampleOutdir}/fastqc/${fc}${sample2}_qc"
    if [ ! -d "${fastQcOutdir}" ]; then
        #qsub -cwd -V -N ${sample2}.qc -o ${sampleOutdir}/fastqc/${sample1}.qc -S /bin/bash -j y -b y -p -500 "mkdir -p ${fastQcOutdir}; fastqc --outdir ${fastQcOutdir} $TMPDIR/${sample2}.fastq.gz"
        mkdir -p ${fastQcOutdir}; fastqc -t ${NSLOTS} --outdir ${fastQcOutdir} $TMPDIR/${sample2}.fastq.gz
    fi
    
    echo
    echo "Histogram of read lengths for R1"
    zcat $TMPDIR/${sample1}.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' | sort -k1,1n
    
    echo
    echo "Histogram of read lengths for R2"
    zcat $TMPDIR/${sample2}.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' | sort -k1,1n
    
    echo
    echo "Histogram of read lengths for unpaired"
    zcat $TMPDIR/${curfile}.unpaired.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' | sort -k1,1n
else
    PErun="FALSE"
    curfile="${fc}${sample1}"
    
    #BUGBUG missing filterNextSeqReadsForPolyG.py for SE data
    #BUGBUG wrong adapter files
    java -XX:ParallelGCThreads=2 org.usadellab.trimmomatic.TrimmomaticSE ${trimmomaticBaseOpts} ${readsFq} $TMPDIR/${sample1}.fastq.gz ${trimmomaticSteps}
    #BUGBUG java doesn't set nonzero exit code on trimmomatic exception
    
    
    if gzip -l $TMPDIR/${sample1}.fastq.gz | awk 'NR==2 {exit($2!=0)}'; then
        echo "No tags passed filtering, quitting successfully"
        exit 0
    fi
    
    echo
    echo "Running fastqc"
    date
    mkdir -p ${sampleOutdir}/fastqc
    fastQcOutdir="${sampleOutdir}/fastqc/${fc}${sample1}_qc"
    if [ ! -d "${fastQcOutdir}" ]; then
        #qsub -cwd -V -N ${sample1}.qc -o ${sampleOutdir}/fastqc/${sample1}.qc -S /bin/bash -j y -b y -p -500 "mkdir -p ${fastQcOutdir}; fastqc --outdir ${fastQcOutdir} $TMPDIR/${sample1}.fastq.gz"
        mkdir -p ${fastQcOutdir}; fastqc -t ${NSLOTS} --outdir ${fastQcOutdir} $TMPDIR/${sample1}.fastq.gz
    fi
    
    echo
    date
    
    echo
    echo "Histogram of read lengths"
    zcat $TMPDIR/${sample1}.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' | sort -k1,1n
fi


echo
echo "Mapping ${readsFq} for ${curfile}"
echo "userAlnOptions=${userAlnOptions}"
echo "Will map to genomes ${genomesToMap}"
date

for curGenome in `echo ${genomesToMap} | perl -pe 's/,/ /g;'`; do
    echo
    echo "Mapping to reference ${curGenome}"
    date
    
    source ${src}/genomeinfo.sh ${curGenome}
    
    if [[ "${processingCommand}" == "mapBwaAln" ]]; then
        #not sure what -R is, making it lower than samse/pe -n reduces mapped PE tags but not SE tags
        #-Y filters sequences with \d+:Y:... after the space in the read name
        #Originally -n 0.04 seemed to allow upto two mismatches at 36 bp (must have rounded up)
        bwaAlnOpts="-n ${permittedMismatches} -l 32 ${userAlnOptions} -t ${NSLOTS} -Y"
        
        #Other options notes:
        #-q 0.20 does soft clip quality-based trimming of 3' end of reads, but only down to 35 bp
        #http://seqanswers.com/forums/showthread.php?t=5628
        #http://seqanswers.com/forums/showthread.php?t=6251
        
        echo "bwa aln ${bwaAlnOpts} ${bwaIndex} ..."
        bwa aln ${bwaAlnOpts} ${bwaIndex} $TMPDIR/${sample1}.fastq.gz > $TMPDIR/${sample1}.${curGenome}.sai
        
        
        #Previously used -n 10 but never really used XA tag and maybe was causing sampe to occasionally truncate last line of output (dropping the tags)
        bwaAlnExtractOpts="-n 3 -r ${readgroup}"
        if [[ "$PErun" == "TRUE" ]] ; then
            echo -e "\nMapping R2 ${reads2fq} for ${sample2}"
            date
            echo
            
            bwa aln ${bwaAlnOpts} ${bwaIndex} $TMPDIR/${sample2}.fastq.gz > $TMPDIR/${sample2}.${curGenome}.sai
            
            #-P didn't have a major effect, but some jobs were ~10-40% faster but takes ~16GB RAM instead of 4GB
            extractcmd="bwa sampe ${bwaAlnExtractOpts} -a ${maxInsertSize} ${bwaIndex} $TMPDIR/${sample1}.${curGenome}.sai $TMPDIR/${sample2}.${curGenome}.sai $TMPDIR/${sample1}.fastq.gz $TMPDIR/${sample2}.fastq.gz"
            
            #Only map unpaired reads if the file nonzero
            if [[ "${unpairedReadLines}" > 0 ]]; then
                echo -e "\nMapping unpaired ${curfile}.unpaired.fastq.gz for ${sample1}"
                date
                echo
                
                bwa aln ${bwaAlnOpts} ${bwaIndex} $TMPDIR/${curfile}.unpaired.fastq.gz > $TMPDIR/${curfile}.unpaired.${curGenome}.sai
                date
                
                echo
                echo "Extracting unpaired reads"
                unpairedReadsSam="$TMPDIR/${curfile}.${curGenome}.unpaired.sam"
                unpairedExtractcmd="bwa samse ${bwaAlnExtractOpts} ${bwaIndex} $TMPDIR/${curfile}.unpaired.${curGenome}.sai $TMPDIR/${curfile}.unpaired.fastq.gz"
                echo -e "unpairedExtractcmd=bwa $unpairedExtractcmd | (...)"
                bwa $unpairedExtractcmd | grep -v "^@" > ${unpairedReadsSam}
                
                extractcmd="${extractcmd} | cat - ${unpairedReadsSam}"
            fi
            #TODO merge headers instead of dropping
        else
            extractcmd="bwa samse ${bwaAlnExtractOpts} ${bwaIndex} $TMPDIR/${sample1}.${curGenome}.sai $TMPDIR/${sample1}.fastq.gz"
        fi
    elif [[ "${processingCommand}" == "mapBwaMem" ]]; then
        #parameters from ira hall consensus pipeline added 2020feb18 along with update to bwa 0.7.17        #https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md
        #I could not find any effect of '-Y'; in some spot checking I hever see bwa output anything with $6~/H/, though I can find some soft-clipped supplementary alignments with or without the option
        bwaMemOptions="-Y -K 100000000"
        
        extractcmd="bwa mem ${bwaMemOptions} ${userAlnOptions} -t ${NSLOTS} -R ${readgroup} ${bwaIndex} $TMPDIR/${sample1}.fastq.gz"
        if [[ "$PErun" == "TRUE" ]] ; then
            extractcmd="${extractcmd} $TMPDIR/${sample2}.fastq.gz"
        fi
        echo "bwa ${extractcmd}"
    elif [[ "${processingCommand}" == "mapMinimap" ]]; then
        #Assumes SE
        #note this ignores the trimmed reads in $TMPDIR/${sample1}.fastq.gz
        extractcmd="minimap2 -a -x map-ont ${referencefasta} ${readsFq}"
        echo "${extractcmd} ..."
    else
        echo "ERROR: impossible unsupported ${processingCommand}"
        exit 6
    fi
    
    date
    echo
    echo "Extracting"
    echo -e "extractcmd=${extractcmd} | (...)"
    date
    echo
    ${extractcmd} > $TMPDIR/${curfile}.${curGenome}.bwaout.bam
    
    
    echo
    echo "Post-processing"
    date
#    #calmd - this is glacially slow for some reason, not nearly as bad when run interactively
#    #Fix NM/MD (bwa aln still seems to set NM/MD wrong sporadically) and precalculate BAQ in the parallel thread to speed subsequent variant calling
#    #See http://www.biostars.org/p/1268
#    #NB redirecting stderr since calmd can be noisy, but you will miss real errors
#    samtools calmd -u -r - ${referencefasta} 2> $TMPDIR/${curfile}.${curGenome}.calmd.log |
    
    
    if [[ "${sampleType}" == "dnase" ]] || [[ "${sampleType}" == "atac" ]] || [[ "${sampleType}" == "chipseq" ]]; then
        unwanted_refs="--failUnwantedRefs --reqFullyAligned"
    else
        unwanted_refs=""
    fi
    
    if [[ "${curGenome}" =~ ^cegsvectors ]]; then
        dropUnmappedReads="--dropUnmappedReads"
        minMAPQ=0
    else
        dropUnmappedReads=""
        if [ "${readsLongEnough}" -eq 1 ]; then
            minMAPQ=20
        else
            minMAPQ=10
        fi
    fi
    
    #Not populating -SPECIES=Human
    #ParallelGCThreads is for java.lang.OutOfMemoryError: unable to create new native thread
    java -XX:ParallelGCThreads=1 -Dpicard.useLegacyParser=false -jar ${PICARDPATH}/picard.jar CreateSequenceDictionary -O=$TMPDIR/${curfile}.${curGenome}.dict -R=${referencefasta} -GENOME_ASSEMBLY=${annotationgenome}
    
    samtools sort -@ $NSLOTS -m 1750M -O bam -T $TMPDIR/${curfile}.sortbyname -l 1 -n $TMPDIR/${curfile}.${curGenome}.bwaout.bam |
    #not much gain other than avoiding disk IO doing this in a pipe as I don't think sort prints any intermediate results. Perhaps sorting fastq before mapping would be faster? https://www.biostars.org/p/15011/
    ${src}/filter_reads.py --reheader $TMPDIR/${curfile}.${curGenome}.dict ${unwanted_refs} ${dropUnmappedReads} --max_mismatches ${permittedMismatches} --min_mapq ${minMAPQ} --max_insert_size ${maxInsertSize} - ${sampleOutdir}/${curfile}.${curGenome}.bam
    
    echo
    date
    
    
    if [[ "${sampleType}" == "amplicon" ]]; then
        echo
        echo "Primer soft clipping"
        samtools view -h ${sampleOutdir}/${curfile}.${curGenome}.bam |
        #TODO hardcoded primer coords for now
        #BUGBUG primerclip leaves no PG record; otherwise would set samtools view --no-PG
        #BUGBUG dumps masterparsefails.log into directory
        #BUGBUG can't stream to stdout; passes log output through stdout?
        primerclip /vol/sars/sequences/wuhCor1/Swift_Amplicons/sarscov2_masterfile.txt /dev/stdin  $TMPDIR/${curfile}.${curGenome}.sam
        samtools view -@ $NSLOTS -O bam -1 $TMPDIR/${curfile}.${curGenome}.sam > ${sampleOutdir}/${curfile}.${curGenome}.bam
        
        echo
        date
    fi
    
    
    #Add MC tag containing mate CIGAR for duplicate calling
    #Why do I need this? samblaster can add this itself but seems to miss some?
    #Needs to be sorted by coordinate
    #java -Xmx2g -Dpicard.useLegacyParser=false -jar ${PICARDPATH}/picard.jar FixMateInformation -INPUT ${TMPDIR}/${curfile}.${curGenome}.bwaout.bam -OUTPUT ${sampleOutdir}/${curfile}.${curGenome}.bam -VERBOSITY ERROR -QUIET TRUE -COMPRESSION_LEVEL 1
    
#    echo
#    echo "Cleanup"
#    date
#    cp ${curfile}.${curGenome}.bam $TMPDIR/${curfile}.${curGenome}.unclean.bam
#    #BUGBUG Should fix the ERROR... read errors, but doesn't do anything to first 100000 lines of test case except increment version to 1.4 "@HD   VN:1.4"
#    java -Xmx2g -jar ${PICARDPATH}/picard.jar/ CleanSam INPUT=${curfile}.${curGenome}.bam OUTPUT=${curfile}.${curGenome}.clean.bam COMPRESSION_LEVEL=1 && mv ${curfile}.${curGenome}.clean.bam ${curfile}.${curGenome}.bam
    
    
    echo
    echo "SAMtools statistics for genome ${curGenome}"
    date
    samtools flagstat ${sampleOutdir}/${curfile}.${curGenome}.bam
    
    
    echo
    echo "QC metrics to be done on un-merged data"
    date
    
    echo
    echo "Mean quality by cycle"
    #BUGBUG performs badly for SRR jobs -- some assumption not met?
    #BUGBUG reports "WARNING   2019-12-17 09:08:12     SinglePassSamProgram    File reports sort order 'queryname', assuming it's coordinate sorted anyway." mainly on cegsvectors? baseq still seems to be output. Not sure why it cares about sort order.
    java -XX:ParallelGCThreads=2 -Xmx3g -Dpicard.useLegacyParser=false -jar ${PICARDPATH}/picard.jar MeanQualityByCycle -INPUT ${sampleOutdir}/${curfile}.${curGenome}.bam -OUTPUT $TMPDIR/${curfile}.baseQ.txt -CHART_OUTPUT $TMPDIR/${curfile}.baseQ.pdf -VALIDATION_STRINGENCY LENIENT
    
    instrument=`bam2instrument ${sampleOutdir}/${curfile}.${curGenome}.bam | uniq | awk 'values[$0] != 1 {print; values[$0]=1}' | perl -pe 's/\n$//g;' | perl -pe 's/\n/;/g;'`
    awk -v instrument=${instrument} -v fc=${fc} -v sample=${curfile} -v bs=${BS} -v genome=${curGenome} -F "\t" 'BEGIN {OFS="\t"} $0!~/^#/ && $0!="" {if($1=="CYCLE") {$0=tolower($0); $(NF+1)="instrument\tfc\tsample\tBS\tgenome"} else {$(NF+1)=instrument "\t" fc "\t" sample "\t" bs "\t" genome;} print}' $TMPDIR/${curfile}.baseQ.txt > ${sampleOutdir}/${curfile}.${curGenome}.baseQ.txt
    
    
    echo
    #Prints positions
    echo -n -e "Histogram of mismatches to reference by position:\t${instrument}\t${fc}\t${curfile}\t${BS}\t${curGenome}\t"
    samflags="-q 20 -F 1548"
    
    #BUGBUG flag 2 not working?
    #TODO slow
    samtools view ${samflags} ${sampleOutdir}/${curfile}.${curGenome}.bam | awk -F "\t" 'BEGIN {OFS="\t"} { \
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
    }' | tee $TMPDIR/${curfile}.mm.txt |
    awk -v minBaseQ=20 -F "\t" 'BEGIN { for (i=0; i<256; i++) { codeFor[sprintf("%c", i)] = i } } codeFor[$4]-33 > minBaseQ {print} ' | 
    cut -f2 | sort -g | uniq -c | sort -k2,2 -g | awk 'BEGIN {ORS="\t"} {print $1}'
    echo
    
    
#    echo
#    echo
#    echo "Gerald's call for mismatched positions"
#    date
#    awk 'BEGIN {OFS="\t"; print "cycle", "A", "C", "G", "T", "N"} {errors[$2,$3]++} END {for(i=1; i<=36; i++) {print i, errors[i, "A"], errors[i, "C"], errors[i, "G"], errors[i, "T"], errors[i, "N"]}}' $TMPDIR/${curfile}.mm.txt
done


echo
echo "Done!"
date
