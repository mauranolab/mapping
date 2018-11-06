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
alias bedmap='bedmap --ec --header'
alias starch='starch --header'
alias closest-features='closest-features --header'


genomesToMap=$1
analysisType=$2
sampleOutdir=$3
DS=$4
src=$5

shift
shift
shift
shift
shift
userAlnOptions=$@


echo "Output directory: ${sampleOutdir}, DS: ${DS}, Genomes: ${genomesToMap}, analysisType:${analysisType}"
date


processingCommand=`echo "${analysisType}" | awk -F "," '{print $1}'`
analysisCommand=`echo "${analysisType}" | awk -F "," '{print $2}'`


permittedMismatches=3
if [[ "${analysisCommand}" == "dnase" ]] || [[ "${analysisCommand}" == "atac" ]] || [[ "${analysisCommand}" == "chipseq" ]]; then
    maxInsertSize=500
else
    maxInsertSize=1000
fi


bam2instrument()
{
    samtools view -F 2304 $* |
    cut -f1 | cut -d ":" -f1 | 
    #Hack to clean up some encodeproject.org data that has underscore in place of colon after sequencer name
    perl -pe 's/_\d+_\d+_\d+_\d+$//g;' | 
    #SRR data
    perl -pe 's/^(SRR\d+)\.\d+$/$1/g;'
}


jobid=$SGE_TASK_ID
readsFq=`cat ${sampleOutdir}/inputs.map.txt | awk -v jobid=$jobid 'NR==jobid'`
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
echo "using $TMPDIR as TMPDIR"

mkdir -p ${sampleOutdir}


echo
echo "Configuring trimming parameters"
#Trimmomatic options
if [[ "${analysisCommand}" == "atac" ]]; then 
    illuminaAdapters="/cm/shared/apps/trimmomatic/0.38/adapters/NexteraPE-PE.fa"
else 
    #TODO Probably need different sequences per barcode. Note this fa file has 2 ident copies of left adapter and none of right adapter (with barcode).
    illuminaAdapters="/cm/shared/apps/trimmomatic/0.38/adapters/TruSeq3-PE-2.fa"
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
#sequencedTags=$(zcat ${readsFq} | awk 'NR%4==2' | wc -l)
#BUGBUG - I think Jesper just took what looks like a readthrough sequence (TCGTATGCCGTCTTC). Not sure why the trimmer isn't dealing with this properly
#if [ `zcat ${readsFq} | awk -v thresh=0.25 -v sequencedTags=$sequencedTags 'NR%4==2 && $1~/TCGTATGCCGTCTTC/ {readsWithDukeSequence+=1} END {if (readsWithDukeSequence/sequencedTags>thresh) {print 1} else {print 0}}'` ]; then
#    echo "More than 25% of reads have DUKE sequence (TCGTATGCCGTCTTC) - Hard clip to 20bp reads"
#    trimmomaticSteps="CROP:20 ${trimmomaticSteps}"
#    minMAPQ=10
#else
#    echo "No DUKE sequence present"


minMAPQ=20
trimmomaticSteps="${trimmomaticSteps} MINLEN:27"
if [[ "${analysisCommand}" == "dnase" ]] || [[ "${analysisCommand}" == "atac" ]] || [[ "${analysisCommand}" == "chipseq" ]]; then
    trimmomaticSteps="${trimmomaticSteps} CROP:36"
fi


#BUGBUG a bit fragile
fc=`readlink -f ${readsFq} | xargs dirname | xargs dirname | xargs dirname | xargs basename`
if [[ ! "${fc}" =~ ^FC ]] ; then
    fc=""
else
    echo "Flowcell ${fc}"
    fc="${fc}."
fi


DS_nosuffix=`echo ${DS} | perl -pe 's/[A-Z]$//g;'`
readgroup="@RG\\tID:${fc}${DS}\\tLB:${DS}\\tSM:${DS_nosuffix}\\tPL:ILLUMINA"
if [[ "${readsFq}" =~ ^\/vol\/mauranolab\/flowcells\/fastq\/ ]]; then
    readgroup_instrument=`awk -F "\t" 'BEGIN {OFS="\t"} $1=="#Instrument" {print $2}' /vol/mauranolab/flowcells/data/${fc/./}/info.txt`
    
    readgroup_date=`awk -F "\t" 'BEGIN {OFS="\t"} $1=="#Load date" {print $2}' /vol/mauranolab/flowcells/data/${fc/./}/info.txt`
    readgroup_bcs=`awk -v ds=${DS} -F "\t" 'BEGIN {OFS="\t"} $0!~/^#/ && 0!="" && $2 $3==ds {split($7, bc1, "_"); split($7, bc2, "_"); print bc1[2] "-" bc1[2]}' /vol/mauranolab/flowcells/data/${fc/./}/info.txt`
    #BUGBUG BC: shows up in bwa command line but at some point disappears from the bam header
    readgroup="${readgroup}\\tDT:${readgroup_date}\\tBC:${readgroup_bcs}\\tPU:${fc/./}-${readgroup_bcs}"
    
    case "${readgroup_instrument}" in
    Balin)
        readgroup="${readgroup}\\tCN:Maurano_Lab\\tPM:NextSeq_500"
        ;;
    GTC_NovaSeq)
        readgroup="${readgroup}\\tCN:NYUMC_GTC\\tPM:NovaSeq_6000"
        ;;
    esac
fi


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
    ${src}/filterNextSeqReadsForPolyG.py --inputfileR1 ${readsFq} --inputfileR2 ${reads2fq} --outputfileR1 $TMPDIR/${sample1}.pretrim.fastq.gz --outputfileR2 $TMPDIR/${sample2}.pretrim.fastq.gz --maxPolyG 75
    
    java org.usadellab.trimmomatic.TrimmomaticPE ${trimmomaticBaseOpts} $TMPDIR/${sample1}.pretrim.fastq.gz $TMPDIR/${sample2}.pretrim.fastq.gz $TMPDIR/${sample1}.fastq $TMPDIR/${sample1}.unpaired.fastq $TMPDIR/${sample2}.fastq $TMPDIR/${sample2}.unpaired.fastq ${trimmomaticSteps}
    #TODO why does this output uncompressed fastq? I think it's just so one can test with -s below
    
    echo -n "Unpaired reads:"
    #Merge anything unpaired from either R1 or R2
    cat $TMPDIR/${sample1}.unpaired.fastq $TMPDIR/${sample2}.unpaired.fastq | tee $TMPDIR/${curfile}.unpaired.fastq | wc -l
    
    
    if [ ! -s "$TMPDIR/${sample1}.fastq" ] && [ ! -s "$TMPDIR/${sample2}.fastq" ]; then
        echo "No tags passed filtering, quitting successfully"
        exit 0
    fi
    
    
    mkdir -p fastqc
    fastQcOutdir="fastqc/${fc}${sample1}_qc"
    if [ ! -d "${fastQcOutdir}" ]; then
        #qsub -cwd -V -N ${sample1}.qc -o fastqc/${sample1}.qc -S /bin/bash -j y -b y -p -500 "mkdir -p ${fastQcOutdir}; fastqc --outdir ${fastQcOutdir} $TMPDIR/${sample1}.fastq"
        mkdir -p ${fastQcOutdir}; fastqc --outdir ${fastQcOutdir} $TMPDIR/${sample1}.fastq
    fi
    
    fastQcOutdir="fastqc/${fc}${sample2}_qc"
    if [ ! -d "${fastQcOutdir}" ]; then
        #qsub -cwd -V -N ${sample2}.qc -o fastqc/${sample1}.qc -S /bin/bash -j y -b y -p -500 "mkdir -p ${fastQcOutdir}; fastqc --outdir ${fastQcOutdir} $TMPDIR/${sample2}.fastq"
        mkdir -p ${fastQcOutdir}; fastqc --outdir ${fastQcOutdir} $TMPDIR/${sample2}.fastq
    fi
    
    echo
    echo "Histogram of read lengths for R1"
    cat $TMPDIR/${sample1}.fastq | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' | sort -k1,1n
    
    echo
    echo "Histogram of read lengths for R2"
    cat $TMPDIR/${sample2}.fastq | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' | sort -k1,1n
    
    echo
    echo "Histogram of read lengths for unpaired"
    cat $TMPDIR/${curfile}.unpaired.fastq | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' | sort -k1,1n
else
    PErun="FALSE"
    curfile="${fc}${sample1}"
    
    #TODO missing filterNextSeqReadsForPolyG.py
    
    #BUGBUG wrong adapter files
    java org.usadellab.trimmomatic.TrimmomaticSE ${trimmomaticBaseOpts} ${readsFq} $TMPDIR/${sample1}.fastq ${trimmomaticSteps}
    
    
    if [ ! -s "$TMPDIR/${sample1}.fastq" ]; then
        echo "No tags passed filtering, quitting successfully"
        exit 0
    fi
    
    
    mkdir -p fastqc
    fastQcOutdir="fastqc/${fc}${sample1}_qc"
    if [ ! -d "${fastQcOutdir}" ]; then
        #qsub -cwd -V -N ${sample1}.qc -o fastqc/${sample1}.qc -S /bin/bash -j y -b y -p -500 "mkdir -p ${fastQcOutdir}; fastqc --outdir ${fastQcOutdir} $TMPDIR/${sample1}.fastq"
        mkdir -p ${fastQcOutdir}; fastqc --outdir ${fastQcOutdir} $TMPDIR/${sample1}.fastq
    fi
    
    echo
    echo "Histogram of read lengths"
    cat $TMPDIR/${sample1}.fastq | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' | sort -k1,1n
fi


echo
echo "Mapping ${readsFq} for ${curfile}"
echo "userAlnOptions=${userAlnOptions}"
echo "Will map to genomes ${genomesToMap}"
date


#On analysis sets
#https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
#https://gatkforums.broadinstitute.org/gatk/discussion/7857/reference-genome-components

#- Neither ENCODE hg38_no_alt_analysis_set nor hg38 have PAR masked
#PAR in hg38 from wikipedia:
#chrY	10001	2781479	PAR1
#chrY	56887903	57217415	PAR2
#mm10 https://www.ncbi.nlm.nih.gov/grc/mouse
#chrY	90745845 	91644698	PAR
#Rat is unclear which end of chrX it's on
#
#- There are no non-ACGTN in mm10, hg19/38 (from UCSC goldenpath, ENCODE analysis set, or the UCSC-provided NCBI analysis set) or rn6 (UCSC), but NCBI hg38 analysis set has 94
#zcat /vol/isg/annotation/fasta/mm10/mm10.fa.gz | grep -v "^>" | perl -pe 's/[ACGTN\n]+//ig;' | wc -c
#
#UCSC source: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/
#NCBI source (has non-N ambiguous bases and decoys):ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
#ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCA_000001405.27_GRCh38.p12
#
#- Decided not to use decoy hs38d1 -- couldn't find info on what exactly is in here but I think has EBV, HHV, etc. and other unmapping bits empirically found in HuRef and GM12878 by Heng Li and/or Broad

#Notes
#- No non-encode analysis set for mm10 seems available
#- None seem to use latest hg38 (p12 as of 10/15/18), can't find good inventory of changes

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
        bwa aln ${bwaAlnOpts} ${bwaIndex} $TMPDIR/${sample1}.fastq > $TMPDIR/${sample1}.${curGenome}.sai
        
        
        #Previously used -n 10 but never really used XA tag and maybe was causing sampe to occasionally truncate last line of output (dropping the tags)
        #TODO perhaps I should add -N to get complete set of hits?
        bwaAlnExtractOpts="-n 3 -r ${readgroup}"
        if [[ "$PErun" == "TRUE" ]] ; then
            echo -e "\nMapping R2 ${reads2fq} for ${sample2}"
            date
            echo
            
            bwa aln ${bwaAlnOpts} ${bwaIndex} $TMPDIR/${sample2}.fastq > $TMPDIR/${sample2}.${curGenome}.sai
            
            #-P didn't have a major effect, but some jobs were ~10-40% faster but takes ~16GB RAM instead of 4GB
            #TODO permit higher insert size for non-DNase?
            extractcmd="sampe ${bwaAlnExtractOpts} -a ${maxInsertSize} ${bwaIndex} $TMPDIR/${sample1}.${curGenome}.sai $TMPDIR/${sample2}.${curGenome}.sai $TMPDIR/${sample1}.fastq $TMPDIR/${sample2}.fastq"
            
            #Only map unpaired reads if the file nonzero
            if [ -s "$TMPDIR/${curfile}.unpaired.fastq" ]; then
                echo -e "\nMapping unpaired ${curfile}.unpaired.fastq for ${sample1}"
                date
                echo
                
                bwa aln ${bwaAlnOpts} ${bwaIndex} $TMPDIR/${curfile}.unpaired.fastq > $TMPDIR/${curfile}.unpaired.${curGenome}.sai
                date
                
                echo
                echo "Extracting unpaired reads"
                unpairedReadsSam="$TMPDIR/${curfile}.${curGenome}.unpaired.sam"
                unpairedExtractcmd="samse ${bwaAlnExtractOpts} ${bwaIndex} $TMPDIR/${curfile}.unpaired.${curGenome}.sai $TMPDIR/${curfile}.unpaired.fastq"
                echo -e "unpairedExtractcmd=bwa $unpairedExtractcmd | (...)"
                bwa $unpairedExtractcmd | grep -v "^@" > ${unpairedReadsSam}
                
                extractcmd="${extractcmd} | cat - ${unpairedReadsSam}"
            fi
            #TODO merge headers instead of dropping
        else
            extractcmd="samse ${bwaAlnExtractOpts} ${bwaIndex} $TMPDIR/${sample1}.${curGenome}.sai $TMPDIR/${sample1}.fastq"
        fi
    elif [[ "${processingCommand}" == "mapBwaMem" ]]; then
        #       -V            output the reference FASTA header in the XR tag
        extractcmd="mem ${userAlnOptions} -t ${NSLOTS} -R ${readgroup} ${bwaIndex} $TMPDIR/${sample1}.fastq"
        if [[ "$PErun" == "TRUE" ]] ; then
            extractcmd="${extractcmd} $TMPDIR/${sample2}.fastq"
        fi
        echo "bwa ${extractcmd}"
    else
        echo "ERROR: impossible unsupported ${processingCommand}"
        exit 6
    fi
    
    date
    echo
    echo "Extracting"
    echo -e "extractcmd=bwa ${extractcmd} | (...)"
    date
    echo
    
    bwa ${extractcmd} > $TMPDIR/${curfile}.${curGenome}.bwaout.bam
    
    #Piping straight to sort seems to increase frequency of spurious errors (stale file handles, etc.)
    
#    #calmd - this is glacially slow for some reason, not nearly as bad when run interactively
#    #Fix NM/MD (bwa aln still seems to set NM/MD wrong sporadically) and precalculate BAQ in the parallel thread to speed subsequent variant calling
#    #See http://www.biostars.org/p/1268
#    #NB redirecting stderr since calmd can be noisy, but you will miss real errors
#    samtools calmd -u -r - ${referencefasta} 2> $TMPDIR/${curfile}.${curGenome}.calmd.log |
    
    samtools sort -@ $NSLOTS -m 1750M -O bam -T $TMPDIR/${curfile}.sortbyname -l 1 -n $TMPDIR/${curfile}.${curGenome}.bwaout.bam > ${sampleOutdir}/${curfile}.${curGenome}.bam
    
    echo
    echo "Post-processing"
    date
    
    
    if [[ "${processingCommand}" == "mapBwaAln" ]]; then
        if [[ "${analysisCommand}" == "dnase" ]] || [[ "${analysisCommand}" == "atac" ]] || [[ "${analysisCommand}" == "chipseq" ]]; then
            unwanted_refs="--failUnwantedRefs"
        else
            unwanted_refs=""
        fi
        
        
        if [[ "${curGenome}" == "cegsvectors" ]]; then
            dropUnmappedReads="--dropUnmappedReads"
        else
            dropUnmappedReads=""
        fi
        
        #BUGBUG filter_reads.py not working with bwa mem supplementary alignments
        #TODO not really any point in piping this as I don't think sort prints any intermediate results. Perhaps sorting fastq before mapping would be faster? https://www.biostars.org/p/15011/
        #TODO should we really be hard-unmapping unscaffolded contigs, etc.?
        ${src}/filter_reads.py --reqFullyAligned ${unwanted_refs} ${dropUnmappedReads} --max_mismatches ${permittedMismatches} --min_mapq ${minMAPQ} --max_insert_size ${maxInsertSize}  ${sampleOutdir}/${curfile}.${curGenome}.bam ${sampleOutdir}/${curfile}.${curGenome}.new.bam && mv ${sampleOutdir}/${curfile}.${curGenome}.new.bam ${sampleOutdir}/${curfile}.${curGenome}.bam
        
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
    #TODO break out as separate job?
    
    
    echo
    echo "Mean quality by cycle"
    date
    #BUGBUG performs badly for SRR jobs -- some assumption not met?
    java -Xmx3g -Dpicard.useLegacyParser=false -jar ${PICARDPATH}/picard.jar MeanQualityByCycle -INPUT ${sampleOutdir}/${curfile}.${curGenome}.bam -OUTPUT $TMPDIR/${curfile}.baseQ.txt -CHART_OUTPUT $TMPDIR/${curfile}.baseQ.pdf -VALIDATION_STRINGENCY LENIENT
    
    instrument=`bam2instrument ${sampleOutdir}/${curfile}.${curGenome}.bam | uniq | awk 'values[$0] != 1 {print; values[$0]=1}' | perl -pe 's/\n$//g;' | perl -pe 's/\n/;/g;'`
    awk -v instrument=${instrument} -v fc=${fc} -v sample=${curfile} -v ds=${DS} -v genome=${curGenome} -F "\t" 'BEGIN {OFS="\t"} $0!~/^#/ && $0!="" {if($1=="CYCLE") {$0=tolower($0); $(NF+1)="instrument\tfc\tsample\tDS\tgenome"} else {$(NF+1)=instrument "\t" fc "\t" sample "\t" ds "\t" genome;} print}' $TMPDIR/${curfile}.baseQ.txt > ${sampleOutdir}/${curfile}.${curGenome}.baseQ.txt
    
    
    echo
    #Prints positions
    echo -n -e "Histogram of mismatches to reference by position:\t${instrument}\t${fc}\t${curfile}\t${DS}\t${curGenome}\t"
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
    
    
    echo
    echo
    echo "Gerald's call for mismatched positions"
    date
    awk 'BEGIN {OFS="\t"; print "cycle", "A", "C", "G", "T", "N"} {errors[$2,$3]++} END {for(i=1; i<=36; i++) {print i, errors[i, "A"], errors[i, "C"], errors[i, "G"], errors[i, "T"], errors[i, "N"]}}' $TMPDIR/${curfile}.mm.txt
done


echo
echo "Done!"
date
