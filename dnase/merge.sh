#!/bin/bash
set -e -o pipefail

name=$1
DS=$2
genome=$3

echo "output name $name"
date

sampleOutdir=$name

#Here we generate a list of expected filenames
files=""
#NB originally checked grep R1 but some old solexa lanes didn't include it
for f1 in `grep $DS inputs.txt | grep -v R2`; do
       #echo "adding $f1"
       
       #First get FC
       #NB This must match exactly what is in runAnalysis.sh
       #BUGBUG a bit fragile
       #BUGBUG misses things like UwStam_CH12-DS22536-FCD0TGK-L002_R1_001.fastq.gz and UwStam_CH12-DS22542-FCD0TDB-L001_R1_002.fastq.gz, but don't see any collision
#       if [[ `basename $f1` =~ "^s_" ]]; then
#              fc=`readlink -f $f1 | perl -pe 's/\/\d\d\d\/s_/\/s_/g;' | xargs dirname | xargs basename | perl -pe 's/_\d+_tag//g;'`
#       else
#              fc=`readlink -f $f1 | xargs dirname | xargs dirname | xargs dirname | xargs basename | perl -pe 's/_\d+_tag//g;'`
#       fi
#       if [[ "$fc" == "." ]] ; then
#              fc=""
#       else
##              echo "Flowcell $fc"
#              fc="${fc}."
#       fi
       fc=""
       
       f2=`echo $f1 | perl -pe 's/_R1(_\d+)?/_R2$1/g;'`
       if echo $f2 | grep -q R2 && grep -q $f2 inputs.txt ; then
              f1=`echo $f1 | perl -pe 's/_R1(_\d+)?/_R1R2\1/g;'`
       fi
       
       curOutputFile=`basename $f1 .fastq.gz`
       curOutputFile="$sampleOutdir/${fc}${curOutputFile}.$genome.bam"
       
       #echo "cur $files"
       
       if [[ -f "$curOutputFile" ]]; then
              files="$files ${curOutputFile}"
              numfiles=$((numfiles+1))
       else
              echo "Error: $curOutputFile doesn't exist!"
              #Don't die to work around weirdness in pipeline FQ files
              #exit 1
       fi
done


if [[ "$numfiles" -eq 0 ]]; then
       echo "Error: No files found to merge!"
       exit 2
fi

echo "Will merge $files"


if [[ "$numfiles" -eq 1 ]]; then
       echo "copying file"
       cp $files $sampleOutdir/$name.bam
else
       echo "merging files"
       samtools merge -l 1 -@ $NSLOTS $sampleOutdir/$name.bam $files
fi

echo "Processing $sampleOutdir/$name.bam"
#TODO AddOrReplaceReadGroups or do bwa -r ""


echo
echo "mark dups"
date
###Obsolete picard version
#Used to need VALIDATION_STRINGENCY=LENIENT to avoid SAM validation error: ERROR...MAPQ should be 0 for unmapped read or CIGAR should have zero elements for unmapped read
#http://seqanswers.com/forums/showthread.php?t=4246
#BUGBUG Can make huge log files despite these options
#Now getting new error: "Ignoring SAM validation error: ERROR: Record 11792763, Read name HISEQ-2500-1:94:C74YHANXX:2:1109:6562:81875, bin field of BAM record does not equal value computed based on alignment start and end, and length of sequence to which read is aligned"
#TODO switch to MarkDuplicatesWithMateCigar for performance (identical output)
#http://sourceforge.net/p/samtools/mailman/message/32910359/
#java -Xmx6g -jar /home/maurano/bin/picard-tools/MarkDuplicates.jar INPUT=$name.bam OUTPUT=$name.markedDups.bam METRICS_FILE=$TMPDIR/$name.picardDups.txt QUIET=TRUE VERBOSITY=ERROR COMPRESSION_LEVEL=9 ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT && mv $name.markedDups.bam $name.bam

###Samblaster is faster
#samblaster used an average of 1GB memory mapping ENCODE DNase data to hg38. 10/889 jobs used >5GB.
samtools sort -@ $NSLOTS -O sam -m 1250M -T $TMPDIR/${name}.sortbyname -l 1 -n $sampleOutdir/$name.bam |
samblaster |
samtools view -Sb - |
samtools sort -@ $NSLOTS -O bam -m 1250M -T $TMPDIR/${name}.sort -l 9 - > $sampleOutdir/$name.markedDups.bam && mv $sampleOutdir/$name.markedDups.bam $sampleOutdir/$name.bam


echo
echo "Indexing"
date
samtools index $sampleOutdir/$name.bam


rm -f $files


echo
echo -e "\nDone!"
date
