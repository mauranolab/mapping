#!/bin/bash
set -eu -o pipefail

module add miller
module add pigz

#newer version of gnu coreutils would permit single-pass split --bytes=1024M --filter='gzip > $FILE.gz' /path/to/input /path/to/output

fc=$1

echo "Running on $HOSTNAME"
echo "Postprocessing ${fc}"


echo
echo "Total read counts"
date
rm -f readcounts.txt; touch readcounts.txt
echo -e "Project\tSample\t#Reads" >> readcounts.txt
for dir in `find .  -maxdepth 1 -name "Project_*" -type d | sort -r | perl -pe 's/^\.\///g;'` .; do
    bs=`find $dir -maxdepth 2 -type f -and -name "*.fastq.gz" -and -printf "%f\n" | cut -d "_" -f1 | sort | uniq`
    for sample1 in ${bs}; do
        echo -e -n "$dir\t$sample1\t" >> readcounts.txt
        files=`find $dir -maxdepth 2 -name "$sample1*.gz"`
        #This is limited by decompression speed
        #Some profiling info in https://bioinformatics.stackexchange.com/questions/935/fast-way-to-count-number-of-reads-and-number-of-bases-in-a-fastq-file
        #Requires sequences/qualities to be on single line
        printf "%'d\n" `pigz -dc $files | awk -F "\t" 'BEGIN {OFS="\t"} END {print NR/4}'` >> readcounts.txt
    done
done


#BUGBUG no commas in total reads
cat readcounts.txt | perl -pe 's/,//g;' | awk -F "\t" 'BEGIN {OFS="\t"} NR>1 {sum+=$3} END {print ""; print "Total_reads", ".", sum}' >> readcounts.txt
cat readcounts.txt


echo
echo "Sorted by BS number"
cat readcounts.txt | awk -F "\t" 'BEGIN {OFS="\t"} $0=="" {exit} $2!="Undetermined" {print}' | mlr --tsv sort -f Sample then cut -f "#Reads"


echo
date
echo "split fastq into 16M read chunks"
splitreads=16000000
for base in `cat readcounts.txt | perl -pe 's/,//g;' | awk -F "\t" 'BEGIN {OFS="\t"} $0=="" {exit} {print}' | awk -v minsizetosplit=0 -F "\t" 'BEGIN {OFS="\t"} NR>1 && $1!="." && $3>minsizetosplit*2' | awk '$1~/Project_CEGS/ || $1~/Project_Maurano/ || $1~/Project_SARS/' | awk -F "\t" 'BEGIN {OFS="\t"} {print $1 ";" $2 ";" $3}'`; do
    project=`echo ${base} | cut -d ";" -f1`
    bs=`echo ${base} | cut -d ";" -f2`
    numreads=`echo ${base} | cut -d ";" -f3`
    
    if [[ "${numreads}" -lt 150000000 ]]; then
        #slightly lower priority so that big jobs start first if possible
        qsubargs="-pe threads 8 -p -1"
    else
        #There can't be too many samples this large per FC, so safe to use a lot of threads
        qsubargs="-pe threads 24"
    fi
    
    echo
    echo -e "${project}\t${bs}"
    #bak.fastq shouldn't exist at this point, but exclude it just in case
    for f1 in `find ${project}/Sample_${bs} -name "*_R1_*.fastq.gz" | grep -v bak.fastq`; do
        f2=`echo $f1 | perl -pe 's/_R1_/_R2_/g;'`
        echo "$f1 $f2"
        #add ${project}/Sample_${bs}/bak.fastq to make backup
        qsub -S /bin/bash ${qsubargs} --time 4:00:00 --mem-per-cpu 2G -j y -N split.${bs}.R1 "/gpfs/data/mauranolab/mapped/src/flowcells/splitfastq.sh ${bs} ${f1} R1 ${splitreads}"
        if [ -s "$f2" ]; then
            qsub -S /bin/bash --time 4:00:00 --mem-per-cpu 2G -j y -N split.${bs}.R2 ${qsubargs} "/gpfs/data/mauranolab/mapped/src/flowcells/splitfastq.sh ${bs} ${f2} R2 ${splitreads}"
        fi
    done
done

echo
echo "Tarballing collaborator projects"
date
basedir=`pwd`
#Rarely in a hurry for this so might as well do in serial to simpify
for d in `find ${basedir}  -maxdepth 1 -name "Project_*" -type d | grep -v Project_CEGS | grep -v Project_Maurano | grep -v Project_SARS`; do
    project=`basename $d`
    echo
    echo "${project}"
    base="${project}_${fc}_"`date +%Y%b%d`
    cd $d
    tar cvf /gpfs/data/isg_sequencing/public_html/${base}_fastq.tar .
    
    echo
    echo "Your data are available at (user:sequencing, password:moretags)"
    echo "https://sequencing:moretags@cascade.isg.med.nyu.edu/flowcells/${base}_fastq.tar"
    echo
    echo "The link to your sequencing data is only guaranteed to work for 30 days,"
    echo "so be sure to download it and store it yourself ASAP."
    
    echo
    awk -v project=${project} -F "\t" 'BEGIN {OFS="\t"} NR==1 || $1==project' ${basedir}/readcounts.txt
done

#echo "Tarballing flowcell data directory"
#base="${fc}_"`date +%Y%b%d`
#cd /gpfs/data/isg_sequencing/data/$fc/
#echo "Your data are available at (user:sequencing, password:moretags)"
#echo "https://sequencing:moretags@cascade.isg.med.nyu.edu/flowcells/${base}_data.tar"
#tar cvf /gpfs/data/isg_sequencing/public_html/${base}_data.tar .


echo
echo "Done!!!"
date
