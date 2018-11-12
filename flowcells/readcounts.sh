#!/bin/bash
set -eu -o pipefail

module add miller

fc=$1


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
        printf "%'d\n" `zcat $files | awk -F "\t" 'BEGIN {OFS="\t"} END {print NR/4}'` >> readcounts.txt
    done
done


#BUGBUG no commas in total reads
cat readcounts.txt | perl -pe 's/,//g;' | awk -F "\t" 'BEGIN {OFS="\t"} NR>1 {sum+=$3} END {print ""; print "Total_reads", ".", sum}' >> readcounts.txt
cat readcounts.txt


echo
echo "Sorted by BS number"
cat readcounts.txt | awk -F "\t" 'BEGIN {OFS="\t"} $0=="" {exit} {print}' | mlr --tsv sort -f Sample


echo
date
echo "split fastq into 16M read chunks"
#TODO parallelize
splitreads=16000000
#bak.fastq shouldn't exist at this point, but exclude it just in case
for base in `cat readcounts.txt | perl -pe 's/,//g;' | awk -F "\t" 'BEGIN {OFS="\t"} $0=="" {exit} {print}' | awk -v splitreads=${splitreads} -F "\t" 'BEGIN {OFS="\t"} NR>1 && $1!="." && $3>splitreads*2' | awk '$1~/Project_CEGS/ || $1~/Project_Maurano/' | awk -F "\t" 'BEGIN {OFS="\t"} {print $1 ";" $2}'`; do
    project=`echo ${base} | cut -d ";" -f1`
    bs=`echo ${base} | cut -d ";" -f2`
    echo
    echo -e "${project}\t${bs}"
    for f1 in `find ${project}/Sample_${bs} -name "*_R1_*.fastq.gz" | grep -v bak.fastq`; do
        f2=`echo $f1 | perl -pe 's/_R1_/_R2_/g;'`
        echo "$f1 $f2"
        #NB assumes PE reads
        /vol/mauranolab/flowcells/src/splitfastq.sh ${bs} ${f1} ${f2} ${splitreads} ${project}/Sample_${bs}/bak.fastq
    done
done
##TODO maybe rename files <16M reads for consistency

echo
echo "Tarballing collaborator projects"
date
basedir=`pwd`
#Rarely in a hurry for this so might as well do in serial to simpify
for d in `find ${basedir}  -maxdepth 1 -name "Project_*" -type d | grep -v Project_CEGS | grep -v Project_Maurano`; do
    project=`basename $d`
    echo
    echo "${project}"
    base="${project}_${fc}"
    cd $d
    tar cvf /vol/mauranolab/flowcells/public_html/${base}_fastq.tar .
    
    echo
    echo "Your data are available at (user:sequencing, password:moretags)"
    echo "https://sequencing:moretags@cascade.isg.med.nyu.edu/flowcells/${base}_fastq.tar"
    
    echo
    awk -v project=${project} -F "\t" 'BEGIN {OFS="\t"} NR==1 || $1==project' ${basedir}/readcounts.txt
done

#echo "Tarballing flowcell data directory"
#base=${fc}
#cd /vol/mauranolab/flowcells/data/$fc/
#echo "Your data are available at (user:sequencing, password:moretags)"
#echo "https://sequencing:moretags@cascade.isg.med.nyu.edu/flowcells/${base}_data.tar"
#tar cvf /vol/mauranolab/flowcells/public_html/${base}_data.tar .


echo
echo "Done!!!"
date
