# DNase pipeline

## Download fastq files info and metadata

1) Visit www.encodeproject.org  
2) In top bar, select data and search  
3) Under Assay --> Select DNase-Seq  
   Under Organism --> Select Homo Sapiens  
   Under Available data --> Select fastq  

Select download  
This results in a file named files.txt, which contains all fastq files for the selected experiments, as well as metadata for each file

4) Copy the first line to the browser, and a metadata file will be downloaded, metadata.tsv. Also downlaod fastq files to a folder called fastq

## Download all JSON files
5) Create a folder to store all JSON files in
6) To download all JSON files per fastq-file run the following command
```
for i in `tail -n +2 "metadata.tsv"|awk '{print $1}'`; do echo https://www.encodeproject.org/files/${i}/?format=json -O JSON/$i;done| xargs -n1 -P20 -L 1 wget
```


## SampleID file
7) This step takes the metadata file as input, the folder with all JSON files, and a output filename  
8) Here it is important to change names to be consistent across all samples. The biosample names chosen here will continue through the whole analysis

```
python FindDSnumber.py metadata.tsv -j JSON/ -o SampleIDs.tsv
```

## Softlinks to fastq files to connect all files from same library
```
mkdir -p renamed
cd renamed
awk  -F'\t' '$1==$6 {print "ln -s ../fastq/" $6 ".fastq.gz " $1 "_" $2 "_R1.fastq.gz"; if($7!="") {print "ln -s ../fastq/" $7 ".fastq.gz " $1 "_" $2 "_R2.fastq.gz"}}' ../SampleIDs.tsv > makelinks.sh  
source makelinks.sh ;
cd ..  
```


## Mapped folder with inputs.txt and submit.sh
```
mkdir -p mapped
cd mapped
awk -v base="/vol/isg/encode/dnase/renamed" -F "\t" 'BEGIN {OFS="\t"} $1==$6 {print base "/" $1 "_" $2 "_R1.fastq.gz", $2; if($7!="") {print base "/" $1 "_" $2 "_R2.fastq.gz", $2}}' ../SampleIDs.tsv |
sort -k2,2 | cut -f1 >inputs.txt 

cat ../SampleIDs.tsv |
perl -pe 's/ /_/g;' -e 's/\-/_/g;' -e "s/[\(\)\'%]//g;" |
awk -F "\t" '{print "/vol/isg/encode/dnase/src/submit.sh hg38 map_dnase " $3  " " $2}' |
sort |
uniq |
sort -k4,3 |
grep -v GroupID >submitJobs.sh

```

## Annotate samples for trackhub
```
cd /vol/isg/encode/dnase/mapped
Rscript --vanilla /vol/isg/encode/dnase/src/trackhub/samplesforTrackhub.R --file /vol/isg/encode/dnase/SampleIDs_20180502_MTM.tsv --out /vol/isg/encode/dnase201805/SamplesForTrackhub.tsv
```
## Create trackhub
This script creates the hub, and genome file at the output location, and creates a subdirectory names hg38 containing the trackhub.txt  
This should be run after all DNase samples have been processed

(requires python 2.7)
```
module remove python
module add python/2.7.10
mkdir -p trackhub/hg38
python /vol/isg/encode/dnase/src/trackhub/DalerTrackhub.py /vol/isg/encode/dnase/SamplesForTrackhub.tsv --output /vol/isg/encode/dnase201805/Encode_DNase | perl -pe 's/^track/\ntrack/g;' > trackhub/hg38/trackDb.txt
```
