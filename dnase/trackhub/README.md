# DNase pipeline - procedure for mapping ENCODE data

```
base=/vol/isg/encode/dnase
src=/vol/mauranolab/mapped/src
```

## Download fastq files info and metadata

1) Visit www.encodeproject.org  
2) In top bar, select Data --> Search  
3) Under Assay --> Select DNase-Seq or ChIP-seq
   Under Organism --> Select Homo Sapiens  
   Under Available data --> Select fastq  
4) Select download to get the link to download files.txt, which contains URLs for the metadata tsv and fastq files for the selected experiments
5) Download the actual files
```
head -1 files.txt | xargs wget -O metadata.tsv
mkdir fastq && cd fastq
cat ../files.txt | xargs -n1 -P20 -L 1 wget

#verify md5 hashes
cat metadata.tsv | mlr --tsv --ocsv --headerless-csv-output --ofs "  " --ors lf rename -g -r ' ,_' then cut -o -f md5sum,File_accession then put '$File_accession=$File_accession . ".fastq.gz"' > /tmp/md5.txt
md5sum -c /tmp/md5.txt
cd ..
```

## Download all JSON files
6) mkdir JSON && cd JSON
7) To download all JSON files per fastq-file run the following command
```
for i in `tail -n +2 ../metadata.tsv | awk '{print $1}'`; do echo https://www.encodeproject.org/files/${i}/?format=json -O $i;done | xargs -n1 -P20 -L 1 wget
cd ..
```


## SampleID file
8) This step takes the metadata file as input, the folder with all JSON files, and a output filename  
9) Here it is important to change names to be consistent across all samples. The biosample names chosen here will continue through the whole analysis

```
python ${src}/trackhub/extractDSfromENCODE_JSON.py metadata.tsv -j JSON/ -o SampleIDs.tsv
```

## Softlinks to fastq files to connect all files from same library
```
mkdir renamed && cd renamed
awk  -F'\t' '$1==$6 {print "ln -s ../fastq/" $6 ".fastq.gz " $1 "_" $2 "_R1.fastq.gz"; if($7!="") {print "ln -s ../fastq/" $7 ".fastq.gz " $1 "_" $2 "_R2.fastq.gz"}}' ../SampleIDs.tsv > makelinks.sh  
source makelinks.sh ;
cd ..  
```


## Mapped folder with inputs.txt and submit.sh
```
mkdir mapped && cd mapped
awk -v base="${base}/renamed" -F "\t" 'BEGIN {OFS="\t"} $1==$6 {print base "/" $1 "_" $2 "_R1.fastq.gz", $2; if($7!="") {print base "/" $1 "_" $2 "_R2.fastq.gz", $2}}' ../SampleIDs.tsv |
sort -k2,2 | cut -f1 > inputs.txt 

cat ../SampleIDs.tsv |
perl -pe 's/ /_/g;' -e 's/\-/_/g;' -e "s/[\(\)\'%]//g;" |
awk -F "\t" '{print "${src}/submit.sh hg38_noalt mapBwaAln,dnase " $3  " " $2}' |
sort |
uniq |
sort -k4,3 |
grep -v GroupID > submitJobs.sh
cd ..
```

## Annotate samples for trackhub
```
cd ${base}/mapped
Rscript --vanilla ${src}/trackhub/samplesforTrackhub.R --file ${base}/SampleIDs_20180502_MTM.tsv --out ${base}/SamplesForTrackhub.tsv
```
## Create trackhub
This script creates the hub, and genome file at the output location, and creates a subdirectory names hg38 containing the trackhub.txt  
This should be run after all DNase samples have been processed

(requires python 2.7)
```
module remove python
module add python/2.7.10
mkdir -p trackhub/hg38
python ${src}/trackhub/DalerTrackhub.py ${base}/SamplesForTrackhub.tsv --output ${base}/Encode_DNase | perl -pe 's/^track/\ntrack/g;' > trackhub/hg38/trackDb.txt
cd ..
```
