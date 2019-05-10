# Trackhub generation from mapping pipeline
This uses the mauranolab fork of daler trackhub: https://daler.github.io/trackhub/

For info on files, hub links, etc. see:
https://cascade.isg.med.nyu.edu/~cadlej01/blog/tag_Hubs.shtml#1550256952


## Procedure for mapping ENCODE data

```
base=/vol/isg/encode/dnase
src=/vol/mauranolab/mapped/src/dnase
```

### Download fastq files info and metadata

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

#Re-download files that are missing (I think wget cleans up incomplete downloads but because the stdout is a mess you don't see the error message)
ls | fgrep -v -f - ../files.txt | xargs -n1 -P20 -L 1 wget

#verify md5 hashes
cat metadata.tsv | mlr --tsv --ocsv --headerless-csv-output --ofs "  " --ors lf rename -g -r ' ,_' then sort -f File_accession then cut -o -f md5sum,File_accession then put '$File_accession=$File_accession . ".fastq.gz"' > metadata.md5sum.txt
bqsub -j y -N md5sum.validate "md5sum -c metadata.md5sum.txt"

cd ..
```

### Download all JSON files
6) mkdir JSON && cd JSON
7) To download all JSON files per fastq-file run the following command
```
for i in `tail -n +2 ../metadata.tsv | awk '{print $1}'`; do echo https://www.encodeproject.org/files/${i}/?format=json -O $i;done | xargs -n1 -P6 -L 1 wget
```

### SampleID file
8) This step takes the metadata file as input, the folder with all JSON files, and a output filename  
```
${src}/trackhub/extractDSfromENCODE_JSON.py metadata.tsv -j JSON/ -o /tmp/SampleIDs.raw.tsv
mlr --tsv sort -f Name,Assay,DS /tmp/SampleIDs.raw.tsv > SampleIDs.raw.tsv
```

9) Here it is important to change names in the Name column to be consistent across all samples. The biosample names chosen here will continue through the whole analysis. Save the edited file in SampleIDs.tsv

```
#Double check that each library ID is unique for a given combination of cell type and histone mark/assay
mlr --tsv cut -f SampleAccession,Name,Assay then sort -f Name,Assay SampleIDs.raw.tsv | uniq | cut -f1 | hist | awk '$1>1'

#miller can parse JSON files for ease of viewing, e.g.:
mlr --json --jvstack head -n 2  JSON/ENCFF001HZO
```

10) Now make softlinks to fastq files to connect all files from same library
```
mkdir renamed && cd renamed
awk  -F'\t' '$1==$7 {print "ln -s ../fastq/" $7 ".fastq.gz " $1 "_" $3 "_R1.fastq.gz"; if($8!="") {print "ln -s ../fastq/" $8 ".fastq.gz " $1 "_" $3 "_R2.fastq.gz"}}' ../SampleIDs.tsv > makelinks.sh  
source makelinks.sh
cd ..  
```


### Mapped folder with inputs.txt and submit.sh
```
mkdir mapped && cd mapped
awk -v base="${base}/renamed" -F "\t" 'BEGIN {OFS="\t"} $1==$7 {print base "/" $1 "_" $3 "_R1.fastq.gz", $3; if($8!="") {print base "/" $1 "_" $3 "_R2.fastq.gz", $2}}' ../SampleIDs.tsv |
sort -k2,2 | cut -f1 > inputs.txt 

tail +2 ../SampleIDs.tsv |
perl -pe 's/ /_/g;' -e 's/\-/_/g;' -e "s/[\(\)\'%]//g;" |
awk -v src=${src} -F "\t" '{print src "/submit.sh hg38_noalt mapBwaAln,dnase " $5  " " $3}' |
#Use this line instead for ChIP-seq
#awk -v src=${src} -F "\t" '{print src "/submit.sh hg38_noalt mapBwaAln,chipseq " $5 "-" $9  " " $3}' |
sort |
uniq |
sort -k4,3 |
grep -v GroupID > submitJobs.sh

cat <<EOF >> submitJobs.sh

rqsub -b y -j y -N analyzeInserts -hold_jid \`cat sgeid.analysis | perl -pe 's/\n/,/g;'\` "$src/analyzeInserts.R 500"
bqsub -j y -N mapped_readcounts -hold_jid \`cat sgeid.analysis | perl -pe 's/\n/,/g;'\` $src/mapped_readcounts.sh
EOF
cd ..
```

### Annotate samples and create trackhub
This should be run after all DNase samples have been processed. 
samplesforTrackhub.R combines metadata from SamplesForTrackhub.tsv with the pipeline output of analysis.sh.
```
#Prep annotation
cd /vol/isg/encode/dnase/mapped/
${src}/trackhub/samplesforTrackhub.R --inputfile ../SampleIDs.tsv --project humanENCODEdnase --out ../SamplesForTrackhub.tsv 2>&1 | grep -E -i "(warning|error)"

cd /vol/isg/encode/chipseq/mapped/
${src}/trackhub/samplesforTrackhub.R --inputfile ../SampleIDs.tsv --project humanENCODEchipseq --out ../SamplesForTrackhub.tsv 2>&1 | grep -E -i "(warning|error)"

cd /vol/isg/encode/mouseencode/mapped/
${src}/trackhub/samplesforTrackhub.R --inputfile ../SampleIDs.tsv --project mouseENCODEdnase --out ../SamplesForTrackhub.tsv 2>&1 | grep -E -i "(warning|error)"

cd /vol/isg/encode/mouseencode_chipseq/mapped/
${src}/trackhub/samplesforTrackhub.R --inputfile ../SampleIDs.tsv --project mouseENCODEchipseq --out ../SamplesForTrackhub.tsv 2>&1 | grep -E -i "(warning|error)"
```
MakeTrackhub.py creates the trackdb file. You will need to create the trackhub/hub.txt and trackhub/genomes.txt landing page.
```
#MakeTrackhub.py requires the forked daler/trackhub from mauranolab. It can be installed locally:
pip install --upgrade --user . 
#Or directly from git
pip install --upgrade --user git+git://github.com/mauranolab/trackhub.git

#Make hub
mkdir -p /vol/isg/encode/trackhub/mm10 /vol/isg/encode/trackhub/hg38
cd /vol/isg/encode/trackhub
${src}/trackhub/MakeTrackhub.py --genome hg38 --supertrack ENCODE_DNase-seq --URLbase https://cascade.isg.med.nyu.edu/mauranolab/encode/dnase/mapped/ --subgroupnames Replicate,Age /vol/isg/encode/dnase/SamplesForTrackhub.tsv > hg38/trackDb.txt
${src}/trackhub/MakeTrackhub.py --genome hg38 --supertrack ENCODE_ChIP-seq --URLbase https://cascade.isg.med.nyu.edu/mauranolab/encode/chipseq/mapped/ --subgroupnames Replicate,Age /vol/isg/encode/chipseq/SamplesForTrackhub.tsv >> hg38/trackDb.txt

${src}/trackhub/MakeTrackhub.py --genome mm10 --supertrack ENCODE_DNase-seq --URLbase https://cascade.isg.med.nyu.edu/mauranolab/encode/mouseencode/mapped/ --subgroupnames Replicate,Age /vol/isg/encode/mouseencode/SamplesForTrackhub.tsv > mm10/trackDb.txt
${src}/trackhub/MakeTrackhub.py --genome mm10 --supertrack ENCODE_ChIP-seq --URLbase https://cascade.isg.med.nyu.edu/mauranolab/encode/mouseencode_chipseq/mapped/ --subgroupnames Replicate,Age /vol/isg/encode/mouseencode_chipseq/SamplesForTrackhub.tsv >> mm10/trackDb.txt
```
For debugging trackhub, add "udcTimeout=1&" to browser URL

You can also debug error messages with:
```
hubCheck -udcDir=/tmp https://cascade.isg.med.nyu.edu/mauranolab/encode/trackhub/hub.txt
```
