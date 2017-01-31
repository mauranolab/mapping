# pipeline



#Locate DSnumber 
JSON files should be downloaded to JSON/
```
py metadata/FindDSnumber.py metadata/metadata.tsv -j JSON/ -o output
```
#Annotate samples for trackhub
```
Rscript --vanilla  trackhub/samplesforTrackhub.R --file "output from FindDSnumber.py" --out output
```
#Create trackhub
This script creates the hub, and genome file at the output location, and creates a subdirectory names hg38 containing the trackhub.txt

```
python src/trackhub/DalerTrackhub.py "trackhub/samplesforTrackhub.R output" --output hub
```
