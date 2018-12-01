#!/bin/bash

echo
echo "Errors in log files"
grep -i error */*.o*


set -eu -o pipefail


grep -h -A 11 "Num_sequenced_reads" *BS*/analysis.*.o* | awk 'BEGIN {skip=0} $0=="--" {skip=0;next} $0=="" {skip=1;next} skip==0' | grep -v -e "syntax[ _]error" | grep -v "Num_SE_" | grep -v "Prop_SE_" | 
#awk '$1!="$1!="Num_hotspots" && SPOT"' |
perl -p -e's/[ \t]\(?([\d\.]+)\%\)?/\t$1/g;' -e's/ /_/g;' -e's/\-BS/\tBS/g;' | 
#perl -pe 's/^(Genomic_coverage\t[^\t]+\t)/\1\t/g;' |
awk -F "\t" 'BEGIN {OFS="\t"; print "Label", "Value", "Cell_type", "BS", "Genome"} $2!="" {print $1, $2, $4, $5, $6} $3!="" {gsub("Num_", "Pct_", $1); if($3=="NA%") {$3="NA"} print $1, $3, $4, $5, $6}' | cut -f1-6 > /tmp/t.txt

R --quiet --no-save << 'EOF'
data <- read("/tmp/t.txt", header=T)
data <- reshape(data, idvar=c("Genome", "Cell_type", "BS"), timevar="Label", direction="wide")
data <- data[,!apply(data, MARGIN=2, FUN=function(x) {all(is.na(x))})]
colnames(data) <- gsub("^Value\\.", "", colnames(data))

colnames(data) <- gsub("^Pct_", "Prop_", colnames(data))
#gives error for columns with NA
data[,grepl("^Prop_", colnames(data))] <- apply(data[,grepl("^Prop.", colnames(data))], MARGIN=c(1,2), FUN=function(x) {as.numeric(x)/100})

data <- data[order(data[,"Genome"], data[,"BS"], decreasing = c(TRUE, FALSE)),]
print.data.frame(data, row.names=F)
write.table(data, row.names=F, col.names=T, quote=F, file="readcounts.summary.txt", append=F, sep="\t")
EOF


echo
echo "Done!!!"
date
