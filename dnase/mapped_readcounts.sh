#!/bin/bash

echo
echo "Errors in log files"
grep -i error */*.o*


set -eu -o pipefail


find . -regextype posix-awk -maxdepth 2 -regex ".+\/[^\/]+(BS|SRR|GSM)[^\/]+\/.+" -and -name analysis.*.o* | xargs grep -h -A 11 "Num_sequenced_reads" | awk 'BEGIN {skip=0} $0=="--" {skip=0;next} $0=="" {skip=1;next} skip==0' | grep -v -e "syntax[ _]error" | grep -v "Num_SE_" | grep -v "Prop_SE_" | 
#awk '$1!="$1!="Num_hotspots" && SPOT"' |
perl -p -e's/[ \t]\(?([\d\.]+)\%\)?/\t$1/g;' -e's/ /_/g;' -e's/\-(BS|SRR|GSM)/\t\1/g;' | 
#perl -pe 's/^(Genomic_coverage\t[^\t]+\t)/\1\t/g;' |
awk -F "\t" 'BEGIN {OFS="\t"; print "Label", "Value", "Cell_type", "BS", "Genome"} $2!="" {print $1, $2, $4, $5, $6} $3!="" {gsub("Num_", "Pct_", $1); if($3=="NA%") {$3="NA"} print $1, $3, $4, $5, $6}' | cut -f1-6 > $TMPDIR/readcounts.summary.long.txt


R --quiet --no-save << 'EOF'
#read() -- smarter wrapper around read.table
#BUGBUG who knows why i'm getting "zcat: stdout: Broken pipe" on first read.table
#TODO seems I need to use gzcat on the mac and zcat on linux
read <- function(filename, nrows=100, header=F, col.classes=NULL, stringsAsFactors=F, sep="\t", comment.char="", ...) {
	#BUGBUG should print example line of file upon failure
	#TODO can't use a pipe or specify which columns to use
	#BUGBUG when the first column is all characters, R tries to use it as row names. But this goofs up my hardcoded colClasses. row.names=NULL doesn't work in the function, though it works interactively.
	
	filename <- Sys.glob(filename)
	
	if(!file.exists(filename)) {
		stop(filename, " does not exist!\n")
	}
	readCmd <- paste("zcat -f \"", filename, "\"", sep="")
#	cat("reading with:", readCmd, "\n")
	if(is.null(col.classes)) {
		#Let R estimate the column classes so we can hardcode it later
		#Hijack nrows parameter for the number of rows upon which to estimate classes (the real number of rows is used for the final read.table)
		
		shortcopy <- read.table(pipe(readCmd), sep=sep, comment.char = comment.char, quote = "", strip.white = TRUE, stringsAsFactors = F, nrows=nrows, header=header, ...)
		guessed.classes <- lapply(shortcopy, typeof)
	} else {
		guessed.classes <- col.classes
	}
	
	#Let UNIX tell us how much memory to allocate
	wcpipe <- pipe(paste(readCmd, " | wc -l | awk '{print $1}' ", sep=""))
	num.lines <- as.integer(readLines(wcpipe))
	close(wcpipe)
	
	return(read.table(pipe(readCmd), sep=sep, comment.char = comment.char, quote = "", strip.white = TRUE, stringsAsFactors = stringsAsFactors, nrows=num.lines, colClasses=guessed.classes, header=header, ...))
}

data <- read(paste(Sys.getenv("TMPDIR"), "/readcounts.summary.long.txt", sep=""), header=T)
data <- reshape(data, idvar=c("Genome", "Cell_type", "BS"), timevar="Label", direction="wide")
data <- data[,!apply(data, MARGIN=2, FUN=function(x) {all(is.na(x))})]
colnames(data) <- gsub("^Value\\.", "", colnames(data))

colnames(data) <- gsub("^Pct_", "Prop_", colnames(data))
#gives error for columns with NA
data[,grepl("^Prop_", colnames(data))] <- apply(data[,grepl("^Prop.", colnames(data))], MARGIN=c(1,2), FUN=function(x) {as.numeric(x)/100})

data <- data[order(data[,"Genome"], data[,"BS"]),]
print.data.frame(data, row.names=F)
write.table(data, row.names=F, col.names=T, quote=F, file="readcounts.summary.txt", append=F, sep="\t")
EOF


echo
echo "Done!!!"
date
