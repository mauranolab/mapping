#!/bin/bash

if [ "$#" -ge 1 ]; then
    dir="${1}"
else
    dir="."
fi

echo
echo "Errors in log files"
find . -regextype posix-awk -maxdepth 2 -regex ".+\/[^\/]+(BS|SRR|GSM|ENCL)[^\/]+\/.+" -and -name "*.o*" | xargs grep -i error

echo
echo "UCSC track links"
projectdir=`pwd | perl -pe 's/^\/vol\/mauranolab\/mapped\///g;'`
projectdir=${projectdir//\//\\/}

echo "Density tracks"
grep -h "Making density track" -A 1 ${dir}/*/analysis.* | grep -v cegsvectors | grep -v "Making density track" | awk '$0!="--"'

echo
echo "Coverage tracks"
grep -h "Making coverage track" -A 1 ${dir}/*/analysis.* | grep -v cegsvectors | grep -v "Making coverage track" | awk '$0!="--"'

echo
echo "variant tracks"
grep -h "variants.bb" ${dir}/*/analysis.* | grep -v cegsvectors

echo
echo "BAM tracks"
grep -h "Making BAM track" -A 1 ${dir}/*/analysis.* | grep -v cegsvectors | grep -v "Making BAM track" | awk '$0!="--"'


set -eu -o pipefail


echo
echo "Read count summaries"
find ${dir} -regextype posix-awk -maxdepth 2 -regex ".+\/[^\/]+\-[^\/]+\/.+" -and -name "analysis.*.o*" |
xargs grep -h -A 11 "Num_sequenced_reads" |
awk 'BEGIN {skip=0} $0=="--" {skip=0;next} $0=="" {skip=1;next} skip==0' |
grep -v -e "syntax[ _]error" | grep -v "Num_SE_" | grep -v "Prop_SE_" |
#Strip % symbols
perl -p -e's/[ \t]\(?([\d\.]+)\%\)?/\t$1/g;' |
perl -p -e's/ /_/g;' |
#Convert the last hyphen in the Sample_Name column to a tab
perl -p -e's/\-([^\-]+\t)/\t\1/g;' |
awk -F "\t" 'BEGIN {OFS="\t"; print "Label", "Value", "Sample_Name", "BS", "Genome"} $2!="" {print $1, $2, $4, $5, $6} $3!="" {gsub("Num_", "Pct_", $1); if($3=="NA%") {$3="NA"} print $1, $3, $4, $5, $6}' |
cut -f1-6 > $TMPDIR/readcounts.summary.long.txt


cd ${dir}
Rscript --quiet --no-save - << 'EOF'
#read() -- smarter wrapper around read.table
#BUGBUG who knows why i'm getting "zcat: stdout: Broken pipe" on first read.table
#TODO seems I need to use gzcat on the mac and zcat on linux
read <- function(filename, nrows=100, header=F, col.classes=NULL, stringsAsFactors=F, sep="\t", comment.char="", quote="", ...) {
	#BUGBUG should print example line of file upon failure
	#TODO can't use a pipe or specify which columns to use
	#BUGBUG when the first column is all characters, R tries to use it as row names. But this goofs up my hardcoded colClasses. row.names=NULL doesn't work in the function, though it works interactively.
	
	filename <- Sys.glob(filename)
	if(length(filename) == 0 || !file.exists(filename)) {
		stop("ERROR read() -- ", filename, " does not exist!\n")
	} else if(length(filename) > 1) {
		stop("ERROR read() -- Multiple files were matched: ",  paste(filename, collapse=", "), "!\n")
	}
	#TODO seems I need to use gzcat on the mac and zcat on linux
	#readCmd <- paste("gzcat -S \"\" -f \"", filename, "\"", sep="")
	readCmd <- paste("zcat -f \"", filename, "\"", sep="")
#	cat("reading with:", readCmd, "\n")

	#Let UNIX tell us how much memory to allocate
	wcpipe <- pipe(paste(readCmd, " | wc -l | awk '{print $1}' ", sep=""))
	num.lines <- as.integer(readLines(wcpipe))
	close(wcpipe)
	
	if(num.lines==0) {
		message("WARNING read() -- input file had zero lines")
		return(NULL)
	} else {
		if(is.null(col.classes)) {
			#Let R estimate the column classes so we can hardcode it later
			#Hijack nrows parameter for the number of rows upon which to estimate classes (the real number of rows is used for the final read.table)
		
			shortcopy <- read.table(pipe(readCmd), sep=sep, comment.char = comment.char, quote = quote, strip.white = TRUE, stringsAsFactors = F, nrows=nrows, header=header, ...)
			guessed.classes <- lapply(shortcopy, typeof)
		} else {
			guessed.classes <- col.classes
		}
		
		return(read.table(pipe(readCmd), sep=sep, comment.char = comment.char, quote = quote, strip.white = TRUE, stringsAsFactors = stringsAsFactors, nrows=num.lines, colClasses=guessed.classes, header=header, ...))
	}
}

data <- read(paste(Sys.getenv("TMPDIR"), "/readcounts.summary.long.txt", sep=""), header=T)
data <- reshape(data, idvar=c("Genome", "Sample_Name", "BS"), timevar="Label", direction="wide")
data <- data[,!apply(data, MARGIN=2, FUN=function(x) {all(is.na(x))})]
colnames(data) <- gsub("^Value\\.", "", colnames(data))

colnames(data) <- gsub("^Pct_", "Prop_", colnames(data))
#fails for columns with NA
data[,grepl("^Prop_", colnames(data))] <- apply(data[,grepl("^Prop.", colnames(data))], MARGIN=c(1,2), FUN=function(x) {as.numeric(x)/100})

data <- data[order(data[,"Genome"], data[,"BS"]),]
print.data.frame(data, row.names=F)
write.table(data, row.names=F, col.names=T, quote=F, file="readcounts.summary.txt", append=F, sep="\t")
EOF
cd - > /dev/null


echo
echo "Done!!!"
date
