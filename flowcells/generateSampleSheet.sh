#!/bin/bash
set -eu -o pipefail


if [ "$#" -ge 1 ]; then
    printPairs="$1"
else
    printPairs="FALSE"
fi

src="/vol/mauranolab/mapped/src"

echo "Paste in a full flowcell entry (including #header lines, followed by a newline and CTRL-D"

#Boilerplate header
cat ${src}/flowcells/SampleSheet.template.txt > $TMPDIR/SampleSheet.withlane.csv


###Parse the sequencing sheet info from STDIN
#Clean up single-line spacer between header and data table
perl -pe 's/^\t+$//g;' |
awk -F "\t" 'BEGIN {OFS="\t"; parse=0} {print} $0=="" && parse==0 {parse=1; print "#Sample Name", "Sample #", "Lab", "Made By", "Sample Type", "Barcode 1 (i7)", "Barcode 2 (i5)", "R1 Trim (P5)", "R2 Trim (P7)", "Sequencing primer R1", "Indexing primer BC1 (i7)", "Indexing primer BC2 (i5)", "Sequencing primer R2", "Library concentration (pM)", "Request Type", "Requested reads (M)", "Read format", "Scale factor", "Relative representation", "Amount put on FC (uL)", "Sequenced reads", "Actual representation"}' |
#Clean up trailing tabs on comment lines
perl -pe 's/^(#.+[^\t])\t+$/\1/g;' |
#Also creates info.txt
tee info.txt |
#NB our sample sheet records RC for BC2/i5, which according to https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/miseq/indexed-sequencing-overview-guide-15057455-04.pdf is valid for iSeq 100, MiniSeq, NextSeq, HiSeq X, HiSeq 4000, or HiSeq 3000. Therefore for runs on NovaSeqTM 6000, MiSeqTM, HiSeq 2500, and HiSeq 2000, BC2 must be RC back to the original sequence.
awk -f ${src}/flowcells/revcomp.awk -F "\t" --source  'BEGIN {OFS=","; split("8,8", bclens, ",")} \
    $1=="#Instrument" { if($2~/NovaSeq/ || $2~/MiSeq/) {doRevComp=1} else {doRevComp=0} } \
    $1=="#Indices" && $2!="" {split($2, bclens, ",")} \
    $0!~/^#/ && $1!="" && $5!="Pool" { \
        if(bclens[1]==0) {$6="_"} \
        if(bclens[2]==0) {$7="_"} \
        split($6, bc1, "_"); \
        split($7, bc2, "_"); \
        if(doRevComp==1) { bc2[2]=revcomp(bc2[2]) } \
        print "Sample_" $2, $2, "", "",  bc1[1], toupper(substr(bc1[2], 0, bclens[1])),  bc2[1], toupper(substr(bc2[2], 0, bclens[2])), "Project_" $3, "", $19; \
    }' >> $TMPDIR/SampleSheet.withlane.csv


echo
echo
lanestatus=`cat $TMPDIR/SampleSheet.withlane.csv | awk -F "," 'BEGIN {foundEmptyLane=0; foundLane=0} $1=="Sample_ID" {parse=1; next} parse==1 {if($11=="") {foundEmptyLane=1} else {foundLane=1}} END {if(foundEmptyLane==1 && foundLane==1) {print "ERROR, incomplete lane info"} else if(foundEmptyLane==1) {print "NOLANE"} else if (foundLane==1) {print "LANE"} else {print "IMPOSSIBLE"}}'`
case "${lanestatus}" in
LANE)
    echo "Samples will be demuxed by lane"
    mv $TMPDIR/SampleSheet.withlane.csv SampleSheet.csv
    ;;
NOLANE)
    cut -d "," -f1-10 $TMPDIR/SampleSheet.withlane.csv > SampleSheet.csv
    ;;
*)
    echo "${lanestatus}"
    exit 1
    ;;
esac


echo
echo
echo "Sample sheet validation"
#validate date format as YYYY-MM-DD
readgroup_date=`awk -F "\t" 'BEGIN {OFS="\t"} $1=="#Load date" {print $2}' info.txt`
if [[ ! "${readgroup_date}" =~ 20[12][0-9]\-[0-2][0-9]\-[0123][0-9] ]]; then
    echo "WARNING: invalid FC load date ${readgroup_date}"
fi


#TODO check read count total
#TODO Better to use pool counts?
awk -F "\t" 'BEGIN {OFS="\t"} $1=="#Expected reads" {expectedReads=$2; programmedReads=0} \
$0!~/^#/ && $1!="" && $5!="Pool" { programmedReads+= $16 } \
END { if(programmedReads != expectedReads) { \
        print "WARNING: Total of " programmedReads " reads requested for a FC that yields " expectedReads; \
    } \
}' info.txt


awk -F "\t" 'BEGIN {OFS="\t"} $1=="#FC kit" { kitcycles="NA"; \
    if($2=="HighOutput_75cycle") {kitcycles=92} \
    else if($2=="MidOutput_150cycle" || $2=="HighOutput_150cycle") {kitcycles=168} \
    else if($2=="MidOutput_300cycle" || $2=="HighOutput_300cycle") {kitcycles=318} \
} \
$1=="#Format" {split($2, fcreadformat, ",")} \
$1=="#Indices" {split($2, fcindices, ",")} \
END {\
    totalcycles=fcreadformat[1]+fcreadformat[2]+fcindices[1]+fcindices[2]; \
    if(totalcycles != kitcycles) { \
        print "WARNING: " totalcycles " cycles programmed while the kit contains reagents for " kitcycles \
    } \
}' info.txt


awk -F "\t" 'BEGIN {OFS="\t"} $1=="#Format" {split($2, fcreadformat, ",")} \
$0!~/^#/ && $1!="" { \
    gsub(/^[SP]E /, "", $17); \
    split($17, sampleReadFormat, ","); \
    if(sampleReadFormat[1]>fcreadformat[1] || sampleReadFormat[2]>fcreadformat[2]) { \
        print "WARNING: Sample " $1 "-" $2 " requires longer read lengths (" sampleReadFormat[1] "," sampleReadFormat[2] ") than programmed" \
    } \
}' info.txt


cat info.txt | awk -F "\t" 'BEGIN {OFS="\t"; inheader=1; lastBSnum=0} \
$1=="#Sample Name" && inheader==1 { inheader=0; next } \
inheader==0 { \
    curBSnum=substr($2, 3, 5); \
    if(curBSnum < lastBSnum) { \
        print "WARNING: Sample " $1 "-" $2 " appears out of order"; \
    } \
    lastBSnum=curBSnum; \
}'


#TODO not sure how to get a single quote inside
#BUGBUG false alarm for - in chipseq samples
cat info.txt | awk -F "\t" 'BEGIN {OFS="\t"; inheader=1; lastBSnum=0} \
$1=="#Sample Name" && inheader==1 { inheader=0; next } \
inheader==0 && $1~/[\-%\(\)\"\/\. ]/ { print "WARNING: Sample name for " $1 "-" $2 " contains invalid characters"; }'


cat info.txt | awk -F "\t" 'BEGIN {OFS="\t"; inheader=1; lastBSnum=0} \
$1=="#Sample Name" && inheader==1 { inheader=0; next } \
inheader==0 && $1!~/\-/ && $5=="ChIP-seq" { print "WARNING: Sample name for " $1 "-" $2 " missing ChIP-seq epitope"; }'


echo "Finished validation"


echo
echo

Rscript --quiet --no-save - ${printPairs} << 'EOF'
suppressPackageStartupMessages(library(reshape))

printPairs <- commandArgs(TRUE)[1]

data <- read.table("SampleSheet.csv", header=T, skip=19, sep=",", comment.char = "", quote = "", strip.white = TRUE, stringsAsFactors = F)

if(nrow(data) == 0) {
	warning("ERROR: No data found!")
	quit(1, save="no", status=0)
}

data[is.na(data[,"index"]),"index"] <- ""
data[is.na(data[,"index2"]),"index2"] <- ""

#I can't find a builtin hamming distance implementation for R
strdist <- function(x,y) {
	if(length(x) != length(y)) {
		stop("ERROR: diff lengths unsupported!")
	}
	xlist <- unlist(strsplit(x, ""))
	ylist <- unlist(strsplit(y, ""))
	#N's count as wildcards and match anything
	length(which(xlist != ylist & xlist!="N" & ylist!="N"))
}


maxBC1len <- max(sapply(data[,"index"], FUN=function(x) {nchar(x)}))
#BUGBUG fails for runs without BC2
maxBC2len <- max(sapply(data[,"index2"], FUN=function(x) {nchar(x)}))
minBClen <- 4

countBCcollisions <- function(data, bc1len=8, bc2len=8) {
    data[,"index"] <- substr(data[,"index"], start=1, stop=bc1len)
    data[,"index2"] <- substr(data[,"index2"], start=1, stop=bc2len)
    
    bcs <- apply(data[, c("index", "index2")], MARGIN=1, FUN=paste, collapse="_")
    #bcs <- data$index
    #bcs <- data$index2
    
    
    numSamplesTooClose <- 0
    results <- matrix(NA, nrow=length(bcs), ncol=length(bcs))
    for(i in 1:length(bcs)) {
        bc1 <- bcs[i]
        for(j in 1:length(bcs)) {
            if(i<j & bc1!="_") {
                bc2 <- bcs[j]
                if(bc2!="_") {
                    results[i,j] = strdist(bc1, bc2)
                    if(results[i,j] <= 2) {
                        if(printPairs) {
                            cat("pair ", i, " (", bc1, ") and ", j, " (", bc2, ") differ by only ", results[i,j], "\n", sep="") 
                        }
                        numSamplesTooClose <- numSamplesTooClose + 1
                    }
                }
            }
        }
    }
    
    #return pairwise matrix of distances
    #return(results)
    
    #return numSamplesTooClose
    return(numSamplesTooClose)
}


if (minBClen <= maxBC1len) {
    bc1range <- c(0,minBClen:maxBC1len)
} else {
    bc1range <- maxBC1len
}

if (minBClen <= maxBC2len) {
    bc2range <- c(0,minBClen:maxBC2len)
} else {
    bc2range <- maxBC2len
}

results <- data.frame()
for(bc1len in bc1range) {
    for(bc2len in bc2range) {
        if(printPairs) {
            cat(paste0("\nbc1len=", bc1len, ";bc2len=", bc2len, ":\n"))
        }
        numSamplesTooClose <- countBCcollisions(data, bc1len, bc2len)
        results <- rbind(results, data.frame(bc1len=bc1len, bc2len=bc2len, numSamplesTooClose=numSamplesTooClose))
    }
}

cat("\nNumber of samples too close by BC sequencing lengths:\n")
if(maxBC2len==0) {
    #Not great but just print results for SE runs
    print.data.frame(results, row.names=F)
} else {
    results.wide <- cast(bc1len~bc2len, value="numSamplesTooClose", data=results, fun.aggregate=sum, add.missing=T, fill=NA)
    #Cast adds bc2 as colnames but not rownames
    rownames(results.wide) <- paste("bc1", results.wide[,1])
    results.wide <- results.wide[,-1]
    print.data.frame(results.wide, row.names=T)
}
EOF


#For plotting pairwise edit distances
#suppressPackageStartupMessages(library(ggplot2))
#old <- theme_set(theme_classic(base_size=7)) #pdf
#old <- theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
#
#
#tmp <- melt(results)
#tmp[!is.na(tmp$value) & tmp$value >=4, "value"] <- NA
#
#ggplot(tmp) + 
#geom_tile(aes(x=X1, y=X2, fill = value), colour = "white") + 
#scale_fill_gradient(low = "white", high = "steelblue")

