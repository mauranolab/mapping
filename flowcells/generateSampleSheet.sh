#!/bin/bash
set -eu -o pipefail

echo "Paste in a full flowcell entry (including #header lines, followed by a newline and CTRL-D"

#Boilerplate header
cat /vol/mauranolab/mapped/src/flowcells/SampleSheet.template.txt > SampleSheet.csv


###Parse the sequencing sheet info from STDIN
#Clean up single-line spacer between header and data table
perl -pe 's/^\t+$//g;' |
awk -F "\t" 'BEGIN {OFS="\t"; parse=0} {print} $0=="" && parse==0 {parse=1; print "#Sample Name", "Sample #", "Lab", "Made By", "Sample Type", "Species", "Barcode 1 (i7)", "Barcode 2 (i5)", "R1 Trim (P5)", "R2 Trim (P7)", "Sequencing primer R1", "Indexing primer BC1 (i7)", "Indexing primer BC2 (i5)", "Sequencing primer R2", "Library concentration (pM)", "Request Type", "Requested reads (M)", "Read format", "Scale factor", "Relative representation", "Amount put on FC (uL)", "Sequenced reads", "Actual representation"}' |
#Clean up trailing tabs on comment lines
perl -pe 's/^(#.+[^\t])\t+$/\1/g;' |
#Also creates info.txt
tee info.txt |
#NB our sample sheet records RC for BC2, which according to https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/miseq/indexed-sequencing-overview-guide-15057455-04.pdf is valid for iSeq 100, MiniSeq, NextSeq, HiSeq X, HiSeq 4000, or HiSeq 3000. Therefore for runs on NovaSeqTM 6000, MiSeqTM, HiSeq 2500, and HiSeq 2000, BC2 must be RC back to the original sequence.
awk -v doRevComp=0 -f ~/lib/revcomp.awk -F "\t" --source  'BEGIN {OFS=","; split("8,8", bclens, ",")} 
    $1=="#Indices" && $2!="" {split($2, bclens, ",")}
    $0!~/^#/ && $1!="" && $5!="Pool" {if(bclens[1]==0) {$7="_"} if(bclens[2]==0) {$8="_"} split($7, bc1, "_"); split($8, bc2, "_"); if(doRevComp==1) {bc2[2]=revcomp(bc2[2])} print "Sample_" $2, $2, "", "",  bc1[1], toupper(substr(bc1[2], 0, bclens[1])),  bc2[1], toupper(substr(bc2[2], 0, bclens[2])), "Project_" $3, "";}' >> SampleSheet.csv

echo
echo
echo "Sample sheet validation"
#validate date format as YYYY-MM-DD
readgroup_date=`awk -F "\t" 'BEGIN {OFS="\t"} $1=="#Load date" {print $2}' info.txt`
if [[ ! "${readgroup_date}" =~ 201[5-9]\-[0-2][0-9]\-[0123][0-9] ]]; then
    echo "WARNING: invalid FC load date ${readgroup_date}"
fi


#TODO check read count total
awk -F "\t" 'BEGIN {OFS="\t"} $1=="#Expected reads" {expectedReads=$2; programmedReads=0} \
$0!~/^#/ && $1!="" { programmedReads+= $17 } \
END { if(programmedReads != expectedReads) { \
        print "WARNING: Total of " programmedReads " reads requested for a FC that yields " expectedReads; \
    } \
}' info.txt


awk -F "\t" 'BEGIN {OFS="\t"} $1=="#FC kit" { kitcycles="NA"; \
    if($2=="HighOutput_75cycle") {kitcycles=92} \
    else if($2=="MidOutput_150cycle" || $2=="HighOutput_150cycle") {kitcycles=168} \
    else if($2=="MidOutput_300cycle" || $2=="HighOutput_300cycle") {kitcycles=347} \
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
    gsub(/^[SP]E /, "", $18); \
    split($18, sampleReadFormat, ","); \
    if(sampleReadFormat[1]>fcreadformat[1] || sampleReadFormat[2]>fcreadformat[2]) { \
        print "WARNING: Sample " $1 "-" $2 " requires longer read lengths (" sampleReadFormat[1] "," sampleReadFormat[2] ") than programmed" \
    } \
}' info.txt


#TODO not sure how to get a single quote inside
cat info.txt | awk -F "\t" 'BEGIN {OFS="\t"; inheader=1; lastBSnum=0} \
$1=="#Sample Name" && inheader==1 { inheader=0; next } \
inheader==0 && $1~/[\-%\(\)\"\/\. ]/ { print "WARNING: Sample name for " $1 "-" $2 " contains invalid characters"; } \
inheader==0 { \
    curBSnum=substr($2, 3, 5); \
    if(curBSnum < lastBSnum) { \
        print "WARNING: Sample " $1 "-" $2 " appears out of order"; \
    } \
    lastBSnum=curBSnum; \
}'

echo "Finished validation"


echo
echo

R --quiet --no-save << 'EOF'
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
                    if(results[i,j] <=2) {
                        cat("pair ", i, " (", bc1, ") and ", j, " (", bc2, ") differ by only ", results[i,j], "\n", sep="")
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
    bc1range <- minBClen:maxBC1len
} else {
    bc1range <- maxBC1len
}

if (minBClen <= maxBC2len) {
    bc2range <- minBClen:maxBC2len
} else {
    bc2range <- maxBC2len
}

results <- data.frame()
for(bc1len in bc1range) {
    for(bc2len in bc2range) {
        cat(paste0("bc1len=", bc1len, ";bc2len=", bc2len, ":\n"))
        numSamplesTooClose <- countBCcollisions(data, bc1len, bc2len)
        results <- rbind(results, data.frame(bc1len=bc1len, bc2len=bc2len, numSamplesTooClose=numSamplesTooClose))
        cat("\n")
    }
}

cat("\nNumber of samples too close by BC sequencing lengths:\n")
print.data.frame(results, row.names=F)
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

