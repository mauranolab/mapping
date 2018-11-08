#!/bin/bash
set -eu -o pipefail

#Boilerplate header
cat /vol/mauranolab/flowcells/src/SampleSheet.template.txt > SampleSheet.csv

#Now parse the sequencing sheet info from STDIN
awk -F "\t" 'BEGIN {OFS="\t"; parse=0} {print} $0!~/^#/ && parse==0 {parse=1; print "#Sample Name", "Sample #", "Lab", "Made By", "Sample Type", "Species", "Barcode 1 (i7)", "Barcode 2 (i5)", "R1 Trim (P5)", "R2 Trim (P7)", "Sequencing primer R1", "Indexing primer BC1 (i7)", "Indexing primer BC2 (i5)", "Sequencing primer R2", "Library concentration (pM)", "Request Type", "Requested reads (M)", "Read format", "Scale factor", "Relative representation", "Amount put on FC (uL)", "Sequenced reads", "Actual representation"}' |
#Also creates info.txt
tee info.txt |
awk -F "\t" 'BEGIN {OFS=","; split("8,8", bclens, ",")} 
    $1=="#Indices" && $2!="" {split($2, bclens, ",")}
    $0!~/^#/ && $1!="" && $1 {if(bclens[1]==0) {$7="_"} if(bclens[2]==0) {$8="_"} split($7, bc1, "_"); split($8, bc2, "_"); print "Sample_" $2, $2, "", "",  bc1[1], toupper(substr(bc1[2], 0, bclens[1])),  bc2[1], toupper(substr(bc2[2], 0, bclens[2])), "Project_" $3, "";}' >> SampleSheet.csv

#validate date format as 
readgroup_date=`awk -F "\t" 'BEGIN {OFS="\t"} $1=="#Load date" {print $2}' info.txt`
if [[ ! "${readgroup_date}" =~ 201[5-9]\-[0-2][0-9]\-[0123][0-9] ]]; then
    echo "WARNING: invalid FC load date ${readgroup_date}"
fi


R --quiet --no-save << 'EOF'
data <- read.table("SampleSheet.csv", header=T, skip=19, sep=",", comment.char = "", quote = "", strip.white = TRUE, stringsAsFactors = F)

#I can't find a builtin hamming distance implementation for R
strdist <- function(x,y) {
	if(length(x)!=length(y)) {
		stop("ERROR: diff lengths unsupported!")
	}
	length(which(unlist(strsplit(x, "")) != unlist(strsplit(y, ""))))
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
                        cat("pair ", i, " (", bc1, ") and ", j, " (", bc2, ") differ by only", results[i,j], "\n", sep="")
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


results <- data.frame()
for(bc1len in minBClen:maxBC1len) {
    for(bc2len in minBClen:maxBC2len) {
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

