#!/bin/env Rscript

print(date())


old <- theme_set(theme_classic(base_size=7)) #pdf
old <- theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Try to work around "unable to start device PNG" on ISG cluster
options(bitmapType="cairo") 

argv <- commandArgs(TRUE)
prefix <- argv[1]
outbase <- argv[2]



cat("\n\nOverlap between libraries by read depth\n")
data <- read(paste0(outbase, ".AllBCs.txt"), header=T)

data$DNA.bin <- cut.pretty(data[,"DNA"], breaks=unique(quantile(data[,"DNA"], probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, pretty.labels="left")
data$DNA.decile <- factor(cut(data[,"DNA"], breaks=unique(quantile(data[,"DNA"], probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))

data$RNA.bin <- cut.pretty(data[,"RNA"], breaks=unique(quantile(data[,"RNA"], probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, pretty.labels="left")
data$RNA.decile <- factor(cut(data[,"RNA"], breaks=unique(quantile(data[,"RNA"], probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))

data$iPCR.bin <- cut.pretty(data[,"iPCR"], breaks=unique(quantile(data[,"iPCR"], probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, pretty.labels="left")
data$iPCR.decile <- factor(cut(data[,"iPCR"], breaks=unique(quantile(data[,"iPCR"], probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))

cat("DNA library overlap with other libraries\n")
DNA.overlaps.byReadDecile <- cbind(sample=prefix, source="DNA", summaryBy(RNA+iPCR~DNA.bin+DNA.decile, FUN=pctTrue, data=subset(data, !is.na(DNA.bin)), keep.names=T))
colnames(DNA.overlaps.byReadDecile)[colnames(DNA.overlaps.byReadDecile)=="DNA.bin"] <- "reads.bin"
colnames(DNA.overlaps.byReadDecile)[colnames(DNA.overlaps.byReadDecile)=="DNA.decile"] <- "reads.decile"
print.data.frame(DNA.overlaps.byReadDecile, row.names=F, right=F)
write.table(melt.data.frame(DNA.overlaps.byReadDecile), row.names=F, col.names=T, quote=F, file=paste0(outbase, ".DNA.intersections.byReadDecile.txt"), append=F, sep="\t")

cat("RNA library overlap with other libraries\n")
RNA.overlaps.byReadDecile <- cbind(sample=prefix, source="RNA", summaryBy(DNA+iPCR~RNA.bin+RNA.decile, FUN=pctTrue, data=subset(data, !is.na(RNA.bin)), keep.names=T))
colnames(RNA.overlaps.byReadDecile)[colnames(RNA.overlaps.byReadDecile)=="RNA.bin"] <- "reads.bin"
colnames(RNA.overlaps.byReadDecile)[colnames(RNA.overlaps.byReadDecile)=="RNA.decile"] <- "reads.decile"
print.data.frame(RNA.overlaps.byReadDecile, row.names=F, right=F)
write.table(melt.data.frame(RNA.overlaps.byReadDecile), row.names=F, col.names=T, quote=F, file=paste0(outbase, ".RNA.intersections.byReadDecile.txt"), append=F, sep="\t")

cat("iPCR library overlap with other libraries\n")
iPCR.overlaps.byReadDecile <- cbind(sample=prefix, source="iPCR", summaryBy(RNA+DNA~iPCR.bin+iPCR.decile, FUN=pctTrue, data=subset(data, !is.na(iPCR.bin)), keep.names=T))
colnames(iPCR.overlaps.byReadDecile)[colnames(iPCR.overlaps.byReadDecile)=="iPCR.bin"] <- "reads.bin"
colnames(iPCR.overlaps.byReadDecile)[colnames(iPCR.overlaps.byReadDecile)=="iPCR.decile"] <- "reads.decile"
print.data.frame(iPCR.overlaps.byReadDecile, row.names=F, right=F)
write.table(melt.data.frame(iPCR.overlaps.byReadDecile), row.names=F, col.names=T, quote=F, file=paste0(outbase, ".iPCR.intersections.byReadDecile.txt"), append=F, sep="\t")



cat("\n\nAnalysis of mapped integration site\n")
sampleData <- read(paste0(outbase, ".txt"), header=T)

###Basic filtering
sampleData <- sampleData[!is.na(sampleData[,"DNA"]),]
sampleData <- sampleData[!is.na(sampleData[,"RNA"]),]
sampleData <- sampleData[!is.na(sampleData[,"iPCR"]),]

#NB throws away sites with the same coords. There's not a ton of these, and I assume they are odd sequencing errors, but some have high read counts
#throwing away sites at same position, opposite strand for now: sampleData[,"strand"]
sampleData <- sampleData[!duplicaterows(paste(sampleData[,"chrom"], sampleData[,"chromStart"])),]


###Thresholding
numDNAreads <- sum(sampleData[,"DNA"])
numRNAreads <- sum(sampleData[,"RNA"])


#This is redundant as barcodes.coords.bed is already thresholded
sampleData  <- subset(sampleData, iPCR >= 2)
#But DNA/RNA counts are not
sampleData  <- subset(sampleData, DNA >= 10)

#Normalize to counts per 1M reads
sampleData[,"DNA"] <- sampleData[,"DNA"] / numDNAreads * 10^6
sampleData[,"RNA"] <- sampleData[,"RNA"] / numRNAreads * 10^6

sampleData$expression <- sampleData[,"RNA"]/sampleData[,"DNA"]
sampleData$zscore <- scale(log(sampleData[,"expression"], base=2))


cat("\nCompare nearby insertions\n")
ins.proximity <- data.frame(subset(sampleData, select=c(chrom, chromStart, expression))[-nrow(sampleData),], subset(sampleData, select=c(chrom, chromStart, expression))[-1,], check.names=T)
ins.proximity <- subset(ins.proximity, chrom==chrom.1)
ins.proximity$dist <- ins.proximity$chromStart.1 - ins.proximity$chromStart
ins.proximity$dist.bin <- cut(ins.proximity$dist, breaks=c(0, 2.5e4, 5e4, 7.5e4, 10^5, Inf), right=F, include.lowest=T) #labels=c("50 bp", "100 bp", "250 bp", "500 bp", "1 kb", "10 kb", "100 kb", "1 Mb")

#png(file="snps.distances.png", width=750, height=750)
#hist(log(ins.proximity$dist, base=10), main="Distances between SNVs (union of all strains)")
#dev.off()

results.proxcor <- NULL
for(curDist in levels(ins.proximity$dist.bin)) {
	curdata <- subset(ins.proximity, dist.bin==curDist)
	if(nrow(curdata)>2) {
		curcor <- with(curdata, cor(expression, expression.1, use="na.or.complete", method="pearson"))
		#cat(curRsquared, ":", curcor, "\n")
		results.proxcor <- rbind(results.proxcor, data.frame(sample=prefix, dist.bin=curDist, meandist=mean.na(curdata$dist), cor=curcor, n=nrow(curdata)))
	}
}
print.data.frame(results.proxcor, row.names=F, right=F)
write.table(results.proxcor, row.names=F, col.names=T, quote=F, file=paste0(outbase, ".proxcor.txt"), append=F, sep="\t")


cat("\nOutput for ucsc\n")
cat("Clipping", length(which(sampleData[,"zscore"] > 3)), "sites z > 3\n")
sampleData[sampleData[,"zscore"] > 3,"zscore"] <- 3
cat("Clipping", length(which(sampleData[,"zscore"] < -3)), "sites z < -3\n")
sampleData[sampleData[,"zscore"] < -3,"zscore"] <- -3
sampleData[,"zscore"] <- (sampleData[,"zscore"] + 3) / 6

write.table(subset(sampleData, select=c("chrom", "chromStart", "chromEnd",  "zscore")), file=paste0(outbase, '.zscore.bed'), quote=F, sep='\t', col.names=F, row.names=F) 


cat("\nDone!!!\n")
print(date())
