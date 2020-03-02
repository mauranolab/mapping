#!/bin/env Rscript

print(date())


old <- theme_set(theme_classic(base_size=7)) #pdf
old <- theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Try to work around "unable to start device PNG" on ISG cluster
options(bitmapType="cairo") 

argv <- commandArgs(TRUE)
prefix <- argv[1]
outbase <- argv[2]


#Basic filtering
filterSampleData <- function(data, useDNA=TRUE) {
	#iPCR filters
	cat("Thresholding iPCR counts\n")
	#The 2-read cutoff is redundant as barcodes.coords.bed was already thresholded
	data  <- subset(data, iPCR >= 2 & !is.na(iPCR))
	
	#DNA filters
	if(useDNA) {
		cat("Thresholding and normalizing DNA counts\n")
		#For normalizing to counts per 1M reads
		numDNAreads <- sum(data[,"DNA"], na.rm=T)
	
		#DNA/RNA counts are not already thresholded
		cat("Removing", length(which(data$DNA < 10 | is.na(data$DNA))), "sites\n")
		data  <- subset(data, DNA >= 10 & !is.na(DNA))
	
		data[,"DNA"] <- data[,"DNA"] / numDNAreads * 10^6
	}
	
	
	#RNA filters
	cat("Normalizing RNA counts\n")
	numRNAreads <- sum(data[,"RNA"], na.rm=T)
	data[,"RNA"] <- data[,"RNA"] / numRNAreads * 10^6
	
	
	#throwing away sites at same position, opposite strand for now (but don't see much of this): data[,"strand"]
	#TODO properly merge these
	cat("Removing", length(which(duplicaterows(paste(data[,"chrom"], data[,"chromStart"])))), "duplicate sites\n")
	data <- data[!duplicaterows(paste(data[,"chrom"], data[,"chromStart"])),]
	
	
	if(useDNA) {
		#Don't think zeroing NAs is so good for 10x data, at least at current coverage
		cat("Zeroing NAs at", length(which(is.na(data[,"RNA"]))), "sites with no RNA reads\n")
		data[is.na(data[,"RNA"]), "RNA"] <- 0
	} else {
		cat("Removing", length(which(is.na(data$RNA))), "sites\n")
		data  <- subset(data, !is.na(RNA))
	}
	
	
	if(useDNA) {
		data$expression <- (data[,"RNA"])/data[,"DNA"]
	} else {
		data$expression <- data[,"RNA"]
	}
	#Add pseudocount. The as.numeric strips attr that messes up further analysis
	data$zscore <- as.numeric(scale(log(data[,"expression"]+1, base=2)))
	
	return(data)
}


cat("\n\nAnalysis of mapped integration site\n")
sampleData <- read(paste0(outbase, ".mapped.txt"), header=T, nrows=20000)
colnames(sampleData)[colnames(sampleData)=="X.chrom"] <- "chrom"

cat("Loaded", nrow(sampleData), "sites\n")


if(all(is.na(sampleData$DNA))) {
	cat("No DNA counts provided\n")
	useDNA <- FALSE
} else {
	cat("Will normalize by DNA counts\n")
	useDNA <- TRUE
}


sampleData <- filterSampleData(sampleData, useDNA)
cat("Finished filtering.\n")
cat(prefix, "\tTotal sites remaining:\t", nrow(sampleData), "\n")


cat("\n\nCompare adjacent insertions\n")
ins.proximity <- data.frame(subset(sampleData, select=c(chrom, chromStart, zscore))[-nrow(sampleData),], subset(sampleData, select=c(chrom, chromStart, zscore))[-1,], check.names=T)
ins.proximity <- subset(ins.proximity, chrom==chrom.1)
ins.proximity$dist <- ins.proximity$chromStart.1 - ins.proximity$chromStart
ins.proximity$dist.bin <- cut.pretty(ins.proximity$dist, breaks=c(1, 1000, 1e4, 1e5, 5e5, Inf), right=F, include.lowest=T, pretty.labels="left") #labels=c("1 bp", "1 kbp", "10 kbp", "20 kbp", "30 kbp", "40 kbp", "50 kbp", "60 kbp", "70 kbp", "80 kbp", "90 kbp", "100 kbp")
#ins.proximity$dist.bin <- cut(ins.proximity$dist, breaks=c(0, 2e4, 4e4, 6e4, 8e4, 10^5, Inf), right=F, include.lowest=T) #labels=c("50 bp", "100 bp", "250 bp", "500 bp", "1 kb", "10 kb", "100 kb", "1 Mb")


results.proxcor <- NULL
for(curDist in levels(ins.proximity$dist.bin)) {
	curdata <- subset(ins.proximity, dist.bin==curDist)
	if(nrow(curdata)>2) {
		curcor <- with(curdata, cor(zscore, zscore.1, use="na.or.complete", method="spearman"))
		#cat(curRsquared, ":", curcor, "\n")
		results.proxcor <- rbind(results.proxcor, data.frame(sample=prefix, dist.bin=curDist, meandist=mean.na(curdata$dist), cor=curcor, n=nrow(curdata)))
	}
}
print.data.frame(results.proxcor, row.names=F, right=F)
write.table(results.proxcor, row.names=F, col.names=T, quote=F, file=paste0(outbase, ".proxcor.txt"), append=F, sep="\t")


#By NumDHS100kb.decile
sampleData$NumDHS100kb.decile <- factor(cut(sampleData$NumDHS100kb, breaks=unique(quantile(sampleData$NumDHS100kb, probs=seq.int(from=0, to=1, length.out=6), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,5,1))

ins.proximity.numDHS <- data.frame(subset(sampleData, select=c(chrom, chromStart, zscore, NumDHS100kb.decile))[-nrow(sampleData),], subset(sampleData, select=c(chrom, chromStart, zscore, NumDHS100kb.decile))[-1,], check.names=T)
ins.proximity.numDHS <- subset(ins.proximity.numDHS, chrom==chrom.1)
ins.proximity.numDHS$dist <- ins.proximity.numDHS$chromStart.1 - ins.proximity.numDHS$chromStart
ins.proximity.numDHS$dist.bin <- cut.pretty(ins.proximity.numDHS$dist, breaks=c(0, 1000, 1e4, 1e5, 5e5, Inf), right=F, include.lowest=T, pretty.labels="left") #labels=c("1 bp", "1 kbp", "10 kbp", "20 kbp", "30 kbp", "40 kbp", "50 kbp", "60 kbp", "70 kbp", "80 kbp", "90 kbp", "100 kbp")
ins.proximity.numDHS$compartments <- apply(ins.proximity.numDHS[,c("NumDHS100kb.decile", "NumDHS100kb.decile.1")], MARGIN=1, FUN=max)

results.proxcor.numDHS <- NULL
for(curDist in levels(ins.proximity.numDHS$dist.bin)) {
	for(curNumDHS in sort(unique(ins.proximity.numDHS$compartments))) {
		curdata <- subset(ins.proximity.numDHS, dist.bin==curDist & compartments==curNumDHS)
		if(nrow(curdata)>2) {
			curcor <- with(curdata, cor(zscore, zscore.1, use="na.or.complete", method="pearson"))
			#cat(curRsquared, ":", curcor, "\n")
			results.proxcor.numDHS <- rbind(results.proxcor.numDHS, data.frame(sample=prefix, dist.bin=curDist, numDHS.decile=curNumDHS, meandist=mean.na(curdata$dist), cor=curcor, n=nrow(curdata)))
		}
	}
}
print.data.frame(results.proxcor.numDHS, row.names=F, right=F)
write.table(results.proxcor.numDHS, row.names=F, col.names=T, quote=F, file=paste0(outbase, ".proxcor.numDHS.txt"), append=F, sep="\t")


#Whether in same TAD
ins.proximity.TAD <- data.frame(subset(sampleData, select=c(chrom, chromStart, zscore, TADid))[-nrow(sampleData),], subset(sampleData, select=c(chrom, chromStart, zscore, TADid))[-1,], check.names=T)
ins.proximity.TAD <- subset(ins.proximity.TAD, chrom==chrom.1)
ins.proximity.TAD$dist <- ins.proximity.TAD$chromStart.1 - ins.proximity.TAD$chromStart
ins.proximity.TAD$dist.bin <- cut.pretty(ins.proximity.TAD$dist, breaks=c(1, 1000, 1e4, 1e5, 5e5, Inf), right=F, include.lowest=T, pretty.labels="left") #labels=c("1 bp", "1 kbp", "10 kbp", "20 kbp", "30 kbp", "40 kbp", "50 kbp", "60 kbp", "70 kbp", "80 kbp", "90 kbp", "100 kbp")
ins.proximity.TAD$sameTAD <- ins.proximity.TAD$TADid==ins.proximity.TAD$TADid.1


results.proxcor.TAD <- NULL
for(curDist in levels(ins.proximity.TAD$dist.bin)) {
	for(curSameTAD in sort(unique(ins.proximity.TAD$sameTAD))) {
		curdata <- subset(ins.proximity.TAD, dist.bin==curDist & sameTAD==curSameTAD)
		if(nrow(curdata)>2) {
			curcor <- with(curdata, cor(zscore, zscore.1, use="na.or.complete", method="pearson"))
			#cat(curRsquared, ":", curcor, "\n")
			results.proxcor.TAD <- rbind(results.proxcor.TAD, data.frame(sample=prefix, dist.bin=curDist, sameTAD=curSameTAD, meandist=mean.na(curdata$dist), cor=curcor, n=nrow(curdata)))
		}
	}
}
print.data.frame(results.proxcor.TAD, row.names=F, right=F)
write.table(results.proxcor.TAD, row.names=F, col.names=T, quote=F, file=paste0(outbase, ".proxcor.TAD.txt"), append=F, sep="\t")


cat("\nlm\n")
#log(abs(DistToTSS)+1)+log(abs(DistToNearestDHSnoCTCF)+1)+log(abs(DistToNearestCTCF)+1)+
fit <- lm(zscore~NumDHS100kb, data=sampleData)
print(summary(fit))

sampleData$fit <- predict(fit, sampleData, type="response")
sampleData$residual <- with(sampleData, zscore-fit)
cat("Clipping", length(which(sampleData[,"residual"] < -2)), "sites z < -2\n")
sampleData[sampleData[,"residual"] < -2, "residual"] <- -2
sampleData[,"residual"] <- (sampleData[,"residual"] + 2) / 4

write.table(subset(sampleData, select=c(chrom, chromStart, chromEnd, BC, residual)), file=paste0(outbase, '.residual.bed'), quote=F, sep='\t', col.names=F, row.names=F) 
#bedops --range 25000 -m T0190_pMH034_BC4_pTR_C1fw_GGlo_PuroP2AGFP_A1rv_HS2.zscore.fithigh.bed > /tmp/fithigh.bed
#awk -F "\t" 'BEGIN {OFS="\t"} {$5= $5-$6; print}' T0190_pMH034_BC4_pTR_C1fw_GGlo_PuroP2AGFP_A1rv_HS2.zscore.fithigh.bed | cut -f1-5 | bedmap --delim "\t" --echo --count --mean /tmp/fithigh.bed - | awk -F "\t" 'BEGIN {OFS="\t"} $4>1' | widen 250000

#bwplot(log(DistToNearestCTCF, base=10)~factor(residual>2), data=sampleData)
#smoothScatter(sampleData$zscore, sampleData$fit)
#xyplot(zscore~fit, data=sampleData)
#with(subset(sampleData, ), smoothScatter(log(abs(DistToNearestCTCF)),residual))


cat("\nOutput for ucsc\n")
#cat("Clipping", length(which(sampleData[,"zscore"] > 3)), "sites z > 3\n")
#sampleData[sampleData[,"zscore"] > 3,"zscore"] <- 3
cat("Clipping", length(which(sampleData[,"zscore"] < -2)), "sites z < -2\n")
sampleData[sampleData[,"zscore"] < -3,"zscore"] <- -2
sampleData[,"zscore"] <- (sampleData[,"zscore"] + 2) / 5

write.table(subset(sampleData, select=c("chrom", "chromStart", "chromEnd", "BC", "zscore")), file=paste0(outbase, '.zscore.bed'), quote=F, sep='\t', col.names=F, row.names=F) 


cat("\n\nActivity related to nearest TSS\n")
sampleData$DistToTSS.bin <- cut.pretty(abs(sampleData$DistToTSS), breaks=c(0, 1000, 1e4, 1e5, 5e5, Inf), right=F, include.lowest=T, pretty.labels="left")
sampleData$DistToTSS.decile <- factor(cut(abs(sampleData$DistToTSS), breaks=unique(quantile(abs(sampleData$DistToTSS), probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))
results.DistToTSS <- cbind(sample=prefix, summaryBy(zscore~DistToTSS.bin, data=sampleData, FUN=list(mean, median, length)))
print.data.frame(results.DistToTSS, row.names=F, right=F)
write.table(results.DistToTSS, file=paste0(outbase, '.DistToTSS.txt'), quote=F, sep='\t', col.names=T, row.names=F) 


cat("\n\nActivity related to nearest DHS\n")
sampleData$DistToNearestDHS.bin <- cut.pretty(abs(sampleData$DistToNearestDHS), breaks=c(0, 1000, 1e4, 1e5, 5e5, Inf), right=F, include.lowest=T, pretty.labels="left")
#sampleData$DistToNearestDHS.decile <- factor(cut(abs(sampleData$DistToNearestDHS), breaks=unique(quantile(abs(sampleData$DistToNearestDHS), probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))
results.DistToNearestDHS <- cbind(sample=prefix, summaryBy(zscore~DistToNearestDHS.bin, data=sampleData, FUN=list(mean, median, length)))
print.data.frame(results.DistToNearestDHS, row.names=F, right=F)
write.table(results.DistToNearestDHS, file=paste0(outbase, '.DistToNearestDHS.txt'), quote=F, sep='\t', col.names=T, row.names=F) 
#summaryBy(DistToNearestDHS~DistToNearestDHS.decile, FUN=function(x) {mean(abs(x))}, data=sampleData)


cat("\n\nActivity related to nearest DHS (no CTCF)\n")
sampleData$DistToNearestDHSnoCTCF.bin <- cut(abs(sampleData$DistToNearestDHSnoCTCF), breaks=c(0, 1000, 1e4, 1e5, 5e5, Inf), labels=c("0 bp", "1 kb", "10 kb", "100 kb", "500 kb"), right=F, include.lowest=T)
#sampleData$DistToNearestDHSnoCTCF.bin <- cut.pretty(abs(sampleData$DistToNearestDHSnoCTCF), breaks=c(0, 1000, 1e4, 1e5, 5e5, Inf), right=F, include.lowest=T, pretty.labels="left")
#sampleData$DistToNearestDHSnoCTCF.decile <- factor(cut(abs(sampleData$DistToNearestDHSnoCTCF), breaks=unique(quantile(abs(sampleData$DistToNearestDHSnoCTCF), probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))
results.DistToNearestDHSnoCTCF <- cbind(sample=prefix, summaryBy(zscore~DistToNearestDHSnoCTCF.bin, data=sampleData, FUN=list(mean, median, length)))
print.data.frame(results.DistToNearestDHSnoCTCF, row.names=F, right=F)
write.table(results.DistToNearestDHSnoCTCF, file=paste0(outbase, '.DistToNearestDHSnoCTCF.txt'), quote=F, sep='\t', col.names=T, row.names=F) 


cat("\n\nActivity related to nearest CTCF site\n")
sampleData$DistToNearestCTCF.bin <- cut(abs(sampleData$DistToNearestCTCF), breaks=c(0, 1000, 1e4, 1e5, 5e5, Inf), labels=c("0 bp", "1 kb", "10 kb", "100 kb", "500 kb"), right=F, include.lowest=T)
#sampleData$DistToNearestCTCF.bin <- cut.pretty(abs(sampleData$DistToNearestCTCF), breaks=c(0, 1000, 1e4, 1e5, 5e5, Inf), right=F, include.lowest=T, pretty.labels="left")
#sampleData$DistToNearestCTCF.bin <- cut(abs(sampleData$DistToNearestCTCF), breaks=c(0, 1e3, 5e3, 1e4,  1e5, Inf), right=F, include.lowest=T)
#sampleData$DistToNearestCTCF.decile <- factor(cut(abs(sampleData$DistToNearestCTCF), breaks=unique(quantile(abs(sampleData$DistToNearestCTCF), probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))
results.DistToNearestCTCF <- cbind(sample=prefix, summaryBy(zscore~DistToNearestCTCF.bin, data=sampleData, FUN=list(mean, median, length)))
print.data.frame(results.DistToNearestCTCF, row.names=F, right=F)
write.table(results.DistToNearestCTCF, file=paste0(outbase, '.DistToNearestCTCF.txt'), quote=F, sep='\t', col.names=T, row.names=F) 


cat("\n\nActivity related to DHS density\n")
sampleData$NumDHS100kb.decile <- factor(cut(sampleData$NumDHS100kb, breaks=unique(quantile(sampleData$NumDHS100kb, probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))
results.NumDHS100kb <- cbind(sample=prefix, summaryBy(zscore~NumDHS100kb.decile, data=sampleData, FUN=list(mean, median, length)))
print.data.frame(results.NumDHS100kb, row.names=F, right=F)
write.table(results.NumDHS100kb, file=paste0(outbase, '.NumDHS100kb.txt'), quote=F, sep='\t', col.names=T, row.names=F) 



cat("\n\nOverlap between libraries by read depth\n")
#Note based on the AllBCs file
data <- read(paste0(outbase, ".AllBCs.txt"), header=T, nrows=5000)

if(useDNA) {
	data$DNA.bin <- cut.pretty(data[,"DNA"], breaks=unique(quantile(data[,"DNA"], probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, pretty.labels="left")
	data$DNA.decile <- factor(cut(data[,"DNA"], breaks=unique(quantile(data[,"DNA"], probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))

	cat("DNA library overlap with other libraries\n")
	DNA.overlaps.byReadDecile <- cbind(sample=prefix, source="DNA", summaryBy(RNA+iPCR~DNA.bin+DNA.decile, FUN=pctTrue, data=subset(data, !is.na(DNA.bin)), keep.names=T))
	colnames(DNA.overlaps.byReadDecile)[colnames(DNA.overlaps.byReadDecile)=="DNA.bin"] <- "reads.bin"
	colnames(DNA.overlaps.byReadDecile)[colnames(DNA.overlaps.byReadDecile)=="DNA.decile"] <- "reads.decile"
	print.data.frame(DNA.overlaps.byReadDecile, row.names=F, right=F)
	DNA.overlaps.byReadDecile.long <- melt.data.frame(DNA.overlaps.byReadDecile)
	colnames(DNA.overlaps.byReadDecile.long)[colnames(DNA.overlaps.byReadDecile.long)=="variable"] <- "target"
	write.table(DNA.overlaps.byReadDecile.long, row.names=F, col.names=T, quote=F, file=paste0(outbase, ".DNA.intersections.byReadDecile.txt"), append=F, sep="\t")
}


data$RNA.bin <- cut.pretty(data[,"RNA"], breaks=unique(quantile(data[,"RNA"], probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, pretty.labels="left")
data$RNA.decile <- factor(cut(data[,"RNA"], breaks=unique(quantile(data[,"RNA"], probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))

cat("RNA library overlap with other libraries\n")
RNA.overlaps.byReadDecile <- cbind(sample=prefix, source="RNA", summaryBy(DNA+iPCR~RNA.bin+RNA.decile, FUN=pctTrue, data=subset(data, !is.na(RNA.bin)), keep.names=T))
colnames(RNA.overlaps.byReadDecile)[colnames(RNA.overlaps.byReadDecile)=="RNA.bin"] <- "reads.bin"
colnames(RNA.overlaps.byReadDecile)[colnames(RNA.overlaps.byReadDecile)=="RNA.decile"] <- "reads.decile"
print.data.frame(RNA.overlaps.byReadDecile, row.names=F, right=F)
RNA.overlaps.byReadDecile.long <- melt.data.frame(RNA.overlaps.byReadDecile)
colnames(RNA.overlaps.byReadDecile.long)[colnames(RNA.overlaps.byReadDecile.long)=="variable"] <- "target"
if(!useDNA) {
	#easier to drop afterwards rather than make two formulae above
	RNA.overlaps.byReadDecile.long <- RNA.overlaps.byReadDecile.long[RNA.overlaps.byReadDecile.long$target!="DNA",]
}
write.table(RNA.overlaps.byReadDecile.long, row.names=F, col.names=T, quote=F, file=paste0(outbase, ".RNA.intersections.byReadDecile.txt"), append=F, sep="\t")


data$iPCR.bin <- cut.pretty(data[,"iPCR"], breaks=unique(quantile(data[,"iPCR"], probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, pretty.labels="left")
data$iPCR.decile <- factor(cut(data[,"iPCR"], breaks=unique(quantile(data[,"iPCR"], probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))

cat("iPCR library overlap with other libraries\n")
iPCR.overlaps.byReadDecile <- cbind(sample=prefix, source="iPCR", summaryBy(RNA+DNA~iPCR.bin+iPCR.decile, FUN=pctTrue, data=subset(data, !is.na(iPCR.bin)), keep.names=T))
colnames(iPCR.overlaps.byReadDecile)[colnames(iPCR.overlaps.byReadDecile)=="iPCR.bin"] <- "reads.bin"
colnames(iPCR.overlaps.byReadDecile)[colnames(iPCR.overlaps.byReadDecile)=="iPCR.decile"] <- "reads.decile"
print.data.frame(iPCR.overlaps.byReadDecile, row.names=F, right=F)
iPCR.overlaps.byReadDecile.long <- melt.data.frame(iPCR.overlaps.byReadDecile)
colnames(iPCR.overlaps.byReadDecile.long)[colnames(iPCR.overlaps.byReadDecile.long)=="variable"] <- "target"
if(!useDNA) {
	#easier to drop afterwards rather than make two formulae above
	RNA.overlaps.byReadDecile.long <- RNA.overlaps.byReadDecile.long[RNA.overlaps.byReadDecile.long$target!="DNA",]
}
write.table(iPCR.overlaps.byReadDecile.long, row.names=F, col.names=T, quote=F, file=paste0(outbase, ".iPCR.intersections.byReadDecile.txt"), append=F, sep="\t")


cat("iPCR insertion density vs. activity\n")
print.data.frame(cbind(prefix, summaryBy(zscore~InsDens,data=sampleData, FUN=list(mean, length))), row.names=F, right=F)

cat("\nDone!!!\n")
print(date())
