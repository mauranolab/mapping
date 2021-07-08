#!/bin/env Rscript

print(date())


old <- theme_set(theme_classic(base_size=8)) #pdf
old <- theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
old <- theme_update(axis.text.x=element_text(size=8), axis.text.y=element_text(size=8))

#Try to work around "unable to start device PNG" on ISG cluster
options(bitmapType="cairo") 

argv <- commandArgs(TRUE)
prefix <- argv[1]
outbase <- argv[2]


#Basic filtering
filterSampleData <- function(data, prefix="merged", useDNA=TRUE) {
	#iPCR filters
	cat("Thresholding iPCR counts\n")
	#The 2-read cutoff is redundant as barcodes.coords.bed was already thresholded
	data  <- subset(data, iPCR >= 2 & !is.na(iPCR))
	
	for(curName in unique(data$Name)) {
		cat("Normalizing", curName, "\n")
		#DNA filters
		if(useDNA) {
			cat("Thresholding and normalizing DNA counts\n")
			#For normalizing to counts per 1M reads
			numDNAreads <- sum(data[data$Name==curName, "DNA"], na.rm=T)
			
			#DNA/RNA counts are not already thresholded
			cat("Removing", length(which(with(data[data$Name==curName,], DNA < 10 | is.na(DNA)))), "sites\n")
			data  <- subset(data, Name!=curName | (DNA >= 10 & !is.na(DNA)))
		
			data[data$Name==curName, "DNA"] <- data[data$Name==curName, "DNA"] / numDNAreads * 10^6
		}
		
		
		#RNA filters
		cat("Normalizing RNA counts\n")
		numRNAreads <- sum(data[data$Name==curName, "RNA"], na.rm=T)
		data[data$Name==curName, "RNA"] <- data[data$Name==curName, "RNA"] / numRNAreads * 10^6
	}
	
	
	#Merge sites that are repeated
	data$mergeIDs <- paste(data[,"chrom"], data[,"chromStart"])
	data.merged <- data[!duplicated(data$mergeIDs),]
	cat("Summarizing", nrow(data)-nrow(data.merged), "duplicate rows\n")
	numericCols <- c("DNA", "RNA", "iPCR") #columns to be averaged rather than copied from first row
	for(i in 1:nrow(data.merged)) {
		curID <- data.merged[i, "mergeIDs"]
		rowsToMerge <- data[,"mergeIDs"]==curID
		if(length(which(rowsToMerge))>1) {
#			cat("Summarizing ", curID, " (", length(which(rowsToMerge)) , " samples)\n", sep="")
			if(length(unique(data[rowsToMerge,"strand"]))!=1) {
				#TODO are averaging some sites at same position but with opposite strand for now (but don't see much of this): data[,"strand"]
				cat("WARNING merging BCs at", curID, "with opposite strand\n")
			}
			data.merged[i,numericCols] <- apply(data[rowsToMerge,numericCols], MARGIN=2, FUN=function(x) { mean(x, na.rm=T) })
		}
	}
	data.merged$Name <- prefix
	
	
	if(useDNA) {
		#Don't think zeroing NAs is so good for 10x data, at least at current coverage
		cat("Zeroing NAs at", length(which(is.na(data.merged[,"RNA"]))), "sites with no RNA reads\n")
		data.merged[is.na(data.merged[,"RNA"]), "RNA"] <- 0
	} else {
		cat("Removing", length(which(is.na(data.merged$RNA))), "sites\n")
		data.merged  <- subset(data.merged, !is.na(RNA))
	}
	
	
	if(useDNA) {
		data.merged$expression <- data.merged[,"RNA"]/data.merged[,"DNA"]
	} else {
		data.merged$expression <- data.merged[,"RNA"]
	}
	#Add pseudocount. The as.numeric strips attr that messes up further analysis
	#Not really a z-score
	data.merged$activity <- as.numeric(log(data.merged[,"expression"]+1, base=2))
#	data.merged$activity <- as.numeric(scale(log(data.merged[,"expression"]+1, base=2)))
	
	return(data.merged)
}


cat("\n\nAnalysis of mapped integration site\n")
sampleData <- read(paste0(outbase, ".mapped.txt"), header=T, nrows=20000)
colnames(sampleData)[colnames(sampleData)=="X.chrom"] <- "chrom"

cat("Loaded", nrow(sampleData), "sites\n")


if(all(is.na(sampleData$DNA))) {
	#For 10x data
	cat("No DNA counts provided\n")
	useDNA <- FALSE
} else {
	cat("Will normalize by DNA counts\n")
	useDNA <- TRUE
}


sampleData <- filterSampleData(sampleData, prefix, useDNA)
cat("Finished filtering.\n")
cat(prefix, "\tTotal sites remaining:\t", nrow(sampleData), "\n")


cat("\n\nCompare adjacent insertions\n")
ins.proximity <- data.frame(subset(sampleData, select=c(chrom, chromStart, activity))[-nrow(sampleData),], subset(sampleData, select=c(chrom, chromStart, activity))[-1,], check.names=T)
ins.proximity <- subset(ins.proximity, chrom==chrom.1)
ins.proximity$dist <- ins.proximity$chromStart.1 - ins.proximity$chromStart
ins.proximity$dist.bin <- cut(ins.proximity$dist, breaks=c(1, 100, 1000, 1e4, 1e5, Inf), labels=c("0+ bp", "100+ bp", "1+ kb", "10+ kb", "100+ kb"), right=F, include.lowest=T, pretty.labels="left")


results.proxcor <- NULL
for(curDist in levels(ins.proximity$dist.bin)) {
	curdata <- subset(ins.proximity, dist.bin==curDist)
	if(nrow(curdata)>2) {
		curcor <- with(curdata, cor(activity, activity.1, use="na.or.complete", method="spearman"))
		#cat(curRsquared, ":", curcor, "\n")
		results.proxcor <- rbind(results.proxcor, data.frame(sample=prefix, dist.bin=curDist, meandist=mean.na(curdata$dist), cor=curcor, meanactivity=mean(c(curdata$activity, curdata$activity.1)), n=nrow(curdata)))
	}
}
print.data.frame(results.proxcor, row.names=F, right=F)
write.table(results.proxcor, row.names=F, col.names=T, quote=F, file=paste0(outbase, ".proxcor.txt"), append=F, sep="\t")


cat("\nBy NumDHS100kb.quintile\n")
sampleData$NumDHS100kb.quintile <- factor(cut(sampleData$NumDHS100kb, breaks=unique(quantile(sampleData$NumDHS100kb, probs=seq.int(from=0, to=1, length.out=6), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,5,1))

ins.proximity.numDHS <- data.frame(subset(sampleData, select=c(chrom, chromStart, activity, NumDHS100kb.quintile))[-nrow(sampleData),], subset(sampleData, select=c(chrom, chromStart, activity, NumDHS100kb.quintile))[-1,], check.names=T)
ins.proximity.numDHS <- subset(ins.proximity.numDHS, chrom==chrom.1)
ins.proximity.numDHS$dist <- ins.proximity.numDHS$chromStart.1 - ins.proximity.numDHS$chromStart
ins.proximity.numDHS$dist.bin <- cut.pretty(ins.proximity.numDHS$dist, breaks=c(0, 1000, 1e4, 1e5, 5e5, Inf), right=F, include.lowest=T, pretty.labels="left") #labels=c("1 bp", "1 kbp", "10 kbp", "20 kbp", "30 kbp", "40 kbp", "50 kbp", "60 kbp", "70 kbp", "80 kbp", "90 kbp", "100 kbp")
ins.proximity.numDHS$compartments <- apply(ins.proximity.numDHS[,c("NumDHS100kb.quintile", "NumDHS100kb.quintile.1")], MARGIN=1, FUN=max)

results.proxcor.numDHS <- NULL
for(curDist in levels(ins.proximity.numDHS$dist.bin)) {
	for(curNumDHS in sort(unique(ins.proximity.numDHS$compartments))) {
		curdata <- subset(ins.proximity.numDHS, dist.bin==curDist & compartments==curNumDHS)
		if(nrow(curdata)>2) {
			curcor <- with(curdata, cor(activity, activity.1, use="na.or.complete", method="pearson"))
			#cat(curRsquared, ":", curcor, "\n")
			results.proxcor.numDHS <- rbind(results.proxcor.numDHS, data.frame(sample=prefix, dist.bin=curDist, numDHS.quintile=curNumDHS, meandist=mean.na(curdata$dist), cor=curcor, n=nrow(curdata)))
		}
	}
}
print.data.frame(results.proxcor.numDHS, row.names=F, right=F)
write.table(results.proxcor.numDHS, row.names=F, col.names=T, quote=F, file=paste0(outbase, ".proxcor.numDHS100kb.txt"), append=F, sep="\t")


cat("\nWhether in same TAD\n")
ins.proximity.TAD <- data.frame(subset(sampleData, select=c(chrom, chromStart, activity, TADid))[-nrow(sampleData),], subset(sampleData, select=c(chrom, chromStart, activity, TADid))[-1,], check.names=T)
ins.proximity.TAD <- subset(ins.proximity.TAD, chrom==chrom.1)
ins.proximity.TAD$dist <- ins.proximity.TAD$chromStart.1 - ins.proximity.TAD$chromStart
ins.proximity.TAD$dist.bin <- cut.pretty(ins.proximity.TAD$dist, breaks=c(1, 1000, 1e4, 1e5, Inf), right=F, include.lowest=T, pretty.labels="left") #labels=c("1 bp", "1 kbp", "10 kbp", "20 kbp", "30 kbp", "40 kbp", "50 kbp", "60 kbp", "70 kbp", "80 kbp", "90 kbp", "100 kbp")
ins.proximity.TAD$sameTAD <- ins.proximity.TAD$TADid==ins.proximity.TAD$TADid.1


results.proxcor.TAD <- NULL
for(curDist in levels(ins.proximity.TAD$dist.bin)) {
	for(curSameTAD in sort(unique(ins.proximity.TAD$sameTAD))) {
		curdata <- subset(ins.proximity.TAD, dist.bin==curDist & sameTAD==curSameTAD)
		if(nrow(curdata)>2) {
			curcor <- with(curdata, cor(activity, activity.1, use="na.or.complete", method="pearson"))
			#cat(curRsquared, ":", curcor, "\n")
			results.proxcor.TAD <- rbind(results.proxcor.TAD, data.frame(sample=prefix, dist.bin=curDist, sameTAD=curSameTAD, meandist=mean.na(curdata$dist), cor=curcor, meanactivity=mean(c(curdata$activity, curdata$activity.1)), n=nrow(curdata)))
		}
	}
}
print.data.frame(results.proxcor.TAD, row.names=F, right=F)
write.table(results.proxcor.TAD, row.names=F, col.names=T, quote=F, file=paste0(outbase, ".proxcor.TAD.txt"), append=F, sep="\t")


cat("\nlm\n")
fit <- lm(activity~scale(NumDHS5kb)+scale(NumDHS100kb), data=sampleData)
print(summary(fit))

sampleData$fit <- predict(fit, sampleData, type="response")
sampleData$residual <- with(sampleData, activity-fit)
cat("Clipping", length(which(sampleData[,"residual"] < -2)), "sites z < -2\n")
sampleData[sampleData[,"residual"] < -2, "residual"] <- -2
sampleData[,"residual"] <- (sampleData[,"residual"] + 2) / 4


cat("\nOutput for ucsc\n")
write.table(subset(sampleData, select=c("chrom", "chromStart", "chromEnd", "BC", "activity")), file=paste0(outbase, '.activity.bed'), quote=F, sep='\t', col.names=F, row.names=F) 


cat("\n\nActivity related to nearest TSS\n")
sampleData$DistToTSS.bin <- cut(abs(sampleData$DistToTSS), breaks=c(0, 500, 5000, 1e4, 5e4, 5e5, Inf), labels=c("0+ bp", "500+ bp", "5+ kb", "10+ kb", "50+ kb", "500+ kb"), right=F, include.lowest=T, pretty.labels="left")
sampleData$DistToTSS.decile <- factor(cut(abs(sampleData$DistToTSS), breaks=unique(quantile(abs(sampleData$DistToTSS), probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))
results.DistToTSS <- cbind(sample=prefix, summaryBy(activity~DistToTSS.bin, data=sampleData, FUN=list(mean, median, length)))
print.data.frame(results.DistToTSS, row.names=F, right=F)
write.table(results.DistToTSS, file=paste0(outbase, '.DistToTSS.txt'), quote=F, sep='\t', col.names=T, row.names=F) 


cat("\n\nActivity related to nearest DHS\n")
sampleData$DistToNearestDHS.bin <- cut(abs(sampleData$DistToNearestDHS), breaks=c(0, 500, 5000, 1e4, 5e4, 5e5, Inf), labels=c("0+ bp", "500+ bp", "5+ kb", "10+ kb", "50+ kb", "500+ kb"), right=F, include.lowest=T, pretty.labels="left")
#sampleData$DistToNearestDHS.decile <- factor(cut(abs(sampleData$DistToNearestDHS), breaks=unique(quantile(abs(sampleData$DistToNearestDHS), probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))
results.DistToNearestDHS <- cbind(sample=prefix, summaryBy(activity~DistToNearestDHS.bin, data=sampleData, FUN=list(mean, median, length)))
print.data.frame(results.DistToNearestDHS, row.names=F, right=F)
write.table(results.DistToNearestDHS, file=paste0(outbase, '.DistToNearestDHS.txt'), quote=F, sep='\t', col.names=T, row.names=F) 


cat("\n\nActivity related to nearest DHS (no CTCF)\n")
sampleData$DistToNearestDHSnoCTCF.bin <- cut(abs(sampleData$DistToNearestDHSnoCTCF), breaks=c(0, 500, 5000, 1e4, 5e4, 5e5, Inf), labels=c("0+ bp", "500+ bp", "5+ kb", "10+ kb", "50+ kb", "500+ kb"), right=F, include.lowest=T)
#sampleData$DistToNearestDHSnoCTCF.decile <- factor(cut(abs(sampleData$DistToNearestDHSnoCTCF), breaks=unique(quantile(abs(sampleData$DistToNearestDHSnoCTCF), probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))
results.DistToNearestDHSnoCTCF <- cbind(sample=prefix, summaryBy(activity~DistToNearestDHSnoCTCF.bin, data=sampleData, FUN=list(mean, median, length)))
print.data.frame(results.DistToNearestDHSnoCTCF, row.names=F, right=F)
write.table(results.DistToNearestDHSnoCTCF, file=paste0(outbase, '.DistToNearestDHSnoCTCF.txt'), quote=F, sep='\t', col.names=T, row.names=F) 


cat("\n\nActivity related to nearest CTCF site\n")
sampleData$DistToNearestCTCF.bin <- cut(abs(sampleData$DistToNearestCTCF), breaks=c(0, 500, 5000, 1e4, 5e4, 5e5, Inf), labels=c("0+ bp", "500+ bp", "5+ kb", "10+ kb", "50+ kb", "500+ kb"), right=F, include.lowest=T)
#sampleData$DistToNearestCTCF.decile <- factor(cut(abs(sampleData$DistToNearestCTCF), breaks=unique(quantile(abs(sampleData$DistToNearestCTCF), probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))
results.DistToNearestCTCF <- cbind(sample=prefix, summaryBy(activity~DistToNearestCTCF.bin, data=sampleData, FUN=list(mean, median, length)))
print.data.frame(results.DistToNearestCTCF, row.names=F, right=F)
write.table(results.DistToNearestCTCF, file=paste0(outbase, '.DistToNearestCTCF.txt'), quote=F, sep='\t', col.names=T, row.names=F) 


cat("\n\nActivity related to DHS density\n")
#using deciles means that bin endpoints are not the same across samples
#sampleData$NumDHS5kb.decile <- factor(cut(sampleData$NumDHS5kb, breaks=unique(quantile(sampleData$NumDHS5kb, probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))
sampleData$NumDHS5kb.bin <- cut.pretty(sampleData$NumDHS5kb, breaks=c(0,1,2,3,4,5,6,Inf), right=F, include.lowest=T, pretty.labels="left")
results.NumDHS5kb <- cbind(sample=prefix, summaryBy(activity+NumDHS5kb~NumDHS5kb.bin, data=sampleData, FUN=list(mean, median, length)))
print.data.frame(results.NumDHS5kb, row.names=F, right=F)
write.table(results.NumDHS5kb, file=paste0(outbase, '.NumDHS5kb.txt'), quote=F, sep='\t', col.names=T, row.names=F)


cat("\n\nActivity related to DHS density\n")
#using deciles means that bin endpoints are not the same across samples
#sampleData$NumDHS100kb.decile <- factor(cut(sampleData$NumDHS100kb, breaks=unique(quantile(sampleData$NumDHS100kb, probs=seq.int(from=0, to=1, length.out=11), na.rm=T)), right=F, include.lowest=T, labels=F), ordered=T, levels=seq.int(1,10,1))
sampleData$NumDHS100kb.bin <- cut.pretty(sampleData$NumDHS100kb, breaks=c(0,5,10,15,20,25,30,35,40,45,50,Inf), right=F, include.lowest=T, pretty.labels="left")
results.NumDHS100kb <- cbind(sample=prefix, summaryBy(activity+NumDHS100kb~NumDHS100kb.bin, data=sampleData, FUN=list(mean, median, length)))
print.data.frame(results.NumDHS100kb, row.names=F, right=F)
write.table(results.NumDHS100kb, file=paste0(outbase, '.NumDHS100kb.txt'), quote=F, sep='\t', col.names=T, row.names=F)


cat("\n\nOverlap between libraries by read depth\n")
#NB based on the AllBCs file
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
print.data.frame(cbind(prefix, summaryBy(activity~InsDens,data=sampleData, FUN=list(mean, length))), row.names=F, right=F)

cat("\nDone!!!\n")
print(date())
