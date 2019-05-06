#!/bin/env Rscript

print(date())

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(directlabels))

old <- theme_set(theme_classic(base_size=7)) #pdf
old <- theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
old <- theme_update(panel.border = element_blank(), strip.background = element_blank()) #remove facet panel border


#Try to work around "unable to start device PNG" on ISG cluster
options(bitmapType="cairo") 

if(length(commandArgs(TRUE)) >= 1) {
	maxlength <- as.integer(commandArgs(TRUE)[1])
} else {
	maxlength <- 300
}
if(length(commandArgs(TRUE)) >= 2) {
	dir <- commandArgs(TRUE)[2]
} else {
	dir <- "."
}

if(maxlength <= 300) {
	breaks <- c(36, 50, 75, 100, 125, 150, 175, seq.int(200, maxlength, 50))
} else {
	breaks <- seq.int(50, maxlength, 50)
}

cat("Loading insert lengths upto", maxlength, "bp\n")
#TODO parallelize?
results <- NULL
for(d in c(list.dirs(path=dir, recursive=F, full.names=T))) {
	for(f in list.files(path=d, pattern=".insertlengths.txt.gz$", recursive=F, full.names=T)) {
		cat("Doing", f, "\n")
		data <- read(f)
		if(ncol(data)==1) {
			colnames(data) <- c("length")
		} else {
			#as of 2019-02-04, name no longer included to save space
			#Backwards compatible for now
			colnames(data) <- c("name", "length")
		}
		dens <- density(subset(data, length<=maxlength & length>=27)$length, bw=10)
		#can also get name from data[1, "name"]; though note mappedgenome was only included inside file from 2018-12-29 on
		results <- rbind(results, data.frame(name=gsub(".insertlengths.txt.gz$", "", basename(f)), x=dens$x, y=dens$y, stringsAsFactors=F))
	}
}

if(is.null(results)) {
	cat("There were no paired end reads.\n")
	cat("insertlengths.RData is not being produced.\n")
	print(date())
	cat("\ndone\n")
	quit(save = "no")
}

results$mappedgenome <- factor(gsub("^.+\\.(hg19|hg38|mm10|rn6)(_full|_noalt|_sacCer3)?$", "\\1\\2", results$name, perl=T))
results$sample <- factor(gsub(".(hg19|hg38|mm10|rn6)(_full|_noalt|_sacCer3)?$", "", results$name))
results$BS <- sapply(as.character(results$sample), function(x) {unlist(strsplit(x, "-"))[2]})
results$sublibrary <- sapply(results$BS, FUN=function(x) {substr(x, 8, 8)})


cat("Saving...\n")
save(list=c("results", "breaks", "maxlength"), file=paste0(dir, "/insertlengths.RData"), compress="bzip2")

#For subsetting
#results <- subset(results, BS %in% c("BS01403A", "BS01403B", "BS01403C", "BS01403D", "BS01409A", "BS01409B", "BS01409C", "BS01409D"))
#results <- subset(results, grepl("Hba_", name))
#results$cleanup <- sapply(results$sample, FUN=function(x) {unlist(strsplit(as.character(x), "_"))[5]})
#results$polymerase <- sapply(results$sample, FUN=function(x) {unlist(strsplit(as.character(x), "_"))[6]})


p <- ggplot(data=subset(results, x>=27 & x<=maxlength), aes(x=x, y=y, group=name, label=BS, color=BS)) +
labs(title="") +
geom_line() +
xlab("Fragment length (bp)") +
ylab("Density") +
scale_x_continuous(breaks=breaks, limits=c(27,maxlength)) +
scale_y_continuous(labels = function(x) {scales::scientific(x, digits=1)}) +
#facet_grid(polymerase~., scales = "free_x", space = "free_x") +
guides(color=F, linetype=F) +
theme(legend.position = c(.85, 0.9)) +
geom_dl(method=list("top.points", cex=0.8))
#theme(plot.margin=unit(c(0,0,0,0), "lines"))) #NB top, rt, bot, left

pdf(paste0(dir, "/insertlengths.pdf"), width=7, height=5)
print(p)
dev.off()


old <- theme_set(theme_classic(base_size=15)) #png
png(paste0(dir, "/insertlengths.png"), width=700, height=500)
print(p)
dev.off()


print(date())
cat("\ndone\n")
