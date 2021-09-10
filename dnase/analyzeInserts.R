#!/bin/env Rscript

print(date())

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(directlabels))

old <- theme_set(theme_classic(base_size=8)) #pdf
old <- theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
old <- theme_update(panel.border = element_blank(), strip.background = element_blank()) #remove facet panel border


#Try to work around "unable to start device PNG" on ISG cluster
options(bitmapType="cairo") 


#read() -- smarter wrapper around read.table
#BUGBUG who knows why i'm getting "zcat: stdout: Broken pipe" on first read.table
#TODO seems I need to use gzcat on the mac and zcat on linux
read <- function(filename, nrows=100, header=F, col.classes=NULL, stringsAsFactors=F, sep="\t", comment.char="", ...) {
	#BUGBUG should print example line of file upon failure
	#TODO can't use a pipe or specify which columns to use
	#BUGBUG when the first column is all characters, R tries to use it as row names. But this goofs up my hardcoded colClasses. row.names=NULL doesn't work in the function, though it works interactively.
	
	filename <- Sys.glob(filename)
	if(length(filename) == 0 || !file.exists(filename)) {
		stop("ERROR read() -- ", filename, " does not exist!\n")
	} else if(length(filename) > 1) {
		stop("ERROR read() -- Multiple files were matched: ",  paste(filename, collapse=", "), "!\n")
	}
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
		
			shortcopy <- read.table(pipe(readCmd), sep=sep, comment.char = comment.char, quote = "", strip.white = TRUE, stringsAsFactors = F, nrows=nrows, header=header, ...)
			guessed.classes <- lapply(shortcopy, typeof)
		} else {
			guessed.classes <- col.classes
		}
		
		return(read.table(pipe(readCmd), sep=sep, comment.char = comment.char, quote = "", strip.white = TRUE, stringsAsFactors = stringsAsFactors, nrows=num.lines, colClasses=guessed.classes, header=header, ...))
	}
}


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
gc_coverage <- NULL
for(d in c(list.dirs(path=dir, recursive=F, full.names=T))) {
	for(f in list.files(path=d, pattern=".insertlengths.txt.gz$", recursive=F, full.names=T)) {
		cat("Doing", f, "\n")
		data <- read(f)
		if(!is.null(data)) {
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
	## Read picardmetrics .gc_bias.detail_metrics
	for (f in list.files(file.path(d, "picardmetrics"), pattern = "(hg19|hg38|mm10|rn6|talOcc4).gc_bias.detail_metrics", recursive=F, full.names=T)) {
		cat("Reading", f, "\n")
		data <- tryCatch(read.delim(f, comment.char = "#", stringsAsFactors = FALSE), error = function(x) NULL)
		if (!is.null(data)) {
			data <- data[, c("GC", "WINDOWS", "READ_STARTS", "MEAN_BASE_QUALITY", "NORMALIZED_COVERAGE", "ERROR_BAR_WIDTH")]
			data$name <- gsub(".gc_bias.detail_metrics", "", basename(f))
			gc_coverage <- rbind(gc_coverage, data)
		}
	}
}

if(is.null(results)) {
	cat("There were no paired end reads.\n")
	cat("insertlengths.RData is not being produced.\n")
	print(date())
	cat("\ndone\n")
	quit(save="no", status=0)
}

results$mappedgenome <- factor(gsub("^.+\\.(hg19|hg38|mm10|rn6|talOcc4)(_full|_noalt|_sacCer3|_full_wuhCor1)?$", "\\1\\2", results$name, perl=T))
results$sample <- factor(gsub(".(hg19|hg38|mm10|rn6|talOcc4)(_full|_noalt|_sacCer3|_full_wuhCor1)?$", "", results$name))
results$BS <- sapply(as.character(results$sample), function(x) {unlist(strsplit(x, "-"))[2]})
results$sublibrary <- sapply(results$BS, FUN=function(x) {substr(x, 8, 8)})

gc_coverage$mappedgenome <- factor(gsub("^.+\\.(hg19|hg38|mm10|rn6|talOcc4)(_full|_noalt|_sacCer3|_full_wuhCor1)?$", "\\1\\2", gc_coverage$name, perl=T))
gc_coverage$sample <- factor(gsub(".(hg19|hg38|mm10|rn6|talOcc4)(_full|_noalt|_sacCer3|_full_wuhCor1)?$", "", gc_coverage$name))
gc_coverage$BS <- sapply(as.character(gc_coverage$sample), function(x) {unlist(strsplit(x, "-"))[2]})
gc_coverage$sublibrary <- sapply(gc_coverage$BS, FUN=function(x) {substr(x, 8, 8)})

cat("Saving...\n")
save(list=c("results", "breaks", "maxlength"), file=paste0(dir, "/insertlengths.RData"), compress="bzip2")
save(list=c("gc_coverage"), file=file.path(dir, "gc_bias.RData"), compress="bzip2")

#For subsetting
#results <- subset(results, BS %in% c("BS01403A", "BS01403B", "BS01403C", "BS01403D", "BS01409A", "BS01409B", "BS01409C", "BS01409D"))
#results <- subset(results, grepl("Hba_", name))
#results$cleanup <- sapply(results$sample, FUN=function(x) {unlist(strsplit(as.character(x), "_"))[5]})
#results$polymerase <- sapply(results$sample, FUN=function(x) {unlist(strsplit(as.character(x), "_"))[6]})
#results <- subset(results, !grepl("^WaterControl", name))

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

## Reset plot theme
theme_set(old)
## Plotting GC graph
for(genome in unique(gc_coverage$mappedgenome)) {
	p <- ggplot(data = subset(gc_coverage, mappedgenome == genome), aes(GC, NORMALIZED_COVERAGE, group = name, label = BS, color = BS)) +
	geom_hline(yintercept = 1, color = "darkgray") +
	geom_line() +
	geom_dl(method=list("maxvar.qp", cex=0.8)) +
	coord_fixed(30, xlim = c(-10, 100), ylim = c(0, 3)) +
	labs(x = "GC content (%)", y = "Normalized Coverage") +
	guides(linetype=F) +
	theme(legend.position = c(.85, 0.9))

	pdf(file.path(dir, sprintf("gc_bias.%s.pdf", genome)), width=7, height=5)
	print(p)
	dev.off()

	png(file.path(dir, sprintf("gc_bias.%s.png", genome)), width=700, height=500)
	print(p + theme_classic(base_size=15))
	dev.off()
}

print(date())
cat("\ndone\n")
