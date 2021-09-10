#!/bin/env Rscript

print(date())

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(directlabels))

old <- theme_set(theme_classic(base_size=8)) #pdf
old <- theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
old <- theme_update(panel.border = element_blank(), strip.background = element_blank()) #remove facet panel border

#Try to work around "unable to start device PNG" on ISG cluster
options(bitmapType="cairo")

arg <- commandArgs(TRUE)
dir <- "."
min_cov <- 200
if (length(arg) >= 1) {
	dir <- arg[1]
}
if (length(arg) >= 2) {
	min_cov <- as.integer(arg[2])
}

## Read picardmetrics' .gc_bias.detail_metrics
gc_bias <- NULL
for(d in c(list.dirs(path = dir, recursive = F, full.names = T))) {
	for (f in list.files(file.path(d, "picardmetrics"), pattern = "(hg19|hg38|mm10|rn6|talOcc4).gc_bias.detail_metrics", recursive=F, full.names=T)) {
		cat("Reading", f, "\n")
		data <- tryCatch(read.delim(f, comment.char = "#", stringsAsFactors = FALSE), error = function(x) NULL)
		if (!is.null(data)) {
			data <- data[, c("GC", "WINDOWS", "READ_STARTS", "MEAN_BASE_QUALITY", "NORMALIZED_COVERAGE", "ERROR_BAR_WIDTH")]
			data$name <- gsub(".gc_bias.detail_metrics", "", basename(f))
			gc_bias <- rbind(gc_bias, data)
		}
	}
}

## Early exit if no file found
if (is.null(gc_bias)) {
	cat("No `.gc_bias.detail_metrics` files were found.\n")
	cat("\ndone\n")
	quit(save = "no", status = 0)
}

## Clean gc_bias and extract fields
gc_bias$mappedgenome <- factor(gsub("^.+\\.(hg19|hg38|mm10|rn6|talOcc4)(_full|_noalt|_sacCer3|_full_wuhCor1)?$", "\\1\\2", gc_bias$name, perl=T))
gc_bias$sample <- factor(gsub(".(hg19|hg38|mm10|rn6|talOcc4)(_full|_noalt|_sacCer3|_full_wuhCor1)?$", "", gc_bias$name))
gc_bias$BS <- sapply(as.character(gc_bias$sample), function(x) { unlist(strsplit(x, "-"))[2] })
gc_bias$sublibrary <- sapply(gc_bias$BS, FUN = function(x) { substr(x, 8, 8) })

## Saving data.frame
cat("Saving...\n")
save(list = c("gc_bias"), file = file.path(dir, "gc_bias.RData"), compress = "bzip2")

## Plotting GC graph
for(genome in unique(gc_bias$mappedgenome)) {
	cat("Plotting", genome, "GC bias...\n")
	d1 <- subset(gc_bias, mappedgenome == genome)
	d2 <- subset(d1, READ_STARTS >= min_cov)
	lost_bs <- setdiff(unique(d1$BS), unique(d2$BS))
	if (length(lost_bs) > 0) {
		cat("Samples lost due to all bin coverage below", min_cov, "reads :", paste(lost_bs, collapse = ","), "\n")
	}

	p <- ggplot(data = d2, aes(GC, NORMALIZED_COVERAGE, group = name, label = BS, color = BS)) +
	geom_hline(yintercept = 1, color = "darkgray") +
	geom_line() +
	geom_dl(method=list("maxvar.points", cex = 0.8)) +
	coord_fixed(30, xlim = c(-10, 100), ylim = c(0, 3)) +
	labs(x = "GC content (%)", y = "Normalized Coverage") +
	guides(color = F, linetype = F) +
	theme(legend.position = c(.85, 0.9))

	ggsave(file.path(dir, sprintf("gc_bias.%s.pdf", genome)), p, width = 7, height = 5)
	ggsave(file.path(dir, sprintf("gc_bias.%s.png", genome)), p + theme_classic(base_size = 15), width = 7, height = 5, dpi = 100)
}

print(date())
cat("\ndone\n")