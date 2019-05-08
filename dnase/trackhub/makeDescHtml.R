#!/bin/env Rscript
# module load R/3.5.0

library(tableHTML)

args <- commandArgs(TRUE)
fname_in <- args[1]
fname_out <- args[2]

# Get the text version of the table
desc_text <- read.table(fname_in, sep="\t", header=TRUE, stringsAsFactors=FALSE)

delete_cols <- c("Genome", "Num.hotspots2", "SPOT2")
desc_text <- desc_text[,!colnames(desc_text) %in% delete_cols]
#for (i in delete_cols) {
#	if(i %in% colnames(desc_text)) {
#		desc_text[, i] <- NULL
#	}
#}

# Add commas to fields with large numeric values.
big_numeric_cols <- c("Num.sequenced.reads", "Num.pass.filter.alignments", "Num.mitochondria.reads", 
						"Num.analyzed.reads", "Num.duplicate.reads", "Genomic.coverage", "Num.hotspots")

for (i in big_numeric_cols) {
	if(i %in% colnames(desc_text)) {
		desc_text[,i] <- format(desc_text[,i], big.mark=",", trim=TRUE)
	}
}

# Get rid of underscores in header lines.
colnames(desc_text) <- chartr(old="_", new=" ", colnames(desc_text))

# Make the html version
desc_HTML <- tableHTML(desc_text, rownames=FALSE) %>% add_css_row(css = list('background-color', 'lightblue'), rows = odd(1:(1+nrow(desc_text)))) %>%
														add_css_row(css = list('background-color', 'white'), rows = even(1:(1+nrow(desc_text))))
write_tableHTML(desc_HTML, file=fname_out, complete_html=FALSE)

