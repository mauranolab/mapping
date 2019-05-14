#!/bin/env Rscript
# module load R/3.5.0

library(tableHTML)

args <- commandArgs(TRUE)
fname_in <- args[1]
fname_out <- args[2]

# Get the text version of the table
desc_text <- read.table(fname_in, sep="\t", header=TRUE, stringsAsFactors=FALSE)

delete_cols <- c("Num.hotspots2", "SPOT2")
if("Genome" %in% colnames(desc_text) && !any(grepl("^cegsvectors_", desc_text$Genome))) {
	delete_cols <- c(delete_cols, "Genome")
} else if("Genome" %in% colnames(desc_text)) {
	desc_text$Genome <- sub('cegsvectors_', '', desc_text$Genome, fixed=TRUE)
}
desc_text <- desc_text[,!colnames(desc_text) %in% delete_cols]

# Add commas to numeric fields.
for (i in 1:ncol(desc_text)) {
	if(is.numeric(desc_text[,i])) {
		desc_text[,i] <- format(desc_text[,i], big.mark=",", trim=TRUE)
	}
}

# Get rid of underscores in header lines.
colnames(desc_text) <- chartr(old="_", new=" ", colnames(desc_text))

# Make the html version
desc_HTML <- tableHTML(desc_text, rownames=FALSE) %>% add_css_row(css = list('background-color', 'lightblue'), rows = odd(1:(1+nrow(desc_text)))) %>%
														add_css_row(css = list('background-color', 'white'), rows = even(1:(1+nrow(desc_text))))
write_tableHTML(desc_HTML, file=fname_out, complete_html=FALSE)

