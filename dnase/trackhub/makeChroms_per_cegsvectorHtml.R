#!/bin/env Rscript

library(tableHTML)   # Documentation:  https://cloud.r-project.org/web/packages/tableHTML/tableHTML.pdf

args <- commandArgs(TRUE)
fname_in <- args[1]
fname_out <- args[2]

# Get the text version of the table
desc_text <- read.table(fname_in, sep="|", header=FALSE, stringsAsFactors=FALSE)

colnames(desc_text) <- c("Reference", "Chromosomes")

# Make the html version
desc_HTML <- tableHTML(desc_text, rownames=FALSE, escape=FALSE) %>% add_css_row(css = list('background-color', 'lightblue'), rows = odd(1:(1+nrow(desc_text)))) %>%
														add_css_row(css = list('background-color', 'white'), rows = even(1:(1+nrow(desc_text))))
write_tableHTML(desc_HTML, file=fname_out, complete_html=FALSE)

