#!/bin/env Rscript

print(date())



suppressPackageStartupMessages(library(directlabels))

old <- theme_set(theme_classic(base_size=7)) #pdf
old <- theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
old <- theme_update(panel.border = element_blank(), strip.background = element_blank()) #remove facet panel border


#Try to work around error "unable to start device PNG" on ISG cluster
options(bitmapType="cairo") 


cat("Loading insert lengths\n")
results <- NULL
for(f in list.files(path=".", pattern=".insertlengths.txt.gz$", recursive=T)) {
       cat("Doing", f, "\n")
       data <- read(f)
       colnames(data) <- c("sample", "length")
       dens <- density(subset(data, length<500 & length>0)$length)
       results <- rbind(results, data.frame(sample=data[1, "sample"], x=dens$x, y=dens$y, stringsAsFactors=F))
}

results$sample <- factor(gsub(".hg19$", "", results$sample))
results$sample <- factor(gsub(".hg38$", "", results$sample))
results$sample <- factor(gsub(".mm10$", "", results$sample))
results$DS <- sapply(as.character(results$sample), function(x) {unlist(strsplit(x, "-"))[2]})


cat("Saving...\n")
save(list=c("results"), file="insertlengths.RData", compress="bzip2")


p <- ggplot(data=subset(results, x>=27 & x<=300), aes(x=x, y=y, group=DS, label=DS, color=DS)) +
labs(title="") +
geom_line() +
xlab("Fragment length (bp)") +
ylab("Density") +
scale_x_continuous(breaks=c(27,36, 125, 200, 300), limits=c(27,300)) +
guides(color=F, linetype=F) +
geom_dl(method=list("top.points", cex=0.8))
#theme(plot.margin=unit(c(0,0,0,0), "lines"))) #NB top, rt, bot, left

pdf("insertlengths.pdf", width=7, height=5)
print(p)
dev.off()


#doesn't work much
old <- theme_set(theme_classic(base_size=10)) #png
png("insertlengths.png", width=700, height=500)
print(p)
dev.off()


print(date())
cat("\ndone\n")
