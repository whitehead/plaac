#!/usr/bin/env Rscript
source("../R/plaac_plot_util.r")
plot_genes("output.tsv", "output.pdf", showSeq=TRUE)
plot_genes("output.tsv", "output.png", showSeq=TRUE)
