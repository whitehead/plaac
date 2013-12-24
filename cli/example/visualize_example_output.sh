#!/usr/bin/env Rscript
source("../R/plaac_plot_util.r")
plot_seqs("output.tsv", "output.pdf", showSeq=TRUE)
plot_seqs("output.tsv", "output.png", showSeq=TRUE)
