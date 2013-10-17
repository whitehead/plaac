#!/usr/bin/env Rscript
source("plaac_plot.r")

plot_genes(infile="prion_candidates_details.tsv",
          outfile="prion_details.pdf",
          showSeq=TRUE
          )

plot_genes(infile="prion_candidates_details.tsv",
          outfile="prion_details.png",
          showSeq=TRUE
          )
