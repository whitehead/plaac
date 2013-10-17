#!/usr/bin/env Rscript
source("../speweyplot.r")

plotgenes(infile="output.tsv",
          outfile="output.pdf",
          seq=TRUE,
          use_pdf=TRUE
          )

plotgenes(infile="output.tsv",
          outfile="output.png",
          seq=TRUE,
          use_pdf=FALSE
          )
