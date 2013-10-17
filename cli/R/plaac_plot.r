#!/usr/bin/env Rscript
## Original program copyright 2009 Whitehead Institute for Biomedial Research
## Additions copyright 2013 University of Massachusetts
## Author: Oliver D. King (oliver.king@umassmed.edu) 
## See LICENSE.TXT for license information.  

args = commandArgs(TRUE)
source("plaac_plot_util.r") 
if (length(args) < 2) {
  cat("usage: Rscript plaac_plot.r plotdata.txt figname [-d]\n")
  cat("   where plotdata.txt is the output of the Java program plaac with the -p option,\n")
  cat("   figname is the name to save the figure as, ending with .pdf or .png.\n")
  cat("   (For png, one file will be output per protein, with a five digit index included in name.)\n")
  cat("   The optional flag -d plots doubly-smoothed rather than singly-smoothed curves.\n")
  cat("See sourcecode of plaac_plot_util.r and plaac_plot.r for other plot options.\n")
} else if (length(args)>2 && args[3]=="-d") {
  plot_genes(args[1], args[2], showSeq=T, showHMMProbs=T, showParses=T, 
             seqUnderline=c("PLAACx2"), seqColor="MAP",
             tracks=c("FIx2","PLAACx2","PAPAx2","THRESH"), w=41)
} else {
  plot_genes(args[1], args[2], showSeq=T, showHMMProbs=T, showParses=F, 
             seqUnderline=c("FIx2","PAPAx2"), seqColor="MAP",
             tracks=c("FI","PLAAC","PAPA","THRESH","HYDRO","CHARGE"), w=21)
}
