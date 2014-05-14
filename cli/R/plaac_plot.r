#!/usr/bin/env Rscript
## Original program copyright 2009 Whitehead Institute for Biomedical Research
## Additions copyright 2013 University of Massachusetts
## Author: Oliver D. King (oliver.king@umassmed.edu) 
## See LICENSE.TXT for license information.  

args = commandArgs(TRUE)
source("plaac_plot_util.r")
if (length(args) < 2) {
  cat("usage: Rscript plaac_plot.r plotdata.txt figname [-d|-c]\n")
  cat("   where plotdata.txt is the output of the Java program plaac with the -p option,\n")
  cat("   figname is the name to save the figure as, ending with .pdf or .png.\n")
  cat("   (For png, one file will be output per protein, with a five digit index included in name.)\n")
  cat("   The optional flag -d plots doubly-smoothed rather than singly-smoothed curves.\n")
  cat("   The optional flag -c plots all sequences con one page, color-coded by AA type.\n")
  cat("   When -c is used, optional integer afterward gives maximum number of AAs to plot per protein.\n")
  cat("See source-code of plaac_plot_util.r and plaac_plot.r for other plot options.\n")
} else if (length(args)>2 && args[3]=="-c") {
  max_n = ifelse(length(args) > 3, as.numeric(args[4]), NA) 
  color_code_seqs(args[1], args[2], max_n=max_n, showLegend=T, showParses=T)
} else if (length(args)>2 && args[3]=="-d") {
  plot_seqs(args[1], args[2], showSeq=T, showHMMProbs=T, showParses=T, showAAColors=T,
             seqUnderline=c("PLAACx2"), seqColor="MAP",
             tracks=c("FIx2","PLAACx2","PAPAx2","THRESH"))
}  else {
  plot_seqs(args[1], args[2], showSeq=T, showHMMProbs=T, showParses=T, showAAColors=T,
             seqUnderline=c("FIx2","PAPAx2"), seqColor="MAP",
             #tracks=c("FI","PLAAC","PAPA","THRESH","HYDRO","CHARGE"), 
             tracks=c("FI","PLAAC","PAPA","THRESH"))
}
