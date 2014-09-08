## Original program copyright 2009 Whitehead Institute for Biomedical Research
## Additions copyright 2013 University of Massachusetts
## Author: Oliver D. King (oliver.king@umassmed.edu) 
## See LICENSE.TXT for license information.  

## split long fasta IDs accross lines in plots
hyphenate = function(myString="--", targetWidth=90) {
  if (nchar(myString) < targetWidth + 5) {return(myString);}
  myChars = strsplit(myString,"")[[1]]
  bestsplit = which(myChars == " ")
  nextbest = which(myChars == ":" | myChars == "-" | myChars == "." | myChars == ";" | myChars == ",")
  targetSplits = seq(targetWidth, nchar(myString), targetWidth); 
  for (k in rev(targetSplits)) { ## from end
    if (min(abs(bestsplit - k)) < 11) {
      myDex = bestsplit[which.min(abs(bestsplit - k))] 
      myString = paste(substr(myString,1,myDex-1),"\n",substr(myString,myDex+1,nchar(myString)), sep="")
    } else if (min(abs(nextbest - k)) < 11) {
      myDex = nextbest[which.min(abs(nextbest - k))]
      myString = paste(substr(myString,1,myDex),"\n",substr(myString,myDex+1,nchar(myString)), sep="")
    } else {
      myDex = k;
      myString= paste(substr(myString,1,myDex),"-\n",substr(myString,myDex+1,nchar(myString)), sep="")
    }          
  }
  return(myString);
}  


## construct coloring of AAs automatically based on order of log-likehood ratios 
## for PrD vs background frequencies in yeast.
autocolors = TRUE;

if (autocolors) {
  aanames = c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y');
  if (1==0) { 
    scer.fg = read.table("prd_freq_scer_29.txt")[2:21,1]
    scer.bg = read.table("scer_bg.txt")[2:21,1]
    scer.bg = scer.bg/sum(scer.bg)
    names(scer.bg)=aanames;
    scer.ll = log2(scer.fg/scer.bg)
    aaordered = aanames[order(-scer.ll)]
    jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    jet100 = rev(jet.colors(100)[1:100]); ## snip off ends?
    ub = max(abs(scer.ll))
    jet.index = round(49*(scer.ll/ub) + 50)
    col.ll = jet100[jet.index];
    names(col.ll)=names(scer.ll);
    col.ll = col.ll[order(-scer.ll)]
    aaColors = col.ll;
  }
  ## hardcoded order based on yeast ordering above
  #aaordered = c("N","Q","Y","G","M","S","P","A","H","T","F","R","V","I","D","L","K","C","W","E")
  
  ## from help for colorRamp {grDevices}
  #jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  #aaColors = jet.colors(22)[2:21]; ## snip off ends, which have very high saturation 
  #names(aaColors) = aaordered
  
  aaColors=c("#002FFF","#0043FF","#08FFF5","#1DFFE1","#32FFCC","#50FFAD","#5BFFA3","#A3FF5B","#ADFF50","#EBFF13", 
             "#FFF300","#FFF300","#FF8100","#FF6200","#FF4E00","#FF4300","#FF1000","#DC0000","#B20000","#7F0000")
  names(aaColors) = c("N","Q","Y","G","M","S","P","A","H","T","F","R","V","I","D","L","K","C","W","E")
  
  
} else { 
  ## user can also specify colors and ordering of residues in legend manually, e.g. to highlight charged residues
  aaColors = unlist(list(
    "Q"="blue3",
    "N"="blue2",
    "Y"="dodgerblue1",
    "G"="cadetblue1",
    "S"="cyan2",
    "T"="cyan1",
    "A"="aquamarine",
    "M"="lawngreen",
    "P"="grey70",
    "V"="greenyellow",
    "I"="yellow1",
    "L"="yellow1",
    "F"="yellow2",
    "W"="yellow2",
    "C"="orange",
    "H"="pink",
    "R"="magenta1",
    "K"="magenta2",
    "D"="red",
    "E"="red3"
    ));
}

### This plots the results (gene parse, hydropathy, FoldIndex) that are produced by plaac.java 
### when using the genelist option. It reads these results from the file "infile" and plots 
### them one gene per page in a pdf specified by "outfile". Setting showSeq=T causes the AA 
### sequence to be shown beneath the plots.

plot_seq = function(dat, 
                      showHMMProbs=TRUE,
                      showParses=FALSE, 
                      showAAColors=FALSE,
                      showSeq=FALSE,  
                      tracks=c("FI","PLAAC","PAPA","THRESH"),
                      seqUnderline=c("FIx2","PLAACx2","PAPAx2"),
                      seqColor="PAPA") {
  
  trackNames = c("HYDRO","CHARGE","FI","FIx2","PLAAC","PLAACx2","PAPA","PAPAx2","THRESH");
  trackColors = c("blue","pink","grey50","grey50","red","red","green3","green3","green3");
  names(trackColors) = trackNames;
  trackCutoff = c(0,0,0,0,0,0,0.05,0.05,0);
  trackScale = c(1,1,1,1,-1,-1,-1,-1,1); ## highlight when scale*score < -cutoff
  names(trackCutoff) = trackNames;
  names(trackScale) = trackNames;
  
  ## datAll = read.table(infile, header=TRUE, sep="\t", stringsAsFactors=F);
  ## test for valid input here?
  
  hmmDex = grep("^HMM", colnames(dat));
  hmmStates = colnames(dat)[hmmDex];
  hmmStates = sub("HMM.", "", hmmStates, fixed=T) 
  numStates = length(hmmDex);
  ## print(paste("numStates: ", numStates))
  
  ## recycles if more than 7 HMM states
  hmmColors = rep(c("black","red","green","blue","cyan", "magenta","yellow"), length.out=numStates);
  
  aaColorVec = aaColors[match(dat$AA,names(aaColors))];
  aaColorVec[which(is.na(aaColorVec))] = "#666666"; # grey
  
  ## resize based on showHMMProbs and showSeq . Input files contains
  ## seqs of different lengths, so single dixed image size won't be ideal 
  ## for all of them when showSeq is T
  
  ## first column gives plot order, which can be different than order in file!
  ## for (dex in sort(unique(datAll$ORDER))) {
  
  # par(mfrow=c(1,1), family="Courier", mar=c(1,4,4,1)+0.1);
  
  n = nrow(dat);
  nn = n; ## width for some horizontal lines 
  ## use pretty() here?
  if (n > 50 & (n %% 50) > 30) {nn = n + 50 - (n %% 50);} # round up to next multiple of 50
  
  ## wasteful to repeat name in every row --- could put at start of output, or in just first row.
  cat(paste("making plots for", dat[1,"SEQid"],"\n"))
  
  titl = hyphenate(paste(dat[1,"SEQid"]), targetWidth=90);
  titSep = which(strsplit(titl,"")[[1]]=="|");
  
  if (length(titSep) > 0) {
    titl = substr(titl, max(titSep)+1, nchar(titl))
  }
  
  ##########################################################################
  ################################### plot HMM output  #####################
  ##########################################################################
  
  if (showHMMProbs) {
    plot(1, xlim=c(-0.01*n, 1.2*n), ylim=c(-0.07/1.2 - 4.5, 3), ##??
         type="n", xlab="", ylab="", axes=F);
    title(titl, cex.main=ifelse(nchar(titl) < 200, 1, 0.7), line=0);
    
    lines(c(-0.02*n,1.02*n),c(2,2),col="grey80");
    lines(c(-0.02*n,1.02*n),c(3,3),col="grey80");
    
    for (i in seq(0, n, by=50)) {lines(c(i,i),c(-1,3),col="grey80")}
    
    axis(2, at=c(-1,0,1));
    axis(2, at=c(2,3),labels=c(0,1));  
    
    vv = 2;
    
    for (i in seq_len(numStates)) {
      lines(seq_len(n), dat[,hmmDex[i]]+vv,type="l",col=hmmColors[i], lwd=2);
    }
    
    legend(1.03*n, 3, hmmStates, col=hmmColors, 
           text.col="black", lwd=2, merge=TRUE,  cex=1.0);
    
    ## only an option if showHMMProbs is TRUE; get rid of vacant space if showParse is FALSE?
   
    if (showAAColors) {
      for (i in seq_along(aaColorVec)) {
        lines(c(i-1,i), rep(-0.25+vv,2), col=aaColorVec[i], lwd=12, lend=1);
      }
    }
    
    if (showParses) {
      vitParse = dat$VIT;
      mapParse = dat$MAP;
        
      diffParse = diff(mapParse)
      segStart = c(0, which(diffParse!=0)) ## starts of runs in parse
      segEnd = c(which(diffParse!=0), length(diffParse)+1) ## ends of runs in parse
      segState = mapParse[segEnd]+1;
      
      for (i in seq_along(segStart)) {
        lines(c(segStart[i],segEnd[i]), rep(-0.18+vv,2), 
              col=hmmColors[segState[i]], bg=hmmColors[segState[i]], lwd=4, lend=1);
      }
      
      diffParse = diff(vitParse)
      segStart = c(0, which(diffParse!=0))
      segEnd = c(which(diffParse!=0), length(diffParse)+1)
      segState = vitParse[segEnd]+1;
      
      for (i in seq_along(segStart)) {
        lines(c(segStart[i],segEnd[i]), rep(-0.32+vv,2), 
              col=hmmColors[segState[i]], bg=hmmColors[segState[i]], lwd=4, lend=1);
      }
      
      text(0,-0.18+vv,"MAP", pos=2);
      text(0,-0.32+vv,"Vit", pos=2);
      
    }
    
  } else if (showSeq) { ## skip HMM output
    plot(1, xlim=c(-0.01*n, 1.2*n), ylim=c(-0.07/1.2 - 2.8, 1.0), 
         type="n", xlab="", ylab="",axes=F);
    title(titl, cex.main=ifelse(nchar(titl) < 200, 1, 0.7), line=3);
    axis(2,at=c(-1,-0.5,0,0.5,1), las=2)  
  } else {  
    plot(1, xlim=c(-0.01*n, 1.2*n), ylim=c(-1.0, 1.0), 
         type="n", xlab="", ylab="",axes=F);
    title(titl, cex.main=ifelse(nchar(titl) < 200, 1, 0.7), line=3);
    axis(2,at=c(-1,-0.5,0,0.5,1), las=2)   
  }
  
  ##########################################################################
  ######################## plot sliding-average tracks #####################
  ##########################################################################
 
  vo = 0;
  
  if (is.element("FI",tracks)) {
    polydex = which(!is.na(dat[,"FI"]));
    polygon(x=c(polydex,rev(polydex)), 
            y=c(dat[polydex,"FI"], rep(0, length(polydex))) +vo, 
            col=trackColors["FI"], border=trackColors["FI"])
  }
  if (is.element("FIx2",tracks)) {
    polydex = which(!is.na(dat[,"FIx2"]));
    polygon(x=c(polydex,rev(polydex)), 
            y=c(dat[polydex,"FIx2"], rep(0, length(polydex))) +vo, 
            col=trackColors["FIx2"], border=trackColors["FIx2"])
  }
  lines(c(-100,nn),c(vo,vo),col="black",lty=1, lwd=1); 
  
  showdex = seq_len(n); ## implicitly restricted to !is.na positions
  ## if protein has length 1, or if only one position is not na, lines don't get plotted.
  ## could plot them as points instead (type=o, pch=16)
  
  if (is.element("CHARGE",tracks)) 
    lines(showdex,abs(dat[showdex,"CHARGE"])+vo, col=trackColors["CHARGE"],lwd=2);
  if (is.element("HYDRO",tracks)) 
    lines(showdex,9/4*(dat[showdex,"HYDRO"]-0.5)+vo, col=trackColors["HYDRO"],lwd=2);
  if (is.element("PLAAC",tracks)) 
    ## originally scaled by 1.3 in sliding average plot, now by 1.386 by using base4 logs
    lines(showdex,-dat[showdex,"PLAAC"]/log(4) + vo,col=trackColors["PLAAC"],lwd=2);        
  if (is.element("PLAACx2",tracks)) 
    lines(showdex,-dat[showdex,"PLAACx2"]/log(4) + vo,col=trackColors["PLAACx2"],lwd=2); 
  
  papaScale = -4;
  papaScore = papaScale*dat[,"PAPA"]
  papaScore2 = papaScale*dat[,"PAPAx2"]
  papaThresh = 0.05*papaScale;
  
  if (is.element("PAPA",tracks)) {
    lines(showdex, pmin(papaScore[showdex], 1)+vo,col=trackColors["PAPA"], lwd=2); 
    lines(c(-0.02*n,1.02*n),c(1,1),col="white",lwd=3); ## cover up capped values 
    ## (for white rather than transparent background!)
  }
  
  if (is.element("PAPAx2",tracks)) {
    lines(showdex, pmin(papaScore2[showdex], 1)+vo,col=trackColors["PAPAx2"], lwd=2); 
    lines(c(-0.02*n,1.02*n),c(1,1),col="white",lwd=3); ## cover up capped values 
    ## (for white rather than transparent background!)
  }
  
  if (is.element("THRESH",tracks)) {
    lines(c(-100,nn),rep(papaThresh,2),col=trackColors["THRESH"],lwd=1.5, lty=2); 
  }
  
  ## should be smarter about axis tick postioning
  if (nn > 100) { 
    axis(3, at=seq(50,nn,by=100), labels=NA, pos=1.0, tcl=-0.25);
    axis(3, at=seq(0,nn,by=100), labels=seq(0,nn,by=100), pos=1.0, tcl=-0.5, cex.axis=1-min(0.5, floor(nn/1000)/10));
  } else if (nn > 10) {
    axis(3, at=seq(5,nn,by=10), labels=NA, pos=1.0, tcl=-0.25);
    axis(3, at=seq(0,nn,by=10), labels=seq(0,nn,by=10), pos=1.0, tcl=-0.5, cex.axis=1-min(0.5, floor(nn/1000)/10));
  } else {
    axis(3, at=seq(0,nn,by=1), labels=seq(0,nn,by=1), pos=1.0, tcl=-0.5, cex.axis=1-min(0.5, floor(nn/1000)/10));
  }
  
  wtracks = which(is.element(c("CHARGE","HYDRO", "FI","PLAAC","PAPA"),sub("x2$","",tracks)));
  
  if (length(wtracks) > 0) {
    legend(1.03*n, 0.7, c("Net charge","Hydropathy", "FoldIndex","-PLAAC","-4*PAPA")[wtracks],
           col = trackColors[c("CHARGE","HYDRO", "FI","PLAAC","PAPA")[wtracks]], 
           text.col="black", lwd=2, merge=TRUE, ,cex=1.0);
  }
  
  ##########################################################################
  ############################# plot protein sequence  #####################
  ##########################################################################
  
  if (showSeq) {
    
    nc = 100; ## number of chars per line
    aaseq = as.character(dat$AA);
    colseq = rep("black", length(aaseq)); ## per-position colors
    
    if (seqColor == "COMBO") {
      colseq[which(dat$FIx2 < 0 & dat$PAPAx2 > 0.05)] = "cyan"
    } else if (seqColor=="VIT") {
      colseq = hmmColors[1+dat$VIT];
    } else if (seqColor=="MAP"){
      colseq = hmmColors[1+dat$MAP];
    } else if (seqColor=="PLAAC") {
      colseq[which(dat$PLAAC > 0)] = "red";
    } else {
      ## add others if desired...
    }
    
    voff = 0.18; ## decrease when n > 1300?
    hoff = 1;       
    vstart = -1.2;
    
    ## check when length(aaseq %% nc)==0 !
    for (i in 1:floor(length(aaseq)/nc + 1)) {
      text(0,vstart-(i-1)*voff, paste((i-1)*nc + 1), pos=2, col="grey40", cex=0.92); 
      seqdex = ((i-1)*nc + 1):(i*nc);
      text((seq_along(seqdex)-1)*1.15*n/nc, vstart-(i-1)*voff, 
           aaseq[seqdex], pos=4, col=colseq[seqdex], cex=0.92);
    }   
    
    stagger = -0.015;
    for (trackName in seqUnderline) {
      ## don't underline first and last w positions? Or assume they are masked with NaN
      if (is.element(trackName, trackNames)) { 
        stagger = stagger + 0.015; ## could adapt based on number of underline Tracks
        ulseq = rep(" ", length(aaseq));
        ulseq[which(trackScale[trackName]*dat[,trackName] < -trackCutoff[trackName])] = "_";
        for (i in 1:floor(length(aaseq)/nc + 1)) {
          seqdex = ((i-1)*nc + 1):(i*nc);
          text((seq_along(seqdex)-1)*1.15*n/nc, vstart-(i-1)*voff-stagger, 
               ulseq[seqdex], pos=4, col=trackColors[trackName], cex=1.1); 
        }   
      }
    }
  } 

}



########################
########################

## infile = "../../oktestplot.txt"

plot_seqs = function(infile="plot_data.txt", 
                      outfile="plaac_plots.pdf",
                      showHMMProbs=TRUE,
                      showParses=FALSE, 
                      showAAColors=FALSE,
                      showSeq=FALSE,  
                      tracks=c("FI","PLAAC","PAPA","THRESH"),
                      seqUnderline=c("FIx2","PLAACx2","PAPAx2"),
                      seqColor="PAPA"){
  
  # datAll = read.table(infile, header=TRUE, sep="\t", stringsAsFactors=F);
  ## changed from default comment.char="#" in order to allow # in fasta ID
  datRaw = scan(infile, what="character", sep="\n", quiet=T, quote="")
  datRaw = datRaw[!grepl("^#", datRaw)] 
  datAll = read.table(text=datRaw, header=TRUE, sep="\t", stringsAsFactors=F, quote="", comment.char="")
  rm(datRaw)
  
  ## test for valid input here?
  
  plotWidth = 12;
  plotHeight = 3 + ifelse(showHMMProbs, 3, 0) + ifelse(showSeq, 4.5, 0)
  
  plotPng = grepl("png$", outfile);
  
  if (plotPng) {
    pngStub = sub(".png$","", outfile); 
    pngCount=1;
  }
  
  # if (grepl("png$", outfile)) {
  #  pngOutfile = sub("\\.png$", "_%05d.png", outfile)
  #  png(file=pngOutfile, bg="transparent", width=plotWidth, height=plotHeight, units="in", 
  #      res=72, pointsize=12, family="Courier"); 
  #   ## On Mac get overplotting for later pngs; fixed with Type=X11, 
  #   ## but then Courier isn't found, and mono font doesn't look good.
  #   ## Could do manualy loop instead...
  #} else {
  #  pdf(file=outfile, width=plotWidth, height=plotHeight, pointsize=12, family="Courier");
  # }
  
  if (!plotPng) {pdf(file=outfile, width=plotWidth, height=plotHeight, pointsize=12, family="Courier");}
  
  ## first column gives plot order, which can be different than order in file!
  for (dex in sort(unique(datAll$ORDER))) {
    
    if (plotPng) {
      png(file=paste(pngStub, "_", sprintf("%05d",pngCount),".png", sep=""), bg="transparent", 
          width=plotWidth, height=plotHeight, units="in", res=72, pointsize=12, family="Courier"); 
      pngCount = pngCount+1;
    }
    
    dat = datAll[which(datAll$ORDER==dex),];
    # dat = datAll[datAll$ORDER==1,];
    
    ## use ellipses ... to avoid explicitly passing the args?
    par(mfrow=c(1,1), family="Courier", mar=c(1,4,4,1)+0.1);
    plot_seq(dat, showHMMProbs=showHMMProbs, showParses=showParses, showAAColors=showAAColors, showSeq=showSeq,  
              tracks=tracks, seqUnderline=seqUnderline, seqColor=seqColor)
    
    if (plotPng) {dev.off();}
  }
  if (!plotPng) {dev.off();}
  
}


#######

## allowing range of n to show would be more flexible than 1:max_n
color_code_seqs = function(infile="plot_data.txt", outfile="color_plots.pdf", max_n=NA, showLegend=T, showParses=T) {
  # datAll = read.table(infile, header=TRUE, sep="\t", stringsAsFactors=F);
  ## changed from default comment.char="#" in order to allow # in fasta ID
  datRaw = scan(infile, what="character", sep="\n", quiet=T, quote="")
  datRaw = datRaw[!grepl("^#", datRaw)] 
  datAll = read.table(text=datRaw, header=TRUE, sep="\t", stringsAsFactors=F, quote="", comment.char="")
  rm(datRaw)
  
  ## use positive max_n to enforce max_length and/or to include padding so different plots will have same scale.
  if (is.na(max_n)) max_n = max(datAll$AANUM);          ## use max_n = NA to put no limit on length, but use no extra padding
  if (max_n < 0) max_n = min(-max_n, max(datAll$AANUM)) ## use negative max_n to impose cap, but avoid padding if longest protein is < -max_n
  if (max_n < max(datAll$AANUM)) {
    print(paste("## Warning: longest protein has length", max(datAll$AANUM)))
    print(paste("## This (and potentially others) will be truncated at", max_n, "residues"))
    print("## You can increase max_n in function call, or set it to NA for no limit.")
  }
  
  hmmColors = rep(c("black","red","green","blue","cyan", "magenta","yellow"), 
                  length.out=length(grep("^HMM", colnames(datAll))));
 
  unique_order = sort(unique(datAll$ORDER));
  num_seq = length(unique_order);
  
  plotPng = grepl("png$", outfile);
  
  if (plotPng) { 
    png(outfile, height=100*(0.6 + 0.3*(num_seq + ifelse(showLegend, 1.5, 0))), width=100*8, units="px", family="Courier", pointsize=16)
  } else{ 
    pdf(outfile, height=0.5 + 0.3*(num_seq + ifelse(showLegend, 1.5, 0)), width=8, family="Courier", pointsize=12)
  }
  
  par(mfrow=c(1,1), mar=c(1,1,1,1)+0.1)
  plot(1, type="n", axes=F, xlab="", ylab="", xaxs="i", yaxs="i", xlim=c(-0.15*max_n, 1.02*(max_n)), 
       ylim=c(0.6, num_seq + 0.4 + ifelse(showLegend,1.6,0)));
  for (k in seq_along(unique_order)) {
    dex = unique_order[k];
    dat = datAll[datAll$ORDER==dex,]; 
    aaColorVec = aaColors[match(dat$AA,names(aaColors))];
    aaColorVec[is.na(aaColorVec)]=="#666666";
    for (i in seq_len(min(length(aaColorVec), max_n))) {
      lines(c(i-1,i), rep(num_seq - k + 1, 2), col=aaColorVec[i], lwd=14, lend=1);
    }
    if (length(aaColorVec) > max_n) {
      print(paste("## truncating sequence", dat$SEQid[1]))
      lines(x=c(max_n, 1.01*max_n, 1.005*max_n, 1.01*max_n, 1.005*max_n, 1.01*max_n, max_n), 
              y=num_seq-k+1 + c(0.2, 0.2, 0.08, 0, -0.08, -0.2, -0.2), col="grey40", lwd=2)
      # arrows(x0=max_n, x1=max_n+3, y0=num_seq-k+1, y1=num_seq-k+1, col="grey", lwd=2, code=2, length=0.1);
    }
    ## truncate name to 10 characters so it can fit in left margin --- this should accommodate monospace fonts in pdf;
    ## for variable width fonts in png output, it's possible that ten-letter names with lots of wide letters may go 
    ## out of bounds, but this may not be an issue in practice.  
    text(0, rep(num_seq - k + 1, 2), substr(dat$SEQid[1],1,10), pos=2, cex=0.8, adj=0)
    if (showParses) {
      myParse = dat$VIT;
      ## myParse = dat$MAP;       
      diffParse = diff(myParse)
      segStart = c(0, which(diffParse!=0))
      segEnd = c(which(diffParse!=0),length(diffParse)+1)
      capDex = max(which(segStart < max_n)) 
      segEnd = pmin(segEnd, max_n);
      segState = myParse[segEnd]+1;
      for (i in seq_len(capDex)) {
        lines(c(segStart[i],segEnd[i]), rep(num_seq - k + 1 - 0.20, 2), pch=22,
              col=hmmColors[segState[i]], bg=hmmColors[segState[i]], lwd=3, lend=1);
        lines(c(segStart[i],segEnd[i]), rep(num_seq - k + 1 + 0.20, 2), pch=22,
              col=hmmColors[segState[i]], bg=hmmColors[segState[i]], lwd=3, lend=1);
      }
    }
  }
  if (showLegend) {
    midx = seq(0.25*max_n, 0.75*max_n, length.out=length(aaColors))
    startx = midx - 0.5*max_n/length(aaColors)/3;
    endx = midx + 0.5*max_n/length(aaColors)/3;
    for (k in seq_along(aaColors)) {
      lines(c(startx[k],endx[k]), rep(num_seq + 1 + 0.1, 2), col=aaColors[k], lwd=12, lend=1);
    }  
    text(midx, rep(num_seq + 1 + 0.15, length(aaColors)), names(aaColors), pos=3, cex=1)
  }
  dev.off()
  
}

#color_seqs("~/Dropbox/plaac/scer_prion_plot.txt", "color_test.pdf", max_n=NA, showlegend=T)
#color_seqs("~/Dropbox/plaac/scer_prion_plot.txt", "color_test.png", showlegend=T)
