## Original program copyright 2009 Whitehead Institute for Biomedial Research
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


### This plots the results (gene parse, hydropathy, FoldIndex) that are produced by spewey.java 
### when using the genelist option. It reads these results from the file "infile" and plots 
### them one gene per page in a pdf specified by "outfile". Setting showSeq=T causes the AA 
### sequence to be shown beneath the plots.

############### Oct 2013
plot_genes = function(infile="plot_data.txt", 
                      outfile="plaac_plots.pdf",
                      showHMMProbs=TRUE,
                      showParses=FALSE, 
                      showSeq=FALSE,  
                      tracks=c("FI","PLAAC","PAPA","THRESH"),
                      seqUnderline=c("FIx2","PLAACx2","PAPAx2"),
                      seqColor="PAPA",
                      maxseq=500, w=21){
  
  trackNames = c("HYDRO","CHARGE","FI","FIx2","PLAAC","PLAACx2","PAPA","PAPAx2","THRESH");
  trackColors = c("blue","pink","grey50","grey50","red","red","green3","green3","green3");
  names(trackColors) = trackNames;
  trackCutoff = c(0,0,0,0,0,0,0.05,0.05,0);
  trackScale = c(1,1,1,1,-1,-1,-1,-1,1); ## highlight when scale*score < -cutoff
  names(trackCutoff) = trackNames;
  names(trackScale) = trackNames;

  datAll = read.table(infile, header=TRUE, sep="\t", stringsAsFactors=F);
  ## test for valid input here?
  
  hmmDex = grep("^HMM", colnames(datAll));
  hmmStates = colnames(datAll)[hmmDex];
  hmmStates = sub("HMM.", "", hmmStates, fixed=T) 
  numStates = length(hmmDex);
  ## print(paste("numStates: ", numStates))
  
  ## recycles if more than 7 HMM states
  hmmColors = rep(c("black","red","green","blue","cyan", "magenta","yellow"), length.out=numStates);
  
  ## resize based on showHMMProbs and showSeq . Input files contains
  ## seqs of different lengths, so single dixed image size won't be ideal 
  ## for all of them when showSeq is T
  
  plotWidth = 12;
  plotHeight = ifelse(showHMMProbs,10,6);
  
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
  

  ## first column gives plot order, which can be different than order in file
  for (dex in sort(unique(datAll$ORDER))) {
    
    if (plotPng) {
      png(file=paste(pngStub, "_", sprintf("%05d",pngCount),".png", sep=""), bg="transparent", 
          width=plotWidth, height=plotHeight, units="in", res=72, pointsize=12, family="Courier"); 
      pngCount = pngCount+1;
    }
    
    par(mfrow=c(1,1), family="Courier", mar=c(1,4,4,1)+0.1);
    
    dat = datAll[datAll$ORDER==dex,];
    n = nrow(dat);
    nn = n; ## width for some horizontal lines 
    if ((n %% 50) > 30) {nn = n + 50 - (n %% 50);} # next multiple of 50
    
    ## wasteful to repeat name in every row --- could put at start of output, or in just first row.
    cat(paste("making plots for", dat[1,"GENE"],"\n"))
    
    titl = hyphenate(paste(dat[1,"GENE"]), targetWidth=90);
    titSep = which(strsplit(titl,"")[[1]]=="|");
    
    if (length(titSep) > 0) {
      titl = substr(titl, max(titSep)+1, nchar(titl))
    }
    
    ##########################################################################
    ################################### plot HMM output  #####################
    ##########################################################################
    
    if (showHMMProbs) {
      plot(1, xlim=c(-0.01*n, 1.2*n), ylim=c(-0.07/1.2 -4, 3), 
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
      if (showParses) {
        vitParse = dat$VIT;
        mapParse = dat$MAP;
        
        diffParse = diff(mapParse)
        segStart = c(1, which(diffParse!=0))
        segEnd = c(which(diffParse!=0),length(diffParse)+1)
        segState = mapParse[segEnd]+1;
        
        for (i in seq_along(segStart)) {
          ## stop a little before?
          lines(c(segStart[i],segEnd[i]-0.5), rep(-0.18+vv,2), pch=22,
                col=hmmColors[segState[i]], bg=hmmColors[segState[i]], lwd=5, lend=2);
        }
        
        diffParse = diff(vitParse)
        segStart = c(1, which(diffParse!=0))
        segEnd = c(which(diffParse!=0), length(diffParse)+1)
        segState = vitParse[segEnd]+1;
        
        for (i in seq_along(segStart)) {
          lines(c(segStart[i],segEnd[i]-0.5), rep(-0.32+vv,2), pch=22,
                col=hmmColors[segState[i]], bg=hmmColors[segState[i]], lwd=5, lend=2);
        }
        
        text(0,-0.18+vv,"MAP", pos=2);
        text(0,-0.32+vv,"Vit", pos=2);
      }
      
    } else { ## skip HMM output
      plot(1, xlim=c(-0.01*n, 1.2*n), ylim=c(-0.07/1.2 -2.3, 1.0), 
           type="n", xlab="", ylab="",axes=F);
      title(titl, cex.main=ifelse(nchar(titl) < 200, 1, 0.7), line=3);
      axis(2,at=c(-1,-0.5,0,0.5,1), las=2) 
      
    }      
    
    ##########################################################################
    ######################## plot sliding-average tracks #####################
    ##########################################################################
    
    vo = 0;
    
    ## make w1 and w2, for singly and double smoothed? Or use NA is output file
    
    if (is.element("FI",tracks)) polygon(x=c(w:(n-w+1),(n-w+1):w), 
                                         y=c(dat[w:(n-w+1),"FI"], 0*(w:(n-w+1)))+vo, 
                                         col=trackColors["FI"], border=trackColors["FI"])
    if (is.element("FIx2",tracks)) polygon(x=c(w:(n-w+1),(n-w+1):w), 
                                         y=c(dat[w:(n-w+1),"FIx2"], 0*(w:(n-w+1)))+vo, 
                                         col=trackColors["FIx2"], border=trackColors["FIx2"])
    lines(c(-100,nn),c(vo,vo),col="black",lty=1, lwd=1); 
    
    if (is.element("CHARGE",tracks)) 
      lines(w:(n-w),abs(dat[w:(n-w),"CHARGE"])+vo, col=trackColors["CHARGE"],lwd=2);
    if (is.element("HYDRO",tracks)) 
      lines(w:(n-w),9/4*(dat[w:(n-w),"HYDRO"]-0.5)+vo, col=trackColors["HYDRO"],lwd=2);
    if (is.element("PLAAC",tracks)) 
      lines(w:(n-w+1),-dat[w:(n-w+1),"PLAAC"]/log(4) + vo,col=trackColors["PLAAC"],lwd=2); ## changed from 1.3 --- now base4 logs
    if (is.element("PLAACx2",tracks)) 
      lines(w:(n-w+1),-dat[w:(n-w+1),"PLAACx2"]/log(4) + vo,col=trackColors["PLAACx2"],lwd=2); ## changed from 1.3 --- now base4 logs
    
    papaScale = -4;
    papaScore = papaScale*dat[,"PAPA"]
    papaScore2 = papaScale*dat[,"PAPAx2"]
    papaThresh = 0.05*papaScale;
    
    if (is.element("PAPA",tracks)) {
      lines(w:(n-w+1), pmin(papaScore[w:(n-w+1)], 1)+vo,col=trackColors["PAPA"], lwd=2); 
      lines(c(-0.02*n,1.02*n),c(1,1),col="white",lwd=3); ## cover up capped values 
      ## (for white rather than transparent background!)
    }
    
    if (is.element("PAPAx2",tracks)) {
      lines(w:(n-w+1), pmin(papaScore2[w:(n-w+1)], 1)+vo,col=trackColors["PAPAx2"], lwd=2); 
      lines(c(-0.02*n,1.02*n),c(1,1),col="white",lwd=3); ## cover up capped values 
      ## (for white rather than transparent background!)
    }
    
    if (is.element("THRESH",tracks)) {
      lines(c(-100,nn),rep(papaThresh,2),col=trackColors["THRESH"],lwd=1.5, lty=2); 
    }
    
    axis(3, at=seq(0,nn,by=50), labels=ifelse(seq(0,nn,by=50)%%100 == 0, seq(0,nn,by=50), NA), pos=1.0);
    
    wtracks = which(is.element(c("CHARGE","HYDRO", "FI","PLAAC","PAPA"),sub("x2$","",tracks)));
    
    if (length(wtracks) > 0) {
      legend(1.03*n, 0.7, c("Net charge","Hydropathy", "FoldIndex","-PLAAC","-4*PAPA")[wtracks],
             col = trackColors[c("CHARGE","HYDRO", "FI","PLAAC","PAPA")[wtracks]], 
             text.col="black", lwd=2, merge=TRUE, ,cex=1.0);
    }
    
    ##########################################################################
    ############################# plot protein sequecne  #####################
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
        ## don't underline first and last 2 positions?
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
    if (plotPng) {dev.off();}
  }
  if (!plotPng) {dev.off();}
  
}


