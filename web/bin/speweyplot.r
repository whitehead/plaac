### This plots the results (gene parse, hydropathy, FoldIndex) that are produced by spewey.java 
### when using the genelist option. It reads these results from the file "infile" and plots 
### them one gene per page in a pdf specified by "outfile". Setting seq=T causes the AA 
### sequence to be shown beneath the plots.

plotgenes <- function(infile="combo.txt", seq=FALSE, use_pdf=TRUE,
                  outfile=paste("Rplot-",strsplit(infile,"\\.")[[1]][1], ".pdf", sep="")){

    if (use_pdf) {
      pdf(file=outfile, width=10, height=10, pointsize=10);
    } else {
      png_files = paste(strsplit(outfile,"\\.")[[1]][1], "%05d.png", sep="")
      png(file=png_files, bg="transparent", width=10, height=10, units="in", res=75, pointsize=10);
    }

    hmmur <- read.table(infile,header=TRUE, sep="\t");
    #hmmx <- hmmur[c(1:2,4:6,11:length(hmmur))];
    hmmx <- hmmur[c(1:2,4:6,12,13)];
    hydrox <- hmmur[c(1,2,4,7:10,11)];

    namevec <- names(hmmx);
    namevec <- namevec[6:length(namevec)];

    for (d in 1:max(hmmx[,1])) {

        dex <-d;

        colorvec = c("black","red","blue","green","yellow", "pink","orange");

        keep <- (hmmx[,1]==dex);
      
        if (sum(keep)<1) {next;}

        hmm <- hmmx[keep,];
        n=dim(hmm)[1];

        titl <- paste(hmm[1,2]);
   
        numstates <- length(hmm)-5;
        par(family="mono");
        plot((-0.01*n):(1.2*n), 7*((-0.01*n):(1.2*n))/(1.2*n)-4,type="n",xlab="",ylab="",axes=F);
       
        abline(h=-1,col="lightgrey");
        abline(h=0,col="lightgrey");
        abline(h=1,col="lightgrey");
        abline(h=2,col="lightgrey");
        abline(h=3,col="lightgrey");
        for (i in 0:(n/50)){ lines(c(50*i,50*i),c(-1,3),col="lightgrey");};

        axis(2,-c(-1,0,0,1));
        axis(3,seq(0,n,by=50),pos=1.0);
        axis(2,c(2,3),c(0,1));  
           
        vv <- 2;

        for (i in 1:numstates) {
          # points(1:n,10*(hmm[,4]==(i-1))-10.30+vv,pch=20,col=colorvec[i],lwd=1);      
          # points(1:n,10*(hmm[,5]==(i-1))-10.16+vv,pch=20,col=colorvec[i],lwd=1); 
          # points(1:n,10*(hmm[,4]==(i-1))-10.30+vv,pch=22,col=colorvec[i],bg=colorvec[i],lwd=1);      
          # points(1:n,10*(hmm[,5]==(i-1))-10.16+vv,pch=22,col=colorvec[i],bg=colorvec[i],lwd=1); 
           
          # lines(1:n,hmm[,i+5]+vv,type="l",col=colorvec[i],lwd=2);
           lines(1:n,hmm[,i+5]+vv,type="l",col=colorvec[i],lwd=2);
        }
        title(titl);
        #text(0,-0.18+vv,"MAP",pos=2);
        #text(0,-0.32+vv,"Vit",pos=2);
        legend(1.03*n, 3, namevec[1:numstates], col=colorvec[1:numstates], 
             text.col= "black", lwd = 3, pch=-1, merge = TRUE, bg='white',cex=1.0);
       
         keep <- (hydrox[,1]==dex);
         hydro <-hydrox[keep,];
        
         vo <- 0;         

         n=dim(hydro)[1]; w = 1;
        
         lines(w:(n-w),pmax(hydro[w:(n-w),6],0)+vo,type="h",col="grey",lwd=3);
         lines(w:(n-w),pmin(hydro[w:(n-w),6],0)+vo,type="h",col="grey",lwd=3);
         lines(w:(n-w),hydro[w:(n-w),6]+vo,col="grey");
         lines(w:(n-w),abs(hydro[w:(n-w),4])+vo,col="pink",lwd=2);
         lines(w:(n-w),9/4*(hydro[w:(n-w),5]-0.5)+vo,col="blue",lwd=2);
         lines(w:(n-w),-hydro[w:(n-w),7]/1.3+vo,col="red",lwd=3);

        # akl: ross1 line
         lines(w:(n-w),-hydro[w:(n-w),8]/1.3+vo,col="black",lwd=2);
        
         legend(1.03*n,0.5,c("net charge","hydropathy","foldindex","-prd score","-ross1"), col = c("pink","blue","grey","red", "black"), 
             text.col= "black", lwd = c(2,2,2,2), pch = c(-1, -1, -1, -1), merge = TRUE, bg='white',cex=1.0);


        if (seq) {
           
              nc <- 100;
              dex <- d;

              colorvec = c("black","red","blue","green","yellow", "pink","orange");
   
              aanames <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*");
              aanums <- hmm[,3];
     
              par(family="mono");
              voff <- 0.16;
              hoff <- 1;

              if (d == 2*round(d/2)) {evenseq <- paste(aanames[aanums]);}
              else {oddseq <- paste(aanames[aanums]);}

              for (targ in 4:4) {

                vstart <- -1.6;
                 
                for (j in 1:numstates) { 
                  ## replace chars by spaces if don't want them colored
                  mask <- !(hmm[,targ]==(j-1));     
                  # str0= aanums;
                  str0 <-paste(aanames[aanums]);
                  for (i in (1:length(hmm[,3]))[mask])  {str0[i]=" ";}
                  str <- paste(str0,collapse='');   
                  for (i in 1:(nchar(str)/nc + 1)) {
                      text(0,vstart-(i-1)*voff, paste((i-1)*nc + 1), pos=2,col="black");
                      text(0,vstart-(i-1)*voff, substr(str,(i-1)*nc + 1, i*nc), pos=4,col=colorvec[j]);
                  }   
              }
              cv <- c("red","red","red","pink","blue","grey","red");

              for (j in 6:7) {
                 mask <- (hydro[,j]< 0.0); 
                 if (j==6) mask = !mask;  
                 str0 <-paste(aanames[aanums]);
                 for (i in (1:length(hydro[,3]))[mask])  {str0[i]=" ";}
                 for (i in (1:length(hydro[,3]))[!mask])  {str0[i]="_";}
                str <- paste(str0,collapse='');   
                for (i in 1:(nchar(str)/nc + 1)) {
                     text(0,vstart-(i-1)*voff, substr(str,(i-1)*nc + 1, i*nc), pos=4, col=cv[j]);
                }
            }
         }      
      }

    }
    dev.off();

}



