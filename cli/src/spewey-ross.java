/** /////////////////////////////////////////////////////////////////////////////
  Updated April 8, 2009, Oliver King (king@bbri.org)
  Compile with: javac spewey.java
  Run with:
  java -Xmx512M spewey fasta-file core-length alpha > output-file
  e.g.
  java -Xmx512M spewey scer.aa 60  0.5 > scer-60-out.txt

  alpha should be between 0 and 1; it controls the bg frequenicies:
  If alpha=1, it uses the S.cer bgfreqs,
  If alpha=0, it uses the bgfreqs from all the sequences in fasta-file.
  In general the bgfreqs are alpha*(scer) + (1-alpha)*(sequences in fasta-file)

  Updated July 8, 2009:
  Can also run with:
  java -Xmx512M spewey fasta-file core-length alpha genelist > output-file
  e.g.
  java -Xmx512M spewey scer.aa 60  0.5 somegenes.txt > combo.txt

  Here genelist is a file that has one geneID per line, with each geneID
  exactly matching a geneID in the fasta-file scer.aa. Each line can optionally
  also contain an alternate name to be used for labeling the plots, after a tab
  This program will output all the information needed to plot graphs showing
  the parses of each gene in genelist, and also hydropathy and charge and FoldIndex.

  For exmple, if one of the entries in the fasta-file scer.aa is
  >YDR172W
  MSDSNQGNNQQNYQQYSQNG...
  then one line in somegenes.txt could be
  YDR172W	SUP35
  with those two names seperated by a tab. Then the common name would be
  shown on the graph rather than the ORF ID.

  The routine to plot these is in the file "plotspewey.r"
  Load this into an R session using
  source("speweyplot.r")
  and then generate a pdf with the plots using
  plotgenes(infile, seq, outfile)
  where infile is the output of spewey above. If seq=T the AA sequence
  will be shown beneath each plot.

  You can first run spewey without the genelist option so that it scores all
  genes. Then by sorting the output of this you can select genes of interest
  and will have geneIDs in the right format to put into somegenes.txt
  so you can then rerun spewey to make graphs.


  java spewey allprd.aa 60 1 > scer-prd-out2.txt
  java spewey ../seqs/scer.aa 60 0 > scer-out-ross.txt


 **/ ////////////////////////////////////////////////////////////////////////


import java.io.*;
import java.util.*;
import java.net.*;
import java.util.regex.*;

class spewey {

    static int [] geneticcode; // length 125

    // 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21
    static char [] aanames =   {'X','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'};
    static char [] aanamesb =  {'-','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'};
    static char [] aanamesss = {'B','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'};
    static String [] aanames2 = {"???","Ala","Cys","Asp","Glu","Phe","Gly","His","Ile","Lys","Leu","Met","Asn",
        "Pro","Gln","Arg","Ser","Thr","Val","Trp","Tyr","***"};

    static char [] ntnames =  {'?','a','c','g','t'};
    static char [] ntnamesb = {'-','a','c','g','t'};
    static String [] codonnames = new String[125];
    static String [] scodonnames = new String[65];

    static Hashtable sys2com;
    static Hashtable com2sys;

    static double log2 = Math.log(2);

    static double [] loglut;
    static int loglutlength = 40*100;

    static double [] aadobson  = {
        0.0, // 'X'
        -3.31, // 'A'
        1.61, // 'C'
        -9.42, // 'D'
        -10.38, // 'E'
        2.80, // 'F'
        -3.96, // 'G'
        -4.31, // 'H'
        0.93, // 'I'
        -9.55, // 'K'
        -0.25, // 'L'
        -1.06, // 'M'
        -6.02, // 'N'
        -11.96, // 'P'
        -6.00, // 'Q'
        -11.93, // 'R'
        -5.08, // 'S'
        -2.12, // 'T'
        0.49, // 'V'
        2.92, // 'W'
        1.03, // 'Y'
        0.0  // '*'
    };

    static double [] aacharge = {
        0, // 'X'
        0, // 'A'
        0, // 'C'
        1, // 'D'
        1, // 'E'
        0, // 'F'
        0, // 'G'
        0, // 'H' -0.5?
        0, // 'I'
        -1, // 'K'
        0, // 'L'
        0, // 'M'
        0, // 'N'
        0, // 'P'
        0, // 'Q'
        -1, // 'R'
        0, // 'S'
        0, // 'T'
        0, // 'V'
        0, // 'W'
        0, // 'Y'
        0  // '*'
    };


    // K-D
    static double [] aahydro = {
        0.0, // 'X'
        1.8, // 'A'
        2.5, // 'C'
        -3.5, // 'D'
        -3.5, // 'E'
        2.8, // 'F'
        -0.4, // 'G'
        -3.2, // 'H'
        4.5, // 'I'
        -3.9, // 'K'
        3.8, // 'L'
        1.9, // 'M'
        -3.5, // 'N'
        -1.6, // 'P'
        -3.5, // 'Q'
        -4.5, // 'R'
        -0.8, // 'S'
        -0.7, // 'T'
        4.2, // 'V'
        -0.9, // 'W'
        -1.3, // 'Y'
        0.0  // '*'
    };

    // From Toombs, McCarty and Ross, MCB 2009
    static double [] fgross1 = {
        0.0, // 'X'
        0.042, // 'A'
        0.033, // 'C'
        0.014, // 'D'
        0.009, // 'E'
        0.075, // 'F'
        0.038, // 'G'
        0.059, // 'H'
        0.102, // 'I'
        0.009, // 'K'
        0.059, // 'L'
        0.038, // 'M'
        0.096, // 'N'
        0.038, // 'P'
        0.024, // 'Q'
        0.054, // 'R'
        0.125, // 'S'
        0.069, // 'T'
        0.102, // 'V'
        0.024, // 'W'
        0.054, // 'Y'
        0.0  // '*'
    };


    static double [] bgross1 = {
        0.0, // 'X'
        0.072, // 'A'
        0.022, // 'C'
        0.051, // 'D'
        0.017, // 'E'
        0.032, // 'F'
        0.040, // 'G'
        0.078, // 'H'
        0.045, // 'I'
        0.045, // 'K'
        0.061, // 'L'
        0.020, // 'M'
        0.089, // 'N'
        0.127, // 'P'
        0.022, // 'Q'
        0.081, // 'R'
        0.109, // 'S'
        0.078, // 'T'
        0.045, // 'V'
        0.012, // 'W'
        0.025, // 'Y'
        0.0  // '*'
    };



    static double [] fgross2 = {
        0.0, // 'X'
        0.057, // 'A'
        0.015, // 'C'
        0.045, // 'D'
        0.023, // 'E'
        0.064, // 'F'
        0.057, // 'G'
        0.080, // 'H'
        0.038, // 'I'
        0.004, // 'K'
        0.045, // 'L'
        0.030, // 'M'
        0.068, // 'N'
        0.072, // 'P'
        0.030, // 'Q'
        0.045, // 'R'
        0.110, // 'S'
        0.087, // 'T'
        0.038, // 'V'
        0.042, // 'W'
        0.049, // 'Y'
        0.0  // '*'
    };

    // D 4% too high, E too low?

    static double [] bgross2 = {
        0.0, // 'X'
        0.064, // 'A'
        0.012, // 'C'
        0.067, // 'D'
        0.024, // 'E'
        0.021, // 'F'
        0.046, // 'G'
        0.070, // 'H'
        0.030, // 'I'
        0.021, // 'K'
        0.052, // 'L'
        0.021, // 'M'
        0.040, // 'N'
        0.095, // 'P'
        0.037, // 'Q'
        0.076, // 'R'
        0.119, // 'S'
        0.095, // 'T'
        0.037, // 'V'
        0.009, // 'W'
        0.064, // 'Y'
        0.0  // '*'
    };
    // F 3% too high, R too low?



    static double [] lodross1;
    static double [] lodross2;
    static double [] llkross1;
    static double [] llkross2;

    static double [] lodross3;

    static double [] odross1 ={
        0.0,
        0.67267686,
        1.5146198,
        0.27887323,
        0.5460614,
        2.313433,
        0.96153843,
        0.75686276,
        2.2562358,
        0.20664589,
        0.9607843,
        1.9615384,
        1.0836071, // N
        0.30196398,// P
        1.0716166, // Q
        0.6664044,
        1.1432927,
        0.8917492,
        2.2562358,
        1.9478673,
        2.1785367,
        0.0
    };

    static double [] odross2 ={
        0.0,
        0.88066554,
        1.2461538,
        0.72039115,// ## 0.66?? [45 vs 67]
        0.77220076,// ## 0.93?? [23 vs 24]
        3.6936572, // ## 3.16?? [64 vs 21] ct 17/264 vs 6/328
        1.2570281,
        1.1460011,
        1.2519685,
        0.17436177,
        0.87114847,
        1.4330357,
        1.7729831, // N
        0.7429888, // P
        0.8229167, // Q
        0.5531136, //  ## 0.58?? [45 vs 76]
        0.9144572,
        0.9143354,
        1.0367454,
        4.710145,
        0.75716186,
        0.0,
    };





    static double pcp = 2.25; // prion charge penalty

    static double [] aaprion = {
        0.0, // 'X'
        1.8, // 'A'
        2.5, // 'C'
        pcp, // 'D'
        pcp, // 'E'
        2.8, // 'F'
        -0.4, // 'G'
        pcp/2, // 'H'
        4.5, // 'I'
        pcp, // 'K'
        3.8, // 'L'
        1.9, // 'M'
        -3.5, // 'N'
        -1.6, // 'P'
        -3.5, // 'Q'
        pcp, // 'R'
        -0.8, // 'S'
        -0.7, // 'T'
        4.2, // 'V'
        -0.9, // 'W'
        -1.3, // 'Y'
        0.0  // '*'
    };


    static double [] gnq = {0,0,0,0,0,0,1.0,0,0,0,0,0,1.0,0,1.0,0,0,0,0,0,0,0};
    static double [] scbg = {0,0.0550,0.0126,0.0586,0.0655,0.0441,0.0498,0.0217,0.0655,0.0735,0.0950,0.0207,
        0.0615,0.0438,0.0396,0.0444,0.0899,0.0592,0.0556,0.0104,0.0337,0};
    static double [] oldfgfreq = {0,0.0488,0.0032,0.0202,0.0234,0.0276,0.1157,0.0149,0.0191,0.0329,0.0456,
        0.0149,0.1444,0.0308,0.2208, 0.0202,0.1008,0.0297,0.0234,0.0064,0.0573,0};
    static double [] newfginvitro = {0,0.05271,0.00091,0.01203,0.00944,0.02854,0.07964,0.01874,0.0219,0.01599,0.03269,
        0.02539,0.19047,0.07032,0.1889,0.02315,0.10565,0.04302,0.02279,0.0011,0.05662,0};
    static double [] newfgsup35 = {0,0.03542,0.0024,0.01701,0.00815,0.02686,0.08151,0.01544,0.02067,0.01774,
        0.02421,0.02599,0.28863,0.047,0.1301,0.02655,0.11715,0.03943,0.01834,0.00191,0.05552,0};
    static double [] newfgred = {0,0.04865,0.00219,0.01638,0.00783,0.02537,0.07603,0.0181,0.02018,0.01641,0.02639,
        0.02975,0.25885,0.05126,0.15178,0.025,0.10988,0.03841,0.01972,0.00157,0.05624,0};
    static double [] unibg = {0.0,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.0};


    static int [] gc2 = {
        9,    // AAA   K
        12,    // AAC	N
        9,    // AAG   K
        12,    // AAT	N
        17,    // ACA	T
        17,    // ACC	T
        17,    // ACG	T
        17,    // ACT	T
        15,    // AGA	R
        16,    // AGC	S
        15,    // AGG	R
        16,    // AGT	S
        8,    // ATA   I
        8,    // ATC   I
        11,    // ATG	M
        8,    // ATT   I
        14,    // CAA	Q
        7,    // CAC   H
        14,    // CAG	Q
        7,    // CAT   H
        13,    // CCA	P
        13,    // CCC	P
        13,    // CCG	P
        13,    // CCT	P
        15,    // CGA	R
        15,    // CGC	R
        15,    // CGG	R
        15,    // CGT	R
        10,    // CTA	L
        10,    // CTC	L
        10,    // CTG	L
        10,    // CTT	L
        4,    // GAA   E
        3,    // GAC   D
        4,    // GAG   E
        3,    // GAT   D
        1,    // GCA   A
        1,    // GCC   A
        1,    // GCG   A
        1,    // GCT   A
        6,    // GGA   G
        6,    // GGC   G
        6,    // GGG   G
        6,    // GGT   G
        18,   // GTA   V
        18,   // GTC   V
        18,   // GTG   V
        18,   // GTT   V
        21,   // TAA   end
        20,   // TAC   Y
        21,   // TAG   end
        20,   // TAT   Y
        16,   // TCA   S
        16,   // TCC   S
        16,   // TCG   S
        16,   // TCT   S
        21,   // TGA   end
        2,   // TGC   C
        19,   // TGG   W
        2,   // TGT   C
        10,   // TTA   L
        5,   // TTC   F
        10,   // TTG   L
        5    // TTT   F
    };


    spewey() {

        // sys2com = readhashtable("sys2com.txt",null);
        // com2sys = readhashtable("com2sys.txt",null);


        geneticcode = new int[125];
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                for (int k=0; k<4; k++) {
                    int codon1 = 16*i + 4*j + k;
                    int codon2 = 25*(i+1) + 5*(j+1) + (k+1);
                    geneticcode[codon2] = gc2[codon1];
                    scodonnames[codon1] =  new String(new char[] {ntnames[i+1], ntnames[j+1], ntnames[k+1]});
                    // System.out.println(gc2[codon1] + ",\t // " + i + j + k);
                }
            }
        }

        scodonnames[64] = "---";
        for (int i=0; i<5; i++) {
            for (int j=0; j<5; j++) {
                for (int k=0; k<5; k++) {
                    int codon2 = 25*i + 5*j + k;
                    codonnames[codon2] =  new String(new char[] {ntnames[i], ntnames[j], ntnames[k]});
                }
            }
        }


        loglut = new double[loglutlength+1];
        for (int i=0; i<=loglutlength; i++) loglut[i] = Math.log(1.0+ Math.exp(-i/100.0));


        // System.out.println(sum(fgross1) + "\n" + sum(bgross1) + "\n" + sum(fgross2) + "\n" + sum(bgross2));
        // double [] oddsratio1 = new double[22];
        // double [] oddsratio2 = new double[22];
        // for (int k=1; k<=20; k++) {
        //    oddsratio1[k]= (fgross1[k]/(1-fgross1[k]))/(bgross1[k]/(1-bgross1[k]));
        //    oddsratio2[k]= (fgross2[k]/(1-fgross2[k]))/(bgross2[k]/(1-bgross2[k]));
        // }
        // printaausage(oddsratio1);
        // printaausage(oddsratio2);
        // Arrays.sort(oddsratio1);
        // Arrays.sort(oddsratio2);
        // printaausage(oddsratio1);
        // printaausage(oddsratio2);

        int [] aac; double [] aaf;


        int boost = 0; // 10;

        aac = computeaafreq(new String("data/rosspsi1.txt"),false,false);
        // System.out.println(intsum(aac));
        aac[0] = 0; aac[21] = 0; aac[12] = aac[12]+boost;
        aaf = normalize(aac);
        // printaausage(aaf);
        fgross1 = aaf;

        aac = computeaafreq(new String("data/rossnai1.txt"),false,false);
        // System.out.println(intsum(aac));
        aac[0] = 0; aac[21] = 0;
        aaf = normalize(aac);
        // printaausage(aaf);
        bgross1 = aaf;

        aac = computeaafreq(new String("data/rosspsi2.txt"),false,false);
        // System.out.println(intsum(aac));
        aac[0] = 0; aac[21] = 0;  aac[12] = aac[12]+boost;
        aaf = normalize(aac);
        // printaausage(aaf);
        fgross2 = aaf;

        aac = computeaafreq(new String("data/rossnai2.txt"),false,false);
        // System.out.println(intsum(aac));
        aac[0] = 0; aac[21] = 0;
        aaf = normalize(aac);
        // printaausage(aaf);
        bgross2 = aaf;

        lodross1 = new double[22];
        lodross2 = new double[22];
        llkross1 = new double[22];
        llkross2 = new double[22];
        lodross3 = new double[22];
        double [] odross1 = new double[22];
        double [] odross2 = new double[22];


        for (int k=1; k<=20; k++) {
            odross1[k]= ((fgross1[k]/(1-fgross1[k]))/(bgross1[k]/(1-bgross1[k])));
            odross2[k]= ((fgross2[k]/(1-fgross2[k]))/(bgross2[k]/(1-bgross2[k])));
            lodross1[k]= Math.log((fgross1[k]/(1-fgross1[k]))/(bgross1[k]/(1-bgross1[k])));
            lodross2[k]= Math.log((fgross2[k]/(1-fgross2[k]))/(bgross2[k]/(1-bgross2[k])));
            llkross1[k]= Math.log(fgross1[k]/bgross1[k]);
            llkross2[k]= Math.log(fgross2[k]/bgross2[k]);
        }
        // printaausage(odross1);
        // printaausage(odross2);
        // printaausage(llkross1);
        // printaausage(llkross2);


    }


    public static void main (String args[]) {

        if (args.length < 3) {
            // digestdisembl2(args[1],3);
            System.out.println("Usage: ./spewey.jar <fastafile> <corelength> <alpha> <bgfastafile> [<gene>...]");
            System.out.println("\t<fastafile> - containing the sequences to process");
            System.out.println("\t<corelength> - ???");
            System.out.println("\t<alpha> - alpha should be between 0 and 1; it controls the bg frequenicies");
            System.out.println("\t\tIf alpha=1, it uses the bgfreqs from <bgfastafile>, by default this is Scer, ");
            System.out.println("\t\tIf alpha=0, it uses the bgfreqs from all the sequences in fasta-file.");
            System.out.println("\t\tIn general the bgfreqs are alpha*(<bgfastafile>) + (1-alpha)*(sequences in <fastafile>)");
            System.out.println("\t<bgfastafile> - a FASTA file containing AA sequences to use as background frequencies\n\t\t(this argument must be present, but contents are currently ignored: defaults to using Scer frequencies)");
            System.out.println("\t<genes> - an optional list of genes to graph");
            System.out.println("(See source code for further explanation of options)");
            System.exit(1);
        }

        spewey ms = new spewey();

        //double [] fgfreq = oldfgfreq; // old AA frequencies in PrDs
        double [] fgfreq = newfgred;	       // AA frequencies in PrDs, updated based on new hits
        double [] bgscer= normalize(scbg);     // overall AA frequencies in S.cer

        int corelength = Integer.parseInt(args[1]);

        // compute global AA frequencies for sequences being scored
        int [] bgi = computeaafreq(args[0],false,false);
        bgi[0] = 0; bgi[21]=0;
        double [] bgthis = normalize(bgi);

        double alpha = Double.parseDouble(args[2]);

	// FIXME, AKL: read in argument, but not used for the moment, will be implemented
	// in  Oliver King's updated version
	String bgfastafile = args[3];
	System.err.println("#: placeholder only: background FASTA file name:" + bgfastafile);

        if (alpha > 1 || alpha < 0) {
            System.out.println("# warning: invalid alpha; using alpha = 0.5");
            alpha = 0.5;
        }

        // use weighted combination of these as BG freqs, alpha*(s.cer) + (1-alpha)*(this org)

        double [] bgcombo = normalize(axpby(alpha,bgscer,(1-alpha),bgthis));

        //printaausage(bgscer);
        //printaausage(bgthis);
        //printaausage(bgcombo);

        if (args.length==4) {
            // can pass in any fg and bg frequencies here, if these are not suitable
            scorefastas(args[0], corelength, fgfreq, bgcombo);
        }
        if (args.length>=5) {
            plotsomefastas(args[0], corelength, fgfreq, bgcombo, args[4]);
        }
    }

    ///////////////////////

    // whole genome
    static void digestdisembl2(String filename, int numss) {
        int [] aa = new int[1];
        int [][] bv = new int[numss][];
        fastareader fs = new fastareader(filename);
        int count = 0;
        System.out.println("%%\tgene-name\tss1-hw\tss2-hw\tss3-hw\tss1-mx\tss2-mx\tss3-mx\tlen");
        while (fs.hasmorefastas()) {
            String name = fs.nextname();
            StringBuffer sb = fs.nextfasta();
            if (count % numss == 0) aa	= string2aa(sb);
            bv[count % numss] = buffer2bits(sb);
            if (count % numss == (numss-1)) {
                int nameend = name.lastIndexOf("_");
                name = name.substring(0,nameend);
                name = name.trim();
                System.out.print(count/numss + "\t" + name +  "\t");
                for (int i=0; i<numss; i++) {
                    System.out.print(hammingweight(bv[i]) + "\t");
                }
                for (int i=0; i<numss; i++) {
                    System.out.print(longestrun(bv[i]) + "\t");
                }
                System.out.println(aa.length);
            }

            count++;
            // fastaprint(name, aa2string(aa1), 100);
        }
    }


    static int hammingweight(int [] arr) {
        int hw = 0;
        for (int i=0; i<arr.length; i++) if (arr[i] != 0) hw++;
        return hw;
    }



    static int hammingdistance(int [] arr1, int [] arr2) {
        int hd = 0;
        for (int i=0; i<arr1.length; i++) if (arr1[i] != arr2[i]) hd++;
        return hd;
    }


    static void plotsomefastas(String filename,	 int corelength, double [] fgfreq, double [] bgfreq, String genelist) {
        HashMap<String,String> ht1 = new HashMap<String,String>(); // for synonym
        HashMap<String,String> ht2 = new HashMap<String,String>(); // for order
        ht1 = readhashtable(genelist,null);
        ht2 = readhashtable(genelist,"incr");
        boolean doall = false; // plot all genes in fasta file

        double [] fg = normalize(fgfreq);
        double [] bg = normalize(bgfreq);

        double [] llk = new double[22];

        for (int i=1; i<21; i++) llk[i] = Math.log(fg[i]/bg[i]);

        // printaausage(fg);
        // printaausage(bg);
        // printaausage(llk);

        hmm myhmm = prionhmm1(fg,bg);

        int ww = 41;
        int ww2 = 41;
        boolean reflect = false; // true;

        fastareader fs = new fastareader(filename);

        System.out.print("n\tgene\taanum\taa\tvit\tmap\tcharge\thydro\tfindex\tllk\tross1");
        for (int i=0; i<myhmm.numclasses; i++) System.out.print("\t" + myhmm.unames[i]);
        System.out.println();

        int genecount = 1;

        char stopcodon = '*';
        while (fs.hasmorefastas()) {
            String name = fs.nextname();
            String nm = name;
            StringBuffer sb = fs.nextfasta();

            if (doall || ht1.containsKey(name))	 {

                String id = new String("" + genecount);

                if  (ht1.containsKey(name)) {nm = (String) (ht1.get(name));}
                if  (ht2.containsKey(name)) {id = (String) (ht2.get(name));}

                if ((sb.charAt(sb.length()-1))==stopcodon){sb.deleteCharAt(sb.length()-1);} // kill stop codon
                int [] seq = string2aa(sb);
                genecount++;

                // make sure name has no spaces here!!
                name = name.replace(' ',  '-');
                name = name.replace('\t', '-');

                myhmm.decodealls(seq);
                // myhmm.printme(genecount + "\t" + name, seq);
                // myhmm.sposteriorl(seq);
                // myhmm.sviterbidecodel(seq);
                disorderreport dr = new disorderreport(seq,ww,ww2,reflect, new double [] {2.78 , -1  ,-1.15}, llk, lodross1);
                //	dr.printme(genecount + "\t" + name);
                for (int i=0; i<seq.length; i++) {
                    System.out.print(id + "\t" + nm + "\t" + (i+1) +"\t" + seq[i] +  "\t" + myhmm.viterbipath[i] + "\t" + myhmm.mappath[i]+"\t"
                            + (float) dr.charge[i] + "\t" + (float) dr.hydro[i] + "\t" + (float) dr.combo[i] + "\t" + (float) dr.other[i] + "\t" + (float) dr.other2[i]);
                    for (int j=0; j<myhmm.postprob.length; j++){ System.out.print("\t" + (float) (myhmm.postprob[j][i]));

                    };
                    System.out.println();
                }

                System.out.println("####### "+ myhmm.lmarginalprob + "\t" + myhmm.lviterbiprob);

            }
        }
    }



    //	MW   vs	  HG   vs   LL	 vs   HMM
    static void scorefastas(String filename,  int corelength, double [] fgfreq, double [] bgfreq) {

        int chomp = 0;

        // int corelength = 50;
        // int targtm = 60;

        double [] fg = normalize(fgfreq);
        double [] bg = normalize(bgfreq);

        double [] llk = new double[22];

        for (int i=1; i<21; i++) llk[i] = Math.log(fg[i]/bg[i]);

        // printaausage(fg);
        // printaausage(bg);
        // printaausage(llk);

        System.out.print("geneID\tMW\tMWstart\tMWend\tMWlen\tLLK\tLLKstart\tLLKend\t");
        System.out.print("LLKlen\tNLLK\tCOREscore\tCOREstart\tCOREend\tCORElen\tPRDstart\tPRDend\tPRDlen\tPROTlen\t");
        System.out.print("HMMmap\tHMMvit\tCOREaa\tSTARTaa\tENDaa\tPRDaa\tPRDscore\tFInumaa\tFImeanhydro\tFImeancharge\tFImeancombo\tFImaxrun\t");
        System.out.print("ROSS1prd\tROSS1fi\tROSS1cen\tROSS1aa\tROSS2prd\tROSS2fi\tROSS2cen\tROSS2aa");
        for (int aac = 1; aac <= 20 ; aac++) {
            System.out.print("\t" + aanames[aac]);
        }
        System.out.println();
        // Hashtable ht = readhashtable(fglist,"1");

        fastareader fs = new fastareader(filename);
        int maxlen = 500;
        hgalg hg = new hgalg(maxlen,bg);
        int genecount=0;

        double [] qnmask = {0,0,0,0,0,0,0,0,0,0,0,0,1.0,0,1.0,0,0,0,0,0,0,0};

        double [][] lops;
        double [] maa1;
        double [] maa2;
        double [] maa3;
        double [] hs1;
        double [] hs2;
        double [] hs3;
        int [] dna;
        int [] rcdna;
        int [] aa;
        String flag;


        int [] aacomp;

        double hmmscore;
        double hmmscorev;

        int hggood = 0;
        int hggood2 = 0;

        hmm fghmm = prionhmm1(fg,bg);
        hmm bghmm = prionhmm111(bg);

        // int ww = 51;
        int ww = 41;
        // int ww2 = 51;
        int ww2 = 41;
        boolean reflect = false; // true;
        char stopcodon = '*';


        while (fs.hasmorefastas()) {
            String name = fs.nextname();
            StringBuffer sb = fs.nextfasta();
            if ((sb.charAt(sb.length()-1))==stopcodon){sb.deleteCharAt(sb.length()-1);} // kill stop codon
            // dna = string2dna(sb);
            aa = string2aa(sb);

            aacomp = aacomposition(aa);

            if (aa.length<12) continue;

            // MW
            maa1 = mapseq(aa,qnmask);
            hs1 = hss2(maa1,80,80);

            // LLK
            maa2 = mapseq(aa,llk);
            hs2 = hss2(maa2,corelength,corelength); // max 500?

            // HMM
            fghmm.decodeall(aa);
            bghmm.decodeall(aa);

            hmmscore = fghmm.lmarginalprob - bghmm.lmarginalprob;
            hmmscorev = fghmm.lviterbiprob - bghmm.lviterbiprob;

            // double [] mixlods = axpby(0.5,lodross1,0.5,lodross1); mixlods[12] = mixlods[12]+1.6;

            //disorderreport  dr = new disorderreport(aa,ww,ww2,reflect, new double [] {2.78 , -1	 ,-1.15}, llk, lodross1);
            // dr.printme();

            //disorderreport dr2 = new disorderreport(aa,ww,ww2,reflect, new double [] {2.78 , -1	 ,-1.15}, llk, lodross2);

            String nm = name.substring(chomp);

            // one-time-use --- change this!!!
            //	 int chopdex = nm.indexOf(" ",	nm.indexOf(" ")+1);
            // nm = nm.substring(0, chopdex);

            int [] seg;
            int [] ppp;

            int [] mp = fghmm.viterbipath;
            // check each map segment for hss;
            maa3 = mapseq(aa,llk);
            for (int i=0; i<maa3.length; i++) {
                if (mp[i]==0) maa3[i] = -100000000;
            }
            hs3 = hss2(maa3,corelength,corelength);
            hs3[2] = Math.max(hs3[2],0.0);
            // expand up;
            // expand down;
            int aastart = (int) hs3[0];
            int aastop	= (int) hs3[1];
            while(aastart>=0 && mp[aastart]==1) aastart--; aastart++;
            while(aastop<mp.length && mp[aastop]==1) aastop++; aastop--;

            int [] prd;
            double prdscore = 0;
            int numcys = 0;

            prd = submatrix(aa, aastart, aastop);
            // score of whole prd?
            for (int kk = 0; kk<prd.length; kk++) {
                prdscore = prdscore + llk[prd[kk]];
                if (prd[kk]==2) numcys++;
            }

            //disorderreport  dr = new disorderreport(aa,ww,ww2,reflect, new double [] {2.78 , -1	 ,-1.15}, llk, lodross1);
            disorderreport  dr = new disorderreport(prd,ww,ww2,reflect, new double [] {2.78 , -1	 ,-1.15}, llk, lodross1);

            //disorderreport dr2 = new disorderreport(aa,ww,ww2,reflect, new double [] {2.78 , -1	 ,-1.15}, llk, lodross2);
            disorderreport dr2 = new disorderreport(prd,ww,ww2,reflect, new double [] {2.78 , -1	 ,-1.15}, llk, lodross2);

            boolean fastaprds = false;
            if (fastaprds) {
                if (prdscore > 0) {
                    System.out.println(">" + nm + "-pprd");
                    System.out.println(aa2string(prd));
                    System.out.println(">" + nm + "-core");
                    System.out.println(aa2string(submatrix(aa, (int) hs3[0], (int) (hs3[1])) ));
                }
            } else{

                if (hs3[2]>-1) {
                    // use function for roundoff for rounding insread!
                    System.out.print(nm + "\t"
                            + (int) hs1[2]   + "\t" + (int) hs1[0] + "\t" + (int) hs1[1] + "\t" + (int) (hs1[1]-hs1[0]+1) + "\t"
                            + (float) (Math.round(1000*hs2[2])/1000.0) + "\t" + (int) hs2[0] + "\t" + (int) hs2[1] + "\t" + (int) (hs2[1]-hs2[0]+1) + "\t"
                            + (float) (Math.round(1000*(hs2[2]/(hs2[1]-hs2[0]+1)))/1000.0) + "\t"
                            + (float) (Math.round(1000*hs3[2])/1000.0)	  + "\t" + (int) hs3[0] + "\t" + (int) hs3[1]
                            + "\t" + (int) (hs3[1]-hs3[0]+1) + "\t"
                            +	aastart + "\t" + aastop + "\t" + (aastop-aastart+1) + "\t" + aa.length
                            + "\t" + (float) (Math.round(1000*hmmscore)/1000.0) + "\t" + (float) (Math.round(1000*hmmscorev)/1000.0) + "\t");
                    System.out.print(aa2string(submatrix(aa, (int) hs3[0], (int) (hs3[1])) ));
                    System.out.print("\t");
                    System.out.print(aa2string(submatrix(aa, aastart, aastart+14)));
                    System.out.print("\t");
                    System.out.print(aa2string(submatrix(aa, aastop-14, aastop)));
                    System.out.print("\t");
                    System.out.print(aa2string(prd));
                    System.out.print("\t" + (float) (Math.round(1000*prdscore)/1000.0));
                    System.out.print("\t" + dr.numdisorderedstrict2 + "\t" + (float) (Math.round(1000*dr.meanhydro)/1000.0) +"\t" + (float) (Math.round(1000*dr.meancharge)/1000.0) + "\t"
                            + (float) (Math.round(1000*dr.meancombo)/1000.0) + "\t" + (int) (dr.maxlen));
                    // System.out.print("\t" + (float) (dr.rossmaxprd) +  "\t" + (float) (dr.rossmaxdis) +	 "\t" + (int) (dr.rossmaxcenter));
                    // System.out.print("\t" + aa2string(submatrix(aa,  dr.rossmaxcenter - ww2/2, dr.rossmaxcenter + ww2/2))); // ok for length < w?
                    System.out.print("\t" + (float) (dr.rossmaxprd) +  "\t" + (float) (dr.rossmaxdis) +	 "\t" + (int) (aastart + dr.rossmaxcenter));
                    System.out.print("\t" + aa2string(submatrix(prd,  dr.rossmaxcenter - ww2/2, dr.rossmaxcenter + ww2/2))); // ok for length < w?
                    // System.out.print("\t" + (float) (dr2.rossmaxprd) + "\t" + (float) (dr2.rossmaxdis) + "\t" + (int) (dr2.rossmaxcenter));
                    // System.out.print("\t" + aa2string(submatrix(aa,  dr2.rossmaxcenter - ww2/2, dr2.rossmaxcenter + ww2/2)));
                    System.out.print("\t" + (float) (dr2.rossmaxprd) + "\t" + (float) (dr2.rossmaxdis) + "\t" + (int) (aastart + dr2.rossmaxcenter));
                    System.out.print("\t" + aa2string(submatrix(prd,  dr2.rossmaxcenter - ww2/2, dr2.rossmaxcenter + ww2/2)));
                    for (int aac = 1; aac <= 20 ; aac++) {
                        System.out.print("\t" + aacomp[aac]);
                    }
                    System.out.println();
                }
            }
        }

    }




    static int [] rc(int [] dna) {
        int n = dna.length;
        int [] rcdna = new int[n];
        for (int i=0; i<n; i++) rcdna[i] = 5-dna[n-1-i];
        return rcdna;

    }

    // static String translate(Hashtable ht, String key) {
    //	if (ht.containsKey(key)) {return (String) ht.get(key);}
    //	else return key;
    // }


    static String translate(HashMap<String,String> ht, String key) {
        if (ht.containsKey(key)) {return (String) ht.get(key);}
        else return key;
    }

    static int [][] trimrows(int [][] mat, int r1, int r2) {
        int [][] newmat = new int [r2-r1+1][];
        for (int i=0; i<r2-r1+1; i++) newmat[i] = mat[i+r1];
        return newmat;

    }


    static float roundoff(double number, int digits) {
        return (float) (Math.round(Math.pow(10,digits)*number)/Math.pow(10,digits));

    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////	       hmms				/////////////////////////



    static hmm prionhmm1(double [] fgfreq, double bgfreq []) {
        // state 0 = normal; state 1 = prion
        double [][] tmat = {{99.9/100, 0.1/100},{2.0/100, 98.0/100}};
        double [] imat = {0.9524, 0.0476};
        double [][] emat = new double[2][];
        emat[0] = normalize(axpb(1,bgfreq,0.00001));
        emat[1] = normalize(axpb(1,fgfreq,0.00001));
        hmm h = new hmm(tmat,emat,imat);
        h.subtrellis = new int [] {1,0};
        h.states = new char[] {'-','+'};
        h.setnames(new String [] {"background","PrD-like"});
        return h;

    }


    static hmm prionhmm111(double bgfreq []) {
        // state 0 = normal; state 1 = prion
        double [][] tmat = {{1, 0},{0, 1}};
        double [] imat = {1, 0};
        double [][] emat = new double[2][];
        emat[0] = normalize(axpb(1,bgfreq,0.00001));
        emat[1] = normalize(axpb(1,bgfreq,0.00001));
        hmm h = new hmm(tmat,emat,imat);
        h.subtrellis = new int [] {1,0};
        h.states = new char[] {'-','+'};
        h.setnames(new String [] {"background","PrD-like"});
        return h;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////


    // computes log(exp(a) + exp(b));
    // okay if a or b is -Infinity, or both
    // a bit slow --- precompute lookup table?
    static double logeapeb2(double a, double b) {
        if (a>b) return (a + Math.log(1+Math.exp(b-a)));
        else if (b>a) return (b + Math.log(1+Math.exp(a-b)));
        else return (a + log2); // takes care of a = b = -Infinity

    }


    static double logeapeb(double a, double b) {
        if (a>b) {
            double c = a-b;
            if (!(c<40)) return a; // takes care of  b = -Inf
            else {
                int dex = (int) Math.floor(100*c);
                //System.out.print("A " + a + "\t" + "\t" + b + "\t" + c + "\t" + dex + "\t"  + loglut[dex+1] + "\t" +	loglut[dex] + "\t");
                return (a + ((100*c-dex)*loglut[dex+1] + (dex+1-100*c)*loglut[dex]));

            }

        }
        else if (b>a) {
            double c = b-a;
            if (!(c<40)) return b; // takes care of  a = -Inf
            else {
                int dex = (int) Math.floor(100*c);
                // System.out.print("B " + a + "\t" + "\t" + b + "\t" + c + "\t" + dex + "\t"  + loglut[dex+1] + "\t" +	 loglut[dex] + "\t");
                return (b + ((100*c-dex)*loglut[dex+1] + (dex+1-100*c)*loglut[dex]));
            }
        }
        else return (a + log2); // takes care of a = b = -Infinity

    }


    static double [] killstops(double [] pi0) {
        double [] pi = copymatrix(pi0);
        pi[48] = 0;
        pi[50] = 0;
        pi[56] = 0;
        // pi = normalize(pi);
        return pi;
    }

    static double [][] killstops(double [][] jscod0) {
        double [][] jscod = copymatrix(jscod0);
        for (int i=0; i<64; i++) {
            jscod[48][i] = 0;
            jscod[50][i] = 0;
            jscod[56][i] = 0;
            jscod[i][48] = 0;
            jscod[i][50] = 0;
            jscod[i][56] = 0;
        }
        return jscod;
    }


    static double [] asdouble(int [] arr) {
        int n = arr.length;
        // System.out.println(n);
        double [] narr = new double[n];
        for (int i=0; i<n; i++) narr[i] = arr[i]+0.0;
        return narr;
    }

    static double [][] asdouble(int [][] arr) {
        int n = arr.length;
        int m = arr[0].length;
        // System.out.println(n);
        double [][] narr = new double[n][m];
        for (int i=0; i<n; i++) for (int j=0; j<m; j++) narr[i][j] = arr[i][j]+0.0;
        return narr;
    }

    // highest-scoring-subsequence
    static double [] hss(double [] seq, int minlength, int maxlength) {
        int n = seq.length;
        if (minlength > n) minlength = n;
        if (maxlength > n) maxlength = n;
        double [] score = new double[3]; // start, end, score
        double bestscore;
        double tempscore;
        double [] psum = new double[n+1];
        int beststart = 0;
        int beststop = minlength;
        psum[0] = 0;
        for (int i=0; i<n; i++) {
            psum[i+1] = psum[i]+seq[i];
        }

        bestscore = psum[minlength]-psum[0];

        for (int i=0; i<n-minlength; i++) {
            for (int j=i+minlength; j<= Math.min(i+maxlength,n); j++) {
                tempscore = psum[j] - psum[i];
                if (tempscore > bestscore) {
                    bestscore = tempscore;
                    beststart = i;
                    beststop = j-1;
                }
            }
        }
        score[0] = beststart;
        score[1] = beststop;
        score[2] = bestscore;
        return score;
    }


    // highest-scoring-subsequence
    // assumes at least one element is positive
    static double [] hss2(double [] seq) {
        int n = seq.length;

        double [] score = new double[3]; // start, end, score
        double bestscore;
        int beststart = 0;
        int beststop = 0;
        int curstart = 0;

        double d = Math.max(seq[0],0);
        bestscore = d;

        for (int i=1; i<n; i++) {
            if (d+seq[i] > 0) {
                d = d + seq[i];
            }
            else {
                d = 0;
                curstart = i+1;
            }
            if (d > bestscore) {
                bestscore = d;
                beststop = i;
                beststart = curstart;
            }
        }

        score[0] = beststart;
        score[1] = beststop;
        score[2] = bestscore;
        return score;
    }



    // highest-scoring-subsequence
    // assumes at least one element is positive
    static double [] hss2(double [] seq, int minlen) {
        int minlength = minlen-1;
        int n = seq.length;
        if (minlength >= n) minlength = n-1;
        double [] score = new double[3]; // start, end, score
        double bestscore;
        int beststart = 0;
        int beststop = minlength;
        int curstart = 0;
        int newstart = 0;

        double [] psum = new double[n+1];

        psum[0] = 0;
        for (int i=0; i<n; i++) {
            psum[i+1] = psum[i]+seq[i];
        }

        double d = psum[minlength+1];
        bestscore = d;

        for (int i = minlength+1; i<n; i++) {
            d = psum[i+1] - psum[curstart];
            newstart = curstart;
            for (int j=curstart+1; j<i-minlength+1; j++) {
                if (psum[i+1]-psum[j] >= d) {
                    d = psum[i+1]-psum[j];
                    newstart = j;
                }
                curstart = newstart;
            }
            if (d > bestscore) {
                bestscore = d;
                beststop = i;
                beststart = curstart;
            }
        }

        score[0] = beststart;
        score[1] = beststop;
        score[2] = bestscore;
        return score;
    }


    // highest-scoring-subsequence
    // assumes at least one element is positive
    static double [] hss2(double [] seq, int minlen, int maxlength) {
        int minlength = minlen-1;
        int n = seq.length;
        if (minlength >= n) minlength = n-1;  //
        if (maxlength > n) maxlength = n;  //  +1?
        double [] score = new double[3];   // start, end, score
        double bestscore;
        int beststart = 0;
        int beststop = minlength;
        int curstart = 0;
        int newstart = 0;

        double [] psum = new double[n+1];

        psum[0] = 0;
        for (int i=0; i<n; i++) {
            psum[i+1] = psum[i]+seq[i];
        }

        double d = psum[minlength+1];
        bestscore = d;

        for (int i = minlength+1; i<n; i++) {
            if ((i-curstart) >= maxlength) curstart++;
            d = psum[i+1] - psum[curstart];
            newstart = curstart;
            for (int j=curstart+1; j<i-minlength+1; j++) {
                if (psum[i+1]-psum[j] >= d) {
                    d = psum[i+1]-psum[j];
                    newstart = j;
                }
                curstart = newstart;
            }
            if (d > bestscore) {
                bestscore = d;
                beststop = i;
                beststart = curstart;
            }
        }

        score[0] = beststart;
        score[1] = beststop;
        score[2] = bestscore;
        return score;
    }



    static int [] firstandlast(int p1, int p2, int [] vec, int skipme) {
        int [] fal = new int[2];
        int dex = 0;
        int hits = 0;
        while ((dex <vec.length) && (hits <= p1)) {
            if (vec[dex] != skipme) hits++;
            dex++;
        }
        fal[0]=dex-1;

        while ((dex <vec.length) && (hits <= p2)) {
            if (vec[dex] != skipme) hits++;
            dex++;
        }
        fal[1]=dex-1;
        return fal;
    }



    static double [][] transpose(double [][] mat) {
        int r = mat.length;
        int c = mat[0].length;
        double [][] tmat = new double[c][r];
        for (int i=0; i<r; i++) {
            for (int j=0; j<c; j++) {
                tmat[j][i] = mat[i][j];
            }
        }
        return tmat;
    }


    static int [][] transpose(int [][] mat) {
        int r = mat.length;
        int c = mat[0].length;
        int [][] tmat = new int[c][r];
        for (int i=0; i<r; i++) {
            for (int j=0; j<c; j++) {
                tmat[j][i] = mat[i][j];
            }
        }
        return tmat;
    }


    static void charcount(String filename) {
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String line;
            int [] counts = new int[256];


            while ((line = in.readLine()) != null) {
                line.trim();
                if (line != null) {
                    for (int i=0; i<line.length(); i++) counts[line.charAt(i)]++;
                }

            }
            for (int j=65; j<91; j++) System.out.println((char) j + "\t" + counts[j] + "\t\t" + (char) (j+32) + "\t" + counts[j+32]);
            System.out.println(intsum(submatrix(counts,65,90)));
            System.out.println(intsum(submatrix(counts,65+32,90+32)));

        }
        catch (IOException e) {
            System.out.println("# Couldn't open");
        }


    }


    static int [] findindices(int [] vec, int target) {
        int n = vec.length;
        int hits = 0;
        for (int i=0; i<n; i++) if (vec[i]==target) hits++;
        int [] dexvec = new int[hits];
        hits = 0;
        for (int i=0; i<n; i++) if (vec[i]==target) {dexvec[hits]=i; hits++;}
        return dexvec;
    }


    static double [] scodf2aaf(double [] scodf) {
        double [] aaf = new double[22];
        for (int i=0; i<64; i++) aaf[gc2[i]] += scodf[i];
        return aaf;
    }



    // e.g. "QN.GY" makes three groups, one with QN, one with GY, one with everything else
    static int [] string2mask(String clumps) {
        int [] mask = new int[22];
        int seg = 1;
        for (int i=0; i<clumps.length(); i++) {
            char c = clumps.charAt(i);
            if (c=='.') seg++;
            else {
                for (int k=0; k<22; k++) {
                    if (aanames[k]==c) {mask[k] = seg; break;}
                }
            }
        }
        return mask;
    }


    // partial sums of rows
    static double [][] cumsum(double [][] mat, boolean byrow) {
        int r = mat.length;
        int c = mat[0].length;
        double [][] cmat = new double[r][c];
        if (byrow) {
            for (int i=0; i<r; i++) {
                double rs = 0;
                for (int j=0; j<c; j++) {
                    rs = rs + mat[i][j];
                    cmat[i][j] = rs;
                }
            }
        }
        else {
            for (int i=0; i<c; i++) {
                double rs = 0;
                for (int j=0; j<r; j++) {
                    rs = rs + mat[j][i];
                    cmat[j][i] = rs;
                }
            }

        }
        return cmat;
    }


    // partial sums of rows
    static int [][] cumsum(int [][] mat, boolean byrow) {
        int r = mat.length;
        int c = mat[0].length;
        int [][] cmat = new int[r][c];
        if (byrow) {
            for (int i=0; i<r; i++) {
                int rs = 0;
                for (int j=0; j<c; j++) {
                    rs = rs + mat[i][j];
                    cmat[i][j] = rs;
                }
            }
        }
        else {
            for (int i=0; i<c; i++) {
                int rs = 0;
                for (int j=0; j<r; j++) {
                    rs = rs + mat[j][i];
                    cmat[j][i] = rs;
                }
            }

        }
        return cmat;
    }


    static int [][] fliplr (int [][] mat) {
        int m = mat.length;
        int n = mat[0].length;
        int [][] mat2 = new int [m][n];
        for (int i=0; i<m; i++) for (int j=0; j<n; j++) mat2[i][j] = mat[i][n-j-1];
        return mat2;


    }


    // partial sum
    static double [] cumsum(double [] arr) {
        int r = arr.length;
        double [] carr = new double[r];
        double cs = 0;
        for (int i=0; i<r; i++) {
            cs = cs + arr[i];
            carr[i]=cs;
        }
        return carr;
    }


    // converts 0-63 to 0-124
    static int scod2cod(int scod) {
        int a1;
        int a2;
        int a3;
        a3 = scod % 4;
        scod = (scod - a3)/4;
        a2 = scod % 4;
        scod = (scod - a2)/4;
        a1 = scod % 4;
        return 25*(a1+1) + 5*(a2+1) + (a3+1);
    }

    // converts 0-124 to 0-63 ---- -1 if not valid codon
    static int cod2scod(int cod) {
        int a1;
        int a2;
        int a3;
        a3 = cod % 5;
        cod = (cod - a3)/5;
        a2 = cod % 5;
        cod = (cod - a2)/5;
        a1 = cod % 5;
        if (a1==0 || a2==0 || a3==0 ) return -1;
        else return 16*(a1-1) + 4*(a2-1) + (a3-1);

    }

    static int [] codon2dna(int [] codon) {
        int m = codon.length;
        int n = 3*m;
        int [] dna = new int[n];
        int a1;
        int a2;
        int a3;
        int cod;
        for (int i=0; i<m; i++) {
            cod = codon[i];
            a3 = cod % 5;
            cod = (cod - a3)/5;
            a2 = cod % 5;
            cod = (cod - a2)/5;
            a1 = cod % 5;
            dna[3*i] = a1;
            dna[3*i+1] = a2;
            dna[3*i+2] = a3;
        }
        return dna;
    }


    static int [] codon2dna(int codon) {
        int m = 1;
        int n = 3*m;
        int [] dna = new int[n];
        int a1;
        int a2;
        int a3;
        int cod;
        for (int i=0; i<m; i++) {
            cod = codon;
            a3 = cod % 5;
            cod = (cod - a3)/5;
            a2 = cod % 5;
            cod = (cod - a2)/5;
            a1 = cod % 5;
            dna[3*i] = a1;
            dna[3*i+1] = a2;
            dna[3*i+2] = a3;
        }
        return dna;
    }


    static int [] scodon2dna(int scodon) {
        return codon2dna(scod2cod(scodon));
    }

    static int []  scod2cod(int [] scod) {
        int [] cod = new int[scod.length];
        for (int i=0; i<scod.length; i++) cod[i] = scod2cod(scod[i]);
        return cod;
    }

    static int []  cod2scod(int [] cod) {
        int [] scod = new int[cod.length];
        for (int i=0; i<cod.length; i++) scod[i] = cod2scod(cod[i]);
        return scod;
    }

    static int [][] submatrix(int [][] matrix, int r1, int r2, int c1, int c2) {
        int m = matrix.length;
        int n = matrix[0].length;
        if (r1<0) r1=0;
        if (r2<r1) r2=r1;
        if (c1<0) c1=0;
        if (c2<c1) c2=c1;
        if (r1>=m) r1=m-1;
        if (r2>=m) r2=m-1;
        if (c1>=n) c1=n-1;
        if (c2>=n) c2=n-1;
        int [][] newmat = new int[r2-r1+1][c2-c1+1];
        for (int i=0; i<r2-r1+1; i++)
            for (int j=0; j<c2-c1+1;j++)
                newmat[i][j] = matrix[r1+i][c1+j];

        return newmat;
    }


    static double  [][] submatrix(double [][] matrix, int r1, int r2, int c1, int c2) {
        int m = matrix.length;
        int n = matrix[0].length;
        if (r1<0) r1=0;
        if (r2<r1) r2=r1;
        if (c1<0) c1=0;
        if (c2<c1) c2=c1;
        if (r1>=m) r1=m-1;
        if (r2>=m) r2=m-1;
        if (c1>=n) c1=n-1;
        if (c2>=n) c2=n-1;
        double [][] newmat = new double[r2-r1+1][c2-c1+1];
        for (int i=0; i<r2-r1+1; i++)
            for (int j=0; j<c2-c1+1;j++)
                newmat[i][j] = matrix[r1+i][c1+j];

        return newmat;
    }

    static int [] submatrix(int [] matrix, int r1, int r2) {
        int m = matrix.length;

        if (r1<0) r1=0;
        if (r2<r1) r2=r1;
        if (r1>=m) r1=m-1;
        if (r2>=m) r2=m-1;

        int [] newmat = new int[r2-r1+1];
        for (int i=0; i<r2-r1+1; i++) newmat[i] = matrix[r1+i];
        return newmat;
    }


    static double [] submatrix(double [] matrix, int r1, int r2) {
        int m = matrix.length;

        if (r1<0) r1=0;
        if (r2<r1) r2=r1;
        if (r1>=m) r1=m-1;
        if (r2>=m) r2=m-1;

        double [] newmat = new double[r2-r1+1];
        for (int i=0; i<r2-r1+1; i++) newmat[i] = matrix[r1+i];
        return newmat;
    }


    // kullback-leibler divergence
    static double kld(double ps[], double qs[]) {
        double h = 0.0;
        double eps = 0.000000001;
        for (int i=0; i<ps.length; i++) {
            if (ps[i] > eps || qs[i] > eps) {
                h = h + ps[i]*plogp(ps[i]/qs[i]);
            }
        }
        return h;
    }



    static double entropy(double [] probs) {
        double h = 0;
        for (int i=0; i<probs.length; i++) h = h - plogp(probs[i]);
        return h;
    }


    static double entropy(double [][] probs) {
        return entropy(flatten(probs));
    }



    // p*log(p)
    static double plogp(double p) {
        if (p==0) return 0; else return p*Math.log(p); // change to base 2?
    }


    static void killlowercase(StringBuffer sb) {
        for (int j=0; j<sb.length(); j++) {
            if (Character.isLowerCase(sb.charAt(j))) sb.setCharAt(j,'-');
        }

    }

    static void printmatrix(double [] mat) {
        for (int i=0; i<mat.length; i++) {
            System.out.println((float) (mat[i]));
        }
        System.out.println();
    }


    static int basetoint(char c) {
        if (c=='A' || c=='a')  return 1;
        else if (c=='C' || c=='c') return 2;
        else if (c=='G' || c=='g') return 3;
        else if (c=='T' || c=='t' || c=='U' || c=='u') return 4;
        else {
            //System.out.println("bogus base " + c);
            return 0;
        }
    }


    static int aatoint(char c) {
        if (c=='A' || c=='a')  return 1;
        else if (c=='C' || c=='c') return 2;
        else if (c=='D' || c=='d') return 3;
        else if (c=='E' || c=='e') return 4;
        else if (c=='F' || c=='f') return 5;
        else if (c=='G' || c=='g') return 6;
        else if (c=='H' || c=='h') return 7;
        else if (c=='I' || c=='i') return 8;
        else if (c=='K' || c=='k') return 9;
        else if (c=='L' || c=='l') return 10;
        else if (c=='M' || c=='m') return 11;
        else if (c=='N' || c=='n') return 12;
        else if (c=='P' || c=='p') return 13;
        else if (c=='Q' || c=='q') return 14;
        else if (c=='R' || c=='r') return 15;
        else if (c=='S' || c=='s') return 16;
        else if (c=='T' || c=='t') return 17;
        else if (c=='V' || c=='v') return 18;
        else if (c=='W' || c=='w') return 19;
        else if (c=='Y' || c=='y') return 20;
        else if (c=='X' || c=='*') return 21; // X could mean something else!
        else {
            //System.out.println("bogus base " + c);
            return 0;
        }
    }



    // use binary search instead?
    static int char2int(char c, char [] arr) {
        for (int i=0; i<arr[i]; i++) {if (arr[i]==c) return i;}
        return -1;

    }


    static int [] string2ints(StringBuffer sb, char [] chararr) {
        int n = sb.length();
        int [] arr = new int[n];
        for (int i=0; i<n; i++) arr[i]=char2int(sb.charAt(i),chararr);
        return arr;

    }


    static int [][] basestoints (String [] sites) {
        int r = sites.length;
        int c = sites[0].length();
        int [][] imat = new int[r][c];
        for (int i=0; i<r; i++) {
            for (int j=0; j<c; j++) {
                imat[i][j] = basetoint(sites[i].charAt(j));
            }
        }
        return imat;
    }


    static double sum(int [] vec) {
        double sm = 0;
        for (int i=0; i<vec.length; i++) sm = sm + vec[i];
        return sm;
    }

    static int intsum(int [] vec) {
        int sm = 0;
        for (int i=0; i<vec.length; i++) sm = sm + vec[i];
        return sm;
    }


    static double sum(double [] vec) {
        double sm = 0;
        for (int i=0; i<vec.length; i++) sm = sm + vec[i];
        return sm;
    }


    static double mean(int [] vec) {
        double mn = 0;
        for (int i=0; i<vec.length; i++) mn = mn + vec[i];
        return ((1.0*mn)/vec.length);
    }


    static double mean(double [] vec) {
        double mn = 0;
        for (int i=0; i<vec.length; i++) mn = mn + vec[i];
        return ((1.0*mn)/vec.length);
    }


    static double variance(double [] vec) {
        double mn = mean(vec);
        double var = 0;
        for (int i=0; i<vec.length; i++) var = var + (vec[i]-mn)*(vec[i]-mn);
        return (1.0*var)/vec.length;
    }

    static double variance(int [] vec) {
        double mn = mean(vec);
        double var = 0;
        for (int i=0; i<vec.length; i++) var = var + (vec[i]-mn)*(vec[i]-mn);
        return (1.0*var)/vec.length;
    }



    // mean, var, median, alpha-trimmed-mean, alpha-trimmed-var
    static double [] arraystats(double [] arr, double alpha) {
        int n = arr.length;
        double [] stats = new double[5];
        stats[0] = mean(arr);
        stats[1] = variance(arr);
        double [] newarr = new double[n];
        System.arraycopy(arr,0,newarr,0,n);
        Arrays.sort(newarr);
        // median
        if ( n%2 == 0) stats[2] = (newarr[n/2] + newarr[n/2-1])/2.0;
        else stats[2] = newarr[(n-1)/2];
        int numtrim = (int) Math.round(alpha*n);
        newarr = submatrix(newarr,numtrim,n-numtrim-1); // +-1?
        stats[3] = mean(newarr);
        stats[4] = variance(newarr);
        return stats;

    }

    static double std(double [] vec) {
        return Math.sqrt(variance(vec));
    }

    static double std(int [] vec) {
        return Math.sqrt(variance(vec));
    }


    static double median(int [] vec) {
        double mn = 0;
        int n = vec.length;
        double [] svec = new double[n];
        for (int i=0; i<n; i++) svec[i] = vec[i];
        java.util.Arrays.sort(svec);
        if (n%2 == 0) return((0.0 + svec[n/2] + svec[n/2-1])/2);
        else return (svec[(n-1)/2]);
    }

    // do this with logs?
    static double binomial(int n, int k) {
        double mul = 1.0;
        int top = n;
        for (int i=1; i<=k; i++) {mul = (mul*top)/k; top--;}
        return mul;
    }


    static void translate(String filename) {
        int [] dna;
        int [] aa;
        fastareader fs = new fastareader(filename);
        while (fs.hasmorefastas()) {
            String name = fs.nextname();
            StringBuffer sb = fs.nextfasta();
            dna	  = string2dna(sb);
            aa	= dna2aa(dna);
            fastaprint(name, aa2string(aa), 100);
        }
    }


    // discretize to 10000 parts? 0.0001 // uppermost cprob should be 1, so index stays in bounds
    static int randomint(Random rg, double [] cprobs) {
        int i = 0;
        double d = rg.nextDouble();
        while (d>cprobs[i]) {
            if (i==21) {System.out.println(d); printmatrix(cprobs);}
            i++;}

            return i;
    }


    static int [] computeaafreq(String filename, boolean isdna, boolean killlc) {
        int [] aacounts = new int[22];
        fastareader fs = new fastareader(filename);
        while (fs.hasmorefastas()) {
            String geneid = fs.nextname();
            StringBuffer sb = fs.nextfasta();
            if (killlc) killlowercase(sb);
            countaas(sb,aacounts, isdna);
        }
        // printmatrix(aacounts);
        return aacounts;
    }


    static int [][] computeaafreq2(String filename, boolean isdna, boolean killlc) {
        int [][] aacounts = new int[22][22];
        fastareader fs = new fastareader(filename);
        while (fs.hasmorefastas()) {
            String geneid = fs.nextname();
            StringBuffer sb = fs.nextfasta();
            if (killlc) killlowercase(sb);
            countaas2(sb, aacounts, isdna);
        }
        // printmatrix(aacounts);
        return aacounts;
    }


    static int [][][] computeaafreq3(String filename, boolean isdna, boolean killlc) {
        int [][][] aacounts = new int[22][22][22];
        fastareader fs = new fastareader(filename);
        while (fs.hasmorefastas()) {
            String geneid = fs.nextname();
            StringBuffer sb = fs.nextfasta();
            if (killlc) killlowercase(sb);
            countaas3(sb, aacounts, isdna);
        }
        // printmatrix(aacounts);
        return aacounts;
    }

    // orphan?
    static int [] computeaafreq(int [][] aas) {
        int [] aacounts = new int[22];
        for (int i=0; i<aas.length; i++) {
            //
        }
        // printmatrix(aacounts);
        return aacounts;
    }


    static void countaas(StringBuffer buffer, int [] aacounts, boolean isdna) {
        int [] aa;
        if (isdna) {
            int [] dna	 = string2dna(buffer);
            int [] codon = dna2codon(dna);
            aa	= codon2aa(codon);
        }
        else aa = string2aa(buffer);
        if (isvalidprotein(aa)) {
            for (int i=0; i<aa.length; i++) {
                aacounts[aa[i]]++;
            }
        }
    }

    static void countaas2(StringBuffer buffer, int [][] aacounts, boolean isdna) {
        int [] aa;
        if (isdna) {
            int [] dna	 = string2dna(buffer);
            int [] codon = dna2codon(dna);
            aa	= codon2aa(codon);
        }
        else aa = string2aa(buffer);
        if (isvalidprotein(aa)) {
            for (int i=0; i<aa.length-1; i++) {
                aacounts[aa[i]][aa[i+1]]++;
            }
        }
    }


    static void countaas3(StringBuffer buffer, int [][][] aacounts, boolean isdna) {
        int [] aa;
        if (isdna) {
            int [] dna	 = string2dna(buffer);
            int [] codon = dna2codon(dna);
            aa	= codon2aa(codon);
        }
        else aa = string2aa(buffer);
        if (isvalidprotein(aa)) {
            for (int i=0; i<aa.length-2; i++) {
                aacounts[aa[i]][aa[i+1]][aa[i+2]]++;
            }
        }
    }


    // no premature stop codons? or gaps
    static boolean isvalidprotein(int [] aa) {
        int m = aa.length;
        // if (aa[0]  != 11) return false;
        // if (aa[m-1]!= 21) return false;
        for (int i=1; i<m-1; i++) {
            if ((aa[i] == 0) || (aa[i] == 21)) return false;
        }
        if (aa[m-1]==0) return false;
        return true;
    }


    static int [] deltavec(int len, int i, int ival) {
        int [] dv = new int[len];
        dv[i] = ival;
        return dv;
    }

    static double [] deltavec(int len, int i, double ival) {
        double [] dv = new double[len];
        dv[i] = ival;
        return dv;
    }


    static int [] string2dna(StringBuffer buffer) {
        int n = buffer.length();
        int [] dna = new int[n];
        for (int i = 0; i<n; i++) {
            dna[i] = basetoint(buffer.charAt(i));
        }
        return dna;
    }


    static int [][] string2dna(StringBuffer [] buffer) {
        int m = buffer.length;
        int [][] dna = new int [m][];
        for (int i=0; i<m; i++) dna[i] = string2dna(buffer[i]);
        return dna;
    }

    static int [][] string2aa(StringBuffer [] buffer) {
        int m = buffer.length;
        int [][] aa = new int [m][];
        for (int i=0; i<m; i++) aa[i] = string2aa(buffer[i]);
        return aa;
    }



    static int [] string2aa(StringBuffer buffer) {
        int m = buffer.length();
        int [] aa = new int[m];
        for (int i = 0; i<m; i++) {
            aa[i] = aatoint(buffer.charAt(i));
        }
        return aa;
    }


    // 1 if upper-case; 0 otherwise
    static int [] buffer2bits(StringBuffer sb) {
        int m = sb.length();
        int [] bv = new int[m];
        for (int i = 0; i<m; i++) {
            if (Character.isUpperCase(sb.charAt(i))) bv[i] = 1;
            else bv[i] = 0;
        }
        return bv;
    }

    // longest run of 1's in 0-1 vec
    static int longestrun(int [] bitvec) {
        int maxlen = 0;
        int n = bitvec.length;
        int i = 0;
        while (i<n) {
            if (bitvec[i] > 0) {
                int startdex = i;
                i++;
                while (i < n && bitvec[i] > 0) {
                    i++;
                }
                int stopdex = i-1;
                int len = stopdex-startdex+1;
                if (len >= maxlen) {
                    maxlen=len;
                }
            }
            else i++;
        }
        return maxlen;
    }

    // longest run of 1's in 0-1 vec
    static int longestrun(double [] bitvec) {
        int maxlen = 0;
        int n = bitvec.length;
        int i = 0;
        while (i<n) {
            if (bitvec[i] > 0) {
                int startdex = i;
                i++;
                while (i < n && bitvec[i] > 0) {
                    i++;
                }
                int stopdex = i-1;
                int len = stopdex-startdex+1;
                if (len >= maxlen) {
                    maxlen=len;
                }
            }
            else i++;
        }
        return maxlen;
    }



    static int [] dna2codon(int [] dna) {
        int n = dna.length;
        int m = n/3;
        int [] codon = new int[m];
        for (int i=0; i<m; i++) {
            codon[i] = 25*dna[3*i] + 5*dna[3*i+1] + dna[3*i+2];
        }
        return codon;
    }

    // "sparse" codon --- base 4 instead of 5
    static int [] dna2scodon(int [] dna) {
        int n = dna.length;
        int m = n/3;
        int [] scodon = new int[m];
        for (int i=0; i<m; i++) {
            scodon[i] = 16*(dna[3*i]-1) + 4*(dna[3*i+1]-1) + (dna[3*i+2]-1);
        }
        return scodon;
    }

    static int [][] dna2scodon(int [][] dna) {
        int [][] sc = new int[dna.length][];
        for (int i=0; i<dna.length; i++) sc[i] = dna2scodon(dna[i]);
        return sc;
    }

    static int [] dna2aa(int [] dna) {
        return codon2aa(dna2codon(dna));
    }


    static int [][] dna2aa(int [][] dna) {
        return codon2aa(dna2codon(dna));
    }

    static int [][] codon2aa(int [][] codon) {
        int [][] aa = new int[codon.length][];
        for (int i=0; i<codon.length; i++) aa[i] = codon2aa(codon[i]);
        return aa;
    }

    static int [][] dna2codon(int [][] dna) {
        int [][] codon = new int[dna.length][];
        for (int i=0; i<dna.length; i++) codon[i] = dna2codon(dna[i]);
        return codon;
    }


    static int tnt2aa(int nt1, int nt2, int nt3) {
        return geneticcode[25*nt1 + 5*nt2 + nt3];
    }


    static int tnt2aa(int [] tnt) {
        return geneticcode[25*tnt[0] + 5*tnt[1] + tnt[2]];
    }


    static int [] codon2aa(int [] codon) {
        int m = codon.length;
        int [] aa = new int[m];
        for (int i=0; i<m; i++) {
            aa[i] = geneticcode[codon[i]];
        }
        return aa;
    }


    static StringBuffer aa2string(int [] aa) {
        int m = aa.length;
        StringBuffer buffer = new StringBuffer(m);
        buffer.setLength(m);
        for (int i=0; i<m; i++) {
            buffer.setCharAt(i,aanames[aa[i]]);
        }
        return buffer;
    }


    static StringBuffer [] aa2string(int [][] aa) {
        int n = aa.length;
        StringBuffer [] buffer = new StringBuffer[n];
        for (int j=0; j<n; j++) {
            int m = aa[j].length;
            buffer[j] = new StringBuffer(m);
            buffer[j].setLength(m);
            for (int i=0; i<m; i++) {
                buffer[j].setCharAt(i,aanames[aa[j][i]]);
            }
        }
        return buffer;
    }



    static StringBuffer dna2string(int [] dna) {
        int n = dna.length;
        StringBuffer buffer = new StringBuffer(n);
        buffer.setLength(n);
        for (int i=0; i<n; i++) {
            buffer.setCharAt(i,ntnames[dna[i]]);
        }
        return buffer;
    }


    static StringBuffer int2string(int [] arr, String [] names) {
        int n = arr.length;
        StringBuffer buffer = new StringBuffer(n);
        // buffer.setLength(n);
        for (int i=0; i<n; i++) {
            buffer.append(names[arr[i]]);
        }
        return buffer;
    }

    static StringBuffer int2string(int [] arr, char [] names) {
        int n = arr.length;
        StringBuffer buffer = new StringBuffer(n);
        buffer.setLength(n);
        for (int i=0; i<n; i++) {
            buffer.setCharAt(i,names[arr[i]]);
        }
        return buffer;
    }



    static double logbinomial(int n, int k, double p){
        // log C(n,k) p^k (1-p)^(n-k)
        double lb = k*Math.log(p) + (n-k)*Math.log(1-p);
        double nn = 1.0*n;
        for (int i=0; i<k; i++) {lb = lb + Math.log((nn-i)/(k-i));}
        return lb;
    }


    static double [][] logbinomiallut(int maxn, double p) {
        double [][] lut = new double[maxn][maxn];
        for (int i=0; i<maxn; i++) {
            for (int j=0; j<=i; j++) {
                lut[i][j] = logbinomial(i,j,p);
            }
        }
        return lut;
    }



    // quicker version
    static double [][] logbinomiallut2(int maxn, double p) {
        double q = 1 - p;
        double [][] lut = new double[maxn+1][maxn+1];

        for (int k=0; k<=maxn; k++) {
            lut[k][k] = k*Math.log(p);
            for (int n=k+1; n<=maxn; n++) {
                lut[n][k] = lut[n-1][k] + Math.log((n+1.0)/(n+1.0-k)*q);
            }
        }
        return lut;
    }


    // brute force
    static double [] lps(int [] seq, double p) {
        // printmatrix(seq);

        double [] sss = {0, 0, 0, 0.0}; // start index, end index, count, logprob

        int besti = 0;
        int bestj = 0;
        double bestscore = 0.0;
        int bestcount = 0;

        double score;

        int n = seq.length;

        for (int i=0; i<n; i++) {
            int ct = 0;
            int len = 0;
            for (int j=i; j<n; j++) {
                len++;
                ct = ct + seq[j];
                score = logbinomial(len,ct,p);
                // System.out.println(i + " " + j + " " + len + " " + ct + " " + (float) score);
                if (score < bestscore) {bestscore = score; besti=i; bestj=j; bestcount=ct;}
            }
        }
        sss[0] = 1.0*besti;
        sss[1] = 1.0*bestj;
        sss[2] = 1.0*bestcount;
        sss[3] = bestscore;
        return sss;
    }



    // brute force, but with look-up-table.
    static double [] lps(int [] seq, double p, double [][] lut, int maxlen) {

        double [] sss = {0, 0, 0, 0.0}; // start index, end index, count, logprob

        // int maxl = lut.length;

        int besti = 0;
        int bestj = 0;
        double bestscore = 0.0;
        int bestcount = 0;

        double score;

        int n = seq.length;


        for (int i=0; i<n; i++) {
            int ct = 0;
            int len = 0;
            for (int j=i; j<Math.min(i+maxlen,n); j++) {
                len++;
                ct = ct + seq[j];
                score = lut[len][ct];
                if (score < bestscore) {bestscore=score; besti=i; bestj=j; bestcount=ct;}

            }
        }
        sss[0] = 1.0*besti;
        sss[1] = 1.0*bestj;
        sss[2] = 1.0*bestcount;
        sss[3] = bestscore;
        return sss;
    }


    // brute force, but with look-up-table. tail=0 both ways; tail=1 overrep; tail=-1 underrep
    static double [] lps(int [] seq, double p, double [][] lut, int minlen, int maxlen, int tail) {

        double [] sss = {0, 0, 0, 0.0}; // start index, end index, count, logprob

        // int maxl = lut.length;

        int besti = 0;
        int bestj = 0;
        double bestscore = 0.0;
        int bestcount = 0;

        double score;

        int n = seq.length;

        boolean legal = true;

        for (int i=0; i<n; i++) {
            int ct = 0;
            int len = 0;
            for (int j=i; j<Math.min(i+maxlen,n); j++) {
                len++;
                ct = ct + seq[j];
                legal = true;
                if (len<minlen) legal = false;
                else if (tail==1 && (ct+0.0)/(len+0.0) < p) legal = false;
                else if (tail==-1 && (ct+0.0)/(len+0.0) > p) legal = false;

                if (legal) {
                    score = lut[len][ct];
                    if (score < bestscore) {bestscore=score; besti=i; bestj=j; bestcount=ct;}
                }

            }
        }
        sss[0] = 1.0*besti;
        sss[1] = 1.0*bestj;
        sss[2] = 1.0*bestcount;
        sss[3] = bestscore;
        return sss;
    }



    //sparse and one-tailed, but with look-up-table
    static double [] slps(int [] seq, double p, double [][] lut, int maxlen) {

        double [] sss = {0, 0, 0, 0.0}; // start index, end index, count, logprob

        // int maxl = lut.length;

        int besti = 0;
        int bestj = 0;
        double bestscore = 0.0;
        int bestcount = 0;

        double score;

        int n = seq.length;


        int [] nzdex = sparsesupport(seq,1);

        int m = nzdex.length;

        int len;
        int ct;

        for (int i=0; i<m; i++) {
            for (int j=i; j<m; j++) {
                len = nzdex[j]-nzdex[i]+1;
                if (len >= maxlen) continue;
                ct = j-i + 1;
                score = lut[len][ct];
                if (score < bestscore) {bestscore=score; besti=nzdex[i]; bestj=nzdex[j]; bestcount=ct;}
            }
        }
        sss[0] = 1.0*besti;
        sss[1] = 1.0*bestj;
        sss[2] = 1.0*bestcount;
        sss[3] = bestscore;
        return sss;
    }


    ////////////////////

    static HashMap<String, String> readhashtable(String filename, String defaultvalue) {
        HashMap<String,String> ht = new HashMap<String,String>();

        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            String line;
            int i = 1;
            while ((line = in.readLine()) != null) {
                line.trim();
                if (line != null) {
                    String [] chunks = line.split("\\t");
                    if (defaultvalue == "incr") {ht.put(chunks[0], new String(""+i));}
                    else if (chunks.length > 1) ht.put(chunks[0], chunks[1]);
                    else if (defaultvalue == null) {ht.put(chunks[0], chunks[0]);}
                    else ht.put(chunks[0], defaultvalue);
                    i++;
                    // System.out.println("#"+chunks[0] + "aaa" + ht.get(chunks[0]) + "bbb");
                }

            }
        }
        catch (IOException e) {
            System.out.println("# Couldn't open " + filename);
        }
        // System.out.println(ht.size());

        return ht;
    }


    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////


    static void fastaprint(String name, StringBuffer seq, int w) {
        System.out.println(">" + name);
        int n = seq.length();
        for (int i=0; i<n/w; i++) {
            System.out.println(seq.substring(i*w,i*w+w));
        }
        if (n % w >0) System.out.println(seq.substring(n/w*w));
    }



    static int [] composition(int [] vec, int max) {
        int [] comp = new int[max+1];
        for (int i=0; i<vec.length; i++) {comp[vec[i]]++;}
        return comp;
    }

    static int [] composition(int [] vec) {
        return composition(vec, maxint(vec));
    }


    static int [] aacomposition(int [] vec) {
        return composition(vec,22);
    }

    static int [] dnacomposition(int [] vec) {
        return composition(vec,5);
    }


    static double [] normalize(int [] arr) {
        int n = arr.length;
        double [] narr = new double [n];
        double sm = 1.0*sum(arr);
        for (int i=0; i<n; i++) narr[i] = arr[i]/sm;
        return narr;
    }

    static double [] normalize(double [] arr) {
        int n = arr.length;
        double [] narr = new double [n];
        double sm = 1.0*sum(arr);
        if (sm < 0.000000000001) sm = 1;
        for (int i=0; i<n; i++) narr[i] = arr[i]/sm;
        return narr;
    }

    // makes each row sum 1
    static double [][] normalizerows(double [][] mat) {
        int n = mat.length;
        double [][] nmat = new double[n][];
        for (int i=0; i<n; i++) nmat[i] = normalize(mat[i]);
        return nmat;
    }

    // makes total sum	1
    static double [][] normalizemat(double [][] mat) {
        int m = mat.length;
        int n = mat[0].length;
        double [][] nmat = new double[m][n];
        double ts = 0.0;
        for (int i=0; i<m; i++) {
            for (int j=0; j<n; j++) {
                ts = ts+mat[i][j];
            }
        }
        if (ts < 0.00000000001) ts = 0.00000000001;


        for (int i=0; i<m; i++) {
            for (int j=0; j<n; j++) {
                nmat[i][j] = mat[i][j]/ts;
            }
        }
        return nmat;
    }


    ///////////////////////////////////////////////////////////////////////////
    /////////////
    /////////////	  matrix functions

    // ax + by
    static double [] axpby(double a, double [] x, double b, double [] y) {
        int n = x.length;
        int m = y.length;
        if (m<n) n=m;
        double [] arr = new double[n];
        for (int i=0; i<n; i++) arr[i] = a*x[i] + b*y[i];
        return arr;
    }


    // ax + by
    static double [][] axpby(double a, double [][] x, double b, double [][] y) {
        int m1 = x.length;
        int m2 = y.length;
        int n1 = x[0].length;
        int n2 = y[0].length;
        if (m1>m2) m1=m2;
        if (n1>n2) n1=n2;
        double [][] mat = new double[m1][n1];
        for (int i=0; i<m1; i++) {
            for (int j=0; j<n1; j++) mat[i][j] = a*x[i][j] + b*y[i][j];
        }
        return mat;
    }


    static double [] flatten(double [][] mat) {
        int m = mat.length;
        int n = mat[0].length;
        double [] arr = new double [m*n];
        for (int i=0; i<m; i++) {
            for (int j=0; j<n; j++) {
                arr[i*m + j] = mat[i][j];
            }
        }
        return arr;

    }

    static double [] rowsums(double [][] mat) {
        int m = mat.length;
        double [] rs = new double[m];
        for (int i=0; i<m; i++) {
            double s = 0;
            for (int j=0; j<mat[i].length; j++) s = s + mat[i][j];
            rs[i] = s;
        }
        return rs;
    }

    static double [] colsums(double [][] mat) {
        int m = mat.length;
        int n = mat[0].length;
        double [] cs = new double[n];
        for (int j=0; j<n; j++) {
            double s = 0;
            for (int i=0; i<m; i++) s = s + mat[i][j];
            cs[j] = s;
        }
        return cs;
    }

    // ax + by + c
    static double [] axpbypc(double a, double [] x, double b, double [] y, double c) {
        int n = x.length;
        int m = y.length;
        if (m<n) n=m;
        double [] arr = new double[n];
        for (int i=0; i<n; i++) arr[i] = a*x[i] + b*y[i] + c;
        return arr;
    }

    // ax + b
    static double [] axpb(double a, double [] x, double b) {
        int n = x.length;
        double [] arr = new double[n];
        for (int i=0; i<n; i++) arr[i] = a*x[i] + b;
        return arr;
    }

    static double [][] axpb(double a, double [][] x, double b) {
        int m = x.length;
        int n = x[0].length;
        double [][] arr = new double[m][n];
        for (int i=0; i<m; i++) {
            for (int j=0; j<n; j++) {
                arr[i][j] = a*x[i][j] + b;
            }
        }
        return arr;
    }


    // |x|
    static double [] absval(double [] x) {
        int n = x.length;
        double [] arr = new double[n];
        for (int i=0; i<n; i++) arr[i] = Math.abs(x[i]);
        return arr;
    }


    static double [][] matmult(double [][] mat1, double [][] mat2){
        int a = mat1.length;
        int b = mat1[0].length;
        int c = mat2.length;
        int d = mat2[0].length;
        if (b!=c)  {
            System.out.println("inner dimensions don't match: "	  + a + " x " + b + ", " +  c + " x " + d);
            return (new double[1][1]);
        }
        double [][] mat3 = new double[a][d];
        double t;
        for (int i=0; i<a; i++) {
            for (int j=0; j<d; j++) {
                t=0;
                for (int k=0; k<b; k++) t = t + mat1[i][k]*mat2[k][j];
                mat3[i][j] = t;
            }
        }
        return mat3;
    }


    static double [] matmult(double [][] mat1, double [] vec){
        int a = mat1.length;
        int b = mat1[0].length;
        int c = vec.length;
        if (b!=c)  {
            System.out.println("inner dimensions don't match: "	  + a + " x " + b + ", " +  c);
            return (new double[1]);
        }
        double [] vec2 = new double[a];
        double t;
        for (int i=0; i<a; i++) {
            t=0;
            for (int k=0; k<b; k++) t = t + mat1[i][k]*vec[k];
            vec2[i] = t;
        }
        return vec2;
    }


    static double [][] outerproduct(double [] vec1, double [] vec2) {
        int m = vec1.length;
        int n = vec2.length;
        double [][] op = new double [m][n];
        for (int i=0; i<m; i++) {
            for (int j=0; j<n; j++) {
                op[i][j] = vec1[i]*vec2[j];
            }
        }
        return op;
    }


    static double innerproduct(double [] vec1, double [] vec2) {
        int m = vec1.length;
        double ip = 0.0;
        for (int i=0; i<m; i++)	 ip = ip + vec1[i]*vec2[i];
        return ip;
    }


    static double [][] diagonal(double [] vec) {
        int n = vec.length;
        double [][] mat = new double [n][n];
        for (int i=0; i<n; i++) mat[i][i] = vec[i];
        return mat;
    }


    static double [] diagonal(double [][] mat) {
        int n = Math.min(mat.length, mat[0].length);
        double [] vec = new double[n];
        for (int i=0; i<n; i++) vec[i] = mat[i][i];
        return vec;
    }


    // pointwise division: vec1./vec2
    static double [] divp(double [] vec1, double [] vec2) {
        int n = vec1.length;
        double [] nvec = new double[n];
        for (int i=0; i<n; i++) nvec[i]=vec1[i]/vec2[i]; // change if both 0?
        return nvec;
    }


    // pointwise division: vec1./vec2
    static double [][] divp(double [][] vec1, double [][] vec2) {
        int m = vec1.length;
        int n = vec1[0].length;
        double [][] nvec = new double[m][n];
        for (int j=0; j<m; j++) {
            for (int i=0; i<n; i++) {
                nvec[j][i]=vec1[j][i]/vec2[j][i]; // change if both 0?
            }
        }
        return nvec;
    }


    static double [][] copymatrix(double [][] mat) {
        int m = mat.length;
        double [][] nmat = new double[m][];
        for (int i=0; i<m; i++) {
            int n = mat[i].length;
            nmat[i] = new double[n];
            for (int j=0; j<n; j++) nmat[i][j]=mat[i][j];
        }
        return nmat;
    }

    static double [] copymatrix(double [] mat) {
        int m = mat.length;
        double [] nmat = new double[m];
        for (int i=0; i<m; i++) {
            nmat[i]=mat[i];
        }
        return nmat;
    }



    static double factorial(int n) {
        double fact = 1.0;
        for (int i=1; i<=n; i++) fact = fact*i;
        return fact;
    }


    static double [][] matsum(double [][] mat1, double [][] mat2){
        int a = mat1.length;
        int b = mat1[0].length;
        int c = mat2.length;
        int d = mat2[0].length;
        if (a!=c || b !=d)  {
            System.out.println("matrix dimensions don't match: "   + a + " x " + b + ", " +  c + " x " + d);
            return (new double[1][1]);
        }
        double [][] mat3 = new double[a][b];
        double t;
        for (int i=0; i<a; i++) {
            for (int j=0; j<b; j++) {
                mat3[i][j] = mat1[i][j] + mat2[i][j];
            }
        }
        return mat3;
    }


    static double [][] scalarmult(double [][] mat, double t) {
        int a = mat.length;
        int b = mat[0].length;
        double [][] mat2 = new double[a][b];
        for (int i=0; i<a; i++) {
            for (int j=0; j<b; j++) {
                mat2[i][j] = t*mat[i][j];
            }
        }
        return mat2;
    }



    static int [][] sparsesupport(double [][] mat, double thresh) {
        int m = mat.length;
        int n = mat[0].length;
        int [][] ss = new int[m][];
        for (int i=0; i<m; i++) {
            int ct = 0;
            for (int j=0; j<n; j++) if (mat[i][j]>thresh) ct++;
            ss[i] = new int[ct];
            ct = 0;
            for (int j=0; j<n; j++) if(mat[i][j]>thresh) {ss[i][ct]=j;ct++;}
        }
        return ss;
    }


    static int [] sparsesupport(int [] vec, int target) {
        int m = vec.length;
        int ct = 0;
        for (int i=0; i<m; i++) {
            if (vec[i]==target) ct++;
        }
        int [] ss = new int[ct];
        ct = 0;
        for (int i=0; i<m; i++) {
            if (vec[i]==target) {ss[ct]=i; ct++;}
        }
        return ss;
    }




    static double [][] logmat(double [][] mat){
        int n = mat.length;
        int m = mat[0].length;
        double [][] lmat = new double[n][m];
        for (int i=0; i<n; i++) {
            for (int j=0; j<m; j++) {
                lmat[i][j] = Math.log(mat[i][j]); // check for 0 or neg?
            }
        }
        return lmat;
    }

    static double [] logmat(double [] mat){
        int n = mat.length;
        double [] lmat = new double[n];
        for (int i=0; i<n; i++) {
            lmat[i] = Math.log(mat[i]);
        }
        return lmat;
    }


    static double [][] eyemat(int r, int c) {
        double [][] mat = new double[r][c];
        for (int i=0; i<Math.min(r,c); i++) {
            mat[i][i] = 1;
        }
        return mat;
    }

    static double [][] zeromat(int r, int c) {
        double [][] mat = new double[r][c];
        return mat;
    }

    static double [][] constantmat(int r, int c, double val) {
        double [][] mat = new double[r][c];
        for (int i=0; i<r; i++) for (int j=0; j<c; j++) mat[i][j] = val;
        return mat;
    }

    static double [] constantmat(int r, double val) {
        double [] mat = new double[r];
        for (int i=0; i<r; i++)	 mat[i] = val;
        return mat;
    }


    // pointwise inverse
    static double [] inverse (double [] vec) {
        int n =	 vec.length;
        double [] iv = new double[n];
        for (int i=0; i<n; i++) iv[i] = 1.0/vec[i];
        return iv;
    }


    static double max(double [] arr){
        int n = arr.length;
        double mx = arr[0];
        int dex = 0;
        for (int i=1; i<n; i++) {
            if (arr[i] > mx) {
                mx = arr[i];
                dex = i;
            }
        }
        return mx;
    }


    static double min(double [] arr){
        int n = arr.length;
        double mx = arr[0];
        int dex = 0;
        for (int i=1; i<n; i++) {
            if (arr[i] < mx) {
                mx = arr[i];
                dex = i;
            }
        }
        return mx;
    }

    // returns val and index
    static double [] max2(double [] arr){
        int n = arr.length;
        double mx = arr[0];
        int dex = 0;
        for (int i=1; i<n; i++) {
            if (arr[i] > mx) {
                mx = arr[i];
                dex = i;
            }
        }
        double [] md = new double[2];
        md[0] = mx; md[1] = dex;
        return md;
    }


    static int maxint(int [] arr){
        int n = arr.length;
        int mx = arr[0];
        int dex = 0;
        for (int i=1; i<n; i++) {
            if (arr[i] > mx) {
                mx = arr[i];
                dex = i;
            }
        }
        return mx;
    }


    // needs to be square here
    static double [][] readmatrix(String filename) {
        BufferedReader in;
        String line;

        double [][] mat = new double[1][1];

        StringTokenizer toker;

        try {
            in = new BufferedReader(new FileReader(filename));
            line = in.readLine();
            toker = new StringTokenizer(line);
            int n = toker.countTokens();
            mat = new double[n][n];
            for (int i=0; i<n; i++) {
                if (i>0) { // first line is already read
                    line = in.readLine();
                    toker = new StringTokenizer(line);
                }
                for (int j=0; j<n; j++) {
                    mat[i][j] = Double.parseDouble(toker.nextToken());
                }
            }
            in.close();

        }
        catch (IOException e) {
            System.out.println("Couldn't open " + filename);
        }

        return mat;

    }


    static void printmatrix(int [][] mat) {
        for (int i=0; i<mat.length; i++) {
            for (int j=0; j<mat[i].length; j++) {
                System.out.print(mat[i][j] + "\t");
            }
            System.out.println();
        }
    }


    static void printmatrix(double [][] mat) {
        for (int i=0; i<mat.length; i++) {
            for (int j=0; j<mat[i].length; j++) {
                System.out.print((float) (mat[i][j]) + "\t");
            }
            System.out.println();
        }
    }


    static void printmatrix(int [] mat) {
        for (int i=0; i<mat.length; i++) {
            System.out.println(mat[i]);
        }
        System.out.println();
    }



    static double [] mapfreq(double [] freq, int [] map) {
        int n = maxint(map) + 1;
        double [] newfreq = new double[n];
        for (int i=0; i<freq.length; i++) newfreq[map[i]] = newfreq[map[i]]+freq[i];
        return newfreq;
    }


    static int [] mapseq(int [] arr, int [] map) {
        int n = arr.length;
        int [] newarr = new int[n];
        for (int i=0; i<n; i++) newarr[i] = map[arr[i]];
        return newarr;
    }


    static double [] mapseq(int [] arr, double [] map) {
        int n = arr.length;
        double [] newarr = new double[n];
        for (int i=0; i<n; i++) newarr[i] = map[arr[i]];
        return newarr;
    }


    static double [] slidingaverage(double [] arr, int ww, boolean reflect) {
        int n = arr.length;
        // System.out.println(n);
        if (n==0) return arr;
        int w = ww/2;
        if (w>=n) w=n-1;
        double [] padarr = new double [n+2*w];
        System.arraycopy(arr,0,padarr,w,n);
        // printmatrix(padarr);
        if (reflect) {
            for (int i=0; i<=w; i++) {
                padarr[i] = arr[w-i];
                padarr[n+2*w-1-i] = arr[n-w-1+i];
            }
        }
        // printmatrix(padarr);
        double [] sa = new double[n];
        double score;
        for (int i=0; i<n; i++) {
            score = 0;
            for (int j=-w;j<=w; j++) {
                score = score+padarr[i+w+j];
            }
            sa[i] = score/(2.0*w+1);
        }

        return sa;
    }


    // if  reflect==F, average over smaller windows
    static double [] slidingaverage2(double [] arr, int ww, boolean reflect) {
        int n = arr.length;
        // System.out.println(n);
        if (n==0) return arr;
        int w = ww/2;
        if (w>=n) w=n-1;
        double [] padarr = new double [n+2*w];
        System.arraycopy(arr,0,padarr,w,n);
        // printmatrix(padarr);
        if (reflect) {
            for (int i=0; i<=w; i++) {
                padarr[i] = arr[w-i];
                padarr[n+2*w-1-i] = arr[n-w-1+i];
            }
        }
        // printmatrix(padarr);
        double [] sa = new double[n];
        double score;
        int inrange;
        for (int i=0; i<n; i++) {
            score = 0;
            inrange = 0;
            for (int j=-w;j<=w; j++) {
                score = score+padarr[i+w+j];
                if ((i+j >=0) &(i+j<n)) inrange++;
            }
            if (reflect) {sa[i] = score/(2.0*w+1);} else {sa[i] = score/(1.0*inrange);}
        }

        return sa;
    }


    static void printaausage(double [] aafreq) {
        for (int i=0; i<22; i++) System.out.println("# " + aanamesb[i] + " " + (float) aafreq[i]);
    }


    static void printaausage(double [][] aafreq) {
        for (int i=0; i<22; i++) {
            System.out.print("# " + aanamesb[i]);
            for (int j=0; j<aafreq.length; j++) System.out.print("\t" + (float) aafreq[j][i]);
            System.out.println();
        }
    }


    static void printaausage(double [][] aafreq, int sd) {
        for (int i=0; i<22; i++) {
            System.out.print("# " + aanamesb[i]);
            for (int j=0; j<aafreq.length; j++) System.out.print("\t" + (float) (Math.round(sd*aafreq[j][i])/(1.0*sd)));
            System.out.println();
        }
    }


}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


class hmm {


    int [] data;  // entries in {0,1,2,..., no-1}
    int ns; // number of states
    int no; // number of output symbols
    double [][] tprob; // (ns x ns) transition probs
    double [][] eprob; // (ns x no) emmision probs
    double []	iprob; // (ns)	    initial probs
    double [][] ctprob;
    double [][] ceprob;
    double [] ciprob;
    double [][] ltprob; // (ns x ns) log transition probs
    double [][] leprob; // (ns x no) log emmision probs
    double []	liprob; // (ns)	     log initial probs


    double [] fprob;	// (ns) prob of transitioning to end state
    double [] lfprob;	// (ns)



    char [] states; // symbols for printing
    String [] names;

    char [] ustates; // unique
    String [] unames;


    double lviterbiprob = -1.0/0;
    double lmarginalprob = -1.0/0;



    int [] subtrellis; // (ns)	bitmask


    int [] classes; // for merging states
    int numclasses;



    double etst = 0.0; // expected time on subtrellis
    double lpst = 0.0; // log prob of never visiting subtrellis

    double lmp = 0.0; // log of marginal prob of output

    double [][] postprob;
    double [] margcollapse;

    int [] viterbipath;
    int [] mappath;



    boolean sparse = false;



    // read from file
    hmm(String filename) {


        tprob = new double[2][2];
        eprob = new double[2][2];
        iprob = new double[2];


        BufferedReader in;
        String line;

        StringTokenizer toker;

        try {
            in = new BufferedReader(new FileReader(filename));
            line = in.readLine();
            toker = new StringTokenizer(line);
            ns = Integer.parseInt(toker.nextToken());
            line = in.readLine();
            toker = new StringTokenizer(line);
            no = Integer.parseInt(toker.nextToken());
            tprob = new double[ns][ns];
            eprob = new double[ns][no];
            iprob = new double[ns];
            line = in.readLine(); // blank
            for (int i=0; i<ns; i++) {
                line = in.readLine();
                toker = new StringTokenizer(line);
                for (int j=0; j<ns; j++) {
                    tprob[i][j] = Double.parseDouble(toker.nextToken());
                }
            }
            line = in.readLine(); // blank
            line = in.readLine();
            toker = new StringTokenizer(line);
            for (int i=0; i<ns; i++) {
                iprob[i] = Double.parseDouble(toker.nextToken());
            }
            line = in.readLine(); // blank
            for (int j=0; j<no; j++) {
                line = in.readLine();
                toker = new StringTokenizer(line);
                for (int i=0; i<ns; i++) {
                    eprob[i][j] = Double.parseDouble(toker.nextToken());
                }
            }
            in.close();

        }
        catch (IOException e) {
            System.out.println("Couldn't open " + filename);
        }

        initialize();

    }


    hmm(double [][] tprob, double [][] eprob, double [] iprob) {
        this.tprob = tprob;
        this.eprob = eprob;
        this.iprob = iprob;
        initialize();
    }


    public void initialize() {
        ctprob = spewey.cumsum(tprob,true);
        ceprob = spewey.cumsum(eprob,true);
        ciprob = spewey.cumsum(iprob);



        ltprob = spewey.logmat(tprob);
        leprob = spewey.logmat(eprob);
        liprob = spewey.logmat(iprob);
        ns = iprob.length;
        no = eprob[0].length;


        double [] rs = spewey.rowsums(tprob);
        fprob = new double[ns];

        boolean freeend = true;

        for (int i=0; i<ns; i++) {
            fprob[i] = Math.max(0.0,1.0-rs[i]);
            if (fprob[i] > 0.0001) freeend = false;
        }
        if (freeend) {
            for (int i=0; i<ns; i++) fprob[i] = 1.0;
        }
        lfprob = spewey.logmat(fprob);


        subtrellis = new int[ns];
        for (int i=0; i<ns; i++) subtrellis[i]=1;
        states = new char[ns];
        for (int i=0; i<ns; i++) states[i]=(char) i;
        names = new String[ns];
        for (int i=0; i<ns; i++) names[i]= new String("s"+i);
        classes = new int[ns];
        for (int i=0; i<ns; i++) classes[i]=i;
        numclasses = ns;
        ustates = states;
        unames = names;
    }


    public void setclasses(int [] classes) {
        this.classes = classes;
        int nc = 0;
        numclasses = 0;
        for (int i=0; i<classes.length; i++) {
            if (classes[i]>nc) nc=classes[i];
        }
        numclasses = nc+1;
        // System.out.println("nc " + numclasses);
        // spewey.printmatrix(classes);
        // System.out.println(numclasses)
        unames = new String[numclasses];
        ustates = new char[numclasses];
        for (int i=classes.length-1; i>=0; i--) {
            unames[classes[i]] = names[i];
            ustates[classes[i]] = states[i];
        }

    }

    public void setnames(String [] names) {
        this.names = names;
        this.unames = names;

    }


    public void print() {
        System.out.println(ns + " \t ## num states");
        System.out.println(no + " \t ## num symbols");
        System.out.println();
        for (int i=0; i<ns; i++) {
            for (int j=0; j<ns; j++) {
                System.out.print((float) tprob[i][j] + "\t");
            }
            System.out.println();
        }
        System.out.println();

        for (int j=0; j<ns; j++) {
            System.out.print((float) iprob[j] + "\t");
        }
        System.out.println();
        System.out.println();

        for (int j=0; j<no; j++) {
            for (int i=0; i<ns; i++) {
                System.out.print((float) eprob[i][j] + "\t");
            }
            System.out.println();
        }
        System.out.println();

    }


    public void write(String filename) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(filename));

            out.write(ns + " \t ## num states");  out.newLine();
            out.write(no + " \t ## num symbols"); out.newLine();
            out.newLine();
            for (int i=0; i<ns; i++) {
                for (int j=0; j<ns; j++) {
                    out.write((float) tprob[i][j] + "\t");
                }
                out.newLine();
            }
            out.newLine();

            for (int j=0; j<ns; j++) {
                out.write((float) iprob[j] + "\t");
            }
            out.newLine();
            out.newLine();

            for (int j=0; j<no; j++) {
                for (int i=0; i<ns; i++) {
                    out.write((float) eprob[i][j] + "\t");
                }
                out.newLine();
            }
            out.newLine();
            out.close();
        }
        catch (IOException e) {
            // do something
        }

    }


    public int [] viterbidecode(int [] seq) {
        int n = seq.length;
        int [] vit = new int[n];

        double [][] s = new double[ns][n];
        int [][] tb = new int[ns][n]; // traceback

        for (int i=0; i<ns; i++) s[i][0] = iprob[i]*eprob[i][seq[0]];
        for (int t=1; t<n; t++) {
            for (int i=0; i<ns; i++) {
                int bestdex = 0;
                double bestscore = tprob[0][i]*s[0][t-1];
                for (int k=1; k<ns; k++) {
                    if (tprob[k][i]*s[k][t-1]  > bestscore) {
                        bestscore =  tprob[k][i]*s[k][t-1];
                        bestdex = k;
                    }
                }
                s[i][t] = bestscore*eprob[i][seq[t]];
                tb[i][t] = bestdex;
            }
        }
        int bestdex = 0;
        double bestscore = s[0][n-1];
        for (int k=1; k<ns; k++) {
            if (s[k][n-1] > bestscore) {bestscore = s[k][n-1]; bestdex = k;}
        }
        vit[n-1] = bestdex;
        for (int t = n-2; t>=0; t--) {vit[t] = tb[vit[t+1]][t+1];}

        spewey.printmatrix(spewey.transpose(s)); System.out.println();

        viterbipath = vit;
        return vit;
    }


    public int [] viterbidecodel(int [] seq) {
        int n = seq.length;
        int [] vit = new int[n];

        double [][] s = new double[ns][n];
        int [][] tb = new int[ns][n];	    // traceback

        for (int i=0; i<ns; i++) s[i][0] = liprob[i]+leprob[i][seq[0]]; // initial probs
        for (int t=1; t<n; t++) {
            for (int i=0; i<ns; i++) {
                int bestdex = 0;
                double bestscore = ltprob[0][i] + s[0][t-1];
                for (int k=1; k<ns; k++) {
                    if (ltprob[k][i] + s[k][t-1]  > bestscore) {
                        bestscore =  ltprob[k][i] + s[k][t-1];
                        bestdex = k;
                    }
                }
                s[i][t] = bestscore + leprob[i][seq[t]];
                tb[i][t] = bestdex;
            }
        }
        int bestdex = 0;
        double bestscore = s[0][n-1] + lfprob[0];
        for (int k=1; k<ns; k++) {
            if (s[k][n-1] + lfprob[k] > bestscore) {bestscore = s[k][n-1] + lfprob[k]; bestdex = k;}
        }
        vit[n-1] = bestdex;

        for (int t = n-2; t>=0; t--) {vit[t] = tb[vit[t+1]][t+1];}

        // spewey.printmatrix(spewey.transpose(s)); System.out.println();

        lviterbiprob = bestscore;
        if (numclasses < ns) vit = spewey.mapseq(vit,classes);
        viterbipath = vit;
        return vit;
    }


    // sparse transmat
    public int [] sviterbidecodel(int [] seq) {
        int n = seq.length;
        int [] vit = new int[n];

        double [][] s = new double[ns][n];
        int [][] tb = new int[ns][n];	    // traceback


        int [][] sps = spewey.sparsesupport(spewey.transpose(tprob),0.000000001);


        for (int i=0; i<ns; i++) s[i][0] = liprob[i]+leprob[i][seq[0]]; // initial probs
        for (int t=1; t<n; t++) {
            for (int i=0; i<ns; i++) {
                int bestdex = 0;
                double bestscore = ltprob[0][i] + s[0][t-1];
                for (int sk=0; sk < sps[i].length; sk++) {
                    int k =sps[i][sk];
                    // for (int k=1; k<ns; k++) {
                    if (ltprob[k][i] + s[k][t-1]  > bestscore) {
                        bestscore =  ltprob[k][i] + s[k][t-1];
                        bestdex = k;
                    }
                }
                s[i][t] = bestscore + leprob[i][seq[t]];
                tb[i][t] = bestdex;
            }
        }
        int bestdex = 0;
        double bestscore = s[0][n-1] + lfprob[0];
        for (int k=1; k<ns; k++) {
            if (s[k][n-1] + lfprob[k] > bestscore) {bestscore = s[k][n-1] + lfprob[k]; bestdex = k;}
        }
        vit[n-1] = bestdex;

        for (int t = n-2; t>=0; t--) {vit[t] = tb[vit[t+1]][t+1];}

        // spewey.printmatrix(spewey.transpose(s)); System.out.println();

        lviterbiprob = bestscore;

        if (numclasses < ns) vit = spewey.mapseq(vit,classes);
        viterbipath = vit;
        return vit;
    }



    // some labels can be known
    public int [] viterbidecodel(int [] seq, int [] lab) {
        int n = seq.length;
        int [] vit = new int[n];

        double [][] s = new double[ns][n];
        int [][] tb = new int[ns][n]; // traceback

        double badscore = -1.0/0;


        if (lab[0]<0) {for (int i=0; i<ns; i++) s[i][0] = liprob[i]+leprob[i][seq[0]];}
        else {
            for (int i=0; i<ns; i++) s[i][0] = badscore;
            s[lab[0]][0] = liprob[lab[0]]+leprob[lab[0]][seq[0]];
        }

        for (int t=1; t<n; t++) {
            if (lab[t]<0) {
                for (int i=0; i<ns; i++) {
                    int bestdex = 0;
                    double bestscore = ltprob[0][i] + s[0][t-1];
                    for (int k=1; k<ns; k++) {
                        if (ltprob[k][i] + s[k][t-1]  > bestscore) {
                            bestscore =	 ltprob[k][i] + s[k][t-1];
                            bestdex = k;
                        }
                    }
                    s[i][t] = bestscore + leprob[i][seq[t]];
                    tb[i][t] = bestdex;
                }
            }
            else {
                for (int i=0; i<ns; i++) s[i][t] = badscore;
                int bestdex = 0;
                double bestscore = ltprob[0][lab[t]] + s[0][t-1];
                for (int k=1; k<ns; k++) {
                    if (ltprob[k][lab[t]] + s[k][t-1]  > bestscore) {
                        bestscore =  ltprob[k][lab[t]] + s[k][t-1];
                        bestdex = k;
                    }
                }
                s[lab[t]][t] = bestscore + leprob[lab[t]][seq[t]];
                tb[lab[t]][t] = bestdex;
            }
        }
        int bestdex = 0;
        double bestscore = s[0][n-1];
        for (int k=1; k<ns; k++) {
            if (s[k][n-1] > bestscore) {bestscore = s[k][n-1]; bestdex = k;}
        }
        vit[n-1] = bestdex;
        for (int t = n-2; t>=0; t--) {vit[t] = tb[vit[t+1]][t+1];}

        // spewey.printmatrix(spewey.transpose(s)); System.out.println();


        viterbipath = vit;
        return vit;
    }


    public void decodealls(int [] seq) {
        if (sparse)  {
            sviterbidecodel(seq);
            smapdecodel(seq);
        }
        else {
            viterbidecodel(seq);
            mapdecodel(seq);
        }
    }


    public void decodeall(int [] seq) {

        viterbidecodel(seq);
        mapdecodel(seq);
        margcollapse = new double[seq.length];
        etst = 0.0; // expected time in subtrellis

        if (ns == numclasses) {
            for (int j=0; j<seq.length; j++) {
                for (int i=0; i<ns; i++) {
                    if (subtrellis[i]==0) margcollapse[j] = margcollapse[j] + postprob[i][j];
                }
            }
            for (int i=0; i<ns; i++) {
                if (subtrellis[i]==0) etst = etst + spewey.sum(postprob[i]);
            }
        }
        else {
            margcollapse = postprob[0];
            etst = spewey.sum(postprob[0]);
        }

        lpst = logprobsubtrellis(seq, subtrellis);


    }


    public double [][] posterior(int [] seq) {
        int n = seq.length;
        double [][] pp = new double[ns][n];


        // forward algorithm
        double [][] a = new double[ns][n];

        for (int i=0; i<ns; i++) a[i][0] = iprob[i]*eprob[i][seq[0]];
        for (int t=1; t<n; t++) {
            for (int i=0; i<ns; i++) {
                double score = 0;
                for (int k=0; k<ns; k++) {
                    score = score + tprob[k][i]*a[k][t-1];
                }
                a[i][t] = score*eprob[i][seq[t]];
            }
        }

        // backward algorithm
        double [][] b = new double[ns][n];

        for (int i=0; i<ns; i++) b[i][n-1] = 1; // iprob[i];
        for (int t=n-2; t>=0; t--) {
            for (int i=0; i<ns; i++) {
                double score = 0;
                for (int k=0; k<ns; k++) {
                    score = score + tprob[i][k]*b[k][t+1]*eprob[k][seq[t+1]];
                }
                b[i][t] = score;
            }
        }

        double pseq = 0;
        for (int i=0; i<ns; i++) pseq = pseq + a[i][n-1]*b[i][n-1];


        // System.out.println("> " + pseq + " " + (a[0][n-1] + a[1][n-1]));


        for (int t=0; t<n; t++) {
            for (int i=0; i<ns; i++) {
                pp[i][t] = a[i][t]*b[i][t]/pseq;
            }
        }

        // System.out.println();System.out.println("a"); spewey.printmatrix(spewey.transpose(a));
        // System.out.println();System.out.println("b"); spewey.printmatrix(spewey.transpose(b));
        // System.out.println();


        if (numclasses<ns) pp = collapseposteriors(pp,classes, numclasses);
        postprob = pp;
        return pp;

    }


    // in log-space
    public double [][] posteriorl(int [] seq) {
        int n = seq.length;
        double [][] pp = new double[ns][n];

        // forward algorithm
        double [][] a = new double[ns][n];
        a = spewey.logmat(a);

        for (int i=0; i<ns; i++) a[i][0] = liprob[i] + leprob[i][seq[0]];
        for (int t=1; t<n; t++) {
            for (int i=0; i<ns; i++) {
                double score = -1.0/0.0;
                for (int k=0; k<ns; k++) {
                    score = spewey.logeapeb(score, ltprob[k][i] + a[k][t-1]);
                }
                a[i][t] = score + leprob[i][seq[t]];
            }
        }
        double ltotprob = -1.0/0.0;
        for (int i=0; i<ns; i++) {
            // a[i][n-1] = a[i][n-1] + lfprob[i];
            ltotprob = spewey.logeapeb(ltotprob,a[i][n-1] + lfprob[i]);
        }

        lmarginalprob = ltotprob;


        // backward algorithm
        double [][] b = new double[ns][n];

        for (int i=0; i<ns; i++) b[i][n-1] = lfprob[i]; // accounts for termination
        for (int t=n-2; t>=0; t--) {
            for (int i=0; i<ns; i++) {
                double score = -1.0/0.0;
                for (int k=0; k<ns; k++) {
                    score = spewey.logeapeb(score, ltprob[i][k] + b[k][t+1] + leprob[k][seq[t+1]]);
                }
                b[i][t] = score;
            }
        }

        double lpseq = -1.0/0.0;
        for (int i=0; i<ns; i++) lpseq = spewey.logeapeb(lpseq, a[i][0] + b[i][0]);


        // lmarginalprob = lpseq;


        // check this!!
        for (int t=0; t<n; t++) {
            for (int i=0; i<ns; i++) {
                pp[i][t] = Math.exp((a[i][t] + b[i][t]) - lpseq);
            }
        }

        if (numclasses<ns) pp = collapseposteriors(pp,classes, numclasses);
        postprob = pp;
        return pp;

    }


    // in log-space, sparse transmat
    public double [][] sposteriorl(int [] seq) {
        int n = seq.length;
        double [][] pp = new double[ns][n];

        // forward algorithm
        double [][] a = new double[ns][n];
        a = spewey.logmat(a);


        int [][] sps = spewey.sparsesupport(spewey.transpose(tprob),0.000000001);
        int [][] sps2 = spewey.sparsesupport(tprob,0.000000001);

        for (int i=0; i<ns; i++) a[i][0] = liprob[i] + leprob[i][seq[0]];
        for (int t=1; t<n; t++) {
            for (int i=0; i<ns; i++) {
                double score = -1.0/0.0;
                for (int sk=0; sk < sps[i].length; sk++) {
                    int k =sps[i][sk];
                    // for (int k=0; k<ns; k++) {
                    score = spewey.logeapeb(score, ltprob[k][i] + a[k][t-1]);
                }
                a[i][t] = score + leprob[i][seq[t]];
            }
        }
        double ltotprob = -1.0/0.0;
        for (int i=0; i<ns; i++) {
            // a[i][n-1] = a[i][n-1] + lfprob[i];
            ltotprob = spewey.logeapeb(ltotprob,a[i][n-1] + lfprob[i]);
        }

        lmarginalprob = ltotprob;


        // backward algorithm
        double [][] b = new double[ns][n];

        for (int i=0; i<ns; i++) b[i][n-1] = lfprob[i]; // accounts for termination
        for (int t=n-2; t>=0; t--) {
            for (int i=0; i<ns; i++) {
                double score = -1.0/0.0;
                for (int sk=0; sk < sps2[i].length; sk++) {
                    int k =sps2[i][sk];
                    // for (int k=0; k<ns; k++) {
                    score = spewey.logeapeb(score, ltprob[i][k] + b[k][t+1] + leprob[k][seq[t+1]]);
                }
                b[i][t] = score;
            }
        }

        double lpseq = -1.0/0.0;
        for (int i=0; i<ns; i++) lpseq = spewey.logeapeb(lpseq, a[i][0] + b[i][0]);


        // lmarginalprob = lpseq;

        // check this!!
        for (int t=0; t<n; t++) {
            for (int i=0; i<ns; i++) {
                pp[i][t] = Math.exp((a[i][t] + b[i][t]) - lpseq);
            }
        }

        if (numclasses<ns) pp = collapseposteriors(pp,classes, numclasses);
        postprob = pp;
        return pp;

    }


    // wrapper
    public double [][] sposteriors(int [] seq) {
        return sposteriors(seq, tprob, eprob, iprob);
    }



    // wrapper
    public double [][] posteriors(int [] seq) {
        return posteriors(seq, tprob, eprob, iprob);
    }



    // use scaling to prevent underflow

    public double [][] posteriors(int [] seq, double [][] tprob, double [][] eprob, double [] iprob) {

        int n = seq.length;
        double [][] pp = new double[ns][n];

        // scaling factor
        double [] sc = new double[n];
        for (int i=0; i<n; i++) sc[i] = 1.0;

        // forward algorithm
        double [][] a = new double[ns][n];

        for (int i=0; i<ns; i++) a[i][0] = sc[0]*iprob[i]*eprob[i][seq[0]];
        for (int t=1; t<n; t++) {
            double sf = 0.0;
            for (int i=0; i<ns; i++) {
                double score = 0;
                for (int k=0; k<ns; k++) {
                    score = score + tprob[k][i]*a[k][t-1];
                }
                a[i][t] = score*eprob[i][seq[t]];
                sf = sf + score*eprob[i][seq[t]];
            }
            sc[t] = 1/sf;
            for (int i=0; i<ns; i++) a[i][t] = a[i][t]*sc[t];
        }

        // backward algorithm
        double [][] b = new double[ns][n];

        for (int i=0; i<ns; i++) b[i][n-1] = sc[n-1]; // 1; // iprob[i];
        for (int t=n-2; t>=0; t--) {
            for (int i=0; i<ns; i++) {
                double score = 0;
                for (int k=0; k<ns; k++) {
                    score = score + tprob[i][k]*b[k][t+1]*eprob[k][seq[t+1]];
                }
                b[i][t] = sc[t]*score;
            }
        }

        double pseq = 0;
        for (int i=0; i<ns; i++) pseq = pseq + a[i][n-1]*b[i][n-1];
        double lpseq = Math.log(pseq);
        for (int t=0; t<n; t++) lpseq = lpseq - Math.log(sc[t]);


        for (int t=0; t<n; t++) {
            double denom = 0;
            for (int i=0; i<ns; i++) {
                denom = denom + a[i][t]*b[i][t];
            }
            for (int i=0; i<ns; i++) {
                pp[i][t] = a[i][t]*b[i][t]/denom;
            }
        }

        System.out.println();System.out.println("a"); spewey.printmatrix(spewey.transpose(a));
        System.out.println();System.out.println("b"); spewey.printmatrix(spewey.transpose(b));
        // System.out.println();


        lmarginalprob = lpseq;
        lmp = lpseq;
        if (numclasses<ns) pp = collapseposteriors(pp,classes, numclasses);
        postprob = pp;
        return pp;

    }


    // sparse transmat; use scaling to prevent underflow

    public double [][] sposteriors(int [] seq, double [][] tprob, double [][] eprob, double [] iprob) {

        int n = seq.length;
        double [][] pp = new double[ns][n];

        // scaling factor
        double [] sc = new double[n];
        for (int i=0; i<n; i++) sc[i] = 1.0;

        // forward algorithm
        double [][] a = new double[ns][n];


        int [][] sps = spewey.sparsesupport(spewey.transpose(tprob),0.000000001);
        int [][] sps2 = spewey.sparsesupport(tprob,0.000000001);

        for (int i=0; i<ns; i++) a[i][0] = sc[0]*iprob[i]*eprob[i][seq[0]];
        for (int t=1; t<n; t++) {
            double sf = 0.0;
            for (int i=0; i<ns; i++) {
                double score = 0;
                for (int sk=0; sk < sps[i].length; sk++) {
                    int k = sps[i][sk];
                    // for (int k=0; k<ns; k++) {
                    score = score + tprob[k][i]*a[k][t-1];
                }
                a[i][t] = score*eprob[i][seq[t]];
                sf = sf + score*eprob[i][seq[t]];
            }
            sc[t] = 1/sf;
            for (int i=0; i<ns; i++) a[i][t] = a[i][t]*sc[t];
        }

        double totprob = 0;

        for (int i=0; i<ns; i++) {
            totprob = totprob + a[i][n-1]*fprob[i];
        }

        double ltotprob = Math.log(totprob);
        for (int t=0; t<n; t++) ltotprob = ltotprob - Math.log(sc[t]);

        lmarginalprob = ltotprob;


        // backward algorithm
        double [][] b = new double[ns][n];

        for (int i=0; i<ns; i++) b[i][n-1] = sc[n-1]; // 1; // iprob[i];
        for (int t=n-2; t>=0; t--) {
            for (int i=0; i<ns; i++) {
                double score = 0;
                for (int sk=0; sk < sps2[i].length; sk++) {
                    // for (int k=0; k<ns; k++) {
                    int k = sps2[i][sk];
                    score = score + tprob[i][k]*b[k][t+1]*eprob[k][seq[t+1]];
                }
                b[i][t] = sc[t]*score;
            }
        }

        double pseq = 0;
        for (int i=0; i<ns; i++) pseq = pseq + a[i][n-1]*b[i][n-1];
        double lpseq = Math.log(pseq);
        for (int t=0; t<n; t++) lpseq = lpseq - Math.log(sc[t]);


        for (int t=0; t<n; t++) {
            double denom = 0;
            for (int i=0; i<ns; i++) {
                denom = denom + a[i][t]*b[i][t];
            }
            for (int i=0; i<ns; i++) {
                pp[i][t] = a[i][t]*b[i][t]/denom;
            }
        }

        // System.out.println();System.out.println("a"); spewey.printmatrix(spewey.transpose(a));
        // System.out.println();System.out.println("b"); spewey.printmatrix(spewey.transpose(b));
        // System.out.println();

        // lmarginalprob = lpseq;
        lmp = lpseq;
        if (numclasses<ns) pp = collapseposteriors(pp,classes, numclasses);
        postprob = pp;
        return pp;

    }


    // also returns forward and backward values
    public double [][][] posteriorsfb(int [] seq, double [][] tprob, double [][] eprob, double [] iprob) {

        int n = seq.length;
        double [][] pp = new double[ns][n];

        // scaling factor
        double [] sc = new double[n];
        for (int i=0; i<n; i++) sc[i] = 1.0;

        // forward algorithm
        double [][] a = new double[ns][n];

        for (int i=0; i<ns; i++) a[i][0] = sc[0]*iprob[i]*eprob[i][seq[0]];
        for (int t=1; t<n; t++) {
            double sf = 0.0;
            for (int i=0; i<ns; i++) {
                double score = 0;
                for (int k=0; k<ns; k++) {
                    score = score + tprob[k][i]*a[k][t-1];
                }
                a[i][t] = score*eprob[i][seq[t]];
                sf = sf + score*eprob[i][seq[t]];
            }
            sc[t] = 1/sf;
            for (int i=0; i<ns; i++) a[i][t] = a[i][t]*sc[t];
        }

        // backward algorithm
        double [][] b = new double[ns][n];

        for (int i=0; i<ns; i++) b[i][n-1] = sc[n-1]; // 1; // iprob[i];
        for (int t=n-2; t>=0; t--) {
            for (int i=0; i<ns; i++) {
                double score = 0;
                for (int k=0; k<ns; k++) {
                    score = score + tprob[i][k]*b[k][t+1]*eprob[k][seq[t+1]];
                }
                b[i][t] = sc[t]*score;
            }
        }

        double pseq = 0;
        for (int i=0; i<ns; i++) pseq = pseq + a[i][n-1]*b[i][n-1];
        double lpseq = Math.log(pseq);
        for (int t=0; t<n; t++) lpseq = lpseq - Math.log(sc[t]);


        for (int t=0; t<n; t++) {
            double denom = 0;
            for (int i=0; i<ns; i++) {
                denom = denom + a[i][t]*b[i][t];
            }
            for (int i=0; i<ns; i++) {
                pp[i][t] = a[i][t]*b[i][t]/denom;
            }
        }

        // System.out.println();System.out.println("a"); spewey.printmatrix(spewey.transpose(a));
        // System.out.println();System.out.println("b"); spewey.printmatrix(spewey.transpose(b));
        // System.out.println();

        lmarginalprob = lpseq;

        lmp = lpseq;


        postprob = pp;
        // return pp;

        double [][][] pfb = new double[3][ns][n];
        pfb[0] = pp;
        pfb[1] = a;
        pfb[2] = b;
        return pfb;

    }


    // allow partial labels
    public double [][] posteriors(int [] seq, int [] lab, double [][] tprob, double [][] eprob, double [] iprob) {

        int n = seq.length;
        double [][] pp = new double[ns][n];

        // scaling factor
        double [] sc = new double[n];
        for (int i=0; i<n; i++) sc[i] = 1.0;

        // forward algorithm
        double [][] a = new double[ns][n];

        int lb, ub;

        if (lab[0]<0) { for (int i=0; i<ns; i++) a[i][0] = sc[0]*iprob[i]*eprob[i][seq[0]]; }
        else a[lab[0]][0] = sc[0]*eprob[lab[0]][seq[0]];

        for (int t=1; t<n; t++) {
            double sf = 0.0;

            if (lab[t]<0) {
                for (int i=0; i<ns; i++) {
                    double score = 0;
                    for (int k=0; k<ns; k++) {
                        score = score + tprob[k][i]*a[k][t-1];
                    }
                    a[i][t] = score*eprob[i][seq[t]];
                    sf = sf + score*eprob[i][seq[t]];
                }
            }
            else {
                double score = 0;
                for (int k=0; k<ns; k++) {
                    score = score + tprob[k][lab[t]]*a[k][t-1];
                }
                a[lab[t]][t] = score*eprob[lab[t]][seq[t]];
                sf = sf + score*eprob[lab[t]][seq[t]];
            }

            sc[t] = 1/sf;
            for (int i=0; i<ns; i++) a[i][t] = a[i][t]*sc[t];
        }

        // backward algorithm
        double [][] b = new double[ns][n];

        for (int i=0; i<ns; i++) b[i][n-1] = sc[n-1]; // 1; // iprob[i];
        for (int t=n-2; t>=0; t--) {

            if (lab[t+1]<0) {
                for (int i=0; i<ns; i++) {
                    double score = 0;
                    for (int k=0; k<ns; k++) {
                        score = score + tprob[i][k]*b[k][t+1]*eprob[k][seq[t+1]];
                    }
                    b[i][t] = sc[t]*score;
                }
            }
            else {
                for (int i=0; i<ns; i++) {
                    double score = 0;
                    score = score + tprob[i][lab[t+1]]*b[lab[t+1]][t+1]*eprob[lab[t+1]][seq[t+1]];
                    b[i][t] = sc[t]*score;
                }
            }


        }

        double pseq = 0;
        for (int i=0; i<ns; i++) pseq = pseq + a[i][n-1]*b[i][n-1];
        double lpseq = Math.log(pseq);
        for (int t=0; t<n; t++) lpseq = lpseq - Math.log(sc[t]);


        for (int t=0; t<n; t++) {
            double denom = 0;
            for (int i=0; i<ns; i++) {
                denom = denom + a[i][t]*b[i][t];
            }
            for (int i=0; i<ns; i++) {
                pp[i][t] = a[i][t]*b[i][t]/denom;
            }
        }

        // System.out.println();System.out.println("a"); spewey.printmatrix(spewey.transpose(a));
        // System.out.println();System.out.println("b"); spewey.printmatrix(spewey.transpose(b));
        // System.out.println();

        lmarginalprob = lpseq;
        lmp = lpseq;

        if (numclasses<ns) pp = collapseposteriors(pp,classes, numclasses);


        postprob = pp;
        return pp;

    }


    public double logprobsubtrellis(int [] seq, int [] st) {

        int n = seq.length;
        double [][] pp = new double[ns][n];

        // scaling factor
        double [] sc = new double[n];
        for (int i=0; i<n; i++) sc[i] = 1.0;


        double [][] a = new double[ns][n];

        // forward algorithm on entire trellis
        for (int i=0; i<ns; i++) {
            if (true) a[i][0] = sc[0]*iprob[i]*eprob[i][seq[0]];
        }
        for (int t=1; t<n; t++) {
            double sf = 0.0;
            for (int i=0; i<ns; i++) {
                if (true) {
                    double score = 0;
                    for (int k=0; k<ns; k++) {
                        score = score + tprob[k][i]*a[k][t-1];
                    }
                    a[i][t] = score*eprob[i][seq[t]];
                    sf = sf + score*eprob[i][seq[t]];
                }
            }
            sc[t] = 1/sf;
            for (int i=0; i<ns; i++) a[i][t] = a[i][t]*sc[t];
        }

        double lpd = 0; // pr(seq)
        for (int t=0; t<n; t++) lpd = lpd - Math.log(sc[t]);

        a = new double[ns][n];

        //  forward algorithm on subtrellis
        for (int i=0; i<ns; i++) {
            if (st[i] == 1) a[i][0] = sc[0]*iprob[i]*eprob[i][seq[0]];
        }
        for (int t=1; t<n; t++) {
            double sf = 0.0;
            for (int i=0; i<ns; i++) {
                if (st[i] == 1) {
                    double score = 0;
                    for (int k=0; k<ns; k++) {
                        score = score + tprob[k][i]*a[k][t-1];
                    }
                    a[i][t] = score*eprob[i][seq[t]];
                    sf = sf + score*eprob[i][seq[t]];
                }
            }
            sc[t] = 1/sf;
            for (int i=0; i<ns; i++) a[i][t] = a[i][t]*sc[t];
        }

        double lpdas = 0; // pr(seq & subtrellis)
        for (int t=0; t<n; t++) lpdas = lpdas - Math.log(sc[t]);


        // not needed
        // prior probability of subtrellis

        /*

           a = new double[ns][n];

           for (int i=0; i<ns; i++) {
           if (st[i] == 1) a[i][0] = sc[0]*iprob[i];
           }
           for (int t=1; t<n; t++) {
           double sf = 0.0;
           for (int i=0; i<ns; i++) {
           if (st[i] == 1) {
           double score = 0;
           for (int k=0; k<ns; k++) {
           score = score + tprob[k][i]*a[k][t-1];
           }
           a[i][t] = score;
           sf = sf + score;
           }
           }
           sc[t] = 1/sf;
           for (int i=0; i<ns; i++) a[i][t] = a[i][t]*sc[t];
           }

           double lps = 0; // pr(subtrellis)
           for (int t=0; t<n; t++) lps = lps - Math.log(sc[t]);
           */


        // System.out.println((float) lpdas + "\t" + (float) lpd + " \t" + (float) Math.exp(lpdas - lpd));
        return (lpdas - lpd);
    }


    public int [] mapdecode(int [] seq) {
        int n = seq.length;
        int [] map = new int[n];
        double [][] pp = posterior(seq);
        postprob = pp;
        for (int i=0; i<n; i++) {
            for (int j=0; j<pp.length; j++) {
                if (pp[j][i] > pp[map[i]][i]) map[i]=j;
            }
        }
        // if (numclasses < ns) map = spewey.mapseq(map,classes);
        mappath = map;
        return map;
    }

    // pp is ns x n
    public static int [] mapdecode(double [][] pp) {
        int n = pp[0].length;
        int ns = pp.length;
        int [] map = new int[n];
        for (int i=0; i<n; i++) {
            for (int j=0; j<pp.length; j++) {
                if (pp[j][i] > pp[map[i]][i]) map[i]=j;
            }
        }
        return map;
    }


    public int [] smapdecodel(int [] seq) {
        int n = seq.length;
        int [] map = new int[n];
        double [][] pp = sposteriorl(seq);
        postprob = pp;
        for (int i=0; i<n; i++) {
            for (int j=0; j<pp.length; j++) {
                if (pp[j][i] > pp[map[i]][i]) map[i]=j;
            }
        }
        // if (numclasses < ns) map = spewey.mapseq(map,classes);
        mappath = map;
        return map;
    }



    public int [] mapdecodel(int [] seq) {
        int n = seq.length;
        int [] map = new int[n];
        double [][] pp = posteriorl(seq);
        postprob = pp;
        for (int i=0; i<n; i++) {
            for (int j=0; j<pp.length; j++) {
                if (pp[j][i] > pp[map[i]][i]) map[i]=j;
            }
        }
        // if (numclasses < ns) map = spewey.mapseq(map,classes);
        mappath = map;
        return map;
    }


    public int [] mapdecodes(int [] seq) {
        int n = seq.length;
        int [] map = new int[n];
        double [][] pp = posteriors(seq);
        postprob = pp;
        for (int i=0; i<n; i++) {
            for (int j=0; j<pp.length; j++) {
                if (pp[j][i] > pp[map[i]][i]) map[i]=j;
            }
        }
        // if (numclasses < ns) map = spewey.mapseq(map,classes);
        mappath = map;
        return map;
    }



    public double marginal(int [] seq) {
        int n = seq.length;
        double mp = 0;


        return mp;
    }


    public static double [][] collapseposteriors(double [][] pp, int [] mask, int newns) {
        double [][] npp = new double[newns][pp[0].length];
        for(int i=0; i<pp.length; i++) {
            for (int j=0; j<pp[0].length; j++) {
                npp[mask[i]][j]+=pp[i][j];
            }
        }
        return npp;


    }


    public int [][] simulate(int n, Random rg) {
        int [] hseq = new int[n];
        int [] oseq = new int[n];
        double r = rg.nextDouble();
        int dex = 0;
        while (r < ciprob[dex]) dex++;
        hseq[0]=dex;
        for (int i=1; i<n; i++) {
            r = rg.nextDouble();
            dex = 0;
            while (r < ctprob[hseq[i-1]][dex]) dex++;
            hseq[i] = dex;
        }
        for (int i=0; i<n; i++) {
            r = rg.nextDouble();
            dex = 0;
            while (r < ceprob[hseq[i]][dex]) dex++;
            oseq[i] = dex;

        }
        int [][] smat = new int[2][n];
        smat[0] = hseq;
        smat[1] = oseq;
        return smat;

    }

    public void printme(String st) {
        int nc = postprob.length;

        for (int i = 0; i<viterbipath.length; i++) {
            System.out.print(st + "\t" + viterbipath[i] + "\t" + mappath[i]);
            for (int j = 0; j<nc; j++) System.out.print("\t" + (float) postprob[j][i]);
            System.out.println();
        }
        System.out.print("# " +(float) etst + "\t" + (float) lpst);
        for (int i=0; i<nc; i++) System.out.print("\t" + subtrellis[i]);
        System.out.println();
        System.out.print("# " +(float) etst + "\t" + (float) lpst);
        for (int i=0; i<nc; i++) System.out.print("\t" + (float) spewey.sum(postprob[i]));
        System.out.println();
    }

    public void printme(String st, int [] seq) {
        int nc = postprob.length;

        for (int i = 0; i<viterbipath.length; i++) {
            System.out.print(st + "\t" + seq[i] + "\t" + viterbipath[i] + "\t" + mappath[i]);
            for (int j = 0; j<nc; j++) System.out.print("\t" + (float) postprob[j][i]);
            System.out.println();
        }
        System.out.print("# " +(float) etst + "\t" + (float) lpst);
        for (int i=0; i<nc; i++) System.out.print("\t" + subtrellis[i]);
        System.out.println();
        System.out.print("# " +(float) etst + "\t" + (float) lpst);
        for (int i=0; i<nc; i++) System.out.print("\t" + (float) spewey.sum(postprob[i]));
        System.out.println();
    }



    public double [][][] scoreseqs(int [][] seqs, int [][] labels, double [][] tprob, double [][] eprob, double [] iprob) {
        // double [][] tprob;
        // double [][] eprob;
        // double [] iprob;
        int ns = iprob.length;

        int n = seqs.length;
        double [][][] scores = new double[n][ns][];
        for (int i=0; i<n; i++) {
            scores[i] = posteriors(seqs[i], labels[i], tprob, eprob, iprob);
        }

        for (int i=0; i<n; i++) {
            spewey.printmatrix(spewey.transpose(scores[i]));
            System.out.println("-------------------------------------------------------------");
        }
        return scores;

    }


    public double [][][] scoreseqs(int [][] seqs, int [][] labels) {
        // double [][] tprob;
        // double [][] eprob;
        // double [] iprob;
        int ns = iprob.length;

        int n = seqs.length;
        double [][][] scores = new double[n][ns][];
        for (int i=0; i<n; i++) {
            scores[i] = posteriors(seqs[i], labels[i], tprob, eprob, iprob);
        }

        for (int i=0; i<n; i++) {
            spewey.printmatrix(spewey.transpose(scores[i]));
            System.out.println("-------------------------------------------------------------");
        }
        return scores;

    }



    public void dottify(boolean record) {
        System.out.println("Digraph G {");
        for (int i=0; i<ns; i++) {
            System.out.println("  n" + i + " [label=\"" + names[i] + "\", shape=ellipse];");
        }
        for (int i=0; i<ns; i++) {
            for (int j=0; j<ns; j++) {
                if (tprob[i][j]>0) {
                    String lab = new String("" + (float) tprob[i][j]);
                    if (lab.length()>7) lab = lab.substring(0,6);
                    System.out.println("  n" + i + " -> n" + j + " [label=\"" + lab + "\", color=gray];");
                }
            }
        }

        if (record) {

            System.out.print("rec [shape=record, label=\"");

            System.out.print("{	  " );
            for (int j=0; j<no; j++) {
                String lab = new String("" + spewey.aanames[j]);
                if (lab.length()>7) lab = lab.substring(0,6);
                System.out.print("|" + lab);
            }
            System.out.print("}");

            for (int i=0; i<ns; i++) {
                System.out.print("|{" + names[i]);
                for (int j=0; j<no; j++) {
                    String lab = new String("" + (float) eprob[i][j]);
                    if (lab.length()>7) lab = lab.substring(0,6);

                    System.out.print("|" + lab);
                }
                System.out.print("}");
            }
            System.out.println("\"];");
        }


        System.out.println("}");
    }




}


///////////////////////////////////////////////////////////////////////////////
// clumsy, since need to call hasmorefastas() before nextname() before nextfasta()

class fastareader {

    String name;
    StringBuffer sb;
    BufferedReader in;
    String line;
    boolean ondeck = false;
    int numread = 0;

    boolean goodfile;



    fastareader(String filename) {
        goodfile = true;
        try {
            in = new BufferedReader(new FileReader(filename));
        }
        catch (IOException e) {
            try {
                in = new BufferedReader(new InputStreamReader(fastareader.class.getResourceAsStream(filename)));
            }
            catch (Exception ie) {
                System.out.println("# Couldn't open " + filename);
                goodfile=false;
            }
        }
    }


    public StringBuffer nextfasta() {
        sb = new StringBuffer("");
        try {
            while ((line = in.readLine()) != null) {
                line.trim();
                if (line.length() == 0) {
                    ondeck = false;
                    return sb;
                }
                else if (line.charAt(0) == '>') {
                    ondeck = true;
                    name = line.substring(1);
                    name.trim();
                    return sb;
                }
                else sb.append(line);
            }
            ondeck = false;
            return sb;
        }
        catch (IOException e) {
            System.out.println("# Couldn't get next fasta");
            return sb;
        }
    }

    public String nextname() {
        numread++;
        return name;
    }

    public boolean hasmorefastas() {
        if (ondeck) return true;
        try {
            while ((line = in.readLine()) != null) {
                if (line.length() >0 && line.charAt(0) == '>') {
                    name = line.trim().substring(1);
                    return true;
                }
            }
            in.close();
            return false;
        }
        catch (IOException e) {
            System.out.println("# Couldn't read lines");
            return false;
        }
    }

}


///////////////////////////////////////////////////
///////////////////////////////// harrison gerstein lps

class hgalg {


    int [][] masks = {
        {0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0},	// QN
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},	// Q
        {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},	// N
        {0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0},	// DERK
        {0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,1,0,0,0},	// VILM
        {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},	// G
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},	// Y
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}	// S
    };

    double [] threshs  = {
        Math.log(1.8e-14),
        Math.log(1e-13),
        Math.log(1e-13),
        Math.log(6.5e-3),
        Math.log(2e-2),
        Math.log(5e-4),
        Math.log(5e-4),
        Math.log(5e-4)
    };

    double [] ps = {
        0.1,
        0.05,
        0.05,
        0.2,
        0.2,
        0.05,
        0.05,
        0.05
    };


    int [] tails = {
        1,
        1,
        1,
        -1,
        -1,
        1,
        1,
        1
    };



    double [][][] luts;


    int maxlen;
    int minlen = 25;



    hgalg(int maxlen, double [] bgfreq) {

        this.maxlen = maxlen;
        if (bgfreq !=null) {
            for (int i=0; i<masks.length; i++) {
                ps[i]=0.0;
                for (int j=0; j<22; j++) {
                    if (masks[i][j]==1) ps[i]+=bgfreq[j];
                }
            }
        }


        luts = new double[masks.length][][]; // need other dims?
        for (int i=0; i<masks.length; i++) {
            luts[i] = spewey.logbinomiallut2(maxlen, ps[i]);
        }
    }



    public double [][] check(int [] seq) {
        double [][] bias = new double [masks.length][4];
        for (int i=0; i<masks.length; i++) {
            bias[i] =  spewey.lps(spewey.mapseq(seq,masks[i]),	 ps[i],	 luts[i],  minlen, maxlen, tails[i]);
        }
        return bias;
    }



    public double [][] check2(int [] seq) {
        double [][] bias = new double [masks.length][4];
        double [] leng = new double[masks.length];
        int lb = 0;
        int ub = 0;
        int minol = 15;
        int longest = 0;
        for (int i=0; i<=2; i++) {
            bias[i] =  spewey.lps(spewey.mapseq(seq,masks[i]),	 ps[i],	 luts[i],   minlen,  maxlen, tails[i]);
            leng[i] = bias[i][1]-bias[i][0]+1;
            if (leng[i]>longest) {longest = (int) leng[i]; lb = (int) bias[i][0]; ub = (int) bias[i][1];}
        }
        // merge these if they overlap by at least 15?
        for (int i=0; i<=2; i++) {
            if (overlap(lb,ub,(int) bias[i][0], (int) bias[i][1]) > minol) {lb=Math.min(lb, (int) bias[i][0]); ub=Math.max(ub, (int) bias[i][1]);};
        }
        // bias[3][0] = lb;
        // bias[3][1] = ub;
        int [] subseq = spewey.submatrix(seq,lb,ub);

        for (int i=3; i< masks.length; i++) {
            bias[i] =  spewey.lps(spewey.mapseq(subseq,masks[i]),   ps[i],  luts[i],  subseq.length,  maxlen, tails[i]);
            bias[i][0] += lb; bias[i][1] += lb;
            leng[i] = bias[i][1]-bias[i][0]+1;
            // if (leng[i]>longest) {longest = (int) leng[i]; lb = (int) bias[i][0]; ub = (int) bias[i][1];}
        }

        // System.out.println(spewey.aa2string(subseq));
        return bias;

    }




    public boolean winner (int [] seq) {
        double [][] bias = new double [masks.length][4];
        double [] leng = new double[masks.length];
        int slb = 0;
        int sub = 0;
        int minol = 15;
        int longest = 0;
        for (int i=0; i<=2; i++) {
            bias[i] =  spewey.lps(spewey.mapseq(seq,masks[i]),	 ps[i],	 luts[i],   minlen,  maxlen, tails[i]);
            leng[i] = bias[i][1]-bias[i][0]+1;
            if (bias[i][3] <= threshs[i] && leng[i]>longest) {longest = (int) leng[i]; slb = (int) bias[i][0]; sub = (int) bias[i][1];}
        }


        if (longest == 0) return false;


        // merge these if they overlap by at least 15?
        int lb = slb;
        int ub = sub;
        for (int i=0; i<=2; i++) {
            if (overlap(lb,ub,(int) bias[i][0], (int) bias[i][1]) > minol) {lb=Math.min(lb, (int) bias[i][0]); ub=Math.max(ub, (int) bias[i][1]);};
        }
        // bias[3][0] = lb;
        // bias[3][1] = ub;
        int [] subseq = spewey.submatrix(seq,slb,sub);

        for (int i=3; i< masks.length; i++) { // uses entire sequence here!!!
            bias[i] =  spewey.lps(spewey.mapseq(subseq,masks[i]),   ps[i],  luts[i],  subseq.length,  maxlen, tails[i]);
            bias[i][0] += slb; bias[i][1] += slb;
            leng[i] = bias[i][1]-bias[i][0]+1;
        }

        boolean good = true;

        if (bias[3][3] >= threshs[3] || bias[4][3] >= threshs[4]) good = false;
        if (bias[5][3] >= threshs[5] && bias[6][3] >= threshs[6] && bias[7][3] >= threshs[7] ) good = false;

        if (good) {
            // spewey.printmatrix(bias);
            return true;
        }

        else {
            int [] seq1 = spewey.submatrix(seq,0,lb-1);
            int [] seq2 = spewey.submatrix(seq,ub+1,seq.length);
            if (seq1.length >= minlen && seq2.length >= minlen) return (winner(seq1) || winner(seq2));
            else if (seq1.length >= minlen) return (winner(seq1));
            else if (seq2.length >= minlen) return (winner(seq2));
            else return false;

        }
    }




    public static int overlap(int l1, int u1, int l2, int u2) {
        if (u2 < l1 || u1 <l2) return 0;
        else if (l2 <= l1 && u2 >=u1) return u1-l1+1;
        else if (l1 <= l2 && u1 >=u2) return u2-l2+1;
        else if (l1 <= l2)  return u1-l2+1;
        else if (l2 <= l1)  return u2-l1+1;
        else return 1;
    }


}

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

class disorderreport {

    // double c1 = 2.785;
    // double c2 = -1.0;
    // double c3 = -1.151;
    double [] ccdef = { 2.785, -1.0, -1.151 };
    double [] cc;
    double [] hydro;
    double [] charge;
    double [] combo;
    double [] other;
    double [] other2;
    int [] aa;
    int ww;
    int ww2;
    boolean reflect;
    double meancharge;
    double meanhydro;
    double meancombo;
    int n;
    int numseg;
    int numdisordered;
    int numdisorderedstrict; // only windows completely contained in seq
    int numdisorderedstrict2; // only windows completely contained in seq and of min length
    int [] startaa;
    int [] stopaa;
    int [] lenaa;
    double [] localmean;
    double [] localsd;
    double [] localhydro;
    double [] localcharge;
    int maxlen;		// for combo
    int minlen = 5;	// for combo
    int minlen2 = 40;	// for combo
    int bestdex;
    double maxlong;
    int targetlen = 40;
    double botl = 0.0; // bestoftargetlength;
    double [] hssr;
    double [] hssr2;

    double [] maa3;
    double [] maa4;

    int minlength = 80; // for hssr2
    int maxlength = 240; // for hssr2


    double [] otherweights;   // for llk
    double [] otherweights2;  // for ross1 or ross2

    double [] agg; // not currently used

    double rossmaxprd;	  // maximum ross score
    double rossmaxdis;	  // corresponding FI
    int rossmaxcenter;	  // index of max above
    //boolean smushpp = false; // give neutral scores to all but first of consectuive prolines.
    boolean smushpp = true;
    // should redo to better handle case where window doesn't include first proline

    disorderreport(int [] aa, int ww, int ww2, boolean reflect, double [] cc, double [] otherweights, double [] otherweights2) {
        this.ww = ww;
        this.ww2 = ww2;
        this.cc = cc;
        this.reflect = reflect;
        this.aa = aa;
        this.otherweights = otherweights;
        this.otherweights = otherweights2;
        n = aa.length;
        double [] maa1 = spewey.mapseq(aa,spewey.axpb(1.0/9,spewey.aahydro,0.5));
        meanhydro = spewey.mean(maa1);
        hydro = spewey.slidingaverage2(maa1,ww,reflect);
        double [] maa2 = spewey.mapseq(aa,spewey.aacharge);
        meancharge = spewey.mean(maa2);
        charge = spewey.slidingaverage2(maa2,ww,false); // don't reflect charges
        meancombo  = cc[2] + cc[1]*Math.abs(meancharge) + cc[0]*meanhydro;
        combo = spewey.axpbypc(cc[0],hydro,cc[1],spewey.absval(charge),cc[2]);
        maa3 = spewey.mapseq(aa,otherweights);
        other = spewey.slidingaverage2(maa3,ww,reflect);
        maa4 = spewey.mapseq(aa,otherweights2);
        if (smushpp) {maa4 = pp2p(maa4,aa);} // gives Ps after first score of 0
        other2 = spewey.slidingaverage2(maa4,ww2,false);
        // average of averages for other2
        other2 = spewey.slidingaverage2(other2,ww2,false);
        // combo = spewey.slidingaverage2(combo,ww,false);

        numdisordered = 0;
        rossmaxprd = -1000;
        rossmaxdis = -1000;
        rossmaxcenter = -1;
        for (int k=0; k<n; k++) {
            if (combo[k]<0) {
                numdisordered++;
            }
        }


        numdisorderedstrict = 0;
        int halfw = (ww-1)/2; // or 0
        if (halfw > n/2) {halfw = n/2;}
        // System.out.println("halfww " + halfw);
        if (combo[halfw]<0) {
            numdisorderedstrict =  numdisorderedstrict + halfw + 1;
        }
        if (combo[n-halfw-1]<0) {
            numdisorderedstrict =  numdisorderedstrict + halfw + 1;
        }
        for (int k = halfw + 1; k < n-halfw - 1; k++) {
            if (combo[k]<0) {
                numdisorderedstrict++;
            }
        }

        // restrict to MAP parse?
        double rossmaxscore = -1000;
        double rossscore = 0;
        for (int k=((ww2-1)/2); k<(n-(ww2-1)/2); k++) {
            rossscore = other2[k]; // Math.min(-1.0*combo[k]-0.0,other2[k]-0.05);
            // rossscore = Math.min(-1.0*combo[k]-0.0,other2[k]-0.05);
            // if ((rossscore > rossmaxscore) & combo[k] < 0) { /// what should score be if no disorder?
            if (rossscore > rossmaxscore) {
                rossmaxscore = rossscore;
                rossmaxprd = other2[k];
                rossmaxdis = combo[k];
                rossmaxcenter = k;
            }

        }

        // maa4 = spewey.mapseq(aa,spewey.aadobson);
        //  agg = dobson.sa(aa); // spewey.slidingaverage(maa4,7,true);

        /////////////////////////////

        hssr = spewey.hss2(maa3);
        if (hssr[1]-hssr[0]+1 < minlength || hssr[1]-hssr[0]+1 > maxlength) {
            hssr2 = spewey.hss2(maa3, minlength, maxlength);
        }
        else hssr2 = hssr;


        /////////////////////////////
        int i = 0 + halfw; // +1?
        startaa = new int[n/minlen];
        stopaa = new int[n/minlen];
        lenaa = new int[n/minlen];
        localmean = new double[n/minlen];
        localsd = new double[n/minlen];
        localhydro = new double[n/minlen];
        localcharge = new double[n/minlen];
        if (localhydro.length>0) localhydro[0] = 0.5;
        numseg = 0;
        numdisorderedstrict2 = 0;
        while (i < n - halfw) {
            if (combo[i] < 0) {
                double sc = combo[i];
                double lc = Math.abs(charge[i]);
                double lh = hydro[i];
                double scsc = combo[i]*combo[i];
                int startdex = i;
                i++;
                while (i < n - halfw && combo[i] <0) {
                    sc = sc + combo[i];
                    lh = lh + hydro[i];
                    lc = lc + Math.abs(charge[i]);
                    scsc = scsc + combo[i]*combo[i];
                    i++;
                }
                int stopdex = i-1;
                int len = stopdex-startdex+1;
                double msc = sc/len;
                double sdsc = Math.sqrt(scsc/len - msc*msc);
                if (startdex == halfw) {startdex = 0;}
                if (stopdex == n - halfw -1 ) {stopdex = n-1;}
                len = stopdex-startdex+1;
                // adjust length if startdex == halfw or stopdex == n - halfww - 1
                if (len >= minlen) {
                    startaa[numseg] = startdex;
                    stopaa[numseg] = stopdex;
                    lenaa[numseg] = len;
                    localmean[numseg] = msc;
                    localsd[numseg] = sdsc;
                    localhydro[numseg] = lh/len;
                    localcharge[numseg] = lc/len;
                    numseg++;
                    numdisorderedstrict2 = numdisorderedstrict2 + len;
                }
                // now trim arrays to length numseg?
            }
            else i++;
        }
        if (lenaa.length>0) maxlen = spewey.maxint(lenaa);
        maxlong = 1;
        bestdex = 0;
        for (i=0; i<numseg; i++) {if (lenaa[i] >= minlen2 && localmean[i]<maxlong) {maxlong = localmean[i]; bestdex = i;};}
    }

    // gives p weight to only first occurance in consecutive run.
    // should redo to better handle case where window doesn't include first proline
    public static double [] pp2p(double [] oldwtvec, int [] aavec) {
        double [] newwtvec = new double[oldwtvec.length];
        newwtvec[0] = oldwtvec[0];
        for (int k=1; k<oldwtvec.length; k++) {
            if ((aavec[k]==13) & (aavec[k-1]==13)) newwtvec[k]=0;
            else newwtvec[k]=oldwtvec[k];
            // System.out.println(aavec[k] + "\t" + newwtvec[k] + "\t" + oldwtvec[k]);
        }
        return newwtvec;
    }


    // new linear combination and threshhold;
    public void reweight(double [] cc) {
        meancombo  = cc[2] + cc[1]*Math.abs(meancharge) + cc[0]*meanhydro;
        combo = spewey.axpbypc(cc[0],hydro,cc[1],spewey.absval(charge),cc[2]);
        int i=0;
        startaa = new int[n/minlen];
        stopaa = new int[n/minlen];
        lenaa = new int[n/minlen];
        localmean = new double[n/minlen];
        localsd = new double[n/minlen];
        localhydro = new double[n/minlen];
        localcharge = new double[n/minlen];
        if (localhydro.length>0) localhydro[0] = 0.5;
        numseg = 0;
        while (i<n) {
            if (combo[i] < 0) {
                double sc = combo[i];
                double lc = Math.abs(charge[i]);
                double lh = hydro[i];
                double scsc = combo[i]*combo[i];
                int startdex = i;
                i++;
                while (i < n && combo[i] <0) {
                    sc = sc + combo[i];
                    lh = lh + hydro[i];
                    lc = lc + Math.abs(charge[i]);
                    scsc = scsc + combo[i]*combo[i];
                    i++;
                }
                int stopdex = i-1;
                int len = stopdex-startdex+1;
                double msc = sc/len;
                double sdsc = Math.sqrt(scsc/len - msc*msc);
                if (len >= minlen) {
                    startaa[numseg] = startdex;
                    stopaa[numseg] = stopdex;
                    lenaa[numseg] = len;
                    localmean[numseg] = msc;
                    localsd[numseg] = sdsc;
                    localhydro[numseg] = lh/len;
                    localcharge[numseg] = lc/len;
                    numseg++;
                }
                // now trim arrays to length numseg?
            }
            else i++;
        }
        if (lenaa.length>0) maxlen = spewey.maxint(lenaa);
        maxlong = 1;
        bestdex = 0;
        for (i=0; i<numseg; i++) {if (lenaa[i] >= minlen2 && localmean[i]<maxlong) {maxlong = localmean[i]; bestdex = i;};}
    }


    public void printme() {
        printme("");
    }


    // feed this to R.
    public void printme(String id){
        System.out.println("### protein length:	   " + n);
        System.out.println("### num disordered(1): " + numdisordered);
        System.out.println("### num disordered(2): " + numdisorderedstrict);
        System.out.println("### num disordered(3): " + numdisorderedstrict2);
        for (int i=0; i<numseg; i++)  {
            System.out.println("### " + startaa[i] + "-" + stopaa[i] + ": " + (float) localmean[i] + " +- " + (float) localsd[i]);
        }
        System.out.println("### " + maxlen + " " + (float) localmean[bestdex] + " " + lenaa[bestdex]);
        System.out.println("# " + n + "\t" + (float) meancharge + "\t" + (float) meanhydro + "\t" + (float) meancombo + "\t" + (float) maxlong);
        System.out.println("# " + (float) hssr[0] + "\t" + (float) hssr[1] + "\t" + (float) hssr[2] + "\t"
                + (float) hssr2[0] + "\t" + (float) hssr2[1] + "\t" + (float) hssr2[2]);
        for (int i=0; i<n; i++) {
            System.out.println(id + "\t" + aa[i] + "\t" + (float) charge[i] + "	 \t"
                    + (float) hydro[i] + "  \t" + (float) combo[i] + "\t" + (float) other[i] + " \t"	 + (float) other2[i]
                    + "\t # ["+(i+1-ww/2)+"-"+(i+1 + ww/2)+"]");
        }
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

