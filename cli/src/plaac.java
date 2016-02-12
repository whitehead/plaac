/** /////////////////////////////////////////////////////////////////////////////
    Original program copyright 2009 Whitehead Institute for Biomedial Research;
    additions copyright 2011 BBRI and copyright 2014 University of Massachusetts Medical School.
    Author: Oliver King (oliver.king@umassmed.edu) 
    
    See LICENSE.TXT for license information.  

    Last updated May 12, 2014
    
    Compile with: javac plaac.java
    To see usage details, run with: java plaac 
   
**/ ////////////////////////////////////////////////////////////////////////

// TODO:
// Compute fg params from file (need pseudocount option and weighting option)
// Check different uses of eps

import java.io.*;
import java.util.*;
import java.net.*;
import java.util.regex.*;

class plaac {

    static char [] aanames =   {'X','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'};
    static String [] aanames2 = {"???","Ala","Cys","Asp","Glu","Phe","Gly","His","Ile","Lys","Leu","Met","Asn",
				 "Pro","Gln","Arg","Ser","Thr","Val","Trp","Tyr","***"};

    static double log2 = Math.log(2);

    // look-up table for logarithms when using logs in HMM
    static double [] loglut;
    static int loglutlength = 40*100;
    
    // charge of each AA residue, used in FoldIndex
    static double [] aacharge = {
	0, // 'X'
	0, // 'A'
	0, // 'C'
	1, // 'D'
	1, // 'E'
	0, // 'F'
	0, // 'G'
	0, // 'H' // allow pH dependency?
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


    // Kyte-Doolittle Hydropathicity: J Mol Biol 157:105-132 1982
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

    // scaled and shifted version of above, used in FoldIndex
    static double [] aahydro2 = axpb(1.0/9.0,aahydro,0.5);

    // From Toombs, McCarty and Ross, MCB 2010
    static double [] fgpapa1 = {
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

    // From Toombs, McCarty and Ross, MCB 2010
    static double [] bgpapa1 = {
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


    // From Toombs, McCarty and Ross, MCB 2010
    static double [] fgpapa2 = {
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

    // From Toombs, McCarty and Ross, MCB 2010
    // D 4% too high, E too low?  F 3% too high, R too low?
    // Note: these published frequencies don't agree with the raw AA counts and don't add up to 1; 
    // Can recompute frequencies from raw AA counts with recompute_papa_parameters()
    static double [] bgpapa2 = {
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
   
    // log of odds-ratios
    static double [] lodpapa1;
    static double [] lodpapa2;
   

    // odds-ratios
    static double [] odpapa1 ={
	0.0,
	0.67267686, // A
	1.5146198,  // C
	0.27887323, // D
	0.5460614,  // E
	2.313433,   // F
	0.96153843, // G
	0.75686276, // H
	2.2562358,  // I
	0.20664589, // K
	0.9607843,  // L
	1.9615384,  // M
	1.0836071,  // N
	0.30196398, // P
	1.0716166,  // Q
	0.6664044,  // R
	1.1432927,  // S
	0.8917492,  // T
	2.2562358,  // V
	1.9478673,  // W
	2.1785367,  // Y
	0.0         // *
    };

    // odds-rations
    // inconsistency causd by bgpapa2?
    static double [] odpapa2 ={
	0.0,        // X
	0.88066554, // A
	1.2461538,  // C
	0.72039115, // D ## 0.66?? [45 vs 67]
	0.77220076, // E ## 0.93?? [23 vs 24]
	3.6936572,  // F ## 3.16?? [64 vs 21] ct 17/264 vs 6/328
	1.2570281,  // G
	1.1460011,  // H
	1.2519685,  // I
	0.17436177, // K
	0.87114847, // L
	1.4330357,  // M
	1.7729831,  // N
	0.7429888,  // P
	0.8229167,  // Q
	0.5531136,  // R  ## 0.58?? [45 vs 76]
	0.9144572,  // S
	0.9143354,  // T
	1.0367454,  // V
	4.710145,   // W
	0.75716186, // Y
	0.0,        // *
    };

    // static double [] gnq = {0,0,0,0,0,0,1.0,0,0,0,0,0,1.0,0,1.0,0,0,0,0,0,0,0};
   
    // global AA frequencies from S. cerevisiae used in Alberti et al Cell 2009. 
    static double [] bg_freq_scer = {0,0.0550,0.0126,0.0586,0.0655,0.0441,0.0498,0.0217,0.0655,0.0735,0.0950,0.0207,
			     0.0615,0.0438,0.0396,0.0444,0.0899,0.0592,0.0556,0.0104,0.0337,0};
   
    // values used in Alberti et al Cell 2009 based on prion domains from Sup3p5, Rnq1p, Ure2p, New1p
    static double [] prd_freq_scer_04 = {0,0.0488,0.0032,0.0202,0.0234,0.0276,0.1157,0.0149,0.0191,0.0329,0.0456,
				  0.0149,0.1444,0.0308,0.2208, 0.0202,0.1008,0.0297,0.0234,0.0064,0.0573,0};
    
    // values derived from 28 experimentally determined prion-like domains in Alberti et al Cell 2009
    static double [] prd_freq_scer_28 = {0,0.04865,0.00219,0.01638,0.00783,0.02537,0.07603,0.0181,0.02018,0.01641,0.02639,
				 0.02975,0.25885,0.05126,0.15178,0.025,0.10988,0.03841,0.01972,0.00157,0.05624,0};
    
    // uniform bg frequencies; not currently used. 
    // static double [] unibg = {0.0,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.0};
    
    static double [][] aastats;

    static boolean verbose = false;

    plaac() {
	
        // used by HMM when doing using log values rather than scaling to prevent underflow
	loglut = new double[loglutlength+1]; 
	for (int i=0; i<=loglutlength; i++) loglut[i] = Math.log(1.0 + Math.exp(-i/100.0));

       	lodpapa1 = new double[22];
	lodpapa2 = new double[22];

	for (int k=1; k<=20; k++) {
	    lodpapa1[k] = Math.log(odpapa1[k]);
	    lodpapa2[k] = Math.log(odpapa2[k]);
	}

       	// alternative: compute lods directly from frequencies in sequence files; 
	// should agree for papa1, not for papa2, but these aren't currently used.
        
	// recompute_papa_parameters();
        

    }


    public static void main (String args[]) {

	plaac ms = new plaac(); 

        // if (1==1) {print_aa_params(newfgred); return;}
	// if (1==1) {print_aa_params(oldfgfreq); return;}

        // store defaults in textfile to be read with -F?
	double [] fgfreq = prd_freq_scer_28;

	double [] bgscer = normalize(bg_freq_scer);     // overall AA frequencies in S.cer
	String inputfile = "";
	String bgfile = "";
       	String bgfreqfile = "";
        String fgfreqfile = "";
	String plotlist = "";
        String hmmdotfile= "";
	int corelength = 60; // for PLACC HMM  
	int ww1 = 41; // window size for sliding averages of FoldIndex
	int ww2 = 41; // window size for sliding averages of PAPA
	int ww3 = 41; // window size for sliding averages of PLAAC LLR
	double alpha = 1.0;
	int hmmtype = 1;
	boolean printheaders = false;   // print the headers at the top of the file (default: false)
       	boolean printparameters = true; // print the parameters at the top of the file (default: true)

        // from PAPA: score only first proline in PP or PXP
       	boolean adjustprolines = true;


        // should check for valid input here!
        int i = 0;

	// to avoid out-of-bounds attempt to grab an option from position i+1 for invalid command line-arguments, 
	// iteration stops at next-to-last position, unless last position is one of the flags that doesn't take options.
        while(i<args.length - 1 || (i<args.length && (args[i].equals("-d") || args[i].equals("-s")))) {
	    if (args[i].equals("-i")) {inputfile = args[i+1]; i++;}
	    else if (args[i].equals("-b")) {bgfile = args[i+1]; i++;}
	    else if (args[i].equals("-B")) {bgfreqfile = args[i+1]; i++;}
            else if (args[i].equals("-F")) {fgfreqfile = args[i+1]; i++;}
	    else if (args[i].equals("-c")) {corelength = Integer.parseInt(args[i+1]); i++;}
	    else if (args[i].equals("-w")) {ww1 = Integer.parseInt(args[i+1]); i++;}
	    else if (args[i].equals("-W")) {ww2 = Integer.parseInt(args[i+1]); i++;}
	    else if (args[i].equals("-a")) {alpha = Double.parseDouble(args[i+1]); i++;}
	    else if (args[i].equals("-m")) {hmmtype = Integer.parseInt(args[i+1]); i++;} // different models. not currently used
	    else if (args[i].equals("-p")) {plotlist = args[i+1]; i++;}
            else if (args[i].equals("-h")) {hmmdotfile = args[i+1]; i++;} 
            else if (args[i].equals("-d")) {printheaders = true; } // don't consume another token, '-d' doesn't take an argument 
            else if (args[i].equals("-s")) {printparameters = false; } // don't consume another token, '-s' doesn't take an argument 
	    else {System.out.println("# skipping unknown option " + args[i]);}
            i++;
	}

	ww3 = ww2; // don't currenty have a separate command-line paramater for ww3

	if (verbose) {
	    System.out.println("## -i -->"+inputfile+"<---");
	    System.out.println("## -b -->"+bgfile+"<---");
	    System.out.println("## -B -->"+bgfreqfile+"<---");
            System.out.println("## -F -->"+fgfreqfile+"<---");
	    System.out.println("## -m -->"+hmmtype+"<---");
	    System.out.println("## -p -->"+plotlist+"<---");
	    System.out.println("## -c -->"+corelength+"<---");
	    System.out.println("## -w -->"+ww1+"<---");
	    System.out.println("## -W -->"+ww2+"<---"); 
	    System.out.println("## -a -->"+alpha+"<---");
            System.out.println("## -h -->"+hmmdotfile+"<---");
            System.out.println("## -d -->"+printheaders+"<---");
            System.out.println("## -s -->"+(!printparameters)+"<---"); // negated, since s means skip printing
	} 

	// compute global AA frequencies for sequences being scored   
        double [] bgf = new double[22];
	double [] bgthis = new double[22];

	if (bgfreqfile.length() > 0) {
            // catch exceptions?
	    bgf = read_aa_params(bgfreqfile);
	} else if (bgfile.length() > 0) {
	    bgf = asdouble(computeaafreq(bgfile,false));
	} else if (inputfile.length() > 0) {
	    bgf = asdouble(computeaafreq(inputfile,false));
	} 

        if (fgfreqfile.length() > 0) {
            // catch exceptions?
	    fgfreq = read_aa_params(bgfreqfile);
	    // FIXME: may want to add pseudocounts here before normalizing --- some AAs may not be seen at all in training samples. 
	    // should check if the read params are counts or frequencies.
	}
	
	// can be redirected to file to use later with -B option
       	if ((bgfile.length() > 0) && (inputfile.length() == 0) ) {
	    print_aa_params(bgf);
	    return;
	}
	
	// checks that bgfreqfile is being read correctly
	if ((bgfreqfile.length() > 0) && (inputfile.length() == 0)) {
	    print_aa_params(bgf);
	    return;
	}
	
	if ((inputfile.length() == 0) && (bgfile.length() == 0)) {
	    System.out.println("------------------------------------------------------------");
	    System.out.println("This program is offered with NO WARRANTY WHATSOEVER; see LICENSE.TXT for details.");
	    System.out.println("------------------------------------------------------------");
	    System.out.println("USAGE: Need to specify an input protein fasta file after the flag -i, e.g.\n  java -jar plaac.jar -i input.fa > output.txt");
	    System.out.println("Note that the example above redirects the output of the program to a tab-delimited text file output.txt");
	    System.out.println(" that can be opened with a spreadsheet program.");
	    System.out.println("Optional arguments:");
	    System.out.println("-c core_length, where the integer core_length is the minimal contiguous prion-like domain length\n  for the HMM parses. Default is 60.");
	    System.out.println("-B bg_freqs.txt, specifying background AA freqs to use for the species, one per line, in the following order:");
            System.out.println("  X, A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y, *");
	    System.out.println(" (Values for X and * will be set to zero and other numbers normalized to add to 1)");
	    // X and * get zeroed out -- could exclude from file spec. Add pseudocounts before normalization?
	    System.out.println("-b background.fa, where background.fa is the name of a protein fasta file used to\n  compute background AA frequencies for the species.");
            System.out.println("  This option is ignored if -B is used, but otherwise if -b is not specified it defaults to the input.fa file.");
	    System.out.println("  If the sequences in input.fa have biased AA composition then a separate background.fa or bg_freqs.txt is recommended.");
            System.out.println("  If -b is specified but -i is not, AA counts for background.fa will be written to standard output, and the program will exit.");
	    System.out.println("  These counts can be redirected to a file (e.g. with > bg_freqs.txt), in a format that can be read by the -B option.");
	    System.out.println("-a alpha, where alpha is a number between 0 and 1 (inclusive) that controls the degree to which the S. cerevisiae") ;
	    System.out.println("  background AA frequencies are mixed with the background AA frequencies from -B, -b, or -i.");
	    System.out.println("  If alpha = 0, just the AA frequencies from the -B, -b, or -i are used, and if alpha = 1 just the\n  S. cerevisiae AA frequencies are used. Default is 1.0.");
            System.out.println("-F fg_freqs.txt, specifying prion-like AA freqs in same format as -B above. Defaults to freqs from 28 S. cerevisiae domains.");
	    System.out.println("-w window_size, the window size for FoldIndex disorder predictions. Default is 41.");
	    System.out.println("-W Window_size, the window size for the PAPA algorithm. Default is 41.");
	    System.out.println("-d, print documentation for headers. If flag is not set, headers will not be printed.");
	    System.out.println("-s, skip printing of run-time parameters at top of file. If flag is not set, run-time parameters will be printed.");
            System.out.println("-h hmm_filename.txt, writes parameters of HMM to hmm_filenmae.txt in dot format, which can be made into a figure with GraphViz.");
	    System.out.println("-p print_list.txt, where print_list.txt has the name of one fasta on each line, and specifies\n  which fastas in input.fa will be plotted");
	    System.out.println("  The names must exactly match those in input.fa, but do need need the > symbol before the name.");
	    System.out.println("  If no print_list.txt is specified the output from the program will be a table of summaries for each protein (one per line) in input.fa;");
	    System.out.println("  If a print_list.txt is specified the output from the program will be a table (one line per residue) that is used");
	    System.out.println("  for making plots for each of the proteins listed in print_list.txt.");
	    System.out.println("  If the option is given as -p all, then plots will be made for all of the proteins in input.fa, \n  which is not advised if input.fa is an entire proteome.");
	    System.out.println("  To make the plots from output that has been redirected to output.txt, at the command-line type type\n  Rscript plaac_plot.r output.txt plotname.pdf.");
	    System.out.println("  This requires that the program R be installed (see http://www.r-project.org/)\n  and will create a file named plotname.pdf, with one plot per page.");
	    System.out.println("  Calling Rscript plaac_plot.r with no file specified will list other options for plotting."); // FIXME: need to add this to plaacplot!
	    return;
	}

	if (alpha > 1 || alpha < 0) {
	    System.out.println("# warning: invalid alpha; using alpha = 1.0");
	    alpha = 1.0;
	} 

        fgfreq[0] = 0; fgfreq[21] = 0; // ignore X and * counts if freq reads from file
	// fgfreq = axpb(1.0, normalize(fgfreq), 0.0001); 
        // fgfreq[0] = 0; fgfreq[21] = 0; 
       	fgfreq = normalize(fgfreq); 

	bgf[0] = 0; bgf[21] = 0; 
        // shouldn't need pseudocounts for bgf, if they are based on whole proteome
	bgthis = normalize(bgf); // empty input file: all zero ends up all zero?
	// use weighted combination of bgscer and bgthis as bg freqs: alpha*(s.cer) + (1-alpha)*(this org)
	double [] bgcombo = normalize(axpby(alpha, bgscer,(1-alpha), bgthis));

	// read table of parameters from file: could do this for all fg and bg AA freqs, etc.
	// double[][] mymat = readmatrix2("util/aa-biases-core.txt",1,1);
	// aastats = transpose(mymat);
	// double [] ogf = aastats[1]; // first row of matrix
	// ogf[0] = 0; ogf[21] = 0; 
	// ogf = normalize(axpb(1.0,ogf,0.0001));
        // printmatrix(mymat);
	// printmatrix(mymat[0]);
        // printaausage(aastats[0]);
	
	if (verbose) {
            System.out.println("## fg-freq");
	    printaausage(bgscer); 
	    System.out.println("## bg-scer");
	    printaausage(bgscer); 
	    System.out.println("## bg-this");
	    printaausage(bgthis); 
	    System.out.println("## bg-mixed");
	    printaausage(bgcombo);
	    System.out.println("## papa-lod");
	    printaausage(lodpapa1); 
	    System.out.println("## aa-hydro");
	    printaausage(aahydro);
	    System.out.println("## aa-hydro2");
	    printaausage(aahydro2);     
	    // System.out.println("## og");
	    // printaausage(ogf); 
	}  

	// include small prob so HMMs won't give -Infinity if there is X or * in protein sequence
        double epsx = 0.00001;
        
	fgfreq[0]  = epsx;
        fgfreq[21] = epsx;
        bgcombo[0]  = epsx;
        bgcombo[21] = epsx;
       	double [] fg = normalize(fgfreq); 
	double [] bg = normalize(bgcombo); 
	double [] llr = new double[22];
	// due to epsx above, LLR would automatically be zero for j=0 (X) and j=21 (*) anyway
	for (int j=1; j<21; j++) llr[j] = Math.log(fg[j]/bg[j]);


        if (printparameters) {
	    System.out.println("############################ parameters at run-time ####################################");
            System.out.println("## alpha="+ alpha +"; corelength="+corelength+"; ww1="+ww1+"; ww2="+ww2+"; ww3="+ww3+"; adjustprolines="+adjustprolines+";");
            // could make these AA params tab-separated so the parameters line up, one AA per column
            System.out.println("## fg_used: {" + aaparams2string(fg) + "}");
            System.out.println("## bg_scer: {" + aaparams2string(bgscer) + "}");
            System.out.println("## bg_input: {" + aaparams2string(bgthis) + "}");
            System.out.println("## bg_used: {" + aaparams2string(bg) + "}");
            System.out.println("## plaac_llr: {" + aaparams2string(llr) + "}");
            System.out.println("## papa_lods: {" + aaparams2string(lodpapa1) + "}");
	    System.out.println("#######################################################################################");
        } 


        hmm hmm1 = prionhmm1(fg,bg);
        hmm hmm0 = prionhmm0(bg);

        if (hmmdotfile.length() > 0) {
	    hmm1.dottify(hmmdotfile, true); 
	}

        if ((inputfile.length() > 0) && (plotlist.length() == 0)) {
	    // can pass in any fg and bg frequencies here, if these are not suitable
	    scoreallfastas(inputfile, corelength, ww1, ww2, ww3, fg, bg, llr, hmm1, hmm0, adjustprolines, printheaders);
	} else if ((inputfile.length() > 0) && (plotlist.length() > 0)) {
	    plotsomefastas(inputfile, corelength, ww1, ww2, ww3, fg, bg, llr, hmm1, hmm0, plotlist, adjustprolines);
	}
    }


   // The frequencies and odds for papa2 scores from Toombs et al MCB 2010 don't appear to be consistent. 
   // Below we compute frequencies from the listed positive and negative sequences.
   // These overwrite hard-coded defaults. 

   static int recompute_papa_parameters() {
	int [] aac; 
	double [] aaf;

	aac = computeaafreq(new String("util/rosspsi1.txt"),false);
	aac[0] = 0; aac[21] = 0;
	aaf = normalize(aac);
	fgpapa1 = aaf;

	aac = computeaafreq(new String("util/rossnai1.txt"),false);
	aac[0] = 0; aac[21] = 0; 
	aaf = normalize(aac);
	bgpapa1 = aaf;

	aac = computeaafreq(new String("util/rosspsi2.txt"),false);
	aac[0] = 0; aac[21] = 0; 
	aaf = normalize(aac);
	fgpapa2 = aaf;

	aac = computeaafreq(new String("util/rossnai2.txt"),false);
	aac[0] = 0; aac[21] = 0; 
	aaf = normalize(aac);
	bgpapa2 = aaf;

	lodpapa1 = new double[22];
	lodpapa2 = new double[22];
	odpapa1 = new double[22];
	odpapa2 = new double[22];

	for (int k=1; k<=20; k++) {
	    odpapa1[k]= ((fgpapa1[k]/(1-fgpapa1[k]))/(bgpapa1[k]/(1-bgpapa1[k])));
	    odpapa2[k]= ((fgpapa2[k]/(1-fgpapa2[k]))/(bgpapa2[k]/(1-bgpapa2[k])));
	    lodpapa1[k]= Math.log((fgpapa1[k]/(1-fgpapa1[k]))/(bgpapa1[k]/(1-bgpapa1[k])));
	    lodpapa2[k]= Math.log((fgpapa2[k]/(1-fgpapa2[k]))/(bgpapa2[k]/(1-bgpapa2[k])));
	}

	if (verbose) {
	    System.out.println("## od-papa1");
	    printaausage(odpapa1);
	    System.out.println("## lod-papa1");
	    printaausage(lodpapa1);
	    System.out.println("## od-papa2");
	    printaausage(odpapa2);
	}

	return 1;
   }


    // Pass PAPA lods? (taken from global variable now). hmm0 not currently used in this function.
    static void plotsomefastas(String filename,	 int corelength, int ww1, int ww2, int ww3, double [] fg, 
			       double [] bg, double [] llr, hmm hmm1, hmm hmm0, String genelist, boolean adjustprolines) {

	HashMap<String,String> ht1 = new HashMap<String,String>(); // for synonyms used as plot titles
	HashMap<String,String> ht2 = new HashMap<String,String>(); // for ordering plots

	boolean plotall = false; // plot all genes in fasta file
	if (genelist.equals("all")) {
	    plotall = true;
	} else {
	    ht1 = readhashtable(genelist, null); 
	    ht2 = readhashtable(genelist,"incr");
	}

	fastareader fs = new fastareader(filename);
	// changed GENE to SEQID
	System.out.print("ORDER\tSEQid\tAANUM\tAA\tVIT\tMAP\tCHARGE\tHYDRO\tFI\tPLAAC\tPAPA\tFIx2\tPLAACx2\tPAPAx2");
	for (int i=0; i<hmm1.numclasses; i++) System.out.print("\tHMM." + hmm1.unames[i]);
	System.out.println();

	int genecount = 1;

	char stopcodon = '*';
	while (fs.hasmorefastas()) {
	    String name = fs.nextname();
	    String nm = name;
	    String nm2 = ">" + nm;
	    StringBuffer sb = fs.nextfasta();

	    // checks if name is on plot list. Should probably keep track of which genes in plot list were not found in fasta.
	    if (plotall || ht1.containsKey(nm) || ht1.containsKey(nm2)) {

		String id = new String("" + genecount);

		if  (ht1.containsKey(name)) nm = (String) (ht1.get(name));
		if  (ht2.containsKey(name)) id = (String) (ht2.get(name));

		if ((sb.charAt(sb.length()-1)) == stopcodon) sb.deleteCharAt(sb.length()-1); // kill terminal stop codon
		int [] seq = string2aa(sb);
		genecount++;

		hmm1.decodealls(seq);
		// bghmm.decodealls(seq);
		// myhmm.printme(genecount + "\t" + name, seq);
		// myhmm.sposteriorl(seq);
		// myhmm.sviterbidecodel(seq);
		disorderreport dr = new disorderreport(seq, ww1, ww2, ww3, new double [] {2.785 , -1  ,-1.151}, llr, lodpapa1, adjustprolines);
		//  dr.printme(genecount + "\t" + name);
		for (int i=0; i<seq.length; i++) {
                    // print i+1 rather than i, for one-based indices
		    System.out.print(id + "\t" + nm + "\t" + (i+1) + "\t" + aanames[seq[i]] + "\t" + hmm1.viterbipath[i] + "\t" + hmm1.mappath[i]+"\t");
		    System.out.format("%.4f\t%.4f\t%.8f\t%.4f\t%.8f\t%.8f\t%.4f\t%.8f", dr.charge[i], dr.hydro[i], dr.fi[i], 
				      dr.plaacllr[i], dr.papa[i], dr.fix2[i], dr.plaacllrx2[i], dr.papax2[i]);
                    // use exp format here? not a big deal for plotting, but if we want to take logs later could be helpful
		    for (int j=0; j<hmm1.postprob.length; j++) System.out.format("\t%.4f", hmm1.postprob[j][i]);
		    System.out.println();
		}
		// System.out.format("#### %.4f\t%.4f\t%.4f\t%.4f\n", hmm1.lmarginalprob, hmm1.lviterbiprob,  hmm1.lviterbiprob,  hmm1.lviterbiprob);
		System.out.println("########################################################");

	    }
	}
    }


    // Pass PAPA lods? (taken from global variable now)
    static void scoreallfastas(String filename,  int corelength, int ww1, int ww2, int ww3, double [] fg, 
			       double [] bg, double [] llr, hmm hmm1, hmm hmm0, boolean adjustprolines, boolean printheaders) {

	// Random rg = new Random(7);
	int chomp = 0;

        if (printheaders) {

	    System.out.println("############################ Description of output columns ############################");	
	    
	    System.out.println("## SEQid: sequence name from fasta file");
	    System.out.println("## MW: Michelitsh-Weissman [PNAS 2000] score --- maximum number of N + Q in window of at most 80 AA");
	    System.out.println("## MWstart: index of start position of window with max MW score.");
	    System.out.println("## MWend: index of end position of window with max MW score");
            System.out.println("## MWlen: length of window used for MW score [smaller of 80 and PROTlen]."); // Redundant with start and end
           
	    System.out.println("## LLR: max sum of PLAAC log-likelihood ratios (base 4) in window of size c [NaN if PROTlen < c]");
            System.out.println("## LLRstart: index of start position of window with max LLR score [-1 if PROTlen < c]");
	    System.out.println("## LLRend: index of end position of window with max LLR score [-1 if PROTlen < c]"); 
           // Redundant with start and c, unless we allow small c for PROlen < c
            System.out.println("## LLRlen: length of window used for LLR score [should be c]");  // Redundant with start and end
	    System.out.println("## NLLR: normalized LLR score, i.e. LLR/LLRlen [NaN if PROTlen < c]"); // Redundant with LLR and LLRlen
            
	    System.out.println("## VITmaxrun: maximum length of consecutive PrD state in Viterbi parse");
	    System.out.println("## COREscore: max sum of PLAAC LLRs in window of size c contained entirely within Viterbi parse [NaN if VITmaxrun < c]");
            System.out.println("## COREstart: index of start position of window with max COREscore [-1 if VITmaxrun < c]");
	    System.out.println("## COREend: index of end position of window with max COREscore [-2 if VITmaxrun < c]"); // Redundant with start and c
            System.out.println("## CORElen: length of window used for CORElen [should be either c or 0]");  // Redundant with start and end
	    // moved PRDScore column to here
            System.out.println("## PRDscore: sum of PLAAC LLRs in full region of Viterbi parse containing CORE region, if any [NaN otherwise].");
            System.out.println("## PRDstart: index of start position of window with PRDscore [-1 if VITmaxrun < c]");
	    System.out.println("## PRDend: index of end position of window with PRDscore [-2 if VITmaxrun < c]");
            System.out.println("## PRDlen: length of window used for PRDscore");  // Redundant with start and end
            System.out.println("## PROTlen: number of AAs in protein, not including terminal stop codon if any.");
	    System.out.println("## HMMall: log-likelihood ratio for sequence under two-state HMM vs. one-state background HMM"); // base e? renamed from HMMmap.
            System.out.println("## HMMvit: log-likelihood ratio for sequence under Viterbi parse of two-state HMM vs. one-state background HMM"); // base e?
	    System.out.println("## COREaa: AA sequence at which COREscore is attained [- if VITmaxrun < c]");
            System.out.println("## STARTaa: first 15 AA of PRDaa [- if VITmaxrun < c]"); 
	    // Redundant with PRDaa --- mostly used to check that the parse starts with a sensible residue, but a few upstream AAs would help for that (if there are any)
            System.out.println("## ENDaa: last 15 AA of PRDaa [- if VITmaxrun < c]"); 
	    // Redundant with PRDaa mostly used to check that the parse starts with a sensible residue, but a few downstream AAs would help for that (if there are any)
	    System.out.println("## PRDaa: AA sequence at which PRDscore is attained [- if VITmaxrun < c]");
            System.out.println("## FInumaa: number of AAs predicted to be disordered by FoldIndex [Prilusky et al, Bioinformatics, 2005] (exludes runs of under 5 AA)"); 
            System.out.println("## FImeanhydro: hydropathy score <H> for entire protein [Uversky et al, Proteins, 2000]");
	    System.out.println("## FImeancharge: mean charge <R> for entire protein [Uversky et al, Proteins, 2000]"); 
            System.out.println("## FImeancombo: disorder score for entire protein: 2.785<H> - |<R>| - 1.151 [Uversky et al, Proteins, 2000]");  
            System.out.println("## FImaxrun: length of longest run of predicted disorder by FoldIndex");
	   
	    // need to explain these better!!
            System.out.println("## PAPAcombo: signed distance to PAPA decision surface [as in King et al Brain Res 2012]"); 
	    // at PAPAcen? Or anywhere? can use this even no region has negative FI score. Need to fix this.
            System.out.println("## PAPAprop: maximal score of PAPA prion propensities (averges of averages) in region with negative FI score [Toombs et al MBC 2012]"); 
	    System.out.println("## PAPAfi: FI score (averages of averages) at PAPAcen");
            System.out.println("## PAPAllr: PLAAC LLR score (average) at PAPAcen");
            System.out.println("## PAPAllr2: PLAAC LLR score (average of averages) at PAPAcen");
	    System.out.println("## PAPAcen: index of center of window at which PAPAprop is obtained"); // move to earlier?
	    System.out.println("## PAPAaa: AA sequence of width W centered at PAPAcen");
	    
	    System.out.println("#######################################################################################");

	}

	System.out.print("SEQid\tMW\tMWstart\tMWend\tMWlen\tLLR\tLLRstart\tLLRend\t");	
	System.out.print("LLRlen\tNLLR\tVITmaxrun\tCOREscore\tCOREstart\tCOREend\tCORElen\tPRDscore\tPRDstart\tPRDend\tPRDlen\tPROTlen\t");				
	System.out.print("HMMall\tHMMvit\tCOREaa\tSTARTaa\tENDaa\tPRDaa\tFInumaa\tFImeanhydro\tFImeancharge\tFImeancombo\tFImaxrun\t");
	System.out.print("PAPAcombo\tPAPAprop\tPAPAfi\tPAPAllr\tPAPAllr2\tPAPAcen\tPAPAaa");
	System.out.println();


	fastareader fs = new fastareader(filename);
	// int maxlen = 500; // not currently used
	// hgalg hg = new hgalg(maxlen,bg);
	int genecount = 0;
        
        //                 {X,A,C,D,E,F,G,H,I,K,L,M,N  ,P,Q,  R,S,T,V,W,Y,*}
	double [] qnmask = {0,0,0,0,0,0,0,0,0,0,0,0,1.0,0,1.0,0,0,0,0,0,0,0};

	double [][] lops;
	double [] maa1;
	double [] maa2;
	double [] maa3;
	double [] hs1;
	double [] hs2;
	double [] hs3;
	
	// test hss vs hss2
	double [] hs1b;
	double [] hs2b;
	double [] hs3b;

	int [] aa;
	String flag;

	double hmmscore;
	double hmmscorev;

	int longestprd;

	char stopcodon = '*';

        boolean debug = false;

	while (fs.hasmorefastas()) {
	    String name = fs.nextname();
	    StringBuffer sb = fs.nextfasta();
	    if ((sb.charAt(sb.length()-1)) == stopcodon) sb.deleteCharAt(sb.length()-1); // kill stop codon
	    aa = string2aa(sb);

	    // something longer here? w=41? c=60? 
	    if (aa.length < 1) continue; // print output row for filtered sequences instead?

	    // MW: Michelitsh-Weissman [PNAS 2000] score
	    // maximum number of N + Q in window of at most 80 AA can allow window to be smaller than 80,
            //  since we're looking for max sum rather than max averge, and sum is monotonic
	    maa1 = mapseq(aa, qnmask);

            int mwsize = 80;
            if (aa.length < 80) mwsize = aa.length;
	    hs1 = hss2(maa1, mwsize, mwsize); //CHANGED: mwsize for short proteins changed here, rather than inside hss2
	    if (debug) {
		hs1b = hss(maa1, mwsize, mwsize); //CHANGED: mwsize for short proteins changed here, rather than inside hss2
		if (!(Double.isNaN(hs1[2]) && Double.isNaN(hs1b[2])) && (hs1[2] != hs1b[2])) {
		    System.out.println("## CHECKX1");
		    System.out.print("## hs1: "); printrowvec(hs1); 
		    System.out.print("## hs1b: "); printrowvec(hs1);
		} 
	    }
	   
	    // LLR
	    maa2 = mapseq(aa, llr);
	    hs2 = hss2(maa2, corelength, corelength); // CHANGED: min length no longer shortened when protlen < corelength
	    if (debug) {
		hs2b = hss(maa2, corelength, corelength);   // CHANGED: min length no longer shortened when protlen < corelength
		if (!(Double.isNaN(hs2[2]) && Double.isNaN(hs2b[2])) && (hs2[2] != hs2b[2])) {
		    System.out.println("## CHECKX2"); 
		    System.out.print("## hs2: "); printrowvec(hs2); 
		    System.out.print("## hsb: "); printrowvec(hs2b);
		} 
	    }

	    // HMM
	    hmm1.decodeall(aa);
	    hmm0.decodeall(aa);

	    hmmscore  = hmm1.lmarginalprob - hmm0.lmarginalprob;
	    hmmscorev = hmm1.lviterbiprob  - hmm0.lviterbiprob;

	    disorderreport  dr = new disorderreport(aa, ww1, ww2, ww3, new double [] {2.785 , -1, -1.151}, llr, lodpapa1, adjustprolines);
	    // dr.printme();
            // System.out.println(dr.papamaxprop + "\t" + dr.papamaxscore);
	   
	    // gets rid of leading ">" from fasta 
	    String nm = name.substring(chomp);

	    // int chopdex = nm.indexOf(" ",	nm.indexOf(" ")+1);
	    // nm = nm.substring(0, chopdex);

	    int [] seg;
	    int [] ppp;

             // could use MAP parse here if desired
	    int [] mp = hmm1.viterbipath;

	    longestprd = longestrun(mp);

	    maa3 = mapseq(aa, llr);
	    double big_neg = -1000000.0;

	    // A bit of a hack here. The idea is to mask off positions that aren't in the PrD parse when finding the hss, 
	    // by giving them a sufficiently strong penalty. The implementation of hss doesn't work with -Infinty since
	    // it uses differences in cumulative sums and the diff between -Inf and -Inf isn't well defined. But masking 
	    // elements by setting them to a large negative number, say -K where K > 2 * max(abs(maa3))*max(corelength) 
	    // should work, since then any segment with a masked residue would have score < -K/2; if the max scoring segment 
	    // has such a score then there is no PrD segment of width >= corelength; 
	    // For fixed rather that variable width windows, an simple implementation that does not use cumulative somes
	    // would do avoid this issue.
	    for (int i=0; i < maa3.length; i++) {
		if (mp[i]==0) maa3[i] = big_neg; 
	    }

	    hs3 = hss2(maa3,corelength,corelength); //CHANGED: no longer relax corelength when corelength > prot length 
	    if (debug) {
		hs3b = hss(maa3,corelength,corelength);
		if (!(Double.isNaN(hs3[2]) && Double.isNaN(hs3b[2])) && (hs3[2] != hs3b[2])) {
		    System.out.println("## CHECKX3"); 
		    System.out.print("## hs3: ");  printrowvec(hs3);
		    System.out.print("## hs3b:");  printrowvec(hs3b);
		} 
		if ((longestprd >= corelength) && (hs3[2] < big_neg/2)) {
		    System.out.println("## CHECKX4 " + longestprd); 
		    System.out.print("## hs3: "); printrowvec(hs3);
		} 
		if ((longestprd < corelength) && (hs3[2] > big_neg/2)) {
		    System.out.println("## CHECKX5 "  + longestprd); 
		    System.out.print("## hs3: "); printrowvec(hs3);
		}
	    }
	   
	    int corestart = (int) hs3[0];
	    int corestop = (int) hs3[1];
	     
	    int aastart = corestart;
	    int aastop	= corestop;

	    int [] prd = {};
	    double prdscore = 0;
	 
	    // has prd of at least corelength
	    if (hs3[2] > big_neg/2) { 
		// expand core up and down within PrD parse;
		while (aastart>=0 && mp[aastart]==1) aastart--; 
		aastart++;
		while (aastop<mp.length && mp[aastop]==1) aastop++;
		aastop--;
	
		prd = submatrix(aa, aastart, aastop);
		// score of whole prd?
		for (int kk = 0; kk<prd.length; kk++) {
		    prdscore = prdscore + llr[prd[kk]];
		}
	    }  else {
		//  if (aastop - aastart + 1 < corelength) {
		hs3[2] = 0.0/0.0; // CHANGED: NaN; used to be zero.
		aastart = -1;
		aastop  = -2;
		corestart = -1;
		corestop = -2;
	    }

            if (debug) {
		if (aastop - aastart + 1 > longestprd) {
		    System.out.println("## CHECKX6 " + longestprd);
		    System.out.print("## hs3: "); printrowvec(hs3);
		}
	    }

	    boolean fastaprds = false;
	    if (fastaprds) {
		if (prdscore > 0) {
		    System.out.println(">" + nm + "-pprd");
		    System.out.println(aa2string(prd));
		    System.out.println(">" + nm + "-core");
		    System.out.println(aa2string(submatrix(aa, (int) hs3[0], (int) (hs3[1])) ));
		}
	    } else {
		// use one-based indices for output table (e.g, first AA is at position 1); zero-based is used internally.
		System.out.format("%s\t%d\t%d\t%d\t%d\t%.3f\t%d\t%d\t%d\t%.3f\t%d\t%.3f\t%d\t%d\t%d\t%.3f\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t",
				  nm,
				  (int) inf2nan(hs1[2]), 
				  (int) (hs1[0]+1), (int) (hs1[1]+1), // one-based
				  (int) (hs1[1]-hs1[0]+1), 
				  inf2nan(hs2[2]), 
				  (int) (hs2[0]+1), (int) (hs2[1]+1), // one-based
				  (int) (hs2[1]-hs2[0]+1) , 
				  inf2nan(hs2[2])/(hs2[1]-hs2[0]+1), longestprd,
				  inf2nan(hs3[2]), 
				  corestart+1, corestop+1, // one-based
				  corestop-corestart+1 , 
				  prdscore, 
				  aastart+1, aastop+1, // one-based
				  aastop-aastart+1 , 
				  aa.length, hmmscore, hmmscorev);
		if (aastop - aastart + 1 >= corelength) {
		    System.out.print(aa2string(submatrix(aa, corestart, corestop)));
		    System.out.print("\t");
		    System.out.print(aa2string(submatrix(aa, aastart, aastart+14)));
		    System.out.print("\t");
		    System.out.print(aa2string(submatrix(aa, aastop-14, aastop)));
		    System.out.print("\t");
		    System.out.print(aa2string(prd));
		} else {
		    System.out.print("-");
		    System.out.print("\t");
		    System.out.print("-");
		    System.out.print("\t");
		    System.out.print("-");
		    System.out.print("\t");
		    System.out.print("-");
		}
		System.out.format("\t%d\t%.3f\t%.3f\t%.3f\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%s",
				  dr.numdisorderedstrict2,
				  dr.meanhydro, 
				  dr.meancharge,
				  dr.meanfi,
				  dr.maxlen,
				  inf2nan(dr.papamaxscore),
				  dr.papamaxprop,
				  dr.papamaxdis,
				  dr.papamaxllr,
				  dr.papamaxllr2,
				  dr.papamaxcenter + 1, //one-based
				  aa2string(submatrix(aa,  dr.papamaxcenter - ww2/2, dr.papamaxcenter + ww2/2))); // ok for length < w? 
		System.out.println();
	    }
	   
	}

    }

    //static String translate(HashMap<String,String> ht, String key) {
    //	if (ht.containsKey(key)) {return (String) ht.get(key);}
    //	else return key;
    //}

    static int [][] trimrows(int [][] mat, int r1, int r2) {
	int [][] newmat = new int [r2-r1+1][];
	for (int i=0; i<r2-r1+1; i++) newmat[i] = mat[i+r1];
	return newmat;

    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////	       hmms				/////////////////////////

    // two-state HMM with differert fg and bg AA freqs
    static hmm prionhmm1(double [] fgfreq, double bgfreq []) {
	// state 0 = normal; state 1 = prion
	double [][] tmat = {{99.9/100, 0.1/100},{2.0/100, 98.0/100}};
	double [] imat = {0.9524, 0.0476}; // frequencies of stationary distribution (scaled principal eigenvector)
	double [][] emat = new double[2][];
        // CHANGED: added 0.00001 to * and X states in advance
	emat[0] = normalize(bgfreq); 
        emat[1] = normalize(fgfreq); 
	hmm h = new hmm(tmat,emat,imat);
	h.subtrellis = new int [] {1,0};
	h.states = new char[] {'-','+'};
	h.setnames(new String [] {"background","PrD-like"});
	return h;
    }

    // Set up as two states, but effectively just one state (bg AA freqs) ---
    // can be used to compute likelihood ratio of two state model vs background model
    // This global LLR isn't great for ranking proteins, though, as very long regions
    // with slight AA biases can out score modest-length regions with strong biases.
    // That is why we score based on max LLR in window of corelength.
    static hmm prionhmm0(double bgfreq []) {
	// state 0 = normal; state 1 = prion
	double [][] tmat = {{1, 0},{0, 1}};
	double [] imat = {1, 0}; // starts and stays in state 1 
	double [][] emat = new double[2][];
         // CHANGED: added 0.00001 to * and X states in advance
	emat[0] = normalize(bgfreq);
	emat[1] = normalize(bgfreq);
	hmm h = new hmm(tmat,emat,imat);
	h.subtrellis = new int [] {1,0};
	h.states = new char[] {'-','+'};
	h.setnames(new String [] {"background","also.background"});
	return h;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    // converts Infinity and -Infinity to NaN, which Excel may have an easier time with.
    static double inf2nan(double x) {
	return (Double.isInfinite(x)) ? 0.0/0.0 : x;
    }
     

    // computes log(exp(a) + exp(b));
    // okay if a or b is -Infinity, or both
    // a bit slow --- precompute lookup table?
    static double logeapeb2(double a, double b) {
	if (a > b) return (a + Math.log(1+Math.exp(b-a)));
	else if (b > a) return (b + Math.log(1+Math.exp(a-b)));
	else return (a + log2); // takes care of a = b = -Infinity

    }

    // computes log(exp(a) + exp(b)) with lookup-table
    static double logeapeb(double a, double b) {
	if (a > b) {
	    double c = a - b;
	    if (!(c < 40)) return a; // takes care of  b = -Inf
	    else {
	    	int dex = (int) Math.floor(100*c);
		//System.out.print("A " + a + "\t" + "\t" + b + "\t" + c + "\t" + dex + "\t"  + loglut[dex+1] + "\t" +	loglut[dex] + "\t");
	    	return (a + ((100*c - dex)*loglut[dex + 1] + (dex + 1 - 100*c)*loglut[dex]));

	    }
	    
	}
	else if (b > a) {
	    double c = b - a;
	    if (!(c < 40)) return b; // takes care of  a = -Inf
	    else {
	    	int dex = (int) Math.floor(100*c);
		// System.out.print("B " + a + "\t" + "\t" + b + "\t" + c + "\t" + dex + "\t"  + loglut[dex+1] + "\t" +	 loglut[dex] + "\t");
	    	return (b + ((100*c - dex)*loglut[dex + 1] + (dex + 1 - 100*c)*loglut[dex]));
	    }
	}
	else return (a + log2); // takes care of a = b = -Infinity

    }


    static double [] asdouble(int [] arr) {
	int n = arr.length;
	// System.out.println(n);
	double [] narr = new double[n];
	for (int i=0; i<n; i++) narr[i] = arr[i] + 0.0;
	return narr;
    }

    static double [][] asdouble(int [][] arr) {
	int n = arr.length;
	int m = arr[0].length;
	// System.out.println(n);
	double [][] narr = new double[n][m];
	for (int i=0; i<n; i++) {
	    for (int j=0; j<m; j++) {
		narr[i][j] = arr[i][j] + 0.0;
	    }
	}
	return narr;
    }

   
     // highest-scoring-subsequence -- brute force, order n*maxlen
    static double [] hss(double [] seq, int minlength, int maxlength) {
	int n = seq.length;

	double [] score = new double[3]; // start, end, score
        if (minlength > n || minlength > maxlength) {
	    score[0] = -1.0;
	    score[1] = -2.0;    // so len = stop - start + 1 = 0;
	    score[2] = -1.0/0;   // - Infinity
	    return score;
	}

       	// if (minlength > n) minlength = n;  // return -Inf instead
	if (maxlength > n) maxlength = n;

	double bestscore;
	double tempscore;
	double [] psum = new double[n+1];
	int beststart = 0;
	int beststop = minlength-1; 
	psum[0] = 0;
	// -Inf will cause probs here; cumsum with initial 0; {1,1,1,2} -> {0,1,2,3,5};
	for (int i=0; i<n; i++) {
	    psum[i+1] = psum[i] + seq[i];
	}

	bestscore = psum[minlength] - psum[0]; // (psum[0]=0, so same as psum[minlength])
	// say n=100, minlength=10; then last valid start index (base 0) would be 90 = 100-10, last segment 90:99
	for (int i=0; i<=n-minlength; i++) { // changed to < to <=; i is start index in seq
	    for (int j=i+minlength; j<=Math.min(i+maxlength,n); j++) { //  j is end index in psum; so j-1 is end index in seq 
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


    // highest-scoring-subsequence or any length
    // assumes at least one element is positive?
    static double [] hss2(double [] seq) {
	int n = seq.length;

	double [] score = new double[3]; // start, end, score
	double bestscore;
	int beststart = 0;
	int beststop = 0;
	int curstart = 0;

	double d = Math.max(seq[0], 0);
	bestscore = d;

	for (int i=1; i<n; i++) {
	    if (d+seq[i] > 0) {
		d = d + seq[i];
	    }
	    else {
		d = 0;
		curstart = i + 1;
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

    // highest-scoring-subsequence with given minimum length
    // assumes at least one element is positive?
    static double [] hss2(double [] seq, int minlength) {	
	int n = seq.length;
	double [] score = new double[3]; // start, end, score
	if (minlength > n ) {
	    score[0] = -1.0;
	    score[1] = -2.0;    // so length = stop - start + 1 = 0
	    score[2] = -1.0/0.0; //  -Infinity
	    return score;
	}

	double bestscore;
	int beststart = 0;
	int beststop = minlength-1;
	int curstart = 0;
	int newstart = 0;

	double [] psum = new double[n+1];

	psum[0] = 0;
	for (int i=0; i<n; i++) {
	    psum[i+1] = psum[i] + seq[i];
	}

	double d = psum[minlength];
	bestscore = d;

	for (int i = minlength; i<n; i++) {
	    d = psum[i+1] - psum[curstart];
	    newstart = curstart;
	    for (int j=curstart+1; j<i-minlength; j++) {
		if (psum[i+1]-psum[j] >= d) {
		    d = psum[i+1] - psum[j];
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


    // highest-scoring-subsequence with given minimum and maximum lengths
    // assumes at least one element is positive?
    static double [] hss2(double [] seq, int minlength, int maxlength) {
	int n = seq.length;
	double [] score = new double[3];   // start, end, score
	
	if (minlength > n || minlength > maxlength) {
	    score[0] = -1.0;
	    score[1] = -2.0; // so len = stop - start + 1 = 0;
	    score[2] = -1.0/0.0; // -Infinity instead?
	    return score;
	}
       
	if (maxlength > n) maxlength = n;  //  +1?

	double bestscore;
	int beststart = 0;
	int beststop = minlength-1;
	int curstart = 0;
	int newstart = 0;

	double [] psum = new double[n+1];

	psum[0] = 0;
	for (int i=0; i<n; i++) {
	    psum[i+1] = psum[i] + seq[i];
	}

	double d = psum[minlength];
	bestscore = d;

	for (int i = minlength; i<n; i++) {
	    if ((i-curstart) >= maxlength) curstart++;
	    d = psum[i+1] - psum[curstart];
	    newstart = curstart;
	    for (int j=curstart+1; j<i-minlength; j++) {
		if (psum[i+1] - psum[j] >= d) {
		    d = psum[i+1] - psum[j];
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


    // indices at which vec takes value target
    static int [] findindices(int [] vec, int target) {
	int n = vec.length;
	int hits = 0;
	for (int i=0; i<n; i++) if (vec[i] == target) hits++;
	int [] dexvec = new int[hits];
	hits = 0;
	for (int i=0; i<n; i++) {
	    if (vec[i] == target) {
		dexvec[hits]=i;
		hits++;
	    }
	}
	return dexvec;
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
		    if (aanames[k]==c) {
			mask[k] = seg; 
			break;
		    }
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
	} else {
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
	} else {
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
	int [][] mat2 = new int[m][n];
	for (int i=0; i<m; i++) {
	    for (int j=0; j<n; j++) { 
		mat2[i][j] = mat[i][n-j-1];
	    }
	}
	return mat2;


    }


    // partial sum
    static double [] cumsum(double [] arr) {
	int r = arr.length;
	double [] carr = new double[r];
	double cs = 0;
	for (int i=0; i<r; i++) {
	    cs = cs + arr[i];
	    carr[i] = cs;
	}
	return carr;
    }


    static int [][] submatrix(int [][] matrix, int r1, int r2, int c1, int c2) {
    	int m = matrix.length;
    	int n = matrix[0].length;
    	if (r1<0) r1 = 0;
    	if (r2<r1) r2 = r1;
    	if (c1<0) c1 = 0;
    	if (c2<c1) c2 = c1;
    	if (r1>=m) r1 = m-1;
    	if (r2>=m) r2 = m-1;
    	if (c1>=n) c1 = n-1;
    	if (c2>=n) c2 = n-1;
    	int [][] newmat = new int[r2-r1+1][c2-c1+1];
    	for (int i=0; i<r2-r1+1; i++) {
	    for (int j=0; j<c2-c1+1;j++) {
		newmat[i][j] = matrix[r1+i][c1+j];
	    }
	}
	return newmat;
    }


    static double  [][] submatrix(double [][] matrix, int r1, int r2, int c1, int c2) {
	int m = matrix.length;
	int n = matrix[0].length;
	if (r1<0) r1 = 0;
	if (r2<r1) r2 = r1;
	if (c1<0) c1 = 0;
	if (c2<c1) c2 = c1;
	if (r1>=m) r1 = m-1;
	if (r2>=m) r2 = m-1;
	if (c1>=n) c1 = n-1;
	if (c2>=n) c2 = n-1;
	double [][] newmat = new double[r2-r1+1][c2-c1+1];
	for (int i=0; i<r2-r1+1; i++) {
	    for (int j=0; j<c2-c1+1;j++) {
		newmat[i][j] = matrix[r1+i][c1+j];
	    }
	}
	return newmat;
    }

    static int [] submatrix(int [] matrix, int r1, int r2) {
	int m = matrix.length;

	if (r1<0) r1 = 0;
	if (r2<r1) r2 = r1;
	if (r1>=m) r1 = m-1;
	if (r2>=m) r2 = m-1;

	int [] newmat = new int[r2-r1+1];
	for (int i=0; i<r2-r1+1; i++) newmat[i] = matrix[r1+i];
	return newmat;
    }


    static double [] submatrix(double [] matrix, int r1, int r2) {
	int m = matrix.length;

	if (r1<0) r1 = 0;
	if (r2<r1) r2 = r1;
	if (r1>=m) r1 = m-1;
	if (r2>=m) r2 = m-1;

	double [] newmat = new double[r2-r1+1];
	for (int i=0; i<r2-r1+1; i++) newmat[i] = matrix[r1+i];
	return newmat;
    }

   

    // p*log(p)
    static double plogp(double p) {
	if (p == 0) {
	    return 0; 
	} else {
	    return p*Math.log(p); // change to base 2?
	}
    }

    // excise them instead?
    static void killlowercase(StringBuffer sb) {
	for (int j=0; j<sb.length(); j++) {
	    if (Character.isLowerCase(sb.charAt(j))) {
		sb.setCharAt(j,'-');
	    }
	}

    }

    static void printmatrix(double [] mat) {
	for (int i=0; i<mat.length; i++) {
	    System.out.println((float) (mat[i]));
	}
	System.out.println();
    }

    static void printrowvec(double [] mat) {
	for (int i=0; i<mat.length; i++) {
	    System.out.print((float) (mat[i]) + " ");
	}
	System.out.println();
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
	else if (c=='*') return 21; // no longer includes X
	else {
	    //System.out.println("bogus base " + c);
	    return 0;
	}
    }



    //  use binary search instead?
    static int char2int(char c, char [] arr) {
	for (int i=0; i<arr[i]; i++) {
	    if (arr[i]==c) return i;
	}
	return -1;

    }


    static int [] string2ints(StringBuffer sb, char [] chararr) {
	int n = sb.length();
	int [] arr = new int[n];
	for (int i=0; i<n; i++) arr[i] = char2int(sb.charAt(i),chararr);
	return arr;

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
	return (1.0*mn)/vec.length;
    }


    static double mean(double [] vec) {
      	double mn = 0;
      	for (int i=0; i<vec.length; i++) mn = mn + vec[i];
	return (1.0*mn)/vec.length;
    }

    // divides by n rather than n-1
    static double variance(double [] vec) {
      	double mn = mean(vec);
      	double var = 0;
      	for (int i=0; i<vec.length; i++) var = var + (vec[i]-mn)*(vec[i]-mn);
	return (1.0*var)/vec.length;
    }

    // divides by n rather than n-1
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
      	System.arraycopy(arr, 0, newarr, 0, n);
      	Arrays.sort(newarr);
	// median
      	if (n%2 == 0) {
	    stats[2] = (newarr[n/2] + newarr[n/2-1])/2.0;
	} else { 
	    stats[2] = newarr[(n-1)/2];
	}
      	int numtrim = (int) Math.round(alpha*n);
	newarr = submatrix(newarr, numtrim , n-numtrim-1); // +-1?
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
	if (n%2 == 0) {
	    return((0.0 + svec[n/2] + svec[n/2-1])/2);
	} else {
	    return (svec[(n-1)/2]);
	}
    }

 
    // singleton AA frequencies from fasta file
    static int [] computeaafreq(String filename, boolean killlc) {
	int [] aacounts = new int[22];
	fastareader fs = new fastareader(filename);
	while (fs.hasmorefastas()) {
	    String geneid = fs.nextname();
	    StringBuffer sb = fs.nextfasta();
	    if (killlc) killlowercase(sb);
	    countaas(sb, aacounts);
	}
	// printmatrix(aacounts);
	return aacounts;
    }

    // digram AA frequencies from fasta file
    static int [][] computeaafreq2(String filename, boolean killlc) {
	int [][] aacounts = new int[22][22];
	fastareader fs = new fastareader(filename);
	while (fs.hasmorefastas()) {
	    String geneid = fs.nextname();
	    StringBuffer sb = fs.nextfasta();
	    if (killlc) killlowercase(sb);
	    countaas2(sb, aacounts);
	}
	// printmatrix(aacounts);
	return aacounts;
    }
  
    // trigram AA frequencies from fasta file
    static int [][][] computeaafreq3(String filename, boolean killlc) {
	int [][][] aacounts = new int[22][22][22];
	fastareader fs = new fastareader(filename);
	while (fs.hasmorefastas()) {
	    String geneid = fs.nextname();
	    StringBuffer sb = fs.nextfasta();
	    if (killlc) killlowercase(sb);
	    countaas3(sb, aacounts);
	}
	// printmatrix(aacounts);
	return aacounts;
    }


    // Singleton AA frequencies for single seq
    static void countaas(StringBuffer buffer, int [] aacounts) {
	int [] aa;
	aa = string2aa(buffer);
	if (isvalidprotein(aa)) {
	    for (int i=0; i<aa.length; i++) {
		aacounts[aa[i]]++;
	    }
	}
    }

    // Digram AA frequencies for single seq
    static void countaas2(StringBuffer buffer, int [][] aacounts) {
	int [] aa;
	aa = string2aa(buffer);
	if (isvalidprotein(aa)) {
	    for (int i=0; i<aa.length-1; i++) {
		aacounts[aa[i]][aa[i+1]]++;
	    }
	}
    }

    // Trigram AA frequencies for single seq
    static void countaas3(StringBuffer buffer, int [][][] aacounts) {
	int [] aa;
	aa = string2aa(buffer);
	if (isvalidprotein(aa)) {
	    for (int i=0; i<aa.length-2; i++) {
		aacounts[aa[i]][aa[i+1]][aa[i+2]]++;
	    }
	}
    }


    // no premature stop codons or unknown residues. terminal stop okay.
    static boolean isvalidprotein(int [] aa) {
	int m = aa.length;
	for (int i=1; i<m-1; i++) {
	    if ((aa[i] == 0) || (aa[i] == 21)) return false;
	}
	if (aa[m-1] == 0) return false;
	return true;
    }

    // delta vector of length len with int value ival at index i
    static int [] deltavec(int len, int i, int ival) {
	int [] dv = new int[len];
	dv[i] = ival;
	return dv;
    }

    // delta vector of length len with double value ival at index i
    static double [] deltavec(int len, int i, double ival) {
	double [] dv = new double[len];
	dv[i] = ival;
	return dv;
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
      	for (int i=0; i<m; i++) aa[i] = aatoint(buffer.charAt(i));
      	return aa;
    }


    // 1 if upper-case; 0 otherwise
    static int [] buffer2bits(StringBuffer sb) {
      	int m = sb.length();
      	int [] bv = new int[m];
      	for (int i = 0; i<m; i++) {
	    if (Character.isUpperCase(sb.charAt(i))) {
		bv[i] = 1;
	    } else {
		bv[i] = 0;
	    }
      	}
      	return bv;
    }

    // longest run of 1 in 0-1 vec
    static int longestrun(int [] bitvec) {
      	int maxlen = 0;
      	int n = bitvec.length;
      	int i = 0; 
      	while (i < n) {
	    if (bitvec[i] > 0) {
		int startdex = i;
		i++;
		while (i < n && bitvec[i] > 0) i++;
		int stopdex = i-1;
		int len = stopdex-startdex+1;		
		if (len >= maxlen) maxlen=len;
	    } else { 
		i++;
	    }
      	}
      	return maxlen;
    }

    // longest run of 1 in 0-1 vec
    static int longestrun(double [] bitvec) {
      	int maxlen = 0;
      	int n = bitvec.length;
      	int i = 0; 
      	while (i < n) {
	    if (bitvec[i] > 0) {
		int startdex = i;
		i++;
		while (i < n && bitvec[i] > 0) i++;
		int stopdex = i-1;
		int len = stopdex - startdex + 1;		
		if (len >= maxlen) maxlen = len;
	    } else {
		i++;
	    }
      	}
      	return maxlen;
    }


    static StringBuffer aa2string(int [] aa) {
       	int m = aa.length;
       	StringBuffer buffer = new StringBuffer(m);
       	buffer.setLength(m);
       	for (int i=0; i<m; i++) buffer.setCharAt(i,aanames[aa[i]]);
       	return buffer;
    }


    static StringBuffer [] aa2string(int [][] aa) {
       	int n = aa.length;
       	StringBuffer [] buffer = new StringBuffer[n];
       	for (int j=0; j<n; j++) {
	    int m = aa[j].length;
	    buffer[j] = new StringBuffer(m);
	    buffer[j].setLength(m);
	    for (int i=0; i<m; i++) buffer[j].setCharAt(i, aanames[aa[j][i]]);   
       	}
       	return buffer;
    }


    static StringBuffer int2string(int [] arr, String [] names) {
       	int n = arr.length;
       	StringBuffer buffer = new StringBuffer(n);
	// buffer.setLength(n);
       	for (int i=0; i<n; i++) buffer.append(names[arr[i]]);
       	return buffer;
    }

    static StringBuffer int2string(int [] arr, char [] names) {
       	int n = arr.length;
       	StringBuffer buffer = new StringBuffer(n);
       	buffer.setLength(n);
       	for (int i=0; i<n; i++) buffer.setCharAt(i, names[arr[i]]);	
       	return buffer;
    }


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
		    if (defaultvalue == "incr") ht.put(chunks[0], new String(""+i));
		    else if (chunks.length > 1) ht.put(chunks[0], chunks[1]);
		    else if (defaultvalue == null) ht.put(chunks[0], chunks[0]);
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
	for (int i=0; i<n/w; i++) System.out.println(seq.substring(i*w, i*w + w));
	if (n%w > 0) System.out.println(seq.substring((n/w)*w)); //CHECKME priority
    }

    static int [] composition(int [] vec, int max) {
	int [] comp = new int[max+1];
	for (int i=0; i<vec.length; i++) comp[vec[i]]++;
	return comp;
    }

    static int [] composition(int [] vec) {
	return composition(vec, maxint(vec));
    }


    static int [] aacomposition(int [] vec) {
	return composition(vec, 22);
    }

  
    static double [] normalize(int [] arr) {
	int n = arr.length;
	double [] narr = new double [n];
	double sm = 1.0*sum(arr);
	for (int i=0; i<n; i++) narr[i] = arr[i]/sm;
	return narr;
    }

    // CHECK eps
    // if all entries in arr are 0, they will still be 0 after normalization (won't sum to zero)
    static double [] normalize(double [] arr) {
	double eps = 0.000000000001;
	int n = arr.length;
	double [] narr = new double [n];
	double sm = 1.0*sum(arr);
	if (sm < eps) sm = 1;
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

    // CHECK eps
    // makes total sum	1
    static double [][] normalizemat(double [][] mat) {
	double eps = 0.000000000001;
	int m = mat.length;
	int n = mat[0].length;
	double [][] nmat = new double[m][n];
	double ts = 0.0;
	for (int i=0; i<m; i++) {
	    for (int j=0; j<n; j++) {
		ts = ts+mat[i][j];
	    }
	}
	if (ts < eps) ts = eps;


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

    // ax + by, x and y vectors
    static double [] axpby(double a, double [] x, double b, double [] y) {
	int n = x.length;
	int m = y.length;
	if (m < n) n=m;
	double [] arr = new double[n];
	for (int i=0; i<n; i++) arr[i] = a*x[i] + b*y[i];
	return arr;
    }


    // ax + by, x and y matrices
    static double [][] axpby(double a, double [][] x, double b, double [][] y) {
	int m1 = x.length;
	int m2 = y.length;
	int n1 = x[0].length;
	int n2 = y[0].length;
	if (m1 > m2) m1 = m2;
	if (n1 > n2) n1 = n2;
	double [][] mat = new double[m1][n1];
	for (int i=0; i<m1; i++) {
	    for (int j=0; j<n1; j++) { 
		mat[i][j] = a*x[i][j] + b*y[i][j];
	    }
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
	    for (int j=0; j<mat[i].length; j++) {
		s = s + mat[i][j];
	    }
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
	    for (int i=0; i<m; i++) { 
		s = s + mat[i][j];
	    }
	    cs[j] = s;
	}
	return cs;
    }

    // ax + by + c
    static double [] axpbypc(double a, double [] x, double b, double [] y, double c) {
	int n = x.length;
	int m = y.length;
	if (m < n) n = m;
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

     // ax + b
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


    // matrix * matrix
    static double [][] matmult(double [][] mat1, double [][] mat2){
	int a = mat1.length;
	int b = mat1[0].length;
	int c = mat2.length;
	int d = mat2[0].length;
	if (b !=c )  {
	    System.out.println("# inner dimensions don't match: " + a + " x " + b + ", " +  c + " x " + d);
	    return (new double[1][1]);
	}
	double [][] mat3 = new double[a][d];
	double t;
	for (int i=0; i<a; i++) {
	    for (int j=0; j<d; j++) {
		t = 0;
		for (int k=0; k<b; k++) {
		    t = t + mat1[i][k]*mat2[k][j];
		}
		mat3[i][j] = t;
	    }
	}
	return mat3;
    }

    // matrix * vector
    static double [] matmult(double [][] mat1, double [] vec){
	int a = mat1.length;
	int b = mat1[0].length;
	int c = vec.length;
	if (b!=c)  {
	    System.out.println("# inner dimensions don't match: " + a + " x " + b + ", " +  c);
	    return (new double[1]);
	}
	double [] vec2 = new double[a];
	double t;
	for (int i=0; i<a; i++) {
	    t = 0;
	    for (int k=0; k<b; k++) {
		t = t + mat1[i][k]*vec[k];
	    }
	    vec2[i] = t;
	}
	return vec2;
    }


    static double [][] outerproduct(double [] vec1, double [] vec2) {
	int m = vec1.length;
	int n = vec2.length;
	double [][] op = new double[m][n];
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
	for (int i=0; i<m; i++)	ip = ip + vec1[i]*vec2[i];
	return ip;
    }


    static double [][] diagonal(double [] vec) {
	int n = vec.length;
	double [][] mat = new double[n][n];
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
	for (int i=0; i<n; i++) nvec[i] = vec1[i]/vec2[i]; // change if both 0?
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
	if (a != c || b != d)  {
	    System.out.println("# matrix dimensions don't match: "  + a + " x " + b + ", " +  c + " x " + d);
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


    // adjacency list of matrix entries > thresh
    static int [][] sparsesupport(double [][] mat, double thresh) {
	int m = mat.length;
	int n = mat[0].length;
	int [][] ss = new int[m][];
	for (int i=0; i<m; i++) {
	    int ct = 0;
	    for (int j=0; j<n; j++) {
		if (mat[i][j]>thresh) {
		    ct++;
		}
	    }
	    ss[i] = new int[ct];
	    ct = 0;
	    for (int j=0; j<n; j++) {
		if (mat[i][j]>thresh) {
		    ss[i][ct] = j;
		    ct++;
		}
	    }
	}
	return ss;
    }


    
    static int [] sparsesupport(int [] vec, int target) {
	int m = vec.length;
	int ct = 0;
	for (int i=0; i<m; i++) {
	    if (vec[i]==target) {
		ct++;
	    }
	}
	int [] ss = new int[ct];
	ct = 0;
	for (int i=0; i<m; i++) {
	    if (vec[i]==target) {
		ss[ct] = i; 
		ct++;
	    }
	}
	return ss;
    }

    // elementwise log
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

    // elementwise log
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
	for (int i=0; i<r; i++) {
	    for (int j=0; j<c; j++) {
		mat[i][j] = val;
	    }
	}
	return mat;
    }

    static double [] constantmat(int r, double val) {
	double [] mat = new double[r];
	for (int i=0; i<r; i++)	 mat[i] = val;
	return mat;
    }


    // elementwise inverse
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

    // returns val and index --- index is double but can be cast to int
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


    // added 08/09/2011
    // need not be square here --- allow row and column names???
    static double [][] readmatrix2(String filename, int skiprow, int skipcol) {
	BufferedReader in;
	String line;
	int n = 0; 
	int m = 0; 

	double [][] mat = new double[1][1];

	StringTokenizer toker;

	try {
	    in = new BufferedReader(new FileReader(filename));
	    for (int i=0; i<skiprow; i++) {
		line = in.readLine();
	    }
	    line = in.readLine();
	    toker = new StringTokenizer(line);
	    n = toker.countTokens() - skipcol;
	    m = 0;
	    while((line != null) && (toker.countTokens() - skipcol == n)) { // make sure toker is non-null?
	    	m++;
	    	line = in.readLine();
	    	if (line != null) toker = new StringTokenizer(line);
	    }
	    in.close();
	}
	catch (IOException e) {
	    System.out.println("Couldn't open " + filename);
	}
        // System.out.println("## dim " + m + " by " + n);


	try {
	    in = new BufferedReader(new FileReader(filename));
	    for (int i=0; i<skiprow; i++) {
		line = in.readLine();
	    }
	    line = in.readLine();
	    toker = new StringTokenizer(line);
	    mat = new double[m][n];
	    for (int i=0; i<m; i++) {
		if (i > 0) { // first line is already read
		    line = in.readLine();
		    toker = new StringTokenizer(line);
		}
		for (int j=0; j<skipcol; j++) {
		    toker.nextToken();
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
	for (int i=0; i<freq.length; i++) {
	    newfreq[map[i]] = newfreq[map[i]]+freq[i];
	}
	return newfreq;
    }


    static int [] mapseq(int [] arr, int [] map) {
	int n = arr.length;
	int [] newarr = new int[n];
	for (int i=0; i<n; i++) {
	    newarr[i] = map[arr[i]];
	}
	return newarr;
    }


    static double [] mapseq(int [] arr, double [] map) {
	int n = arr.length;
	double [] newarr = new double[n];
	for (int i=0; i<n; i++) {
	    newarr[i] = map[arr[i]];
	}
	return newarr;
    }

 
    // window size ww should be odd. 
    static double [] slidingaverage(double [] arr, int ww, boolean shrink, boolean weight) {
	int n = arr.length;
	if (n == 0) return arr;
	int w = ww/2; // floor if odd
	if ( w>= n) w=n-1; // necessary ??
	double [] sa = new double[n];
	double score;
	double denom;
        double wt;
        int mini, maxi;
        if (shrink) {
	    mini = 0;
	    maxi = n-1;
	} else {
	    mini = w;
	    maxi = n-w-1;
	    for (int i=0; i<mini; i++) sa[i] = 0.0/0.0;
	    for (int i=maxi+1; i<n; i++) sa[i] = 0.0/0.0;
	}
	for (int i=mini; i<=maxi; i++) {
	    score = 0.0;
	    denom = 0.0;
	    for (int j=-w; j<=w; j++) { // should tighten bounds to avoid next test
        	if ((i + j >= 0) && (i + j < n)) { 
		    wt = 1.0;
                    if (weight) {
			wt = 1.0 + Math.min(i+j, w) + Math.min(n-i-j-1, w);
			//System.out.println(n + "\t" + ww + "\t" + i + "\t" + j + "\t" + wt);
		    }
		    // if (weight) wt = Math.max(i+j, w) + Math.max(n-i-j-1, w);
		    denom = denom + wt;
		    score = score + wt*arr[i + j];
		}
	    }
	    sa[i] = score/denom;   
	}
	return sa;
    }

    
    //  if mergeme=Z, scores only first Z in ZZ or Z.Z in each window
    static double [] slidingaverage(double [] arr, int ww, boolean shrink, boolean weight, int mergeme, int [] seq) {
	int n = arr.length;
	if (n == 0) return arr;
      	int w = ww/2; // floor if odd
	if ( w >= n) w=n-1; // necessary ??
	double [] sa = new double[n];
	double score;
	double denom;
	double wt;
	int mini, maxi;
	if (shrink) {
	    mini = 0;
	    maxi = n-1;
	} else {
	    mini = w;
	    maxi = n-w-1;
	    for (int i=0; i<mini; i++) sa[i] = 0.0/0.0;
	    for (int i=maxi+1; i<n; i++) sa[i] = 0.0/0.0;
	}
	for (int i=mini; i<=maxi; i++) {
	    score = 0.0;
	    denom = 0.0;
	    for (int j=-w; j<=w; j++) { // should tighten bounds to avoid next test
		if ((i + j >= 0) && (i + j < n)) {
		    wt = 1.0;
		    if (weight)	wt = 1.0 + Math.min(i+j, w) + Math.min(n-i-j-1, w);
		    denom = denom + wt;
		    if (!((seq[i + j] == mergeme) && (i + j >= 1) && (seq[i + j - 1] == mergeme))
			  && !((seq[i + j] == mergeme) && (i + j >= 2) && (seq[i + j - 2] == mergeme))) {
			score = score + wt*arr[i + j];
		    } 
		} 
	    } 
	    sa[i] = score/denom;
	}
	return sa;
    }

    // companion to read_aa_params
    static void print_aa_params(double [] aaparams) {
	for (int i=0; i<22; i++) {
	    // System.out.println((float) aaparams[i] + " # " + aanames[i]);
	    System.out.format("%.6f # %s\n", aaparams[i], aanames[i]);
	}
    }

    static void print_aa_params(int [] aaparams) {
	for (int i=0; i<22; i++) {
            System.out.format("%.6f # %s\n", aaparams[i], aanames[i]);
	}
    }
   
    
    // companion to print_aa_params.
    // file should be 22 numbers, one per line, in order of string array aanames. 
    // {'X','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'};
    // optionally each line can be: param # name, where if present the first letter of the name 
    // on the ith line will be checked for agreement with aanames[i] (counting lines from index zero). 
    static double [] read_aa_params(String filename) {
	BufferedReader in;
	String line;
	double [] vec = new double[22];
	StringTokenizer toker;
	String name; 

	try {
	    in = new BufferedReader(new FileReader(filename));
	    for (int i=0; i<22; i++) {
	       line = in.readLine();
	       toker = new StringTokenizer(line);
	       int n = toker.countTokens();
	       vec[i] = Double.parseDouble(toker.nextToken());
	       if (n > 2) {
		   toker.nextToken(); // should be #
		   name = toker.nextToken();
		   if (name.charAt(0) != aanames[i]) {
		       System.out.println("# warning: " + filename + " does not have expected name in line" + (i+1));
		   }
	       }
	    }
	    in.close();
	}
	catch (IOException e) {
	    System.out.println("# Couldn't open " + filename);
	}
	// FIXME: catch other sorts of exceptions for invalid formats, not enough lines, ...
	return vec;
    }


 
    static String aaparams2string(double [] aafreq) {
        String aastring = new String("");
	for (int i=0; i<22; i++) {
            aastring = aastring + String.format("%s=%.5f;", aanames[i], aafreq[i]);
	}
        return aastring;
     }     



    // starts with comment char so can be ignored by R
    static void printaausage(double [] aafreq) {
	for (int i=0; i<22; i++) {
	    // System.out.println("# " + aanames[i] + " " + (float) aafreq[i]);
	    System.out.format("# %s %.6f\n", aanames[i], aafreq[i]);
	}
    }

    static void printaausage(double [][] aafreq) {
	for (int i=0; i<22; i++) {
	    System.out.print("# " + aanames[i]);
	    for (int j=0; j<aafreq.length; j++) {
		//System.out.print("\t" + (float) aafreq[j][i]);
               	System.out.format("\t%.6f",aafreq[j][i]);
	    }
	    System.out.println();
	}
    }

    // static void printaausage(double [][] aafreq, int sd) {
    //	for (int i=0; i<22; i++) {
    //	    System.out.print("# " + aanames[i]);
    //	    for (int j=0; j<aafreq.length; j++) {
    //		System.out.print("\t" + (float) (Math.round(sd*aafreq[j][i])/(1.0*sd)));
    //	    }
    //	    System.out.println();
    //	}
    //}


}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// This class implements Viterbi and MAP decoding of hidden markov models (HMMs).
// For an introduction to HMMS and descriptions of the algorithms, see e.g. Chapter 3 of
//   Durbin R, Eddy S, Krogh A, Mitchison G. 
//   "Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids",
//   Cambridge University Press, 1998.
// As noted in 3.6, straightforward implementations of decoding algorithms are prone to
// numerical underflow, which can be addressed by using log-probabilties or scaling parameters.
// Versions below that opeate in log-space end in l, and those that use scaling end in s 
// (e.g. naive implementation posterio vs posteriorl vs posteriors).
// It's not needed for the two state HMMs, but for HMMs with many states and relatively 
// few direct transitions allowed, it can be faster to use a sparse matrix representations 
// of the transition matrix. Versions that use sparse representations start with s, e.g. 
// sposterior.
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
    double []   ciprob;
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

    // 0-1 mask of length ns
    int [] subtrellis; 

    // multiple states can belong to same class, with probabilities combined in post-processing
    int [] classes; 
    int numclasses;

    double etst = 0.0; // expected time on subtrellis
    double lpst = 0.0; // log prob of never visiting subtrellis

    double lmp = 0.0; // log of marginal prob of output

    double [][] postprob;
    double [] margcollapse;

    int [] viterbipath;
    int [] mappath;

    boolean sparse = false;

    double sparseeps =  0.000000001; // in converting to sparse matrix, ignore transition probs smaller than this. 
                                   // depending on application may want to change this.

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
	ctprob = plaac.cumsum(tprob,true);
	ceprob = plaac.cumsum(eprob,true);
	ciprob = plaac.cumsum(iprob);

	ltprob = plaac.logmat(tprob);
	leprob = plaac.logmat(eprob);
	liprob = plaac.logmat(iprob);
	ns = iprob.length;
	no = eprob[0].length;

	double [] rs = plaac.rowsums(tprob);
	fprob = new double[ns];

        // if freeend = true, don't model sequence length (transition to end state is free);
        // otherwise, sequence length is modeled as exponential with decay given by 1 - sum of 
        // transition probabilities

	boolean freeend = true; 

	for (int i=0; i<ns; i++) {
	    fprob[i] = Math.max(0.0, 1.0 - rs[i]);
	    if (fprob[i] > 0.0001) freeend = false; // may want to change this
	}
	if (freeend) {
	    for (int i=0; i<ns; i++) {
		fprob[i] = 1.0;
	    }
	}
	lfprob = plaac.logmat(fprob);

	subtrellis = new int[ns];
	for (int i=0; i<ns; i++) subtrellis[i] = 1;
	states = new char[ns];
	for (int i=0; i<ns; i++) states[i] = (char) i;
	names = new String[ns];
	for (int i=0; i<ns; i++) names[i] = new String("s"+i);
	classes = new int[ns];
	for (int i=0; i<ns; i++) classes[i] = i;
	numclasses = ns;
	ustates = states;
	unames = names;
    }


    public void setclasses(int [] classes) {
	this.classes = classes;
	int nc = 0;
	numclasses = 0;
	for (int i=0; i<classes.length; i++) {
	    if (classes[i]>nc) nc = classes[i];
	}
	numclasses = nc+1;
	// System.out.println("nc " + numclasses);
	// plaac.printmatrix(classes);
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

    // should add state names etc to output file, modify hmm(String) to read them
    public void print() {
	System.out.println(ns + " \t ## num states");
	System.out.println(no + " \t ## num symbols");
	System.out.println();
	for (int i=0; i<ns; i++) {
	    for (int j=0; j<ns; j++) {
		System.out.format("%.6f\t", tprob[i][j]);
	    }
	    System.out.println();
	}
	System.out.println();

	for (int j=0; j<ns; j++) {
	    System.out.format("%.6f\t", iprob[j]);
	}
	System.out.println();
	System.out.println();

	for (int j=0; j<no; j++) {
	    for (int i=0; i<ns; i++) {
		System.out.format("%.6f\t",eprob[i][j]); 
	    }
	    System.out.println();
	}
	System.out.println();

    }

    // should add state names etc to output file, modify hmm(String) to read them
    public void write(String filename) {
	try {
	    BufferedWriter out = new BufferedWriter(new FileWriter(filename));

	    out.write(ns + " \t ## num states");  
	    out.newLine();
	    out.write(no + " \t ## num symbols"); 
	    out.newLine();
	    out.newLine();
	    for (int i=0; i<ns; i++) {
		for (int j=0; j<ns; j++) {
		    out.write(String.format("%.6f\t", tprob[i][j])); // omit last tab?
		}
		out.newLine();
	    }
	    out.newLine();

	    for (int j=0; j<ns; j++) {
		out.write(String.format("%.6f\t", iprob[j]));
	    }
	    out.newLine();
	    out.newLine();

	    for (int j=0; j<no; j++) {
		for (int i=0; i<ns; i++) {
                    out.write(String.format("%.6f\t",eprob[i][j])); 
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

	for (int i=0; i<ns; i++) {
	    s[i][0] = iprob[i]*eprob[i][seq[0]];
	}
	for (int t=1; t<n; t++) {
	    for (int i=0; i<ns; i++) {
		int bestdex = 0;
		double bestscore = tprob[0][i]*s[0][t-1];
		for (int k=1; k<ns; k++) {
		    if (tprob[k][i]*s[k][t-1] > bestscore) {
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
	    if (s[k][n-1] > bestscore) {
		bestscore = s[k][n-1]; 
		bestdex = k;
	    }
	}
	vit[n-1] = bestdex;
	for (int t=n-2; t>=0; t--) {
	    vit[t] = tb[vit[t+1]][t+1];
	}

	plaac.printmatrix(plaac.transpose(s)); 
	System.out.println();

	viterbipath = vit;
	return vit;
    }


    public int [] viterbidecodel(int [] seq) {
	int n = seq.length;
	int [] vit = new int[n];

	double [][] s = new double[ns][n];
	int [][] tb = new int[ns][n];	    // traceback

	for (int i=0; i<ns; i++) {
	    s[i][0] = liprob[i]+leprob[i][seq[0]];
	} // initial probs
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
	    if (s[k][n-1] + lfprob[k] > bestscore) {
		bestscore = s[k][n-1] + lfprob[k];
		bestdex = k;
	    }
	}
	vit[n-1] = bestdex;

	for (int t=n-2; t>=0; t--) {
	    vit[t] = tb[vit[t+1]][t+1];
	}

	// plaac.printmatrix(plaac.transpose(s)); System.out.println();

	lviterbiprob = bestscore;
	if (numclasses < ns) vit = plaac.mapseq(vit, classes);
	viterbipath = vit;
	return vit;
    }

 
    // sparse transmat
    public int [] sviterbidecodel(int [] seq) {
	double eps = sparseeps;
	int n = seq.length;
	int [] vit = new int[n];

	double [][] s = new double[ns][n];
	int [][] tb = new int[ns][n];	    // traceback

	int [][] sps = plaac.sparsesupport(plaac.transpose(tprob), eps);


	for (int i=0; i<ns; i++) {
	    s[i][0] = liprob[i]+leprob[i][seq[0]];
	} // initial probs
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
	    if (s[k][n-1] + lfprob[k] > bestscore) {
		bestscore = s[k][n-1] + lfprob[k]; 
		bestdex = k;
	    }
	}
	vit[n-1] = bestdex;

	for (int t=n-2; t>=0; t--) {
	    vit[t] = tb[vit[t+1]][t+1];
	}

	// plaac.printmatrix(plaac.transpose(s)); System.out.println();

	lviterbiprob = bestscore;

	if (numclasses < ns) vit = plaac.mapseq(vit, classes);
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

	if (lab[0] < 0) {
	    for (int i=0; i<ns; i++) {
		s[i][0] = liprob[i]+leprob[i][seq[0]];
	    }
	} else {
	    for (int i=0; i<ns; i++) {
		s[i][0] = badscore;
	    }
	    s[lab[0]][0] = liprob[lab[0]]+leprob[lab[0]][seq[0]];
	}

	for (int t=1; t<n; t++) {
	    if (lab[t] < 0) {
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
	    if (s[k][n-1] > bestscore) {
		bestscore = s[k][n-1]; 
		bestdex = k;
	    }
	}
	vit[n-1] = bestdex;
	for (int t=n-2; t>=0; t--) {
	    vit[t] = tb[vit[t+1]][t+1];
	}

	// plaac.printmatrix(plaac.transpose(s)); System.out.println();

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
		if (subtrellis[i]==0) etst = etst + plaac.sum(postprob[i]);
	    }
	}
	else {
	    margcollapse = postprob[0];
	    etst = plaac.sum(postprob[0]);
	}

	lpst = logprobsubtrellis(seq, subtrellis);

    }


    public double [][] posterior(int [] seq) {
	int n = seq.length;
	double [][] pp = new double[ns][n];

	// forward algorithm
	double [][] a = new double[ns][n];

	for (int i=0; i<ns; i++) {
	    a[i][0] = iprob[i]*eprob[i][seq[0]];
	}
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

	for (int i=0; i<ns; i++) {
	    b[i][n-1] = 1;
	} // iprob[i];
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
	for (int i=0; i<ns; i++) {
	    pseq = pseq + a[i][n-1]*b[i][n-1];
	}

	// System.out.println("> " + pseq + " " + (a[0][n-1] + a[1][n-1]));


	for (int t=0; t<n; t++) {
	    for (int i=0; i<ns; i++) {
		pp[i][t] = a[i][t]*b[i][t]/pseq;
	    }
	}

	// System.out.println();System.out.println("a"); plaac.printmatrix(plaac.transpose(a));
	// System.out.println();System.out.println("b"); plaac.printmatrix(plaac.transpose(b));
	// System.out.println();


	if (numclasses<ns) pp = collapseposteriors(pp, classes, numclasses);
	postprob = pp;
	return pp;

    }

    // in log-space
    public double [][] posteriorl(int [] seq) {
	int n = seq.length;
	double [][] pp = new double[ns][n];

	// forward algorithm
	double [][] a = new double[ns][n];
	a = plaac.logmat(a);

	for (int i=0; i<ns; i++) {
	    a[i][0] = liprob[i] + leprob[i][seq[0]];
	}
	for (int t=1; t<n; t++) {
	    for (int i=0; i<ns; i++) {
		double score = -1.0/0.0;
		for (int k=0; k<ns; k++) {
		    score = plaac.logeapeb(score, ltprob[k][i] + a[k][t-1]);
		}
		a[i][t] = score + leprob[i][seq[t]];
	    }
	}
	double ltotprob = -1.0/0.0;
	for (int i=0; i<ns; i++) {
	    // a[i][n-1] = a[i][n-1] + lfprob[i];
	    ltotprob = plaac.logeapeb(ltotprob,a[i][n-1] + lfprob[i]);
	}

	lmarginalprob = ltotprob;

	// backward algorithm
	double [][] b = new double[ns][n];

	for (int i=0; i<ns; i++) {
	    b[i][n-1] = lfprob[i];
	} // accounts for termination
	for (int t=n-2; t>=0; t--) {
	    for (int i=0; i<ns; i++) {
		double score = -1.0/0.0;
		for (int k=0; k<ns; k++) {
		    score = plaac.logeapeb(score, ltprob[i][k] + b[k][t+1] + leprob[k][seq[t+1]]);
		}
		b[i][t] = score;
	    }
	}

	double lpseq = -1.0/0.0;
	for (int i=0; i<ns; i++) {
	    lpseq = plaac.logeapeb(lpseq, a[i][0] + b[i][0]);
	}

	// lmarginalprob = lpseq;

	//CHECKME
	for (int t=0; t<n; t++) {
	    for (int i=0; i<ns; i++) {
		pp[i][t] = Math.exp((a[i][t] + b[i][t]) - lpseq);
	    }
	}

	if (numclasses<ns) pp = collapseposteriors(pp, classes, numclasses);
	postprob = pp;
	return pp;

    }

    // in log-space, sparse transmat
    public double [][] sposteriorl(int [] seq) {
	double eps = sparseeps;
	int n = seq.length;
	double [][] pp = new double[ns][n];

	// forward algorithm
	double [][] a = new double[ns][n];
	a = plaac.logmat(a);

	int [][] sps = plaac.sparsesupport(plaac.transpose(tprob), eps);
	int [][] sps2 = plaac.sparsesupport(tprob, eps);

	for (int i=0; i<ns; i++) {
	    a[i][0] = liprob[i] + leprob[i][seq[0]];
	}
	for (int t=1; t<n; t++) {
	    for (int i=0; i<ns; i++) {
		double score = -1.0/0.0;
		for (int sk=0; sk < sps[i].length; sk++) {
		    int k =sps[i][sk];
		    // for (int k=0; k<ns; k++) {
		    score = plaac.logeapeb(score, ltprob[k][i] + a[k][t-1]);
		}
		a[i][t] = score + leprob[i][seq[t]];
	    }
	}
	double ltotprob = -1.0/0.0;
	for (int i=0; i<ns; i++) {
	    // a[i][n-1] = a[i][n-1] + lfprob[i];
	    ltotprob = plaac.logeapeb(ltotprob,a[i][n-1] + lfprob[i]);
	}

	lmarginalprob = ltotprob;


	// backward algorithm
	double [][] b = new double[ns][n];

	for (int i=0; i<ns; i++) {
	    b[i][n-1] = lfprob[i];
	} // accounts for termination
	for (int t=n-2; t>=0; t--) {
	    for (int i=0; i<ns; i++) {
		double score = -1.0/0.0;
		for (int sk=0; sk < sps2[i].length; sk++) {
		    int k =sps2[i][sk];
		    // for (int k=0; k<ns; k++) {
		    score = plaac.logeapeb(score, ltprob[i][k] + b[k][t+1] + leprob[k][seq[t+1]]);
		}
		b[i][t] = score;
	    }
	}

	double lpseq = -1.0/0.0;
	for (int i=0; i<ns; i++) {
	    lpseq = plaac.logeapeb(lpseq, a[i][0] + b[i][0]);
	}

	// lmarginalprob = lpseq;

	// CHECKME
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

	for (int i=0; i<ns; i++) {
	    a[i][0] = sc[0]*iprob[i]*eprob[i][seq[0]];
	}
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
	    for (int i=0; i<ns; i++) {
		a[i][t] = a[i][t]*sc[t];
	    }
	}

	// backward algorithm
	double [][] b = new double[ns][n];

	for (int i=0; i<ns; i++) {
	    b[i][n-1] = sc[n-1];
	} // 1; // iprob[i];
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
	for (int i=0; i<ns; i++) {
	    pseq = pseq + a[i][n-1]*b[i][n-1];
	}
	double lpseq = Math.log(pseq);
	for (int t=0; t<n; t++) {
	    lpseq = lpseq - Math.log(sc[t]);
	}

	for (int t=0; t<n; t++) {
	    double denom = 0;
	    for (int i=0; i<ns; i++) {
		denom = denom + a[i][t]*b[i][t];
	    }
	    for (int i=0; i<ns; i++) {
		pp[i][t] = a[i][t]*b[i][t]/denom;
	    }
	}

	System.out.println();
	System.out.println("a"); 
	plaac.printmatrix(plaac.transpose(a));
	System.out.println();
	System.out.println("b"); 
	plaac.printmatrix(plaac.transpose(b));

	lmarginalprob = lpseq;
	lmp = lpseq;
	if (numclasses<ns) pp = collapseposteriors(pp, classes, numclasses);
	postprob = pp;
	return pp;

    }

    // sparse transmat; use scaling to prevent underflow
    public double [][] sposteriors(int [] seq, double [][] tprob, double [][] eprob, double [] iprob) {
	
	double eps = sparseeps;

	int n = seq.length;
	double [][] pp = new double[ns][n];

	// scaling factor
	double [] sc = new double[n];
	for (int i=0; i<n; i++) sc[i] = 1.0;

	// forward algorithm
	double [][] a = new double[ns][n];

	int [][] sps = plaac.sparsesupport(plaac.transpose(tprob), eps);
	int [][] sps2 = plaac.sparsesupport(tprob, eps);

	for (int i=0; i<ns; i++) {
	    a[i][0] = sc[0]*iprob[i]*eprob[i][seq[0]];
	}
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
	    for (int i=0; i<ns; i++) {
		a[i][t] = a[i][t]*sc[t];
	    }
	}

	double totprob = 0;

	for (int i=0; i<ns; i++) {
	    totprob = totprob + a[i][n-1]*fprob[i];
	}

	double ltotprob = Math.log(totprob);
	for (int t=0; t<n; t++) {
	    ltotprob = ltotprob - Math.log(sc[t]);
	}

	lmarginalprob = ltotprob;


	// backward algorithm
	double [][] b = new double[ns][n];

	for (int i=0; i<ns; i++) {
	    b[i][n-1] = sc[n-1];
	} // 1; // iprob[i];
	for (int t=n-2; t>=0; t--) {
	    for (int i=0; i<ns; i++) {
		double score = 0;
		for (int sk=0; sk<sps2[i].length; sk++) {
		    // for (int k=0; k<ns; k++) {
		    int k = sps2[i][sk];
		    score = score + tprob[i][k]*b[k][t+1]*eprob[k][seq[t+1]];
		}
		b[i][t] = sc[t]*score;
	    }
	}

	double pseq = 0;
	for (int i=0; i<ns; i++) {
	    pseq = pseq + a[i][n-1]*b[i][n-1];
	}
	double lpseq = Math.log(pseq);
	for (int t=0; t<n; t++) {
	    lpseq = lpseq - Math.log(sc[t]);
	}


	for (int t=0; t<n; t++) {
	    double denom = 0;
	    for (int i=0; i<ns; i++) {
		denom = denom + a[i][t]*b[i][t];
	    }
	    for (int i=0; i<ns; i++) {
		pp[i][t] = a[i][t]*b[i][t]/denom;
	    }
	}

	// System.out.println();System.out.println("a"); plaac.printmatrix(plaac.transpose(a));
	// System.out.println();System.out.println("b"); plaac.printmatrix(plaac.transpose(b));
	// System.out.println();

	// lmarginalprob = lpseq;
	lmp = lpseq;
	if (numclasses<ns) pp = collapseposteriors(pp, classes, numclasses);
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

	for (int i=0; i<ns; i++) {
	    a[i][0] = sc[0]*iprob[i]*eprob[i][seq[0]];
	}
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
	    for (int i=0; i<ns; i++) { 
		a[i][t] = a[i][t]*sc[t];}
	}

	// backward algorithm
	double [][] b = new double[ns][n];

	for (int i=0; i<ns; i++) {
	    b[i][n-1] = sc[n-1];
	} // 1; // iprob[i];
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
	for (int i=0; i<ns; i++) {
	    pseq = pseq + a[i][n-1]*b[i][n-1];
	}
	double lpseq = Math.log(pseq);
	for (int t=0; t<n; t++) {
	    lpseq = lpseq - Math.log(sc[t]);
	}


	for (int t=0; t<n; t++) {
	    double denom = 0;
	    for (int i=0; i<ns; i++) {
		denom = denom + a[i][t]*b[i][t];
	    }
	    for (int i=0; i<ns; i++) {
		pp[i][t] = a[i][t]*b[i][t]/denom;
	    }
	}

	// System.out.println();System.out.println("a"); plaac.printmatrix(plaac.transpose(a));
	// System.out.println();System.out.println("b"); plaac.printmatrix(plaac.transpose(b));
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


    // // allow partial labels --- not tested thoroughly
    // public double [][] posteriors(int [] seq, int [] lab, double [][] tprob, double [][] eprob, double [] iprob) {

    // 	int n = seq.length;
    // 	double [][] pp = new double[ns][n];

    // 	// scaling factor
    // 	double [] sc = new double[n];
    // 	for (int i=0; i<n; i++) sc[i] = 1.0;

    // 	// forward algorithm
    // 	double [][] a = new double[ns][n];

    // 	int lb, ub;

    // 	if (lab[0]<0) { 
    // 	    for (int i=0; i<ns; i++) {
    // 		a[i][0] = sc[0]*iprob[i]*eprob[i][seq[0]]; 
    // 	    }
    // 	} else {
    // 	    a[lab[0]][0] = sc[0]*eprob[lab[0]][seq[0]];
    // 	}

    // 	for (int t=1; t<n; t++) {
    // 	    double sf = 0.0;

    // 	    if (lab[t]<0) {
    // 		for (int i=0; i<ns; i++) {
    // 		    double score = 0;
    // 		    for (int k=0; k<ns; k++) {
    // 			score = score + tprob[k][i]*a[k][t-1];
    // 		    }
    // 		    a[i][t] = score*eprob[i][seq[t]];
    // 		    sf = sf + score*eprob[i][seq[t]];
    // 		}
    // 	    } else {
    // 		double score = 0;
    // 		for (int k=0; k<ns; k++) {
    // 		    score = score + tprob[k][lab[t]]*a[k][t-1];
    // 		}
    // 		a[lab[t]][t] = score*eprob[lab[t]][seq[t]];
    // 		sf = sf + score*eprob[lab[t]][seq[t]];
    // 	    }

    // 	    sc[t] = 1/sf;
    // 	    for (int i=0; i<ns; i++) {
    // 		a[i][t] = a[i][t]*sc[t];
    // 	    }
    // 	}

    // 	// backward algorithm
    // 	double [][] b = new double[ns][n];

    // 	for (int i=0; i<ns; i++) {
    // 	    b[i][n-1] = sc[n-1];
    // 	} // 1; // iprob[i];
    // 	for (int t=n-2; t>=0; t--) {

    // 	    if (lab[t+1]<0) {
    // 		for (int i=0; i<ns; i++) {
    // 		    double score = 0;
    // 		    for (int k=0; k<ns; k++) {
    // 			score = score + tprob[i][k]*b[k][t+1]*eprob[k][seq[t+1]];
    // 		    }
    // 		    b[i][t] = sc[t]*score;
    // 		}
    // 	    } else {
    // 		for (int i=0; i<ns; i++) {
    // 		    double score = 0;
    // 		    score = score + tprob[i][lab[t+1]]*b[lab[t+1]][t+1]*eprob[lab[t+1]][seq[t+1]];
    // 		    b[i][t] = sc[t]*score;
    // 		}
    // 	    }


    // 	}

    // 	double pseq = 0;
    // 	for (int i=0; i<ns; i++) {
    // 	    pseq = pseq + a[i][n-1]*b[i][n-1];
    // 	}
    // 	double lpseq = Math.log(pseq);
    // 	for (int t=0; t<n; t++) {
    // 	    lpseq = lpseq - Math.log(sc[t]);
    // 	}


    // 	for (int t=0; t<n; t++) {
    // 	    double denom = 0;
    // 	    for (int i=0; i<ns; i++) {
    // 		denom = denom + a[i][t]*b[i][t];
    // 	    }
    // 	    for (int i=0; i<ns; i++) {
    // 		pp[i][t] = a[i][t]*b[i][t]/denom;
    // 	    }
    // 	}

    // 	// System.out.println();System.out.println("a"); plaac.printmatrix(plaac.transpose(a));
    // 	// System.out.println();System.out.println("b"); plaac.printmatrix(plaac.transpose(b));
    // 	// System.out.println();

    // 	lmarginalprob = lpseq;
    // 	lmp = lpseq;

    // 	if (numclasses<ns) pp = collapseposteriors(pp, classes, numclasses);


    // 	postprob = pp;
    // 	return pp;

    // }


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
	    for (int i=0; i<ns; i++) {
		a[i][t] = a[i][t]*sc[t];
	    }
	}

	double lpd = 0; // pr(seq)
	for (int t=0; t<n; t++) {
	    lpd = lpd - Math.log(sc[t]);
	}

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
	    for (int i=0; i<ns; i++) {
		a[i][t] = a[i][t]*sc[t];
	    }
	}

	double lpdas = 0; // pr(seq & subtrellis)
	for (int t=0; t<n; t++) {
	    lpdas = lpdas - Math.log(sc[t]);
	}


	// // not needed --- prior probability of subtrellis

	//    a = new double[ns][n];

	//    for (int i=0; i<ns; i++) {
	//    if (st[i] == 1) a[i][0] = sc[0]*iprob[i];
	//    }
	//    for (int t=1; t<n; t++) {
	//    double sf = 0.0;
	//    for (int i=0; i<ns; i++) {
	//    if (st[i] == 1) {
	//    double score = 0;
	//    for (int k=0; k<ns; k++) {
	//    score = score + tprob[k][i]*a[k][t-1];
	//    }
	//    a[i][t] = score;
	//    sf = sf + score;
	//    }
	//    }
	//    sc[t] = 1/sf;
	//    for (int i=0; i<ns; i++) a[i][t] = a[i][t]*sc[t];
	//    }

	//    double lps = 0; // pr(subtrellis)
	//    for (int t=0; t<n; t++) lps = lps - Math.log(sc[t]);

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
		if (pp[j][i] > pp[map[i]][i]) map[i] = j;
	    }
	}
	// if (numclasses < ns) map = plaac.mapseq(map,classes);
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
		if (pp[j][i] > pp[map[i]][i]) map[i] = j;
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
		if (pp[j][i] > pp[map[i]][i]) map[i] = j;
	    }
	}
	// if (numclasses < ns) map = plaac.mapseq(map,classes);
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
		if (pp[j][i] > pp[map[i]][i]) map[i] = j;
	    }
	}
	// if (numclasses < ns) map = plaac.mapseq(map,classes);
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
		if (pp[j][i] > pp[map[i]][i]) map[i] = j;
	    }
	}
	// if (numclasses < ns) map = plaac.mapseq(map,classes);
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
	for (int i=0; i<pp.length; i++) {
	    for (int j=0; j<pp[0].length; j++) {
		npp[mask[i]][j]+=pp[i][j];
	    }
	}
	return npp;


    }

 
    // public int [][] simulate(int n, Random rg) {
    // 	int [] hseq = new int[n];
    // 	int [] oseq = new int[n];
    // 	double r = rg.nextDouble();
    // 	int dex = 0;
    // 	while (r < ciprob[dex]) dex++;
    // 	hseq[0]=dex;
    // 	for (int i=1; i<n; i++) {
    // 	    r = rg.nextDouble();
    // 	    dex = 0;
    // 	    while (r < ctprob[hseq[i-1]][dex]) dex++;
    // 	    hseq[i] = dex;
    // 	}
    // 	for (int i=0; i<n; i++) {
    // 	    r = rg.nextDouble();
    // 	    dex = 0;
    // 	    while (r < ceprob[hseq[i]][dex]) dex++;
    // 	    oseq[i] = dex;

    // 	}
    // 	int [][] smat = new int[2][n];
    // 	smat[0] = hseq;
    // 	smat[1] = oseq;
    // 	return smat;

    // }
   

    public void printme(String st) {
	int nc = postprob.length;

	for (int i=0; i<viterbipath.length; i++) {
	    System.out.print(st + "\t" + viterbipath[i] + "\t" + mappath[i]);
	    for (int j=0; j<nc; j++) System.out.print("\t" + (float) postprob[j][i]);
	    System.out.println();
	}
	System.out.print("# " + (float) etst + "\t" + (float) lpst);
	for (int i=0; i<nc; i++) System.out.print("\t" + subtrellis[i]);
	System.out.println();
	System.out.print("# " +( float) etst + "\t" + (float) lpst);
	for (int i=0; i<nc; i++) System.out.print("\t" + (float) plaac.sum(postprob[i]));
	System.out.println();
    }

    public void printme(String st, int [] seq) {
	int nc = postprob.length;

	for (int i=0; i<viterbipath.length; i++) {
	    System.out.print(st + "\t" + seq[i] + "\t" + viterbipath[i] + "\t" + mappath[i]);
	    for (int j=0; j<nc; j++) System.out.print("\t" + (float) postprob[j][i]);
	    System.out.println();
	}
	System.out.print("# " + (float) etst + "\t" + (float) lpst);
	for (int i=0; i<nc; i++) System.out.print("\t" + subtrellis[i]);
	System.out.println();
	System.out.print("# " + (float) etst + "\t" + (float) lpst);
	for (int i=0; i<nc; i++) System.out.print("\t" + (float) plaac.sum(postprob[i]));
	System.out.println();
    }


    // // allows partial labels known 
    // public double [][][] scoreseqs(int [][] seqs, int [][] labels, double [][] tprob, double [][] eprob, double [] iprob) {
    // 	// double [][] tprob;
    // 	// double [][] eprob;
    // 	// double [] iprob;
    // 	int ns = iprob.length;

    // 	int n = seqs.length;
    // 	double [][][] scores = new double[n][ns][];
    // 	for (int i=0; i<n; i++) {
    // 	    scores[i] = posteriors(seqs[i], labels[i], tprob, eprob, iprob);
    // 	}

    // 	for (int i=0; i<n; i++) {
    // 	    plaac.printmatrix(plaac.transpose(scores[i]));
    // 	    System.out.println("-------------------------------------------------------------");
    // 	}
    // 	return scores;

    // }

    // // allows partial labels known 
    // public double [][][] scoreseqs(int [][] seqs, int [][] labels) {
    // 	// double [][] tprob;
    // 	// double [][] eprob;
    // 	// double [] iprob;
    // 	int ns = iprob.length;

    // 	int n = seqs.length;
    // 	double [][][] scores = new double[n][ns][];
    // 	for (int i=0; i<n; i++) {
    // 	    scores[i] = posteriors(seqs[i], labels[i], tprob, eprob, iprob);
    // 	}

    // 	for (int i=0; i<n; i++) {
    // 	    plaac.printmatrix(plaac.transpose(scores[i]));
    // 	    System.out.println("-------------------------------------------------------------");
    // 	}
    // 	return scores;

    // }


    // dur = expected run length; ini = relative num of segments of each type
    public static double [][] autotransmat(double [] ini, double [] dur) {
	int n = ini.length;
	double [][] tm = new double[n][n];
	double rs = plaac.sum(ini);
	for (int i=0; i<n; i++) {
	    double lrs = rs - ini[i];
	    for (int j=0; j<n; j++) {
		tm[i][j] = (ini[j]/lrs)*(1.0/dur[i]);
	    }
	    tm[i][i] = 1.0 - 1.0/dur[i]; // fix diagonal
	}
	return tm;
    }


    // HMM representation suitable for graphing with graphviz/dot (http://www.graphviz.org/)
    // Some of the layout is tuned for 2 state HMMs, may need to change for others
    // if record = true, includes emission probabilities, otherwise just state transitions.
    // dot -Tpng filename > hmmimage.png


    public void dottify(String filename, boolean showemissionprobs) {
	
	int maxlab = 7; // truncate labels --- not currently used
	// pad shorter names with whitespace so symbol have same size and table columns have same width?
	
	try {
	    BufferedWriter out = new BufferedWriter(new FileWriter(filename));
	    String lab;
	    out.write("Digraph G {"); out.newLine();
            out.write("edge [fontname=Courier, fontsize=8, labelfontname=Courier,labelfontsize=8];");  out.newLine();
	    out.write(" node [fontname=Courier, fontsize=10]"); out.newLine();
	    out.write(" start [label=start, shape=circle, height=0.25, style=filled, color=grey, rank=source];");  out.newLine();
	    for (int i=0; i<ns; i++) {
		out.write("  n" + i + " [label=\"" + names[i] + "\", shape=circle, height=1.2];"); out.newLine();
		lab = String.format("%.3f", iprob[i]);
               	out.write("  start -> n" + i + " [label=\"" + lab + "\", color=gray];"); out.newLine();
	    }
	    for (int i=0; i<ns; i++) {
		for (int j=0; j<ns; j++) {
		    if (tprob[i][j] > 0) {
			lab = String.format("%.3f", tprob[i][j]);
		        if (i==j) { 
			    if (i % 2 == 0) {
				// invisible inner self-edge makes visible outer self-edge a little bigger
                                out.write("  n" + i + ":w -> n" + j + ":w [label=\"" + "\", color=gray, style=invis];"); out.newLine();
				out.write("  n" + i + ":w -> n" + j + ":w [label=\"" + lab + "\", color=gray];"); out.newLine();
			    } else {
				out.write("  n" + i + ":e -> n" + j + ":e [label=\"" + "\", color=gray, style=invis];"); out.newLine();
				out.write("  n" + i + ":e -> n" + j + ":e [label=\"" + lab + "\", color=gray];"); out.newLine();
			    }
			} else {
			    if (i < j) {
				// gives nicer layout, at least for ns=2: longish (but not shown) label makes nodes further apart, 
				// and takes the more direct route between nodes, displacing real transition edges to curve above and below
				out.write("  n" + i + " -> n" + j + " [label=\"" + "spacerlabel" + "\", color=gray, constraint=false, style=invis];"); out.newLine();
			    }
			    out.write("  n" + i + " -> n" + j + " [label=\"" + lab + "\", color=gray, constraint=false];"); out.newLine();
			}
		    }
		}
	    }

	    if (showemissionprobs) {
		//  hack: option to skip first and last output states, for X and *
		boolean skipfirstandlast = true;
		int starto = 0;
		int endo = no;
                if (skipfirstandlast) {starto++; endo--;}
		for (int i=0; i<ns; i++) {
		    out.write("rec" + i + " [shape=record, label=\"");
		    out.write("{ <fs> AA");
		    for (int j=starto; j<endo; j++) {
			lab = new String("" + plaac.aanames[j]);
			// if (lab.length() > maxlab) lab = lab.substring(0, maxlab-1);
			out.write("|" + lab);
		    }
		    out.write("}");
		    out.write("|{ <f" + i + "> " + "prob");
		    for (int j=starto; j<endo; j++) {
			// may need to increase to %.5f or more
			lab = String.format("%.4f", eprob[i][j]);
			out.write("|" + lab);
		    }
		    out.write("}");
		    out.write("\"];"); out.newLine();
		}
		for (int i=0; i<ns; i++) {
		    // out.write("  n" + i + " -> rec" + i + ":f" + i + " [style=dashed];"); out.newLine();
		    out.write("  n" + i + " -> rec" + i + " [style=dashed];"); out.newLine();
		}
	    }
	    out.write("}");
            out.newLine();
	    out.close();
	}
	catch (IOException e) {
	    System.out.println("## problem writing to dotfile");
	}
    }

}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// This class reads sequences one at a time from a fasta file. It's a bit clumsy,
// since one needs to call hasmorefastas() before nextname() before nextfasta().
// Could use peek or lookahead instead? Or at least keep track of which of these
// has been called most recently. 
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

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
	    System.out.println("# Couldn't open " + filename);
	    goodfile = false;
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
	if (!goodfile) return false; // some error message, too??
	if (ondeck) return true;
	try {
	    while ((line = in.readLine()) != null) {
		if (line.length() > 0 && line.charAt(0) == '>') {
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

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
// This class implements the "lowest probability subsequences" algorithm in
//   Harrison PM and Gerstein M (2003)
//   A method to assess compositional bias in biological sequences and its application to prion-like 
//   glutamine/asparagine-rich domains in eukaryotic proteomes. Genome Biol 4(6):R40
// It doesn't inlcude the pre-filtering step described in that paper for speedup, but can use lookup
// tables of binomial probabilities for speedup. This class is commented out, as it has not been 
// tested thoroughly or recently. It also was not clear to me why one should look for subsequences with 
// the lowest probability under binomial distributions --- why not subsequences with the lowest tail 
// probability under binomial distributions? E.g., compute probability of seeing at least k Qs out of
// n contiguous residues rather than the probabililty of seeing exactly k Qs out of n contiguous residues.
// As n gets very large, the probability of seeing any particular value of k becomes small, even if there 
// is no bias whatsoever.
/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

/*
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

    // direction of bias to look for: too many or too few
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

    // look-up table
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
	    luts[i] = logbinomiallut2(maxlen, ps[i]);
	}
    }

    static double binomial(int n, int k) {
	double mul = 1.0;
	int top = n;
	for (int i=1; i<=k; i++) {
	    mul = (mul*top)/k; 
	    top--;
	}
	return mul;
    }
     
    static double logbinomial(int n, int k, double p){
	// log C(n,k) p^k (1-p)^(n-k)
       	double lb = k*Math.log(p) + (n - k)*Math.log(1 - p);
       	double nn = 1.0*n;
       	for (int i=0; i<k; i++) lb = lb + Math.log((nn - i)/(k - i));
	return lb;
    }

    static double [][] logbinomiallut(int maxn, double p) {
       	double [][] lut = new double[maxn][maxn];
       	for (int i=0; i<maxn; i++) {
	    for (int j=0; j<=i; j++) {
		lut[i][j] = logbinomial(i, j, p);
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
		lut[n][k] = lut[n-1][k] + Math.log((n + 1.0)/(n + 1.0 - k)*q);
	    }
       	}
       	return lut;
    }


    public double [][] check(int [] seq) {
	double [][] bias = new double [masks.length][4];
	for (int i=0; i<masks.length; i++) {
	    bias[i] = lps(plaac.mapseq(seq,masks[i]), ps[i], luts[i], minlen, maxlen, tails[i]);
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
	    bias[i] =  lps(plaac.mapseq(seq,masks[i]), ps[i], luts[i], minlen, maxlen, tails[i]);
	    leng[i] = bias[i][1]-bias[i][0]+1;
	    if (leng[i]>longest) {longest = (int) leng[i]; lb = (int) bias[i][0]; ub = (int) bias[i][1];}
	}
	// merge these if they overlap by at least 15?
	for (int i=0; i<=2; i++) {
	    if (overlap(lb,ub,(int) bias[i][0], (int) bias[i][1]) > minol) {lb=Math.min(lb, (int) bias[i][0]); ub=Math.max(ub, (int) bias[i][1]);};
	}
	// bias[3][0] = lb;
	// bias[3][1] = ub;
	int [] subseq = plaac.submatrix(seq,lb,ub);

	for (int i=3; i< masks.length; i++) {
	    bias[i] =  lps(plaac.mapseq(subseq,masks[i]), ps[i], luts[i], subseq.length, maxlen, tails[i]);
	    bias[i][0] += lb; bias[i][1] += lb;
	    leng[i] = bias[i][1]-bias[i][0]+1;
	    // if (leng[i]>longest) {longest = (int) leng[i]; lb = (int) bias[i][0]; ub = (int) bias[i][1];}
	}

	// System.out.println(plaac.aa2string(subseq));
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
	    bias[i] = lps(plaac.mapseq(seq,masks[i]), ps[i],	luts[i], minlen, maxlen, tails[i]);
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
	int [] subseq = plaac.submatrix(seq,slb,sub);

	for (int i=3; i< masks.length; i++) { // uses entire sequence here!!!
	    bias[i] =  lps(plaac.mapseq(subseq,masks[i]), ps[i], luts[i], subseq.length,  maxlen, tails[i]);
	    bias[i][0] += slb; bias[i][1] += slb;
	    leng[i] = bias[i][1]-bias[i][0]+1;
	}

	boolean good = true;

	if (bias[3][3] >= threshs[3] || bias[4][3] >= threshs[4]) good = false;
	if (bias[5][3] >= threshs[5] && bias[6][3] >= threshs[6] && bias[7][3] >= threshs[7] ) good = false;

	if (good) {
	    // plaac.printmatrix(bias);
	    return true;
	}

	else {
	    int [] seq1 = plaac.submatrix(seq,0,lb-1);
	    int [] seq2 = plaac.submatrix(seq,ub+1,seq.length);
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
    		score = logbinomial(len, ct, p);
    		// System.out.println(i + " " + j + " " + len + " " + ct + " " + (float) score);
    		if (score < bestscore) {
    		    bestscore = score;
    		    besti = i; 
    		    bestj = j; 
    		    bestcount = ct;
    		}
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
    		if (score < bestscore) {
    		    bestscore = score; 
    		    besti = i; 
    		    bestj = j; 
    		    bestcount = ct;
    		}
    
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
    		    if (score < bestscore) {
    			bestscore = score; 
    			besti = i; 
    			bestj = j; 
    			bestcount = ct;
    		    }
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
    
    
    	int [] nzdex = plaac.sparsesupport(seq,1);
    
    	int m = nzdex.length;
    
    	int len;
    	int ct;
    
    	for (int i=0; i<m; i++) {
    	    for (int j=i; j<m; j++) {
    		len = nzdex[j]-nzdex[i]+1;
    		if (len >= maxlen) continue;
    		ct = j-i + 1;
    		score = lut[len][ct];
    		if (score < bestscore) {
    		    bestscore = score; 
    		    besti = nzdex[i];
    		    bestj = nzdex[j]; 
    		    bestcount = ct;
    		}
    	    }
    	}
    	sss[0] = 1.0*besti;
    	sss[1] = 1.0*bestj;
    	sss[2] = 1.0*bestcount;
    	sss[3] = bestscore;
    	return sss;
    }
    
}
*/


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
// This class implements FoldIndex for predicting protein disorder, following 
//   Prilusky J, Felder C, Zeev-Ben-Mordehai T, Rydberg EH, Man O, Beckmann JS, Silman I, Sussman JL (2005) 
//   FoldIndex: a simple tool to predict whether a given protein sequence is intrinsically unfolded. 
//   Bioinformatics 21:3435-3438.] 
// and PAPA for prion-forming propensity, following 
//  Toombs JA, McCarty BR, Ross ED. (2010)
//  Compositional determinants of prion formation in yeast. Mol Cell Biol 30:319-332.
// and  
//   Toombs JA, Petri M, Paul KR, Kan GY, Ben-Hur A, Ross ED (2012) 
//   De novo design of synthetic prion domains. Proc Natl Acad Sci USA 109:6519-6524.
// It also computed the smoothed profiles of log-likelohood ratios from PLAAC.
/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

class disorderreport {

    // double [] ccdef = {2.785, -1.0, -1.151}; // weights from FoldIndex; now passed as cc 
    double [] cc;    
    double [] hydro;
    double [] charge;
    double [] fi;  
    double [] plaacllr; 
    double [] papa; 
    double [] papax2;      // x2 are twice-smoothed versions
    double [] plaacllrx2;  // x2 are twice-smoothed versions  
    double [] fix2;        // x2 are twice-smoothed versions
    int [] aa;
    int ww1;
    int ww2;
    int ww3;
    // boolean reflect;
    double meancharge;
    double meanhydro;
    double meanfi;
    int n;
    int numseg;
    int numdisordered;
    int numdisorderedstrict;  // only windows completely contained in seq 
    int numdisorderedstrict2; // only windows completely contained in seq and of min length
    int [] startaa;
    int [] stopaa;
    int [] lenaa;
    double [] localmean;
    double [] localsd;
    double [] localhydro;
    double [] localcharge;
    int maxlen;		// for fi
    int minlen = 5;	// for fi
    int minlen2 = 40;	// for fi
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
    
    double [] plaacweights;   
    double [] papaweights;  
    
    double papascore;     // temp
    double papamaxprop;	  // maximum papa propensity
    double papamaxscore;  // tradeoff between prop and FI
    double papamaxdis;	  // corresponding FI
    double papamaxllr;	  // corresponding LLR
    double papamaxllr2;	  // corresponding LLR2
    int papamaxcenter;	  // index of max above
   

    disorderreport(int [] aa, int ww1, int ww2, int ww3, double [] cc, 
		   double [] plaacweights, double [] papaweights, boolean adjustprolines) {
    	this.ww1 = ww1;
    	this.ww2 = ww2;
	this.ww3 = ww3;
    	this.cc = cc;
    	this.aa = aa;
    	this.plaacweights = plaacweights;
    	this.papaweights = papaweights;
    	n = aa.length;
    
	double [] maa1 = plaac.mapseq(aa, plaac.aahydro2);
    	meanhydro = plaac.mean(maa1);
    	hydro = plaac.slidingaverage(maa1,ww1,true,false);
    	
        double [] maa2 = plaac.mapseq(aa,plaac.aacharge);
    	meancharge = plaac.mean(maa2);
	charge = plaac.slidingaverage(maa2,ww1,true,false);

	meanfi  = cc[2] + cc[1]*Math.abs(meancharge) + cc[0]*meanhydro;
	// note absolute value of local net charge here!
	fi = plaac.axpbypc(cc[0],hydro,cc[1],plaac.absval(charge),cc[2]);

	maa3 = plaac.mapseq(aa,plaacweights);
	plaacllr = plaac.slidingaverage(maa3,ww3,true,false);
	maa4 = plaac.mapseq(aa,papaweights);

	if (adjustprolines) {
	    papa = plaac.slidingaverage(maa4,ww2,true,false,13,aa);
	} else { 
	    papa = plaac.slidingaverage(maa4,ww2,true, false);
	}

	// Average of averages. Scores for first and last ww residues depend on how boundaries are handled.
        // Here we use the approach used in PAPA, for consistency: 
        // In first round of averages, smaller windows are used when sliding window would extend out-of-bounds;
	// In second round of averages, scores are weighted based on number of residues in each round-one window. 
	papax2 = plaac.slidingaverage(papa,ww2,false,true);
	plaacllrx2 = plaac.slidingaverage(plaacllr,ww3,false,true);
	fix2 = plaac.slidingaverage(fi,ww1,false,true);
	

        numdisordered = 0;
        for (int k=0; k<n; k++) {
	    if (fi[k]<0) { 
		numdisordered++;
	    }
        } 

        numdisorderedstrict = 0;
	int halfw = (ww1-1)/2; // or 0
	if (halfw > n/2) halfw = n/2;
	// System.out.println("halfww " + halfw);
	if (fi[halfw]<0) { 
	    numdisorderedstrict =  numdisorderedstrict + halfw + 1;
	}
	if (fi[n-halfw-1]<0) { 
	    numdisorderedstrict =  numdisorderedstrict + halfw + 1;
	}
	for (int k = halfw + 1; k < n-halfw - 1; k++) {
	    if (fi[k]<0) { 
		numdisorderedstrict++;
	    }
	} 

	// restrict to MAP parse?
	papamaxscore = -1.0/0; // -Infinity
        papamaxcenter = -1;	


	int papamode = 1; // 1: maximum PAPA score with FI < 0 
                          // 2: maximum PAPA score, even if FI > 0
                          // 3: tradeoff between PAPA score and FI: maximum of min(PAPA - 0.05, -FI) 
                          // 4: tradeoff between PAPA score, FI, and LLR: maximum of min(PAPA - 0.05, -FI, LLR) 
	
	if (papamode == 1) {
	    for (int k=(ww2-1)/2; k<n-(ww2-1)/2; k++) {
		papascore = papax2[k];
		if ((papascore > papamaxscore) & (fix2[k] < 0)) { /// what should score be if no disorder? -Inf? NA?
		    papamaxcenter = k;
		    papamaxscore = papascore; 	   
		} 
	    } 
	} else if (papamode == 2) {
	    for (int k=(ww2-1)/2; k<n-(ww2-1)/2; k++) {
		papascore = papax2[k];
		if (papascore > papamaxscore) { 
		    papamaxcenter = k;
		    papamaxscore = papascore;
		} 
	    } 
	} else if (papamode == 3) {
	    for (int k=(ww2-1)/2; k<n-(ww2-1)/2; k++) {
                // if both positive, find distance to decision boundary. weight dimensions differently?
		if ((fix2[k] < 0) && (papax2[k] > 0.05)) {
		    papascore = Math.min(-1*fix2[k], (papax2[k] - 0.05));
		} else {  
		    papascore = -1*Math.sqrt((Math.min(0,-1.0*fix2[k])*Math.min(0,-1.0*fix2[k]) 
					      + Math.min(0,papax2[k]-0.05)*Math.min(0,papax2[k]-0.05)));
		}
		if ((papascore > papamaxscore)) { 
		    papamaxcenter = k;
		    papamaxscore = papascore;	    
		} 
	    } 
	} else if (papamode == 4) { // improve this!
	    for (int k=(ww2-1)/2; k<n-(ww2-1)/2; k++) {
	        papascore = Math.min(-1.0*fix2[k]-0.0,papax2[k]-0.05); // weight dimensions differently?
                papascore = Math.min(papascore, plaacllrx2[k]-0.0);    // use singly-smoothed LLR instead?
                if ((papascore > papamaxscore)
		    || ((papascore == papamaxscore)                      
			&& ((papax2[k] > papax2[papamaxcenter]) || (fix2[k] < fix2[papamaxcenter]) 
			    || (plaacllrx2[k] > plaacllrx2[papamaxcenter])))) { 
		    papamaxcenter = k;
		    papamaxscore = papascore;
		} 
	    } 
	}

        // undefined if no valid center
	papamaxprop = 0.0/0.0; 
	papamaxdis = 0.0/0.0;  
        papamaxllr = 0.0/0.0;  
        papamaxllr2 = 0.0/0.0; 
        // papamaxscore = 0.0/0.0; // already set 

	if (papamaxcenter >= 0) {
	    papamaxprop = papax2[papamaxcenter];
	    papamaxdis = fix2[papamaxcenter];
	    papamaxllr2 = plaacllrx2[papamaxcenter];
	    papamaxllr = plaacllr[papamaxcenter];
	}

	// System.out.println(papamaxprop + "\t" + papamaxscore);
	/////////////////////////////

	hssr = plaac.hss2(maa3);
	if (hssr[1]-hssr[0]+1 < minlength || hssr[1]-hssr[0]+1 > maxlength) {
	    hssr2 = plaac.hss2(maa3, minlength, maxlength);
	} else {
	    hssr2 = hssr;
	}

	/////////////////////////////
	int i = 0 + halfw; // +1?
	startaa = new int[n/minlen];
	stopaa = new int[n/minlen];
	lenaa = new int[n/minlen];
	localmean = new double[n/minlen];
	localsd = new double[n/minlen];
	localhydro = new double[n/minlen];
	localcharge = new double[n/minlen];
	if (localhydro.length > 0) localhydro[0] = 0.5; //?
	numseg = 0;
	numdisorderedstrict2 = 0;
	while (i < n - halfw) {
	    if (fi[i] < 0) {
		double sc = fi[i];
		double lc = Math.abs(charge[i]);
		double lh = hydro[i];
		double scsc = fi[i]*fi[i];
		int startdex = i;
		i++;
		while (i < n - halfw && fi[i] <0) {
		    sc = sc + fi[i];
		    lh = lh + hydro[i];
		    lc = lc + Math.abs(charge[i]);
		    scsc = scsc + fi[i]*fi[i];
		    i++;
		}
		int stopdex = i - 1;
		int len = stopdex - startdex + 1;
		double msc = sc/len;
		double sdsc = Math.sqrt(scsc/len - msc*msc);
		if (startdex == halfw) startdex = 0;
		if (stopdex == n - halfw -1 ) stopdex = n-1;
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
	if (lenaa.length>0) maxlen = plaac.maxint(lenaa);
	maxlong = 1;
	bestdex = 0;
	for (i=0; i<numseg; i++) {
	    if (lenaa[i] >= minlen2 && localmean[i] < maxlong) {
		maxlong = localmean[i]; 
		bestdex = i;
	    }
	}
    }

   
    public void printme() {
	printme("");
    }

    // FoldIndex summary
    public void printme(String id){
	System.out.println("### protein length:	   " + n);
	System.out.println("### num disordered(1): " + numdisordered);
	System.out.println("### num disordered(2): " + numdisorderedstrict);
	System.out.println("### num disordered(3): " + numdisorderedstrict2);
	for (int i=0; i<numseg; i++)  {
	    System.out.println("### " + startaa[i] + "-" + stopaa[i] + ": " + (float) localmean[i] 
			       + " +- " + (float) localsd[i]);
	}
	System.out.println("### " + maxlen + " " + (float) localmean[bestdex] + " " + lenaa[bestdex]);
	System.out.println("# " + n + "\t" + (float) meancharge + "\t" + (float) meanhydro + "\t" 
			   + (float) meanfi + "\t" + (float) maxlong);
	System.out.println("# " + (float) hssr[0] + "\t" + (float) hssr[1] + "\t" + (float) hssr[2] + "\t"
			   + (float) hssr2[0] + "\t" + (float) hssr2[1] + "\t" + (float) hssr2[2]);
	for (int i=0; i<n; i++) {
	    System.out.println(id + "\t" + aa[i] + "\t" + (float) charge[i] + "	 \t"
			       + (float) hydro[i] + "  \t" + (float) fi[i] + "\t" + (float) plaacllr[i] 
			       + " \t"	 + (float) papa[i]
			       + "\t # ["+(i+1-ww1/2)+"-"+(i+1 + ww1/2)+"]");
	}
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

