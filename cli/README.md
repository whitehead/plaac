PLAAC command-line application
==============================

PLAAC is a Java command-line application for identifying proteins with Prion-Like Amino Acid Composition.

For most users, the recommended way to obtain PLAAC is to download a release package containing the precompiled `plaac.jar`. Building from
source is only necessary if you need an unreleased development version or want to modify the PLAAC source code.

Installing a released version
-----------------------------

Download the PLAAC command-line ZIP package from the [GitHub Releases page](https://github.com/whitehead/plaac/releases) and unpack it.

The package contains the precompiled `plaac.jar` together with the ancillary files needed by the command-line application.

> PLAAC requires Java Runtime Environment 11 or later. Check the
> installed Java version with:

    java -version

> If java is not installed, it can be downloaded from http://www.java.com.  

> **Note as of PLAAC 1.1.0, Java 11 is the minimum supported version of Java**

To display the PLAAC command-line options, run:

    java -jar plaac.jar

The release package is the preferred way to install a stable version of PLAAC. It does not require a Java compiler or a Git checkout.

Building from source
--------------------

> The source build requires a Java Development Kit (JDK) capable of
> compiling Java 11 bytecode. The resulting `plaac.jar` is compatible
> with Java 11 and later.

To build the command-line application from the Git repository, first
clone the repository and change into the `cli` directory:

    git clone https://github.com/whitehead/plaac.git
    cd plaac/cli

Source code and supporting files for the algorithm are located in the ```src``` subdirectory, R code for visualization is located in ```R``` subdirectory.  

Now run the provided shell script:

    ./build_plaac.sh
	
The script compiles the Java source and produces:

    target/plaac.jar

It also generates files used by the web application and documentation.

The build uses the Git repository to determine the PLAAC version. Release
versions are identified by Git tags; development builds include the
nearest release version together with development and Git commit
information.

Usage
-----

### Basic usage

To run PLAAC and see the command-line options (using either the pre-compiled package `.jar` file, or the `.jar` that you compiled from
source), open a terminal or command prompt and for complete usage options, type:

    java -jar plaac.jar

with no additional arguments.

### Use case 1:

Compute per-sequence scores for every sequence in a protein fasta file:

    java -jar plaac.jar -i input.fa > output.txt

This will score every sequence in the input protein fasta file ```input.fa```, and write the results 
(one line per protein) as a table to output.txt. The header lines of ```output.txt``` start with the 
symbol ```#``` and include descriptions of the columns.

### Use-case 2: 

Compute per-residue scores, which are suitable for plotting, 
for [a subset of] sequences in a protein fasta file:

    java -jar plaac.jar -i input.fa -p print_list.txt > plotdata.txt

where ```print_lists.txt``` gives the names of the subset of sequences in ```input.fa``` to print per-residue 
scores for; or 

    java -jar plaac.jar -i input.fa -p all > plotdata.txt

which will give per-residue scores for all the sequences in ```input.fa``` - this is not recommended 
if ```input.fa``` is an entire proteome, as the plotting routines are not optimized for very large files.

    java -jar plaac.jar -b background.fa -i input.fa -p all > plotdata.txt

### Background frequencies

The main complication to both use-cases above is specifying the background amino-acid frequencies 
to use for plaac. As the prion-like AA frequencies were derived from *S. cerevisiae*, one option 
(the default) is to use background frequencies from *S. cerevisiae* as well (i.e., treat the per-AA likelihood-ratios 
are species independent).  This is done with the command-line option ```-a 1``` (the default).
The other extreme is to use use the background frequencies of the species being scores 
(i.e., treat the prion-like AA frequencies as species-independent). 
This is done with the command-line option ```-a 0```. One can also linearly interpolate between these 
two extremes with any value of a between 0 and 1, e.g. use ```-a 0.5``` for an average of the background 
frequencies of *S. cerevisiae* and the species being scored. 

The background AA frequencies for *S. cerevisiae* are built in to ```plaac.java```. For any other species, 
```plaac.java``` needs to know what frequencies to use. If a < 1, the default is to compute the background
AA frequencies for the input species from file given by ```-i input.fa```, and then interpolate 
between these with the background frequencies from *S. cerevisiae* (degree of interpolation given by ```-a```).

But if ```-i input.fa``` is not an entire proteome, this may not be a good estimate of background frequencies;  
in particular if ```-i input.fa``` is pre-selected to consist of proteins with prion-like domains (e.g. candidates 
for plotting) then this is a bad idea, as the background AA frequency can be strongly skewed. 
To avoid this problem, instead one can include the entire proteome in ```-i input.fa``` and just give the 
IDs of the proteins to plot in print\_list.txt (one per line, exactly matching sequence names in the input.fa); 
or one can use just a subset of proteins in ```-i input.fa``` but tell plaac to compute the background frequencies 
for the species from a different fasta file (with ```-b background.fa```) or from a table (with ```-B bg_freqs.txt```). 
In both cases the frequencies will still be interpolated with the *S. cerevisiae* background AA frequencies 
(unless ```a``` = 0).

So, **CAUTION**: 
Unless you are using ```a```=1 (just *S. cerevisiae* background frequencies), if ```input.fa``` consists of just a small 
number of sequences it is probably *not* a good idea to use a command like:

    java -jar plaac.jar -i input.fa -a 0.5 > output.txt

to score each of the sequences at the protein level, or  

    java -jar plaac.jar -i input.fa -a 0.5 -p all > plotdata.txt

to score each of these sequences at the residue level.

Instead, you could use:

    java -jar plaac.jar -b background.fa -i input.fa -a 0.5 > output.txt

or 

    java -jar plaac.jar -b background.fa -i input.fa -a 0.5 -p all > plotdata.txt

or include the whole proteome as the input with ```-i``` and use ```print_list.txt``` to give names of sequences to plot. 

Plotting results
----------------

> The plotting routines are written in R, and to call them from the
> command-line requires `R` and the front-end ```Rscript```, which is
> included with many installations of R.  To check whether it is
> already installed, open a terminal or command prompt and type:

    Rscript --version
   
> If `Rscript` is not installed, it (along with `R`) can be downloaded from http://www.r-project.org

For usage options of the plotting routines, change into the ```R``` subdirectory of ```cli``` (i.e. ```cd R```) and type:

    Rscript plaac_plot.r

with no additional arguments. 

The basic usage is:

    Rscript plaac_plot.r  plotdata.txt figname[.pdf|.png] [-d]

where ```plotdata.txt``` is the output from plaac with the ```-p``` option above, 
and ```figname.pdf``` (in which case plots will be one sequence per page in a pdf)
or with ```figname.png``` (in which plots will be one sequence per png image).
If the optional argument ```-d``` is omitted, plots show sliding averages of per-residue scores; 
if ```-d``` is included, plots show sliding averages of sliding averages (```-d``` for doubly-smoothed). 

Individual tracks can be disabled or enabled by editing the arguments to ```plot_seqs()``` 
in the file ```plaac_plot.r```; the lower-level plotting functions are defined in the file ```plac_plot_util.r```, 
and can be called directly from within R if desired.

Notes
------

If FASTA sequences contain characters other the 20 AAs, or internal stop codons (```*```), the log-likelihood for 
those positions is set to zero, the charge is set to zero, and the shifted-and-scaled hydropathy score ```(1/9)*(KW+0.5)``` is set to zero. Results may not agree with ```PAPA``` or ```FoldIndex``` in these cases.

For sequences of length less than defined window sizes (```w``` or ```W``` or core length ```c```) some of the results may be automatically 
set to NA (undefined); others are scored based on the maximum attainable window size (e.g. for maximum number of Q+N in window of 
length 80 for the MW score). The idea is to allow shorter windows only if the scores take the form 
of sums rather than averages or per-residue scores, and if all residues have non-negative scores, 
so that scores are monotonic with window length. If other behavior is desired, check the protein length column and filter results as needed. 

The file [```src/scer_fg_28.fasta```](https://github.com/whitehead/PLAAC/blob/master/cli/src/scer_fg_28.fasta) contains the names of the yeast proteins and their prion-like domains (PrLDs) used to compute the foreground frequencies for the algorithm.

