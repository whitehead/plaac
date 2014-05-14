Quick start from github
------------------------


Source code and supporting files for the algorithm are located in the ```src``` subdirectory, R code for visualization is located in ```R``` subdirectory.  Compilation can be performed using the provided shell script:

    build_plaac.sh

This results in a jar file ```plaac.jar``` located in ```target```.  For usage details, this jar file can be run by changing to the ```target``` subdirectory ```cd target``:

    java -jar plaac.jar

Detailed installation
---------------------

PLAAC is written in java and requires that the Java Runtime Environment is installed. To check whether it is already installed, open a terminal or command prompt and type:

    java -version

If java is not installed, it can be downloaded from http://www.java.com

The plotting routines are written in R, and to call them from the command-line requires the R scripting front-end ```Rscript```, which is included with many installations of R.  To check whether it is already installed, open a terminal or command prompt and type:

    Rscript --version
   
If Rscript is not installed, it (along with R) can be downloaded from http://www.r-project.org


A precompiled version of the java program PLAAC is available for download in ```web/bin/plaac.jar```.  

To run PLAAC, download this file, open a terminal or command prompt and for complete usage options, type:

    java -jar plaac.jar

with no additional arguments.


Usage
-----

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
