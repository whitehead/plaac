PLAAC
=====

PLAAC (Prion-Like Amino Acid Composition) searches protein sequences
to identify candidate prion subsequences using a hidden-Markov model
(HMM) algorithm.  The PLAAC website is located at: 

 http://plaac.wi.mit.edu/

This README.md file describes how to download, install and run PLAAC, as
well as developer instructions on how to compile from source and set
up the web framework.  PLLAC is released under the open-source MIT
license.

Citing PLAAC
------------

To cite use of either the PLAAC website or underlying code, in any
publication, please cite Lancaster et al. (2014), which describes
additional extensions and implementation details of the algorithm.
Details of the PLAAC algorithm are described in Alberti et al. (2009)
and also at http://plaac.wi.mit.edu/details.


* Alberti, S., Halfmann, R., King, O., Kapila, A., Lindquist, S. (2009) 
[A systematic survey identifies prions and illuminates sequence
features of prionogenic proteins](http://www.sciencedirect.com/science/article/pii/S0092867409002669).
*Cell* 137, 146–58.

* Lancaster, A.K., Nutter-Upham, A., Lindquist, King, O.D. (2014)
[PLAAC: a web and command-line application to identify proteins with
Prion-Like Amino Acid Composition](http://bioinformatics.oxfordjournals.org/content/early/2014/05/13/bioinformatics.btu310.abstract) *Bioinformatics* doi:10.1093/bioinformatics/btu310


Installation
------------

### Command-line application (cli)

#### Released package

For most users, the recommended way to install PLAAC is to download the latest release package from the [GitHub Releases page](https://github.com/whitehead/plaac/releases).

Download the PLAAC command-line ZIP package and unpack it. The package contains the precompiled `plaac.jar` and the ancillary files needed to
run PLAAC.

> PLAAC requires Java 11 or later. To check whether Java is installed,
> open a terminal or command prompt and type: `java -version`.  If
> Java is not installed, it can be downloaded from
> http://www.java.com/.

The release package contains a precompiled version of PLAAC, so compiling the Java source code is not necessary for normal use.

To display the available command-line options, run:

    java -jar plaac.jar

#### Building from source

Building PLAAC from source is primarily intended for developers and
users who need to work with unreleased code or modify the PLAAC source.

The command-line application can be built from a Git checkout using:

    cd cli
    ./build_plaac.sh

This will build a ```plaac.jar``` file , as well as column outputs for the website in the ```_plaac_headers.haml``` both in the ```cli/target``` subdirectory.  

PLAAC versions are derived from Git tags. Release builds use the Git/GitHub
tag as the source of the version number, while development builds include
information identifying the development version and Git commit.

More  detailed instructions for the use of the ```plaac``` cli are found in [```cli/README.md```](https://github.com/whitehead/plaac/blob/master/cli/README.md).

### Web-application

The web application is intended primarily for developers maintaining
the PLAAC website.


The web application uses the same compiled Java command-line application
as the standalone CLI. To build the web application from source, first
follow the CLI build instructions above.

The generated files can then be installed into the web application as
described in `web/README.md`.

The web application also requires the R plotting scripts. See
`web/README.md` for the remaining installation and development
instructions.

The remaining installation steps are detailed in the [```web/README.md```](https://github.com/whitehead/plaac/blob/master/web/README.md) (note that ```Rscript``` 
should be installed).
