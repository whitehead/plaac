PLAAC
=====

PLAAC (Prion-Like Amino Acid Composition) searches protein sequences
to identify candidate prion subsequences using a hidden-Markov model
(HMM) algorithm.  The PLAAC website is located at: 

 http://plaac.wi.mit.edu/

This README.md file describes how to download, install and run
standalone PLAAC, as well as developer instructions on how to compile
from source and set up the web framework.  PLAAC is released under the
open-source MIT license.

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


Command-line application (cli)
------------------------------

### Installing a released package

Released versions of the PLAAC command-line application are available from the
[GitHub Releases page](https://github.com/whitehead/plaac/releases).

Each release provides two packages:

- `plaac-v<VERSION>.zip` — a self-contained command-line distribution containing
  the precompiled `plaac.jar`, example input, R plotting scripts, and documentation.
- `plaac-<VERSION>-py3-none-any.whl` — a Python package that provides a command-line
  interface to the same PLAAC Java application.

For most users who want to run PLAAC from the command line, the ZIP package is
the simplest option.

#### Quick start: command-line ZIP package

1. Download the `plaac-v<VERSION>.zip` file from the
   [GitHub Releases page](https://github.com/whitehead/plaac/releases).

2. Unpack the ZIP file:

       unzip plaac-v<VERSION>.zip

   This creates a directory named `plaac-v<VERSION>` containing the PLAAC
   application and supporting files.

3. Change into the unpacked directory:

       cd plaac-v<VERSION>

4. Check that Java 11 or later is installed:

       java -version

   If Java is not installed, it can be downloaded from
   [java.com](http://www.java.com/).

5. Run PLAAC using the example FASTA file included in the distribution:

       java -jar plaac.jar -i four_classic_prions.fasta

   This runs the command-line application using the included example input
   and writes the results to standard output.

To display the available command-line options, run:

       java -jar plaac.jar

The release package contains a precompiled version of PLAAC, so compiling the
Java source code is not necessary for normal use.

More detailed instructions for the use of the `plaac` cli are found in
[`cli/README.md`](cli/README.md).

#### Quick start: Python wheel

The Python wheel provides a Python package and command-line interface for PLAAC.

1. Download the `plaac-<VERSION>-py3-none-any.whl` file from the
   [GitHub Releases page](https://github.com/whitehead/plaac/releases).

2. Install the downloaded wheel:

       python -m pip install plaac-<VERSION>-py3-none-any.whl

   You can also provide the path to the downloaded file if it is in another
   directory.

3. Use the PLAAC API wrapper as described in the
   [Python wrapper documentation](cli/python/README.md).

The Python package also requires Java 11 or later because the PLAAC application
includes the Java implementation.

### Building from source

Building PLAAC from source is primarily intended for developers and
users who need to work with unreleased code or modify the PLAAC
source. See [`cli/README.md`](cli/README.md#building-from-source) for
step-by-step details.

> **Note**: PLAAC versions are derived from Git tags. Release builds use the
> Git/GitHub tag as the source of the version number, while
> development builds include information identifying the development
> version and Git commit.

Web-application
----------------

The web application is intended primarily for developers maintaining
the PLAAC website.

The web application uses the same compiled Java command-line application
as the standalone CLI. To build the web application from source, first
follow the CLI source build instructions above to obtain a `plaac.jar` file.

The generated `plaac.jar` can then be installed into the web
application.  The application also requires the R plotting scripts (in
particular `Rscript`). See [`web/README.md`](web/README.md) for the
remaining installation and development instructions.

