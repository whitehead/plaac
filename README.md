PLAAC
=====

PLAAC (Prion-Like Amino Acid Composition) searches protein sequences
to identify candidate prion subsequences using a hidden-Markov model
(HMM) algorithm.  The PLAAC website is located at: 

 http://plaac.wi.mit.edu/

This README file is for developers and others interested in accessing
the Java code and the web framework, which is released under the
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


Installation
------------

### Command-line application (cli)

    cd cli
    ./build_plaac.sh

This will build a ```plaac.jar``` file , as well as column outputs for the website in the ```_plaac_headers.haml``` both in the ```cli/target``` subdirectory.  More  detailed instructions for the use of the ```plaac``` cli are found in [```cli/README.md```](https://github.com/whitehead/plaac/blob/master/cli/README.md).


### Web-application

Perform the same steps as for the cli, above, then:

    cp target/plaac.jar ../web/bin/
    cp R/plaac_plot.r R/plaac_plot_util.r ../web/bin/

If the column outputs from cli version have changed (i.e. in ```cli/src/plaac.java```), you will also need to update the HAML code for the website:

    cp target/_plaac_headers.haml ../web/views/

The remaining installation steps are detailed in the [```web/README.md```](https://github.com/whitehead/plaac/blob/master/web/README.md) (note that ```Rscript``` 
should be installed).
