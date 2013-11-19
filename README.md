PLAAC
=====

Prion-Like Amino Acid Composition


Installation
============

Command-line application (cli)
------------------------------

```
cd cli
./build_plaac.sh
```

This will build a jar file in the ```cli/target``` subdirectory.  Instructions for use of the
```plaac``` cli are found in ```cli/README```

Web-application
---------------

Perform the same steps as for the cli, above, then:

```
cp target/plaac.jar ../web/bin/
cp R/plaac_plot.r R/plaac_plot_util.r ../web/bin/
```

The remaining installation steps are detailed in the ```web/README``` (note that ```Rscript``` 
should be installed).
