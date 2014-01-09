Source code and supporting files for the algorithm are located in the
```src``` subdirectory, R code for visualization is located in ```R```
subdirectory.  Compilation can be performed using the provided shell
script:

```$ build_plaac.sh```

This results in a jar file ```plaac.jar``` located in ```target```.  
For usage details, this jar file can be run using:

```$ java -jar target/plaac.jar```

The file ```src/scer_fg_28.fasta``` contains the names of the yeast
proteins and their prion-like domains (PrLDs) used to compute the
foreground frequencies for the algorithm.
