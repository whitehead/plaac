#!/usr/bin/env python
# Author: Alex Lancaster
# Whitehead Institute for Biomedical Research, 2013
# Modified by Siddhant Sharma (18-10-25) - ISTA Austria
#
# Build documentation from print headers output of plaac.jar
# Note this assumes that build_plaac.sh has already
# been run.

import os, sys, html

with open(sys.argv[1], 'r') as doc:
    for line in doc:
        if line.startswith("## "):
            col, desc = line.split(":", 1)  # only take the first split
            col = col[3:]
            desc = desc.strip()

            print ("%%li\n  %%strong %s\n  %s" % (col, html.escape(desc)))
