#!/usr/bin/env python
# Author: Alex Lancaster
# Whitehead Institute for Biomedical Research, 2013
#
# Build documentation from print headers output of plaac.jar
# Note this assumes that build_plaac.sh has already
# been run.

import os, sys, cgi

with open(sys.argv[1], 'r') as doc:
    for line in doc:
        if line.startswith("## "):
            col, desc = line.split(":", 1)  # only take the first split
            col = col[3:]
            desc = desc.strip()

            print "%%li\n  %%strong %s\n  %s" % (col, cgi.escape(desc))
