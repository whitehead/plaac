#!/usr/bin/env python
# Author: Alex Lancaster
# Whitehead Institute for Biomedical Research, 2013
#
# Downloads proteome FASTA files from Uniprot
# builds background frequency data from plaac.jar
# and metadata and optionally copies to web area
# Note this assumes that build_plaac.sh has already
# been run

import urllib2
import re
import itertools
import datetime
import gzip
import urllib
import os, os.path
import StringIO
import string
import shutil


UNIPROT_BASEURL="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/proteomes/"
BACKGROUND_BUILD_DIR = "bg_freqs_build"
BACKGROUND_WEB_DIR = "../web/bg_freqs"
METADATA_FILENAME = "background_metadata.txt"

msg = 'Copy final background frequency files and metadata to' + BACKGROUND_WEB_DIR + '?'
copy_files_flag = True if raw_input("%s (y/N) " % msg).lower() == 'y' else False


# create output directory if it doesn't already exist
if not os.path.isdir(BACKGROUND_BUILD_DIR):
    os.mkdir(BACKGROUND_BUILD_DIR)

# get current date
download_date = str(datetime.date.today())

# parse the release notes
response = urllib2.urlopen(UNIPROT_BASEURL + "relnotes.txt")

release = re.compile("Release ([0-9]+_[0-9]+)")
species = re.compile("([a-zA-Z]+ [a-zA-Z]+) \(([A-Z]+)\)")

# build list of proteomes to download 
proteomes = {}
for line in response.readlines():
    rel = re.match(release, line)
    spec = re.match(species, line)
    if rel: uniprot_version = rel.group(1)
    if spec:
        species_name, id_name = spec.group(1, 2)
        proteomes[id_name] = species_name

# write metadata header
metadata_path = os.path.join(BACKGROUND_BUILD_DIR, METADATA_FILENAME)
metadata = open(metadata_path, 'w')
header = string.join(["uniprot_version", "id_name", "species_name", "fasta_url", "download_date", "aa_freq_filename", "\n"], "\t")
metadata.write(header)

for id_name, species_name in proteomes.iteritems():

    # generate filenames
    fasta_filename = id_name + '.fasta'
    fasta_filename_gzip = fasta_filename + '.gz'
    fasta_url = UNIPROT_BASEURL + fasta_filename_gzip
    aa_freq_filename = "bg_freqs_" + id_name + '.txt'

    # write metadata line 
    output_species = string.join([uniprot_version, id_name, species_name, fasta_url, download_date, aa_freq_filename, "\n"], "\t")
    metadata.write(output_species)
    
    # downloading and uncompress the file
    fasta_outpath = os.path.join(BACKGROUND_BUILD_DIR, fasta_filename)
    aa_freq_outpath = os.path.join(BACKGROUND_BUILD_DIR, aa_freq_filename)
    
    if not os.path.isfile(fasta_outpath):
        response = urllib2.urlopen(fasta_url)
        compressedFile = StringIO.StringIO(response.read())
        decompressedFile = gzip.GzipFile(fileobj=compressedFile)
        with open(fasta_outpath, 'w') as fasta_out:
            fasta_out.write(decompressedFile.read())

    cmd_line = "java -jar target/plaac.jar -b " + fasta_outpath + " > " + aa_freq_outpath
    print cmd_line
    os.system(cmd_line)

    # copy to web area
    if copy_files_flag:
        shutil.copyfile(aa_freq_outpath,
                        os.path.join(BACKGROUND_WEB_DIR, aa_freq_filename))

# close metadata file
metadata.close()
if copy_files_flag:
    shutil.copyfile(metadata_path,
                    os.path.join(BACKGROUND_WEB_DIR, METADATA_FILENAME))


    
