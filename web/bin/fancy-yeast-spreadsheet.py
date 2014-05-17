#!/usr/bin/env python
# Author: Alex Lancaster <alexl@wi.mit.edu>
# Whitehead Institute for Biomedical Research
# Date: 2014-05-17

# Usage:
# fancy-yeast-spreadsheet.py inputfile-from-plaac [-|numhits] output.xls [gff3-file] [list-of-genes-to-include.txt] 

# FIXME/TODO:
# - add original ranking from Alberti & Halfmann
# - add experimental status (if there is one) from Table S2 in A&H
# - get human homolog gene name
# - get mouse/Drosophila?
# - add category (?) of gene

import os, sys, csv, string, math
import pybedtools

import urllib
from xlwt import Workbook, Style, Formula, easyxf

debug = False

def get_width(num_characters):
    return int((1+num_characters) * 256)

def move_to_back(li, item):
    try:
        li.remove(item)
    except ValueError:
        print item, li
    li.append(item)
    return li

def set_col_width(colname, header, ws, width):
    header_idx = header.index(colname)
    ws.col(header_idx).width = get_width(width)

# import csv file
infilename = sys.argv[1]

# by default print out all hits
allhits = True
include = False

if len(sys.argv) > 2:
    if sys.argv[2] == '-':
        allhits = True
    else:
        numhits =  int(sys.argv[2])
        allhits = False

    if len(sys.argv) > 3:
        outputfilename = sys.argv[3]

        if len(sys.argv) > 4:

            # use GFF3 file to get data for spreadsheet
            gff_filename = sys.argv[4]
            # sgd_gff = pybedtools.BedTool(gff_filename).remove_invalid().saveas()
            sgd_gff = pybedtools.BedTool(gff_filename).remove_invalid().saveas()

            # filter out only the ORFs ('gene'  in GFF), also skip the 2-micron plasmid
            # also check for transposable_element_gene (needed for YAR009C)
            all_orfs = sgd_gff.filter(lambda f: (f[2] == 'gene' or f[2] == 'transposable_element_gene') and f.chrom != '2-micron')

            orf_map = {}
            common_map = {}
            for orf in all_orfs:
                systematic = orf.attrs['ID']
                orf_map[systematic] = orf.attrs
                try:
                    gene = orf.attrs['gene']
                    common_map[gene] = systematic
                except KeyError:
                    pass
                # also map back to itself
                common_map[systematic] = systematic

            if len(sys.argv) > 5:
                filterfile = open(sys.argv[5], 'r')
                genes = [li.rstrip().upper() for li in filterfile.readlines()]
                include = []

                candidate_list = open("%s-candidate-list" % outputfilename, 'w')
                for gene in genes:
                    # convert all genes to systematic names
                    systematic = common_map[gene]
                    include.append(systematic)
                    candidate_list.write("%s\t%s\n" % (systematic, gene))
                candidate_list.close()


                

tsvfile = open(infilename, 'r')
while True:
    # store file pos, just before reading current line
    filepos = tsvfile.tell()
    line = tsvfile.readline()

    # skip comment lines
    if line[0] != '#':
        # found the true header, keep to preserve order
        header = line.rstrip().split('\t')
        break

# wind file back
tsvfile.seek(filepos)
# create reader
reader = csv.DictReader(tsvfile, delimiter='\t')
# save input
savedinput = []
for inputrow in reader:
    savedinput.append(inputrow)

# sort entries by COREscore, then by LLR
# need to take care of "NaN"
sortedinput = sorted(savedinput, key=lambda d: (float('inf') if math.isnan(float(d['COREscore'])) else -float(d['COREscore']),
                                                float('inf') if math.isnan(float(d['LLR'])) else -float(d['LLR'])))
for order, inputrow in enumerate(sortedinput):
    inputrow["rank_by_corescore"] = order + 1

def write_sheet(sortedinput, sheetnum, header):

    ws = wb.add_sheet('prion-scores-%d' % sheetnum)

    header_style = easyxf('font: bold true, color black;')
    url_style = easyxf('font: underline single, color blue;')
    corescore_pos_style = easyxf('font: bold true; pattern: pattern solid, fore_color light_orange;')
    corescore_zero_style =  easyxf('pattern: pattern solid_fill, fore_color light_orange;')
    llk_pos_style = easyxf('font: bold true; pattern: pattern solid, fore_color light_yellow;')
    llk_zero_style =  easyxf('pattern: pattern solid, fore_color light_yellow;')
    rossprd_high_style = easyxf('font: bold true; pattern: pattern solid, fore_color light_turquoise;')
    rossprd_low_style =  easyxf('pattern: pattern solid, fore_color light_turquoise;')
    shrink_col_style = easyxf('pattern: pattern solid, fore_color gray25; alignment: horizontal center, vertical center, wrap on;')
    row_style = easyxf('alignment: horizontal center, vertical center;')
    #default_style = easyxf('')


    # append new columns to header
    start_header = header[0:1]
    rest = header[2:]
    header = []
    header.append("rank_by_corescore")
    header.extend(start_header)
    header.append("gene_name")
    header.append("gene_aliases")
    header.append("orf_status")
    header.append("description")

    # rearrange columns in rest: push lesser-used cols to end of list
    for item in ["PRDstart", "PRDend", "PRDlen", "PROTlen", "COREaa", "STARTaa", "ENDaa", "PRDaa", "PRDscore", \
                     "MWstart", "MWend", "MWlen", "HMMall", "HMMvit", "FInumaa", "FImeanhydro", "FImeancharge", "FImeancombo", "FImaxrun"]:
        move_to_back(rest, item)

    # now add them back
    header.extend(rest)

    # initialize row number
    row = 0
    # write header
    for col in range(len(header)):
        ws.row(0).write(col, header[col], header_style)

    # parse fields and generate output Excel file

    for inputrow in sortedinput:
        orf_id = inputrow['SEQid']

        # if we are filtering, only print out specified hits
        if include:
            if orf_id in include:
                if debug: print "include:", orf_id, "rank:", inputrow["rank_by_corescore"]
            else:
                #print "skip:", orf_id, "rank:", inputrow["rank_by_corescore"]
                continue
        else:
            if debug: print orf_id, "rank:", inputrow["rank_by_corescore"]

        # loop through input cols
        for col in range(len(header)):

            # get key
            key = header[col]

            # note: offset row by one

            if key == "rank_by_corescore":
                # use the current row as the rank
                value = inputrow[key]
                if debug: print orf_id, value
                ws.row(row+1).write(col, value)
            # hyperlink ORF name to SGD
            elif key == "SEQid":
                try:
                    value = inputrow[key]
                    formula = 'HYPERLINK("http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=%s";"%s")' % (value, value)
                    link = Formula(formula)
                    ws.row(row+1).write(col, link, url_style)
                except:
                    ws.row(row+1).write(col, value)
            # get human readable name from GFF3 file
            elif key == "gene_name":
                try:
                    common_name = orf_map[orf_id]['gene']
                except (NameError, KeyError) as e:
                    # default back to using the ORF id
                    common_name = orf_id
                ws.row(row+1).write(col, common_name)
            elif key == "gene_aliases":
                try:
                    aliases = orf_map[orf_id]['Alias']
                except (NameError, KeyError) as e:
                    # if no aliases leave blank
                    aliases = ""
                if aliases == common_name:
                    aliases = ""
                if aliases.startswith(common_name):
                    aliases = string.replace(aliases, common_name+',', '')
                # also remove URL quoting, if it exists
                ws.row(row+1).write(col, urllib.unquote(aliases))
            elif key == "orf_status":
                try:
                    orf_status = orf_map[orf_id]['orf_classification']
                except (NameError, KeyError) as e:
                    orf_status = "NA"
                ws.row(row+1).write(col, orf_status)
            elif key == "description":
                try:
                    display = orf_map[orf_id]['display']
                except (NameError, KeyError) as e:
                    display = orf_id
                ws.row(row+1).write(col, urllib.unquote(display), shrink_col_style)
            elif key == "COREscore":
                value = inputrow[key]
                if float(value) > 0:
                    ws.row(row+1).write(col, value, corescore_pos_style)
                else:
                    ws.row(row+1).write(col, value, corescore_zero_style)
            elif key == "LLR":
                value = inputrow[key]
                if float(value) > 0:
                    ws.row(row+1).write(col, value, llk_pos_style)
                else:
                    ws.row(row+1).write(col, value, llk_zero_style)
            elif key == "PAPAprd":
                value = inputrow[key]
                if float(value) > 0.05:
                    ws.row(row+1).write(col, value, rossprd_high_style)
                else:
                    ws.row(row+1).write(col, value, rossprd_low_style)
            elif key == "PAPAfi":
                value = inputrow[key]
                if float(value) < 0.0:
                    ws.row(row+1).write(col, value, rossprd_high_style)
                else:
                    ws.row(row+1).write(col, value, rossprd_low_style)
            else:
                # get value at key
                value = inputrow[key]
                ws.row(row+1).write(col, value)

        ws.row(row+1).set_style(row_style)
        ws.row(row+1).height_mismatch = 1
        ws.row(row+1).height = 600

        row += 1
        if not(allhits):
            if row >= numhits:
                break

    # set column widths after everything generated
    set_col_width("description", header, ws, 25)

    set_col_width("LLR", header, ws, 5)
    set_col_width("LLRstart", header, ws, 7)
    set_col_width("LLRend", header, ws, 7)
    set_col_width("LLRlen", header, ws, 6)
    set_col_width("NLLR", header, ws, 6)

    set_col_width("COREscore", header, ws, 8)
    set_col_width("COREstart", header, ws, 8)
    set_col_width("CORElen", header, ws, 8)
    set_col_width("COREend", header, ws, 6)

    for colname in ["PRDstart", "PRDend", "PRDend", "PRDlen", "PROTlen"]:
        set_col_width(colname, header, ws, 8)

# create Excel writer
wb = Workbook()

# each sheet can have a maximum of 65536 rows, so need to split each
# sheet

total_rows = len(sortedinput)
max_rows_per_sheet = 65536 - 1

sheetnum = 0
for row_start in xrange(0, total_rows, max_rows_per_sheet):
    start = row_start
    end = row_start + max_rows_per_sheet if row_start + max_rows_per_sheet < total_rows else total_rows
    print "doing sheet:", sheetnum, "rows: %d to %d" % (start, end)
    write_sheet(sortedinput[start:end], sheetnum, header)
    sheetnum += 1

wb.save(outputfilename)
tsvfile.close()        
