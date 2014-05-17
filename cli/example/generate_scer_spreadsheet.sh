#!/bin/bash
SCRIPT_DIR=`dirname $0`

#
# Build plaac.jar
#
$SCRIPT_DIR/../build_plaac.sh

#
# Run plaac.jar
#

# params
fasta=$SCRIPT_DIR/Scer.fasta
core_length=60
alpha=1
output_file=Scer-all-proteins-$(date +%Y-%m-%d).out
output_xls=Scer-all-proteins-$(date +%Y-%m-%d).xls

# run plaac to get hits
cmd="java -jar $SCRIPT_DIR/../target/plaac.jar -i $fasta -c $core_length -a $alpha > $output_file"
echo $cmd
eval $cmd

# make fancy spreadsheet
cmd="$SCRIPT_DIR/../../web/bin/fancy-yeast-spreadsheet.py  $output_file - $output_xls ${fasta}.list saccharomyces_cerevisiae.gff"
echo $cmd
eval $cmd
