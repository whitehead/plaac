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
fasta=$SCRIPT_DIR/TAIR10_pep_20101214.fasta
core_length=60
alpha=1
output_file=output.tsv

java -jar $SCRIPT_DIR/../target/plaac.jar -i $fasta -c $core_length -a $alpha > $output_file
