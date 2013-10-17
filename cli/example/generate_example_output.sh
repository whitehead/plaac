#!/bin/bash
SCRIPT_DIR=`dirname $0`

#
# Build spewey.jar
#
$SCRIPT_DIR/../build.sh

#
# Run spewey.jar
#

# params
fasta=$SCRIPT_DIR/TAIR10_pep_20101214.fasta
core_length=60
alpha=1
output_file=output.tsv

java -jar $SCRIPT_DIR/../target/spewey.jar $fasta $core_length $alpha > $output_file
