#!/usr/bin/env bash

# Tests for checking regression changes in output for MOT3 against a gold output
#
# Usage:
#   cli/tests/test_mot3.sh            # builds plaac.jar if missing
#   cli/tests/test_mot3.sh --no-build # require an existing jar


source "$(dirname "$0")/test_common.sh"

echo "Running MOT3 regression CLI test..."

run -i "$CLI_DIR/example/MOT3.fasta"

# normalize by removing version
sed -E 's/^(## plaac_version=).*/\1<version>;/' "$CLI_DIR/example/MOT3-candidates-plaac-1.1.0.tsv" > gold.normalized.tsv
sed -E 's/^(## plaac_version=).*/\1<version>;/' $OUT > output.normalized.tsv

diff -u gold.normalized.tsv output.normalized.tsv || {
    echo "FAIL: MOT3 output differs from expected output"
    exit 1
}

echo "PASS: MOT3 output matches expected output"


