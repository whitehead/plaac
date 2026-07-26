#!/usr/bin/env bash
#
# Regenerate the committed golden output. Run this ONLY when a change to
# plaac.java is expected to change the numeric output, and after you have
# reviewed the diff and confirmed the new numbers are correct. The golden file
# is a guard, not ground truth -- regenerating it silently defeats the test.
#
# Usage: cli/tests/regenerate_golden.sh

set -euo pipefail
TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLI_DIR="$(cd "$TESTS_DIR/.." && pwd)"
JAR="$CLI_DIR/target/plaac.jar"
EXAMPLE="$CLI_DIR/example/four_classic_prions.fasta"
EXPECTED="$TESTS_DIR/expected/four_classic_prions.a1.c60.txt"

[[ -f "$JAR" ]] || ( cd "$CLI_DIR" && ./build_plaac.sh )

export LC_ALL=C
java -Duser.language=en -Duser.country=US -jar "$JAR" \
    -i "$EXAMPLE" -a 1 -c 60 -d > "$EXPECTED"

echo "Wrote $EXPECTED"
echo "Review the diff carefully before committing:"
echo "  git diff -- $EXPECTED"
