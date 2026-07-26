#!/usr/bin/env bash
#
# Golden / regression tests for the PLAAC command-line application.
#
# PLAAC's scoring is deterministic: for a fixed input FASTA, alpha, and core
# length its output is a pure function of the sequences and the built-in
# frequency tables (no randomness, no timestamps, no host-specific paths). That
# makes byte-for-byte golden comparison a valid guard against unintended changes
# to the ~5k-line plaac.java.
#
# Two tests are run:
#   1. golden     -- full output matches a committed expected file exactly.
#   2. recovery   -- the four classic yeast prions each get a called PrLD
#                    (PRDlen > 0). This is a semantic check that survives a
#                    deliberate, reviewed regeneration of the golden file.
#
# Usage:
#   cli/tests/run_tests.sh            # builds plaac.jar if missing, then tests
#   cli/tests/run_tests.sh --no-build # require an existing target/plaac.jar
#
# Exits non-zero if any test fails.

set -u

# --- locate directories (this script lives in cli/tests) --------------------
TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLI_DIR="$(cd "$TESTS_DIR/.." && pwd)"
JAR="$CLI_DIR/target/plaac.jar"
EXAMPLE="$CLI_DIR/example/four_classic_prions.fasta"
EXPECTED="$TESTS_DIR/expected/four_classic_prions.a1.c60.txt"

# Pin the locale so number formatting (decimal separator) is reproducible on
# CI hosts configured with a comma-decimal locale.
export LC_ALL=C
JAVA_LOCALE_FLAGS=(-Duser.language=en -Duser.country=US)

BUILD=1
[[ "${1:-}" == "--no-build" ]] && BUILD=0

fail=0
pass_msg() { printf '  \033[32mPASS\033[0m  %s\n' "$1"; }
fail_msg() { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; fail=1; }

# --- ensure the jar exists --------------------------------------------------
if [[ ! -f "$JAR" ]]; then
    if [[ "$BUILD" -eq 1 ]]; then
        echo "plaac.jar not found; building..."
        ( cd "$CLI_DIR" && ./build_plaac.sh ) || { echo "build failed"; exit 2; }
    else
        echo "plaac.jar not found at $JAR (run without --no-build to build it)"; exit 2
    fi
fi

command -v java >/dev/null 2>&1 || { echo "java not found on PATH"; exit 2; }

echo "Running PLAAC golden tests..."
# Written to a stable (gitignored) path rather than a temp file so CI can upload
# it as an artifact -- useful for diagnosing a cross-platform golden mismatch.
ACTUAL="$TESTS_DIR/actual_output.txt"

java "${JAVA_LOCALE_FLAGS[@]}" -jar "$JAR" -i "$EXAMPLE" -a 1 -c 60 -d > "$ACTUAL" 2>/dev/null

# --- test 1: byte-for-byte golden -------------------------------------------
if diff -u "$EXPECTED" "$ACTUAL" > /tmp/plaac_golden.diff 2>&1; then
    pass_msg "golden output matches expected (four_classic_prions, -a 1 -c 60)"
else
    fail_msg "golden output differs from expected:"
    sed 's/^/        /' /tmp/plaac_golden.diff | head -40
    echo "        (if this change is intentional, regenerate with tests/regenerate_golden.sh)"
fi

# --- test 2: known-positive recovery (PRDlen > 0 for all four prions) -------
# Find the PRDlen column index from the header row so the check does not depend
# on a hard-coded column position.
recovery_report="$(
    awk -F'\t' '
        /^#/ { next }
        $1 == "SEQid" {
            for (i = 1; i <= NF; i++) if ($i == "PRDlen") col = i
            next
        }
        NF > 1 {
            n++
            if (col == 0)            { print "NO_PRDLEN_COLUMN"; bad=1; next }
            if ($col + 0 <= 0)       { print "NOCALL " $1 " PRDlen=" $col; bad=1 }
            else                     { print "OK " $1 " PRDlen=" $col }
        }
        END { if (n == 0) print "NO_DATA_ROWS"; else print "COUNT " n }
    ' "$ACTUAL"
)"

if echo "$recovery_report" | grep -qE 'NO_PRDLEN_COLUMN|NO_DATA_ROWS|NOCALL'; then
    fail_msg "PrLD recovery check failed:"
    echo "$recovery_report" | sed 's/^/        /'
elif [[ "$(echo "$recovery_report" | grep -c '^OK ')" -eq 4 ]]; then
    pass_msg "all four classic prions call a PrLD (PRDlen > 0)"
    echo "$recovery_report" | grep '^OK ' | sed 's/^/        /'
else
    fail_msg "expected 4 PrLD-positive proteins; got:"
    echo "$recovery_report" | sed 's/^/        /'
fi

echo
if [[ "$fail" -eq 0 ]]; then
    echo "All tests passed."
else
    echo "Some tests FAILED."
fi
exit "$fail"
