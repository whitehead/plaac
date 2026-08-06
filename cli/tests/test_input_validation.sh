#!/usr/bin/env bash
#
# Tests for command-line argument validation and graceful failure.
#
# These lock in the behavior added alongside them:
#   - a non-numeric value for a numeric option reports a usage error on stderr
#     and exits non-zero, instead of throwing a raw Java stack trace;
#   - an unknown option is reported on stderr, not stdout (so it never pollutes
#     the TSV output);
#   - a value-taking flag with no value (e.g. a trailing "-p") is warned about
#     rather than silently ignored, but the run still completes.
#
# Usage:
#   cli/tests/test_input_validation.sh            # builds plaac.jar if missing
#   cli/tests/test_input_validation.sh --no-build # require an existing jar
#
# Exits non-zero if any test fails.

set -u

TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLI_DIR="$(cd "$TESTS_DIR/.." && pwd)"
JAR="$CLI_DIR/target/plaac.jar"
EXAMPLE="$CLI_DIR/example/four_classic_prions.fasta"

export LC_ALL=C

BUILD=1
[[ "${1:-}" == "--no-build" ]] && BUILD=0

fail=0
pass_msg() { printf '  \033[32mPASS\033[0m  %s\n' "$1"; }
fail_msg() { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; fail=1; }

if [[ ! -f "$JAR" ]]; then
    if [[ "$BUILD" -eq 1 ]]; then
        echo "plaac.jar not found; building..."
        ( cd "$CLI_DIR" && ./build_plaac.sh ) || { echo "build failed"; exit 2; }
    else
        echo "plaac.jar not found at $JAR (run without --no-build to build it)"; exit 2
    fi
fi
command -v java >/dev/null 2>&1 || { echo "java not found on PATH"; exit 2; }

OUT="$(mktemp)"; ERR="$(mktemp)"
trap 'rm -f "$OUT" "$ERR"' EXIT

run() { # run <args...>; sets $rc and fills $OUT/$ERR
    java -jar "$JAR" "$@" >"$OUT" 2>"$ERR"; rc=$?
}

echo "Running PLAAC input-validation tests..."

# 1. sanity: a valid invocation still succeeds
run -i "$EXAMPLE" -a 1 -c 60
if [[ "$rc" -eq 0 && "$(grep -cE '^(Sup35p|Ure2p|Rnq1p|Mot3p)' "$OUT")" -eq 4 ]]; then
    pass_msg "valid invocation exits 0 and prints all four proteins"
else
    fail_msg "valid invocation regressed (exit=$rc, rows=$(grep -cE '^(Sup35p|Ure2p|Rnq1p|Mot3p)' "$OUT"))"
fi

# 2. non-numeric alpha: usage error on stderr, non-zero exit, no stack trace
run -i "$EXAMPLE" -a foo -c 60
if [[ "$rc" -ne 0 ]] && grep -qi "requires a number" "$ERR" && ! grep -qi "Exception in thread" "$ERR"; then
    pass_msg "non-numeric -a exits non-zero with a usage error (no stack trace)"
else
    fail_msg "non-numeric -a not handled cleanly (exit=$rc)"; sed 's/^/        /' "$ERR" | head -5
fi

# 3. non-numeric core length: usage error, non-zero exit
run -i "$EXAMPLE" -c foo
if [[ "$rc" -ne 0 ]] && grep -qi "requires an integer" "$ERR" && ! grep -qi "Exception in thread" "$ERR"; then
    pass_msg "non-numeric -c exits non-zero with a usage error (no stack trace)"
else
    fail_msg "non-numeric -c not handled cleanly (exit=$rc)"; sed 's/^/        /' "$ERR" | head -5
fi

# 4. unknown option: reported on stderr, never on stdout
run -i "$EXAMPLE" -Z -a 1 -c 60
if ! grep -qi "skipping" "$OUT" && grep -qi "skipping unknown option" "$ERR"; then
    pass_msg "unknown option is reported on stderr, not stdout"
else
    fail_msg "unknown option leaked to stdout or was not reported on stderr"
fi

# 5. trailing value-flag with no value: warned, but run still completes
run -i "$EXAMPLE" -a 1 -c 60 -p
if [[ "$rc" -eq 0 ]] && grep -qi "requires a value but was given none" "$ERR" \
   && [[ "$(grep -cE '^(Sup35p|Ure2p|Rnq1p|Mot3p)' "$OUT")" -eq 4 ]]; then
    pass_msg "trailing value-flag warns on stderr but the run completes"
else
    fail_msg "trailing value-flag not handled as expected (exit=$rc)"; sed 's/^/        /' "$ERR" | head -5
fi

echo
if [[ "$fail" -eq 0 ]]; then
    echo "All input-validation tests passed."
else
    echo "Some input-validation tests FAILED."
fi
exit "$fail"
