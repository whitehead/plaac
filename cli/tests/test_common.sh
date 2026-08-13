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
