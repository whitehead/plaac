#!/usr/bin/env bash
#
# Tests that SEQid always occupies exactly one field of the tab-delimited output.
#
# Sequence names are taken verbatim from the FASTA header, so before the fix that
# accompanies these tests a tab inside a header split SEQid across several fields and
# shifted every column after it. Nothing errored: the header row kept its 41 columns
# while data rows gained one per tab, so anything reading the table by column index
# silently read the wrong values. MitoCarta and several other curated FASTA
# distributions ship tab-separated headers, which is how this is met in practice.
#
# A second, related inconsistency: the name captured in fastareader.nextfasta()
# discarded the result of its trim(), so the first record (whose name comes from
# hasmorefastas(), which does trim) lost trailing whitespace while every later record
# kept it. The same file could therefore produce a different column count for its
# first record than for the rest.
#
# Usage:
#   cli/tests/test_seqid_field.sh            # builds plaac.jar if missing
#   cli/tests/test_seqid_field.sh --no-build # require an existing jar
#
# Exits non-zero if any test fails.

source "$(dirname "$0")/test_common.sh"

WORK="$(mktemp -d)"
trap 'rm -f "$OUT" "$ERR"; rm -rf "$WORK"' EXIT

SEQ=MQNSNQSQNQGQFQQNNMQQQQQQQQQQNQFQQNMPMHQFNMQNQGQFQQNGMQPQFHQQ

# field counts of the header row and of every data row, ignoring "#" comment lines
header_cols() { grep -v '^#' "$1" | awk -F'\t' 'NR==1 {print NF; exit}'; }
data_cols()   { grep -v '^#' "$1" | awk -F'\t' 'NR>1 {print NF}' | sort -u | tr '\n' ' '; }

echo "Running SEQid field tests..."

# 1. a tab inside the header must not add columns
printf '>sp|P00001|TEST_HUMAN\tGeneID:1\tMYGENE\n%s\n' "$SEQ" > "$WORK/tab_header.fasta"
run -i "$WORK/tab_header.fasta"
h="$(header_cols "$OUT")"; d="$(data_cols "$OUT")"
if [[ "$rc" -eq 0 && "$d" == "$h " ]]; then
    pass_msg "tab in header keeps the data row at $h columns"
else
    fail_msg "tab in header changed the column count (header=$h data=$d exit=$rc)"
fi

# 2. the tab-separated and space-separated forms of one header agree in every field
printf '>sp|P00001|TEST_HUMAN GeneID:1 MYGENE\n%s\n' "$SEQ" > "$WORK/space_header.fasta"
run -i "$WORK/space_header.fasta"; grep -v '^#' "$OUT" > "$WORK/space.tsv"
run -i "$WORK/tab_header.fasta";   grep -v '^#' "$OUT" > "$WORK/tab.tsv"
if diff -q "$WORK/space.tsv" "$WORK/tab.tsv" >/dev/null; then
    pass_msg "tab-separated and space-separated headers give identical tables"
else
    fail_msg "tab-separated header does not match the space-separated form"
    diff "$WORK/space.tsv" "$WORK/tab.tsv" | sed 's/^/        /' | head -6
fi

# 3. every record of a multi-record file gets the same column count, whatever
#    trailing whitespace its header carries
printf '>REC_ONE\t\n%s\n>REC_TWO\t\n%s\n' "$SEQ" "$SEQ" > "$WORK/trailing_tab.fasta"
run -i "$WORK/trailing_tab.fasta"
h="$(header_cols "$OUT")"; d="$(data_cols "$OUT")"
if [[ "$rc" -eq 0 && "$d" == "$h " ]]; then
    pass_msg "trailing tabs give every record $h columns"
else
    fail_msg "trailing tab produced uneven rows (header=$h data=$d exit=$rc)"
fi

printf '>REC_ONE   \n%s\n>REC_TWO   \n%s\n' "$SEQ" "$SEQ" > "$WORK/trailing_space.fasta"
run -i "$WORK/trailing_space.fasta"
names="$(grep -v '^#' "$OUT" | awk -F'\t' 'NR>1 {print $1}' | tr '\n' ',')"
if [[ "$names" == "REC_ONE,REC_TWO," ]]; then
    pass_msg "trailing spaces are trimmed consistently across records"
else
    fail_msg "trailing spaces handled inconsistently (SEQids: $names)"
fi

# 4. names without control characters must be untouched, embedded spaces and all
run -i "$CLI_DIR/example/four_classic_prions.fasta"
if [[ "$rc" -eq 0 && "$(grep -cE '^(Sup35p|Ure2p|Rnq1p|Mot3p)' "$OUT")" -eq 4 ]]; then
    pass_msg "ordinary names are unchanged"
else
    fail_msg "ordinary names regressed (exit=$rc)"
fi

# 5. the per-residue table (-p) puts SEQid in one field as well
run -i "$WORK/space_header.fasta" -p all; grep -v '^#' "$OUT" | grep -v '^###' > "$WORK/space_res.tsv"; rc_s=$rc
run -i "$WORK/tab_header.fasta"   -p all; grep -v '^#' "$OUT" | grep -v '^###' > "$WORK/tab_res.tsv"
if [[ "$rc" -eq 0 && "$rc_s" -eq 0 ]] && diff -q "$WORK/space_res.tsv" "$WORK/tab_res.tsv" >/dev/null; then
    pass_msg "per-residue table matches between tab- and space-separated headers"
else
    fail_msg "per-residue table differs between tab- and space-separated headers (exit=$rc)"
    diff "$WORK/space_res.tsv" "$WORK/tab_res.tsv" | sed 's/^/        /' | head -4
fi

echo
if [[ "$fail" -eq 0 ]]; then
    echo "All SEQid field tests passed."
else
    echo "Some SEQid field tests FAILED."
fi
exit "$fail"
