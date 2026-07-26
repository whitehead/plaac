PLAAC command-line tests
========================

Regression tests for the PLAAC CLI. There is no external test framework
dependency — the harness is plain `bash` + `java`, the same tools already
needed to build and run PLAAC.

Running
-------

From the repository root (or anywhere):

    cli/tests/run_tests.sh

This builds `cli/target/plaac.jar` if it is missing, then runs the tests. Pass
`--no-build` to require a pre-built jar instead.

What is tested
--------------

1. **Golden output** — PLAAC is run on the bundled
   `cli/example/four_classic_prions.fasta` with `-a 1 -c 60 -d`, and its full
   output is compared byte-for-byte against
   [`expected/four_classic_prions.a1.c60.txt`](expected/four_classic_prions.a1.c60.txt).
   PLAAC's scoring is deterministic (Viterbi parse over fixed frequency tables,
   no randomness, no timestamps or host paths in the output), so any difference
   is a real change in behavior. This locks the scientific output — including
   the built-in frequency tables and the column set the web UI's HAML is
   generated from — against unintended changes to `plaac.java`.

2. **Known-positive recovery** — Sup35p, Ure2p, Rnq1p, and Mot3p are bona fide
   yeast prions; the harness asserts each is assigned a called PrLD
   (`PRDlen > 0`). This is a semantic check that still holds if the golden file
   is deliberately regenerated after a reviewed numeric change.

Reproducibility note: the harness pins the JVM locale
(`-Duser.language=en -Duser.country=US`, `LC_ALL=C`) so number formatting is
identical regardless of the host's default locale.

Regenerating the golden file
----------------------------

The golden file is a guard, not ground truth. Regenerate it **only** when a
change to `plaac.java` is *expected* to change the output, and only after
reviewing the diff:

    cli/tests/regenerate_golden.sh
    git diff -- cli/tests/expected/    # review before committing
