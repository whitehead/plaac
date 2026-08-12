PLAAC Python API
================

A small Python wrapper around the PLAAC command-line jar, for driving PLAAC from
scripts and pipelines instead of shelling out by hand. It locates the jar, runs
it with reproducible settings, checks for failure, and parses the per-protein
output into structured rows.

The core module is **standard-library only**. pandas is optional and used only by
`PlaacRun.to_dataframe()`.

Prerequisites
-------------

- A Java runtime on `PATH` (to run the jar).
- A `plaac.jar`. The published wheel **bundles a precompiled jar**, so an
  installed package is self-contained. When working from a source checkout,
  build it once:

      cd cli && ./build_plaac.sh      # produces cli/target/plaac.jar

  The wrapper finds the jar automatically (see *Locating the jar* below).

Install
-------

Install the wheel from a GitHub release (it bundles `plaac.jar`), or build a
self-contained wheel from a source checkout after building the jar:

    cd cli && ./build_plaac.sh
    cp target/plaac.jar python/plaac/plaac.jar     # bundle the built jar
    pip install ./python                            # or: pip install "./python[pandas]"

The package version is derived from the Git tag via `setuptools-scm` — the same
single source of truth the CLI build uses — so the CLI and the Python package in
a release report the same version.

For development, add the directory to `PYTHONPATH` / `sys.path` and `import
plaac`; the wrapper will use a jar built in `cli/target/`.

Usage
-----

```python
import plaac

# Per-protein (summary) mode: run.rows is a list of dicts, one per protein.
run = plaac.run("cli/example/four_classic_prions.fasta", alpha=1, core_length=60)

for row in plaac.called_prlds(run.rows):        # proteins with PRDlen > 0
    print(row["SEQid"], row["PRDscore"], row["PRDlen"])

df = run.to_dataframe()                          # optional; requires pandas

# Per-residue mode (for plotting): run.rows is None, data is in run.stdout.
per_res = plaac.run("input.fasta", print_list="all")
```

Key parameters of `plaac.run()` (they mirror the CLI):

| Parameter | CLI flag | Notes |
|---|---|---|
| `alpha` | `-a` | background-frequency blend (default 1.0 = S. cerevisiae) |
| `core_length` | `-c` | minimum called-domain length (default 60) |
| `background_fasta` | `-b` | compute background from another FASTA |
| `background_freqs` | `-B` | background from a frequency table |
| `print_list` | `-p` | `"all"` or a path → per-residue output |
| `heap` | `-Xmx` | e.g. `"4g"`; reproducible memory instead of a RAM fraction |
| `timeout` | — | seconds; PLAAC prints no progress, so this bounds hung runs |
| `check` | — | raise `PlaacError` on non-zero exit (default `True`) |
| `extra_args` | — | appended verbatim for options not modeled here |

Locating the jar
----------------

`plaac.run()` (via `plaac.find_jar()`) resolves the jar as follows:

1. an explicit `jar=...` argument, or
2. the `PLAAC_JAR` environment variable.

   Either is treated as an instruction: if given but no file exists there,
   resolution **fails** rather than silently using a different jar.
3. If neither is given, it looks for a jar **bundled inside the installed
   package** (`plaac/plaac.jar`, present in the released wheel), then falls back
   to `cli/target/plaac.jar` (built) in a source checkout.

(The web app's `web/bin/plaac.jar` is deliberately *not* searched — it is
managed separately and would be the most likely to differ in version.)

Version matching
----------------

The Java implementation and the Python package must be the same version. By
default `plaac.run()` checks the jar's version (via `plaac.jar -V`) against
`plaac.__version__` and raises on a mismatch — so a stale or mismatched jar is
surfaced, not silently used. A **release wheel bundles a matching jar** by
construction, so this never fires there. In a source checkout the installed
package version may be `"unknown"`, in which case the check is a no-op; if you
install the package *and* a jar built from a different checkout is found, the
check tells you to rebuild. Pass `check_jar_version=False` to bypass it.

Reproducibility
---------------

The wrapper pins the JVM locale (`-Duser.language=en -Duser.country=US`) so
number formatting does not depend on the host locale — the same reason the test
suite in `cli/tests/` does.

Testing
-------

    pip install pytest      # or: pip install "./cli/python[test]"
    pytest cli/python -v

Pure-parsing tests run anywhere; the end-to-end tests self-skip if `plaac.jar`
or Java is unavailable.
