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

- A Java runtime on `PATH` (same as the jar itself).
- A built `plaac.jar`. Build it once from the repo root:

      cd cli && ./build_plaac.sh      # produces cli/target/plaac.jar

  The wrapper finds the jar automatically (see *Locating the jar* below).

Install
-------

    pip install ./cli/python          # or: pip install "./cli/python[pandas]"

Or, without installing, add the directory to `PYTHONPATH` / `sys.path` and
`import plaac`.

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
3. If neither is given, it auto-discovers `cli/target/plaac.jar` (built) and
   `web/bin/plaac.jar` (shipped), relative to the module.

Reproducibility
---------------

The wrapper pins the JVM locale (`-Duser.language=en -Duser.country=US`) so
number formatting does not depend on the host locale — the same reason the test
suite in `cli/tests/` does.

Testing
-------

    python3 -m unittest discover -s cli/python -v

Pure-parsing tests run anywhere; the end-to-end tests self-skip if `plaac.jar`
or Java is unavailable.
