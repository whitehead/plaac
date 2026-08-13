"""
plaac.py -- a small Python API for the PLAAC command-line tool.

PLAAC (https://github.com/whitehead/plaac) is distributed as a Java jar. This
module wraps it so PLAAC can be driven from Python without hand-rolling a
subprocess call every time: it locates the jar, runs it with reproducible
settings, checks for failure, and parses the tab-separated per-protein output
into structured rows.

Design notes:
  * Standard library only. ``PlaacRun.to_dataframe()`` uses pandas if it is
    installed, but pandas is not required for anything else.
  * The JVM locale is pinned so number formatting is identical regardless of the
    host locale (a comma-decimal locale would otherwise change the output text).
  * ``timeout`` is exposed because PLAAC prints no progress: on a large proteome
    a slow run is indistinguishable from a hung one except by this bound.
  * ``heap`` sets ``-Xmx`` explicitly so memory use is reproducible across
    machines rather than defaulting to a fraction of physical RAM.

Example:
    >>> import plaac
    >>> run = plaac.run("cli/example/four_classic_prions.fasta", alpha=1, core_length=60)
    >>> [r["SEQid"] for r in plaac.called_prlds(run.rows)]
    ['Sup35p', 'Ure2p', 'Rnq1p', 'Mot3p']
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Optional, Sequence, Union

# Package version, derived from the Git tag via setuptools-scm at build time and
# read here from the installed distribution metadata (the same single source of
# truth the Java CLI uses). "unknown" when running from an unbuilt source tree.
try:
    from importlib.metadata import version as _pkg_version, PackageNotFoundError
    try:
        __version__ = _pkg_version("plaac")
    except PackageNotFoundError:
        __version__ = "unknown"
except ImportError:  # pragma: no cover -- importlib.metadata is in 3.8+
    __version__ = "unknown"

DEFAULT_ALPHA = 1.0
DEFAULT_CORE_LENGTH = 60

# Pin number formatting (decimal separator) so output is host-locale independent.
_LOCALE_FLAGS = ["-Duser.language=en", "-Duser.country=US"]

# Placeholders PLAAC writes for undefined values.
_NA_TOKENS = {"NA", "NaN", "nan", ""}

# Columns that hold names or amino-acid sequences; kept as strings (never coerced
# to numbers, and their "-" placeholder is left intact rather than made None).
_STRING_COLUMNS = {"SEQid", "COREaa", "STARTaa", "ENDaa", "PRDaa", "PAPAaa",
                   "PRDall", "PRDothers"}

PathLike = Union[str, "os.PathLike[str]"]


class PlaacError(RuntimeError):
    """Raised when the PLAAC jar cannot be found or the run fails."""


def find_jar(explicit: Optional[PathLike] = None) -> Path:
    """Locate ``plaac.jar``.

    An explicit ``jar=`` path or a ``$PLAAC_JAR`` setting is treated as an
    instruction, not a hint: if one is given but no file exists there, this
    raises rather than silently falling back to a different jar. Otherwise it
    looks for a jar bundled inside the installed package (so a pip-installed
    wheel is self-contained), then falls back to the jar's usual built and
    shipped locations in a source checkout.
    """
    if explicit is not None:
        p = Path(explicit)
        if p.is_file():
            return p.resolve()
        raise PlaacError(f"plaac.jar not found at the path given (jar=): {explicit}")

    env = os.environ.get("PLAAC_JAR")
    if env:
        p = Path(env)
        if p.is_file():
            return p.resolve()
        raise PlaacError(f"$PLAAC_JAR is set but no file exists there: {env}")

    here = Path(__file__).resolve().parent          # .../plaac (src/plaac in a checkout)
    cli_dir = here.parents[2]                         # .../cli in a source checkout
    # web/bin/plaac.jar is deliberately NOT searched: it is a separately-managed
    # artifact for the web app and is the most likely to differ in version from
    # the Python package, which would trip the version-skew check (check_version).
    discovered = [
        here / "plaac.jar",                          # bundled in the wheel (self-contained)
        cli_dir / "target" / "plaac.jar",            # built in a source checkout
    ]
    for c in discovered:
        if c.is_file():
            return c.resolve()
    searched = "\n  ".join(str(c) for c in discovered)
    raise PlaacError(
        "Could not find plaac.jar. Install the wheel (which bundles it), build "
        "it with cli/build_plaac.sh, or pass jar=... / set $PLAAC_JAR. "
        "Searched:\n  " + searched
    )


@lru_cache(maxsize=None)
def get_jar_version(jar_path: PathLike, java: str = "java") -> str:
    """Return the version a PLAAC jar reports via its ``-V`` flag.

    Cached per (jar, java) so it costs at most one JVM start per jar in a
    process. Raises :class:`PlaacError` if the version cannot be parsed.
    """
    result = subprocess.run(
        [java, "-jar", str(jar_path), "-V"],
        check=True, capture_output=True, text=True,
    )
    output = result.stdout.strip()
    prefix = "PLAAC Version: "                        # human-readable prefix to strip
    if not output.startswith(prefix):
        raise PlaacError(f"Could not determine PLAAC jar version from: {output!r}")
    return output[len(prefix):].strip()


def check_version(jar_path: PathLike, java: str = "java") -> str:
    """Raise if the jar's version differs from the installed package version.

    The jar and the Python bindings can drift out of sync -- especially if a
    prebuilt jar is picked up -- so we require them to match exactly. This is a
    no-op when the package version is ``"unknown"`` (e.g. running from a source
    tree that was never installed), so it only enforces a match for installed
    packages; a release wheel's bundled jar always matches by construction.
    Returns the jar version.
    """
    jar_version = get_jar_version(jar_path, java)
    if __version__ != "unknown" and jar_version != __version__:
        raise PlaacError(
            f"PLAAC version mismatch: the Python package is {__version__} but the "
            f"jar at {jar_path} reports {jar_version}. Rebuild the jar from the "
            f"same checkout (cli/build_plaac.sh), or pass jar= / set $PLAAC_JAR."
        )
    return jar_version


def _coerce(value: str) -> Any:
    """Convert an output cell to int/float/None; leave anything else as str."""
    if value in _NA_TOKENS:
        return None
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def parse_summary(text: str) -> list[dict[str, Any]]:
    """Parse PLAAC per-protein (summary) output into a list of dicts.

    Comment lines (starting with ``#``) are skipped; the column header is the
    line starting with ``SEQid``. Numeric cells are coerced to int/float and
    PLAAC's undefined placeholders (``NA``/``NaN``) become ``None``. Name and
    amino-acid-sequence columns are kept as strings.

    Raises :class:`PlaacError` if the output contains data but no column header
    (i.e. it is malformed or from an unexpected mode), or if a data row has more
    fields than the header (a sign the output drifted) -- both are surfaced
    rather than silently returning empty or dropping columns. Output that is
    empty or all comments returns ``[]``.
    """
    header: Optional[list[str]] = None
    rows: list[dict[str, Any]] = []
    saw_content = False
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        saw_content = True
        fields = line.split("\t")
        if header is None:
            if fields[0] == "SEQid":
                header = fields
            continue
        if len(fields) > len(header):
            raise PlaacError(
                "malformed PLAAC output: a data row has more fields "
                f"({len(fields)}) than the header ({len(header)}); "
                f"row starts with {fields[0]!r}"
            )
        row: dict[str, Any] = {}
        for i, name in enumerate(header):
            raw = fields[i] if i < len(fields) else ""
            row[name] = raw if name in _STRING_COLUMNS else _coerce(raw)
        rows.append(row)
    if header is None and saw_content:
        raise PlaacError(
            "could not find the PLAAC column header (a line starting with "
            "'SEQid') in the output; it may be malformed or from an unexpected "
            "mode"
        )
    return rows


def called_prlds(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    """Return the rows with a called prion-like domain (``PRDlen > 0``)."""
    out = []
    for r in rows:
        prdlen = r.get("PRDlen")
        if isinstance(prdlen, (int, float)) and prdlen > 0:
            out.append(r)
    return out


@dataclass
class PlaacRun:
    """The result of a single PLAAC invocation."""

    command: list[str]
    returncode: int
    stdout: str
    stderr: str
    rows: Optional[list[dict[str, Any]]]  # parsed summary rows; None in per-residue mode

    def to_dataframe(self):
        """Return the summary rows as a pandas DataFrame (pandas required)."""
        if self.rows is None:
            raise PlaacError(
                "No summary rows to convert (per-residue mode). Use .stdout instead."
            )
        import pandas as pd  # optional dependency, imported lazily
        return pd.DataFrame(self.rows)


def run(
    input_fasta: PathLike,
    *,
    jar: Optional[PathLike] = None,
    alpha: float = DEFAULT_ALPHA,
    core_length: int = DEFAULT_CORE_LENGTH,
    background_fasta: Optional[PathLike] = None,
    background_freqs: Optional[PathLike] = None,
    print_list: Optional[str] = None,
    java: str = "java",
    heap: Optional[str] = None,
    timeout: Optional[float] = None,
    check: bool = True,
    check_jar_version: bool = True,
    extra_args: Optional[Sequence[str]] = None,
) -> PlaacRun:
    """Run PLAAC on a protein FASTA and return a :class:`PlaacRun`.

    In the default (summary) mode ``run.rows`` holds one parsed dict per protein.
    Pass ``print_list="all"`` or a path to a list of sequence names to get
    per-residue output instead; then ``run.rows`` is ``None`` and the data is in
    ``run.stdout``.

    Parameters mirror the PLAAC CLI: ``alpha`` (``-a``), ``core_length``
    (``-c``), ``background_fasta`` (``-b``), ``background_freqs`` (``-B``),
    ``print_list`` (``-p``). ``heap`` sets ``-Xmx`` (e.g. ``"4g"``); ``timeout``
    is in seconds; ``check=True`` raises on a non-zero exit; ``extra_args`` is
    appended verbatim for options not modeled here. ``check_jar_version=True``
    verifies the jar's version matches this package's (skip with ``False`` for a
    deliberately mismatched jar).
    """
    if isinstance(extra_args, str):
        raise TypeError(
            "extra_args must be a sequence of strings, not a single str "
            f"(got {extra_args!r}); pass e.g. ['-s'], not '-s'"
        )
    jar_path = find_jar(jar)
    if shutil.which(java) is None:
        raise PlaacError(f"'{java}' not found on PATH. Install a Java runtime (JRE).")
    if check_jar_version:
        check_version(jar_path, java)
    input_path = Path(input_fasta)
    if not input_path.is_file():
        raise PlaacError(f"Input FASTA not found: {input_path}")

    cmd = [java, *_LOCALE_FLAGS]
    if heap:
        cmd.append(f"-Xmx{heap}")
    cmd += ["-jar", str(jar_path), "-i", str(input_path),
            "-a", str(alpha), "-c", str(core_length)]
    if background_fasta:
        cmd += ["-b", str(background_fasta)]
    if background_freqs:
        cmd += ["-B", str(background_freqs)]
    if print_list is not None:
        cmd += ["-p", str(print_list)]
    if extra_args is not None:
        cmd += list(extra_args)

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except FileNotFoundError as e:  # java disappeared between the which() check and now
        raise PlaacError(f"Failed to launch java: {e}") from e
    except subprocess.TimeoutExpired as e:
        raise PlaacError(
            f"PLAAC timed out after {timeout}s. PLAAC prints no progress, so a large "
            f"proteome can look hung; raise timeout=... for big inputs."
        ) from e

    if check and proc.returncode != 0:
        raise PlaacError(
            f"PLAAC exited with code {proc.returncode}.\n"
            f"Command: {' '.join(cmd)}\nstderr:\n{proc.stderr}"
        )

    # Only parse in summary mode on a successful run; a failed run (reachable
    # only with check=False) has no meaningful rows.
    rows = None
    if print_list is None and proc.returncode == 0:
        rows = parse_summary(proc.stdout)
    return PlaacRun(
        command=cmd,
        returncode=proc.returncode,
        stdout=proc.stdout,
        stderr=proc.stderr,
        rows=rows,
    )
