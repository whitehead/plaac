"""
Tests for the plaac Python wrapper (pytest).

The package uses a src layout, so `import plaac` resolves to the *installed*
package rather than the in-tree source -- these tests therefore exercise the
built/installed wheel. Install the package first (`pip install ./cli/python`, or
`-e` for development). Pure-parsing tests run anywhere; the end-to-end tests are
skipped if a built plaac.jar or a Java runtime is unavailable.

Run:  pip install ./cli/python && pytest cli/python
"""

import os
import subprocess
import tempfile

import pytest

import plaac

# tests/ lives at cli/python/tests/, so the CLI dir is three levels up.
_CLI_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_EXAMPLE = os.path.join(_CLI_DIR, "example", "four_classic_prions.fasta")
_CLASSIC = ["Sup35p", "Ure2p", "Rnq1p", "Mot3p"]


def _jar_available():
    import shutil
    try:
        plaac.find_jar()
    except plaac.PlaacError:
        return False
    return shutil.which("java") is not None


requires_jar = pytest.mark.skipif(
    not _jar_available(), reason="plaac.jar / java not available"
)

# A small, realistic snippet of PLAAC summary output for pure-parsing tests.
_SAMPLE = (
    "## alpha=1.0; corelength=60;\n"
    "SEQid\tNLLR\tPRDstart\tPRDlen\tPROTlen\tCOREaa\n"
    "Sup35p\t0.854\t1\t133\t685\tNQGNNQQNY\n"
    "ShortSeq\tNaN\t-1\t0\t20\t-\n"
)


# --- parsing (no jar needed) -------------------------------------------------

def test_parse_types_and_na():
    rows = plaac.parse_summary(_SAMPLE)
    assert len(rows) == 2
    r0 = rows[0]
    assert r0["SEQid"] == "Sup35p"          # name kept as str
    assert isinstance(r0["NLLR"], float)    # float coercion
    assert isinstance(r0["PRDlen"], int)    # int coercion
    assert r0["PRDlen"] == 133
    assert r0["COREaa"] == "NQGNNQQNY"      # AA column kept as str
    assert rows[1]["NLLR"] is None          # NaN -> None
    assert rows[1]["COREaa"] == "-"         # '-' left intact in AA col


def test_comment_lines_skipped():
    rows = plaac.parse_summary(_SAMPLE)
    assert all("SEQid" in r for r in rows)


def test_called_prlds():
    called = plaac.called_prlds(plaac.parse_summary(_SAMPLE))
    assert [r["SEQid"] for r in called] == ["Sup35p"]


def test_empty_or_comments_only_returns_empty():
    assert plaac.parse_summary("") == []
    assert plaac.parse_summary("## params\n# note\n") == []


def test_header_only_returns_empty():
    assert plaac.parse_summary("SEQid\tPRDlen\n") == []


def test_missing_header_raises():
    # Data present but no 'SEQid' header -> malformed, must not silently return []
    with pytest.raises(plaac.PlaacError):
        plaac.parse_summary("Sup35p\t0.5\t100\n")


def test_overlong_row_raises():
    # A data row with more fields than the header signals format drift
    with pytest.raises(plaac.PlaacError):
        plaac.parse_summary("SEQid\tPRDlen\nfoo\t10\textra\n")


def test_prd_structured_columns_kept_as_strings():
    txt = ("SEQid\tPRDcount\tPRDall\tPRDothers\n"
           "S\t1\t{[1-133 (51.2 @ 5)];}\t{}\n")
    row = plaac.parse_summary(txt)[0]
    assert row["PRDcount"] == 1                          # count stays int
    assert row["PRDall"] == "{[1-133 (51.2 @ 5)];}"      # structured -> str
    assert row["PRDothers"] == "{}"


# --- errors / misuse (no jar needed) ----------------------------------------

def test_missing_jar_raises():
    with pytest.raises(plaac.PlaacError):
        plaac.find_jar("/nonexistent/definitely/not/plaac.jar")


def test_extra_args_str_raises_typeerror():
    # A bare string would explode into single characters via list(); reject it
    # (raised before any jar/env work, so no jar needed).
    with pytest.raises(TypeError):
        plaac.run(_EXAMPLE, extra_args="-s")


@requires_jar
def test_missing_input_raises():
    with pytest.raises(plaac.PlaacError):
        plaac.run("/nonexistent/input.fasta")


# --- end-to-end (require jar + java) ----------------------------------------

@requires_jar
def test_four_classic_prions():
    run = plaac.run(_EXAMPLE, alpha=1, core_length=60)
    assert run.returncode == 0
    assert run.rows is not None
    assert sorted(r["SEQid"] for r in run.rows) == sorted(_CLASSIC)


@requires_jar
def test_all_four_are_prd_positive():
    run = plaac.run(_EXAMPLE, alpha=1, core_length=60)
    called = plaac.called_prlds(run.rows)
    assert sorted(r["SEQid"] for r in called) == sorted(_CLASSIC)


@requires_jar
def test_per_residue_mode_returns_text_not_rows():
    run = plaac.run(_EXAMPLE, alpha=1, core_length=60, print_list="all")
    assert run.rows is None
    assert len(run.stdout) > 0


@requires_jar
def test_failed_run_yields_no_rows():
    # A bad numeric arg makes the jar exit non-zero; with check=False the wrapper
    # must not fabricate rows from failed output.
    run = plaac.run(_EXAMPLE, extra_args=["-c", "not_a_number"], check=False)
    assert run.returncode != 0
    assert run.rows is None


@requires_jar
def test_to_dataframe_raises_in_per_residue_mode():
    run = plaac.run(_EXAMPLE, print_list="all")
    with pytest.raises(plaac.PlaacError):
        run.to_dataframe()


# --- jar/package version-skew check (require jar) ---------------------------

@requires_jar
def test_get_jar_version_returns_string():
    v = plaac.get_jar_version(plaac.find_jar())
    assert isinstance(v, str) and v


@requires_jar
def test_check_version_ok_when_unknown_or_matching():
    # From a source checkout __version__ is "unknown", so check_version is a
    # no-op and returns the jar version without raising.
    v = plaac.check_version(plaac.find_jar())
    assert isinstance(v, str) and v


@requires_jar
def test_check_version_raises_on_mismatch(monkeypatch):
    monkeypatch.setattr(plaac, "__version__", "9.9.9-does-not-match")
    with pytest.raises(plaac.PlaacError):
        plaac.check_version(plaac.find_jar())


@requires_jar
def test_run_skips_version_check_when_disabled(monkeypatch):
    # A deliberately mismatched version must not raise when the check is off.
    monkeypatch.setattr(plaac, "__version__", "9.9.9-does-not-match")
    run = plaac.run(_EXAMPLE, alpha=1, core_length=60, check_jar_version=False)
    assert run.returncode == 0


# --- option pass-throughs (require jar) -------------------------------------

@requires_jar
def test_heap_option_sets_xmx():
    run = plaac.run(_EXAMPLE, alpha=1, core_length=60, heap="256m")
    assert "-Xmx256m" in run.command
    assert run.returncode == 0
    assert len(run.rows) == 4


@requires_jar
def test_background_fasta_option():
    # alpha < 1 so the -b background actually influences scoring.
    run = plaac.run(_EXAMPLE, alpha=0.5, core_length=60, background_fasta=_EXAMPLE)
    assert "-b" in run.command
    assert _EXAMPLE in run.command
    assert run.returncode == 0
    assert len(run.rows) == 4


@requires_jar
def test_background_freqs_option():
    # Generate a frequency table the documented way ('-b FASTA' with no -i prints
    # one), then feed it back via -B.
    jar = str(plaac.find_jar())
    freqs = subprocess.run(
        ["java", "-jar", jar, "-b", _EXAMPLE],
        capture_output=True, text=True, check=True,
    ).stdout
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write(freqs)
        freqs_path = fh.name
    try:
        run = plaac.run(_EXAMPLE, alpha=0.5, core_length=60,
                        background_freqs=freqs_path)
        assert "-B" in run.command
        assert freqs_path in run.command
        assert run.returncode == 0
        assert len(run.rows) == 4
    finally:
        os.unlink(freqs_path)
