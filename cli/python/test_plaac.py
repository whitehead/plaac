"""
Tests for the plaac Python wrapper. Standard-library unittest only (no pytest).

Requires a built cli/target/plaac.jar (CI builds it before running these) and a
Java runtime on PATH; the end-to-end tests self-skip if either is missing so the
pure-parsing tests still run anywhere.

Run:  python3 -m unittest discover -s cli/python -v
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plaac  # noqa: E402

_CLI_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_EXAMPLE = os.path.join(_CLI_DIR, "example", "four_classic_prions.fasta")
_CLASSIC = ["Sup35p", "Ure2p", "Rnq1p", "Mot3p"]


def _jar_available():
    try:
        plaac.find_jar()
    except plaac.PlaacError:
        return False
    return shutil.which("java") is not None


# A small, realistic snippet of PLAAC summary output for pure-parsing tests.
_SAMPLE = (
    "## alpha=1.0; corelength=60;\n"
    "SEQid\tNLLR\tPRDstart\tPRDlen\tPROTlen\tCOREaa\n"
    "Sup35p\t0.854\t1\t133\t685\tNQGNNQQNY\n"
    "ShortSeq\tNaN\t-1\t0\t20\t-\n"
)


class TestParsing(unittest.TestCase):
    def test_parse_types_and_na(self):
        rows = plaac.parse_summary(_SAMPLE)
        self.assertEqual(len(rows), 2)
        r0 = rows[0]
        self.assertEqual(r0["SEQid"], "Sup35p")          # name kept as str
        self.assertIsInstance(r0["NLLR"], float)         # float coercion
        self.assertIsInstance(r0["PRDlen"], int)         # int coercion
        self.assertEqual(r0["PRDlen"], 133)
        self.assertEqual(r0["COREaa"], "NQGNNQQNY")      # AA column kept as str
        self.assertIsNone(rows[1]["NLLR"])               # NaN -> None
        self.assertEqual(rows[1]["COREaa"], "-")         # '-' left intact in AA col

    def test_comment_lines_skipped(self):
        rows = plaac.parse_summary(_SAMPLE)
        self.assertTrue(all("SEQid" in r for r in rows))

    def test_called_prlds(self):
        rows = plaac.parse_summary(_SAMPLE)
        called = plaac.called_prlds(rows)
        self.assertEqual([r["SEQid"] for r in called], ["Sup35p"])

    def test_empty_or_comments_only_returns_empty(self):
        self.assertEqual(plaac.parse_summary(""), [])
        self.assertEqual(plaac.parse_summary("## params\n# note\n"), [])

    def test_header_only_returns_empty(self):
        self.assertEqual(plaac.parse_summary("SEQid\tPRDlen\n"), [])

    def test_missing_header_raises(self):
        # Data present but no 'SEQid' header -> malformed, must not silently []
        with self.assertRaises(plaac.PlaacError):
            plaac.parse_summary("Sup35p\t0.5\t100\n")

    def test_overlong_row_raises(self):
        # A data row with more fields than the header signals format drift
        with self.assertRaises(plaac.PlaacError):
            plaac.parse_summary("SEQid\tPRDlen\nfoo\t10\textra\n")

    def test_prd_structured_columns_kept_as_strings(self):
        txt = ("SEQid\tPRDcount\tPRDall\tPRDothers\n"
               "S\t1\t{[1-133 (51.2 @ 5)];}\t{}\n")
        row = plaac.parse_summary(txt)[0]
        self.assertEqual(row["PRDcount"], 1)                 # count stays int
        self.assertEqual(row["PRDall"], "{[1-133 (51.2 @ 5)];}")  # structured -> str
        self.assertEqual(row["PRDothers"], "{}")


class TestErrors(unittest.TestCase):
    def test_missing_jar_raises(self):
        with self.assertRaises(plaac.PlaacError):
            plaac.find_jar("/nonexistent/definitely/not/plaac.jar")

    def test_extra_args_str_raises_typeerror(self):
        # A bare string would explode into single characters via list(); reject it
        # (raised before any jar/env work, so no jar needed).
        with self.assertRaises(TypeError):
            plaac.run(_EXAMPLE, extra_args="-s")

    @unittest.skipUnless(_jar_available(), "plaac.jar / java not available")
    def test_missing_input_raises(self):
        with self.assertRaises(plaac.PlaacError):
            plaac.run("/nonexistent/input.fasta")


@unittest.skipUnless(_jar_available(), "plaac.jar / java not available")
class TestEndToEnd(unittest.TestCase):
    def test_four_classic_prions(self):
        run = plaac.run(_EXAMPLE, alpha=1, core_length=60)
        self.assertEqual(run.returncode, 0)
        self.assertIsNotNone(run.rows)
        names = [r["SEQid"] for r in run.rows]
        self.assertEqual(sorted(names), sorted(_CLASSIC))

    def test_all_four_are_prd_positive(self):
        run = plaac.run(_EXAMPLE, alpha=1, core_length=60)
        called = plaac.called_prlds(run.rows)
        self.assertEqual(sorted(r["SEQid"] for r in called), sorted(_CLASSIC))

    def test_per_residue_mode_returns_text_not_rows(self):
        run = plaac.run(_EXAMPLE, alpha=1, core_length=60, print_list="all")
        self.assertIsNone(run.rows)
        self.assertTrue(len(run.stdout) > 0)

    def test_failed_run_yields_no_rows(self):
        # A bad numeric arg makes the jar exit non-zero; with check=False the
        # wrapper must not fabricate rows from failed output.
        run = plaac.run(_EXAMPLE, extra_args=["-c", "not_a_number"], check=False)
        self.assertNotEqual(run.returncode, 0)
        self.assertIsNone(run.rows)

    def test_to_dataframe_raises_in_per_residue_mode(self):
        run = plaac.run(_EXAMPLE, print_list="all")
        with self.assertRaises(plaac.PlaacError):
            run.to_dataframe()


@unittest.skipUnless(_jar_available(), "plaac.jar / java not available")
class TestPassThroughOptions(unittest.TestCase):
    """The thin option pass-throughs: verify each flag reaches the jar's command
    line and the run still succeeds."""

    def test_heap_option_sets_xmx(self):
        run = plaac.run(_EXAMPLE, alpha=1, core_length=60, heap="256m")
        self.assertIn("-Xmx256m", run.command)
        self.assertEqual(run.returncode, 0)
        self.assertEqual(len(run.rows), 4)

    def test_background_fasta_option(self):
        # alpha < 1 so the -b background actually influences scoring.
        run = plaac.run(_EXAMPLE, alpha=0.5, core_length=60, background_fasta=_EXAMPLE)
        self.assertIn("-b", run.command)
        self.assertIn(_EXAMPLE, run.command)
        self.assertEqual(run.returncode, 0)
        self.assertEqual(len(run.rows), 4)

    def test_background_freqs_option(self):
        # Generate a frequency table the documented way ('-b FASTA' with no -i
        # prints one), then feed it back via -B.
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
            self.assertIn("-B", run.command)
            self.assertIn(freqs_path, run.command)
            self.assertEqual(run.returncode, 0)
            self.assertEqual(len(run.rows), 4)
        finally:
            os.unlink(freqs_path)


if __name__ == "__main__":
    unittest.main()
