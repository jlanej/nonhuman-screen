"""Tests for the CLI's failure wiring.

The point of these is narrow but important: a kraken2 failure must reach the
shell as a distinct non-zero exit code, never as exit 0 with a zero-filled
("clean-looking") summary.  kraken2 itself is not needed — ``shutil.which`` is
stubbed so the binary check passes.
"""

from unittest import mock

from nonhuman_screen import cli
from nonhuman_screen.engine import Kraken2Error


def _args(*extra):
    return ["classify", "--bam", "s.bam", "--kraken2-db", "db", *extra]


class TestExitCodes:
    def test_missing_kraken2_binary_exits_2(self):
        with mock.patch.object(cli.shutil, "which", return_value=None):
            assert cli.main(_args()) == 2

    def test_kraken2_failure_exits_3_whole_bam(self, capsys):
        with mock.patch.object(cli.shutil, "which", return_value="/bin/kraken2"), \
             mock.patch.object(
                 cli, "_classify_all", side_effect=Kraken2Error("db not found"),
             ):
            assert cli.main(_args()) == 3
        assert "db not found" in capsys.readouterr().err

    def test_kraken2_failure_exits_3_variant_mode(self, capsys):
        with mock.patch.object(cli.shutil, "which", return_value="/bin/kraken2"), \
             mock.patch.object(
                 cli, "_classify_variants",
                 side_effect=Kraken2Error("exited with code 1"),
             ):
            assert cli.main(_args("--variants", "v.vcf.gz")) == 3
        assert "exited with code 1" in capsys.readouterr().err

    def test_success_exits_0(self):
        with mock.patch.object(cli.shutil, "which", return_value="/bin/kraken2"), \
             mock.patch.object(cli, "_classify_all", return_value=0):
            assert cli.main(_args()) == 0
