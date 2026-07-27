"""Tests for the CLI's failure wiring.

The point of these is narrow but important: a kraken2 failure must reach the
shell as a distinct non-zero exit code, never as exit 0 with a zero-filled
("clean-looking") summary.  kraken2 itself is not needed — ``shutil.which`` is
stubbed so the binary check passes.
"""

from unittest import mock

import pytest

from nonhuman_screen import cli
from nonhuman_screen.engine import Kraken2Error
from nonhuman_screen.result import TaxonomicFractions, VariantNHF


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


def _variant(taxonomy_available):
    return VariantNHF(
        chrom="chr1", pos=99, ref="A", alt="T",
        supporting_read_names=frozenset({"r1", "r2"}),
        fractions=TaxonomicFractions(
            total=2, nonhuman=0.5, bacterial=0.5, archaeal=0.0, fungal=0.0,
            protist=0.0, viral=0.0, univec_core=0.0, human_lineage=0.5,
            unclassified=0.0,
        ),
        taxonomy_available=taxonomy_available,
    )


class TestVariantTsv:
    """The --variants table is the only NHF surface some consumers read."""

    def _tsv(self, tmp_path, taxonomy_available):
        args = mock.Mock(
            variants="v.vcf", bam="s.bam", kraken2_db="db", ref_fasta=None,
            min_baseq=0, confidence=0.0, threads=1, memory_mapping=False,
            out_prefix=str(tmp_path / "out"),
        )
        with mock.patch.object(cli, "_read_vcf_variants", return_value=[]), \
             mock.patch(
                 "nonhuman_screen.bam.classify_variants_alt_reads",
                 return_value=[_variant(taxonomy_available)],
             ):
            assert cli._classify_variants(args) == 0
        return (tmp_path / "out.variant_nhf.tsv").read_text().splitlines()

    @pytest.mark.parametrize("available", [True, False])
    def test_taxonomy_available_column_round_trips(self, tmp_path, available):
        header, row = self._tsv(tmp_path, available)
        cols = dict(zip(header.split("\t"), row.split("\t")))
        assert cols["taxonomy_available"] == str(available)
        # Consumers join on these; adding a column must not disturb them.
        assert cols["variant_key"] == "chr1:99:A:T"
        assert cols["supporting_reads"] == "2"
        assert cols["nonhuman_fraction"] == "0.5"

    def test_header_and_row_widths_match(self, tmp_path):
        """A header/row width mismatch would silently misalign every column."""
        header, row = self._tsv(tmp_path, True)
        assert len(header.split("\t")) == len(cli._TSV_COLUMNS)
        assert len(row.split("\t")) == len(cli._TSV_COLUMNS)

    def test_taxonomy_available_is_last_column(self, tmp_path):
        """Appended last, so a positional reader sees no existing column move."""
        header, _ = self._tsv(tmp_path, True)
        assert header.split("\t")[-1] == "taxonomy_available"
