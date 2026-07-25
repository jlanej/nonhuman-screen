# Changelog

All notable changes to `nonhuman-screen` are documented here. The format is
based on [Keep a Changelog](https://keepachangelog.com/), and this project
adheres to [Semantic Versioning](https://semver.org/).

Continuous builds published from `main` carry the version `0.1.<run-number>`
(see the README's *Releasing* section); the entries below are grouped by the
`<major>.<minor>` base they were released under.

## Unreleased

### Changed
- **Breaking:** a kraken2 run that exits non-zero, or whose per-read output does
  not account for every input read, now raises `Kraken2Error` instead of
  returning a zero-filled result. Previously both cases surfaced as
  `nonhuman_fraction = 0.0`, i.e. a failed screen indistinguishable from a clean
  sample. Pass `strict=False` to `Kraken2Runner` / the BAM helpers to restore
  the old behaviour, which now also sets
  `ClassificationResult.classification_failed`.
- `nonhuman-screen classify` exits **3** on classification failure and writes no
  summary file.

### Added
- `Kraken2Error`, exported at the package top level.
- `ClassificationResult.classification_failed`.
- A warning when the database taxonomy contains no human node (taxid 9606), in
  which case the human lineage/clade exclusions are empty and the
  human-homology guard cannot fire.
- Functional failure control: a broken `--kraken2-db` must exit 3 and leave no
  summary behind.

### Fixed
- `tests/test_integration.py` now builds its mini database single-threaded
  (`--threads 1` + `OMP_NUM_THREADS=1`), matching
  `scripts/build_test_fixtures.py`; the previous two-thread build died with
  SIGPIPE on macOS.

### Documentation
- methodology: the per-domain columns are **not** a partition of `nonhuman` —
  `protist` includes `Eukaryota` (2759) and `Opisthokonta` (33154), which are
  human ancestors excluded from `nonhuman`, so `protist` over-reports; and
  whether the domains sum to `nonhuman` is a property of the database (they do
  for PrackenDB).
- methodology/database: the database must contain the human genome; both
  conservatism mechanisms are defined relative to human.
- methodology: documented exactly where `taxonomy_available` is observable — it
  is **not** reachable through `VariantNHF` or `classify --variants`.
- cli: added exit code 3; corrected the whole-BAM scope (one sequence per unique
  read name, not per record) and the per-read semantics of `--min-baseq`.
- testing: documented what the controls do *not* establish — the mini database
  has no human sequence, so the human-homology guard is never exercised
  end-to-end and the per-read false-positive rate is structurally zero.
- Removed stale references to the source repository from shipped docstrings
  (`JellyfishKmerQuery`, "the companion BED file") and from
  `download_kraken2_db.sh`, which now also prints the database directory to
  pass to `--kraken2-db`.

## [0.1.0] — 2026-07-10

Initial release. Extracted from
[`kmer-denovo-filter`](https://github.com/jlanej/kmer_denovo_filter), where the
methodology was originally developed.

### Added
- `Kraken2Runner` engine: kraken2 LCA classification of reads into a non-human
  fraction plus per-domain breakdowns, with a human-homology guard and
  UniVec-Core exclusion (stdlib-only core).
- `ClassificationResult`, `TaxonomicFractions`, `VariantNHF` result types;
  `ClassificationResult.taxonomy_available` flag for fail-closed consumers.
- Allele-based non-human fraction (`[bam]` extra): `classify_variant_alt_reads`,
  `classify_variants_alt_reads` (batched, single kraken2 invocation),
  `classify_reads_from_bam`, `reads_supporting_alt`, `read_supports_alt`.
- `parse_kmer_votes` helper.
- `nonhuman-screen classify` CLI (requires the `[bam]` extra).
- Docs: methodology, database setup, CLI reference; `download_kraken2_db.sh`;
  Dockerfile pinning kraken2 v2.17.1; CI + PyPI-publish workflows.
