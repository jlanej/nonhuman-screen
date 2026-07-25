# Command-line interface

```
nonhuman-screen classify --bam BAM --kraken2-db DIR [options]
```

Requires the `[bam]` extra (pysam) and a `kraken2` binary on `PATH`.

## Options

| Flag | Description |
|---|---|
| `--bam` | Input BAM/CRAM (required). |
| `--kraken2-db` | kraken2 database directory (required). |
| `--ref-fasta` | Reference FASTA (required for CRAM). |
| `--variants` | VCF/BCF. When given, compute the allele-based NHF for every concrete ALT allele instead of a whole-BAM summary. |
| `--min-baseq` | Minimum base quality over the variant span for a **read** to count as ALT support (default: 0). Applied per read, not per base: if *any* base the read contributes across the variant's reference span is below the threshold, the whole read is dropped from ALT support. |
| `--out-prefix` | Write outputs to files instead of stdout. With `--variants`: `<prefix>.variant_nhf.tsv` **and** `<prefix>.summary.json`. Without `--variants` (whole-BAM mode): only `<prefix>.summary.json`. |
| `--threads` | kraken2 threads (default: 1). |
| `--confidence` | kraken2 confidence threshold 0.0–1.0 (default: 0.0). |
| `--memory-mapping` | Pass `--memory-mapping` to kraken2 to reduce RAM. |

Symbolic and breakend ALT alleles (`<DEL>`, `N[chr2:321[`, `*`) are skipped.
VCF positions are converted to the 0-based internal convention automatically.

## Examples

Per-variant allele NHF for candidate calls:

```bash
nonhuman-screen classify \
    --bam sample.bam --kraken2-db kraken2_db --ref-fasta ref.fa \
    --variants candidates.vcf.gz --out-prefix sample_contam
```

Produces `sample_contam.variant_nhf.tsv`:

```
variant_key       supporting_reads  nonhuman_fraction  bacterial  archaeal  ...
chr1:12344:A:T    18                0.0                0.0        0.0       ...
chr8:40119:G:GAC  12                0.83               0.83       0.0       ...
```

and `sample_contam.summary.json` (the full per-variant breakdown including
per-domain fractions and supporting-read counts).

Whole-BAM summary:

```bash
nonhuman-screen classify --bam sample.bam --kraken2-db kraken2_db
```

emits a JSON summary with `reads`, `classified`, `taxonomy_available`,
`nonhuman_fraction`, and the per-domain `fractions`.

Every record in the file is scanned, unmapped ones included, but reads are
de-duplicated by name: **one sequence per unique read name**, first occurrence
wins. So `reads` counts *fragments*, not records — paired mates and
secondary/supplementary alignments of the same fragment collapse to one. This
is the only mode that reports `taxonomy_available` (see
[methodology.md §8](methodology.md)).

> **Memory:** whole-BAM mode holds one sequence per fragment in memory before
> invoking kraken2, so it suits targeted or small BAMs. For a whole-genome BAM
> use `--variants` (which classifies only ALT-supporting reads), or drive the
> API with your own read selection.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Success. |
| 2 | `kraken2` not found on `PATH`. |
| 3 | kraken2 classification failed — non-zero exit, or per-read output that did not account for every input read. **No summary file is written**, because a failed screen must never be recorded as `nonhuman_fraction = 0.0`. See [methodology.md §7](methodology.md). |
