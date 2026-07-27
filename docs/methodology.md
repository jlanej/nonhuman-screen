# Methodology

`nonhuman-screen` classifies reads by non-human taxonomic content using
kraken2's lowest-common-ancestor (LCA) algorithm and reduces the per-read
verdicts to a conservative **non-human fraction (NHF)** plus per-domain
breakdowns.

## 1. Classification

Reads are written to a temporary FASTQ (with dummy qualities) and classified by
the `kraken2` binary. kraken2 breaks each read into k-mers, looks up the LCA of
each k-mer against the database, and assigns the read a single taxid.

The `--confidence` threshold (0.0–1.0, default 0.0) controls how strict the LCA
must be: a value of 0.2 requires ≥20% of a read's k-mers to agree on the
assigned clade. Higher values reduce spurious assignments at the cost of
sensitivity.

## 2. Lineage-aware domain assignment

The database's `nodes.dmp` taxonomy is parsed into a child→parent map, and
descendant sets are computed for each domain from its NCBI root taxid:

| Domain | Root taxid | Definition |
|---|---|---|
| Bacteria | 2 | all descendants |
| Archaea | 2157 | all descendants |
| Fungi | 4751 | all descendants |
| Protist | — | `eukaryota(2759) − metazoa(33208) − fungi − viridiplantae(33090)` |
| Viruses | 10239 | all descendants |
| UniVec-Core | 81077 | synthetic vector/adapter sequences (see §4) |
| Human clade | 9606 | *Homo sapiens* and any sub-taxa |
| Human lineage | — | ancestors of 9606 up to root (Eukaryota, Root, …) |

A read is **non-human** iff it is classified **and** its taxid is *not* on the
human→root lineage, *not* in the human clade, and *not* under UniVec-Core.

**This test is negative, so its meaning depends on the database.** Anything
classified that is not explicitly excluded counts as non-human — including taxa
outside the five reported domains, such as `other sequences` (28384, the parent
of UniVec-Core) or `unclassified sequences` (12908). An LCA can only land on
such a node if the database holds sequence under it.

For [PrackenDB](database.md) it follows from the documented composition
(bacteria, archaea, protists, fungi, human, RefSeq viral, UniVec-Core) that
nothing sits under 28384 except UniVec-Core itself and nothing sits under
12908, so every LCA that could fall outside the five domains instead lands on
`cellular organisms` (131567), `root` (1), or `Eukaryota` (2759) — all on the
human lineage, all excluded. On that reasoning the PrackenDB non-human fraction
*is* the union of the five domains, with no unlabelled residual. It is an
inference from the composition, not a measurement: to confirm it for a given
database, check whether `nonhuman` exceeds the sum of the five domain columns on
a real run. Do that before treating the two definitions as interchangeable on
any other database (`nt`, custom builds, anything with siblings under 28384 or
12908).

**The database's taxonomy must contain the human node (taxid 9606).** Without
it both `human_lineage` and `human_clade` resolve to empty sets, so *every*
classified taxid — root included — passes the non-human test, and the
human-homology guard below can never fire. The engine logs a warning when it
loads such a taxonomy; non-human fractions from it are not comparable to
fractions from a human-containing database.

## 3. Human-homology guard (HHG)

Some non-human references share k-mers with the human genome (notably
integrating viruses — endogenous retroviruses, HBV, HPV). To avoid over-calling,
the per-read k-mer detail string is inspected: **if any k-mer voted for human
(taxid 9606), the read is removed from every non-human numerator**, regardless
of its LCA assignment. Such reads are marked with guard status `HHG`.

## 4. UniVec-Core exclusion

UniVec-Core (taxid 81077) is a curated set of synthetic sequencing-vector and
adapter sequences. These are laboratory artefacts, not biological contamination,
and can share k-mers with human DNA. Reads assigned to UniVec-Core are tracked
separately (`univec_core_read_names`) and **unconditionally excluded** from the
non-human fraction. This is a second safety net applied even when the HHG does
not fire.

## 5. The non-human fraction and its partition

Over any set of reads, the classified+unclassified reads are partitioned into
four disjoint sets whose fractions sum to 1.0 (modulo rounding):

```
nonhuman + univec_core + human_lineage + unclassified = 1.0
```

- **nonhuman** — definitively outside the human lineage (after HHG + UniVec).
- **univec_core** — synthetic vectors.
- **human_lineage** — human clade, HHG-cleared reads, and reads assigned to
  ranks too broad to call (Root, Eukaryota, …).
- **unclassified** — kraken2 returned no assignment.

### Relationship between the domain columns and `nonhuman`

The per-domain breakdowns (`bacterial`, `archaeal`, `fungal`, `protist`,
`viral`) are **not** a partition of `nonhuman`, and should not be read as one:

- `bacterial`, `archaeal`, `fungal` and `viral` are subsets of `nonhuman`.
- **`protist` is not.** It is computed as
  `eukaryota − metazoa − fungi − viridiplantae`, and that set still contains
  `Eukaryota` (2759) and `Opisthokonta` (33154) — both human ancestors, both
  correctly excluded from `nonhuman`. A read whose LCA is Eukaryota (a routine
  outcome when a database holds both fungi and protists) is therefore reported
  as `protist > 0` while `nonhuman = 0`, and is simultaneously counted in
  `human_lineage`. The `protist` column over-reports; it never causes
  `nonhuman` to over- or under-report.
- Whether the domains *sum* to `nonhuman` is a property of the database, not of
  the method — see the note in §2. They do for PrackenDB.

Only the four-way partition above is guaranteed. Use `nonhuman` for decisions
and the domain columns for attribution.

## 6. Allele-based NHF

For variant-level screening, the reads supporting a variant's **ALT allele** are
identified by `read_supports_alt`, which extracts the read sequence aligned to
the variant's reference span (via pysam's aligned pairs) and compares it
strictly to the ALT — handling SNPs, MNPs, insertions, deletions, and complex
indels. Only ALT-supporting reads are classified, and the NHF is computed over
that set. When many variants are screened, the union of their ALT-supporting
reads is classified in a single kraken2 invocation and split back per variant.

**Caveat — whole-read scope:** kraken2 classifies the *entire* read, not the
sub-region overlapping the variant. A read is judged non-human as a whole; the
allele-based NHF therefore measures "how many ALT-supporting reads are non-human
reads," which is the intended signal for contamination screening but should not
be read as locus-level taxonomy.

## 7. Failure semantics: a broken screen must not look like a clean one

A *failed run* fails towards **understating** non-human content — a false
negative that reads as a clean sample — so failures raise rather than return a
number. (Note that the two *degraded-database* modes go the other way: a
missing `nodes.dmp` (§8) and a taxonomy without the human node (§2) both make
the non-human fraction over-count. Those still return a usable tally, so they
warn rather than raise.)

**A failed kraken2 run raises.** `Kraken2Runner.classify_sequences` raises
`Kraken2Error` when the subprocess exits non-zero, or when its per-read output
does not account for every input read. The second case matters as much as the
first: unaccounted reads still count towards the per-variant denominator but
towards no numerator, so a truncated output silently deflates every fraction.
Both would otherwise surface as `nonhuman_fraction = 0.0`.

Pass `strict=False` to downgrade both to a warning plus
`ClassificationResult.classification_failed = True`. If you do, you must check
that flag — an unchecked flag is the failure mode this design exists to
prevent.

The CLI is always strict and exits **3** on failure, writing no summary file.

**Zero support is not zero contamination.** `nonhuman_fraction == 0.0` with
`supporting_reads == 0` (or `total == 0`) means *no evidence*, not *clean*.
Always read the count alongside the fraction.

## 8. When taxonomy is unavailable

If `nodes.dmp` cannot be read, classification falls back to **exact-taxid
matching only** (no lineage walk), and the results become unreliable in *both*
directions:

- The **per-domain breakdowns** (`bacterial`, `viral`, …) **under-count**: only
  reads assigned to the exact domain-root taxid match, and `protist` cannot be
  computed at all, so all descendant taxa are missed.
- The **consolidated non-human fraction over-counts**: the fallback treats a
  read as non-human unless its taxid is exactly human (9606), root (1), or
  UniVec-Core, so human-lineage ancestors (e.g. genus *Homo*, Eukaryota) and
  non-9606 human subspecies fall through and are mis-counted as non-human.

In this mode the engine logs a warning and sets
`ClassificationResult.taxonomy_available = False`. **Consumers that gate
decisions on non-human content should treat `taxonomy_available is False` as
"unknown" — neither "clean" nor "contaminated"** — since a database missing its
taxonomy dumps otherwise silently corrupts the signal. Always verify your
database includes `nodes.dmp` (and `names.dmp` for taxon names); see
[database.md](database.md).

### Where `taxonomy_available` is observable

Everywhere. It originates on the engine `ClassificationResult` and is copied
onto each `VariantNHF`, so it is reachable from:

- `Kraken2Runner.classify_sequences(...)` and `classify_reads_from_bam(...)` —
  on the returned result.
- `classify_variant_alt_reads(...)` / `classify_variants_alt_reads(...)` — as
  `VariantNHF.taxonomy_available`.
- `nonhuman-screen classify` — in `summary.json` in whole-BAM mode, and with
  `--variants` as the last column of `<prefix>.variant_nhf.tsv` plus a field
  per entry in `<prefix>.summary.json`.

**Why it is worth gating on.** Unlike a failed run (§7), a taxonomy-less
database still produces a plausible-looking tally, so it cannot raise — but
that tally *over*-counts non-human content. A pipeline that down-ranks variants
above some NHF threshold will therefore reject real calls, with nothing in the
numbers to indicate why. Checking this column is what turns that into a
detectable condition.

One kraken2 invocation loads the taxonomy once, so the value is the same on
every variant of a batch. It describes the classification run, which means it
is only meaningful once a run happened: a batch in which no variant had any
ALT-supporting read never invokes kraken2 and reports the default `True`
alongside `supporting_reads == 0`. Read it together with the read count.
