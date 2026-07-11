# Testing

`nonhuman-screen` has three layers of tests:

| Layer | Files | Needs | What it covers |
|---|---|---|---|
| **Unit** | `tests/test_engine.py`, `test_alleles.py`, `test_bam.py` | nothing (kraken2 is mocked) | classification logic, allele support, fraction math |
| **Integration** | `tests/test_integration.py` | `kraken2` + `kraken2-build` | builds a tiny DB at runtime and checks `Kraken2Runner` classifies bacterial/human/unclassified reads |
| **Functional (CLI controls)** | `tests/test_functional_cli.py` | `kraken2` binary only | runs the real `nonhuman-screen classify` **command** end-to-end on real BAMs, using committed fixtures |

This document describes the **functional positive/negative controls**: their
rationale, how the fixtures are created, and how the tests assert on them.

---

## 1. Rationale — why positive and negative controls

The unit and integration tests exercise the Python API. The functional tests
answer a different question: **does the shipped command-line tool actually detect
non-human contamination on a real BAM, and stay quiet when there is none?**

They do that with two matched controls that differ only by the presence of
non-human reads:

- **Negative control** — the unmodified GIAB human BAMs. A correct tool must
  report **no** non-human content on clean human data (no false positives).
- **Positive control** — the same child BAM with a known number of bacterial
  (*Streptococcus*) reads spliced in. A correct tool must **detect** them (no
  false negatives).

Because the two inputs are identical except for the injected reads, a passing
run proves the signal comes from the contamination, not from anything else in
the pipeline.

---

## 2. The committed fixtures

Everything the functional tests need is generated once and committed under
`tests/data/`, so CI needs only the `kraken2` **binary** (not `kraken2-build`)
and never rebuilds the database or the BAMs:

```
tests/data/
├── strep_ref.fna                  # synthetic Streptococcus reference (50 kb)
├── kraken2_mini_strep/            # the mini kraken2 database (~134 KB)
│   ├── hash.k2d  opts.k2d  taxo.k2d
│   └── taxonomy/{nodes.dmp, names.dmp}
└── bam/
    ├── HG002_child.bam(.bai)      # negative control (GIAB human reads)
    ├── HG003_father.bam(.bai)     # negative control
    ├── HG004_mother.bam(.bai)     # negative control
    └── HG002_child.strep.bam(.bai)# positive control (child + injected strep)
```

These files live only in the git repository — they are **excluded from the
built sdist/wheel**, so they don't bloat the PyPI release (the package only
ships `src/`).

---

## 3. Creation method — `scripts/build_test_fixtures.py`

All fixtures are produced deterministically by one script. Regenerate with:

```bash
# First generation (needs the source human BAMs):
python scripts/build_test_fixtures.py --giab-src /path/to/kmer_denovo_filter/tests/data/giab
# Re-generating the DB / injected BAM from the already-committed negatives:
python scripts/build_test_fixtures.py
```

### 3.1 The synthetic *Streptococcus* reference

`strep_ref.fna` is a **deterministic pseudo-random 50 kb DNA sequence**
(`random.Random(1314)`), written with a kraken2 taxid tag in its header:

```
>strep_ref|kraken:taxid|1314 Streptococcus pyogenes (synthetic test sequence)
```

It is deliberately **synthetic rather than a real *S. pyogenes* genome**:

- **Deterministic & self-contained** — no network download, reproducible byte-
  for-byte from the seed.
- **No accidental human homology** — a real bacterial genome can share the odd
  conserved 35-mer with the human reference, which would muddy the controls. A
  random sequence shares no k-mers with the human reads, so the negative control
  is *exactly* zero.

The sequence carries the **real NCBI taxid 1314 (*Streptococcus pyogenes*)**, so
kraken2 assigns matching reads to that species and `nonhuman-screen` resolves it
— via the taxonomy — as bacterial / non-human. In other words: real taxonomy,
synthetic sequence.

### 3.2 The minimal taxonomy

`taxonomy/nodes.dmp` and `names.dmp` encode a small lineage that places 1314
under Bacteria (taxid 2):

```
1314 Streptococcus pyogenes → 1301 Streptococcus → 1300 Streptococcaceae →
186826 Lactobacillales → 91061 Bacilli → 1239 Bacillota → 2 Bacteria →
131567 cellular organisms → 1 root
```

This is **minimal but Bacteria-rooted**, not the full NCBI lineage (e.g. it
attaches Bacillota directly to Bacteria, skipping the Terrabacteria group). That
is sufficient: `nonhuman-screen` computes the bacterial domain as *all
descendants of taxid 2*, and 1314's lineage passes through 2, so it is counted
as bacterial. (Verified: `Kraken2Runner._load_all_taxid_sets(db)['bacterial']`
contains 1314.) No human/eukaryote nodes are included — see §4.

### 3.3 Building and pruning the database

The DB is built with the same recipe the integration test uses, at the kraken2
default **k = 35** (minimizer 31):

```
kraken2-build --add-to-library strep_ref.fna --db <db> --no-masking
kraken2-build --build --db <db> --no-masking --kmer-len 35 --minimizer-len 31 --threads 1
```

- `--no-masking` avoids the `dustmasker` (BLAST+) dependency.
- `--threads 1` **plus** `OMP_NUM_THREADS=1` are required on macOS: the
  multi-threaded `build_db` pipeline (`find | xargs cat | build_db`) raises a
  SIGPIPE there. The DB is tiny, so single-threaded costs nothing. Linux CI is
  unaffected either way.

Only the files needed at classification time are kept (the rest of the
`kraken2-build` scratch — `library/`, `seqid2taxid.map`, intermediate maps — is
pruned): `hash.k2d`, `opts.k2d`, `taxo.k2d` (used by the `kraken2` classifier)
and `taxonomy/nodes.dmp`, `taxonomy/names.dmp` (used by `nonhuman-screen` for
lineage-aware domain assignment and taxon names). Result: ~134 KB.

The DB is built with **kraken2 2.17.1** to match the version CI installs, so
there is no on-disk format skew. The `.k2d` format is little-endian and
platform-independent — a database built on macOS/arm64 is read correctly by
kraken2 on linux/x86_64 (confirmed in CI).

### 3.4 Injecting the positive control

`HG002_child.strep.bam` is `HG002_child.bam` with **500 synthetic strep reads**
appended (`random.Random(20240711)`):

- Each injected read is a **150 bp window** sampled from `strep_ref.fna`, so its
  k-mers all match the database and it classifies confidently as taxid 1314.
- The reads are added as **unmapped** records (SAM flag 4, `reference_id = -1`)
  with names `strep_inject_0000…0499`. Unmapped reads sort after all mapped
  reads, so the BAM stays coordinate-sorted and is re-indexed.
- Unmapped is deliberate: the whole-BAM `classify` mode reads every record
  (`fetch(until_eof=True)`), so the reads are picked up without needing to
  fabricate an alignment (see §6 for the scope this implies).

---

## 4. How the controls actually work

The database contains **only** the strep reference — no human sequence. So:

- **Negative control (human BAM):** none of the real human reads share k-mers
  with the synthetic strep sequence, so kraken2 leaves them **unclassified**
  (taxid 0). Unclassified reads are not counted as non-human → non-human
  fraction = **0.0** (unclassified fraction = 1.0).
- **Positive control (injected BAM):** the 500 injected reads are pure strep
  sequence → classified as taxid 1314 → a descendant of Bacteria (2) → counted
  as bacterial and non-human.

Note the CLI's whole-BAM mode counts **unique read (fragment) names**, so paired
mates and secondary/supplementary alignments collapse to one fragment. This is
why the record counts and the reported `reads` differ.

Observed values (from the committed fixtures, kraken2 2.17.1):

| Control | BAM records | `reads` (fragments) | `classified` | `nonhuman_fraction` | `bacterial` |
|---|---:|---:|---:|---:|---:|
| Negative (`HG002_child.bam`) | 11,097 | 5,924 | 0 | **0.0000** | 0.0000 |
| Positive (`HG002_child.strep.bam`) | 11,597 | 6,424 | 500 | **0.0778** | 0.0778 |

(`0.0778 = 500 / 6,424`.) The taxonomy loaded correctly in both
(`taxonomy_available = true`).

---

## 5. Testing method — `tests/test_functional_cli.py`

The tests invoke the **actual CLI** (not the Python API) as a subprocess —
preferring the installed `nonhuman-screen` console script, falling back to
`python -m nonhuman_screen.cli` when the bin dir isn't on `PATH`:

```
nonhuman-screen classify --bam <fixture> --kraken2-db tests/data/kraken2_mini_strep --out-prefix <tmp>/out
```

then parse `<tmp>/out.summary.json`. Assertions:

- **Negatives** (parametrized over all three human BAMs):
  `taxonomy_available is True`, `classified == 0`, `nonhuman_fraction == 0.0`,
  `bacterial == 0.0`.
- **Positive:** `classified >= 400` (~500 injected), `nonhuman_fraction > 0.02`,
  `bacterial > 0.02`, and `bacterial == nonhuman_fraction` (the injected reads
  are bacterial, nothing else).
- **Separation:** the positive's `nonhuman_fraction` exceeds every negative's.

**Why these thresholds are robust:** the negatives are *exactly* 0.0 (not merely
small), because the DB shares no k-mers with human reads. Any positive margin
therefore separates the classes cleanly, so `> 0.02` (against an observed
~0.078) and `>= 400` (against 500) leave a wide safety band and are not brittle.

The whole module is `skipif`-gated on the `kraken2` binary and the committed DB
being present, so it silently skips in environments without kraken2 (e.g. the
plain unit-test matrix) and runs for real where kraken2 is installed.

---

## 6. CI

The `integration` job in `.github/workflows/ci.yml`:

1. installs **kraken2 2.17.1** from source (matching the DB's builder version),
2. runs `tests/test_integration.py` (which builds its own DB), and
3. runs `tests/test_functional_cli.py` against the **committed** DB + BAMs.

Because the fixtures are committed, the functional job needs only the classifier
binary, and the positive/negative controls run on every push.

---

## 7. Scope and caveats

- **Synthetic, not a real genome.** The strep reference is a random sequence
  assigned the real taxid 1314. It proves the *pipeline* (kraken2 → taxonomy →
  non-human fraction) detects reads assigned to a bacterial taxon; it is not a
  test of kraken2's ability to recognise a real *S. pyogenes* genome.
- **Whole-BAM mode only.** The injected reads are unmapped, so these controls
  exercise `classify` (whole-BAM). They do **not** exercise the allele-based
  `--variants` mode, which would need mapped reads supporting an ALT allele at a
  locus. That is a natural future addition.
- **DB is version-coupled.** The committed `.k2d` files are tied to the kraken2
  on-disk format; CI pins 2.17.1 to match. If you bump kraken2, rebuild the
  fixtures with the matching version (`scripts/build_test_fixtures.py`).
