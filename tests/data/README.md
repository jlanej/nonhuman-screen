# Functional-test fixtures

Committed inputs for the end-to-end CLI controls in
[`tests/test_functional_cli.py`](../test_functional_cli.py). See
[`docs/testing.md`](../../docs/testing.md) for the full rationale, creation
method, and assertions — this file is just an inventory.

| Path | What it is |
|---|---|
| `strep_ref.fna` | Synthetic 50 kb *Streptococcus pyogenes* (taxid 1314) reference (deterministic; assigned the real taxid, random sequence). |
| `kraken2_mini_strep/` | Mini kraken2 database built from `strep_ref.fna` with kraken2 2.17.1 (k=35), pruned to `hash/opts/taxo.k2d` + `taxonomy/{nodes,names}.dmp` (~134 KB). |
| `bam/HG002_child.bam`, `HG003_father.bam`, `HG004_mother.bam` | **Negative controls** — GIAB human BAMs; expected non-human fraction `0.0`. |
| `bam/HG002_child.strep.bam` | **Positive control** — `HG002_child.bam` with 500 synthetic strep reads spliced in (unmapped); expected non-human fraction `> 0`. |

Regenerate everything deterministically:

```bash
python scripts/build_test_fixtures.py            # rebuild DB + positive from the committed negatives
python scripts/build_test_fixtures.py --giab-src /path/to/giab   # also (re)copy the source human BAMs
```

These files are **git-only**: they are excluded from the built sdist/wheel.
