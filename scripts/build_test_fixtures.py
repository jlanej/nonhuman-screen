#!/usr/bin/env python3
"""Generate the committed functional-test fixtures for nonhuman-screen.

Produces, under ``tests/data/``:

* ``kraken2_mini_strep/`` — a tiny kraken2 database containing a single
  synthetic *Streptococcus pyogenes* (taxid 1314) reference, pruned to just the
  files needed at classification time (hash/opts/taxo.k2d + taxonomy dumps).
* ``bam/HG00{2,3,4}_*.bam`` — the GIAB human BAMs, copied verbatim (the
  **negative** controls: real human reads, expected non-human fraction ~0).
* ``bam/HG002_child.strep.bam`` — the child BAM with synthetic strep reads
  spliced in (the **positive** control: expected non-human fraction well above 0).
* ``strep_ref.fna`` — the synthetic strep reference (kept so the DB and the
  injected reads can be regenerated identically).

These are committed so CI only needs the ``kraken2`` binary (not
``kraken2-build``) to run the functional tests. Regenerate with:

    python scripts/build_test_fixtures.py --giab-src /path/to/kmer_denovo_filter/tests/data/giab

Requires: kraken2-build, samtools, pysam, and (for the first run) the GIAB BAMs.
The database is built with kraken2 2.17.1 to match the pinned CI version.
"""

from __future__ import annotations

import argparse
import os
import random
import shutil
import subprocess
import sys

import pysam

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "tests", "data")
DB_DIR = os.path.join(DATA, "kraken2_mini_strep")
BAM_DIR = os.path.join(DATA, "bam")
STREP_FA = os.path.join(DATA, "strep_ref.fna")

# Streptococcus pyogenes and its (minimal, but Bacteria-rooted) lineage.
STREP_TAXID = 1314
STREP_REF_LEN = 50_000
STREP_SEED = 1314
N_INJECTED_READS = 500
INJECT_READ_LEN = 150
INJECT_SEED = 20240711

GIAB_BAMS = ["HG002_child.bam", "HG003_father.bam", "HG004_mother.bam"]
POSITIVE_SRC = "HG002_child.bam"  # the BAM we inject strep reads into

# (taxid, parent, rank) — root -> Bacteria -> ... -> S. pyogenes
_NODES = [
    (1, 1, "no rank"),
    (131567, 1, "no rank"),
    (2, 131567, "superkingdom"),
    (1239, 2, "phylum"),
    (91061, 1239, "class"),
    (186826, 91061, "order"),
    (1300, 186826, "family"),
    (1301, 1300, "genus"),
    (STREP_TAXID, 1301, "species"),
]
_NAMES = {
    1: "root", 131567: "cellular organisms", 2: "Bacteria",
    1239: "Bacillota", 91061: "Bacilli", 186826: "Lactobacillales",
    1300: "Streptococcaceae", 1301: "Streptococcus",
    STREP_TAXID: "Streptococcus pyogenes",
}


def _gen_seq(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def write_strep_reference() -> str:
    seq = _gen_seq(STREP_REF_LEN, STREP_SEED)
    with open(STREP_FA, "w") as fh:
        fh.write(
            f">strep_ref|kraken:taxid|{STREP_TAXID} "
            "Streptococcus pyogenes (synthetic test sequence)\n"
        )
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i + 70] + "\n")
    print(f"[fixtures] wrote {STREP_FA} ({STREP_REF_LEN} bp, taxid {STREP_TAXID})")
    return seq


def build_database() -> None:
    if os.path.isdir(DB_DIR):
        shutil.rmtree(DB_DIR)
    tax = os.path.join(DB_DIR, "taxonomy")
    os.makedirs(tax, exist_ok=True)

    with open(os.path.join(tax, "nodes.dmp"), "w") as fh:
        for taxid, parent, rank in _NODES:
            fh.write(
                f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\t0\t|\t0"
                f"\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n"
            )
    with open(os.path.join(tax, "names.dmp"), "w") as fh:
        for taxid, name in sorted(_NAMES.items()):
            fh.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")

    # Single-threaded: kraken2-build's multi-threaded build_db pipeline
    # (find | xargs cat | build_db) is flaky on macOS (SIGPIPE); the DB is tiny
    # so threading buys nothing. OMP_NUM_THREADS=1 keeps build_db from warning.
    env = {**os.environ, "OMP_NUM_THREADS": "1"}
    subprocess.run(
        ["kraken2-build", "--add-to-library", STREP_FA, "--db", DB_DIR,
         "--no-masking"],
        check=True, capture_output=True, env=env,
    )
    subprocess.run(
        ["kraken2-build", "--build", "--db", DB_DIR, "--threads", "1",
         "--no-masking", "--kmer-len", "35", "--minimizer-len", "31"],
        check=True, capture_output=True, env=env,
    )

    # Prune to the minimal set needed at classification time:
    #   kraken2 classify -> hash.k2d, opts.k2d, taxo.k2d
    #   nonhuman-screen  -> taxonomy/nodes.dmp, taxonomy/names.dmp
    keep = {"hash.k2d", "opts.k2d", "taxo.k2d"}
    for name in os.listdir(DB_DIR):
        path = os.path.join(DB_DIR, name)
        if name in keep:
            continue
        if name == "taxonomy":
            for t in os.listdir(path):
                if t not in {"nodes.dmp", "names.dmp"}:
                    os.remove(os.path.join(path, t))
            continue
        (shutil.rmtree if os.path.isdir(path) else os.remove)(path)

    size = sum(
        os.path.getsize(os.path.join(dp, f))
        for dp, _, fs in os.walk(DB_DIR) for f in fs
    )
    print(f"[fixtures] built {DB_DIR} ({size / 1e6:.2f} MB, pruned)")


def copy_negatives(giab_src: str) -> None:
    os.makedirs(BAM_DIR, exist_ok=True)
    for bam in GIAB_BAMS:
        for ext in ("", ".bai"):
            src = os.path.join(giab_src, bam + ext)
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join(BAM_DIR, bam + ext))
    print(f"[fixtures] copied {len(GIAB_BAMS)} negative-control BAMs -> {BAM_DIR}")


def inject_positive(strep_seq: str) -> None:
    src = os.path.join(BAM_DIR, POSITIVE_SRC)
    dst = os.path.join(BAM_DIR, POSITIVE_SRC.replace(".bam", ".strep.bam"))
    rng = random.Random(INJECT_SEED)
    with pysam.AlignmentFile(src, "rb") as fin:
        n_orig = 0
        with pysam.AlignmentFile(dst, "wb", header=fin.header) as fout:
            for read in fin:
                fout.write(read)
                n_orig += 1
            for i in range(N_INJECTED_READS):
                start = rng.randint(0, len(strep_seq) - INJECT_READ_LEN)
                a = pysam.AlignedSegment(fout.header)
                a.query_name = f"strep_inject_{i:04d}"
                a.query_sequence = strep_seq[start:start + INJECT_READ_LEN]
                a.flag = 4  # unmapped
                a.reference_id = -1
                a.reference_start = -1
                a.mapping_quality = 0
                a.query_qualities = pysam.qualitystring_to_array("I" * INJECT_READ_LEN)
                fout.write(a)
    pysam.index(dst)
    print(
        f"[fixtures] wrote {dst}: {n_orig} human reads "
        f"+ {N_INJECTED_READS} injected strep reads"
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--giab-src",
        default=os.path.join(BAM_DIR),
        help="Directory holding the source GIAB BAMs (default: the already-"
             "committed tests/data/bam). Point at kmer_denovo_filter/tests/data/"
             "giab for the very first generation.",
    )
    args = ap.parse_args(argv)

    for tool in ("kraken2-build", "samtools"):
        if shutil.which(tool) is None:
            print(f"error: {tool} not found on PATH", file=sys.stderr)
            return 2

    os.makedirs(DATA, exist_ok=True)
    os.makedirs(BAM_DIR, exist_ok=True)
    strep_seq = write_strep_reference()
    if os.path.abspath(args.giab_src) != os.path.abspath(BAM_DIR):
        copy_negatives(args.giab_src)
    build_database()
    inject_positive(strep_seq)
    print("[fixtures] done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
