# kraken2 binary and database setup

`nonhuman-screen` shells out to the `kraken2` binary and reads a kraken2
database directory. Neither is bundled.

## The kraken2 binary

Install kraken2 and make sure `kraken2` is on `PATH`. This package is developed
and tested against **kraken2 v2.17.1** (pinned in the package `Dockerfile`).

```bash
git clone --depth 1 --branch v2.17.1 https://github.com/DerrickWood/kraken2.git
cd kraken2 && ./install_kraken2.sh /usr/local/bin
kraken2 --version
```

Or use the package's Docker image, which installs the pinned version.

## The database

### Requirement: the database must contain the human genome

This is a methodological requirement, not a convenience. Both of the method's
conservatism mechanisms are defined relative to human:

- the **human-lineage / human-clade exclusion** needs taxid 9606 present in
  `nodes.dmp`, or it resolves to an empty set and *every* classified read —
  root included — is counted as non-human;
- the **human-homology guard** needs human *sequence* in the index, or no
  k-mer can ever vote for 9606 and the guard never fires.

A microbial-only or otherwise custom database therefore produces non-human
fractions that are not comparable to those from a human-containing database.
The engine logs a warning when it loads a taxonomy with no human node, but it
cannot detect a missing human *genome* in the index — verify that yourself.
PrackenDB (below) includes the human genome and satisfies both.

### Files

A kraken2 database directory must contain the index files `hash.k2d`,
`opts.k2d`, and `taxo.k2d`. It should **also** contain the taxonomy dumps
`nodes.dmp` and `names.dmp` (either under `taxonomy/` or at the database root):
these are strongly recommended but technically optional — without them the
engine still runs but falls back to exact-taxid matching, which corrupts the
non-human signal in both directions (see [methodology.md §8](methodology.md)).
`nodes.dmp` drives lineage-aware classification and `names.dmp` supplies taxon
names.

### PrackenDB (recommended)

PrackenDB is a curated, pre-built kraken2 database of NCBI reference assemblies
across all domains, built with kraken2's default k-mer length (**35**) and
shipping the taxonomy dumps. Download and validate it with the bundled script:

```bash
scripts/download_kraken2_db.sh --db /path/to/kraken_db
# override the source with --url <tarball>
```

The script downloads the tarball (`k2_NCBI_reference_20251007.tar.gz` by
default), extracts it, and validates that `hash.k2d`, `opts.k2d`, and
`taxo.k2d` are present — failing if any is missing, and warning separately if
`nodes.dmp` or `names.dmp` is absent.

Some pre-built tarballs extract into a **versioned subdirectory** (e.g.
`k2_NCBI_reference_20251007/`) rather than into `--db` itself. The script
detects that and prints the directory that actually holds the index; pass
**that** path to `--kraken2-db`, since kraken2 resolves `--db` to the directory
containing `hash.k2d` and `nonhuman-screen` looks for `nodes.dmp` there or in
its `taxonomy/` subdirectory.

### Confirming the database k-mer length

```python
from nonhuman_screen import Kraken2Runner
print(Kraken2Runner.read_kmer_length("/path/to/kraken_db"))  # e.g. 35
```

This reads the `k` field from the database's `opts.k2d` (an `IndexOptions`
struct). It is a quick way to confirm the database is intact and which k-mer
size kraken2 will use.

## Memory

Large databases are memory-hungry because kraken2 loads the hash table into
RAM. Pass `memory_mapping=True` (Python) or `--memory-mapping` (CLI) to memory-map
the index files instead, trading speed for a much smaller resident footprint.
