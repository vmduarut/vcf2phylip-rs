# Compatibility, provenance and future semantics

## Reference baseline

The CLI vocabulary is informed by **vcf2phylip 2.9**, Edgardo M. Ortiz,
[upstream commit 85a260e655b5f3c5dc2cb9915eb24552a002357d](https://github.com/edgardomortiz/vcf2phylip/tree/85a260e655b5f3c5dc2cb9915eb24552a002357d).
Please cite upstream when appropriate: **DOI
[10.5281/zenodo.2540861](https://doi.org/10.5281/zenodo.2540861)**.
This project is not endorsed by upstream. The Rust code and synthetic tests
were independently designed, not a translation of upstream Python. The
reference source was consulted for argument spellings and behavioral caveats.

Upstream is GPL-3.0. This project conservatively uses **GPL-3.0-only** as well;
see LICENSE. No upstream authorship is claimed for the new Rust implementation.
Dependencies retain their own respective licenses. This provenance note is not
legal advice or a declaration of permission to relicense upstream material.

## Argument vocabulary retained

| Argument | Milestone 3 status |
| --- | --- |
| `-i`, `--input` | Required input path |
| `-m`, `--min-samples-locus` | Positive integer, default 4, capped with warning; independently applied to DNA/binary |
| `-o`, `--outgroup` | Exact validated name, repeatable; stable outgroup-first output order |
| `--output-folder`, `--output-prefix` | Output directory and safe basename; defaults cwd and input basename |
| `-p`, `--phylip-disable` | Disable default relaxed PHYLIP; requires FASTA |
| `-f`, `--fasta` | Add FASTA output |
| `-n`, `--nexus` | Reserved, explicitly unsupported |
| `-b`, `--nexus-binary` | Reserved, explicitly unsupported |
| `-r`, `--resolve-IUPAC` | Seeded allele-copy resolution; explicit `--seed` required |
| `--seed` | u64, permitted only with `-r` |
| `--partial-missing` | `missing` (default) or `called-alleles`, nucleotide calls only |
| `-w`, `--write-used-sites` | Retained nucleotide coordinates and called-sample count |
| `-v`, `--version` | This Rust package's version, not upstream's version |
| `--inspect` | Diagnostics only; incompatible with output flags |

Names, errors, output text, compression detection and multiple outgroups
intentionally differ from upstream. In particular, input sample names are not
cleaned by removing `./`. Malformed records fail rather than being skipped.
There is no promise of exact Python argparse behavior, output filenames,
alignment bytes, counters, RNG sequence or performance equivalence.

## Known reference bugs are not compatibility promises

The baseline has behavior that a future scientific converter must not inherit:

1. Missingness counted from the first character can treat `0/.` and `./0`
   asymmetrically during filtering.
2. Character-based GT interpretation can split multidigit allele indices.
3. Replacing `<NON_REF>` with REF can manufacture reference support from a
   symbolic/nonreference allele.
4. Binary-only selection can diverge from nucleotide-path accepted-site
   counts and the used-site list.

Milestone 1 parses genotype indices numerically and preserves missingness and
symbolic alleles. Milestone 2 interprets canonical SNP genotypes independently,
with symmetric partial missingness and independent DNA/binary thresholds.
Milestone 3 emits nucleotide FASTA/relaxed PHYLIP and used-site output with
hand-computed synthetic tests, not upstream output-equivalence claims. See the
[scientific contract](../README.md#scientific-contract) and
[output safety contract](../README.md#conversion-and-output-safety).

## Planned, not shipped

M4 will add nucleotide NEXUS/binary NEXUS/SNAPP; M5 will add bounded disk-backed
matrix storage. M3 currently uses unbounded RAM for one sequence per sample.
Scientific inspection remains a separate no-write mode.
The independently designed implementation and synthetic tests do not copy
upstream Python algorithms or sequence fixtures.
