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

| Argument | Milestone 1 status |
| --- | --- |
| `-i`, `--input` | Required input path |
| `-m`, `--min-samples-locus` | Positive integer, default 4, capped with warning; not applied |
| `-o`, `--outgroup` | Exact validated name, repeatable; not reordered yet |
| `--output-folder`, `--output-prefix` | Reserved, explicitly unsupported |
| `-p`, `--phylip-disable` | Reserved, explicitly unsupported |
| `-f`, `--fasta` | Reserved, explicitly unsupported |
| `-n`, `--nexus` | Reserved, explicitly unsupported |
| `-b`, `--nexus-binary` | Reserved, explicitly unsupported |
| `-r`, `--resolve-IUPAC` | Reserved, explicitly unsupported |
| `-w`, `--write-used-sites` | Reserved, explicitly unsupported |
| `-v`, `--version` | This Rust package's version, not upstream's version |
| `--inspect` | New, mandatory for functional diagnostic path |

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
symbolic alleles. It has **no filtering, alignment or used-site counter** and
does not claim to have implemented or validated fixes to conversion outputs.

## Planned, not shipped

Later milestones may add genotype interpretation, IUPAC conversion, filtering,
sample ordering, matrix storage/transposition, PHYLIP/FASTA/NEXUS writers and
scientific statistics. Only `cli` and `reader` are implemented now; there is no
empty scaffolding pretending the other components exist.

The intended conservative SNP policy is canonical single-base A/C/G/T alleles,
with explicit handling/rejection of symbolic alleles, indels, spanning deletions
and ambiguous reference bases rather than silently coercing them. Partial
missingness must be order-independent: `0/.` and `./0` have identical filtering
status, with a documented policy separating unknown calls from fully called
ones. The exact filtering/conversion contract, binary eligibility, used-site
accounting and deterministic ambiguity handling still need implementation and
scientific tests. None of those future rules is applied by `--inspect`.
