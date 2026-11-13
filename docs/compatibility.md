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

| Argument | Current status |
| --- | --- |
| `-i`, `--input` | Required input path |
| `-m`, `--min-samples-locus` | Optional absolute positive integer minimum, capped with warning; mutually exclusive with explicit call rate |
| `--min-call-rate` | Finite fraction in [0,1], default 0.9 (90%); without `-m`, requires max(1, ceil(rate × actual input sample count)) |
| `-o`, `--outgroup` | Exact validated name, repeatable; stable outgroup-first output order |
| `--output-folder`, `--output-prefix` | Output directory and safe basename; defaults cwd and input basename |
| `-p`, `--phylip-disable` | Disable default relaxed PHYLIP; requires FASTA |
| `-f`, `--fasta` | Add FASTA output |
| `-r`, `--resolve-IUPAC` | Seeded allele-copy resolution; explicit `--seed` required |
| `--seed` | u64, permitted only with `-r` |
| `--partial-missing` | `missing` (default) or `called-alleles`, nucleotide calls only |
| `-w`, `--write-used-sites` | Retained nucleotide coordinates and called-sample count |
| `-v`, `--version` | This Rust package's version, not upstream's version |
| `--inspect` | Diagnostics only; incompatible with output flags |
| `--storage` | `memory` (unbounded default) or `disk` (bounded cell buffer) |
| `--memory-limit` | Disk-only cell-vector capacity budget; default 64 MiB; bytes/B/KiB/MiB/GiB |
| `--temp-dir` | Disk-only existing matrix spool directory; default OS temp |

Names, errors, output text, compression detection and multiple outgroups
intentionally differ from upstream. In particular, input sample names are not
cleaned by removing `./`. Malformed records fail rather than being skipped.
There is no promise of exact Python argparse behavior, output filenames,
alignment bytes, counters, RNG sequence or performance equivalence.

The minimum-name spelling is retained, but filtering defaults intentionally
change from an absolute four samples to 90% of input samples with an absolute
floor of one when `-m` is absent. Explicit `-m` and explicit `--min-call-rate`
are mutually exclusive, even if the rate is 0.9, in either argument order.
`-m` alone selects count-only filtering; the implicit default rate does not apply.
Use `-m 2` for count-only filtering, or `-m 4` for this
project's prior default. Without `-m`, rate 1 requires all samples and rate 0
requires at least one; outgroups remain part of the input denominator.
The effective threshold applies to nucleotide-called sample counts. Output `.min<N>` uses
the effective count. Reports show only the active policy. Rates use finite
`f64` values, with decimal
integer threshold arithmetic to avoid floating multiplication boundary errors.
Public `Cli::min_samples_locus` is `Option<usize>` and `Cli::min_call_rate` changes
from `f64` to `Option<f64>`. Use `Some(n)`/`None` for count-only filtering,
`None`/`Some(rate)` for an explicit rate, or `None`/`None` for the default 0.9 rate.
Supplying both as `Some` fails public validation. The low-level
`VariantProcessor::new` count-only API retains its legacy `None` behavior.

## Known reference bugs are not compatibility promises

The baseline has behavior that a future scientific converter must not inherit:

1. Missingness counted from the first character can treat `0/.` and `./0`
   asymmetrically during filtering.
2. Character-based GT interpretation can split multidigit allele indices.
3. Replacing `<NON_REF>` with REF can manufacture reference support from a
   symbolic/nonreference allele.

Milestone 1 parses genotype indices numerically and preserves missingness and
symbolic alleles. Milestone 2 interprets canonical SNP genotypes independently,
with symmetric partial missingness and nucleotide called-sample thresholds.
Milestone 3 emits nucleotide FASTA/relaxed PHYLIP and used-site output with
hand-computed synthetic tests, not upstream output-equivalence claims. See the
[scientific contract](../README.md#scientific-contract) and
[output safety contract](../README.md#conversion-and-output-safety).

## Deferred and current behavior

Optional bounded-cell-buffer disk-backed storage is available; the default
memory backend remains unbounded in loci. The same processing and writer paths
produce identical backend outputs; this is not upstream byte compatibility.
Scientific inspection remains a separate no-write mode.
The independently designed implementation and synthetic tests do not copy
upstream Python algorithms or sequence fixtures.

## Nucleotide-only API change

This is an intentional public library API break, including for struct literals:

- Removed `genotype::BinaryCall` (including its `symbol` method) and
  `genotype::binary_call`.
- Removed `processing::ProcessedLocus::{binary,binary_called_samples,binary_accepted}`.
- Removed `processing::Stats::{binary_eligible,binary_sites_accepted,binary_missing_calls,binary_unsupported_ploidy_calls}`.

Callers should use `NucleotideCall`, `nucleotide_call`, the nucleotide fields of
`ProcessedLocus`, and the remaining nucleotide/record counters in `Stats`.
No replacement dosage encoding is provided. Nucleotide eligibility, thresholds,
partial missingness, arbitrary positive ploidy, IUPAC and seeded copy-weighted
resolution are unchanged, as are FASTA, relaxed PHYLIP and used-site bytes.

Inspection stdout now contains only nucleotide and input diagnostics: auxiliary
accepted/eligible/missing/unsupported-ploidy counts are gone. Its threshold line
is `Minimum called samples: N (effective)`, the acceptance line is
`Accepted nucleotide loci: N; rejected non-SNP: N`, and the missing-call line is
`Nucleotide missing calls: N`. Consumers of inspection text must update their
parsers. CLI help likewise describes only nucleotide calls; unsupported options
remain unrecognized rather than being accepted or silently ignored.
