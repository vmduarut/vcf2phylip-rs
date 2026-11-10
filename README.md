# vcf2phylip-rs — Milestone 2

An independently designed Rust library and thin CLI for **streaming VCF/BCF
inspection and scientific genotype processing**. Alignment conversion is **not implemented**. This is not yet a
replacement for vcf2phylip, nor a claim of scientific-output compatibility.

## Build and inspect

```sh
cargo build --locked
cargo run --locked -- --input variants.vcf.gz --inspect
cargo run --locked -- -i variants.bcf --inspect -m 4 -o 'exact sample name'
cargo run --locked -- -i variants.vcf --inspect --partial-missing called-alleles -r --seed 42
cargo run --locked -- --help
```

Tested with Rust/Cargo **1.99.0** on Linux; `rust-version = "1.99"` is the
conservative tested minimum, not a claim that earlier toolchains cannot work.
Direct dependencies are exactly pinned and `Cargo.lock` is included. All
requested versions resolved: noodles 0.119.0, clap 4.6.7, thiserror 2.0.21,
anyhow 1.0.104, flate2 1.1.11, and dev dependency tempfile 3.27.0. No native
libraries or toolchain installation is needed. Seeded resolution uses exactly
pinned rand_chacha 0.10.0 (rand_core 0.10.1 in the lockfile), without the rand
convenience API. There is no application parallelism; noodles brings rayon
transitively, and tempfile brings its own test-only randomness dependencies.

## Shipped behavior

- Content-based detection (extensions are ignored): plain VCF, gzip VCF,
  concatenated/multimember gzip VCF, BGZF VCF, raw or gzip/BGZF BCF 2.2.
- Reads one record at a time with noodles adapters; memory scales with header
  and current record, not the number of loci. No random access/index required.
- Requires a valid header, at least one sample, and `FORMAT/GT` declared as
  `Number=1,Type=String`. Empty/duplicate sample names are errors. Names remain
  unchanged, including `./`, punctuation and spaces.
- Header-only input is valid inspection input. FORMAT can omit GT at a locus;
  this is carried as absent GT rather than assumed reference. A missing VCF
  GT value `.` has unknown ploidy; explicit `./.` retains two missing alleles.
  BCF's decoded allele representation is preserved (a missing encoded allele
  can be a one-allele genotype). Diagnostic absent-GT counts count `None`
  values, **not** all missing/partially missing genotypes.
- Present genotypes preserve integer indices, allele missingness, ploidy and
  per-allele phase. Haploid, diploid, polyploid and multidigit indices are
  supported. Out-of-range allele indices are errors, not reference calls.
- Strict failure on malformed records, wrong sample/field counts, invalid GT,
  invalid FORMAT layout or damaged compressed input. No unsafe skip mode.
  Errors include record ordinal/successful count and contig/position when
  available. Success counts are printed only after the stream is drained.
- The default minimum is 4, capped to the sample count with a warning; positive
  user values are similarly capped, zero/negative values rejected. This is
  applied independently to nucleotide-called and binary-encodable sample counts.
- `-o/--outgroup` requires an exact name; repeat it for multiple names, in
  requested order. Unknown or duplicate requests are errors. No comma/semicolon
  splitting, renaming, or actual sample reordering occurs during inspection.
- Inspection prints diagnostic counts, never creates output directories,
  matrices or used-site files. Without `--inspect`, conversion fails clearly.
  Reserved output flags also fail with `--inspect`, so they cannot silently be
  ignored. `--help`, `--version` and upstream-style `-v` work without input.

This is structural/genotype validation, not a complete VCF specification
validator or scientific quality-control tool. It does not impose sorted loci,
canonical nucleotides or a SNP-only policy in the parser. The scientific layer
rejects unsupported biological alleles without treating them as malformed.
Symbolic alleles and indels remain preserved in reader records. For untrusted input there are no
explicit header/record size limits. BGZF is consumed as concatenated gzip;
gzip checksums/truncated members are checked, but a missing optional BGZF EOF
marker alone is not an error and BGZF indexing metadata is not validated.

## Library

`reader::VariantReader::{open,from_reader,next_variant}` exposes the parsed
header, detected format, compression status, successful record count and a
`Variant` containing a noodles record plus optional genotype buffers in sample
order. An error terminates the reader; callers must handle that error, not
interpret a subsequent `None` as successful validation.
`filter::ProcessingConfig::new` validates sample dimensions, positive minimum
(default 4, capped), partial policy and allele resolution. Its `process` method
returns a `ProcessedLocus` with nucleotide calls, optional binary calls,
independent called counts and acceptance booleans. It never accumulates a matrix.
`filter::Stats::observe` accumulates diagnostic counts; absent-GT diagnostics
remain separate. `genotype` exposes call types/policies and single-call functions;
`iupac::BaseMask` distinguishes known four-base ambiguity from missing data.
`cli::Cli` validates mode, sample names and seed pairing.

## Scientific contract

- An eligible SNP has one canonical A/C/G/T REF and at least one distinct
  canonical single-base ALT, case-insensitively. All REF/ALT bases must be
  distinct. ALT `.`, N, `*`, symbolic alleles, indels, MNPs and mixed records are
  rejected biologically. Malformed GT, including out-of-range indices, is still
  a strict error, even at a rejected site.
- One locus per record; no splitting or sorting. Nucleotide consensus supports
  arbitrary positive ploidy, ignoring phase and dosage. A=1/C=2/G=4/T=8 union
  masks map to all 15 standard uppercase IUPAC codes. Missing calls emit N;
  fully known four-base ambiguity also emits N but **counts as called**.
- Default `--partial-missing missing` makes any partially missing genotype N
  and excludes it from the nucleotide-called count. `called-alleles` uses known
  copies and counts a genotype with at least one known copy. `0/.` and `./0`
  always have identical called status. Fully missing/absent GT is never REF.
  No lowercase uncertainty encoding is used.
- Binary eligibility requires exactly REF + one distinct canonical ALT.
  Only fully called diploids encode: 0/0=0, 0/1 or 1/0=1, 1/1=2. Partial or
  missing calls and unsupported ploidies emit `?`; counters separate missing
  diploid/unknown-ploidy calls from unsupported known ploidies (ploidy takes
  precedence, including all-missing non-diploids). Only encodable calls count
  toward the binary threshold. Nucleotide policies never alter binary calls.
- Both eligibility/threshold paths always run, independent of future output
  selections. Multiallelic SNPs can pass DNA but never binary. Inspection prints
  accepted locus counts, not alignment files. Missing-call counters include
  eligible sites even if their threshold fails, never biologically rejected sites.
- `-r`/`--resolve-IUPAC` requires an explicit `--seed <u64>`; a seed without `-r`
  is rejected. Partial policy is applied first, then one known allele **copy**
  is sampled (dosage-weighted, not uniform over unique bases). Homozygotes are
  unchanged. Acceptance counts do not change, but inspection executes resolution
  and reports its policy/seed; resolved symbols are available through the library.

### Deterministic resolution v1

Each call gets a fresh ChaCha20 RNG: 32-byte key = seed (8 little-endian bytes),
original zero-based input record ordinal (8 little-endian bytes), then 16 zeros.
The ChaCha stream ID is the original zero-based sample index, checked to fit u64;
word position starts at zero. For n known copies in original GT order, read
`next_u64`, reject words below `2^64 mod n`, then use `word mod n` as the copy
index. This explicitly specified rejection sampling avoids modulo bias.
No shared RNG or standard-library hash is used. Calls are independent of
scheduling, output selection and the number of earlier RNG draws. Future
parallel callers must retain **original** record/sample keys, including rejected
records. Reordering/filtering samples or renumbering records can change results;
no invariance under those operations is promised. Golden vectors lock this v1
scheme, not upstream Python RNG equivalence.

## Milestones

- M1 complete: streaming readers and strict structural/genotype validation.
- M2 complete: genotype/IUPAC policies, independent filters, keyed seeded
  resolution and scientific inspection counters.
- M3/M4 future: sample ordering, matrix storage/transposition and alignment/site
  writers. No matrix, PHYLIP/FASTA/NEXUS output or used-site files exist yet.

## Development checks

```sh
cargo check
cargo fmt --all -- --check
cargo test --locked
cargo clippy --locked --all-targets -- -D warnings
```

Tests synthesize their own VCF, gzip, BGZF and BCF data. There are no upstream
sequence fixtures or copied conversion algorithms. See
[compatibility and provenance](docs/compatibility.md) and [GPL-3.0](LICENSE).
