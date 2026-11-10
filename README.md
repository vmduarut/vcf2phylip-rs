# vcf2phylip-rs — Milestone 1

An independently designed Rust library and thin CLI for **streaming VCF/BCF
inspection**. Alignment conversion is **not implemented**. This is not yet a
replacement for vcf2phylip, nor a claim of scientific-output compatibility.

## Build and inspect

```sh
cargo build --locked
cargo run --locked -- --input variants.vcf.gz --inspect
cargo run --locked -- -i variants.bcf --inspect -m 4 -o 'exact sample name'
cargo run --locked -- --help
```

Tested with Rust/Cargo **1.99.0** on Linux; `rust-version = "1.99"` is the
conservative tested minimum, not a claim that earlier toolchains cannot work.
Direct dependencies are exactly pinned and `Cargo.lock` is included. All
requested versions resolved: noodles 0.119.0, clap 4.6.7, thiserror 2.0.21,
anyhow 1.0.104, flate2 1.1.11, and dev dependency tempfile 3.27.0. No native
libraries or toolchain installation is needed. There is no application RNG or
parallelism; noodles brings rayon transitively, and tempfile brings its own
test-only randomness dependencies.

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
  configuration validation only: **no filtering is performed**.
- `-o/--outgroup` requires an exact name; repeat it for multiple names, in
  requested order. Unknown or duplicate requests are errors. No comma/semicolon
  splitting, renaming, or actual sample reordering occurs during inspection.
- Inspection prints diagnostic counts, never creates output directories,
  matrices or used-site files. Without `--inspect`, conversion fails clearly.
  Reserved output flags also fail with `--inspect`, so they cannot silently be
  ignored. `--help`, `--version` and upstream-style `-v` work without input.

This is structural/genotype validation, not a complete VCF specification
validator or scientific quality-control tool. It does not impose sorted loci,
canonical nucleotides or a SNP-only policy. Symbolic alleles and indels are
preserved for inspection, not converted. For untrusted input there are no
explicit header/record size limits. BGZF is consumed as concatenated gzip;
gzip checksums/truncated members are checked, but a missing optional BGZF EOF
marker alone is not an error and BGZF indexing metadata is not validated.

## Library

`reader::VariantReader::{open,from_reader,next_variant}` exposes the parsed
header, detected format, compression status, successful record count and a
`Variant` containing a noodles record plus optional genotype buffers in sample
order. An error terminates the reader; callers must handle that error, not
interpret a subsequent `None` as successful validation.
`cli::Cli` and its validation methods provide the current configuration contract.

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
