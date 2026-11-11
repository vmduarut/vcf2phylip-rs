# Installation and worked examples

## Local source (available now)

Requirements: Rust/Cargo 1.99.0 or later, a native linker/C toolchain, and network
access to the locked crates on the first build (or an existing Cargo cache).
No external compression library is required. The current tested platform is
Linux; Linux/macOS CI is supplied but has not been executed remotely here.

From the repository root:

```sh
cargo build --release --locked
./target/release/vcf2phylip-rs --version
./target/release/vcf2phylip-rs --help
# Optional user-local installation; this writes into Cargo's installation root:
cargo install --path . --locked
```

The executable is `vcf2phylip-rs` (hyphens); the library crate is
`vcf2phylip_rs` (underscores). Put Cargo's bin directory on PATH if installing.
Version 0.1.0 has not been published. Do not assume `cargo install vcf2phylip-rs`,
a BioConda package, a container image, or a downloadable binary exists. Only
after an explicitly authorized publication should instructions name its real
URL/version.

## Reproducible tiny alignment

`examples/example.vcf` is independently generated, with four samples and three
canonical SNP records (including one multiallelic SNP). The expected files are
hand-computed, not copied from upstream. Default IUPAC calls retain all three
sites at an explicit four-sample, count-only threshold (also the effective
default for this four-sample input). For example, alpha's genotypes
0/0, 0/1 and 1/2 produce `A`, `Y` and `W` respectively.

Run from the repository root; `mktemp` avoids output-name collisions:

```sh
out=$(mktemp -d)
cargo run --locked -- -i examples/example.vcf --inspect -m 4 --min-call-rate 0
cargo run --locked -- -i examples/example.vcf -m 4 --min-call-rate 0 -f -w --output-folder "$out/memory"
cargo run --locked -- -i examples/example.vcf -m 4 --min-call-rate 0 -f -w --storage disk \
  --memory-limit 4B --temp-dir "$out" --output-folder "$out/disk"
for backend in memory disk; do
  diff examples/expected.phy "$out/$backend/example.min4.phy"
  diff examples/expected.fasta "$out/$backend/example.min4.fasta"
  diff examples/expected.used_sites.tsv "$out/$backend/example.min4.used_sites.tsv"
done
```

Inspection reports `Samples: 4; records: 3`, `Accepted nucleotide loci: 3`
and `accepted binary loci: 2`. The multiallelic site is eligible for the
nucleotide alignment but not the binary path; binary counts are diagnostics,
not an available binary output. Each `diff` is silent on a match. The FASTA
contains `alpha: AYW`, `beta: RTR`, `gamma: GCT` and `delta: ATG` (one record
per sample); the site report lists chr1 positions 1, 2 and 3, each with four
called samples. The two backends should produce the same files.

Lifecycle/status messages and conversion success summaries go to stderr, not
stdout (the conversion summary destination intentionally changed). Inspection's
diagnostic stdout retains existing labels with additive configured minimum/rate
lines and an effective-threshold annotation. Use `--quiet` to suppress informational
stderr while preserving warnings/errors, or `--verbose` for extra configuration
and missing-call diagnostics; these flags conflict. `-v` still prints version.
For example, capture inspection with `--inspect --quiet > inspection.txt`, or
save a conversion log with `2> run.log`. Successful conversion summaries list
record filtering totals, persisted output paths, matrix payload (not peak RSS),
and elapsed wall time; failed runs never print successful completion.

Four bytes is intentionally tiny to exercise tile boundaries, not a recommended
production budget. It bounds the cell buffer only, not RSS. See the README's
resource contract and [benchmarking](benchmarking.md) before sizing large jobs.
Outputs are refused by default; use a new directory/prefix, or explicitly rerun
with `--force` to replace only selected regular output files:

```sh
cargo run --locked -- -i examples/example.vcf -m 4 --min-call-rate 0 -f -w --output-folder "$out/memory" --force
```

Existing outputs survive input/validation/staged-writing errors. Publication is
atomic per file, not across all outputs: a later publication failure can leave
earlier files replaced, with no rollback. Checks refuse observed directories,
symlinks and aliases of the input (including hardlinks), even with `--force`.
Checks are not atomic with replacement: racing entries, including symlinks, can
be replaced, but rename does not follow a target symlink. Observed regular-file
permissions are copied to the staged replacement (failure aborts publication),
not ownership or ACLs. See the README's output safety contract for concurrency
and durability limitations.

Other supported modes (replace paths with your own; scratch must exist):

```sh
vcf2phylip-rs -i variants.bcf --inspect --partial-missing called-alleles
vcf2phylip-rs -i variants.vcf.gz -p -f -r --seed 42 --output-prefix resolved
vcf2phylip-rs -i variants.vcf.gz -f -o outgroup --storage disk \
  --memory-limit 64MiB --temp-dir /existing/scratch --output-folder alignments
```

Plain/gzip/BGZF VCF and BCF are detected by content. `-p -f` is FASTA-only;
without `-p`, relaxed PHYLIP remains enabled.
Outgroups must match exact sample names. Seeded calls use this project's defined
algorithm, not upstream Python RNG equivalence. Refer to the README and
[compatibility/provenance](compatibility.md) for the full scientific contract.

## Interpreting an alignment

`-m/--min-samples-locus` keeps its upstream name but now defaults to 1.
`--min-call-rate` defaults to 0.9 (90%), accepting finite values in [0,1].
Both conditions must hold: the effective minimum is the larger of `-m` (capped
to input sample count with a warning) and `ceil(rate × actual input samples)`.
Three samples therefore require three calls by default, or two with
`--min-call-rate 0.66`. **`-m 2` alone does not disable the default 90% condition.**
Use `-m 2 --min-call-rate 0` for upstream-style count-only filtering;
`-m 4 --min-call-rate 0` restores this project's previous default. Rate 1
requires all samples. Inspection, memory and disk conversion share the same
threshold, reflected in the output `.min<N>` suffix and reports.

IUPAC symbols encode the set of known bases, not phase or allele dosage: both
A/A/G and A/G/G become `R`. `N` can mean a missing call or a called four-base
ambiguity; only the ambiguity counts as called. By default, a partially
missing genotype becomes `N` and is not counted; `--partial-missing
called-alleles` instead uses its known copies. For example, A/C/G/T/. then
emits `N` as a four-base ambiguity and counts as called. These choices affect site
retention, so record the policy and threshold alongside downstream analyses.
With `-r --seed`, a known allele copy is sampled after applying that policy;
dosage can affect the result and the output no longer preserves ambiguity.

Only canonical SNP records meeting the called-sample threshold enter the
alignment. Indels, symbolic alleles and other unsupported biological alleles
are not converted; malformed genotypes still cause errors. This converter is
not a substitute for variant quality control or a complete VCF validator.
Outgroup ordering does not root a tree. Consult the downstream phylogenetic
tool's treatment of ambiguity and missing data before interpreting results.
