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
Do not assume `cargo install vcf2phylip-rs`, a BioConda package, or a downloadable
binary exists: registry/release publication is not part of M6. Only after an
explicitly authorized publication should instructions name its real URL/version.

## Reproducible tiny alignment

`examples/example.vcf` is independently generated, with four samples and three
canonical SNP records (including one multiallelic SNP). The expected files are
hand-computed, not copied from upstream. Default IUPAC calls retain all three
sites at the default four-sample threshold. For example, alpha's genotypes
0/0, 0/1 and 1/2 produce `A`, `Y` and `W` respectively.

Run from the repository root; `mktemp` avoids output-name collisions:

```sh
out=$(mktemp -d)
cargo run --locked -- -i examples/example.vcf --inspect
cargo run --locked -- -i examples/example.vcf -f -w --output-folder "$out/memory"
cargo run --locked -- -i examples/example.vcf -f -w --storage disk \
  --memory-limit 4B --temp-dir "$out" --output-folder "$out/disk"
for backend in memory disk; do
  diff examples/expected.phy "$out/$backend/example.min4.phy"
  diff examples/expected.fasta "$out/$backend/example.min4.fasta"
  diff examples/expected.used_sites.tsv "$out/$backend/example.min4.used_sites.tsv"
done
```

Four bytes is intentionally tiny to exercise tile boundaries, not a recommended
production budget. It bounds the cell buffer only, not RSS. See the README's
resource contract and [benchmarking](benchmarking.md) before sizing large jobs.
Outputs are never overwritten; use a new directory or prefix for each run.

Other supported modes (replace paths with your own; scratch must exist):

```sh
vcf2phylip-rs -i variants.bcf --inspect --partial-missing called-alleles
vcf2phylip-rs -i variants.vcf.gz -p -f -r --seed 42 --output-prefix resolved
vcf2phylip-rs -i variants.vcf.gz -f -o outgroup --storage disk \
  --memory-limit 64MiB --temp-dir /existing/scratch --output-folder alignments
```

Plain/gzip/BGZF VCF and BCF are detected by content. `-p -f` is FASTA-only;
without `-p`, relaxed PHYLIP remains enabled. NEXUS and SNAPP are **not** shipped.
Outgroups must match exact sample names. Seeded calls use this project's defined
algorithm, not upstream Python RNG equivalence. Refer to the README and
[compatibility/provenance](compatibility.md) for the full scientific contract.
