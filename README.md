# vcf2phylip-rs

Convert VCF/BCF genotypes into nucleotide alignments for phylogenetic analysis,
or inspect input and filtering statistics without writing an alignment. This
independently designed Rust CLI and library streams input records and offers
in-memory or bounded-cell-buffer disk storage for conversion.

**Available:** FASTA and relaxed PHYLIP. NEXUS/SNAPP output is deferred for now.
Exact upstream output compatibility is not claimed.

Version 0.1.0 has not been published; no registry package, container image or
downloadable release is claimed to exist. Build from local source for now.
See the [documentation guide](docs/README.md), [installation and examples](docs/install.md),
[compatibility and provenance](docs/compatibility.md),
[benchmarking](docs/benchmarking.md), [packaging templates](docs/packaging.md),
and the [release checklist](docs/releases.md).

## Source install and quickstart

From the repository root, with Rust/Cargo 1.99.0 or later and a native linker:

```sh
cargo install --path . --locked
out=$(mktemp -d)
vcf2phylip-rs -i examples/example.vcf --inspect
vcf2phylip-rs -i examples/example.vcf -m 4 -f -w --output-folder "$out"
```

Ensure Cargo's bin directory is on your PATH. The conversion writes
`example.min4.phy`, `example.min4.fasta`, and `example.min4.used_sites.tsv` in
the new directory; they can be compared with the corresponding `examples/expected.*`
files. `-f` adds FASTA to the default PHYLIP output, and `-w` records retained
sites. Existing outputs are refused unless `--force` is supplied. `--locked` requires Cargo to use
the checked-in lockfile without changing dependency resolution.

For more options (replace input and scratch paths with your own):

```sh
cargo build --locked
cargo run --locked -- --input variants.vcf.gz --inspect
cargo run --locked -- -i variants.bcf --inspect -m 4 -o 'exact sample name'
cargo run --locked -- -i variants.vcf --inspect --partial-missing called-alleles -r --seed 42
cargo run --locked -- -i variants.vcf.gz --output-folder alignments
cargo run --locked -- -i variants.vcf.gz --output-folder alignments --force
cargo run --locked -- -i variants.bcf -f -w -o outgroup --output-prefix study
cargo run --locked -- -i variants.vcf -p -f -r --seed 42
cargo run --locked -- -i variants.vcf --storage disk --memory-limit 64MiB --temp-dir /existing/scratch -f -w
cargo run --locked -- --help
```

Tested with Rust/Cargo **1.99.0** on Linux; `rust-version = "1.99"` is the
conservative tested minimum, not a claim that earlier toolchains cannot work.
Direct dependencies are exactly pinned and `Cargo.lock` is included. All
requested versions resolved: noodles 0.119.0, clap 4.6.7, thiserror 2.0.21,
anyhow 1.0.104, flate2 1.1.11, tempfile 3.27.0, and same-file 1.0.6. No external compression
library is needed. Seeded resolution uses exactly
pinned rand_chacha 0.10.0 (rand_core 0.10.1 in the lockfile), without the rand
convenience API. There is no application parallelism; noodles brings rayon
transitively, and tempfile brings randomness dependencies for staging filenames.

## Shipped behavior

- Content-based detection (extensions are ignored): plain VCF, gzip VCF,
  concatenated/multimember gzip VCF, BGZF VCF, raw or gzip/BGZF BCF 2.2.
- Reads one record at a time with noodles adapters. Inspection memory scales
  with header and current record; conversion uses the selected matrix backend
  (see resource limits below). No input index required.
- Requires a valid header, at least one sample, and `FORMAT/GT` declared as
  `Number=1,Type=String`. Empty/duplicate sample names are errors. Names remain
  unchanged, including `./`, punctuation and spaces in inspection. Conversion
  rejects whitespace/control-containing names rather than silently changing them.
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
- `-m/--min-samples-locus` is an optional absolute called-sample minimum,
  capped to the input sample count with a warning; zero/negative values fail.
  `--min-call-rate` defaults to **0.9 (90%)**, accepting finite values in [0,1].
  **Choose one policy: explicit `-m` and explicit `--min-call-rate` cannot be
  combined**, even if the supplied rate is 0.9. `-m 1` alone requires only one
  called sample and `-m 2` requires two; the implicit default rate does not apply.
  Without `-m`, the effective minimum is
  `max(1, ceil(rate × input sample count))`, using all header samples (including
  outgroups). The effective count applies independently to nucleotide-called and
  binary-encodable sample counts. `-m 4` reproduces the previous count-only default.
  Neither option supplied means rate-based filtering at 0.9.
- `-o/--outgroup` requires an exact name; repeat it for multiple names, in
  requested order. Unknown or duplicate requests are errors. No comma/semicolon
  splitting, renaming, or actual sample reordering occurs during inspection.
- Inspection prints diagnostic counts, never creates output directories,
  matrices or used-site files. Output flags fail with `--inspect`, so they cannot silently be
  ignored. `--help`, `--version` and upstream-style `-v` work without input.

### Logging and reports

Final summaries go to **stderr**, only after success: input path, actual header sample count and
the active threshold policy (absolute minimum or call rate/percentage) and effective required sample
count, then successful filter totals, persisted output paths and
elapsed wall time. Conversion summaries previously went to stdout; this is an
intentional destination change. Conversion stdout is now empty; inspection's
existing diagnostic stdout labels remain, with the active configured minimum/rate
lines and an effective-threshold annotation.

`--quiet` suppresses informational messages, never warnings or errors. `--verbose`
adds configuration and missing-call diagnostics; it conflicts with `--quiet`.
`-v` still means version, not verbose. For example:

```sh
vcf2phylip-rs -i examples/example.vcf --inspect --quiet > inspection.txt
vcf2phylip-rs -i examples/example.vcf --inspect --verbose 2> run.log
```

Record totals partition into rejected non-SNP, below minimum called samples,
and retained nucleotide loci. The broad non-SNP category includes absent ALT,
noncanonical alleles (including indels/symbolic alleles), and duplicate alleles.
Below-minimum counts only SNP-eligible records; missing **sample calls** are not
rejected loci. No progress percentages, extra input scans, or per-sample stages
are inferred. Completion is reported only after successful processing and, for
conversion, persistence of all selected outputs. On partial publication failure,
the command fails without a success summary (see output safety below).

## Conversion and output safety

Without `--inspect`, relaxed sequential PHYLIP is the default. `-f` adds FASTA;
`-p -f` selects FASTA only. `-p` alone is an error. The reserved `-n` and `-b`
flags fail before output effects. Both writers use the same nucleotide matrix.
PHYLIP starts with `S L`, then full names separated from unwrapped sequences by
two spaces. FASTA wraps sequence lines at 80 bases. Identifiers retain their
exact UTF-8 bytes: no truncation, cleaning, escaping or remapping. Empty/duplicate
identifiers and any Unicode whitespace/control characters are rejected for both
formats, preventing ambiguous FASTA first-token identifiers and PHYLIP labels.
Non-whitespace punctuation and long names are preserved. Some downstream tools
may impose stricter ASCII or identifier conventions than these relaxed formats.

Repeated `-o` names come first in requested order, followed by all other samples
in original header order; this is only an ordering convention, not tree rooting.
Random calls are computed with original record/sample indices before ordering.
`-w` writes nucleotide-retained sites only, in input order including repeated
coordinates: tab-separated `#CHROM POS NUM_SAMPLES`, with original 1-based POS
and the called-sample count under the selected partial policy.

Outputs default to the current directory, or `--output-folder` (created if
missing). Names are `<prefix>.min<N>.phy`, `.fasta`, `.used_sites.tsv`, with
effective minimum N (after capping an explicit `-m`, or deriving it from the rate).
The default prefix is the input basename with a
terminal `.vcf`/`.bcf` and optional `.gz`/`.bgz`/`.bgzf` removed case-insensitively;
arbitrary middle `.vcf` substrings remain. Explicit prefixes must be nonempty
basenames, not `.`/`..`, absolute paths or slash/backslash-containing paths.

All selected targets are preflighted before records are processed. By default,
existing targets (including dangling symlinks) are refused. `--force` permits
replacement of existing regular files; checks refuse observed directories,
symlinks, other non-regular targets, and aliases of the input (including hardlinks). Only outputs
selected for this run are replaced; other files are untouched. `--inspect`
rejects `--force`, like other output flags. Targets are checked again immediately
before each publication, but these checks are not atomic with replacement: a
racing entry, including a symlink, can be replaced. Rename replaces the target
directory entry without following a target symlink. Observed existing regular
file permissions are applied to the staged file before replacement; failure to
apply them aborts that publication. Ownership, ACLs and other metadata are not
copied; newly created outputs retain the temporary file's restrictive permissions.
Validated flags, sample labels and outgroups precede output directory creation.
Files are staged with `NamedTempFile` in the target directory; all input must
validate through EOF and every writer must flush successfully before publication.
Header-only or zero-retained-locus conversion is an actionable error with no
final alignment or site files. Ordinary errors clean temporary staging files;
an empty newly created output directory may remain.

Final publication uses `persist_noclobber` by default, or `persist` with `--force`:
atomic **per file**, not a transaction across all selected files. Existing outputs
remain unchanged on input, validation or staged-writing errors. A filesystem
failure or racing creator during publication can leave earlier outputs finalized
(and their old contents replaced with `--force`); later outputs remain unchanged.
Diagnostics state this, and no rollback deletes or restores published files.
No durability across power loss,
cleanup after SIGKILL, or protection from a hostile directory owner is promised.

### Storage and resource limits

The memory backend stores one growable `Vec<u8>` per sample, no second transposed matrix
or binary matrix. Payload is S × L bytes, with checked dimensions; allocator
capacity, sample names, header/current record, metadata and I/O buffers add RAM.
For example 1,000 samples × 1,000,000 retained loci need at least 1 GB just for
sequence payload. Success reports payload bytes, **not peak RAM**. Memory is
unbounded in the number of loci; allocation failures are reported where possible,
but process/OS OOM termination is possible. Use datasets that fit available RAM.

`--storage memory` is the compatibility default and is **not bounded**.
`--storage disk` uses a single RAII-managed, restrictive-permission temporary
file; `--temp-dir` selects an **existing** directory, otherwise the OS temporary
directory is used. Output/site staging still lives in the output directory,
not `--temp-dir`. Ordinary errors remove both spool and staging files.
`--memory-limit` accepts positive integer bytes (optional `B`) or binary
`KiB`, `MiB`, `GiB` suffixes (case sensitive, no spaces). Default disk budget:
64 MiB. Explicit budget/temp options require disk mode; inspection rejects them.
There is no automatic backend switch or silent fallback.

The disk budget covers the **actual capacity of the one cell/tile vector**,
including the reused output-read buffer, not entire RSS. With S samples its
minimum is S bytes. Tile width T = min(floor(budget/S), 65,536) loci. A tile
allocates exactly S×T requested bytes (allocator-reported capacity is checked
against budget); no other cell vectors or full-sample reads are allocated by
storage. Dimensions/offsets are checked. Header/sample names, parsing and genotype
buffers for the current record, O(S) output-order validation, fixed writer/site
buffers and allocator/runtime overhead are **outside this budget**. Large headers,
single records/ploidies and OS page cache can still consume substantial memory.
This is a bounded matrix path, **not a strict RSS/OOM guarantee**.

Tiles flush sample-contiguous segments in original sample order. Uniform tile
offsets and the retained-site count determine every read; there is no growing
tile index or retained-record list. The final partial tile is packed on EOF.
Disk spool logical size is S×L bytes; allocated blocks may differ. Each writer
reads only each sample's own tile segments (O(S×L) bytes per format, with
O(S×ceil(L/T)) seeks), not the entire spool once per sample. Tiny budgets can
cause many seeks/writes. FASTA wrapping persists across chunks; PHYLIP streams
unwrapped rows. `-w` streams directly into its staged sidecar while processing,
so it adds disk usage but no per-locus memory list. Multiple output files add
staging disk space and one read pass each. Success payload bytes are neither
temporary-disk high-water nor peak RSS. See [benchmarking](docs/benchmarking.md).

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
`processing::VariantProcessor::new` is a low-level count-only constructor validating
sample dimensions, positive minimum (legacy `None` default 4, capped), partial
policy and allele resolution. CLI inspection/conversion instead use
`Cli::validate_samples` to compute the active threshold once and expose the
processor as `ValidatedConfig::processor`. Public `Cli::min_samples_locus` is now
`Option<usize>`: use `Some(n)` for count-only filtering, with `min_call_rate: None`.
`Cli::min_call_rate` is now `Option<f64>`: use `Some(rate)` for an explicit rate,
with `min_samples_locus: None`. Both `None` select the default 0.9 rate; both
`Some` fail public validation, just as supplying both flags fails CLI parsing.
`VariantProcessor::process`
returns a `ProcessedLocus` with nucleotide calls, optional binary calls,
independent `nucleotide_called_samples`/`binary_called_samples` counts and
acceptance booleans. It never accumulates a matrix.
`processing::Stats::observe` accumulates diagnostic counts, including
`records_processed`, `nucleotide_sites_accepted` and `binary_sites_accepted`;
absent-GT diagnostics remain separate. `genotype` exposes call types/policies,
`NucleotideResolution` and single-call functions;
`iupac::BaseMask` distinguishes known four-base ambiguity from missing data.
`cli::Cli` validates mode, sample names and seed pairing.
`conversion::convert(&Cli)` runs the staged conversion (no-clobber unless `Cli::force`) and returns
statistics/output paths.
`convert` never prints and has no observer/callback API. Its `ConversionReport`
includes `sample_count`, `matrix_payload_bytes`, input format/compression,
effective minimum and warnings captured during conversion, without reopening the
input. Callers own warning
reporting. The CLI prints each returned warning once, including in quiet mode,
then its final summary only after successful conversion. Failed conversions
return only the error, so configuration warnings are not reported on failure.
`writers::{fasta,phylip}::write` accept generic `Write`, a matrix and a sample
index permutation, and propagate buffered write/flush failures.
Both writers now accept any `matrix::Alignment` (existing `&Matrix` calls remain
valid). Its `samples`, `sites`, and fallible `row_chunks` callback expose bounded
reads without requiring a full row. `storage::DiskMatrix::{new,append,finish}`
implements the same trait; finish must precede reading. `Storage` selects either
backend. Library and CLI conversion share the same private processing loop for
both backends.

## Scientific contract

The default requires at least 90% of input samples called at each retained locus:
three samples require three calls; `--min-call-rate 0.66` requires two. Alternatively,
`-m 1` selects count-only filtering and retains an otherwise eligible locus with
at least one called sample. Explicit count and rate options are mutually exclusive.
Without `-m`, rate 0
still requires one call and rate 1 requires every sample. Rates parse as finite
`f64`; the
threshold calculation uses integer arithmetic on its canonical shortest decimal
representation, avoiding binary multiplication/ceil artifacts such as
`0.07 × 100` spuriously requiring eight rather than seven calls. Nucleotide and
diagnostic binary paths use the same effective count, but their own called
criteria below; neither is based merely on non-`N` characters in output.

By default, a genotype with more than one distinct known base is represented
by its IUPAC ambiguity code: for example, A/G becomes `R` and C/T becomes `Y`.
This also applies to higher-ploidy calls: the code represents the set of
distinct bases, not their dosage or phase. Thus A/A/G and A/G/G both become
`R`. A fully known A/C/G/T call and a missing call both appear as `N` in an
alignment, but the former counts as called and the latter does not. Partial
calls depend on the policy below; `--resolve-IUPAC` changes this representation.

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
- Both eligibility/threshold paths always run, independent of output
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

## Development checks

```sh
cargo check --locked
cargo fmt --all -- --check
cargo test --locked
cargo clippy --locked --all-targets -- -D warnings
python3 -m unittest discover -s scripts -p 'test_*.py'
cargo package --locked # add --allow-dirty only for local packaging review
```

Tests synthesize their own VCF, gzip, BGZF and BCF data. There are no upstream
sequence fixtures or copied conversion algorithms. See
[compatibility and provenance](docs/compatibility.md) and [GPL-3.0](LICENSE).
