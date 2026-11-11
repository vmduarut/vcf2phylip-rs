# Reproducible, opt-in benchmarking

No baseline is downloaded, tool installed, cache dropped, remote executed, or
build triggered by the stdlib-only Python 3 harness. Build explicitly first:

```sh
cargo build --release --locked
python3 -m unittest discover -s scripts -p 'test_*.py'
python3 scripts/benchmark.py --samples 8 --sites 200 --memory-limit 64KiB
```

Defaults are deliberately a **tiny smoke**, not a representative performance
study. Large datasets/full-scale runs require separate approval. The harness
creates unique retained sandboxes under ignored `target/bench/`; it never
overwrites outputs. `--sandbox-parent /tmp/opencode` is also supported. Remove
only your completed benchmark sandboxes manually when no longer needed.

For an independently obtained, locally supplied upstream **vcf2phylip 2.9**
script, expected provenance commit
`85a260e655b5f3c5dc2cb9915eb24552a002357d`:

```sh
python3 scripts/benchmark.py --python-tool /absolute/local/path/vcf2phylip.py \
  --samples 8 --sites 200 --seed 42 --repeats 2
```

The baseline runs with `--python-executable` (default: current Python 3). “2.9”
is the converter version, not a Python interpreter version. Supplying the
script authorizes executing that local file; it must be trusted. The report
records its SHA256 and expected provenance; **a hash alone does not verify its
commit identity**. Do not substitute an unrelated baseline and claim that commit.
Without this option the result is explicitly **Rust-only**, with no comparative
Python speedup claim. Neither mode prints automatic speedup conclusions.

## Correctness before comparability

The generator is deterministic for the recorded Python version/seed/settings.
It emits canonical A/C biallelic SNPs and fully called diploids or explicit
`./.`. Sample S0 is always called, so upstream `-m 1` retains exactly all
generated sites. Rust runs explicitly add `--min-call-rate 0` to match that
count-only filter rather than applying the new default 90% requirement.
This avoids known upstream partial-missing/multidigit/symbolic-allele bugs and
cross-language RNG differences. No randomized IUPAC resolution is requested.
All runs request PHYLIP, FASTA and retained-site reports. Validation compares
sample order/names, dimensions, per-sample SHA256 of normalized sequences and
every retained coordinate/called count with independently generated expectations.
Raw file hashes and captured stdout/stderr are retained too. FASTA wrapping and
PHYLIP separators may differ; only whitespace/case are normalized for sequences.
Unsupported baseline output formats fail validation rather than silently
claiming equivalence. Stats stdout is preserved; accepted dimensions and site
counts are checked, not inferred from exit success. Any failed command or
validation leaves `validated: false` and diagnostic files; no CSV is published.

## Measurement and reproducibility caveats

`report.json` and `runs.csv` contain argv, repeats/settings, wall time and
per-child peak RSS. Linux/macOS use `os.wait4`, **not cumulative RUSAGE_CHILDREN**:
Linux ru_maxrss KiB and macOS bytes normalize to bytes. Other platforms report
RSS as null with an unsupported marker. Wall time includes process startup,
not input generation, validation, hashing or build. Child CPU/library versions,
Rust/Cargo/tool/Python versions, binary/input/lock hashes, source commit/dirty
status, platform/uname/processor and CPU count are recorded. A binary hash does
not prove it was built from the recorded source; rebuild explicitly. Fixed
memory→disk→Python order and uncontrolled warm page caches are potential biases.
Use repeated runs and controlled hardware for a genuine performance study.
For very short commands, wait4's child high-water may include the launcher's
inherited pre-exec memory footprint; identical RSS floors are not evidence of
identical steady-state converter memory. The harness does not subtract this or
present tiny-run RSS as a precise matrix-memory measurement.

A polling thread samples **logical file bytes and allocated blocks × 512** in
two separate trees: matrix/temp spool and output/staging. Linux/macOS support
allocated-block accounting; the fallback on platforms without `st_blocks` is
zero (unavailable, not zero allocation). The default interval is 10 ms,
configurable with `--poll-seconds`; short-lived files and peaks can be missed.
Renames can also race directory scans. Independent logical/allocated and tree
peaks need not occur simultaneously and must not be summed as an exact total.
The output/staging high-water **includes final outputs**; final-output usage is
reported separately. Baseline temporary files obeying TMPDIR/TMP/TEMP land in
the spool tree; files it creates elsewhere cannot be measured by this sampler.
Polling overhead may perturb tiny runs. There is no disk quota or RSS cap.

The small local milestone validation ran Rust-only; no upstream baseline was
available/supplied, so no relative performance conclusion was made. Actual
numbers live in local ignored reports, not as hardcoded documentation claims.
