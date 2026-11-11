# BCF decoder boundaries

The BCF 2.2 adapter preflights decompressed headers and one record at a time
before calling the pinned noodles 0.119.0 (bcf 0.94.0 / vcf 0.96.0) decoder.
That is also the latest published noodles release checked during this fix.
No dependency panic catcher is used as a validation mechanism.

Default limits, independent of `--memory-limit` and alignment storage mode:

| Resource | Limit |
| --- | --- |
| Decompressed BCF header | 16 MiB |
| Combined shared + individual record bytes | 64 MiB |
| Slots in each string/contig dictionary | 1,000,000 |
| Samples × max(FORMAT fields, 1) per record | 1,000,000 |
| Typed elements, including keys, per record | 4,000,000 |

Header and record declarations are checked before buffering their payloads.
Payloads are read incrementally using fallible buffer reservations, so a short
file declaring a large payload does not allocate that declared size upfront.
Parsed header IDX values are checked (including checked addition) **before**
noodles inserts into its sparse dictionary. Implicit indices are bounded too.
The dictionary slot accounting is conservative for repeated declarations.

Typed descriptors are checked without recursive length parsing; extended
lengths must be scalar integers. Payload byte counts use checked arithmetic and
must fit in the frame. Sample counts must match the header before decoding.
FORMAT types are checked against their definitions; null/zero-length FORMAT
values, reserved numeric sentinels, scalar end-of-vector sentinels, empty
Character elements, duplicate keys, and misplaced GT return errors rather than
reaching known dependency panic paths. INFO end-of-vector sentinels are rejected.

GT currently supports Int8 encoding only. Int16 and Int32 GT (which can be valid
BCF) return an explicit unsupported-encoding error; use VCF for these inputs.
Float/null GT encodings and invalid negative GT codes are rejected. Non-GT
Int16/Int32 integer fields remain supported. These checks should be revisited
when upgrading noodles, particularly its sample and INFO value decoders.
Like the pinned decoder, the adapter tolerates uninterpreted trailing individual
block bytes (the noodles writer can emit excess mixed-ploidy GT padding).

Ordinary gzip, concatenated gzip, and BGZF remain streaming: only the bounded
header/current record is buffered, not the whole decompressed file. Library
callers processing trusted, unusually large inputs can set `reader::BcfLimits`
with `VariantReader::from_reader_with_bcf_limits`. There is no CLI override.

These are **not a total process-memory guarantee**, a general-purpose sandbox,
or a bound on VCF text parsing, decompression CPU, header object overhead, or
alignment/output storage. Noodles still makes ordinary infallible allocations
within these limits, and decoded records and temporary frame copies consume
more memory than their encoded size. Run untrusted workloads under suitable
OS memory/CPU/disk limits as well. Increasing the limits weakens these bounds.
