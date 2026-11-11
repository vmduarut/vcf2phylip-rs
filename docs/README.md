# Documentation guide

## Using vcf2phylip-rs

- [Installation and worked examples](install.md): prerequisites, inspection,
  conversion and comparison against the included expected outputs.
- [Compatibility and provenance](compatibility.md): supported options,
  differences from upstream and deferred formats.
- [Scientific contract](../README.md#scientific-contract) and
  [output safety and resource limits](../README.md#conversion-and-output-safety):
  interpretation, filtering, storage and opt-in output replacement behavior.
- [Benchmarking](benchmarking.md): opt-in local measurements and their limits;
  it does not establish a published speedup.

## Maintainers and distributors

- [Packaging preparation](packaging.md): container and BioConda templates,
  validation status and outstanding work.
- [Release checklist](releases.md): toolchain policy, CI, candidate artifacts and
  steps requiring separate authorization before publication.

Version 0.1.0 is source-only for now; these templates and candidate workflows
are not published packages or releases. See [installation](install.md).
