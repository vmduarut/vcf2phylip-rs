# Toolchain policy and release checklist

## Version and MSRV policy

Current package version: **0.1.0**, pre-1.0. Milestone numbers are development
checkpoints, not published version numbers. No release has been published by M6.
Use semantic versioning: before 1.0, incompatible public API/CLI/scientific
contract changes require a minor bump; compatible fixes use patch bumps. After
1.0, incompatible changes require a major bump. Describe any alignment-byte,
seeded-resolution, filtering, error, or storage-contract change explicitly in
release notes; it is not merely a performance detail.

`Cargo.toml` declares **Rust 1.99** as the minimum supported Rust version (MSRV),
tested locally with rustc/Cargo 1.99.0. This conservative policy does not claim
earlier versions fail. Do not lower it based only on dependencies' rust-version
metadata: first test the entire locked graph and all targets with that toolchain.
CI explicitly tests 1.99.0 and moving stable on Linux/macOS. An MSRV increase
requires release notes, explicit review, and synchronized Cargo, CI, recipes and
documentation changes. Direct dependency pins and Cargo.lock are intentional;
update them together and rerun the full checks, including packaging verification.

The official 1.99.0 channel manifest and the container tags/digests in the recipes
were reachable during M6 preparation. That is availability evidence, not proof of
remote CI success or a container build. No nightly fallback is used. Action SHAs
were resolved from the upstream checkout v4, upload-artifact v4 and
dtolnay/rust-toolchain stable refs via GitHub's API; review and refresh pins
deliberately instead of allowing an unreviewed moving action dependency.

## CI and manual candidate builds

`.github/workflows/ci.yml` runs fmt, locked check/clippy/tests, the Python stdlib
tests and Cargo's package verification build on Linux/macOS for MSRV and stable.
No publishing permissions or secrets are required. Local checks do not substitute
for the first successful hosted run, especially macOS.

`.github/workflows/release.yml` runs **only** by `workflow_dispatch`, with a
default-false acknowledgement gate. Select a reviewed ref and set acknowledgement
true to build candidates from the dispatch's exact commit. It never responds to
pushes/tags, creates releases, publishes crates/images, or writes repository
contents. Artifacts expire after 14 days. There is no publication job to approve.
Repository maintainers must separately authorize any eventual public distribution.

Each native host produces a tar.gz and SHA-256 sidecar containing the executable,
LICENSE, README, docs, examples, verified Cargo source package (.crate), and
BUILD.txt with commit/toolchain/host. The host triple is discovered, not guessed
from the runner label. These are unsigned native-runner candidates: no cross-build,
bit-for-bit reproducibility, universal macOS support, old-glibc compatibility,
notarization, signing, SBOM, or public availability is promised. Determine and
document runtime OS/architecture requirements on real supported hosts.

## Checklist before any separately authorized publication

1. Review a clean commit; confirm version/lockfile, scientific changes, MSRV,
   changelog/release notes, attribution and license. Keep deferred features clear.
2. Run the README development checks. Run `cargo package --locked --list` and
   `cargo package --locked`; review the archive for LICENSE, README, core source,
   tests and examples, and absence of caches/benchmarks/build outputs. The package
   allowlist in Cargo.toml intentionally excludes `.git`, `.tgrep`, and `target`.
   `--allow-dirty` is for local review only, never final release provenance.
3. Require successful Linux/macOS MSRV/stable CI. Manually build candidates and
   compare their BUILD.txt commit with the reviewed commit. Verify the SHA-256
   after downloading (e.g. `sha256sum -c NAME.tar.gz.sha256` on Linux or
   `shasum -a 256 -c NAME.tar.gz.sha256` on macOS). Run --version/--help and the
   example conversion on every supported host. Review linkage and permissions.
4. Check source-package rebuilds using the locked graph and supported toolchain.
   `.crate` is a gzip tar archive; unpack it and run `cargo build --release --locked`.
   Cargo package verification itself builds the library/binary, not every test;
   run the unpacked package's tests too. Crates are fetched unless already cached.
5. Before binary/container distribution, review GPL-3.0 corresponding-source and
   third-party license obligations. Supply the matching source, lockfile, build
   instructions and required dependency sources/notices with an appropriate durable
   source offer/access mechanism; a short-lived Actions artifact alone is not a
   long-term source distribution plan. Preserve [upstream attribution](compatibility.md).
6. Container builds and BioConda rendering/build/lint are separate outstanding
   validations; follow [packaging](packaging.md). Do not label templates as tested
   releases. Confirm image architectures, runtime libraries and nonroot behavior.
7. Only after explicit authorization: assign real repository/release URLs (none
   invented in Cargo metadata), review the final source/binary checksums, arrange
   signing/provenance as required, and perform the chosen publication manually.
   Crates.io, GitHub releases, BioConda, Docker/OCI registries and remote pushes
   are all outside the supplied workflows. Update install docs only after the
   published artifacts actually exist; retain a rollback/advisory plan.
