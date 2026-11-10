# Packaging preparation (not published packages)

## Validation status

Cargo package verification and the synthetic examples are exercised locally.
Docker is present, but **no container build/run was performed**. Apptainer and
Conda/BioConda tooling are unavailable here, so those recipes are **unvalidated**.
Neither images nor conda packages were published. Base-image tag/digest registry
availability was checked read-only; that does not validate build/runtime behavior.
Do not install tooling or publish anything merely to follow this document without
separate approval. Commands below are future validation instructions.

## Docker

The root Dockerfile builds with locked Cargo dependencies, using pinned Rust
1.99.0/bookworm and Debian bookworm-slim multiarchitecture image digests. Build
and runtime use the same Debian generation to avoid an obvious glibc mismatch.
The final image contains the executable and GPL license, not the compiler; input
and output default to `/data`. It runs as numeric nonroot UID/GID 65532.
`.dockerignore` allowlists only the build inputs, omitting repository/cache/output
data. Native architecture only; no multiarchitecture release was tested.

From the repository root, on an authorized Docker-capable host:

```sh
docker build -t vcf2phylip-rs:local .
docker run --rm vcf2phylip-rs:local --version
docker run --rm vcf2phylip-rs:local --help
out=$(mktemp -d)
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD/examples:/input:ro" -v "$out:/data" \
  vcf2phylip-rs:local -i /input/example.vcf -f -w
diff examples/expected.phy "$out/example.min4.phy"
```

Match the invoking UID/GID for bind-mount writes; ensure input is readable and
output/scratch directories writable by that identity. Do not solve permissions
with world-writable data or privileged/root containers. `/tmp` is the default
matrix scratch; for large disk jobs bind a writable scratch directory and pass
`--storage disk --memory-limit 64MiB --temp-dir /scratch`. Output staging still
uses the output directory. Host mount/SELinux/rootless Docker policies may require
site-specific adjustments. Networks are needed during build, not conversion.

`RUST_VERSION` defaults to 1.99.0. To change toolchains, **also** override
`RUST_IMAGE` with a verified compatible tag/digest (the pinned digest otherwise
still selects the old image). `RUNTIME_IMAGE` is separately configurable. Review
security updates to base images and dependency/license inventory before any OCI
publication. The recipes do not perform distribution/source-compliance work for you.

## Apptainer

`packaging/apptainer.def` is an independent two-stage recipe using the same pinned
base images, local source files and `cargo build --release --locked`. It does not
refer to a nonexistent published vcf2phylip image. Run from the repository root so
the `%files` paths resolve. Example on a host with configured fakeroot support:

```sh
apptainer build --fakeroot vcf2phylip-rs.sif packaging/apptainer.def
apptainer test vcf2phylip-rs.sif
apptainer run vcf2phylip-rs.sif --version
out=$(mktemp -d)
apptainer run --bind "$PWD/examples:/input:ro,$out:/output" \
  vcf2phylip-rs.sif -i /input/example.vcf -f -w --output-folder /output
diff examples/expected.phy "$out/example.min4.phy"
```

Building requires network access and a working Apptainer build environment;
fakeroot availability depends on site policy. Do not change system configuration
to enable it as part of M6. Runtime uses the invoking user's identity and host
bind/scratch policies. For HPC disk conversion, explicitly bind existing writable
scratch and select it with `--temp-dir`. Building the SIF, testing dependency
linkage, permissions and outputs on the intended cluster are outstanding.

## BioConda template and blockers

`packaging/bioconda/meta.yaml.in` is deliberately **not renderable as-is**. Every
`@...@` token is an unresolved submission prerequisite, not a real release URL or
checksum. `build.sh` must accompany the rendered `meta.yaml`. Do not submit until:

1. A reviewed public project URL, maintainer handle and immutable source version
   exist. Review BioConda/conda-forge policy and compiler availability on intended
   Linux/macOS targets; `compiler('rust')` must resolve to Rust >=1.99. The C
   compiler/stdlib selectors support the native linker and platform ABI. Determine
   runtime requirements from actual build/linkage results rather than guessing.
2. Prepare an authorized immutable **vendored** source tarball. From a clean
   matching source tree, run `cargo vendor --locked vendor > vendor-config.toml`
   (downloads into Cargo cache/vendor), preserve Cargo.lock, LICENSE, examples and
   generated vendor config, and include all dependency sources/licenses. The
   normal Cargo `.crate` archive is not vendored and is not this source tarball.
   Record its real SHA-256 and real published URL in the rendered recipe. No
   public tarball or hash is fabricated here.
3. `build.sh` copies that vendor config into `.cargo/config.toml` and builds
   `--release --locked --offline`. Verify the generated config uses relative
   `vendor` paths, inspect the source contents for secrets/local config, and test
   in a clean network-disabled sandbox. Do not rely on a developer's Cargo cache.
   The script installs only the executable into `$PREFIX/bin`; the conda recipe
   records LICENSE and exercises --help, --version and an exact example conversion.
4. Render, lint and build with current BioConda tooling and conda-forge/bioconda
   channels under their documented channel priority, then run package tests on
   every supported platform. Verify compiler version, relocation, runtime library
   dependencies, license metadata and source checksums. These checks have **not**
   been run. Only then request separate approval for a BioConda submission.

Vendoring, actual builds, registry uploads and submissions are not automated by
these assets. See the [release checklist](releases.md) for publication boundaries.
