#!/usr/bin/env bash
set -euo pipefail

# Source preparation must run `cargo vendor --locked vendor > vendor-config.toml`
# first. Do not download crates from a conda build sandbox.
test -d vendor
test -s vendor-config.toml
mkdir -p .cargo
cp vendor-config.toml .cargo/config.toml
cargo build --release --locked --offline
target_dir="${CARGO_TARGET_DIR:-target}"
if [[ -n "${CARGO_BUILD_TARGET:-}" ]]; then
    target_dir="${target_dir}/${CARGO_BUILD_TARGET}"
fi
mkdir -p "${PREFIX}/bin"
install -m755 "${target_dir}/release/vcf2phylip-rs" "${PREFIX}/bin/vcf2phylip-rs"
