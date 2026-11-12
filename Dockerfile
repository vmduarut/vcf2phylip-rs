# Registry tags were checked during M6; digests pin the actual multiarch images.
ARG RUST_VERSION=1.99.0
ARG RUST_IMAGE=rust:${RUST_VERSION}-bookworm@sha256:114c7a4425406451c2866b6aafe69fe29b1b298832db1277d411ac73c82d04d6
ARG RUNTIME_IMAGE=debian:bookworm-slim@sha256:734a72b3d6142a71993b4bc580ca01d8ca687f77a2c422d96a961d3376b179bc
FROM ${RUST_IMAGE} AS build
WORKDIR /build
COPY Cargo.toml Cargo.lock LICENSE README.md ./
COPY src ./src
RUN cargo build --release --locked

FROM ${RUNTIME_IMAGE}
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /build/target/release/vcf2phylip-rs /usr/local/bin/vcf2phylip-rs
COPY LICENSE /usr/share/licenses/vcf2phylip-rs/LICENSE
RUN mkdir /data && chown 65532:65532 /data
USER 65532:65532
WORKDIR /data
CMD ["/usr/local/bin/vcf2phylip-rs", "--help"]
