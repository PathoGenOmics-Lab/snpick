# syntax=docker/dockerfile:1

# ---- Build stage: static musl binary (snpick has no C dependencies) ----
FROM rust:1-slim AS build
RUN apt-get update \
    && apt-get install -y --no-install-recommends musl-tools \
    && rm -rf /var/lib/apt/lists/*
RUN rustup target add x86_64-unknown-linux-musl
WORKDIR /src
COPY . .
RUN cargo build --release --target x86_64-unknown-linux-musl --bin snpick \
    && strip target/x86_64-unknown-linux-musl/release/snpick

# ---- Runtime stage: minimal scratch image ----
FROM scratch
LABEL org.opencontainers.image.title="snpick" \
      org.opencontainers.image.description="Fast extraction of variable sites from FASTA alignments" \
      org.opencontainers.image.source="https://github.com/PathoGenOmics-Lab/snpick" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"
COPY --from=build /src/target/x86_64-unknown-linux-musl/release/snpick /usr/local/bin/snpick
ENTRYPOINT ["/usr/local/bin/snpick"]
