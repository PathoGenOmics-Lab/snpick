# Installation

## Bioconda (recommended)

[![Bioconda version](https://anaconda.org/bioconda/snpick/badges/version.svg)](https://anaconda.org/bioconda/snpick)

```bash
conda install -c bioconda snpick
# or
mamba install -c bioconda snpick
```

## Pre-built binary

Grab the binary for your platform from the
[latest release](https://github.com/PathoGenOmics-Lab/snpick/releases/latest) — Linux
(`x86_64`, `aarch64`) and macOS (`x86_64`, `aarch64`), with `SHA256SUMS.txt` published for
verification:

```bash
# choose one:
#   snpick-linux-x86_64 | snpick-linux-aarch64 | snpick-macos-x86_64 | snpick-macos-aarch64
curl -LO https://github.com/PathoGenOmics-Lab/snpick/releases/latest/download/snpick-linux-x86_64
chmod +x snpick-linux-x86_64
./snpick-linux-x86_64 --help
```

!!! tip "Verify the download"
    Fetch `SHA256SUMS.txt` from the same release and compare the checksum of your binary
    against the matching line.

## From source

Requires a [Rust toolchain](https://rustup.rs/) (edition 2021).

```bash
git clone https://github.com/PathoGenOmics-Lab/snpick.git
cd snpick
cargo build --release
# Binary at target/release/snpick
```

The release profile enables fat LTO and a single codegen unit for maximum throughput, so a
release build takes noticeably longer than a debug build — this is expected.

## Check the install

```bash
snpick --version
snpick --help
```
