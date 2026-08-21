---
tags:
  - install
---

# Installation

=== ":simple-anaconda: Bioconda"

    [![Bioconda version](https://anaconda.org/bioconda/snpick/badges/version.svg)](https://anaconda.org/bioconda/snpick)

    The recommended way to install SNPick:

    ```bash
    conda install -c bioconda snpick
    # or
    mamba install -c bioconda snpick
    ```

=== ":material-download: Pre-built binary"

    Grab the binary for your platform from the
    [latest release](https://github.com/PathoGenOmics-Lab/snpick/releases/latest): Linux
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

=== ":material-language-rust: From source"

    Requires a [Rust toolchain](https://rustup.rs/) (edition 2021, Rust 1.85 or newer).

    ```bash
    git clone https://github.com/PathoGenOmics-Lab/snpick.git
    cd snpick
    cargo build --release
    # Binary at target/release/snpick
    ```

    !!! note "Release builds are slow on purpose"
        The release profile enables fat LTO and a single codegen unit for maximum throughput, so
        a release build takes noticeably longer than a debug build.

=== ":material-docker: Container"

    ```bash
    docker run --rm -v "$PWD:/data" ghcr.io/pathogenomics-lab/snpick \
      -f /data/alignment.fasta -o /data/snps.fasta
    ```

    An Apptainer/Singularity image is available too; see the
    [`snpick.def`](https://github.com/PathoGenOmics-Lab/snpick/blob/main/snpick.def) recipe.

## Library & bindings

snpick is also a Rust **library crate** (these docs' API surface), with **Python**
([maturin](https://www.maturin.rs)/[pyo3](https://pyo3.rs), under `bindings/python`) and
**WebAssembly** (under `bindings/wasm`) bindings for notebooks, pipelines and the browser.

## Check the install

```bash
snpick --version
snpick --help
```
