# snpick (Python bindings)

Python bindings for [snpick](https://github.com/PathoGenOmics-Lab/snpick), built on the Rust
library crate with [pyo3](https://pyo3.rs) and [maturin](https://www.maturin.rs).

## Build

```bash
cd bindings/python
pip install maturin
maturin develop --release   # or: maturin build --release
```

## Usage

```python
import snpick

s = snpick.stats("alignment.fasta")          # gzip and "-" (stdin) accepted
print(s["variable_sites"], s["fconst"])
```

`stats(path, include_gaps=False)` returns a dict with `sequences`, `alignment_length`,
`variable_sites`, `constant_sites`, `ambiguous_sites`, and the `fconst` tuple `(A, C, G, T)`.
