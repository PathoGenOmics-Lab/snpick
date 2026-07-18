---
tags:
  - performance
---

# Benchmarks

Benchmarks use simulated *M. tuberculosis*-like genomes (4.4 Mbp, ~65% GC, 3.6% variable sites).

## Scaling by number of sequences

<p align="center">
  <img src="../assets/benchmark.png" alt="Benchmark: sequence scaling" width="760" />
</p>

## Scaling by sequence length

<p align="center">
  <img src="../assets/benchmark_length.png" alt="Benchmark: length scaling" width="760" />
</p>

SNPick maintains **O(L)** memory regardless of sequence count, while snp-sites requires
**O(N × L)** — it holds the full matrix in memory and is eventually killed on large inputs.

| Dataset | SNPick | snp-sites |
|---|---|---|
| 250 seqs × 4.4 Mbp | **0.9 s**, 105 MB | 9.5 s, 520 MB |
| 1000 seqs × 4.4 Mbp | **~3 s**, ~140 MB | >26 min (killed), 3+ GB |

!!! tip "Reproducing"
    Wall-clock depends on core count; pin it with [`--threads`](usage.md#control-threads-hpc-reproducibility)
    for comparable runs. The extracted sites and VCF are identical regardless of thread count.
