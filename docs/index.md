---
tags:
  - overview
---

# SNPick

<p align="center" markdown>
  ![SNPick logo](assets/logo.png){ width="600" }
</p>

**Fast, memory-efficient extraction of variable sites from FASTA alignments.**

[Get started :material-rocket-launch:](installation.md){ .md-button .md-button--primary }
[Usage reference :material-console:](usage.md){ .md-button }
[View on GitHub :fontawesome-brands-github:](https://github.com/PathoGenOmics-Lab/snpick){ .md-button }

SNPick extracts variable (SNP) sites from whole-genome FASTA alignments. It produces reduced
alignments ready for phylogenetic inference with ascertainment-bias correction (ASC) in
**IQ-TREE** and **RAxML**, and optionally generates VCF files.

!!! question "Why not snp-sites?"
    snp-sites works well for small datasets but struggles with large alignments — it loads the
    whole matrix into memory and scales poorly. SNPick uses a zero-copy, memory-mapped
    architecture that handles thousands of genomes in seconds with minimal RAM.

## What you get

<div class="grid cards" markdown>

- :material-dna:{ .lg .middle } __Variable sites only__

    ---

    A reduced FASTA with just the informative columns — a drop-in, much smaller input for
    phylogenetics.

- :material-tree:{ .lg .middle } __ASC-ready__

    ---

    Constant-site counts (`fconst`) printed for IQ-TREE's `+ASC` models, so branch lengths stay
    unbiased.

- :material-file-table:{ .lg .middle } __Optional VCF__

    ---

    VCF v4.2 with per-sample genotypes, a configurable contig name, alignment-column `POS`
    (or ungapped reference coordinates with `--ref-coords`), and gap/ambiguity handling.

- :material-filter:{ .lg .middle } __Filter & mask__

    ---

    Per-site missingness / MAC / MAF / allele filters, BED region masking and sample
    selection — without breaking `fconst`.

- :material-pipe:{ .lg .middle } __Pipeline-native__

    ---

    gzip & stdin/stdout streaming, PHYLIP/NEXUS output, and a `--stats-json` sidecar so nothing
    scrapes stderr.

- :material-lightning-bolt:{ .lg .middle } __Built for scale__

    ---

    Zero-copy mmap + parallel, auto-vectorized scan: **O(L)** memory, thousands of genomes in
    seconds.

</div>

## SNPick vs snp-sites

| | **SNPick** | **snp-sites** |
|---|---|---|
| Architecture | Zero-copy mmap, parallel scan | Full matrix in memory |
| 250 seqs × 4.4 Mbp | **0.9 s**, 105 MB | 9.5 s, 520 MB |
| 1000 seqs × 4.4 Mbp | **~3 s**, ~140 MB | >26 min (killed), 3+ GB |
| ASC `fconst` output | :material-check:{ .snpick-yes } Built-in | :material-close: Not supported |
| VCF output | :material-check:{ .snpick-yes } Optional | :material-check: Default |
| Gap handling | :material-check:{ .snpick-yes } Optional (`-g`) | :material-check: Default |
| IUPAC ambiguous | :material-check:{ .snpick-yes } Tracked as ambiguous | :material-alert: Treated as variant |

## Quick start

```bash
# Install from Bioconda
conda install -c bioconda snpick

# Extract variable sites
snpick -f alignment.fasta -o snps.fasta

# With a VCF, on 8 threads, quietly
snpick -f alignment.fasta -o snps.fasta --vcf -t 8 -q
```

## Citation

If you use SNPick in your research, please cite:

```bibtex
@software{snpick,
  author  = {Ruiz-Rodriguez, Paula and Coscolla, Mireia},
  title   = {SNPick: Fast extraction of variable sites from FASTA alignments},
  url     = {https://github.com/PathoGenOmics-Lab/snpick},
  doi     = {10.5281/zenodo.14191809},
  license = {GPL-3.0}
}
```

<sub>Paula Ruiz-Rodriguez and Mireia Coscolla — Institute for Integrative Systems Biology,
I²SysBio, University of Valencia-CSIC, Valencia, Spain.</sub>
