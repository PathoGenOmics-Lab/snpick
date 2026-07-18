# SNPick

<p align="center">
  <img src="assets/logo.png" alt="SNPick logo" width="620" />
</p>

**Fast, memory-efficient extraction of variable sites from FASTA alignments.**

SNPick extracts variable (SNP) sites from whole-genome FASTA alignments. It produces reduced
alignments ready for phylogenetic inference with ascertainment-bias correction (ASC) in
**IQ-TREE** and **RAxML**, and optionally generates VCF files.

!!! question "Why not snp-sites?"
    snp-sites works well for small datasets but struggles with large alignments — it loads the
    whole matrix into memory and scales poorly. SNPick uses a zero-copy, memory-mapped
    architecture that handles thousands of genomes in seconds with minimal RAM.

## SNPick vs snp-sites

| | **SNPick** | **snp-sites** |
|---|---|---|
| Architecture | Zero-copy mmap, parallel scan | Full matrix in memory |
| 250 seqs × 4.4 Mbp | **0.9 s**, 105 MB | 9.5 s, 520 MB |
| 1000 seqs × 4.4 Mbp | **~3 s**, ~140 MB | >26 min (killed), 3+ GB |
| ASC `fconst` output | :material-check: Built-in | :material-close: Not supported |
| VCF output | :material-check: Optional | :material-check: Default |
| Gap handling | :material-check: Optional (`-g`) | :material-check: Default |
| IUPAC ambiguous | :material-check: Tracked as ambiguous | :material-alert: Treated as variant |

## Quick start

```bash
# Install from Bioconda
conda install -c bioconda snpick

# Extract variable sites
snpick -f alignment.fasta -o snps.fasta

# With a VCF, on 8 threads, quietly
snpick -f alignment.fasta -o snps.fasta --vcf -t 8 -q
```

<div class="grid cards" markdown>

- :material-download: **[Installation](installation.md)** — Bioconda, source, or a pre-built binary
- :material-console: **[Usage](usage.md)** — every flag, with examples
- :material-file-document: **[Output formats](output.md)** — reduced FASTA, VCF v4.2, ASC `fconst`
- :material-chart-line: **[Benchmarks](benchmarks.md)** — scaling by sequences and length

</div>

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
