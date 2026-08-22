<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="logo/SNPick-dark.png" />
    <img src="logo/SNPick.png" alt="SNPick logo" width="720" />
  </picture>
</p>

<div align="center">

**Fast, memory-efficient extraction of variable sites from FASTA alignments.**

[![Documentation](https://img.shields.io/badge/docs-online-007a6c?style=flat-square)](https://pathogenomics-lab.github.io/snpick/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/snpick?style=flat-square&label=bioconda&color=007a6c)](https://anaconda.org/bioconda/snpick)
[![Downloads](https://img.shields.io/conda/dn/bioconda/snpick?style=flat-square&label=downloads)](https://anaconda.org/bioconda/snpick)
[![License: GPL v3](https://img.shields.io/badge/license-GPL%20v3-007a6c?style=flat-square)](LICENSE)
[![DOI](https://img.shields.io/badge/doi-10.5281%2Fzenodo.14191809-007a6c?style=flat-square)](https://doi.org/10.5281/zenodo.14191809)

**[📖 Read the documentation →](https://pathogenomics-lab.github.io/snpick/)**

<sub>Paula Ruiz-Rodriguez and Mireia Coscolla, Institute for Integrative Systems Biology, I<sup>2</sup>SysBio, University of Valencia-CSIC, Valencia, Spain</sub>

</div>

---

SNPick pulls the variable (SNP) sites out of a whole-genome FASTA alignment and writes a reduced alignment ready for phylogenetics: ascertainment-bias correction (ASC `fconst`) for **IQ-TREE** and **RAxML**, plus an optional VCF. A zero-copy, memory-mapped scan handles thousands of genomes in seconds on a laptop, where `snp-sites` runs out of memory (250 sequences in 1.7 s; at 1000 it is killed).

## Quick start

```bash
# Install from Bioconda
conda install -c bioconda snpick

# Extract the variable sites
snpick -f alignment.fasta -o snps.fasta

# ...with a VCF, on 8 threads
snpick -f alignment.fasta.gz -o snps.fasta --vcf -t 8
```

That is the core of it. The flags for per-site filtering, region masking, VCF coordinates, alternative output formats (PHYLIP / NEXUS) and QC all live in the documentation, and `snpick --help` prints the authoritative list.

## Documentation

Everything is at **[pathogenomics-lab.github.io/snpick](https://pathogenomics-lab.github.io/snpick/)**:

- **[Installation](https://pathogenomics-lab.github.io/snpick/installation/)**: Bioconda, pre-built binaries, from source, Docker / Apptainer
- **[Usage](https://pathogenomics-lab.github.io/snpick/usage/)**: the complete flag reference and common recipes
- **[Tutorial](https://pathogenomics-lab.github.io/snpick/tutorials/tutorial/)**: a hands-on walkthrough
- **[Output formats](https://pathogenomics-lab.github.io/snpick/output/)**: FASTA, PHYLIP, NEXUS and the VCF layout
- **[Benchmarks](https://pathogenomics-lab.github.io/snpick/benchmarks/)**: how it scales against snp-sites
- **[Architecture](https://pathogenomics-lab.github.io/snpick/architecture/)**: the two-pass memory-mapped scan

snpick is also a Rust **library crate**, with **Python** (`bindings/python`) and **WebAssembly** (`bindings/wasm`) bindings.

## Citation

If SNPick helps your research, please cite it:

```bibtex
@software{snpick,
  author  = {Ruiz-Rodriguez, Paula and Coscolla, Mireia},
  title   = {SNPick: Fast extraction of variable sites from FASTA alignments},
  url     = {https://github.com/PathoGenOmics-Lab/snpick},
  doi     = {10.5281/zenodo.14191809},
  license = {GPL-3.0}
}
```

<h2 align="center">✨ Contributors</h2>

<div align="center">
<table>
  <tr>
    <td align="center">
      <a href="https://github.com/paururo">
        <img src="https://avatars.githubusercontent.com/u/50167687?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Paula Ruiz-Rodriguez</b></sub>
      </a>
      <br />
      <a href="" title="Code">💻</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Data">🔣</a>
      <a href="" title="Design">🎨</a>
      <a href="" title="Tool">🔧</a>
    </td>
    <td align="center">
      <a href="https://github.com/mireiacoscolla">
        <img src="https://avatars.githubusercontent.com/u/29301737?v=4&s=100" width="100px;" alt=""/>
        <br />
        <sub><b>Mireia Coscolla</b></sub>
      </a>
      <br />
      <a href="" title="Funding">🔍</a>
      <a href="" title="Ideas">🤔</a>
      <a href="" title="Mentoring">🧑‍🏫</a>
      <a href="" title="Research">🔬</a>
      <a href="" title="User Testing">📓</a>
    </td>
  </tr>
</table>

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification ([emoji key](https://allcontributors.org/docs/en/emoji-key)).

</div>

<div align="center">
<sub>Released under the <a href="LICENSE">GPL-3.0</a> license.</sub>
</div>
