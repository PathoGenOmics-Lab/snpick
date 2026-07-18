# Changelog

All notable changes to SNPick are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.2] - 2026-07-18

### Added
- `-t, --threads <N>` to pin the Rayon thread pool (deterministic wall-clock on
  shared HPC/SLURM nodes). Thread count never changes the output, only speed.
- `--chrom <NAME>` to set the VCF `CHROM` / `##contig` name (default `1`, e.g.
  `NC_000962.3`).
- `-q, --quiet` to silence the `[snpick]` progress logs (errors still print).
- A valid header-only VCF is now written when `--vcf` is requested but the
  alignment has no variable sites, so pipelines that declare the `.vcf` as an
  output no longer break.

### Fixed
- Silent SNP loss when a single-line FASTA record contained a blank line: such
  records were read with misaligned byte offsets. They now use the
  newline-skipping scanner.
- Gap reference bases produced an invalid VCF `REF` of `-` under `--include-gaps`;
  they are now rendered as `*`, matching the `ALT` encoding.
- All-gap columns were tallied in no category under `--include-gaps`, breaking the
  reported `variable + constant + ambiguous == length` invariant.
- Bare output filenames (`-o snps.fasta`, no directory) were rejected during path
  resolution, breaking the documented quickstart.
- The release workflow's "tag already exists" guard was never read, so merges
  without a version bump re-published over existing tags.
- The README's pre-built binary download command pointed at an asset name the
  release never publishes (404).

### Performance
- VCF data rows are assembled in a reused byte buffer instead of one formatted
  write per genotype field — a large speedup on many-sample inputs. Output is
  byte-identical.

## [1.0.1]

- First Bioconda release.
- Cross-platform Build & Release CI workflow (Linux/macOS, x86_64/aarch64).
- README overhaul and benchmarks.

## [1.0.0]

- Initial release: zero-copy memory-mapped extraction of variable sites from
  FASTA alignments, optional VCF v4.2 output, and ASC `fconst` reporting.

[1.0.2]: https://github.com/PathoGenOmics-Lab/snpick/compare/1.0.1...1.0.2
[1.0.1]: https://github.com/PathoGenOmics-Lab/snpick/compare/1.0.0...1.0.1
[1.0.0]: https://github.com/PathoGenOmics-Lab/snpick/releases/tag/1.0.0
