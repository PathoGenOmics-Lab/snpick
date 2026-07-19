# Changelog

All notable changes to SNPick are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.2] - 2026-07-18

### Added

**Filtering & selection**
- Per-site filtering: `--max-missing`, `--mac`, `--maf`, `--min-samples`,
  `--max-alleles`. Filtered sites stay variable (never re-enter `fconst`), so ASC
  correction remains valid.
- `--keep-samples` / `--exclude-samples` (comma-separated IDs or `@file`) subset
  the panel before the scan, so `fconst` is recomputed for the retained samples.
- `--mask <BED>` (and `--mask-ref`) excludes regions from both the output and `fconst`.

**Coordinates & references**
- `--reference <ID>` chooses the REF/polarity sequence and the VCF `##reference`.
- `--ref-coords` writes VCF `POS` as ungapped reference positions.
- `--sites-output <TSV>` maps each site's alignment column to its reference position.

**I/O & formats**
- Transparent gzip/bgzip input, plus stdin (`-f -`) and stdout (`-o -`) streaming.
- `--format {fasta,phylip,nexus}` for the reduced alignment.
- `--stats-json <FILE>` (or `-`) writes a typed JSON run summary with the `fconst` array.

**Robustness & QC**
- `--dry-run` (stats only) and `--check` (composition audit) exit without writing.
- `--on-invalid {ignore,warn,error}` guards against non-nucleotide input.
- `--iupac-mode resolve` optionally resolves IUPAC codes to their bases.
- Empty/duplicate sequence IDs are rejected (`--allow-dup-ids` to permit).
- `-v`/`-vv` diagnostics; distinct exit codes (`2` bad input, `1` I/O).

**Threading & VCF**
- `-t, --threads <N>` pins the Rayon thread pool (deterministic wall-clock; output
  is unaffected by thread count).
- `--chrom <NAME>` sets the VCF `CHROM` / `##contig` name (default `1`).
- `-q, --quiet` silences the `[snpick]` progress logs (errors still print).
- A valid header-only VCF is written when `--vcf` is requested but there are no
  variable sites, so pipelines declaring the `.vcf` output no longer break.

**Packaging**
- Published as a reusable Rust library crate, with Python (pyo3/maturin) and
  WebAssembly bindings, and Docker/Apptainer container recipes.

### Fixed
- Silent SNP loss when a single-line FASTA record contained a blank line: such
  records were read with misaligned byte offsets. They now use the
  newline-skipping scanner.
- Gap reference bases produced an invalid VCF `REF` of `-` under `--include-gaps`;
  the `REF` is now rendered as `N` (a valid VCF base). ALT gaps stay `*` (the snp-sites
  convention), so `REF` and `ALT` gap encodings differ by design.
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

## [1.0.1] - 2026-03-31

- First Bioconda release.
- Cross-platform Build & Release CI workflow (Linux/macOS, x86_64/aarch64).
- README overhaul and benchmarks.

## [1.0.0] - 2024-11-16

- Initial release: zero-copy memory-mapped extraction of variable sites from
  FASTA alignments, optional VCF v4.2 output, and ASC `fconst` reporting.

[1.0.2]: https://github.com/PathoGenOmics-Lab/snpick/compare/1.0.1...1.0.2
[1.0.1]: https://github.com/PathoGenOmics-Lab/snpick/compare/1.0.0...1.0.1
[1.0.0]: https://github.com/PathoGenOmics-Lab/snpick/releases/tag/1.0.0
