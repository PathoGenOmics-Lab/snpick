---
tags:
  - cli
  - usage
---

# Usage

```
snpick [OPTIONS] --fasta <FASTA>
```

`-f/--fasta` is always required. `-o/--output` is required **unless** `--dry-run` or `--check`
is given. Use `-` in place of a path for `--fasta`, `--output` or `--stats-json` to read from
stdin / write to stdout. Run `snpick --help` for the authoritative, always-current list.

## Core

| Argument | Description |
|---|---|
| `-f, --fasta <FILE>` | Input FASTA alignment (`-` = stdin; gzip/bgzip auto-detected) |
| `-o, --output <FILE>` | Output FASTA of variable sites (`-` = stdout) |
| `-g, --include-gaps` | Treat gaps (`-`) as a 5th character |
| `--format <FMT>` | Reduced-alignment format: `fasta` (default), `phylip`, `nexus` |
| `-t, --threads <N>` | Threads for the parallel scan (default: all logical cores) |
| `-q, --quiet` | Silence progress logs (errors still shown) |
| `-v`, `-vv` | Increase detail (diagnostics; allele-class histogram) |

## VCF & coordinates

| Argument | Description |
|---|---|
| `--vcf` | Write a VCF named after the output (`snps.fasta` → `snps.vcf`) |
| `--vcf-output <FILE>` | Write the VCF to a custom path (implies `--vcf`) |
| `--chrom <NAME>` | `CHROM` / `##contig` name (default `1`, e.g. `NC_000962.3`) |
| `--reference <ID>` | Sequence ID used for REF/polarity (default: first) |
| `--ref-coords` | Write `POS` as the ungapped reference position, not the alignment column |
| `--sites-output <TSV>` | Map each site: `alignment_pos`, `ref_pos`, `ref`, `alt` |

## Filtering & masking

| Argument | Description |
|---|---|
| `--max-missing <f>` | Drop sites with more than fraction `f` of missing genotypes |
| `--mac <n>` | Drop sites with minor-allele count below `n` (e.g. `2` drops singletons) |
| `--maf <f>` | Drop sites with minor-allele frequency below `f` |
| `--min-samples <n>` | Drop sites with fewer than `n` samples with data |
| `--max-alleles <n>` | Drop sites with more than `n` distinct alleles (`2` = biallelic) |
| `--keep-samples <IDs>` | Keep only these samples (comma-separated, or `@file`) |
| `--exclude-samples <IDs>` | Drop these samples (mutually exclusive with `--keep-samples`) |
| `--mask <BED>` | Mask BED regions (excluded from output **and** `fconst`) |
| `--mask-ref` | Interpret `--mask` coordinates as reference positions |

!!! info "Filtering keeps ASC valid"
    Per-site filters remove **variable** sites from the output only; they are never reclassified
    as constant, so the `fconst` counts are untouched. Sample selection filters *before* the
    scan, so `fconst` is recomputed for exactly the retained samples.

## QC, reporting & robustness

| Argument | Description |
|---|---|
| `--stats-json <FILE>` | Machine-readable JSON summary (`-` = stdout) |
| `--dry-run` | Report statistics without writing any FASTA/VCF |
| `--check` | Audit the alignment composition and exit |
| `--on-invalid <MODE>` | `ignore` (default), `warn`, or `error` on out-of-alphabet bytes |
| `--iupac-mode <MODE>` | `missing` (default) or `resolve` IUPAC codes to bases |
| `--allow-dup-ids` | Permit duplicate sequence IDs instead of erroring |

**Exit codes:** `0` success · `1` I/O error · `2` bad input/data.

## Examples

### Basic extraction

```bash
snpick -f alignment.fasta -o snps.fasta
```

### VCF with reference coordinates and a stats sidecar

```bash
snpick -f alignment.fasta.gz -o snps.fasta \
       --vcf --chrom NC_000962.3 --reference H37Rv --ref-coords \
       --stats-json stats.json    # (1)!
```

1. `--ref-coords` sets `POS` (and the contig length) to ungapped reference positions;
   `--reference` pins REF polarity to H37Rv; `--stats-json` gives IQ-TREE the `fconst` array
   without scraping stderr.

### Filter and mask

```bash
# Drop singletons, keep biallelic sites, mask repetitive regions (reference coords)
snpick -f alignment.fasta -o snps.fasta \
       --mac 2 --max-alleles 2 --mask exclude.bed --mask-ref
```

### Subset samples

```bash
snpick -f alignment.fasta -o clade.fasta --keep-samples @clade_ids.txt
```

### Alternative formats and piping

```bash
snpick -f alignment.fasta -o snps.nex --format nexus
gzip -dc big.fasta.gz | snpick -f - -o - --vcf-output snps.vcf > snps.fasta
```

### Dry run and composition audit

```bash
snpick -f alignment.fasta --dry-run --stats-json -   # stats only, JSON to stdout
snpick -f alignment.fasta --check                    # composition breakdown, then exit
```
