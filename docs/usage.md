---
tags:
  - cli
  - usage
---

# Usage

```
snpick [OPTIONS] --fasta <FASTA> --output <OUTPUT>
```

## Options

| Argument | Required | Description |
|---|:---:|---|
| `-f, --fasta <FILE>` | :material-check: | Input FASTA alignment (all sequences must have equal length) |
| `-o, --output <FILE>` | :material-check: | Output FASTA containing only the variable sites |
| `-g, --include-gaps` | | Treat gaps (`-`) as a 5th character instead of ignoring them |
| `--vcf` | | Also write a VCF, named after the output (`snps.fasta` → `snps.vcf`) |
| `--vcf-output <FILE>` | | Write the VCF to a custom path (implies `--vcf`) |
| `-t, --threads <N>` | | Threads for the parallel scan (default: all logical cores) |
| `--chrom <NAME>` | | `CHROM` / contig name written to the VCF (default: `1`) |
| `-q, --quiet` | | Silence progress logs on stderr (errors are still reported) |
| `-h, --help` | | Print help |
| `-V, --version` | | Print version |

## Examples

### Basic extraction

```bash
snpick -f alignment.fasta -o snps.fasta
```

**Input** (`alignment.fasta`):

```text
>sequence1
ATGCTAGCTAGCTAGCTA
>sequence2
ATGCTAGCTGGCTAGCTA
>sequence3
ATGCTAGCTAGCTAGCTA
```

**Output** (`snps.fasta`):

```text
>sequence1
A
>sequence2
G
>sequence3
A
```

**stderr:**

```text
[snpick] Mapped 63 bytes. 3 sequences × 18 positions.
[snpick] 1 variable, 17 constant (A:4 C:4 G:4 T:5), 0 ambiguous-only, 18 total.
[snpick] ASC fconst: 4,4,4,5
[snpick] Done in 0.00s. 1 vars from 3 seqs × 18 pos.
```

### With a VCF

```bash
snpick -f alignment.fasta -o snps.fasta --vcf
# or a custom path:
snpick -f alignment.fasta -o snps.fasta --vcf-output variants.vcf
```

See [Output formats](output.md) for the VCF layout.

### Set the VCF contig name

By default the VCF `CHROM` and `##contig` are `1`. Match your reference so downstream tools
(bcftools, GATK, IGV) line up without post-processing:

```bash
snpick -f alignment.fasta -o snps.fasta --vcf --chrom NC_000962.3 # (1)!
```

1. `--chrom` sets **both** the `##contig` header ID and the per-row `CHROM` column, and is
   rejected if it contains whitespace (which would break the tab-delimited columns).

### Include gaps

```bash
snpick -f alignment.fasta -o snps.fasta -g
```

Without `-g`, gap columns never make a site variable. With `-g`, a gap is treated as a 5th
allele (rendered as `*` in the VCF). See [Gaps & ambiguity](output.md#gaps-and-ambiguous-bases).

### Control threads (HPC / reproducibility)

The scan is parallelised with Rayon and, by default, uses every logical core. On a shared
SLURM node, pin it to your allocation:

```bash
snpick -f alignment.fasta -o snps.fasta -t 8
```

!!! note "Determinism"
    The thread count **never changes the output** — the per-position bitmask is merged with a
    commutative OR — only the wall-clock time.

### Quiet mode for pipelines

```bash
snpick -f alignment.fasta -o snps.fasta --vcf -q
```

All progress goes to **stderr** and `stdout` stays clean, so `-q` is only needed to silence the
`[snpick]` chatter; errors are always printed.

## Exit codes

| Code | Meaning |
|:---:|---|
| `0` | Success (including the "no variable sites" case) |
| `1` | Error — bad input, unequal sequence lengths, unwritable output, etc. |
