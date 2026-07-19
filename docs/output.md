---
tags:
  - vcf
  - output
  - phylogenetics
---

# Output formats

## Reduced FASTA

The primary output (`-o`) is a FASTA where every sequence keeps only the **variable columns**,
in alignment order. Headers (ID and description) are preserved verbatim; bases are upper-cased.

Constant and ambiguous-only columns are dropped, so the reduced alignment is a drop-in input for
phylogenetic tools while staying orders of magnitude smaller.

!!! info "No variable sites"
    If the alignment has no variable columns, SNPick writes a valid FASTA with each record's
    header and an empty sequence line, and exits `0`.

### Alternative formats

`--format` writes the reduced alignment as something other than FASTA:

=== "PHYLIP"

    ```bash
    snpick -f alignment.fasta -o snps.phy --format phylip
    ```

    Relaxed PHYLIP (`ntax nchar` header, then `name  sequence`), read by IQ-TREE and RAxML.

=== "NEXUS"

    ```bash
    snpick -f alignment.fasta -o snps.nex --format nexus
    ```

    A NEXUS `DATA` block (`DIMENSIONS` / `FORMAT DATATYPE=DNA` / `MATRIX`) for MrBayes, PAUP*
    and SplitsTree.

## ASC `fconst` for ascertainment-bias correction

Removing invariant sites biases branch lengths unless the model is told how many constant sites
of each base were dropped. SNPick reports these counts on stderr, ready for IQ-TREE's `+ASC`
models:

```text
[snpick] ASC fconst: 744123,1382922,1382180,743556
```

The four numbers are the constant-site counts for **A, C, G, T**. Feed them to IQ-TREE:

```bash
iqtree2 -s snps.fasta -m GTR+ASC -fconst 744123,1382922,1382180,743556
```

## VCF v4.2

With `--vcf` (or `--vcf-output`), SNPick also writes a VCF v4.2 file with one row per variable
site and per-sample genotypes.

```text
##fileformat=VCFv4.2
##source=snpick v1.0.2
##reference=ref
##contig=<ID=1,length=8>
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ref	s1	s2
1	2	.	T	C	.	PASS	NS=3	GT	0	0	1
1	4	.	C	T	.	PASS	NS=3	GT	0	1	0
```

| Field | Meaning |
|---|---|
| `CHROM` | Contig name — `1` by default, override with [`--chrom`](usage.md#vcf-coordinates) |
| `POS` | **1-based alignment column** by default, or the ungapped reference position with `--ref-coords` |
| `REF` | Base of the reference sequence (the first, or the one set by `--reference`); if that base is ambiguous, the first observed base in A, C, G, T order; if it is a gap (under `-g`), `N` |
| `ALT` | The other observed alleles, comma-separated |
| `INFO=NS` | Number of samples with data (a called base; gaps count only under `-g`) |
| `FORMAT=GT` | Per-sample allele index: `0` = REF, `1..` = the *n*-th ALT, `.` = missing/ambiguous |

The `##reference` header carries the **ID** of the reference sequence (the first record, or the
one chosen with `--reference <ID>`), and `##contig` its length.

!!! warning "POS is an alignment coordinate by default"
    Without `--ref-coords`, `POS` is the column index in the alignment, so when the reference
    contains gaps it diverges from the true genomic position and `##contig` length is the
    alignment length. `--ref-coords` maps `POS` (and the contig length) onto ungapped reference
    positions; `--reference <ID>` picks which sequence is that reference.

    Under `--ref-coords`, columns where the reference has a **gap** (insertions relative to the
    reference) have no reference coordinate of their own, so they take the position of the
    preceding reference base. Several such variable columns therefore share one `POS` (and any
    leading-gap column clamps to `POS 1`), producing duplicate-`POS` records — expected for
    insertions, but some downstream tools may need `bcftools norm` or a sort.

### Genotype matrix guard

The genotype matrix is `variants × samples` bytes. To avoid accidental multi-gigabyte VCFs,
SNPick refuses to emit a VCF whose matrix would exceed **4 GB** and tells you to drop `--vcf` or
reduce the input. The reduced FASTA is unaffected by this guard.

### Header-only VCF

If `--vcf` is requested but there are no variable sites, SNPick still writes a valid VCF
containing just the header and sample columns (no data rows), so a Snakemake/Nextflow rule that
declares the `.vcf` as an output does not break.

## Gaps and ambiguous bases

- **Ambiguous bases** (N, R, Y, …) are treated as **missing data**, never as alleles. A column
  is variable only if it has ≥2 standard bases (A, C, G, T). Ambiguous genotypes are written as
  `.` in the VCF.
- **Gaps** (`-`) are **ignored by default**. With `-g` a gap becomes a 5th allele. In the VCF a
  gap in **ALT** is written as `*`, but a gap in **REF** (the reference sequence has a gap at a
  variable column) is written as **`N`**, since VCF v4.2 forbids `-`/`*` in the REF field. The
  sample whose base is that gap still gets genotype `0`, so `REF=N` there marks a gap, not an N.

!!! note "The `*` gap encoding"
    Writing ALT gaps as `*` is an alignment convention shared with snp-sites. Note that in strict
    VCF v4.2, `*` denotes a spanning deletion, so some downstream tools may interpret gap sites
    accordingly. A gap in REF is written as `N` instead (`*` is not a legal REF allele).

With `--iupac-mode resolve`, IUPAC codes (R = A|G, …) are resolved to their bases when
classifying, so a column that is all-A plus one `R` becomes variable. Resolution applies only to
the presence scan; an `R` is still excluded from `NS` and genotyped as missing.

## Machine-readable summary (`--stats-json`)

`--stats-json <FILE>` (or `-` for stdout) writes a flat JSON run summary — so pipelines get the
`fconst` array as a typed field instead of scraping the `[snpick] ASC fconst:` stderr line:

```json
{
  "snpick_version": "1.0.2",
  "input": "alignment.fasta",
  "reference": "H37Rv",
  "sequences": 250,
  "alignment_length": 4411532,
  "variable_sites": 158934,
  "written_sites": 141002,
  "constant_sites": 4252598,
  "constant_by_base": { "A": 744123, "C": 1382922, "G": 1382180, "T": 743556 },
  "ambiguous_sites": 0,
  "fconst": [744123, 1382922, 1382180, 743556],
  "include_gaps": false,
  "threads": 8
}
```

`variable_sites` is the classified count; `written_sites` is what remains after per-site
filtering. Pair with `--dry-run` to compute the summary without writing any alignment.

## Variable-site coordinate map (`--sites-output`)

`--sites-output <TSV>` writes one row per variable site mapping its alignment column to the
ungapped reference position and alleles — useful for cross-referencing gene coordinates or a
masking BED:

```text
alignment_pos	ref_pos	ref	alt
2	1	A	T
8	7	C	A
```

## Composition audit (`--check`)

`--check` prints an A/C/G/T, N, gap, IUPAC and invalid-byte breakdown, then exits — a quick
pre-flight to catch a mis-supplied protein alignment or truncated download (see also
`--on-invalid`).
