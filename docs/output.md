# Output formats

## Reduced FASTA

The primary output (`-o`) is a FASTA where every sequence keeps only the **variable columns**,
in alignment order. Headers (ID and description) are preserved verbatim; bases are upper-cased.

Constant and ambiguous-only columns are dropped, so the reduced alignment is a drop-in input for
phylogenetic tools while staying orders of magnitude smaller.

!!! info "No variable sites"
    If the alignment has no variable columns, SNPick writes a valid FASTA with each record's
    header and an empty sequence line, and exits `0`.

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
##reference=first_sequence
##contig=<ID=1,length=8>
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ref	s1	s2
1	2	.	T	C	.	PASS	NS=3	GT	0	0	1
1	4	.	C	T	.	PASS	NS=3	GT	0	1	0
```

| Field | Meaning |
|---|---|
| `CHROM` | Contig name — `1` by default, override with [`--chrom`](usage.md#set-the-vcf-contig-name) |
| `POS` | **1-based alignment column** (not an ungapped reference coordinate) |
| `REF` | Base of the first sequence; if that base is ambiguous, the first observed base in A, C, G, T order |
| `ALT` | The other observed alleles, comma-separated |
| `INFO=NS` | Number of samples with data (a called base; gaps count only under `-g`) |
| `FORMAT=GT` | Per-sample allele index: `0` = REF, `1..` = the *n*-th ALT, `.` = missing/ambiguous |

!!! warning "POS is an alignment coordinate"
    `POS` is the column index in the alignment, so when the reference sequence contains gaps it
    diverges from the true genomic position, and `##contig` length is the alignment length.

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
- **Gaps** (`-`) are **ignored by default**. With `-g` a gap becomes a 5th allele and, in the
  VCF, is rendered as `*`.

!!! note "The `*` gap encoding"
    Writing gaps as `*` is an alignment convention shared with snp-sites. Note that in strict
    VCF v4.2, `*` denotes a spanning deletion, so some downstream tools may interpret gap sites
    accordingly.
