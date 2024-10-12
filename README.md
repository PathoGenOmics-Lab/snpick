<p align="center">
  <img src="logo/SNPick.png" title="SNPick.png logo" style="width:750px; height: auto;">
</p>

</div>

__Paula Ruiz-Rodriguez<sup>1</sup>__ 
__and Mireia Coscolla<sup>1</sup>__
<br>
<sub> 1. Institute for Integrative Systems Biology, I<sup>2</sup>SysBio, University of Valencia-CSIC, Valencia, Spain </sub>  

# SNPick

SNPick is a fast and memory-efficient tool designed to extract SNP (Single Nucleotide Polymorphism) sites from large-scale FASTA alignments, taking into consideration IUPAC ambiguous nucleotides. SNP analysis is critical for understanding genetic variation, phylogenetic relationships, and evolutionary biology.

## Problematic

Tools like **snp-sites** have been widely used to extract SNP positions from FASTA alignments. **snp-sites** is a valuable tool for smaller datasets, providing an effective solution for extracting SNPs in many research scenarios. However, it faces significant challenges when processing large alignments or a high number of sequences.

When dealing with thousands of sequences or large sequence lengths, **snp-sites** tends to become inefficient both in terms of runtime and memory usage. For instance, the tool may struggle or even fail entirely when processing datasets with hundreds of thousands of sequences or sequences that are several megabases in length. This limitation is particularly relevant for researchers who need to work with whole-genome alignments from many individuals, especially in epidemiological studies or population genomics.

## SNPick Advantage

**SNPick** was developed to overcome these limitations, providing a highly scalable and optimized approach to SNP extraction. Unlike **snp-sites**, SNPick employs parallel processing and an efficient memory management strategy, making it suitable for very large datasets. This scalability allows users to work with massive alignments, extracting SNPs in a fraction of the time and with significantly lower memory requirements compared to traditional tools.

The key benefits of SNPick include:
- **Scalability**: Able to handle datasets with tens of thousands of sequences and large genome alignments without a significant increase in memory consumption or processing time.
- **Efficiency**: Optimized for speed, taking advantage of multi-threading and efficient memory allocation.
- **Versatility**: Supports the extraction of SNPs while handling IUPAC ambiguous nucleotides, making it highly adaptable to different types of genomic data.

## Usage

SNPick can be easily integrated into existing pipelines for genomic data analysis. It is suitable for use in environments where memory efficiency and processing speed are crucial, such as analyzing large microbial datasets or conducting comparative genomic studies.

### Example

To extract SNPs from a given FASTA file:

```sh
snpick --fasta input_alignment.fasta --output snp_alignment.fasta --threads 8
```

**Input**: The `--fasta` parameter specifies the input alignment file in FASTA format. This file contains the aligned sequences from which SNPs will be extracted.

**Output**: The `--output` parameter specifies the output file, which will contain only the SNP sites extracted from the input alignment.

This command will generate an output file containing only the SNP sites, optimized for downstream analysis.

### Input and Output Example

**Input File (`input_alignment.fasta`)**:
```fasta
>sequence1
ATGCTAGCTAGCTAGCTA
>sequence2
ATGCTAGCTGGCTAGCTA
>sequence3
ATGCTAGCTAGCTAGCTA
```

**Output File (`snp_alignment.fasta`)**:
```fasta
>sequence1
T
>sequence2
G
>sequence3
T
```

In this example, the SNP site is at the position where `sequence2` has a `G` while the others have a `T`.

### Another Example

If you have a large dataset and want to utilize all available CPU cores for faster processing, you can use the following command:

```sh
snpick --fasta large_dataset.fasta --output snp_large_output.fasta --threads 16
```

In this example, `large_dataset.fasta` is the input file, `snp_large_output.fasta` is the output file, and `--threads 16` specifies the use of 16 CPU threads to speed up the SNP extraction process.

