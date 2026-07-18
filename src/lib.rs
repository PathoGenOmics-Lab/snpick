//! # SNPick
//!
//! Fast, memory-efficient extraction of variable (SNP) sites from FASTA alignments.
//!
//! This library exposes the core of the `snpick` command-line tool as reusable modules:
//!
//! - [`fasta`] — zero-copy FASTA indexing over memory-mapped data.
//! - [`scan`] — pass 1: the parallel bitmask scan and site classification ([`scan::analyze`]).
//! - [`extract`] — pass 2: reduced FASTA / genotype-matrix extraction.
//! - [`vcf`] — VCF v4.2 writer.
//! - [`filter`] — per-site allele/missingness counting and filtering.
//! - [`coords`] — reference-anchored coordinates and BED masking.
//! - [`audit`] — alignment composition audit.
//! - [`input`] — transparent gzip/stdin input mapping.
//! - [`types`] — shared constants, lookup tables and data structures.
//!
//! # Example
//!
//! ```no_run
//! use snpick::fasta::index_fasta;
//! use snpick::scan::{analyze, pass1_scan};
//! use snpick::types::build_lookup;
//!
//! let data: &[u8] = b">a\nACGT\n>b\nACGA\n";
//! let lookup = build_lookup(false);
//! let (records, seq_len, layout) = index_fasta(data).unwrap();
//! let bitmask = pass1_scan(data, &records, seq_len, layout, &lookup);
//! let ref_seq: Vec<u8> = data[3..7].to_vec();
//! let (variable_sites, counts) = analyze(&bitmask, &ref_seq, &lookup, false);
//! assert_eq!(variable_sites.len(), 1);
//! assert_eq!(counts.constant.total(), 3);
//! ```

pub mod audit;
pub mod coords;
pub mod extract;
pub mod fasta;
pub mod filter;
pub mod input;
pub mod scan;
pub mod types;
pub mod vcf;
