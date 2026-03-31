use std::fs::File;
use std::io::{self, BufWriter, Write};

use crate::fasta::FastaRecord;
use crate::types::*;

/// Parameters for variable site extraction (pass 2).
pub struct ExtractParams<'a> {
    pub records: &'a [FastaRecord<'a>],
    pub output: &'a str,
    pub collect_vcf: bool,
    pub lookup: &'a [u8; 256],
    pub upper: &'a [u8; 256],
    pub layout: SeqLayout,
}

/// Pass 2: extract variable sites from alignment and write output FASTA.
///
/// For single-line FASTA: direct byte access via `data[seq_offset + pos]`.
/// For multi-line: linear scan per record, skipping newlines.
///
/// Returns VCF genotype matrix if `collect_vcf` is true.
pub fn pass2_extract(
    data: &[u8], var_positions: &mut [VariablePosition], params: &ExtractParams<'_>,
) -> io::Result<Option<Vec<u8>>> {
    let ExtractParams { records, output, collect_vcf, lookup, upper, layout } = params;
    let collect_vcf = *collect_vcf;
    let layout = *layout;
    let num_var = var_positions.len();
    let num_samples = records.len();
    let pos_indices: Vec<usize> = var_positions.iter().map(|v| v.index).collect();

    let out_file = File::create(output).map_err(|e| io::Error::new(e.kind(),
        format!("Cannot create output '{}': {}", output, e)))?;
    let mut writer = BufWriter::with_capacity(IO_BUF, out_file);

    let mut vcf_geno: Vec<u8> = if collect_vcf { vec![0u8; num_var * num_samples] } else { Vec::new() };
    let mut ns_counts: Vec<usize> = if collect_vcf { vec![0usize; num_var] } else { Vec::new() };
    let mut var_buf = vec![0u8; num_var];

    for (si, rec) in records.iter().enumerate() {
        if layout.single_line {
            let base = rec.seq_offset;
            for (vi, &p) in pos_indices.iter().enumerate() {
                var_buf[vi] = upper[data[base + p] as usize];
            }
        } else {
            let mut pos = rec.seq_offset;
            let end = data.len();
            let mut base_idx = 0usize;
            let mut var_idx = 0usize;
            while var_idx < num_var && pos < end {
                let b = data[pos];
                pos += 1;
                if b == b'\n' || b == b'\r' { continue; }
                if base_idx == pos_indices[var_idx] {
                    var_buf[var_idx] = upper[b as usize];
                    var_idx += 1;
                }
                base_idx += 1;
            }
        }

        writer.write_all(b">")?;
        writer.write_all(rec.id)?;
        if !rec.desc.is_empty() {
            writer.write_all(b" ")?;
            writer.write_all(rec.desc)?;
        }
        writer.write_all(b"\n")?;
        writer.write_all(&var_buf)?;
        writer.write_all(b"\n")?;

        if collect_vcf {
            for (vi, &nuc) in var_buf.iter().enumerate() {
                vcf_geno[vi * num_samples + si] = nuc;
                if lookup[nuc as usize] != 0 {
                    ns_counts[vi] += 1;
                }
            }
        }
    }

    writer.flush()?;

    if collect_vcf {
        for (vi, vp) in var_positions.iter_mut().enumerate() {
            vp.ns = ns_counts[vi];
        }
    }

    eprintln!("[snpick] Pass 2: Wrote {} sequences to {}.", num_samples, output);
    if collect_vcf { Ok(Some(vcf_geno)) } else { Ok(None) }
}
