//! Pass 2: variable site extraction and output FASTA generation.
//!
//! Reads only the variable positions from each sequence (sparse access)
//! and writes a reduced FASTA. Optionally collects a genotype matrix for VCF.

use std::fs::File;
use std::io::{self, BufWriter, Write};

use crate::fasta::FastaRecord;
use crate::types::*;

/// Write a NEXUS taxon label, single-quoting (and doubling embedded quotes) when it contains
/// characters NEXUS treats as punctuation/whitespace, so labels like `iso[London]` survive.
fn write_nexus_label(w: &mut impl Write, id: &[u8]) -> io::Result<()> {
    let safe = !id.is_empty()
        && id.iter().all(|&b| b.is_ascii_alphanumeric() || b == b'_' || b == b'.');
    if safe {
        return w.write_all(id);
    }
    w.write_all(b"'")?;
    for &b in id {
        if b == b'\'' {
            w.write_all(b"''")?;
        } else {
            w.write_all(&[b])?;
        }
    }
    w.write_all(b"'")
}

/// Open an output sink: a file, or stdout when the path is "-".
pub fn open_sink(path: &str) -> io::Result<Box<dyn Write>> {
    if path == "-" {
        Ok(Box::new(io::stdout().lock()))
    } else {
        let f = File::create(path).map_err(|e| io::Error::new(e.kind(),
            format!("Cannot create output '{}': {}", path, e)))?;
        Ok(Box::new(f))
    }
}

/// Parameters for variable site extraction (pass 2).
pub struct ExtractParams<'a> {
    pub records: &'a [FastaRecord<'a>],
    pub output: &'a str,
    pub collect_vcf: bool,
    pub lookup: &'a [u8; 256],
    pub upper: &'a [u8; 256],
    pub layout: SeqLayout,
    pub format: OutputFormat,
    pub progress: bool,
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
    let ExtractParams { records, output, collect_vcf, lookup, upper, layout, format, progress } = params;
    let collect_vcf = *collect_vcf;
    let layout = *layout;
    let format = *format;
    let progress = *progress;
    let num_var = var_positions.len();
    let num_samples = records.len();
    let pos_indices: Vec<usize> = var_positions.iter().map(|v| v.index).collect();

    let mut writer = BufWriter::with_capacity(IO_BUF, open_sink(output)?);

    // Format preamble.
    match format {
        OutputFormat::Fasta => {}
        OutputFormat::Phylip => writeln!(writer, "{} {}", num_samples, num_var)?,
        OutputFormat::Nexus => {
            writeln!(writer, "#NEXUS")?;
            writeln!(writer, "BEGIN DATA;")?;
            writeln!(writer, "  DIMENSIONS NTAX={} NCHAR={};", num_samples, num_var)?;
            writeln!(writer, "  FORMAT DATATYPE=DNA MISSING=N GAP=-;")?;
            writeln!(writer, "  MATRIX")?;
        }
    }

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

        match format {
            OutputFormat::Fasta => {
                writer.write_all(b">")?;
                writer.write_all(rec.id)?;
                if !rec.desc.is_empty() {
                    writer.write_all(b" ")?;
                    writer.write_all(rec.desc)?;
                }
                writer.write_all(b"\n")?;
                writer.write_all(&var_buf)?;
                writer.write_all(b"\n")?;
            }
            OutputFormat::Phylip => {
                writer.write_all(rec.id)?;
                writer.write_all(b"  ")?;
                writer.write_all(&var_buf)?;
                writer.write_all(b"\n")?;
            }
            OutputFormat::Nexus => {
                writer.write_all(b"    ")?;
                write_nexus_label(&mut writer, rec.id)?;
                writer.write_all(b"  ")?;
                writer.write_all(&var_buf)?;
                writer.write_all(b"\n")?;
            }
        }

        if collect_vcf {
            for (vi, &nuc) in var_buf.iter().enumerate() {
                vcf_geno[vi * num_samples + si] = nuc;
                if lookup[nuc as usize] != 0 {
                    ns_counts[vi] += 1;
                }
            }
        }

        if progress && (si % 256 == 0) {
            eprint!("\r[snpick] Extracting sequence {}/{}", si + 1, num_samples);
        }
    }
    if progress {
        eprintln!("\r[snpick] Extracted {} sequences.            ", num_samples);
    }

    if let OutputFormat::Nexus = format {
        writeln!(writer, "  ;")?;
        writeln!(writer, "END;")?;
    }

    writer.flush()?;

    if collect_vcf {
        for (vi, vp) in var_positions.iter_mut().enumerate() {
            vp.ns = ns_counts[vi];
        }
    }

    if collect_vcf { Ok(Some(vcf_geno)) } else { Ok(None) }
}
