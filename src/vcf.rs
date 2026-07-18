//! VCF output writer.
//!
//! Generates VCF v4.2 from the genotype matrix built during pass 2.
//! Uses a per-position lookup table for O(1) allele → index mapping.

use std::fs::File;
use std::io::{self, BufWriter, Write};

use crate::fasta::FastaRecord;
use crate::types::{VariablePosition, IO_BUF};

/// Write VCF output from genotype matrix and variable positions.
pub fn write_vcf(
    vcf_geno: &[u8], num_samples: usize, var_positions: &[VariablePosition],
    vcf_path: &str, records: &[FastaRecord], seq_length: usize, chrom: &str,
) -> io::Result<()> {
    let out = File::create(vcf_path).map_err(|e| io::Error::new(e.kind(),
        format!("Cannot create VCF '{}': {}", vcf_path, e)))?;
    let mut w = BufWriter::with_capacity(IO_BUF, out);

    // Header
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##source=snpick v{}", env!("CARGO_PKG_VERSION"))?;
    writeln!(w, "##reference=first_sequence")?;
    writeln!(w, "##contig=<ID={},length={}>", chrom, seq_length)?;
    writeln!(w, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">")?;
    writeln!(w, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    write!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for rec in records {
        write!(w, "\t")?;
        w.write_all(rec.id)?;
    }
    writeln!(w)?;

    // Data rows. Each row is assembled in a reused byte buffer, so the
    // num_var × num_samples genotypes are single byte pushes rather than one
    // formatted write!() per field. Output bytes are identical to the naive form.
    let mut lut = [255u8; 256];
    let mut row: Vec<u8> = Vec::with_capacity(64 + num_samples * 2);
    let mut alt: Vec<u8> = Vec::with_capacity(8);
    for (vi, vp) in var_positions.iter().enumerate() {
        // ALT alleles, gap → '*'.
        alt.clear();
        for (i, &b) in vp.alt_bases.iter().enumerate() {
            if i > 0 { alt.push(b','); }
            alt.push(if b == b'-' { b'*' } else { b });
        }

        // Build allele → index LUT for this position.
        lut[vp.ref_base as usize] = 0;
        for (i, &ab) in vp.alt_bases.iter().enumerate() {
            lut[ab as usize] = (i + 1) as u8;
        }

        // A gap reference renders as '*' (VCF v4.2 REF must be A/C/G/T/N, never '-'),
        // consistent with the '-'→'*' mapping applied to ALT.
        let ref_byte = if vp.ref_base == b'-' { b'*' } else { vp.ref_base };

        row.clear();
        write!(row, "{}\t{}\t.\t", chrom, vp.index + 1)?;
        row.push(ref_byte);
        row.push(b'\t');
        row.extend_from_slice(&alt);
        write!(row, "\t.\tPASS\tNS={}\tGT", vp.ns)?;

        let base = vi * num_samples;
        for si in 0..num_samples {
            row.push(b'\t');
            let idx = lut[vcf_geno[base + si] as usize];
            // idx is 0..=4 (ref + at most 4 alts), so a single ASCII digit.
            if idx == 255 { row.push(b'.'); } else { row.push(b'0' + idx); }
        }
        row.push(b'\n');
        w.write_all(&row)?;

        // Reset LUT entries.
        lut[vp.ref_base as usize] = 255;
        for &ab in &vp.alt_bases { lut[ab as usize] = 255; }
    }

    w.flush()
}
