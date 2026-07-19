//! VCF output writer.
//!
//! Generates VCF v4.2 from the genotype matrix built during pass 2.
//! Uses a per-position lookup table for O(1) allele → index mapping.

use std::io::{self, BufWriter, Write};

use crate::extract::open_sink;
use crate::fasta::FastaRecord;
use crate::types::{VariablePosition, IO_BUF};

/// Write VCF output from genotype matrix and variable positions.
#[allow(clippy::too_many_arguments)]
pub fn write_vcf(
    vcf_geno: &[u8], num_samples: usize, var_positions: &[VariablePosition],
    vcf_path: &str, records: &[FastaRecord], seq_length: usize, chrom: &str, reference: &str,
    pos_map: Option<&[u32]>,
) -> io::Result<()> {
    // `-` streams to stdout, like -o / --sites-output / --stats-json; otherwise a file.
    let out = open_sink(vcf_path)?;
    let mut w = BufWriter::with_capacity(IO_BUF, out);

    // Header
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##source=snpick v{}", env!("CARGO_PKG_VERSION"))?;
    writeln!(w, "##reference={}", reference)?;
    // With reference coordinates the contig length is the ungapped reference length.
    let contig_len = match pos_map {
        Some(rp) => *rp.last().unwrap_or(&(seq_length as u32)) as usize,
        None => seq_length,
    };
    writeln!(w, "##contig=<ID={},length={}>", chrom, contig_len)?;
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

        // VCF v4.2 REF must be A/C/G/T/N — never '-' and never '*' (which is an ALT-only
        // spanning-deletion allele that bcftools/htslib reject in REF). Render a gap reference
        // as 'N'. (ALT gaps stay '*', the snp-sites convention.)
        let ref_byte = if vp.ref_base == b'-' { b'N' } else { vp.ref_base };

        row.clear();
        let pos = match pos_map {
            Some(rp) => rp[vp.index].max(1),
            None => (vp.index + 1) as u32,
        };
        write!(row, "{}\t{}\t.\t", chrom, pos)?;
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
