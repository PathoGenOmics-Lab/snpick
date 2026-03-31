use std::fs::File;
use std::io::{self, BufWriter, Write};

use crate::fasta::FastaRecord;
use crate::types::VariablePosition;

/// Write VCF output from genotype matrix and variable positions.
pub fn write_vcf(
    vcf_geno: &[u8], num_samples: usize, var_positions: &[VariablePosition],
    vcf_path: &str, records: &[FastaRecord], seq_length: usize,
) -> io::Result<()> {
    let out = File::create(vcf_path).map_err(|e| io::Error::new(e.kind(),
        format!("Cannot create VCF '{}': {}", vcf_path, e)))?;
    let mut w = BufWriter::with_capacity(4 * 1024 * 1024, out);

    // Header
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##source=snpick v{}", env!("CARGO_PKG_VERSION"))?;
    writeln!(w, "##reference=first_sequence")?;
    writeln!(w, "##contig=<ID=1,length={}>", seq_length)?;
    writeln!(w, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">")?;
    writeln!(w, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    write!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for rec in records {
        write!(w, "\t")?;
        w.write_all(rec.id)?;
    }
    writeln!(w)?;

    // Data rows
    let mut lut = [255u8; 256];
    for (vi, vp) in var_positions.iter().enumerate() {
        let alt: String = vp.alt_bases.iter()
            .map(|&b| if b == b'-' { "*".to_string() } else { (b as char).to_string() })
            .collect::<Vec<_>>().join(",");

        // Build allele → index LUT for this position
        lut[vp.ref_base as usize] = 0;
        for (i, &ab) in vp.alt_bases.iter().enumerate() {
            lut[ab as usize] = (i + 1) as u8;
        }

        write!(w, "1\t{}\t.\t{}\t{}\t.\tPASS\tNS={}\tGT",
            vp.index + 1, vp.ref_base as char, alt, vp.ns)?;

        let row = vi * num_samples;
        for si in 0..num_samples {
            let idx = lut[vcf_geno[row + si] as usize];
            if idx == 255 { write!(w, "\t.")?; } else { write!(w, "\t{}", idx)?; }
        }
        writeln!(w)?;

        // Reset LUT entries
        lut[vp.ref_base as usize] = 255;
        for &ab in &vp.alt_bases { lut[ab as usize] = 255; }
    }

    w.flush()
}
