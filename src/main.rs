use bio::io::fasta::{self, FastaRead};
use clap::Parser;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

// =============================================================================
// CLI
// =============================================================================

#[derive(Parser, Debug)]
#[command(
    name = "snpick",
    version = env!("CARGO_PKG_VERSION"),
    author = "Paula Ruiz-Rodriguez <paula.ruiz.rodriguez@csic.es>",
    about = "A fast, memory-efficient tool for extracting variable sites from FASTA alignments."
)]
struct Args {
    #[arg(short, long)] fasta: String,
    #[arg(short, long)] output: String,
    #[arg(short = 'g', long)] include_gaps: bool,
    #[arg(long)] vcf: bool,
    #[arg(long)] vcf_output: Option<String>,
}

// =============================================================================
// Bitmask
// =============================================================================

const BIT_A: u8 = 0b00001;
const BIT_C: u8 = 0b00010;
const BIT_G: u8 = 0b00100;
const BIT_T: u8 = 0b01000;
const BIT_GAP: u8 = 0b10000;
const MAX_SEQ_LENGTH: usize = 10_000_000_000;
const IO_BUF: usize = 16 * 1024 * 1024; // 16 MB I/O buffers

fn build_lookup(include_gaps: bool) -> [u8; 256] {
    let mut t = [0u8; 256];
    t[b'A' as usize] = BIT_A; t[b'a' as usize] = BIT_A;
    t[b'C' as usize] = BIT_C; t[b'c' as usize] = BIT_C;
    t[b'G' as usize] = BIT_G; t[b'g' as usize] = BIT_G;
    t[b'T' as usize] = BIT_T; t[b't' as usize] = BIT_T;
    if include_gaps { t[b'-' as usize] = BIT_GAP; }
    t
}

fn bits_to_bases(bits: u8, include_gaps: bool) -> Vec<u8> {
    let mut v = Vec::with_capacity(5);
    if bits & BIT_A != 0 { v.push(b'A'); }
    if bits & BIT_C != 0 { v.push(b'C'); }
    if bits & BIT_G != 0 { v.push(b'G'); }
    if bits & BIT_T != 0 { v.push(b'T'); }
    if include_gaps && bits & BIT_GAP != 0 { v.push(b'-'); }
    v
}

fn check_paths_differ(input: &str, output: &str) -> io::Result<()> {
    let inp = std::fs::canonicalize(input)?;
    let out_path = Path::new(output);
    let out_abs = if out_path.exists() { std::fs::canonicalize(output)? }
    else {
        let parent = out_path.parent().unwrap_or(Path::new("."));
        std::fs::canonicalize(parent).unwrap_or(parent.to_path_buf())
            .join(out_path.file_name().unwrap_or_default())
    };
    if inp == out_abs {
        return Err(io::Error::new(io::ErrorKind::InvalidInput,
            format!("Input and output resolve to same file: {}", inp.display())));
    }
    Ok(())
}

// =============================================================================
// Data structures
// =============================================================================

struct VariablePosition { index: usize, ref_base: u8, alt_bases: Vec<u8> }
struct ConstantSiteCounts { a: usize, c: usize, g: usize, t: usize }
impl std::fmt::Display for ConstantSiteCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "A:{} C:{} G:{} T:{}", self.a, self.c, self.g, self.t)
    }
}
impl ConstantSiteCounts {
    fn total(&self) -> usize { self.a + self.c + self.g + self.t }
    fn fconst(&self) -> String { format!("{},{},{},{}", self.a, self.c, self.g, self.t) }
}

// =============================================================================
// Pass 1: streaming bitmask — O(L) memory, bio reader
// =============================================================================

struct ScanResult {
    bitmask: Vec<u8>,
    ref_seq: Vec<u8>,
    names: Vec<String>,
    seq_length: usize,
}

fn pass1_scan(path: &str, lookup: &[u8; 256]) -> io::Result<ScanResult> {
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(IO_BUF, file);
    let mut fasta = fasta::Reader::new(reader);
    let mut record = fasta::Record::new();

    let mut bitmask: Vec<u8> = Vec::new();
    let mut ref_seq: Vec<u8> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut seq_length = 0usize;
    let mut count = 0usize;

    loop {
        fasta.read(&mut record)?;
        if record.is_empty() { break; }
        let seq = record.seq();

        if count == 0 {
            seq_length = seq.len();
            if seq_length == 0 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "First sequence has length 0."));
            }
            if seq_length > MAX_SEQ_LENGTH {
                return Err(io::Error::new(io::ErrorKind::InvalidData,
                    format!("Sequence length {} exceeds maximum.", seq_length)));
            }
            bitmask = vec![0u8; seq_length];
            ref_seq = seq.to_vec();
        } else if seq.len() != seq_length {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Sequence '{}' (#{}) has length {} but expected {}.",
                    record.id(), count + 1, seq.len(), seq_length)));
        }

        names.push(record.id().to_string());

        // Hot loop: OR bitmask via lookup table
        let bm = &mut bitmask;
        for i in 0..seq_length {
            unsafe { *bm.get_unchecked_mut(i) |= lookup[*seq.get_unchecked(i) as usize]; }
        }

        count += 1;
    }

    if count == 0 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Input FASTA is empty."));
    }

    eprintln!("[snpick] Pass 1: {} sequences × {} positions.", count, seq_length);
    Ok(ScanResult { bitmask, ref_seq, names, seq_length })
}

fn analyze(bitmask: &[u8], ref_seq: &[u8], lookup: &[u8; 256], include_gaps: bool)
-> (Vec<VariablePosition>, ConstantSiteCounts) {
    let mut vars = Vec::new();
    let mut cs = ConstantSiteCounts { a: 0, c: 0, g: 0, t: 0 };
    for (pos, &bits) in bitmask.iter().enumerate() {
        let ones = bits.count_ones();
        if ones > 1 {
            let rb = ref_seq[pos].to_ascii_uppercase();
            let ref_base = if lookup[rb as usize] != 0 { rb } else { bits_to_bases(bits, include_gaps)[0] };
            let alt_bases: Vec<u8> = bits_to_bases(bits, include_gaps).into_iter().filter(|&b| b != ref_base).collect();
            vars.push(VariablePosition { index: pos, ref_base, alt_bases });
        } else if ones == 1 {
            if bits & BIT_A != 0 { cs.a += 1; }
            else if bits & BIT_C != 0 { cs.c += 1; }
            else if bits & BIT_G != 0 { cs.g += 1; }
            else if bits & BIT_T != 0 { cs.t += 1; }
        }
    }
    (vars, cs)
}

// =============================================================================
// Pass 2: streaming extract — O(V) per seq
// =============================================================================

fn pass2_extract(
    path: &str, output: &str, var_positions: &[VariablePosition],
    num_samples: usize, seq_length: usize, collect_vcf: bool,
) -> io::Result<Option<Vec<u8>>> {
    let num_var = var_positions.len();
    let pos_indices: Vec<usize> = var_positions.iter().map(|v| v.index).collect();

    let file = File::open(path)?;
    let reader = BufReader::with_capacity(IO_BUF, file);
    let mut fasta = fasta::Reader::new(reader);
    let mut record = fasta::Record::new();

    let out_file = File::create(output)?;
    let mut writer = BufWriter::with_capacity(IO_BUF, out_file);

    let mut vcf_geno: Vec<u8> = if collect_vcf { vec![0u8; num_var * num_samples] } else { Vec::new() };
    let mut var_buf = vec![0u8; num_var];
    let mut si = 0usize;

    loop {
        fasta.read(&mut record)?;
        if record.is_empty() { break; }
        let seq = record.seq();

        if seq.len() != seq_length {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Pass 2: '{}' length {} ≠ {}. File modified?",
                    record.id(), seq.len(), seq_length)));
        }

        // Extract variable positions — O(V) direct index
        for (vi, &p) in pos_indices.iter().enumerate() {
            unsafe { *var_buf.get_unchecked_mut(vi) = *seq.get_unchecked(p); }
        }

        writer.write_all(b">")?;
        writer.write_all(record.id().as_bytes())?;
        writer.write_all(b"\n")?;
        writer.write_all(&var_buf)?;
        writer.write_all(b"\n")?;

        if collect_vcf {
            for (vi, &nuc) in var_buf.iter().enumerate() {
                vcf_geno[vi * num_samples + si] = nuc.to_ascii_uppercase();
            }
        }
        si += 1;
    }

    writer.flush()?;

    if si != num_samples {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            format!("Pass 2: Expected {} sequences but got {}.", num_samples, si)));
    }

    eprintln!("[snpick] Pass 2: Wrote {} sequences to {}.", si, output);
    if collect_vcf { Ok(Some(vcf_geno)) } else { Ok(None) }
}

// =============================================================================
// VCF writer
// =============================================================================

fn write_vcf(
    vcf_geno: &[u8], num_samples: usize, var_positions: &[VariablePosition],
    vcf_path: &str, names: &[String], seq_length: usize,
) -> io::Result<()> {
    let out = File::create(vcf_path)?;
    let mut w = BufWriter::with_capacity(4 * 1024 * 1024, out);
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##source=snpick v{}", env!("CARGO_PKG_VERSION"))?;
    writeln!(w, "##reference=first_sequence")?;
    writeln!(w, "##contig=<ID=1,length={}>", seq_length)?;
    writeln!(w, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">")?;
    writeln!(w, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")?;
    write!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for n in names { write!(w, "\t{}", n)?; }
    writeln!(w)?;

    let mut lut = [255u8; 256];
    for (vi, vp) in var_positions.iter().enumerate() {
        let alt: String = vp.alt_bases.iter()
            .map(|&b| if b == b'-' { "*".to_string() } else { (b as char).to_string() })
            .collect::<Vec<_>>().join(",");
        lut[vp.ref_base as usize] = 0;
        for (i, &ab) in vp.alt_bases.iter().enumerate() { lut[ab as usize] = (i + 1) as u8; }
        write!(w, "1\t{}\t.\t{}\t{}\t.\tPASS\tNS={}\tGT", vp.index + 1, vp.ref_base as char, alt, num_samples)?;
        let row = vi * num_samples;
        for si in 0..num_samples {
            let idx = lut[vcf_geno[row + si] as usize];
            if idx == 255 { write!(w, "\t.")?; } else { write!(w, "\t{}", idx)?; }
        }
        writeln!(w)?;
        lut[vp.ref_base as usize] = 255;
        for &ab in &vp.alt_bases { lut[ab as usize] = 255; }
    }
    w.flush()
}

// =============================================================================
// main
// =============================================================================

fn main() -> io::Result<()> {
    let args = Args::parse();
    let start = Instant::now();
    let lookup = build_lookup(args.include_gaps);

    check_paths_differ(&args.fasta, &args.output)?;
    if let Some(ref vp) = args.vcf_output { check_paths_differ(&args.fasta, vp)?; }

    // Pass 1: scan bitmask
    let scan = pass1_scan(&args.fasta, &lookup)?;
    let (var_positions, constant) = analyze(&scan.bitmask, &scan.ref_seq, &lookup, args.include_gaps);
    let num_var = var_positions.len();
    let num_samples = scan.names.len();
    let seq_length = scan.seq_length;

    eprintln!("[snpick] {} variable, {} constant ({}).", num_var, constant.total(), constant);
    eprintln!("[snpick] ASC fconst: {}", constant.fconst());
    if num_var == 0 { eprintln!("[snpick] No variable positions."); return Ok(()); }

    // Free scan data
    drop(scan.bitmask);
    drop(scan.ref_seq);

    // Pass 2: extract
    let vcf_geno = pass2_extract(&args.fasta, &args.output, &var_positions, num_samples, seq_length, args.vcf)?;

    if let Some(ref geno) = vcf_geno {
        let vp = args.vcf_output.unwrap_or_else(|| "output.vcf".to_string());
        write_vcf(geno, num_samples, &var_positions, &vp, &scan.names, seq_length)?;
        eprintln!("[snpick] VCF written to {}.", vp);
    }

    eprintln!("[snpick] Done in {:.2}s. {} vars from {} seqs × {} pos.",
        start.elapsed().as_secs_f64(), num_var, num_samples, seq_length);
    Ok(())
}

// =============================================================================
// Tests
// =============================================================================
#[cfg(test)]
mod tests {
    use super::*;
    fn tmp(name: &str, c: &str) -> String {
        let p = format!("/tmp/snpick_t_{}.fa", name);
        std::fs::write(&p, c).unwrap(); p
    }

    #[test] fn test_lookup() {
        let lk = build_lookup(false);
        assert_eq!(lk[b'A' as usize], BIT_A); assert_eq!(lk[b'N' as usize], 0);
        assert_eq!(build_lookup(true)[b'-' as usize], BIT_GAP);
    }
    #[test] fn test_pass1() {
        let p = tmp("p1b", ">s1\nATGC\n>s2\nATCC\n");
        let lk = build_lookup(false);
        let r = pass1_scan(&p, &lk).unwrap();
        assert_eq!(r.names.len(), 2); assert_eq!(r.seq_length, 4);
        assert_eq!(r.bitmask[2], BIT_G | BIT_C);
        std::fs::remove_file(&p).ok();
    }
    #[test] fn test_e2e() {
        let p = tmp("e2eb", ">ref\nATGCATGC\n>s1\nATGTATGC\n>s2\nACGCATGC\n");
        let o = "/tmp/snpick_t_e2eb_out.fa";
        let lk = build_lookup(false);
        let r = pass1_scan(&p, &lk).unwrap();
        let (v, _) = analyze(&r.bitmask, &r.ref_seq, &lk, false);
        pass2_extract(&p, o, &v, r.names.len(), r.seq_length, false).unwrap();
        let c = std::fs::read_to_string(o).unwrap();
        let l: Vec<&str> = c.lines().collect();
        assert_eq!(l[1], "TC"); assert_eq!(l[3], "TT"); assert_eq!(l[5], "CC");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }
    #[test] fn test_vcf() {
        let p = tmp("vcfb", ">ref\nATGC\n>s1\nCTGC\n>s2\nGTGC\n>s3\nTTGC\n");
        let fo = "/tmp/snpick_t_vcfb_out.fa"; let vo = "/tmp/snpick_t_vcfb.vcf";
        let lk = build_lookup(false);
        let r = pass1_scan(&p, &lk).unwrap();
        let (v, _) = analyze(&r.bitmask, &r.ref_seq, &lk, false);
        let g = pass2_extract(&p, fo, &v, r.names.len(), r.seq_length, true).unwrap().unwrap();
        write_vcf(&g, r.names.len(), &v, vo, &r.names, r.seq_length).unwrap();
        let c = std::fs::read_to_string(vo).unwrap();
        let dl: Vec<&str> = c.lines().filter(|l| !l.starts_with('#')).collect();
        let f: Vec<&str> = dl[0].split('\t').collect();
        assert_eq!(f[3], "A"); assert_eq!(f[4], "C,G,T"); assert_eq!(f[9], "0"); assert_eq!(f[12], "3");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(fo).ok(); std::fs::remove_file(vo).ok();
    }
    #[test] fn test_ambig() {
        let p = tmp("ambb", ">ref\nATGC\n>s1\nCTGC\n>s2\nNTGC\n");
        let fo = "/tmp/snpick_t_ambb_out.fa"; let vo = "/tmp/snpick_t_ambb.vcf";
        let lk = build_lookup(false);
        let r = pass1_scan(&p, &lk).unwrap();
        let (v, _) = analyze(&r.bitmask, &r.ref_seq, &lk, false);
        let g = pass2_extract(&p, fo, &v, r.names.len(), r.seq_length, true).unwrap().unwrap();
        write_vcf(&g, r.names.len(), &v, vo, &r.names, r.seq_length).unwrap();
        let c = std::fs::read_to_string(vo).unwrap();
        let f: Vec<&str> = c.lines().filter(|l|!l.starts_with('#')).collect::<Vec<_>>()[0].split('\t').collect();
        assert_eq!(f[11], ".");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(fo).ok(); std::fs::remove_file(vo).ok();
    }
    #[test] fn test_gaps() {
        let p = tmp("gapb", ">ref\nATGC\n>s1\nA-GC\n");
        let lk_no = build_lookup(false); let lk_yes = build_lookup(true);
        let r1 = pass1_scan(&p, &lk_no).unwrap();
        let (v1, _) = analyze(&r1.bitmask, &r1.ref_seq, &lk_no, false);
        assert!(v1.is_empty());
        let r2 = pass1_scan(&p, &lk_yes).unwrap();
        let (v2, _) = analyze(&r2.bitmask, &r2.ref_seq, &lk_yes, true);
        assert_eq!(v2.len(), 1);
        std::fs::remove_file(&p).ok();
    }
    #[test] fn test_paths() {
        let p = tmp("pdb", ">x\nA\n");
        assert!(check_paths_differ(&p, "/tmp/snpick_t_pdb2.fa").is_ok());
        assert!(check_paths_differ(&p, &p).is_err());
        std::fs::remove_file(&p).ok();
    }
    #[test] fn test_integrity() {
        let p = tmp("intb", ">s1\nATGC\n>s2\nATCC\n");
        let o = "/tmp/snpick_t_intb.fa";
        let lk = build_lookup(false);
        let r = pass1_scan(&p, &lk).unwrap();
        let (v, _) = analyze(&r.bitmask, &r.ref_seq, &lk, false);
        std::fs::write(&p, ">s1\nATGC\n").unwrap(); // truncate
        let result = pass2_extract(&p, o, &v, 2, r.seq_length, false);
        assert!(result.is_err());
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }
}
