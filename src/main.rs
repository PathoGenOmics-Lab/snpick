use clap::Parser;
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
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
// Streaming FASTA reader — zero-alloc per-record, reads into reusable buffers
// =============================================================================




// =============================================================================
// Streaming FASTA reader v2 — chunk-based, zero-alloc per base
// =============================================================================

struct FastaReader2 {
    buf: Vec<u8>,
    pos: usize,
    len: usize,
    file: File,
    eof: bool,
    // Reusable buffers
    name: Vec<u8>,
    seq: Vec<u8>,
}

const READ_BUF_SIZE: usize = 16 * 1024 * 1024; // 16 MB

impl FastaReader2 {
    fn open(path: &str) -> io::Result<Self> {
        let file = File::open(path)?;
        let mut r = Self {
            buf: vec![0u8; READ_BUF_SIZE],
            pos: 0, len: 0, file, eof: false,
            name: Vec::with_capacity(256),
            seq: Vec::with_capacity(8 * 1024 * 1024),
        };
        r.fill()?;
        Ok(r)
    }

    fn fill(&mut self) -> io::Result<()> {
        if self.eof { return Ok(()); }
        // Move remaining to start
        let remaining = self.len - self.pos;
        if remaining > 0 {
            self.buf.copy_within(self.pos..self.len, 0);
        }
        self.len = remaining;
        self.pos = 0;
        // Read more
        let n = self.file.read(&mut self.buf[self.len..])?;
        if n == 0 { self.eof = true; }
        self.len += n;
        Ok(())
    }

    #[inline(always)]
    fn peek(&mut self) -> io::Result<Option<u8>> {
        if self.pos >= self.len {
            self.fill()?;
            if self.pos >= self.len { return Ok(None); }
        }
        Ok(Some(self.buf[self.pos]))
    }

    #[inline(always)]
    fn advance(&mut self) { self.pos += 1; }

    /// Read next FASTA record into self.name and self.seq.
    /// Returns false at EOF.
    fn next(&mut self) -> io::Result<bool> {
        self.name.clear();
        self.seq.clear();

        // Find '>'
        loop {
            match self.peek()? {
                None => return Ok(false),
                Some(b'>') => { self.advance(); break; }
                _ => { self.advance(); }
            }
        }

        // Read name until newline
        loop {
            match self.peek()? {
                None => return Ok(!self.name.is_empty()),
                Some(b'\n') | Some(b'\r') => { self.advance(); break; }
                Some(b) => { self.name.push(b); self.advance(); }
            }
        }
        // Skip extra \n or \r
        while let Some(b'\n') | Some(b'\r') = self.peek()? {
            self.advance();
        }

        // Read sequence bases — stop at '>' or EOF
        loop {
            if self.pos >= self.len {
                self.fill()?;
                if self.pos >= self.len { break; }
            }
            // Process chunk: scan for '>' in remaining buffer
            let start = self.pos;
            let mut end = start;
            while end < self.len {
                let b = self.buf[end];
                if b == b'>' { break; }
                end += 1;
            }
            // Append non-whitespace bytes from buf[start..end]
            for &b in &self.buf[start..end] {
                if b != b'\n' && b != b'\r' {
                    self.seq.push(b);
                }
            }
            self.pos = end;
            if end < self.len { break; } // hit '>'
        }

        Ok(true)
    }
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
// Pass 1: streaming bitmask scan — O(L) memory
// =============================================================================

type ScanResult = (Vec<u8>, Vec<u8>, Vec<String>, usize); // (bitmask, ref_seq, names, seq_length)

fn pass1_scan(path: &str, lookup: &[u8; 256]) -> io::Result<ScanResult> {
    let mut reader = FastaReader2::open(path)?;
    let mut bitmask: Vec<u8> = Vec::new();
    let mut ref_seq: Vec<u8> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut seq_length = 0usize;
    let mut count = 0usize;

    while reader.next()? {
        if count == 0 {
            seq_length = reader.seq.len();
            if seq_length == 0 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "First sequence has length 0."));
            }
            if seq_length > MAX_SEQ_LENGTH {
                return Err(io::Error::new(io::ErrorKind::InvalidData,
                    format!("Sequence length {} exceeds maximum.", seq_length)));
            }
            bitmask = vec![0u8; seq_length];
            ref_seq = reader.seq.clone();
        } else if reader.seq.len() != seq_length {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Sequence '{}' (#{}) has length {} but expected {}.",
                    String::from_utf8_lossy(&reader.name), count + 1, reader.seq.len(), seq_length)));
        }

        names.push(String::from_utf8_lossy(&reader.name).into_owned());

        // Hot loop: OR bitmask
        let seq = &reader.seq;
        for i in 0..seq_length {
            unsafe { *bitmask.get_unchecked_mut(i) |= lookup[*seq.get_unchecked(i) as usize]; }
        }

        count += 1;
    }

    if count == 0 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Input FASTA is empty."));
    }

    eprintln!("[snpick] Pass 1: {} sequences × {} positions.", count, seq_length);
    Ok((bitmask, ref_seq, names, seq_length))
}

fn analyze(bitmask: &[u8], ref_seq: &[u8], lookup: &[u8; 256], include_gaps: bool
) -> (Vec<VariablePosition>, ConstantSiteCounts) {
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
// Pass 2: streaming extract — O(V) per sequence, write FASTA + VCF geno
// =============================================================================

fn pass2_extract(
    path: &str, output: &str, var_positions: &[VariablePosition],
    num_samples: usize, seq_length: usize, collect_vcf: bool,
) -> io::Result<Option<Vec<u8>>> {
    let num_var = var_positions.len();
    let pos_indices: Vec<usize> = var_positions.iter().map(|v| v.index).collect();

    let mut reader = FastaReader2::open(path)?;
    let out = File::create(output)?;
    let mut writer = BufWriter::with_capacity(16 * 1024 * 1024, out);

    let mut vcf_geno: Vec<u8> = if collect_vcf { vec![0u8; num_var * num_samples] } else { Vec::new() };
    let mut var_buf = vec![0u8; num_var];
    let mut si = 0usize;

    while reader.next()? {
        if reader.seq.len() != seq_length {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Pass 2: '{}' has length {} but expected {}. File modified?",
                    String::from_utf8_lossy(&reader.name), reader.seq.len(), seq_length)));
        }

        // Extract variable positions — O(V)
        let seq = &reader.seq;
        for (vi, &p) in pos_indices.iter().enumerate() {
            var_buf[vi] = seq[p];
        }

        writer.write_all(b">")?;
        writer.write_all(&reader.name)?;
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
// Write VCF
// =============================================================================

fn write_vcf(
    vcf_geno: &[u8], num_samples: usize, var_positions: &[VariablePosition],
    vcf_path: &str, sample_names: &[String], seq_length: usize,
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
    for n in sample_names { write!(w, "\t{}", n)?; }
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

    let (bitmask, ref_seq, names, seq_length) = pass1_scan(&args.fasta, &lookup)?;
    let num_samples = names.len();

    let (var_positions, constant) = analyze(&bitmask, &ref_seq, &lookup, args.include_gaps);
    drop(bitmask); drop(ref_seq);

    let num_var = var_positions.len();
    eprintln!("[snpick] {} variable, {} constant ({}).", num_var, constant.total(), constant);
    eprintln!("[snpick] ASC fconst: {}", constant.fconst());
    if num_var == 0 { eprintln!("[snpick] No variable positions."); return Ok(()); }

    let vcf_geno = pass2_extract(&args.fasta, &args.output, &var_positions, num_samples, seq_length, args.vcf)?;

    if let Some(ref geno) = vcf_geno {
        let vp = args.vcf_output.unwrap_or_else(|| "output.vcf".to_string());
        write_vcf(geno, num_samples, &var_positions, &vp, &names, seq_length)?;
        eprintln!("[snpick] VCF written to {}.", vp);
    }

    eprintln!("[snpick] Done in {:.2}s. {} vars from {} seqs × {} pos.",
        start.elapsed().as_secs_f64(), num_var, num_samples, seq_length);
    Ok(())
}

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
    }
    #[test] fn test_pass1() {
        let p = tmp("p1", ">s1\nATGC\n>s2\nATCC\n");
        let lk = build_lookup(false);
        let (bm, rs, names, sl) = pass1_scan(&p, &lk).unwrap();
        assert_eq!(names.len(), 2); assert_eq!(sl, 4); assert_eq!(rs, b"ATGC");
        assert_eq!(bm[2], BIT_G | BIT_C);
        std::fs::remove_file(&p).ok();
    }
    #[test] fn test_e2e() {
        let p = tmp("e2e6", ">ref\nATGCATGC\n>s1\nATGTATGC\n>s2\nACGCATGC\n");
        let o = "/tmp/snpick_t_e2e6.fa";
        let lk = build_lookup(false);
        let (bm, rs, names, sl) = pass1_scan(&p, &lk).unwrap();
        let (v, _) = analyze(&bm, &rs, &lk, false);
        pass2_extract(&p, o, &v, names.len(), sl, false).unwrap();
        let c = std::fs::read_to_string(o).unwrap();
        let l: Vec<&str> = c.lines().collect();
        assert_eq!(l[1], "TC"); assert_eq!(l[3], "TT"); assert_eq!(l[5], "CC");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }
    #[test] fn test_vcf() {
        let p = tmp("vcf6", ">ref\nATGC\n>s1\nCTGC\n>s2\nGTGC\n>s3\nTTGC\n");
        let fo = "/tmp/snpick_t_vcf6.fa"; let vo = "/tmp/snpick_t_vcf6.vcf";
        let lk = build_lookup(false);
        let (bm, rs, names, sl) = pass1_scan(&p, &lk).unwrap();
        let (v, _) = analyze(&bm, &rs, &lk, false);
        let g = pass2_extract(&p, fo, &v, names.len(), sl, true).unwrap().unwrap();
        write_vcf(&g, names.len(), &v, vo, &names, sl).unwrap();
        let c = std::fs::read_to_string(vo).unwrap();
        let dl: Vec<&str> = c.lines().filter(|l| !l.starts_with('#')).collect();
        let f: Vec<&str> = dl[0].split('\t').collect();
        assert_eq!(f[3], "A"); assert_eq!(f[4], "C,G,T"); assert_eq!(f[9], "0"); assert_eq!(f[12], "3");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(fo).ok(); std::fs::remove_file(vo).ok();
    }
    #[test] fn test_ambig() {
        let p = tmp("amb6", ">ref\nATGC\n>s1\nCTGC\n>s2\nNTGC\n");
        let fo = "/tmp/snpick_t_amb6.fa"; let vo = "/tmp/snpick_t_amb6.vcf";
        let lk = build_lookup(false);
        let (bm, rs, names, sl) = pass1_scan(&p, &lk).unwrap();
        let (v, _) = analyze(&bm, &rs, &lk, false);
        let g = pass2_extract(&p, fo, &v, names.len(), sl, true).unwrap().unwrap();
        write_vcf(&g, names.len(), &v, vo, &names, sl).unwrap();
        let c = std::fs::read_to_string(vo).unwrap();
        let f: Vec<&str> = c.lines().filter(|l|!l.starts_with('#')).collect::<Vec<_>>()[0].split('\t').collect();
        assert_eq!(f[11], ".");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(fo).ok(); std::fs::remove_file(vo).ok();
    }
    #[test] fn test_gaps() {
        let p = tmp("gap6", ">ref\nATGC\n>s1\nA-GC\n");
        let lk_no = build_lookup(false);
        let (bm, rs, _, _) = pass1_scan(&p, &lk_no).unwrap();
        let (v1, _) = analyze(&bm, &rs, &lk_no, false);
        assert!(v1.is_empty());
        let lk_yes = build_lookup(true);
        let (bm2, rs2, _, _) = pass1_scan(&p, &lk_yes).unwrap();
        let (v2, _) = analyze(&bm2, &rs2, &lk_yes, true);
        assert_eq!(v2.len(), 1);
        std::fs::remove_file(&p).ok();
    }
    #[test] fn test_paths() {
        let p = tmp("pd6", ">x\nA\n");
        assert!(check_paths_differ(&p, "/tmp/snpick_t_pd62.fa").is_ok());
        assert!(check_paths_differ(&p, &p).is_err());
        std::fs::remove_file(&p).ok();
    }
}
