use bio::io::fasta::{self, FastaRead};
use clap::Parser;
use memmap2::Mmap;
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
// Constants & lookup tables
// =============================================================================

const BIT_A: u8 = 0b00001;
const BIT_C: u8 = 0b00010;
const BIT_G: u8 = 0b00100;
const BIT_T: u8 = 0b01000;
const BIT_GAP: u8 = 0b10000;
const MAX_SEQ_LENGTH: usize = 10_000_000_000;
const IO_BUF: usize = 16 * 1024 * 1024;

fn build_lookup(include_gaps: bool) -> [u8; 256] {
    let mut t = [0u8; 256];
    t[b'A' as usize] = BIT_A; t[b'a' as usize] = BIT_A;
    t[b'C' as usize] = BIT_C; t[b'c' as usize] = BIT_C;
    t[b'G' as usize] = BIT_G; t[b'g' as usize] = BIT_G;
    t[b'T' as usize] = BIT_T; t[b't' as usize] = BIT_T;
    if include_gaps { t[b'-' as usize] = BIT_GAP; }
    t
}

fn build_upper() -> [u8; 256] {
    let mut t = [0u8; 256];
    for (i, v) in t.iter_mut().enumerate() { *v = i as u8; }
    for c in b'a'..=b'z' { t[c as usize] = c - 32; }
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

struct VariablePosition { index: usize, ref_base: u8, alt_bases: Vec<u8>, ns: usize }
struct ConstantSiteCounts { a: usize, c: usize, g: usize, t: usize }
struct SiteCounts { constant: ConstantSiteCounts, variable: usize, ambiguous: usize }

impl std::fmt::Display for ConstantSiteCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "A:{} C:{} G:{} T:{}", self.a, self.c, self.g, self.t)
    }
}
impl ConstantSiteCounts {
    fn total(&self) -> usize { self.a + self.c + self.g + self.t }
    fn fconst(&self) -> String { format!("{},{},{},{}", self.a, self.c, self.g, self.t) }
}

struct ScanResult {
    bitmask: Vec<u8>,
    ref_seq: Vec<u8>,
    names: Vec<String>,
    descs: Vec<String>,
    seq_length: usize,
}

// =============================================================================
// Pass 1: streaming bitmask — O(L) memory, bio reader (fastest for sequential)
// =============================================================================

fn pass1_scan(path: &str, lookup: &[u8; 256]) -> io::Result<ScanResult> {
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(IO_BUF, file);
    let mut fasta = fasta::Reader::new(reader);
    let mut record = fasta::Record::new();

    let mut bitmask: Vec<u8> = Vec::new();
    let mut ref_seq: Vec<u8> = Vec::new();
    let mut names: Vec<String> = Vec::new();
    let mut descs: Vec<String> = Vec::new();
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
        descs.push(record.desc().unwrap_or("").to_string());

        // Hot loop: OR bitmask via lookup table — process in 8-byte chunks for ILP
        let bm = &mut bitmask;
        let s = seq;
        let chunks = seq_length / 8;
        let remainder = seq_length % 8;
        for c in 0..chunks {
            let base = c * 8;
            bm[base]     |= lookup[s[base]     as usize];
            bm[base + 1] |= lookup[s[base + 1] as usize];
            bm[base + 2] |= lookup[s[base + 2] as usize];
            bm[base + 3] |= lookup[s[base + 3] as usize];
            bm[base + 4] |= lookup[s[base + 4] as usize];
            bm[base + 5] |= lookup[s[base + 5] as usize];
            bm[base + 6] |= lookup[s[base + 6] as usize];
            bm[base + 7] |= lookup[s[base + 7] as usize];
        }
        let base = chunks * 8;
        for i in 0..remainder {
            bm[base + i] |= lookup[s[base + i] as usize];
        }

        count += 1;
    }

    if count == 0 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Input FASTA is empty."));
    }

    eprintln!("[snpick] Pass 1: {} sequences × {} positions.", count, seq_length);
    Ok(ScanResult { bitmask, ref_seq, names, descs, seq_length })
}

fn analyze(bitmask: &[u8], ref_seq: &[u8], lookup: &[u8; 256], include_gaps: bool)
-> (Vec<VariablePosition>, SiteCounts) {
    let mut vars = Vec::new();
    let mut cs = ConstantSiteCounts { a: 0, c: 0, g: 0, t: 0 };
    let mut ambiguous = 0usize;
    for (pos, &bits) in bitmask.iter().enumerate() {
        let ones = bits.count_ones();
        if ones > 1 {
            let rb = ref_seq[pos].to_ascii_uppercase();
            let ref_base = if lookup[rb as usize] != 0 { rb } else { bits_to_bases(bits, include_gaps)[0] };
            let alt_bases: Vec<u8> = bits_to_bases(bits, include_gaps).into_iter()
                .filter(|&b| b != ref_base).collect();
            vars.push(VariablePosition { index: pos, ref_base, alt_bases, ns: 0 });
        } else if ones == 1 {
            if bits & BIT_A != 0 { cs.a += 1; }
            else if bits & BIT_C != 0 { cs.c += 1; }
            else if bits & BIT_G != 0 { cs.g += 1; }
            else if bits & BIT_T != 0 { cs.t += 1; }
        } else {
            ambiguous += 1;
        }
    }
    let num_variable = vars.len();
    (vars, SiteCounts { constant: cs, variable: num_variable, ambiguous })
}

// =============================================================================
// Pass 2: mmap + direct index for variable sites extraction
// =============================================================================

/// Index sequence start offsets in the mmap for pass 2 direct access.
/// Returns vec of (seq_data_offset, is_single_line) per sequence.
fn index_seq_offsets(data: &[u8], num_seqs: usize, seq_length: usize, names: &[String])
-> io::Result<Vec<usize>> {
    let mut offsets = Vec::with_capacity(num_seqs);
    let mut pos = 0;
    let len = data.len();

    for (si, expected_name) in names.iter().enumerate() {
        // Find '>'
        while pos < len && data[pos] != b'>' { pos += 1; }
        if pos >= len {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Pass 2: expected {} sequences, found only {}.", num_seqs, si)));
        }
        pos += 1; // skip '>'

        // Verify name matches
        let name = expected_name.as_bytes();
        let header_start = pos;
        while pos < len && data[pos] != b'\n' && data[pos] != b'\r'
            && data[pos] != b' ' && data[pos] != b'\t' { pos += 1; }
        let found_id = &data[header_start..pos];
        if found_id != name {
            let found = std::str::from_utf8(found_id).unwrap_or("?");
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Pass 2: sequence #{} is '{}' but pass 1 had '{}'. File modified?",
                    si + 1, found, expected_name)));
        }

        // Skip to end of header line
        while pos < len && data[pos] != b'\n' { pos += 1; }
        if pos < len { pos += 1; } // skip \n

        // This is where sequence data starts
        offsets.push(pos);

        // Skip past sequence data
        pos += seq_length; // for single-line, this is exact
        // Handle potential newlines in multi-line FASTA
        while pos < len && data[pos] != b'>' { pos += 1; }
    }

    Ok(offsets)
}

struct ExtractParams<'a> {
    names: &'a [String],
    descs: &'a [String],
    output: &'a str,
    collect_vcf: bool,
    lookup: &'a [u8; 256],
    upper: &'a [u8; 256],
}

fn pass2_extract(
    data: &[u8], offsets: &[usize], var_positions: &mut [VariablePosition],
    params: &ExtractParams<'_>,
) -> io::Result<Option<Vec<u8>>> {
    let ExtractParams { names, descs, output, collect_vcf, lookup, upper } = params;
    let collect_vcf = *collect_vcf;
    let num_var = var_positions.len();
    let num_samples = names.len();
    let pos_indices: Vec<usize> = var_positions.iter().map(|v| v.index).collect();

    let out_file = File::create(output).map_err(|e| io::Error::new(e.kind(),
        format!("Cannot create output '{}': {}", output, e)))?;
    let mut writer = BufWriter::with_capacity(IO_BUF, out_file);

    let mut vcf_geno: Vec<u8> = if collect_vcf { vec![0u8; num_var * num_samples] } else { Vec::new() };
    let mut ns_counts: Vec<usize> = if collect_vcf { vec![0usize; num_var] } else { Vec::new() };
    let mut var_buf = vec![0u8; num_var];

    for (si, &base) in offsets.iter().enumerate() {
        // Direct byte access — O(V) random reads from mmap, OS handles page cache
        for (vi, &p) in pos_indices.iter().enumerate() {
            var_buf[vi] = upper[data[base + p] as usize];
        }

        // Write output
        writer.write_all(b">")?;
        writer.write_all(names[si].as_bytes())?;
        if !descs[si].is_empty() {
            writer.write_all(b" ")?;
            writer.write_all(descs[si].as_bytes())?;
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

// =============================================================================
// VCF writer
// =============================================================================

fn write_vcf(
    vcf_geno: &[u8], num_samples: usize, var_positions: &[VariablePosition],
    vcf_path: &str, names: &[String], seq_length: usize,
) -> io::Result<()> {
    let out = File::create(vcf_path).map_err(|e| io::Error::new(e.kind(),
        format!("Cannot create VCF '{}': {}", vcf_path, e)))?;
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
        write!(w, "1\t{}\t.\t{}\t{}\t.\tPASS\tNS={}\tGT", vp.index + 1, vp.ref_base as char, alt, vp.ns)?;
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

const MAX_VCF_GENO_BYTES: usize = 4_000_000_000;

fn run() -> io::Result<()> {
    let args = Args::parse();
    let start = Instant::now();
    let lookup = build_lookup(args.include_gaps);
    let upper = build_upper();

    let do_vcf = args.vcf || args.vcf_output.is_some();

    check_paths_differ(&args.fasta, &args.output)?;
    let vcf_path = if do_vcf {
        let vp = args.vcf_output.unwrap_or_else(|| {
            let out = Path::new(&args.output);
            let stem = out.file_stem().and_then(|s| s.to_str()).unwrap_or("output");
            let parent = out.parent().unwrap_or(Path::new("."));
            parent.join(format!("{}.vcf", stem)).to_string_lossy().into_owned()
        });
        check_paths_differ(&args.fasta, &vp)?;
        check_paths_differ(&args.output, &vp)?;
        Some(vp)
    } else { None };

    // Pass 1: streaming bitmask via bio (optimal for sequential scan)
    let scan = pass1_scan(&args.fasta, &lookup)?;
    let t1 = start.elapsed().as_secs_f64();
    let (mut var_positions, site_counts) = analyze(&scan.bitmask, &scan.ref_seq, &lookup, args.include_gaps);
    let num_var = var_positions.len();
    let num_samples = scan.names.len();
    let seq_length = scan.seq_length;

    let ScanResult { names, descs, .. } = scan;

    eprintln!("[snpick] {} variable, {} constant ({}), {} ambiguous-only, {} total.",
        site_counts.variable, site_counts.constant.total(), site_counts.constant,
        site_counts.ambiguous, seq_length);
    eprintln!("[snpick] ASC fconst: {}", site_counts.constant.fconst());
    eprintln!("[snpick] Pass 1 took {:.2}s.", t1);

    if num_var == 0 {
        eprintln!("[snpick] No variable positions — writing empty output.");
        let out = File::create(&args.output)?;
        let mut w = BufWriter::new(out);
        for (i, name) in names.iter().enumerate() {
            write!(w, ">{}", name)?;
            if !descs[i].is_empty() { write!(w, " {}", descs[i])?; }
            writeln!(w)?; writeln!(w)?;
        }
        w.flush()?;
        return Ok(());
    }

    if do_vcf {
        let geno_bytes = num_var.saturating_mul(num_samples);
        if geno_bytes > MAX_VCF_GENO_BYTES {
            return Err(io::Error::new(io::ErrorKind::InvalidInput,
                format!("VCF genotype matrix would require {} GB ({} vars × {} samples). \
                    Use without --vcf or reduce input.",
                    geno_bytes / 1_000_000_000, num_var, num_samples)));
        }
    }

    // Pass 2: mmap for O(V) random access per sequence (no full re-read)
    let file = File::open(&args.fasta)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let data = &mmap[..];

    let offsets = index_seq_offsets(data, num_samples, seq_length, &names)?;
    let ep = ExtractParams { names: &names, descs: &descs, output: &args.output,
        collect_vcf: do_vcf, lookup: &lookup, upper: &upper };
    let vcf_geno = pass2_extract(data, &offsets, &mut var_positions, &ep)?;

    if let (Some(ref geno), Some(ref vp)) = (&vcf_geno, &vcf_path) {
        write_vcf(geno, num_samples, &var_positions, vp, &names, seq_length)?;
        eprintln!("[snpick] VCF written to {}.", vp);
    }

    eprintln!("[snpick] Done in {:.2}s. {} vars from {} seqs × {} pos.",
        start.elapsed().as_secs_f64(), num_var, num_samples, seq_length);
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("[snpick] Error: {}", e);
        std::process::exit(1);
    }
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
        let p = tmp("p1d", ">s1\nATGC\n>s2\nATCC\n");
        let lk = build_lookup(false);
        let r = pass1_scan(&p, &lk).unwrap();
        assert_eq!(r.names.len(), 2); assert_eq!(r.seq_length, 4);
        assert_eq!(r.bitmask[2], BIT_G | BIT_C);
        std::fs::remove_file(&p).ok();
    }
    #[test] fn test_e2e() {
        let p = tmp("e2ed", ">ref\nATGCATGC\n>s1\nATGTATGC\n>s2\nACGCATGC\n");
        let o = "/tmp/snpick_t_e2ed_out.fa";
        let lk = build_lookup(false);
        let up = build_upper();
        let r = pass1_scan(&p, &lk).unwrap();
        let (mut v, _) = analyze(&r.bitmask, &r.ref_seq, &lk, false);
        let f = File::open(&p).unwrap();
        let m = unsafe { Mmap::map(&f).unwrap() };
        let offsets = index_seq_offsets(&m, r.names.len(), r.seq_length, &r.names).unwrap();
        pass2_extract(&m, &offsets, &mut v, &ExtractParams { names: &r.names, descs: &r.descs, output: o, collect_vcf: false, lookup: &lk, upper: &up }).unwrap();
        let c = std::fs::read_to_string(o).unwrap();
        let l: Vec<&str> = c.lines().collect();
        assert_eq!(l[1], "TC"); assert_eq!(l[3], "TT"); assert_eq!(l[5], "CC");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }
    #[test] fn test_vcf() {
        let p = tmp("vcfd", ">ref\nATGC\n>s1\nCTGC\n>s2\nGTGC\n>s3\nTTGC\n");
        let fo = "/tmp/snpick_t_vcfd_out.fa"; let vo = "/tmp/snpick_t_vcfd.vcf";
        let lk = build_lookup(false);
        let up = build_upper();
        let r = pass1_scan(&p, &lk).unwrap();
        let (mut v, _) = analyze(&r.bitmask, &r.ref_seq, &lk, false);
        let f = File::open(&p).unwrap();
        let m = unsafe { Mmap::map(&f).unwrap() };
        let offsets = index_seq_offsets(&m, r.names.len(), r.seq_length, &r.names).unwrap();
        let g = pass2_extract(&m, &offsets, &mut v, &ExtractParams { names: &r.names, descs: &r.descs, output: fo, collect_vcf: true, lookup: &lk, upper: &up }).unwrap().unwrap();
        write_vcf(&g, r.names.len(), &v, vo, &r.names, r.seq_length).unwrap();
        let c = std::fs::read_to_string(vo).unwrap();
        let dl: Vec<&str> = c.lines().filter(|l| !l.starts_with('#')).collect();
        let f: Vec<&str> = dl[0].split('\t').collect();
        assert_eq!(f[3], "A"); assert_eq!(f[4], "C,G,T"); assert_eq!(f[9], "0");
        assert_eq!(f[7], "NS=4"); assert_eq!(f[12], "3");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(fo).ok(); std::fs::remove_file(vo).ok();
    }
    #[test] fn test_ambig() {
        let p = tmp("ambd", ">ref\nATGC\n>s1\nCTGC\n>s2\nNTGC\n");
        let fo = "/tmp/snpick_t_ambd_out.fa"; let vo = "/tmp/snpick_t_ambd.vcf";
        let lk = build_lookup(false);
        let up = build_upper();
        let r = pass1_scan(&p, &lk).unwrap();
        let (mut v, _) = analyze(&r.bitmask, &r.ref_seq, &lk, false);
        let f = File::open(&p).unwrap();
        let m = unsafe { Mmap::map(&f).unwrap() };
        let offsets = index_seq_offsets(&m, r.names.len(), r.seq_length, &r.names).unwrap();
        let g = pass2_extract(&m, &offsets, &mut v, &ExtractParams { names: &r.names, descs: &r.descs, output: fo, collect_vcf: true, lookup: &lk, upper: &up }).unwrap().unwrap();
        write_vcf(&g, r.names.len(), &v, vo, &r.names, r.seq_length).unwrap();
        let c = std::fs::read_to_string(vo).unwrap();
        let dl: Vec<&str> = c.lines().filter(|l|!l.starts_with('#')).collect();
        let f: Vec<&str> = dl[0].split('\t').collect();
        assert_eq!(f[7], "NS=2"); assert_eq!(f[11], ".");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(fo).ok(); std::fs::remove_file(vo).ok();
    }
    #[test] fn test_desc_preserved() {
        let p = tmp("descd", ">s1 some description\nATGC\n>s2 another desc\nATCC\n");
        let o = "/tmp/snpick_t_descd_out.fa";
        let lk = build_lookup(false);
        let up = build_upper();
        let r = pass1_scan(&p, &lk).unwrap();
        let (mut v, _) = analyze(&r.bitmask, &r.ref_seq, &lk, false);
        let f = File::open(&p).unwrap();
        let m = unsafe { Mmap::map(&f).unwrap() };
        let offsets = index_seq_offsets(&m, r.names.len(), r.seq_length, &r.names).unwrap();
        pass2_extract(&m, &offsets, &mut v, &ExtractParams { names: &r.names, descs: &r.descs, output: o, collect_vcf: false, lookup: &lk, upper: &up }).unwrap();
        let c = std::fs::read_to_string(o).unwrap();
        assert!(c.contains(">s1 some description"));
        assert!(c.contains(">s2 another desc"));
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }
    #[test] fn test_gaps() {
        let p = tmp("gapd", ">ref\nATGC\n>s1\nA-GC\n");
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
        let p = tmp("pdd", ">x\nA\n");
        assert!(check_paths_differ(&p, "/tmp/snpick_t_pdd2.fa").is_ok());
        assert!(check_paths_differ(&p, &p).is_err());
        std::fs::remove_file(&p).ok();
    }
    #[test] fn test_integrity() {
        let p = tmp("intd", ">s1\nATGC\n>s2\nATCC\n");
        let lk = build_lookup(false);
        let r = pass1_scan(&p, &lk).unwrap();
        let (mut v, _) = analyze(&r.bitmask, &r.ref_seq, &lk, false);
        std::fs::write(&p, ">s1\nATGC\n").unwrap();
        let f = File::open(&p).unwrap();
        let m = unsafe { Mmap::map(&f).unwrap() };
        let result = index_seq_offsets(&m, r.names.len(), r.seq_length, &r.names);
        assert!(result.is_err());
        std::fs::remove_file(&p).ok();
    }
    #[test] fn test_name_mismatch() {
        let p = tmp("nmd", ">s1\nATGC\n>s2\nATCC\n");
        let lk = build_lookup(false);
        let r = pass1_scan(&p, &lk).unwrap();
        std::fs::write(&p, ">s2\nATGC\n>s1\nATCC\n").unwrap();
        let f = File::open(&p).unwrap();
        let m = unsafe { Mmap::map(&f).unwrap() };
        let result = index_seq_offsets(&m, r.names.len(), r.seq_length, &r.names);
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("s2") && err.contains("s1"));
        std::fs::remove_file(&p).ok();
    }
    #[test] fn test_empty() {
        let p = tmp("empd", "");
        let lk = build_lookup(false);
        assert!(pass1_scan(&p, &lk).is_err());
        std::fs::remove_file(&p).ok();
    }
}
