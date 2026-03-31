use clap::Parser;
use memmap2::Mmap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
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

fn resolve_path(p: &str) -> io::Result<std::path::PathBuf> {
    let path = Path::new(p);
    if path.exists() {
        return std::fs::canonicalize(path);
    }
    let parent = path.parent().unwrap_or(Path::new("."));
    let parent_abs = std::fs::canonicalize(parent).map_err(|e| {
        io::Error::new(e.kind(), format!("Cannot resolve parent of '{}': {}", p, e))
    })?;
    Ok(parent_abs.join(path.file_name().unwrap_or_default()))
}

fn check_paths_differ(a: &str, b: &str) -> io::Result<()> {
    let pa = resolve_path(a)?;
    let pb = resolve_path(b)?;
    if pa == pb {
        return Err(io::Error::new(io::ErrorKind::InvalidInput,
            format!("Paths resolve to same file: {}", pa.display())));
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

// =============================================================================
// FASTA layout info
// =============================================================================

/// Whether all sequences are single-line (no embedded newlines).
/// When true, `data[seq_offset + pos]` gives the base at `pos` directly.
/// When false, we must build a newline-skip map per record for pass 2.
#[derive(Clone, Copy)]
struct SeqLayout {
    single_line: bool,
}

// =============================================================================
// FASTA record: zero-copy view into mmap
// =============================================================================

struct FastaRecord<'a> {
    id: &'a [u8],
    desc: &'a [u8],
    seq_offset: usize,
}

/// Index FASTA records from mmap. Zero-copy: stores &[u8] slices.
/// Returns (records, seq_length, layout).
fn index_fasta(data: &[u8]) -> io::Result<(Vec<FastaRecord<'_>>, usize, SeqLayout)> {
    if data.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Input FASTA is empty (0 bytes)."));
    }

    let mut records = Vec::new();
    let mut pos = 0;
    let len = data.len();
    let mut seq_length = 0usize;
    let mut is_single_line = true;

    while pos < len {
        // Skip whitespace before '>'
        while pos < len && (data[pos] == b'\n' || data[pos] == b'\r') { pos += 1; }
        if pos >= len { break; }

        if data[pos] != b'>' {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Expected '>' at byte {}, got '{}'.", pos, data[pos] as char)));
        }
        pos += 1;

        // Parse header line
        let header_start = pos;
        while pos < len && data[pos] != b'\n' && data[pos] != b'\r' { pos += 1; }
        let header = &data[header_start..pos];
        if pos < len && data[pos] == b'\r' { pos += 1; }
        if pos < len && data[pos] == b'\n' { pos += 1; }

        // Split id and desc
        let (id, desc) = if let Some(sp) = header.iter().position(|&b| b == b' ' || b == b'\t') {
            (&header[..sp], &header[sp + 1..])
        } else {
            (header, &data[0..0])
        };

        let seq_offset = pos;

        // Scan sequence: count bases
        let mut seq_len = 0usize;
        let mut line_count = 0usize;
        while pos < len && data[pos] != b'>' {
            let line_start = pos;
            while pos < len && data[pos] != b'\n' && data[pos] != b'\r' { pos += 1; }
            let line_len = pos - line_start;
            if line_len > 0 {
                seq_len += line_len;
                line_count += 1;
            }
            if pos < len && data[pos] == b'\r' { pos += 1; }
            if pos < len && data[pos] == b'\n' { pos += 1; }
        }

        if line_count > 1 { is_single_line = false; }

        if records.is_empty() {
            if seq_len == 0 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "First sequence has length 0."));
            }
            if seq_len > MAX_SEQ_LENGTH {
                return Err(io::Error::new(io::ErrorKind::InvalidData,
                    format!("Sequence length {} exceeds maximum.", seq_len)));
            }
            seq_length = seq_len;
        } else if seq_len != seq_length {
            let id_str = std::str::from_utf8(id).unwrap_or("?");
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Sequence '{}' (#{}) has length {} but expected {}.",
                    id_str, records.len() + 1, seq_len, seq_length)));
        }

        records.push(FastaRecord { id, desc, seq_offset });
    }

    if records.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "Input FASTA is empty."));
    }

    Ok((records, seq_length, SeqLayout { single_line: is_single_line }))
}

// =============================================================================
// Pass 1: bitmask scan — zero-copy over mmap, no allocations
// =============================================================================

fn pass1_scan(data: &[u8], records: &[FastaRecord], seq_length: usize,
              layout: SeqLayout, lookup: &[u8; 256]) -> Vec<u8> {
    let mut bitmask = vec![0u8; seq_length];

    if layout.single_line {
        // Fast path: L1-friendly chunks, direct slice access
        const CHUNK: usize = 32 * 1024;
        let mut chunk_start = 0;
        while chunk_start < seq_length {
            let chunk_end = (chunk_start + CHUNK).min(seq_length);
            let bm_chunk = &mut bitmask[chunk_start..chunk_end];
            for rec in records {
                let seq_chunk = &data[rec.seq_offset + chunk_start..rec.seq_offset + chunk_end];
                for (bm_byte, &seq_byte) in bm_chunk.iter_mut().zip(seq_chunk.iter()) {
                    *bm_byte |= lookup[seq_byte as usize];
                }
            }
            chunk_start = chunk_end;
        }
    } else {
        // Multi-line: byte-by-byte scan skipping newlines
        for rec in records {
            let mut pos = rec.seq_offset;
            let end = data.len();
            let mut i = 0;
            while i < seq_length && pos < end {
                let b = data[pos];
                pos += 1;
                if b == b'\n' || b == b'\r' { continue; }
                bitmask[i] |= lookup[b as usize];
                i += 1;
            }
        }
    }

    bitmask
}

/// Extract ref_seq from first record
fn get_ref_seq(data: &[u8], rec: &FastaRecord, seq_length: usize, layout: SeqLayout) -> Vec<u8> {
    if layout.single_line {
        data[rec.seq_offset..rec.seq_offset + seq_length].to_vec()
    } else {
        let mut seq = Vec::with_capacity(seq_length);
        let mut pos = rec.seq_offset;
        let end = data.len();
        while seq.len() < seq_length && pos < end {
            let b = data[pos];
            pos += 1;
            if b != b'\n' && b != b'\r' { seq.push(b); }
        }
        seq
    }
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
// Pass 2: extract variable sites — zero-copy, layout-aware
// =============================================================================

struct ExtractParams<'a> {
    records: &'a [FastaRecord<'a>],
    output: &'a str,
    collect_vcf: bool,
    lookup: &'a [u8; 256],
    upper: &'a [u8; 256],
    layout: SeqLayout,
}

fn pass2_extract(
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
            // Fast path: direct byte access
            let base = rec.seq_offset;
            for (vi, &p) in pos_indices.iter().enumerate() {
                var_buf[vi] = upper[data[base + p] as usize];
            }
        } else {
            // Multi-line: scan to build position→byte mapping, then extract
            // We only need the variable positions, so scan linearly and pick them
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

        // Write output
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

// =============================================================================
// VCF writer
// =============================================================================

fn write_vcf(
    vcf_geno: &[u8], num_samples: usize, var_positions: &[VariablePosition],
    vcf_path: &str, records: &[FastaRecord], seq_length: usize,
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
    for rec in records {
        write!(w, "\t")?;
        w.write_all(rec.id)?;
    }
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

    // Memory-map input file (single mapping for both passes)
    let file = File::open(&args.fasta).map_err(|e| io::Error::new(e.kind(),
        format!("Cannot open '{}': {}", args.fasta, e)))?;
    let file_len = file.metadata()?.len();
    if file_len == 0 {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            format!("Input file '{}' is empty (0 bytes).", args.fasta)));
    }
    let mmap = unsafe { Mmap::map(&file).map_err(|e| io::Error::new(e.kind(),
        format!("Cannot memory-map '{}': {}", args.fasta, e)))? };
    let data = &mmap[..];

    // Index all records (zero-copy slices into mmap)
    let (records, seq_length, layout) = index_fasta(data)?;
    let num_samples = records.len();

    eprintln!("[snpick] Mapped {} bytes. {} sequences × {} positions.{}",
        data.len(), num_samples, seq_length,
        if layout.single_line { "" } else { " (multi-line FASTA)" });

    // Pass 1: bitmask scan (zero-copy, direct over mmap)
    let bitmask = pass1_scan(data, &records, seq_length, layout, &lookup);
    let ref_seq = get_ref_seq(data, &records[0], seq_length, layout);
    let t1 = start.elapsed().as_secs_f64();

    let (mut var_positions, site_counts) = analyze(&bitmask, &ref_seq, &lookup, args.include_gaps);
    let num_var = var_positions.len();

    drop(bitmask);
    drop(ref_seq);

    eprintln!("[snpick] {} variable, {} constant ({}), {} ambiguous-only, {} total.",
        site_counts.variable, site_counts.constant.total(), site_counts.constant,
        site_counts.ambiguous, seq_length);
    eprintln!("[snpick] ASC fconst: {}", site_counts.constant.fconst());
    eprintln!("[snpick] Pass 1 took {:.2}s.", t1);

    if num_var == 0 {
        eprintln!("[snpick] No variable positions — writing empty output.");
        let out = File::create(&args.output)?;
        let mut w = BufWriter::new(out);
        for rec in &records {
            w.write_all(b">")?;
            w.write_all(rec.id)?;
            if !rec.desc.is_empty() { w.write_all(b" ")?; w.write_all(rec.desc)?; }
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

    // Pass 2: extract variable sites (same mmap, layout-aware byte access)
    let ep = ExtractParams {
        records: &records, output: &args.output,
        collect_vcf: do_vcf, lookup: &lookup, upper: &upper, layout,
    };
    let vcf_geno = pass2_extract(data, &mut var_positions, &ep)?;

    if let (Some(ref geno), Some(ref vp)) = (&vcf_geno, &vcf_path) {
        write_vcf(geno, num_samples, &var_positions, vp, &records, seq_length)?;
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
    fn setup(path: &str) -> Mmap {
        let f = File::open(path).unwrap();
        unsafe { Mmap::map(&f).unwrap() }
    }

    #[test] fn test_lookup() {
        let lk = build_lookup(false);
        assert_eq!(lk[b'A' as usize], BIT_A); assert_eq!(lk[b'N' as usize], 0);
        assert_eq!(build_lookup(true)[b'-' as usize], BIT_GAP);
    }
    #[test] fn test_index() {
        let p = tmp("idx3", ">s1\nATGC\n>s2\nATCC\n");
        let m = setup(&p);
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        assert_eq!(recs.len(), 2); assert_eq!(sl, 4);
        assert!(layout.single_line);
        assert_eq!(recs[0].id, b"s1"); assert_eq!(recs[1].id, b"s2");
        std::fs::remove_file(&p).ok();
    }
    #[test] fn test_pass1() {
        let p = tmp("p1f", ">s1\nATGC\n>s2\nATCC\n");
        let m = setup(&p);
        let lk = build_lookup(false);
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        assert_eq!(bm[2], BIT_G | BIT_C);
        std::fs::remove_file(&p).ok();
    }
    #[test] fn test_e2e() {
        let p = tmp("e2ef", ">ref\nATGCATGC\n>s1\nATGTATGC\n>s2\nACGCATGC\n");
        let o = "/tmp/snpick_t_e2ef_out.fa";
        let m = setup(&p);
        let lk = build_lookup(false);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, false);
        let ep = ExtractParams { records: &recs, output: o, collect_vcf: false, lookup: &lk, upper: &up, layout };
        pass2_extract(&m, &mut v, &ep).unwrap();
        let c = std::fs::read_to_string(o).unwrap();
        let l: Vec<&str> = c.lines().collect();
        assert_eq!(l[1], "TC"); assert_eq!(l[3], "TT"); assert_eq!(l[5], "CC");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }
    #[test] fn test_vcf() {
        let p = tmp("vcff", ">ref\nATGC\n>s1\nCTGC\n>s2\nGTGC\n>s3\nTTGC\n");
        let fo = "/tmp/snpick_t_vcff_out.fa"; let vo = "/tmp/snpick_t_vcff.vcf";
        let m = setup(&p);
        let lk = build_lookup(false);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, false);
        let ep = ExtractParams { records: &recs, output: fo, collect_vcf: true, lookup: &lk, upper: &up, layout };
        let g = pass2_extract(&m, &mut v, &ep).unwrap().unwrap();
        write_vcf(&g, recs.len(), &v, vo, &recs, sl).unwrap();
        let c = std::fs::read_to_string(vo).unwrap();
        let dl: Vec<&str> = c.lines().filter(|l| !l.starts_with('#')).collect();
        let f: Vec<&str> = dl[0].split('\t').collect();
        assert_eq!(f[3], "A"); assert_eq!(f[4], "C,G,T"); assert_eq!(f[9], "0");
        assert_eq!(f[7], "NS=4"); assert_eq!(f[12], "3");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(fo).ok(); std::fs::remove_file(vo).ok();
    }
    #[test] fn test_ambig() {
        let p = tmp("ambf", ">ref\nATGC\n>s1\nCTGC\n>s2\nNTGC\n");
        let fo = "/tmp/snpick_t_ambf_out.fa"; let vo = "/tmp/snpick_t_ambf.vcf";
        let m = setup(&p);
        let lk = build_lookup(false);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, false);
        let ep = ExtractParams { records: &recs, output: fo, collect_vcf: true, lookup: &lk, upper: &up, layout };
        let g = pass2_extract(&m, &mut v, &ep).unwrap().unwrap();
        write_vcf(&g, recs.len(), &v, vo, &recs, sl).unwrap();
        let c = std::fs::read_to_string(vo).unwrap();
        let dl: Vec<&str> = c.lines().filter(|l| !l.starts_with('#')).collect();
        let f: Vec<&str> = dl[0].split('\t').collect();
        assert_eq!(f[7], "NS=2"); assert_eq!(f[11], ".");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(fo).ok(); std::fs::remove_file(vo).ok();
    }
    #[test] fn test_desc_preserved() {
        let p = tmp("descf", ">s1 some description\nATGC\n>s2 another desc\nATCC\n");
        let o = "/tmp/snpick_t_descf_out.fa";
        let m = setup(&p);
        let lk = build_lookup(false);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, false);
        let ep = ExtractParams { records: &recs, output: o, collect_vcf: false, lookup: &lk, upper: &up, layout };
        pass2_extract(&m, &mut v, &ep).unwrap();
        let c = std::fs::read_to_string(o).unwrap();
        assert!(c.contains(">s1 some description"));
        assert!(c.contains(">s2 another desc"));
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }
    #[test] fn test_gaps() {
        let p = tmp("gapf", ">ref\nATGC\n>s1\nA-GC\n");
        let m = setup(&p);
        let lk_no = build_lookup(false); let lk_yes = build_lookup(true);
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        let bm1 = pass1_scan(&m, &recs, sl, layout, &lk_no);
        let rs1 = get_ref_seq(&m, &recs[0], sl, layout);
        let (v1, _) = analyze(&bm1, &rs1, &lk_no, false);
        assert!(v1.is_empty());
        let bm2 = pass1_scan(&m, &recs, sl, layout, &lk_yes);
        let (v2, _) = analyze(&bm2, &rs1, &lk_yes, true);
        assert_eq!(v2.len(), 1);
        std::fs::remove_file(&p).ok();
    }
    #[test] fn test_paths() {
        let p = tmp("pdf", ">x\nA\n");
        assert!(check_paths_differ(&p, "/tmp/snpick_t_pdf2.fa").is_ok());
        assert!(check_paths_differ(&p, &p).is_err());
        std::fs::remove_file(&p).ok();
    }
    #[test] fn test_multiline_pass1() {
        // Multi-line FASTA: pass 1 correctness
        let p = tmp("mlf1", ">s1\nAT\nGC\n>s2\nAT\nCC\n");
        let m = setup(&p);
        let lk = build_lookup(false);
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        assert_eq!(sl, 4);
        assert!(!layout.single_line);
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (v, _) = analyze(&bm, &rs, &lk, false);
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].index, 2); // G vs C at position 2
        std::fs::remove_file(&p).ok();
    }
    #[test] fn test_multiline_pass2() {
        // Multi-line FASTA: pass 2 correctness (B1 regression test)
        let p = tmp("mlf2", ">s1\nAT\nGC\n>s2\nAT\nCC\n");
        let o = "/tmp/snpick_t_mlf2_out.fa";
        let m = setup(&p);
        let lk = build_lookup(false);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        assert!(!layout.single_line);
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, false);
        let ep = ExtractParams { records: &recs, output: o, collect_vcf: false, lookup: &lk, upper: &up, layout };
        pass2_extract(&m, &mut v, &ep).unwrap();
        let c = std::fs::read_to_string(o).unwrap();
        let l: Vec<&str> = c.lines().collect();
        assert_eq!(l[1], "G"); // s1 position 2 = G
        assert_eq!(l[3], "C"); // s2 position 2 = C
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }
    #[test] fn test_multiline_vcf() {
        // Multi-line FASTA with VCF output (full pipeline)
        let p = tmp("mlf3", ">ref\nAT\nGC\n>s1\nAT\nCC\n>s2\nCT\nGC\n");
        let fo = "/tmp/snpick_t_mlf3_out.fa"; let vo = "/tmp/snpick_t_mlf3.vcf";
        let m = setup(&p);
        let lk = build_lookup(false);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        assert!(!layout.single_line);
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, false);
        let ep = ExtractParams { records: &recs, output: fo, collect_vcf: true, lookup: &lk, upper: &up, layout };
        let g = pass2_extract(&m, &mut v, &ep).unwrap().unwrap();
        write_vcf(&g, recs.len(), &v, vo, &recs, sl).unwrap();
        // Check FASTA output
        let c = std::fs::read_to_string(fo).unwrap();
        let l: Vec<&str> = c.lines().collect();
        assert_eq!(l[1], "AG");  // ref: pos0=A, pos2=G
        assert_eq!(l[3], "AC");  // s1: pos0=A, pos2=C
        assert_eq!(l[5], "CG");  // s2: pos0=C, pos2=G
        // Check VCF
        let vc = std::fs::read_to_string(vo).unwrap();
        let dl: Vec<&str> = vc.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(dl.len(), 2); // 2 variable positions
        std::fs::remove_file(&p).ok(); std::fs::remove_file(fo).ok(); std::fs::remove_file(vo).ok();
    }
    #[test] fn test_multiline_wide() {
        // 3-char line width, 7 bases total (last line shorter)
        let p = tmp("mlfw", ">s1\nATG\nCAT\nG\n>s2\nATG\nCCT\nG\n");
        let m = setup(&p);
        let lk = build_lookup(false);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        assert_eq!(sl, 7);
        assert!(!layout.single_line);
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, false);
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].index, 4); // position 4: A vs C
        let o = "/tmp/snpick_t_mlfw_out.fa";
        let ep = ExtractParams { records: &recs, output: o, collect_vcf: false, lookup: &lk, upper: &up, layout };
        pass2_extract(&m, &mut v, &ep).unwrap();
        let c = std::fs::read_to_string(o).unwrap();
        let l: Vec<&str> = c.lines().collect();
        assert_eq!(l[1], "A"); // s1 pos 4 = A
        assert_eq!(l[3], "C"); // s2 pos 4 = C
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }
    #[test] fn test_mixed_wrapping() {
        // s1 is single-line, s2 is multi-line — different wrapping per record
        let p = tmp("mxwr", ">s1\nATGCATGC\n>s2\nATGC\nATCC\n");
        let o = "/tmp/snpick_t_mxwr_out.fa";
        let m = setup(&p);
        let lk = build_lookup(false);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        assert_eq!(sl, 8);
        assert!(!layout.single_line); // s2 is multi-line
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, false);
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].index, 6); // pos 6: G vs C
        let ep = ExtractParams { records: &recs, output: o, collect_vcf: false, lookup: &lk, upper: &up, layout };
        pass2_extract(&m, &mut v, &ep).unwrap();
        let c = std::fs::read_to_string(o).unwrap();
        let l: Vec<&str> = c.lines().collect();
        assert_eq!(l[1], "G"); // s1 pos 6 = G
        assert_eq!(l[3], "C"); // s2 pos 6 = C
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }
    #[test] fn test_empty() {
        let p = tmp("empf", "");
        let m = setup(&p);
        assert!(index_fasta(&m).is_err());
        std::fs::remove_file(&p).ok();
    }
}
