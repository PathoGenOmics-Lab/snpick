mod extract;
mod fasta;
mod scan;
mod types;
mod vcf;

use clap::Parser;
use memmap2::Mmap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use crate::extract::{pass2_extract, ExtractParams};
use crate::fasta::{get_ref_seq, index_fasta};
use crate::scan::{analyze, pass1_scan};
use crate::types::*;
use crate::vcf::write_vcf;

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
// Path validation
// =============================================================================

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
// Pipeline
// =============================================================================

fn run() -> io::Result<()> {
    let args = Args::parse();
    let start = Instant::now();
    let lookup = build_lookup(args.include_gaps);
    let upper = build_upper();

    let do_vcf = args.vcf || args.vcf_output.is_some();

    // Validate paths
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

    // Memory-map input
    let file = File::open(&args.fasta).map_err(|e| io::Error::new(e.kind(),
        format!("Cannot open '{}': {}", args.fasta, e)))?;
    let file_len = file.metadata()?.len();
    if file_len == 0 {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            format!("Input file '{}' is empty (0 bytes).", args.fasta)));
    }
    let mmap = unsafe { Mmap::map(&file).map_err(|e| io::Error::new(e.kind(),
        format!("Cannot memory-map '{}': {}", args.fasta, e)))? };
    // Hint: pass 1 reads sequentially; OS can prefetch and release pages eagerly
    mmap.advise(memmap2::Advice::Sequential).ok();
    let data = &mmap[..];

    // Index records
    let (records, seq_length, layout) = index_fasta(data)?;
    let num_samples = records.len();

    eprintln!("[snpick] Mapped {} bytes. {} sequences × {} positions.{}",
        data.len(), num_samples, seq_length,
        if layout.single_line { "" } else { " (multi-line FASTA)" });

    // Pass 1: bitmask scan
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

    // Handle zero-variant case
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

    // VCF size guard
    if do_vcf {
        let geno_bytes = num_var.saturating_mul(num_samples);
        if geno_bytes > MAX_VCF_GENO_BYTES {
            return Err(io::Error::new(io::ErrorKind::InvalidInput,
                format!("VCF genotype matrix would require {} GB ({} vars × {} samples). \
                    Use without --vcf or reduce input.",
                    geno_bytes / 1_000_000_000, num_var, num_samples)));
        }
    }

    // Pass 2: extract variable sites
    let ep = ExtractParams {
        records: &records, output: &args.output,
        collect_vcf: do_vcf, lookup: &lookup, upper: &upper, layout,
    };
    let vcf_geno = pass2_extract(data, &mut var_positions, &ep)?;

    // Write VCF
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
    use crate::extract::ExtractParams;
    use crate::fasta::{get_ref_seq, index_fasta};
    use crate::scan::{analyze, pass1_scan};
    use crate::vcf::write_vcf;

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
        assert_eq!(lk[b'A' as usize], BIT_A);
        assert_eq!(lk[b'N' as usize], 0);
        assert_eq!(build_lookup(true)[b'-' as usize], BIT_GAP);
    }

    #[test] fn test_index() {
        let p = tmp("idx3", ">s1\nATGC\n>s2\nATCC\n");
        let m = setup(&p);
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(sl, 4);
        assert!(layout.single_line);
        assert_eq!(recs[0].id, b"s1");
        assert_eq!(recs[1].id, b"s2");
        std::fs::remove_file(&p).ok();
    }

    #[test] fn test_pass1() {
        let p = tmp("p1g", ">s1\nATGC\n>s2\nATCC\n");
        let m = setup(&p);
        let lk = build_lookup(false);
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        assert_eq!(bm[2], BIT_G | BIT_C);
        std::fs::remove_file(&p).ok();
    }

    #[test] fn test_e2e() {
        let p = tmp("e2eg", ">ref\nATGCATGC\n>s1\nATGTATGC\n>s2\nACGCATGC\n");
        let o = "/tmp/snpick_t_e2eg_out.fa";
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
        let p = tmp("vcfg", ">ref\nATGC\n>s1\nCTGC\n>s2\nGTGC\n>s3\nTTGC\n");
        let fo = "/tmp/snpick_t_vcfg_out.fa"; let vo = "/tmp/snpick_t_vcfg.vcf";
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
        assert_eq!(f[3], "A"); assert_eq!(f[4], "C,G,T");
        assert_eq!(f[7], "NS=4"); assert_eq!(f[9], "0"); assert_eq!(f[12], "3");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(fo).ok(); std::fs::remove_file(vo).ok();
    }

    #[test] fn test_ambig() {
        let p = tmp("ambg", ">ref\nATGC\n>s1\nCTGC\n>s2\nNTGC\n");
        let fo = "/tmp/snpick_t_ambg_out.fa"; let vo = "/tmp/snpick_t_ambg.vcf";
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
        let p = tmp("descg", ">s1 some description\nATGC\n>s2 another desc\nATCC\n");
        let o = "/tmp/snpick_t_descg_out.fa";
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
        let p = tmp("gapg", ">ref\nATGC\n>s1\nA-GC\n");
        let m = setup(&p);
        let lk_no = build_lookup(false);
        let lk_yes = build_lookup(true);
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        let bm1 = pass1_scan(&m, &recs, sl, layout, &lk_no);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (v1, _) = analyze(&bm1, &rs, &lk_no, false);
        assert!(v1.is_empty());
        let bm2 = pass1_scan(&m, &recs, sl, layout, &lk_yes);
        let (v2, _) = analyze(&bm2, &rs, &lk_yes, true);
        assert_eq!(v2.len(), 1);
        std::fs::remove_file(&p).ok();
    }

    #[test] fn test_paths() {
        let p = tmp("pdg", ">x\nA\n");
        assert!(check_paths_differ(&p, "/tmp/snpick_t_pdg2.fa").is_ok());
        assert!(check_paths_differ(&p, &p).is_err());
        std::fs::remove_file(&p).ok();
    }

    #[test] fn test_multiline_pass1() {
        let p = tmp("mlg1", ">s1\nAT\nGC\n>s2\nAT\nCC\n");
        let m = setup(&p);
        let lk = build_lookup(false);
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        assert_eq!(sl, 4);
        assert!(!layout.single_line);
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (v, _) = analyze(&bm, &rs, &lk, false);
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].index, 2);
        std::fs::remove_file(&p).ok();
    }

    #[test] fn test_multiline_pass2() {
        let p = tmp("mlg2", ">s1\nAT\nGC\n>s2\nAT\nCC\n");
        let o = "/tmp/snpick_t_mlg2_out.fa";
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
        assert_eq!(l[1], "G"); assert_eq!(l[3], "C");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }

    #[test] fn test_multiline_vcf() {
        let p = tmp("mlg3", ">ref\nAT\nGC\n>s1\nAT\nCC\n>s2\nCT\nGC\n");
        let fo = "/tmp/snpick_t_mlg3_out.fa"; let vo = "/tmp/snpick_t_mlg3.vcf";
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
        let c = std::fs::read_to_string(fo).unwrap();
        let l: Vec<&str> = c.lines().collect();
        assert_eq!(l[1], "AG"); assert_eq!(l[3], "AC"); assert_eq!(l[5], "CG");
        let vc = std::fs::read_to_string(vo).unwrap();
        let dl: Vec<&str> = vc.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(dl.len(), 2);
        std::fs::remove_file(&p).ok(); std::fs::remove_file(fo).ok(); std::fs::remove_file(vo).ok();
    }

    #[test] fn test_multiline_wide() {
        let p = tmp("mlgw", ">s1\nATG\nCAT\nG\n>s2\nATG\nCCT\nG\n");
        let o = "/tmp/snpick_t_mlgw_out.fa";
        let m = setup(&p);
        let lk = build_lookup(false);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        assert_eq!(sl, 7); assert!(!layout.single_line);
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, false);
        assert_eq!(v.len(), 1); assert_eq!(v[0].index, 4);
        let ep = ExtractParams { records: &recs, output: o, collect_vcf: false, lookup: &lk, upper: &up, layout };
        pass2_extract(&m, &mut v, &ep).unwrap();
        let c = std::fs::read_to_string(o).unwrap();
        let l: Vec<&str> = c.lines().collect();
        assert_eq!(l[1], "A"); assert_eq!(l[3], "C");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }

    #[test] fn test_mixed_wrapping() {
        let p = tmp("mxwg", ">s1\nATGCATGC\n>s2\nATGC\nATCC\n");
        let o = "/tmp/snpick_t_mxwg_out.fa";
        let m = setup(&p);
        let lk = build_lookup(false);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        assert_eq!(sl, 8); assert!(!layout.single_line);
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, false);
        assert_eq!(v.len(), 1); assert_eq!(v[0].index, 6);
        let ep = ExtractParams { records: &recs, output: o, collect_vcf: false, lookup: &lk, upper: &up, layout };
        pass2_extract(&m, &mut v, &ep).unwrap();
        let c = std::fs::read_to_string(o).unwrap();
        let l: Vec<&str> = c.lines().collect();
        assert_eq!(l[1], "G"); assert_eq!(l[3], "C");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }

    #[test] fn test_crlf_multiline() {
        // Windows line endings (\r\n) with multi-line wrapping
        let p = tmp("crlfml", "");
        std::fs::write(&p, b">s1\r\nAT\r\nGC\r\n>s2\r\nAT\r\nCC\r\n").unwrap();
        let m = setup(&p);
        let lk = build_lookup(false);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        assert_eq!(sl, 4);
        assert!(!layout.single_line);
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, false);
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].index, 2);
        let o = "/tmp/snpick_t_crlfml_out.fa";
        let ep = ExtractParams { records: &recs, output: o, collect_vcf: false, lookup: &lk, upper: &up, layout };
        pass2_extract(&m, &mut v, &ep).unwrap();
        let c = std::fs::read_to_string(o).unwrap();
        let l: Vec<&str> = c.lines().collect();
        assert_eq!(l[1], "G"); assert_eq!(l[3], "C");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }

    #[test] fn test_single_sequence() {
        // Single sequence should produce 0 variable sites
        let p = tmp("sing", ">s1\nATGC\n");
        let m = setup(&p);
        let lk = build_lookup(false);
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        assert_eq!(recs.len(), 1);
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (v, sc) = analyze(&bm, &rs, &lk, false);
        assert!(v.is_empty());
        assert_eq!(sc.constant.total(), 4);
        assert_eq!(sc.variable, 0);
        std::fs::remove_file(&p).ok();
    }

    #[test] fn test_no_trailing_newline() {
        let p = tmp("noeof", "");
        std::fs::write(&p, b">s1\nATGC\n>s2\nATCC").unwrap();
        let m = setup(&p);
        let lk = build_lookup(false);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        assert_eq!(recs.len(), 2);
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, false);
        assert_eq!(v.len(), 1);
        let o = "/tmp/snpick_t_noeof_out.fa";
        let ep = ExtractParams { records: &recs, output: o, collect_vcf: false, lookup: &lk, upper: &up, layout };
        pass2_extract(&m, &mut v, &ep).unwrap();
        let c = std::fs::read_to_string(o).unwrap();
        let l: Vec<&str> = c.lines().collect();
        assert_eq!(l[1], "G"); assert_eq!(l[3], "C");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }

    #[test] fn test_empty() {
        let p = tmp("empg", "");
        let m = setup(&p);
        assert!(index_fasta(&m).is_err());
        std::fs::remove_file(&p).ok();
    }
}
