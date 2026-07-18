mod audit;
mod coords;
mod extract;
mod fasta;
mod filter;
mod input;
mod scan;
mod types;
mod vcf;

use clap::Parser;
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use crate::extract::{pass2_extract, ExtractParams};
use crate::fasta::{get_ref_seq, index_fasta, FastaRecord};
use crate::filter::{count_sites, SiteFilters};
use crate::scan::{analyze, pass1_scan};
use crate::types::*;
use crate::vcf::write_vcf;

/// Emit a `[snpick]` progress line to stderr unless `--quiet` was passed.
/// Errors are always printed; only progress chatter is gated.
macro_rules! progress {
    ($quiet:expr, $($arg:tt)*) => {
        if !$quiet { eprintln!($($arg)*); }
    };
}

// =============================================================================
// CLI
// =============================================================================

/// How to treat IUPAC ambiguity codes (R, Y, S, ...).
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq, Default)]
enum IupacMode {
    /// Treat ambiguity codes as missing data (default).
    #[default]
    Missing,
    /// Resolve ambiguity codes to their constituent bases when classifying sites.
    Resolve,
}

/// Policy for out-of-alphabet bytes in the input.
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq, Default)]
enum OnInvalid {
    /// Skip the check entirely (no cost).
    #[default]
    Ignore,
    /// Warn about out-of-alphabet bytes but continue.
    Warn,
    /// Error out if any out-of-alphabet byte is present.
    Error,
}

#[derive(Parser, Debug)]
#[command(
    name = "snpick",
    version = env!("CARGO_PKG_VERSION"),
    author = "Paula Ruiz-Rodriguez <paula.ruiz.rodriguez@csic.es>",
    about = "Fast, memory-efficient extraction of variable sites from FASTA alignments.",
    long_about = "snpick extracts variable (SNP) sites from a whole-genome FASTA alignment, \
producing a reduced alignment ready for phylogenetic inference (optionally with a VCF and the \
ASC fconst constant-site counts for IQ-TREE / RAxML). It uses a zero-copy, memory-mapped, \
parallel two-pass scan that scales to thousands of genomes with minimal RAM.",
    after_help = "NOTES:\n  \
- All input sequences must have the same length (an alignment). The first sequence is the\n    \
reference for REF/ALT polarity.\n  \
- fconst is printed as A,C,G,T (the order IQ-TREE's -fconst expects).\n  \
- In the VCF, POS is the 1-based ALIGNMENT column, not an ungapped reference coordinate.\n  \
- IUPAC ambiguous bases (N, R, Y, ...) are treated as missing data, never as alleles.\n  \
- Gaps ('-') are ignored unless -g is given, where they become a 5th allele ('*' in the VCF).\n\n\
EXAMPLES:\n  \
snpick -f aln.fasta -o snps.fasta\n  \
snpick -f aln.fasta -o snps.fasta --vcf --chrom NC_000962.3\n  \
snpick -f aln.fasta -o snps.fasta -g -t 8 -q"
)]
struct Args {
    /// Input FASTA alignment (all sequences must have equal length).
    #[arg(short, long)] fasta: String,
    /// Output FASTA containing only the variable sites (not needed with --dry-run/--check).
    #[arg(short, long, required_unless_present_any = ["dry_run", "check"])] output: Option<String>,
    /// Treat gaps ('-') as a 5th character instead of ignoring them.
    #[arg(short = 'g', long)] include_gaps: bool,
    /// Also write a VCF, named after the output (snps.fasta -> snps.vcf).
    #[arg(long)] vcf: bool,
    /// Write the VCF to a custom path (implies --vcf).
    #[arg(long)] vcf_output: Option<String>,
    /// Silence progress logs on stderr (errors are still reported).
    #[arg(short = 'q', long)] quiet: bool,
    /// Number of threads for the parallel scan (default: all logical cores).
    #[arg(short = 't', long, value_parser = parse_threads)]
    threads: Option<usize>,
    /// CHROM / contig name written to the VCF (e.g. NC_000962.3).
    #[arg(long, default_value = "1")]
    chrom: String,
    /// Sequence ID to use as the REF/polarity reference (default: first sequence).
    #[arg(long)] reference: Option<String>,
    /// Permit duplicate sequence IDs instead of erroring.
    #[arg(long)] allow_dup_ids: bool,
    /// Write a machine-readable JSON run summary to this path ('-' for stdout).
    #[arg(long)] stats_json: Option<String>,
    /// Report site statistics without writing any FASTA/VCF output.
    #[arg(long)] dry_run: bool,
    /// Keep only these sample IDs (comma-separated, or @file with one ID per line).
    #[arg(long)] keep_samples: Option<String>,
    /// Drop these sample IDs (comma-separated, or @file). Excludes --keep-samples.
    #[arg(long, conflicts_with = "keep_samples")] exclude_samples: Option<String>,
    /// Drop sites whose fraction of missing genotypes exceeds this (0.0-1.0).
    #[arg(long)] max_missing: Option<f64>,
    /// Drop sites whose minor-allele count is below this.
    #[arg(long)] mac: Option<u32>,
    /// Drop sites whose minor-allele frequency is below this (0.0-1.0).
    #[arg(long)] maf: Option<f64>,
    /// Drop sites with fewer than this many samples with data.
    #[arg(long)] min_samples: Option<u32>,
    /// Drop sites with more than this many distinct alleles (e.g. 2 = biallelic).
    #[arg(long)] max_alleles: Option<u32>,
    /// BED file of regions to mask out (excluded from output AND fconst).
    #[arg(long)] mask: Option<String>,
    /// Interpret --mask coordinates as reference positions, not alignment columns.
    #[arg(long, requires = "mask")] mask_ref: bool,
    /// Use ungapped reference positions for VCF POS instead of the alignment column.
    #[arg(long)] ref_coords: bool,
    /// Write a TSV mapping each variable site (alignment_pos, ref_pos, ref, alt) to this path.
    #[arg(long)] sites_output: Option<String>,
    /// What to do about out-of-alphabet bytes: ignore (default), warn, or error.
    #[arg(long, value_enum, default_value_t = OnInvalid::Ignore)] on_invalid: OnInvalid,
    /// Audit the alignment composition and exit without writing output.
    #[arg(long)] check: bool,
    /// How to treat IUPAC ambiguity codes: missing (default) or resolve to bases.
    #[arg(long, value_enum, default_value_t = IupacMode::Missing)] iupac_mode: IupacMode,
}

/// Parse a positive thread count (Rayon treats 0 as "use default", which would
/// be a confusing silent no-op, so reject it explicitly).
fn parse_threads(s: &str) -> Result<usize, String> {
    let n: usize = s.parse().map_err(|_| format!("'{}' is not a valid thread count", s))?;
    if n == 0 {
        return Err("thread count must be at least 1".to_string());
    }
    Ok(n)
}

// =============================================================================
// Path validation
// =============================================================================

fn resolve_path(p: &str) -> io::Result<std::path::PathBuf> {
    let path = Path::new(p);
    if path.exists() {
        return std::fs::canonicalize(path);
    }
    // A bare filename has an *empty* parent (Some("")), not None, so canonicalize
    // would fail with ENOENT. Treat that as the current directory.
    let parent = match path.parent() {
        Some(p) if !p.as_os_str().is_empty() => p,
        _ => Path::new("."),
    };
    let parent_abs = std::fs::canonicalize(parent).map_err(|e| {
        io::Error::new(e.kind(), format!("Cannot resolve parent of '{}': {}", p, e))
    })?;
    Ok(parent_abs.join(path.file_name().unwrap_or_default()))
}

fn check_paths_differ(a: &str, b: &str) -> io::Result<()> {
    if a == "-" || b == "-" {
        return Ok(()); // stdin/stdout can't collide with a file
    }
    let pa = resolve_path(a)?;
    let pb = resolve_path(b)?;
    if pa == pb {
        return Err(io::Error::new(io::ErrorKind::InvalidInput,
            format!("Paths resolve to same file: {}", pa.display())));
    }
    Ok(())
}

/// Reject empty sequence IDs, and duplicate IDs unless `allow_dups`.
/// Duplicate sample names yield spec-invalid VCFs and trees downstream tools reject.
fn validate_ids(records: &[FastaRecord], allow_dups: bool) -> io::Result<()> {
    let mut seen: HashSet<&[u8]> = HashSet::with_capacity(records.len());
    for (i, rec) in records.iter().enumerate() {
        if rec.id.is_empty() {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Record #{} has an empty sequence ID.", i + 1)));
        }
        if !allow_dups && !seen.insert(rec.id) {
            let id = std::str::from_utf8(rec.id).unwrap_or("?");
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("Duplicate sequence ID '{}' (record #{}). Use --allow-dup-ids to permit.",
                    id, i + 1)));
        }
    }
    Ok(())
}

/// Emit a flat JSON run summary to `path` ('-' = stdout). Hand-written (no serde);
/// only the string fields need escaping.
#[allow(clippy::too_many_arguments)]
fn write_stats_json(
    path: &str, input: &str, reference: &str, seq_length: usize, num_samples: usize,
    sc: &SiteCounts, written: usize, include_gaps: bool, threads: usize,
) -> io::Result<()> {
    let esc = |s: &str| s.replace('\\', "\\\\").replace('"', "\\\"");
    let c = &sc.constant;
    let json = format!(
        "{{\n  \"snpick_version\": \"{}\",\n  \"input\": \"{}\",\n  \"reference\": \"{}\",\n  \
\"sequences\": {},\n  \"alignment_length\": {},\n  \"variable_sites\": {},\n  \
\"written_sites\": {},\n  \
\"constant_sites\": {},\n  \"constant_by_base\": {{ \"A\": {}, \"C\": {}, \"G\": {}, \"T\": {} }},\n  \
\"ambiguous_sites\": {},\n  \"fconst\": [{}, {}, {}, {}],\n  \
\"include_gaps\": {},\n  \"threads\": {}\n}}\n",
        env!("CARGO_PKG_VERSION"), esc(input), esc(reference),
        num_samples, seq_length, sc.variable, written,
        c.total(), c.a, c.c, c.g, c.t, sc.ambiguous, c.a, c.c, c.g, c.t,
        include_gaps, threads,
    );
    if path == "-" {
        io::stdout().write_all(json.as_bytes())
    } else {
        std::fs::write(path, json)
    }
}

/// Parse a sample-ID selector: a comma-separated list, or `@path` to read one ID per line
/// (blank lines and `#` comments ignored).
fn parse_id_set(spec: &str) -> io::Result<HashSet<Vec<u8>>> {
    let mut set = HashSet::new();
    if let Some(path) = spec.strip_prefix('@') {
        let content = std::fs::read_to_string(path).map_err(|e| io::Error::new(e.kind(),
            format!("Cannot read sample list '{}': {}", path, e)))?;
        for line in content.lines() {
            let line = line.trim();
            if !line.is_empty() && !line.starts_with('#') {
                set.insert(line.as_bytes().to_vec());
            }
        }
    } else {
        for id in spec.split(',') {
            let id = id.trim();
            if !id.is_empty() { set.insert(id.as_bytes().to_vec()); }
        }
    }
    Ok(set)
}

/// Apply --keep-samples / --exclude-samples in place. Filtering happens before the scan, so
/// site classification and fconst are computed for exactly the retained samples. Unknown IDs
/// are an error (catches typos).
fn apply_sample_filter(
    records: &mut Vec<FastaRecord>, keep: &Option<String>, exclude: &Option<String>,
) -> io::Result<()> {
    let (spec, keeping) = match (keep, exclude) {
        (Some(s), _) => (s, true),
        (_, Some(s)) => (s, false),
        _ => return Ok(()),
    };
    let set = parse_id_set(spec)?;
    let present: HashSet<&[u8]> = records.iter().map(|r| r.id).collect();
    for id in &set {
        if !present.contains(id.as_slice()) {
            return Err(io::Error::new(io::ErrorKind::InvalidInput,
                format!("Sample '{}' not found in the alignment.", String::from_utf8_lossy(id))));
        }
    }
    records.retain(|r| set.contains(r.id) == keeping);
    if records.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidInput,
            "No sequences remain after sample filtering."));
    }
    Ok(())
}

/// Write a TSV mapping each variable site to its alignment column, reference position,
/// REF and ALT alleles (gaps rendered as '*').
fn write_sites_tsv(path: &str, var: &[VariablePosition], ref_pos: Option<&[u32]>) -> io::Result<()> {
    let mut w = BufWriter::new(File::create(path).map_err(|e| io::Error::new(e.kind(),
        format!("Cannot create sites TSV '{}': {}", path, e)))?);
    writeln!(w, "alignment_pos\tref_pos\tref\talt")?;
    for vp in var {
        let rp = ref_pos.map(|r| r[vp.index]).unwrap_or((vp.index + 1) as u32);
        let refc = if vp.ref_base == b'-' { '*' } else { vp.ref_base as char };
        let alt: String = vp.alt_bases.iter()
            .map(|&b| if b == b'-' { "*".to_string() } else { (b as char).to_string() })
            .collect::<Vec<_>>().join(",");
        writeln!(w, "{}\t{}\t{}\t{}", vp.index + 1, rp, refc, alt)?;
    }
    w.flush()
}

/// Index of the record to use as the REF/polarity reference.
fn reference_index(records: &[FastaRecord], reference: &Option<String>) -> io::Result<usize> {
    match reference {
        Some(id) => {
            let target = id.as_bytes();
            records.iter().position(|r| r.id == target).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidInput,
                    format!("--reference '{}' not found among sequence IDs.", id))
            })
        }
        None => Ok(0),
    }
}

// =============================================================================
// Pipeline
// =============================================================================

fn run() -> io::Result<()> {
    let args = Args::parse();
    let quiet = args.quiet;

    // Configure the Rayon pool up front. Absent, Rayon uses all logical cores;
    // pinning it keeps wall-clock deterministic on shared HPC/SLURM nodes.
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().map_err(|e| {
            io::Error::other(format!("Cannot configure {} thread(s): {}", n, e))
        })?;
    }

    let start = Instant::now();
    let lookup = build_lookup(args.include_gaps);
    let upper = build_upper();

    let dry_run = args.dry_run;
    // --dry-run reports statistics only, so it writes no FASTA/VCF.
    let do_vcf = (args.vcf || args.vcf_output.is_some()) && !dry_run;

    // A CHROM with whitespace would break the tab-delimited VCF columns.
    if do_vcf && (args.chrom.is_empty() || args.chrom.bytes().any(|b| b.is_ascii_whitespace())) {
        return Err(io::Error::new(io::ErrorKind::InvalidInput,
            "--chrom must be non-empty and contain no whitespace."));
    }

    // Output path (clap guarantees it is present unless --dry-run).
    let out_path: Option<String> = if dry_run { None } else { args.output.clone() };

    // Validate paths (skip when writing nothing).
    if let Some(ref out) = out_path {
        check_paths_differ(&args.fasta, out)?;
    }
    let vcf_path = if do_vcf {
        let out = out_path.as_deref().unwrap_or("output");
        let vp = args.vcf_output.clone().unwrap_or_else(|| {
            let o = Path::new(out);
            let stem = o.file_stem().and_then(|s| s.to_str()).unwrap_or("output");
            let parent = o.parent().unwrap_or(Path::new("."));
            parent.join(format!("{}.vcf", stem)).to_string_lossy().into_owned()
        });
        check_paths_differ(&args.fasta, &vp)?;
        check_paths_differ(out, &vp)?;
        Some(vp)
    } else { None };

    // Map the input, transparently decompressing gzip/bgzip and reading stdin ("-") as needed.
    let input = input::map_input(&args.fasta)?;
    if input.mmap.is_empty() {
        return Err(io::Error::new(io::ErrorKind::InvalidData,
            format!("Input '{}' is empty (0 bytes).", args.fasta)));
    }
    // Hint: pass 1 reads sequentially; OS can prefetch and release pages eagerly
    input.mmap.advise(memmap2::Advice::Sequential).ok();
    let data = &input.mmap[..];

    // Index records
    let (mut records, seq_length, layout) = index_fasta(data)?;

    validate_ids(&records, args.allow_dup_ids)?;
    apply_sample_filter(&mut records, &args.keep_samples, &args.exclude_samples)?;
    let num_samples = records.len();

    let ref_idx = reference_index(&records, &args.reference)?;
    let ref_name = std::str::from_utf8(records[ref_idx].id).unwrap_or("reference").to_string();

    progress!(quiet, "[snpick] Mapped {} bytes. {} sequences × {} positions.{}",
        data.len(), num_samples, seq_length,
        if layout.single_line { "" } else { " (multi-line FASTA)" });

    // Composition audit for --on-invalid / --check (a full extra pass, only on demand).
    if args.check || args.on_invalid != OnInvalid::Ignore {
        let hist = audit::composition(data, &records, seq_length, layout);
        let invalid = audit::invalid_count(&hist);
        if args.on_invalid == OnInvalid::Error && invalid > 0 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, format!(
                "{} out-of-alphabet byte(s) found — is this a nucleotide alignment? \
                 (use --on-invalid ignore to allow).", invalid)));
        }
        if args.on_invalid == OnInvalid::Warn && invalid > 0 {
            eprintln!("[snpick] Warning: {} out-of-alphabet byte(s) in the alignment.", invalid);
        }
        if args.check {
            progress!(quiet, "[snpick] Composition: {}", audit::summary(&hist));
            return Ok(());
        }
    }

    // Pass 1: bitmask scan. Under --iupac-mode resolve, ambiguity codes are expanded to their
    // bases *only here* (for classification); NS counting and REF selection keep the strict table.
    let scan_lookup = if args.iupac_mode == IupacMode::Resolve {
        build_iupac_lookup(args.include_gaps)
    } else {
        lookup
    };
    let mut bitmask = pass1_scan(data, &records, seq_length, layout, &scan_lookup);
    let ref_seq = get_ref_seq(data, &records[ref_idx], seq_length, layout);
    let t1 = start.elapsed().as_secs_f64();

    // Ungapped reference positions (only computed when a feature needs them).
    let ref_pos: Option<Vec<u32>> =
        if args.ref_coords || args.sites_output.is_some() || args.mask_ref {
            Some(coords::ref_positions(&ref_seq))
        } else {
            None
        };

    // Mask out BED regions before classification: masked columns are zeroed, so they classify
    // as ambiguous and never enter the output or the fconst constant counts.
    if let Some(bed) = &args.mask {
        let ivs = coords::parse_bed(bed)?;
        let rc = if args.mask_ref { ref_pos.as_deref() } else { None };
        let mask = coords::build_mask(&ivs, seq_length, rc);
        let masked = mask.iter().filter(|&&m| m).count();
        for (col, &m) in mask.iter().enumerate() {
            if m { bitmask[col] = 0; }
        }
        progress!(quiet, "[snpick] Masked {} columns from {} BED region(s).", masked, ivs.len());
    }

    let (mut var_positions, site_counts) = analyze(&bitmask, &ref_seq, &lookup, args.include_gaps);

    drop(bitmask);
    drop(ref_seq);

    progress!(quiet, "[snpick] {} variable, {} constant ({}), {} ambiguous-only, {} total.",
        site_counts.variable, site_counts.constant.total(), site_counts.constant,
        site_counts.ambiguous, seq_length);
    progress!(quiet, "[snpick] ASC fconst: {}", site_counts.constant.fconst());
    progress!(quiet, "[snpick] Pass 1 took {:.2}s.", t1);

    // Per-site filtering (opt-in). Removes variable sites from the OUTPUT only; the fconst
    // constant-site counts are unchanged, so ASC stays valid.
    let filters = SiteFilters {
        max_missing: args.max_missing, mac: args.mac, maf: args.maf,
        min_samples: args.min_samples, max_alleles: args.max_alleles,
    };
    if filters.active() && !var_positions.is_empty() {
        let pos_indices: Vec<usize> = var_positions.iter().map(|v| v.index).collect();
        let stats = count_sites(data, &records, &pos_indices, layout, &lookup);
        let ns_total = num_samples as u32;
        let before = var_positions.len();
        let mut i = 0;
        var_positions.retain(|_| {
            let keep = stats[i].passes(&filters, ns_total, args.include_gaps);
            i += 1;
            keep
        });
        progress!(quiet, "[snpick] Filtered {} of {} variable sites ({} remain).",
            before - var_positions.len(), before, var_positions.len());
    }
    let num_var = var_positions.len();

    // Machine-readable stats sidecar (also emitted for dry-run and zero-variant).
    if let Some(ref sp) = args.stats_json {
        write_stats_json(sp, &args.fasta, &ref_name, seq_length, num_samples,
            &site_counts, num_var, args.include_gaps, rayon::current_num_threads())?;
        progress!(quiet, "[snpick] Stats JSON written to {}.",
            if sp == "-" { "stdout" } else { sp });
    }

    if dry_run {
        progress!(quiet, "[snpick] Dry run — no output written.");
        return Ok(());
    }

    // Past the dry-run gate, an output path is guaranteed (clap-enforced).
    let out = out_path.as_deref().expect("output is required unless --dry-run");
    let pos_map = if args.ref_coords { ref_pos.as_deref() } else { None };

    // Optional variable-site coordinate map (alignment column -> reference position).
    if let Some(sp) = &args.sites_output {
        write_sites_tsv(sp, &var_positions, ref_pos.as_deref())?;
        progress!(quiet, "[snpick] Sites TSV written to {}.", sp);
    }

    // Handle zero-variant case
    if num_var == 0 {
        progress!(quiet, "[snpick] No variable positions — writing empty output.");
        let mut w = BufWriter::new(crate::extract::open_sink(out)?);
        for rec in &records {
            w.write_all(b">")?;
            w.write_all(rec.id)?;
            if !rec.desc.is_empty() { w.write_all(b" ")?; w.write_all(rec.desc)?; }
            writeln!(w)?; writeln!(w)?;
        }
        w.flush()?;
        // Still honour a requested VCF: emit a valid header-only file so a
        // pipeline that declares the .vcf as an output doesn't break.
        if let Some(ref vp) = vcf_path {
            write_vcf(&[], num_samples, &[], vp, &records, seq_length, &args.chrom, &ref_name, pos_map)?;
            progress!(quiet, "[snpick] VCF written to {} (header only — no variable sites).", vp);
        }
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
        records: &records, output: out,
        collect_vcf: do_vcf, lookup: &lookup, upper: &upper, layout,
    };
    let vcf_geno = pass2_extract(data, &mut var_positions, &ep)?;
    progress!(quiet, "[snpick] Pass 2: Wrote {} sequences to {}.", num_samples, out);

    // Write VCF
    if let (Some(ref geno), Some(ref vp)) = (&vcf_geno, &vcf_path) {
        write_vcf(geno, num_samples, &var_positions, vp, &records, seq_length, &args.chrom, &ref_name, pos_map)?;
        progress!(quiet, "[snpick] VCF written to {}.", vp);
    }

    progress!(quiet, "[snpick] Done in {:.2}s. {} vars from {} seqs × {} pos.",
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
    use memmap2::Mmap;
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
        write_vcf(&g, recs.len(), &v, vo, &recs, sl, "1", "ref", None).unwrap();
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
        write_vcf(&g, recs.len(), &v, vo, &recs, sl, "1", "ref", None).unwrap();
        let c = std::fs::read_to_string(vo).unwrap();
        let dl: Vec<&str> = c.lines().filter(|l| !l.starts_with('#')).collect();
        let f: Vec<&str> = dl[0].split('\t').collect();
        assert_eq!(f[7], "NS=2"); assert_eq!(f[11], ".");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(fo).ok(); std::fs::remove_file(vo).ok();
    }

    #[test] fn test_vcf_gap_ref() {
        // Gap in the reference at a variable site: REF must be '*' (valid VCF v4.2),
        // never a bare '-'.
        let p = tmp("vgref", ">ref\nA-GC\n>s1\nATGC\n>s2\nA-GT\n");
        let fo = "/tmp/snpick_t_vgref_out.fa"; let vo = "/tmp/snpick_t_vgref.vcf";
        let m = setup(&p);
        let lk = build_lookup(true);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, true);
        let ep = ExtractParams { records: &recs, output: fo, collect_vcf: true, lookup: &lk, upper: &up, layout };
        let g = pass2_extract(&m, &mut v, &ep).unwrap().unwrap();
        write_vcf(&g, recs.len(), &v, vo, &recs, sl, "1", "ref", None).unwrap();
        let c = std::fs::read_to_string(vo).unwrap();
        let dl: Vec<&str> = c.lines().filter(|l| !l.starts_with('#')).collect();
        let f: Vec<&str> = dl[0].split('\t').collect();
        assert_eq!(f[1], "2");            // POS (gap-in-ref site)
        assert_eq!(f[3], "*");            // REF: gap → '*', not '-'
        assert_eq!(f[4], "T");            // ALT
        assert!(!c.contains("\t-\t"));    // no bare hyphen field anywhere
        std::fs::remove_file(&p).ok(); std::fs::remove_file(fo).ok(); std::fs::remove_file(vo).ok();
    }

    #[test] fn test_vcf_gap_alt() {
        // Gap sample at a multi-allelic site renders as the '*' ALT allele.
        let p = tmp("vgalt", ">ref\nATGC\n>s1\nTTGC\n>s2\n-TGC\n");
        let fo = "/tmp/snpick_t_vgalt_out.fa"; let vo = "/tmp/snpick_t_vgalt.vcf";
        let m = setup(&p);
        let lk = build_lookup(true);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, true);
        let ep = ExtractParams { records: &recs, output: fo, collect_vcf: true, lookup: &lk, upper: &up, layout };
        let g = pass2_extract(&m, &mut v, &ep).unwrap().unwrap();
        write_vcf(&g, recs.len(), &v, vo, &recs, sl, "1", "ref", None).unwrap();
        let c = std::fs::read_to_string(vo).unwrap();
        let dl: Vec<&str> = c.lines().filter(|l| !l.starts_with('#')).collect();
        let f: Vec<&str> = dl[0].split('\t').collect();
        assert_eq!(f[1], "1");            // POS
        assert_eq!(f[3], "A");            // REF
        assert_eq!(f[4], "T,*");          // ALT: gap sample → '*'
        std::fs::remove_file(&p).ok(); std::fs::remove_file(fo).ok(); std::fs::remove_file(vo).ok();
    }

    #[test] fn test_vcf_header_only() {
        // Zero variable sites but VCF requested: a valid header-only VCF that
        // still lists every sample, so downstream pipelines find the file.
        let p = tmp("hdr", ">s1\nATGC\n>s2\nATGC\n");
        let vo = "/tmp/snpick_t_hdr.vcf";
        let m = setup(&p);
        let (recs, sl, _layout) = index_fasta(&m).unwrap();
        write_vcf(&[], recs.len(), &[], vo, &recs, sl, "1", "ref", None).unwrap();
        let c = std::fs::read_to_string(vo).unwrap();
        let data: Vec<&str> = c.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(data.len(), 0);
        assert!(c.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2"));
        std::fs::remove_file(&p).ok(); std::fs::remove_file(vo).ok();
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

    #[test] fn test_validate_ids() {
        let p = tmp("ids", ">a\nAC\n>a\nAT\n>b\nAG\n");
        let m = setup(&p);
        let (recs, _, _) = index_fasta(&m).unwrap();
        assert!(validate_ids(&recs, false).is_err());   // duplicate 'a'
        assert!(validate_ids(&recs, true).is_ok());      // allowed
        std::fs::remove_file(&p).ok();
    }

    #[test] fn test_stats_json() {
        let sc = SiteCounts {
            constant: ConstantSiteCounts { a: 7, c: 7, g: 7, t: 4 },
            variable: 5,
            ambiguous: 1,
        };
        let p = "/tmp/snpick_t_stats.json";
        write_stats_json(p, "aln.fasta", "H37Rv", 30, 6, &sc, 3, false, 8).unwrap();
        let c = std::fs::read_to_string(p).unwrap();
        assert!(c.contains("\"variable_sites\": 5"));
        assert!(c.contains("\"written_sites\": 3"));
        assert!(c.contains("\"fconst\": [7, 7, 7, 4]"));
        assert!(c.contains("\"reference\": \"H37Rv\""));
        assert!(c.contains("\"sequences\": 6"));
        std::fs::remove_file(p).ok();
    }

    #[test] fn test_count_sites() {
        // pos0 singleton (A5 T1), pos1 common (A3 C3), pos2 triallelic (A2 C2 G2).
        let p = tmp("cnt", ">s1\nAAA\n>s2\nAAA\n>s3\nAAC\n>s4\nACC\n>s5\nACG\n>s6\nTCG\n");
        let m = setup(&p);
        let lk = build_lookup(false);
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (v, _) = analyze(&bm, &rs, &lk, false);
        let idx: Vec<usize> = v.iter().map(|x| x.index).collect();
        let stats = crate::filter::count_sites(&m, &recs, &idx, layout, &lk);
        assert_eq!(stats.len(), 3);
        assert_eq!(stats[0].counts, [5, 0, 0, 1]);
        assert_eq!(stats[1].counts, [3, 3, 0, 0]);
        assert_eq!(stats[2].counts, [2, 2, 2, 0]);
        std::fs::remove_file(&p).ok();
    }

    #[test] fn test_exclude_samples_reclassifies() {
        // A site variable only because of one sample becomes constant (and enters
        // fconst) once that sample is excluded — the correctness snp-sites misses.
        let p = tmp("excl", ">a\nAA\n>b\nAA\n>c\nAT\n");
        let m = setup(&p);
        let lk = build_lookup(false);
        let (mut recs, sl, layout) = index_fasta(&m).unwrap();
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (v, sc) = analyze(&bm, &rs, &lk, false);
        assert_eq!(v.len(), 1);
        assert_eq!(sc.constant.total(), 1);
        apply_sample_filter(&mut recs, &None, &Some("c".to_string())).unwrap();
        assert_eq!(recs.len(), 2);
        let bm2 = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs2 = get_ref_seq(&m, &recs[0], sl, layout);
        let (v2, sc2) = analyze(&bm2, &rs2, &lk, false);
        assert_eq!(v2.len(), 0);
        assert_eq!(sc2.constant.total(), 2);
        std::fs::remove_file(&p).ok();
    }

    #[test] fn test_keep_samples_and_unknown() {
        let p = tmp("keep", ">a\nAT\n>b\nCT\n>c\nGT\n");
        let m = setup(&p);
        let (mut recs, _, _) = index_fasta(&m).unwrap();
        apply_sample_filter(&mut recs, &Some("a,c".to_string()), &None).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, b"a");
        assert_eq!(recs[1].id, b"c");
        let (mut recs2, _, _) = index_fasta(&m).unwrap();
        assert!(apply_sample_filter(&mut recs2, &Some("a,zzz".to_string()), &None).is_err());
        std::fs::remove_file(&p).ok();
    }

    #[test] fn test_reference_index() {
        let p = tmp("refi", ">a\nAC\n>b\nAT\n>c\nAG\n");
        let m = setup(&p);
        let (recs, _, _) = index_fasta(&m).unwrap();
        assert_eq!(reference_index(&recs, &None).unwrap(), 0);
        assert_eq!(reference_index(&recs, &Some("b".to_string())).unwrap(), 1);
        assert!(reference_index(&recs, &Some("zzz".to_string())).is_err());
        std::fs::remove_file(&p).ok();
    }

    #[test] fn test_reference_polarity() {
        // Choosing a different reference flips which base is REF vs ALT.
        let p = tmp("refp", ">a\nAT\n>b\nCT\n");
        let m = setup(&p);
        let lk = build_lookup(false);
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        // default reference = record 0 ('a') -> REF A
        let rs0 = get_ref_seq(&m, &recs[0], sl, layout);
        let (v0, _) = analyze(&bm, &rs0, &lk, false);
        assert_eq!(v0[0].ref_base, b'A');
        assert_eq!(v0[0].alt_bases, vec![b'C']);
        // reference = record 1 ('b') -> REF C
        let rs1 = get_ref_seq(&m, &recs[1], sl, layout);
        let (v1, _) = analyze(&bm, &rs1, &lk, false);
        assert_eq!(v1[0].ref_base, b'C');
        assert_eq!(v1[0].alt_bases, vec![b'A']);
        std::fs::remove_file(&p).ok();
    }

    #[test] fn test_all_gap_column_counted() {
        // Under --include-gaps an all-gap column is tallied as ambiguous, so
        // variable + constant + ambiguous still equals seq_length.
        let p = tmp("allgap", ">s1\nA-\n>s2\nA-\n");
        let m = setup(&p);
        let lk = build_lookup(true);
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (v, sc) = analyze(&bm, &rs, &lk, true);
        assert_eq!(v.len(), 0);
        assert_eq!(sc.constant.total(), 1);
        assert_eq!(sc.ambiguous, 1);
        assert_eq!(sc.variable + sc.constant.total() + sc.ambiguous, sl);
        std::fs::remove_file(&p).ok();
    }

    #[test] fn test_paths() {
        let p = tmp("pdg", ">x\nA\n");
        assert!(check_paths_differ(&p, "/tmp/snpick_t_pdg2.fa").is_ok());
        assert!(check_paths_differ(&p, &p).is_err());
        // A bare (directory-less) output name must resolve to the cwd, not error.
        assert!(check_paths_differ(&p, "snpick_t_pdg_bare_out.fa").is_ok());
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
        write_vcf(&g, recs.len(), &v, vo, &recs, sl, "1", "ref", None).unwrap();
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

    #[test] fn test_single_line_blank_line() {
        // Blank line after the header in an otherwise single-line file: must not
        // drop the SNP. The record is forced onto the newline-skipping path.
        let p = tmp("blank", ">s1\n\nAAAA\n>s2\n\nAAAT\n");
        let o = "/tmp/snpick_t_blank_out.fa";
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
        assert_eq!(v[0].index, 3);
        let ep = ExtractParams { records: &recs, output: o, collect_vcf: false, lookup: &lk, upper: &up, layout };
        pass2_extract(&m, &mut v, &ep).unwrap();
        let c = std::fs::read_to_string(o).unwrap();
        let l: Vec<&str> = c.lines().collect();
        assert_eq!(l[1], "A"); assert_eq!(l[3], "T");
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

    #[test] fn test_ambiguous_ref_base() {
        // Reference (first sequence) is N at a variable site: REF falls back to
        // the first observed base in A,C,G,T order and the ref genotype is '.'.
        let p = tmp("nref", ">ref\nNTGC\n>s1\nATGC\n>s2\nCTGC\n");
        let fo = "/tmp/snpick_t_nref_out.fa"; let vo = "/tmp/snpick_t_nref.vcf";
        let m = setup(&p);
        let lk = build_lookup(false);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, false);
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].ref_base, b'A');
        assert_eq!(v[0].alt_bases, vec![b'C']);
        let ep = ExtractParams { records: &recs, output: fo, collect_vcf: true, lookup: &lk, upper: &up, layout };
        let g = pass2_extract(&m, &mut v, &ep).unwrap().unwrap();
        write_vcf(&g, recs.len(), &v, vo, &recs, sl, "1", "ref", None).unwrap();
        let c = std::fs::read_to_string(vo).unwrap();
        let dl: Vec<&str> = c.lines().filter(|l| !l.starts_with('#')).collect();
        let f: Vec<&str> = dl[0].split('\t').collect();
        assert_eq!(f[3], "A");        // REF fallback
        assert_eq!(f[4], "C");        // ALT
        assert_eq!(f[7], "NS=2");     // ref N excluded
        assert_eq!(f[9], ".");        // ref sample genotype missing
        assert_eq!(f[10], "0");       // s1 = A = ref
        assert_eq!(f[11], "1");       // s2 = C = alt
        std::fs::remove_file(&p).ok(); std::fs::remove_file(fo).ok(); std::fs::remove_file(vo).ok();
    }

    #[test] fn test_lowercase_normalization() {
        // Soft-masked (lowercase) input classifies case-insensitively and is
        // written uppercase in the reduced FASTA.
        let p = tmp("lower", ">ref\natgc\n>s1\natgt\n");
        let o = "/tmp/snpick_t_lower_out.fa";
        let m = setup(&p);
        let lk = build_lookup(false);
        let up = build_upper();
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (mut v, _) = analyze(&bm, &rs, &lk, false);
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].index, 3);
        let ep = ExtractParams { records: &recs, output: o, collect_vcf: false, lookup: &lk, upper: &up, layout };
        pass2_extract(&m, &mut v, &ep).unwrap();
        let c = std::fs::read_to_string(o).unwrap();
        let l: Vec<&str> = c.lines().collect();
        assert_eq!(l[1], "C"); assert_eq!(l[3], "T");
        std::fs::remove_file(&p).ok(); std::fs::remove_file(o).ok();
    }

    #[test] fn test_inconsistent_length_error() {
        // Sequences of differing lengths are a malformed alignment → error.
        let p = tmp("badlen", ">s1\nATGC\n>s2\nATG\n");
        let m = setup(&p);
        assert!(index_fasta(&m).is_err());
        std::fs::remove_file(&p).ok();
    }

    #[test] fn test_crlf_single_line() {
        // CRLF line endings on single-line sequences use the fast path.
        let p = tmp("crlf1", "");
        std::fs::write(&p, b">s1\r\nATGC\r\n>s2\r\nATCC\r\n").unwrap();
        let m = setup(&p);
        let lk = build_lookup(false);
        let (recs, sl, layout) = index_fasta(&m).unwrap();
        assert_eq!(sl, 4);
        assert!(layout.single_line);
        let bm = pass1_scan(&m, &recs, sl, layout, &lk);
        let rs = get_ref_seq(&m, &recs[0], sl, layout);
        let (v, _) = analyze(&bm, &rs, &lk, false);
        assert_eq!(v.len(), 1);
        assert_eq!(v[0].index, 2);
        std::fs::remove_file(&p).ok();
    }

    #[test] fn test_empty() {
        let p = tmp("empg", "");
        let m = setup(&p);
        assert!(index_fasta(&m).is_err());
        std::fs::remove_file(&p).ok();
    }
}
