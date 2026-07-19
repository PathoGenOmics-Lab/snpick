//! Criterion micro-benchmark for the pass-1 bitmask scan. Advisory only — run locally with
//! `cargo bench`; not a CI gate (shared runners produce noisy thresholds).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use snpick::fasta::index_fasta;
use snpick::scan::pass1_scan;
use snpick::types::build_lookup;

fn make_alignment(nseq: usize, len: usize) -> Vec<u8> {
    let mut fa = Vec::new();
    for i in 0..nseq {
        fa.extend_from_slice(format!(">s{}\n", i).as_bytes());
        let mut seq = vec![b'A'; len];
        seq[i % len] = b'C';
        seq[(i * 7) % len] = b'G';
        fa.extend_from_slice(&seq);
        fa.push(b'\n');
    }
    fa
}

fn bench_scan(c: &mut Criterion) {
    let fa = make_alignment(500, 20_000);
    let lookup = build_lookup(false);
    let (recs, sl, layout) = index_fasta(&fa).unwrap();
    c.bench_function("pass1_scan_500x20000", |b| {
        b.iter(|| pass1_scan(black_box(&fa), &recs, sl, layout, &lookup))
    });
}

criterion_group!(benches, bench_scan);
criterion_main!(benches);
