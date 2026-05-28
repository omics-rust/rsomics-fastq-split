use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

fn bench_fastq_split(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-fastq-split");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fq = manifest.join("tests/golden/se5.fastq");
    c.bench_function("rsomics-fastq-split golden", |b| {
        b.iter(|| {
            let tmp = TempDir::new().unwrap();
            let out_prefix = tmp.path().join("split");
            let out = Command::new(black_box(bin))
                .args([
                    "--in1",
                    fq.to_str().unwrap(),
                    "--out1",
                    out_prefix.to_str().unwrap(),
                    "-s",
                    "2",
                ])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_fastq_split);
criterion_main!(benches);
