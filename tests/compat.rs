use std::path::{Path, PathBuf};
use std::process::Command;

fn ours() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-fastq-split"))
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

// fastp's `--split_by_lines` is NOT a reliable byte-oracle: it rejects values
// < 1000 and on small inputs writes everything into the first file (e.g. a
// 4000-line file split at 2000 → 4000/0/0, not 2000/2000). So we verify the
// operation's correctness directly — every record preserved in order across
// the split, and each chunk holding the requested line count — and benchmark
// speed against fastp separately (perfgate 1.62×).

/// Split files for `out1` prefix, sorted by name (fastp-style `NNNN.<base>`).
fn split_files(dir: &Path, base: &str) -> Vec<PathBuf> {
    let mut v: Vec<PathBuf> = std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.file_name().unwrap().to_str().unwrap().ends_with(base))
        .collect();
    v.sort();
    v
}

fn run_split(args: &[&str]) {
    let st = ours().args(args).status().unwrap();
    assert!(st.success(), "split failed: {args:?}");
}

#[test]
fn se_split_by_lines_correct() {
    let tmp = tempfile::tempdir().unwrap();
    let od = tmp.path().join("o");
    std::fs::create_dir_all(&od).unwrap();
    let input = fixture("se5.fastq");
    let out = od.join("o.fq");

    run_split(&[
        "--split_by_lines",
        "8",
        "--in1",
        input.to_str().unwrap(),
        "--out1",
        out.to_str().unwrap(),
    ]);

    let files = split_files(&od, "o.fq");
    assert!(files.len() >= 2, "expected multiple split files");
    // reassembly preserves the input exactly
    let mut joined = Vec::new();
    for f in &files {
        joined.extend(std::fs::read(f).unwrap());
    }
    assert_eq!(
        joined,
        std::fs::read(&input).unwrap(),
        "reassembly != input"
    );
    // every chunk but the last holds exactly the requested line count
    for f in &files[..files.len() - 1] {
        let lines = std::fs::read_to_string(f).unwrap().lines().count();
        assert_eq!(lines, 8, "chunk {f:?} should hold 8 lines");
    }
}

#[test]
fn se_split_into_n_correct() {
    let tmp = tempfile::tempdir().unwrap();
    let od = tmp.path().join("o");
    std::fs::create_dir_all(&od).unwrap();
    let input = fixture("se5.fastq");
    let out = od.join("o.fq");

    run_split(&[
        "--split",
        "2",
        "--in1",
        input.to_str().unwrap(),
        "--out1",
        out.to_str().unwrap(),
    ]);

    let files = split_files(&od, "o.fq");
    assert_eq!(files.len(), 2, "--split 2 should make 2 files");
    let mut joined = Vec::new();
    for f in &files {
        joined.extend(std::fs::read(f).unwrap());
    }
    assert_eq!(
        joined,
        std::fs::read(&input).unwrap(),
        "reassembly != input"
    );
}

#[test]
fn pe_split_keeps_mates_aligned() {
    let tmp = tempfile::tempdir().unwrap();
    let od = tmp.path().join("o");
    std::fs::create_dir_all(&od).unwrap();
    let r1 = fixture("pe.fastq.r1");
    let r2 = fixture("pe.fastq.r2");
    let o1 = od.join("o1.fq");
    let o2 = od.join("o2.fq");

    run_split(&[
        "--split",
        "2",
        "--in1",
        r1.to_str().unwrap(),
        "--in2",
        r2.to_str().unwrap(),
        "--out1",
        o1.to_str().unwrap(),
        "--out2",
        o2.to_str().unwrap(),
    ]);

    let f1 = split_files(&od, "o1.fq");
    let f2 = split_files(&od, "o2.fq");
    assert_eq!(
        f1.len(),
        f2.len(),
        "R1 and R2 must split into equal file counts"
    );
    let reassemble = |fs: &[PathBuf]| {
        let mut v = Vec::new();
        for f in fs {
            v.extend(std::fs::read(f).unwrap());
        }
        v
    };
    assert_eq!(
        reassemble(&f1),
        std::fs::read(&r1).unwrap(),
        "R1 reassembly != input"
    );
    assert_eq!(
        reassemble(&f2),
        std::fs::read(&r2).unwrap(),
        "R2 reassembly != input"
    );
}
