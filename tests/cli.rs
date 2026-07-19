//! End-to-end CLI tests that drive the built `snpick` binary. Cargo sets
//! `CARGO_BIN_EXE_snpick` for integration tests, so these exercise argument
//! parsing and output routing exactly as a user would.

use std::fs;
use std::path::PathBuf;
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_snpick");

/// A unique temp directory per test, removed on drop.
struct TempDir(PathBuf);
impl TempDir {
    fn new(tag: &str) -> Self {
        let p = std::env::temp_dir().join(format!("snpick-cli-{}-{}", tag, std::process::id()));
        let _ = fs::remove_dir_all(&p);
        fs::create_dir_all(&p).unwrap();
        TempDir(p)
    }
    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
    fn write(&self, name: &str, contents: &str) -> PathBuf {
        let p = self.path(name);
        fs::write(&p, contents).unwrap();
        p
    }
}
impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

const ALN: &str = ">ref\nATGCAT\n>s1\nATGTAT\n>s2\nACGCAT\n>s3\nATGCAG\n";

#[test]
fn vcf_output_dash_streams_to_stdout() {
    let d = TempDir::new("vcfdash");
    let fa = d.write("aln.fa", ALN);
    let out = d.path("out.fa");
    let o = Command::new(BIN)
        .current_dir(&d.0)
        .args(["-f", fa.to_str().unwrap(), "-o", out.to_str().unwrap(), "--vcf-output", "-", "-q"])
        .output()
        .unwrap();
    assert!(o.status.success(), "stderr: {}", String::from_utf8_lossy(&o.stderr));
    // The VCF must arrive on stdout, not in a file literally named "-".
    let stdout = String::from_utf8_lossy(&o.stdout);
    assert!(stdout.starts_with("##fileformat=VCFv4.2"), "stdout was: {:?}", stdout);
    assert!(!d.path("-").exists(), "a file named '-' should NOT have been created");
    // The reduced FASTA still went to its file.
    assert!(out.exists());
}

#[test]
fn check_does_not_reject_output_path() {
    // --check writes nothing, so pointing -o at the input must not trip the collision guard,
    // matching --dry-run's behavior.
    let d = TempDir::new("checkpath");
    let fa = d.write("aln.fa", ALN);
    for mode in ["--check", "--dry-run"] {
        let o = Command::new(BIN)
            .args(["-f", fa.to_str().unwrap(), mode, "-o", fa.to_str().unwrap(), "-q"])
            .output()
            .unwrap();
        assert!(o.status.success(), "{mode} -o <input> should succeed; stderr: {}",
            String::from_utf8_lossy(&o.stderr));
    }
}

#[test]
fn multiple_stdout_outputs_rejected() {
    // Two sidecars both aimed at stdout would concatenate JSON+TSV into one unparseable stream.
    let d = TempDir::new("multistdout");
    let fa = d.write("aln.fa", ALN);
    let out = d.path("out.fa");
    let o = Command::new(BIN)
        .args(["-f", fa.to_str().unwrap(), "-o", out.to_str().unwrap(),
               "--stats-json", "-", "--sites-output", "-", "-q"])
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(2), "stdout: {}", String::from_utf8_lossy(&o.stdout));
    assert!(String::from_utf8_lossy(&o.stderr).contains("stdout"));
    // Exactly one sidecar to stdout is still fine.
    let ok = Command::new(BIN)
        .args(["-f", fa.to_str().unwrap(), "-o", out.to_str().unwrap(), "--stats-json", "-", "-q"])
        .output()
        .unwrap();
    assert!(ok.status.success());
    assert!(String::from_utf8_lossy(&ok.stdout).starts_with('{'));
}

#[test]
fn vcf_derive_rejects_stdout_alignment() {
    // `-o - --vcf` cannot derive a sensible VCF filename; it must error, not write "-.vcf".
    let d = TempDir::new("vcfderive");
    let fa = d.write("aln.fa", ALN);
    let o = Command::new(BIN)
        .current_dir(&d.0)
        .args(["-f", fa.to_str().unwrap(), "-o", "-", "--vcf", "-q"])
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(2), "stderr: {}", String::from_utf8_lossy(&o.stderr));
    assert!(!d.path("-.vcf").exists(), "no '-.vcf' file should be created");
}
