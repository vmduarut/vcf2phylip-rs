use std::{fs, process::Command};

const VCF: &str = "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\nchr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n";

fn assert_safe(output: &std::process::Output) {
    for stream in [&output.stdout, &output.stderr] {
        assert!(!stream.iter().any(|&b| b == 0x1b || b == b'\r'));
    }
}

#[test]
fn missing_path_is_escaped_even_when_quiet() {
    for quiet in [false, true] {
        let dir = tempfile::tempdir().unwrap();
        let mut command = Command::new(env!("CARGO_BIN_EXE_vcf2phylip-rs"));
        command.current_dir(dir.path()).args([
            "-i",
            "bad\x1b\r\nname.vcf",
            "--output-prefix",
            "safe",
        ]);
        if quiet {
            command.arg("--quiet");
        }
        let output = command.output().unwrap();
        assert!(!output.status.success());
        assert_safe(&output);
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert!(stderr.contains("bad\\u{1b}\\r\\nname.vcf"), "{stderr}");
        assert!(stderr.contains("Error:"));
        assert!(!stderr.contains("\nname.vcf"));
        if quiet {
            assert!(!stderr.contains("Converting"));
        }
    }
}

#[test]
fn successful_path_logs_and_inspection_identifier_errors_are_safe() {
    let dir = tempfile::tempdir().unwrap();
    let name = "sample\x1b\r\nname.vcf";
    fs::write(dir.path().join(name), VCF).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_vcf2phylip-rs"))
        .current_dir(dir.path())
        .args([
            "-i",
            name,
            "-m",
            "1",
            "--output-prefix",
            "safe",
            "--output-folder",
            "out\x1b\r\ndir",
        ])
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    assert_safe(&output);
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("sample\\u{1b}\\r\\nname"), "{stderr}");
    assert!(dir.path().join("out\x1b\r\ndir/safe.min1.phy").exists());

    let output = Command::new(env!("CARGO_BIN_EXE_vcf2phylip-rs"))
        .current_dir(dir.path())
        .args([
            "-i",
            name,
            "--inspect",
            "--outgroup",
            "bad\x1b\r\nname",
            "--quiet",
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_safe(&output);
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("bad\\u{1b}\\r\\nname"), "{stderr}");
}

#[test]
fn clap_errors_keep_layout_and_escape_untrusted_values() {
    for arg in ["--bogus", "--bad\x1b\r\nINJECTED"] {
        let output = Command::new(env!("CARGO_BIN_EXE_vcf2phylip-rs"))
            .arg(arg)
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(2));
        assert_safe(&output);
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert!(stderr.contains("\n\nUsage:"), "{stderr}");
        assert!(!stderr.contains("\\n\\nUsage:"), "{stderr}");
        if arg != "--bogus" {
            assert!(stderr.contains("--bad\\u{1b}\\r\\nINJECTED"), "{stderr}");
            assert!(!stderr.contains("\nINJECTED"), "{stderr}");
        }
    }
}
