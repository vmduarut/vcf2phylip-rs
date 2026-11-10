use std::{fs, process::Command};

#[test]
fn documented_example_matches_hand_computed_outputs_in_both_backends() {
    let root = tempfile::tempdir().unwrap();
    let input = root.path().join("example.vcf");
    fs::write(&input, include_str!("../examples/example.vcf")).unwrap();
    for storage in ["memory", "disk"] {
        let output = root.path().join(storage);
        let mut command = Command::new(env!("CARGO_BIN_EXE_vcf2phylip-rs"));
        command
            .arg("-i")
            .arg(&input)
            .args(["--storage", storage, "-f", "-w", "--output-folder"])
            .arg(&output);
        if storage == "disk" {
            command.args(["--memory-limit", "4B", "--temp-dir"]);
            command.arg(root.path());
        }
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        for (suffix, expected) in [
            ("phy", include_str!("../examples/expected.phy")),
            ("fasta", include_str!("../examples/expected.fasta")),
            (
                "used_sites.tsv",
                include_str!("../examples/expected.used_sites.tsv"),
            ),
        ] {
            assert_eq!(
                fs::read_to_string(output.join(format!("example.min4.{suffix}"))).unwrap(),
                expected,
                "{storage}: {suffix}"
            );
        }
    }
}
