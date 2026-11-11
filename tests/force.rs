use std::{fs, path::Path, process::Command};

const INPUT: &str = include_str!("../examples/example.vcf");
const OUTPUTS: [(&str, &str); 3] = [
    ("phy", include_str!("../examples/expected.phy")),
    ("fasta", include_str!("../examples/expected.fasta")),
    (
        "used_sites.tsv",
        include_str!("../examples/expected.used_sites.tsv"),
    ),
];

#[test]
fn help_documents_long_force_option() {
    let result = Command::new(env!("CARGO_BIN_EXE_vcf2phylip-rs"))
        .arg("--help")
        .output()
        .unwrap();
    assert!(result.status.success());
    let help = String::from_utf8(result.stdout).unwrap();
    assert!(help.contains("--force"));
    assert!(help.contains("Replace existing regular output files"));
}

fn command(input: &Path, folder: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_vcf2phylip-rs"));
    command
        .arg("-i")
        .arg(input)
        .args(["-f", "-w", "--output-prefix", "example", "--output-folder"])
        .arg(folder);
    command
}

fn seed_outputs(folder: &Path) {
    for (suffix, _) in OUTPUTS {
        fs::write(folder.join(format!("example.min4.{suffix}")), suffix).unwrap();
    }
}

fn assert_originals(folder: &Path) {
    for (suffix, _) in OUTPUTS {
        assert_eq!(
            fs::read_to_string(folder.join(format!("example.min4.{suffix}"))).unwrap(),
            suffix
        );
    }
}

#[test]
fn default_refuses_but_force_replaces_all_selected_outputs_exactly() {
    for storage in ["memory", "disk"] {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("input.vcf");
        fs::write(&input, INPUT).unwrap();
        seed_outputs(dir.path());
        let unrelated = dir.path().join("unselected.min4.fasta");
        fs::write(&unrelated, "untouched").unwrap();
        let result = command(&input, dir.path()).output().unwrap();
        assert!(!result.status.success());
        assert!(String::from_utf8_lossy(&result.stderr).contains("existing output"));
        assert_originals(dir.path());

        let result = command(&input, dir.path())
            .args(["--force", "--storage", storage, "--verbose"])
            .output()
            .unwrap();
        assert!(result.status.success(), "{result:?}");
        assert!(
            String::from_utf8_lossy(&result.stderr).contains("Replace existing output files: true")
        );
        for (suffix, expected) in OUTPUTS {
            assert_eq!(
                fs::read_to_string(dir.path().join(format!("example.min4.{suffix}"))).unwrap(),
                expected
            );
        }
        assert_eq!(fs::read_to_string(unrelated).unwrap(), "untouched");
        assert_eq!(fs::read_to_string(input).unwrap(), INPUT);
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 5);
    }
}

#[test]
fn force_creates_absent_outputs_and_does_not_replace_unselected_formats() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.vcf");
    fs::write(&input, INPUT).unwrap();
    let result = command(&input, dir.path()).arg("--force").output().unwrap();
    assert!(result.status.success(), "{result:?}");
    for (suffix, expected) in OUTPUTS {
        assert_eq!(
            fs::read_to_string(dir.path().join(format!("example.min4.{suffix}"))).unwrap(),
            expected
        );
    }
    seed_outputs(dir.path());
    let result = Command::new(env!("CARGO_BIN_EXE_vcf2phylip-rs"))
        .arg("-i")
        .arg(&input)
        .args(["--force", "--output-prefix", "example", "--output-folder"])
        .arg(dir.path())
        .output()
        .unwrap();
    assert!(result.status.success(), "{result:?}");
    for suffix in ["fasta", "used_sites.tsv"] {
        assert_eq!(
            fs::read_to_string(dir.path().join(format!("example.min4.{suffix}"))).unwrap(),
            suffix
        );
    }
    assert_eq!(
        fs::read_to_string(dir.path().join("example.min4.phy")).unwrap(),
        OUTPUTS[0].1
    );
}

#[test]
fn force_preserves_outputs_on_late_parse_zero_loci_and_validation_errors() {
    let header = INPUT
        .lines()
        .filter(|line| line.starts_with('#'))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    for invalid in [format!("{INPUT}broken\n"), header] {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("input.vcf");
        fs::write(&input, invalid).unwrap();
        seed_outputs(dir.path());
        let result = command(&input, dir.path()).arg("--force").output().unwrap();
        assert!(!result.status.success());
        assert_originals(dir.path());
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 4);
    }
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.vcf");
    fs::write(&input, INPUT).unwrap();
    seed_outputs(dir.path());
    let result = command(&input, dir.path())
        .args(["--force", "--outgroup", "not-a-sample"])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_originals(dir.path());
}

#[test]
fn force_refuses_directory_and_hardlinked_input_before_changing_any_outputs() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.vcf");
    fs::write(&input, INPUT).unwrap();
    seed_outputs(dir.path());
    let target = dir.path().join("example.min4.fasta");
    fs::remove_file(&target).unwrap();
    fs::create_dir(&target).unwrap();
    let result = command(&input, dir.path()).arg("--force").output().unwrap();
    assert!(!result.status.success());
    assert!(target.is_dir());
    assert_eq!(
        fs::read_to_string(dir.path().join("example.min4.phy")).unwrap(),
        "phy"
    );
    fs::remove_dir(&target).unwrap();
    fs::hard_link(&input, &target).unwrap();
    let result = command(&input, dir.path()).arg("--force").output().unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("overwrite input"));
    assert_eq!(fs::read_to_string(&input).unwrap(), INPUT);
    assert_eq!(fs::read_to_string(&target).unwrap(), INPUT);
    assert_eq!(
        fs::read_to_string(dir.path().join("example.min4.phy")).unwrap(),
        "phy"
    );
}

#[test]
fn force_refuses_direct_input_target_and_inspection() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("example.min4.phy");
    fs::write(&input, INPUT).unwrap();
    let result = command(&input, dir.path()).arg("--force").output().unwrap();
    assert!(!result.status.success());
    assert_eq!(fs::read_to_string(input).unwrap(), INPUT);
    let result = Command::new(env!("CARGO_BIN_EXE_vcf2phylip-rs"))
        .args(["-i", "does-not-exist", "--inspect", "--force"])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("cannot be used with --inspect"));
}

#[cfg(unix)]
#[test]
fn force_preserves_existing_regular_file_permissions() {
    use std::os::unix::fs::PermissionsExt;

    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.vcf");
    fs::write(&input, INPUT).unwrap();
    seed_outputs(dir.path());
    let modes = [0o644, 0o640, 0o600];
    for ((suffix, _), mode) in OUTPUTS.into_iter().zip(modes) {
        fs::set_permissions(
            dir.path().join(format!("example.min4.{suffix}")),
            fs::Permissions::from_mode(mode),
        )
        .unwrap();
    }
    let result = command(&input, dir.path()).arg("--force").output().unwrap();
    assert!(result.status.success(), "{result:?}");
    for ((suffix, expected), mode) in OUTPUTS.into_iter().zip(modes) {
        let target = dir.path().join(format!("example.min4.{suffix}"));
        assert_eq!(
            fs::metadata(&target).unwrap().permissions().mode() & 0o7777,
            mode
        );
        assert_eq!(fs::read_to_string(target).unwrap(), expected);
    }
}

#[cfg(unix)]
#[test]
fn force_refuses_symlinks_including_dangling_and_input_aliases() {
    use std::os::unix::fs::symlink;
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.vcf");
    let outside = dir.path().join("outside");
    fs::write(&input, INPUT).unwrap();
    fs::write(&outside, "outside contents").unwrap();
    let target = dir.path().join("example.min4.phy");
    for destination in [&input, &outside, &dir.path().join("absent")] {
        symlink(destination, &target).unwrap();
        let result = command(&input, dir.path()).arg("--force").output().unwrap();
        assert!(!result.status.success());
        assert!(fs::symlink_metadata(&target).unwrap().is_symlink());
        assert_eq!(fs::read_to_string(&input).unwrap(), INPUT);
        assert_eq!(fs::read_to_string(&outside).unwrap(), "outside contents");
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 3);
        fs::remove_file(&target).unwrap();
    }
    // An input symlink pointing at an otherwise regular output is also an alias.
    fs::write(&target, INPUT).unwrap();
    let input_alias = dir.path().join("alias.vcf");
    symlink(&target, &input_alias).unwrap();
    let result = command(&input_alias, dir.path())
        .arg("--force")
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("overwrite input"));
    assert_eq!(fs::read_to_string(target).unwrap(), INPUT);
}
