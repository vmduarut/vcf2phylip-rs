use std::{fs, process::Command};

use clap::{Parser, error::ErrorKind};
use vcf2phylip_rs::{
    cli::{Cli, ConfigError},
    conversion::convert,
};

fn configured(samples: usize, flags: &[&str]) -> vcf2phylip_rs::cli::ValidatedConfig {
    let cli = Cli::try_parse_from(
        ["app", "-i", "unused"]
            .into_iter()
            .chain(flags.iter().copied()),
    )
    .unwrap();
    let names = (0..samples).map(|i| format!("S{i}")).collect::<Vec<_>>();
    let config = cli.validate_samples(&names).unwrap();
    assert_eq!(config.minimum_samples, config.processing.minimum_samples());
    config
}

#[test]
fn selected_count_or_rate_determines_threshold() {
    for (samples, flags, expected) in [
        (3, vec![], 3),
        (3, vec!["-m", "1"], 1),
        (3, vec!["-m", "2"], 2),
        (3, vec!["--min-samples-locus", "1"], 1),
        (3, vec!["--min-call-rate", "0.66"], 2),
        (3, vec!["--min-call-rate", "0.9"], 3),
        (3, vec!["--min-call-rate", "0"], 1),
        (3, vec!["--min-call-rate", "1"], 3),
        (100, vec!["--min-call-rate", "0.07"], 7),
        (100, vec!["--min-call-rate", "0.07000000000000002"], 8),
        (100, vec!["--min-call-rate", "7e-2"], 7),
        (100, vec!["--min-call-rate", "0.06999999999999999"], 7),
        (3, vec!["--min-call-rate", "5e-324"], 1),
        (3, vec!["--min-call-rate=-0"], 1),
    ] {
        assert_eq!(
            configured(samples, &flags).minimum_samples,
            expected,
            "{flags:?}"
        );
    }
    let capped = configured(3, &["-m", "4"]);
    assert_eq!(capped.minimum_samples, 3);
    assert_eq!(
        capped.warnings,
        ["minimum samples capped from 4 to 3 (input sample count)"]
    );
    assert!(configured(3, &[]).warnings.is_empty());
    assert!(configured(3, &["-m", "1"]).warnings.is_empty());
}

#[test]
fn explicit_threshold_options_conflict_in_both_orders_and_all_modes() {
    let dir = tempfile::tempdir().unwrap();
    for mode in [
        vec!["--inspect"],
        vec!["--storage", "memory"],
        vec!["--storage", "disk"],
    ] {
        for count in ["-m", "--min-samples-locus"] {
            for rate in ["0", "0.9", "1"] {
                for flags in [
                    vec![count, "1", "--min-call-rate", rate],
                    vec!["--min-call-rate", rate, count, "1"],
                ] {
                    let error = Cli::try_parse_from(
                        ["app", "-i", "missing.vcf"]
                            .into_iter()
                            .chain(mode.iter().copied())
                            .chain(flags.iter().copied()),
                    )
                    .unwrap_err();
                    assert_eq!(
                        error.kind(),
                        ErrorKind::ArgumentConflict,
                        "{mode:?} {flags:?}"
                    );
                    let result = Command::new(env!("CARGO_BIN_EXE_vcf2phylip-rs"))
                        .args(["-i", "missing.vcf"])
                        .args(&mode)
                        .args(&flags)
                        .current_dir(dir.path())
                        .output()
                        .unwrap();
                    assert_eq!(result.status.code(), Some(2));
                    assert!(result.stdout.is_empty());
                    let stderr = String::from_utf8(result.stderr).unwrap();
                    assert!(stderr.contains("cannot be used with"), "{stderr}");
                }
            }
        }
        let mut cli = Cli::try_parse_from(
            ["app", "-i", "missing.vcf"]
                .into_iter()
                .chain(mode.iter().copied()),
        )
        .unwrap();
        for rate in [0.0, 0.9, 1.0] {
            cli.min_samples_locus = Some(1);
            cli.min_call_rate = Some(rate);
            assert!(matches!(
                cli.validate_mode(),
                Err(ConfigError::ConflictingThresholds)
            ));
            assert!(matches!(
                cli.validate_samples(&["A".to_owned()]),
                Err(ConfigError::ConflictingThresholds)
            ));
            let error = convert(&cli).unwrap_err();
            assert!(matches!(
                error.downcast_ref::<ConfigError>(),
                Some(ConfigError::ConflictingThresholds)
            ));
        }
    }
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
}

#[test]
fn invalid_rates_rejected_by_parser_and_public_validation() {
    for rate in [
        "NaN", "inf", "-inf", "Infinity", "-0.01", "1.01", "1e309", "garbage",
    ] {
        for count_flags in [vec![], vec!["-m", "1"]] {
            assert!(
                Cli::try_parse_from(
                    ["app", "-i", "x", &format!("--min-call-rate={rate}")]
                        .into_iter()
                        .chain(count_flags)
                )
                .is_err(),
                "{rate}"
            );
        }
    }
    let mut cli = Cli::try_parse_from(["app", "-i", "x"]).unwrap();
    for rate in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.01, 1.01] {
        cli.min_call_rate = Some(rate);
        cli.min_samples_locus = None;
        assert!(cli.validate_mode().is_err());
        cli.min_samples_locus = Some(1);
        assert!(cli.validate_mode().is_err());
    }
    cli.min_call_rate = None;
    cli.min_samples_locus = None;
    assert!(cli.validate_samples(&[]).is_err());
}

#[test]
fn invalid_counts_rejected_by_parser_and_public_validation() {
    for count in [
        "0",
        "-1",
        "1.5",
        "garbage",
        "999999999999999999999999999999999999999",
    ] {
        assert!(
            Cli::try_parse_from(["app", "-i", "x", &format!("--min-samples-locus={count}")])
                .is_err()
        );
    }
    let mut cli = Cli::try_parse_from(["app", "-i", "x"]).unwrap();
    cli.min_samples_locus = Some(0);
    assert!(cli.validate_mode().is_err());
    assert!(cli.validate_samples(&["A".to_owned()]).is_err());
}

#[test]
fn inspection_memory_disk_and_filenames_share_threshold_and_call_policy() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("cohort.vcf");
    fs::write(&input, concat!(
        "##fileformat=VCFv4.3\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n",
        "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n",
        "chr1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1\t./.\n",
        "chr1\t3\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/.\t1\n",
        "chr1\t4\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t./.\t./.\n",
        "chr1\t5\t.\tA\tC\t.\tPASS\t.\tGT\t./.\t./.\t./.\n",
        "chr1\t6\t.\tA\tAT\t.\tPASS\t.\tGT\t0/0\t./.\t./.\n",
    )).unwrap();
    for (flags, minimum, dna, binary, positions) in [
        (vec![], 3, 1, 1, vec![1]),
        (vec!["-m", "1"], 1, 4, 4, vec![1, 2, 3, 4]),
        (vec!["-m", "2"], 2, 3, 2, vec![1, 2, 3]),
        (vec!["--min-samples-locus", "1"], 1, 4, 4, vec![1, 2, 3, 4]),
        (vec!["--min-call-rate", "0.9"], 3, 1, 1, vec![1]),
        (vec!["-m", "4"], 3, 1, 1, vec![1]),
        (vec!["--min-call-rate", "0"], 1, 4, 4, vec![1, 2, 3, 4]),
        (vec!["--min-call-rate", "0.66"], 2, 3, 2, vec![1, 2, 3]),
        (
            vec![
                "--min-call-rate",
                "1",
                "--partial-missing",
                "called-alleles",
            ],
            3,
            2,
            1,
            vec![1, 3],
        ),
    ] {
        let mut expected = None;
        for storage in ["memory", "disk"] {
            let output = tempfile::tempdir().unwrap();
            let cli = Cli::try_parse_from(
                [
                    "app",
                    "-i",
                    input.to_str().unwrap(),
                    "--storage",
                    storage,
                    "--output-folder",
                    output.path().to_str().unwrap(),
                    "-f",
                    "-w",
                ]
                .into_iter()
                .chain(flags.iter().copied()),
            )
            .unwrap();
            let report = convert(&cli).unwrap();
            assert_eq!(report.stats.nucleotide_accepted, dna);
            assert_eq!(report.stats.binary_accepted, binary);
            let bytes = report
                .outputs
                .iter()
                .map(|path| {
                    assert!(
                        path.file_name()
                            .unwrap()
                            .to_str()
                            .unwrap()
                            .contains(&format!(".min{minimum}."))
                    );
                    fs::read(path).unwrap()
                })
                .collect::<Vec<_>>();
            let sites = String::from_utf8(bytes[2].clone()).unwrap();
            let actual_positions = sites
                .lines()
                .skip(1)
                .map(|line| line.split('\t').nth(1).unwrap().parse::<usize>().unwrap())
                .collect::<Vec<_>>();
            assert_eq!(actual_positions, positions);
            if let Some(expected) = &expected {
                assert_eq!(expected, &bytes);
            } else {
                expected = Some(bytes);
            }
        }
        let result = Command::new(env!("CARGO_BIN_EXE_vcf2phylip-rs"))
            .arg("-i")
            .arg(&input)
            .arg("--inspect")
            .args(&flags)
            .output()
            .unwrap();
        assert!(result.status.success(), "{result:?}");
        let stdout = String::from_utf8(result.stdout).unwrap();
        let count_policy = flags.contains(&"-m") || flags.contains(&"--min-samples-locus");
        assert_eq!(
            stdout.contains("Configured absolute minimum called samples:"),
            count_policy,
            "{stdout}"
        );
        assert_eq!(
            stdout.contains("Minimum call rate:"),
            !count_policy,
            "{stdout}"
        );
        assert!(!stdout.contains("ignored"), "{stdout}");
        assert!(
            stdout.contains(&format!("Minimum called samples: {minimum} (effective")),
            "{stdout}"
        );
        assert!(
            stdout.contains(&format!(
                "Accepted nucleotide loci: {dna}; accepted binary loci: {binary}"
            )),
            "{stdout}"
        );
    }
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
}
