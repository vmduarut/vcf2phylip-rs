use std::{fs, process::Command};

use clap::Parser;
use vcf2phylip_rs::{cli::Cli, conversion::convert};

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
fn effective_threshold_combines_rate_and_absolute_count() {
    for (samples, flags, expected) in [
        (3, vec![], 3),
        (3, vec!["-m", "2"], 3), // Explicit -m does not disable the default rate.
        (3, vec!["--min-call-rate", "0.66"], 2),
        (3, vec!["--min-call-rate", "0.66", "-m", "3"], 3),
        (3, vec!["--min-call-rate", "0", "-m", "2"], 2),
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
}

#[test]
fn invalid_rates_rejected_by_parser_and_public_validation() {
    for rate in [
        "NaN", "inf", "-inf", "Infinity", "-0.01", "1.01", "1e309", "garbage",
    ] {
        assert!(
            Cli::try_parse_from(["app", "-i", "x", &format!("--min-call-rate={rate}")]).is_err(),
            "{rate}"
        );
    }
    let mut cli = Cli::try_parse_from(["app", "-i", "x"]).unwrap();
    for rate in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.01, 1.01] {
        cli.min_call_rate = rate;
        assert!(cli.validate_mode().is_err());
    }
    assert!(cli.validate_samples(&[]).is_err());
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
    )).unwrap();
    for (flags, minimum, dna, binary, positions) in [
        (vec![], 3, 1, 1, vec![1]),
        (vec!["-m", "2"], 3, 1, 1, vec![1]),
        (vec!["--min-call-rate", "0.66"], 2, 3, 2, vec![1, 2, 3]),
        (
            vec!["--min-call-rate", "0", "-m", "2"],
            2,
            3,
            2,
            vec![1, 2, 3],
        ),
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
