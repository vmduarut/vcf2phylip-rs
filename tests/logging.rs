use std::{fs, process::Command};

use clap::Parser;
use vcf2phylip_rs::{cli::Cli, conversion::convert_with_observer};

const HEADER: &str = "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n";

fn fixture() -> String {
    let mut text = HEADER.to_owned();
    for (i, (reference, alt, a, b)) in [
        ("A", "C", "0/1", "1/1"),
        ("A", "C,G,T", "0/1/2/3", "0/0"), // called N, not missing
        ("A", "AT", "0/0", "0/1"),
        ("A", ".", "0/0", "0/0"),
        ("A", "A", "0/0", "0/1"),
        ("A", "C", ".", "./."),
        ("A", "C", "0/.", "0/0"),
    ]
    .into_iter()
    .enumerate()
    {
        text.push_str(&format!(
            "chr1\t{}\t.\t{reference}\t{alt}\t.\tPASS\t.\tGT\t{a}\t{b}\n",
            i + 1
        ));
    }
    text
}

fn run(dir: &std::path::Path, flags: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_vcf2phylip-rs"))
        .args(["-i", "input.vcf"])
        .args(flags)
        .current_dir(dir)
        .output()
        .unwrap()
}

#[test]
fn conversion_levels_streams_and_accurate_partition() {
    for (prefix, flags) in [
        ("normal", vec![]),
        ("quiet", vec!["--quiet"]),
        ("verbose", vec!["--verbose"]),
    ] {
        let dir = tempfile::tempdir().unwrap();
        fs::write(dir.path().join("input.vcf"), fixture()).unwrap();
        let mut flags = flags;
        flags.extend(["-m", "4"]); // Explicitly exercise upstream count-cap warning.
        let output = run(dir.path(), &flags);
        assert!(output.status.success(), "{output:?}");
        assert!(output.stdout.is_empty());
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert_eq!(stderr.matches("warning:").count(), 1, "{stderr}");
        assert!(stderr.contains("capped from 4 to 2"));
        assert_eq!(stderr.contains("Configuration:"), prefix == "verbose");
        if prefix == "quiet" {
            assert_eq!(stderr.lines().count(), 1);
        } else {
            let (body, elapsed) = stderr.split_once("Conversion elapsed time: ").unwrap();
            assert!(elapsed.ends_with("s\nDone.\n"));
            assert!(
                elapsed
                    .trim_end_matches("s\nDone.\n")
                    .parse::<f64>()
                    .is_ok()
            );
            let verbose = if prefix == "verbose" {
                "Configuration:\nStorage: Memory\nDisk cell-buffer budget: 67108864 bytes (disk only)\nPartial missing policy: Missing\nRandom seed: None\nNumber of outgroups: 0\n"
            } else {
                ""
            };
            let missing = if prefix == "verbose" {
                "Nucleotide missing sample calls at SNP-eligible records: 3 (not rejected loci).\n"
            } else {
                ""
            };
            assert_eq!(
                body,
                format!(
                    concat!(
                        "Converting file 'input.vcf':\n",
                        "warning: minimum samples capped from 4 to 2 (input sample count)\n",
                        "Number of samples in input: 2\n",
                        "Configured absolute minimum called samples: 4\nMinimum call rate: 0.9 (90.00%)\n",
                        "Minimum called samples per locus: 2 (effective)\n",
                        "Input format: Vcf\nInput gzip/BGZF compressed: false\n",
                        "{}Total records processed: 7\n",
                        "Records excluded as non-SNP: 3 (includes absent ALT, noncanonical and duplicate alleles)\n",
                        "SNP loci excluded for insufficient called samples: 2\nSNP loci retained: 2\n",
                        "{}Matrix payload: 4 bytes (Memory, not peak RSS)\n",
                        "PHYLIP matrix saved to ./input.min2.phy\n"
                    ),
                    verbose, missing
                )
            );
            assert!(!stderr.contains(';'));
        }
        assert_eq!(
            fs::read_to_string(dir.path().join("input.min2.phy")).unwrap(),
            "2 2\nA  MN\nB  CA\n"
        );
        // Repeat fails preflight, including in quiet mode, and never claims success.
        let failed = run(dir.path(), &flags);
        assert!(!failed.status.success());
        assert!(failed.stdout.is_empty());
        let error = String::from_utf8(failed.stderr).unwrap();
        assert!(error.contains("existing output"));
        assert!(!error.contains("Done."));
        assert!(!error.contains("saved to"));
        assert_eq!(error.matches("warning:").count(), 1);
    }
}

#[test]
fn quiet_without_warning_is_silent_and_inspection_stdout_is_unchanged() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("input.vcf"), fixture()).unwrap();
    let converted = run(dir.path(), &["--quiet", "-m", "2"]);
    assert!(converted.status.success());
    assert!(converted.stdout.is_empty());
    assert!(converted.stderr.is_empty());
    let baseline = run(dir.path(), &["--inspect", "-m", "2"]);
    assert!(baseline.status.success());
    for level in ["--quiet", "--verbose"] {
        let output = run(dir.path(), &["--inspect", "-m", "2", level]);
        assert!(output.status.success());
        assert_eq!(output.stdout, baseline.stdout);
        if level == "--quiet" {
            assert!(output.stderr.is_empty());
        } else {
            let stderr = String::from_utf8(output.stderr).unwrap();
            assert!(stderr.contains("Configuration:"));
            assert!(stderr.contains("Nucleotide missing sample calls at SNP-eligible records: 3"));
            assert!(stderr.contains("Inspection elapsed time:"));
            assert!(stderr.ends_with("Done.\n"));
            assert!(!stderr.contains(';'));
        }
    }
    let capped = run(dir.path(), &["--inspect", "--quiet", "-m", "4"]);
    assert!(capped.status.success());
    let stderr = String::from_utf8(capped.stderr).unwrap();
    assert_eq!(stderr.lines().count(), 1);
    assert_eq!(stderr.matches("warning:").count(), 1);
}

#[test]
fn output_labels_match_selected_formats() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(dir.path().join("input.vcf"), fixture()).unwrap();
    let output = run(dir.path(), &["-m", "2", "--fasta", "--write-used-sites"]);
    assert!(output.status.success(), "{output:?}");
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains(concat!(
        "PHYLIP matrix saved to ./input.min2.phy\n",
        "FASTA matrix saved to ./input.min2.fasta\n",
        "Used sites saved to ./input.min2.used_sites.tsv\n"
    )));
    assert!(!stderr.contains(';'));
}

#[test]
fn conflicts_and_errors_never_claim_success_or_paths() {
    let dir = tempfile::tempdir().unwrap();
    for flags in [vec!["--quiet", "--verbose"], vec!["--verbose", "--quiet"]] {
        let output = run(dir.path(), &flags);
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        assert!(
            String::from_utf8(output.stderr)
                .unwrap()
                .contains("cannot be used with")
        );
    }
    for text in [
        None,
        Some(HEADER.to_owned()),
        Some(format!("{}bad late record\n", fixture())),
    ] {
        if let Some(text) = text {
            fs::write(dir.path().join("input.vcf"), text).unwrap();
        }
        for level in ["--quiet", "--verbose"] {
            let output = run(dir.path(), &[level]);
            assert!(!output.status.success());
            assert!(output.stdout.is_empty());
            let stderr = String::from_utf8(output.stderr).unwrap();
            assert!(stderr.contains("Error:"));
            assert!(!stderr.contains("Done."));
            assert!(!stderr.contains("saved to"));
            assert!(!dir.path().join("input.min2.phy").exists());
        }
    }
}

#[test]
fn observer_reuses_header_once_and_stats_partition_records() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.vcf");
    fs::write(&input, fixture()).unwrap();
    let cli = Cli::try_parse_from([
        "app",
        "-i",
        input.to_str().unwrap(),
        "--output-folder",
        dir.path().to_str().unwrap(),
    ])
    .unwrap();
    let mut calls = 0;
    let report = convert_with_observer(&cli, |reader, config| {
        calls += 1;
        assert_eq!(reader.records_read(), 0);
        assert_eq!(reader.header().sample_names().len(), 2);
        assert!(config.warnings.is_empty());
    })
    .unwrap();
    assert_eq!(calls, 1);
    let stats = report.stats;
    assert_eq!(stats.nucleotide_below_minimum(), 2);
    assert_eq!(stats.nucleotide_missing_calls, 3);
    assert_eq!(
        stats.records,
        stats.rejected_non_snp + stats.nucleotide_below_minimum() + stats.nucleotide_accepted
    );
    assert!(report.outputs.iter().all(|path| path.is_file()));
}

#[test]
fn public_struct_literals_remain_compatible() {
    use vcf2phylip_rs::{filter::Stats, genotype::PartialMissingPolicy, storage::StorageMode};

    // Deliberately no struct update syntax: track the public shape, including
    // the new rate field; logging itself still adds no required library fields.
    let cli = Cli {
        input: "input.vcf".into(),
        inspect: false,
        storage: StorageMode::Memory,
        memory_limit: None,
        temp_dir: None,
        min_samples_locus: 4,
        min_call_rate: 0.9,
        outgroup: vec![],
        output_folder: None,
        output_prefix: None,
        force: false,
        phylip_disable: false,
        fasta: false,
        resolve_iupac: false,
        seed: None,
        partial_missing: PartialMissingPolicy::Missing,
        write_used_sites: false,
        version: (),
    };
    assert!(cli.validate_mode().is_ok());
    let stats = Stats {
        records: 7,
        rejected_non_snp: 3,
        nucleotide_accepted: 2,
        binary_eligible: 0,
        binary_accepted: 0,
        nucleotide_missing_calls: 0,
        binary_missing_calls: 0,
        binary_unsupported_ploidy_calls: 0,
    };
    assert_eq!(stats.nucleotide_below_minimum(), 2);
    assert_eq!(Stats::default().nucleotide_below_minimum(), 0);
}

#[test]
fn binary_wrapper_preserves_help_and_version() {
    for flag in ["--help", "--version", "-v"] {
        let output = Command::new(env!("CARGO_BIN_EXE_vcf2phylip-rs"))
            .arg(flag)
            .output()
            .unwrap();
        assert!(output.status.success(), "{output:?}");
        assert!(output.stderr.is_empty());
        let stdout = String::from_utf8(output.stdout).unwrap();
        if flag == "--help" {
            for option in ["--quiet", "--verbose", "--inspect", "--version", "--input"] {
                assert!(stdout.contains(option), "{stdout}");
            }
        } else {
            assert_eq!(
                stdout.trim(),
                concat!("vcf2phylip-rs ", env!("CARGO_PKG_VERSION"))
            );
        }
    }
}
