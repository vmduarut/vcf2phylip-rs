use clap::{CommandFactory, Parser, error::ErrorKind};
use vcf2phylip_rs::cli::Cli;

#[test]
fn nexus_options_are_not_in_help() {
    let help = Cli::command().render_long_help().to_string();
    for option in ["-n", "--nexus", "-b", "--nexus-binary"] {
        assert!(
            !help
                .split_whitespace()
                .any(|word| word.trim_end_matches(',') == option),
            "help contains {option}"
        );
    }
}

#[test]
fn nexus_options_are_unknown() {
    for option in ["-n", "--nexus", "-b", "--nexus-binary"] {
        let error = Cli::try_parse_from(["vcf2phylip", "-i", "input.vcf", option])
            .expect_err("removed option should be rejected");
        assert_eq!(error.kind(), ErrorKind::UnknownArgument, "option {option}");
    }
}

#[test]
fn threshold_defaults_and_help_preserve_upstream_absolute_name() {
    let cli = Cli::try_parse_from(["app", "-i", "x"]).unwrap();
    assert_eq!(cli.min_samples_locus, None);
    assert_eq!(cli.min_call_rate, None);
    let help = Cli::command().render_long_help().to_string();
    for text in [
        "--min-samples-locus",
        "--min-call-rate",
        "Defaults to 0.9 (90%) when neither threshold option is supplied",
        "cannot be combined with --min-call-rate",
        "Cannot be combined with -m/--min-samples-locus",
        "ceil(rate * input samples)",
    ] {
        assert!(help.contains(text), "{help}");
    }
    assert!(!help.contains("--min-called-samples"));
    assert!(!help.contains("[default: 1]"));
    let cli = Cli::try_parse_from(["app", "-i", "x", "--min-samples-locus", "1"]).unwrap();
    assert_eq!(cli.min_samples_locus, Some(1));
}
