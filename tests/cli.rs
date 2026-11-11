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
