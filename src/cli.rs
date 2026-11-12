//! Compatibility arguments and validation for inspection and conversion.
use std::{collections::HashSet, path::PathBuf};

use crate::{
    filter::ProcessingConfig,
    genotype::{PartialMissingPolicy, Resolution},
    storage::StorageMode,
};
use clap::{ArgAction, Parser};
use thiserror::Error;

/// Call-rate policy used when neither threshold option is supplied.
pub const DEFAULT_CALL_RATE: f64 = 0.9;

#[derive(Debug, Parser)]
#[command(
    version,
    disable_version_flag = true,
    about = "Convert VCF/BCF to FASTA and relaxed PHYLIP alignments"
)]
pub struct Cli {
    #[arg(short, long)]
    pub input: PathBuf,
    /// Stream and validate input; print diagnostic counts only, never an alignment.
    #[arg(long)]
    pub inspect: bool,
    /// Matrix backend; memory is unbounded, disk uses a bounded cell buffer.
    #[arg(long, value_enum, default_value = "memory")]
    pub storage: StorageMode,
    /// Disk cell-buffer budget (default 64 MiB). Bytes, KiB, MiB, GiB supported.
    #[arg(long, value_parser = parse_bytes)]
    pub memory_limit: Option<usize>,
    /// Existing directory for the disk matrix spool (default OS temporary dir).
    #[arg(long)]
    pub temp_dir: Option<PathBuf>,
    /// Absolute minimum called samples per locus (capped to input sample count).
    /// Count-only filtering; cannot be combined with --min-call-rate.
    #[arg(short = 'm', long, value_parser = positive_usize, conflicts_with = "min_call_rate")]
    pub min_samples_locus: Option<usize>,
    /// Minimum fraction of all input samples called per locus (finite 0..=1).
    /// Requires ceil(rate * input samples), independently for nucleotide and
    /// diagnostic binary calls, with at least one called sample required.
    /// Defaults to 0.9 (90%) when neither threshold option is supplied.
    /// Cannot be combined with -m/--min-samples-locus.
    #[arg(long, value_parser = parse_call_rate)]
    pub min_call_rate: Option<f64>,
    /// Exact sample name; repeat for multiple outgroups (no delimiter splitting).
    #[arg(short, long)]
    pub outgroup: Vec<String>,
    /// Output directory (defaults to the current directory).
    #[arg(long, default_value = ".", conflicts_with = "inspect")]
    pub output_folder: PathBuf,
    /// Output basename (defaults to input basename without known format suffixes).
    #[arg(long)]
    pub output_prefix: Option<String>,
    /// Replace existing regular output files after successful processing (per-file atomic).
    #[arg(long)]
    pub force: bool,
    /// Disable the default relaxed PHYLIP output.
    #[arg(short = 'p', long)]
    pub phylip_disable: bool,
    /// Also write FASTA (80 columns per sequence line).
    #[arg(short = 'f', long)]
    pub fasta: bool,
    /// Resolve nucleotide calls by sampling allele copies; requires --seed.
    #[arg(short = 'r', long = "resolve-IUPAC")]
    pub resolve_iupac: bool,
    /// Explicit reproducible seed for --resolve-IUPAC.
    #[arg(long)]
    pub seed: Option<u64>,
    /// Treatment of partially missing nucleotide genotypes (binary is always strict).
    #[arg(long, value_enum, default_value = "missing")]
    pub partial_missing: PartialMissingPolicy,
    /// Write retained nucleotide site coordinates and called sample counts.
    #[arg(short = 'w', long)]
    pub write_used_sites: bool,
    #[arg(short = 'v', long, action = ArgAction::Version)]
    pub version: (),
}

#[derive(Debug, Error)]
pub enum ConfigError {
    #[error("no alignment output selected; use --fasta with --phylip-disable")]
    NoOutputs,
    #[error("conversion/output flags cannot be used with --inspect")]
    OutputUnsupported,
    #[error("minimum samples must be positive")]
    InvalidMinimum,
    #[error("minimum call rate must be finite and between 0 and 1 inclusive")]
    InvalidCallRate,
    #[error(
        "--min-samples-locus and --min-call-rate cannot be used together; select one threshold policy"
    )]
    ConflictingThresholds,
    #[error("outgroup is not an exact input sample name: {0:?}")]
    UnknownOutgroup(String),
    #[error("outgroup specified more than once: {0:?}")]
    DuplicateOutgroup(String),
    #[error("input has no samples")]
    NoSamples,
    #[error("--resolve-IUPAC requires --seed; --seed requires --resolve-IUPAC")]
    InvalidSeed,
    #[error("--memory-limit and --temp-dir require --storage disk")]
    StorageOptions,
}

#[derive(Debug)]
pub struct ValidatedConfig {
    pub minimum_samples: usize,
    pub outgroup_indices: Vec<usize>,
    pub warnings: Vec<String>,
    pub processing: ProcessingConfig,
}

fn positive_usize(value: &str) -> Result<usize, String> {
    match value.parse::<usize>() {
        Ok(n) if n > 0 => Ok(n),
        _ => Err("minimum samples must be a positive integer".into()),
    }
}

fn parse_call_rate(value: &str) -> Result<f64, String> {
    value
        .parse::<f64>()
        .ok()
        .filter(|rate| rate.is_finite() && (0.0..=1.0).contains(rate))
        .ok_or_else(|| ConfigError::InvalidCallRate.to_string())
}

/// Multiply the validated rate's canonical shortest decimal representation by
/// the sample count with integer decimal arithmetic. This avoids binary-float
/// multiplication/ceil artifacts (e.g. 0.07 * 100 > 7), without an epsilon that
/// could erase a genuinely above-boundary decimal. u128 holds a usize * digit
/// plus carry on all supported (32/64-bit) platforms. Even subnormal rates work.
fn rate_minimum(rate: f64, samples: usize) -> usize {
    if rate == 1.0 {
        return samples;
    }
    let decimal = rate.to_string();
    let Some((_, fraction)) = decimal.split_once('.') else {
        return 0; // zero (including negative zero)
    };
    let mut carry = 0_u128;
    let mut remainder = false;
    for digit in fraction.bytes().rev() {
        let product = u128::from(digit - b'0') * samples as u128 + carry;
        remainder |= !product.is_multiple_of(10);
        carry = product / 10;
    }
    carry as usize + usize::from(remainder)
}

fn parse_bytes(value: &str) -> Result<usize, String> {
    let (digits, multiplier) = [
        ("KiB", 1024usize),
        ("MiB", 1024 * 1024),
        ("GiB", 1024 * 1024 * 1024),
        ("B", 1),
    ]
    .into_iter()
    .find_map(|(suffix, n)| value.strip_suffix(suffix).map(|s| (s, n)))
    .unwrap_or((value, 1));
    digits
        .parse::<usize>()
        .ok()
        .and_then(|n| n.checked_mul(multiplier))
        .filter(|&n| n > 0)
        .ok_or_else(|| {
            "expected positive bytes or integer KiB/MiB/GiB (no spaces), within platform usize"
                .into()
        })
}

impl Cli {
    pub fn validate_mode(&self) -> Result<(), ConfigError> {
        if self.storage == StorageMode::Memory
            && (self.memory_limit.is_some() || self.temp_dir.is_some())
        {
            return Err(ConfigError::StorageOptions);
        }
        if self.inspect
            && (self.storage != StorageMode::Memory
                || self.memory_limit.is_some()
                || self.temp_dir.is_some()
                || self.output_folder != std::path::Path::new(".")
                || self.output_prefix.is_some()
                || self.force
                || self.phylip_disable
                || self.fasta
                || self.write_used_sites)
        {
            return Err(ConfigError::OutputUnsupported);
        }
        if !self.inspect && self.phylip_disable && !self.fasta {
            return Err(ConfigError::NoOutputs);
        }
        if self.min_samples_locus == Some(0) {
            return Err(ConfigError::InvalidMinimum);
        }
        if self
            .min_call_rate
            .is_some_and(|rate| !rate.is_finite() || !(0.0..=1.0).contains(&rate))
        {
            return Err(ConfigError::InvalidCallRate);
        }
        if self.min_samples_locus.is_some() && self.min_call_rate.is_some() {
            return Err(ConfigError::ConflictingThresholds);
        }
        if self.resolve_iupac != self.seed.is_some() {
            return Err(ConfigError::InvalidSeed);
        }
        Ok(())
    }

    pub fn validate_samples(&self, samples: &[String]) -> Result<ValidatedConfig, ConfigError> {
        self.validate_mode()?;
        if samples.is_empty() {
            return Err(ConfigError::NoSamples);
        }
        let mut seen = HashSet::new();
        let mut outgroup_indices = Vec::new();
        for name in &self.outgroup {
            if !seen.insert(name) {
                return Err(ConfigError::DuplicateOutgroup(name.clone()));
            }
            outgroup_indices.push(
                samples
                    .iter()
                    .position(|s| s == name)
                    .ok_or_else(|| ConfigError::UnknownOutgroup(name.clone()))?,
            );
        }
        let mut warnings = Vec::new();
        let minimum_samples = match self.min_samples_locus {
            Some(requested) => {
                let minimum = requested.min(samples.len());
                if minimum != requested {
                    warnings.push(format!(
                        "minimum samples capped from {requested} to {minimum} (input sample count)"
                    ));
                }
                minimum
            }
            None => rate_minimum(
                self.min_call_rate.unwrap_or(DEFAULT_CALL_RATE),
                samples.len(),
            )
            .max(1),
        };
        Ok(ValidatedConfig {
            minimum_samples,
            outgroup_indices,
            warnings,
            processing: ProcessingConfig::new(
                samples.len(),
                Some(minimum_samples),
                self.partial_missing,
                self.seed
                    .map_or(Resolution::Iupac, |seed| Resolution::Random { seed }),
            )
            .map_err(|_| ConfigError::InvalidMinimum)?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::{Cli, rate_minimum};
    use clap::Parser;
    use std::path::PathBuf;

    #[test]
    fn output_folder_defaults_to_current_directory_and_accepts_explicit_path() {
        let default = Cli::try_parse_from(["app", "-i", "input.vcf"]).unwrap();
        assert_eq!(default.output_folder, PathBuf::from("."));
        let explicit =
            Cli::try_parse_from(["app", "-i", "input.vcf", "--output-folder", "out"]).unwrap();
        assert_eq!(explicit.output_folder, PathBuf::from("out"));
        assert!(
            Cli::try_parse_from([
                "app",
                "-i",
                "input.vcf",
                "--inspect",
                "--output-folder",
                "."
            ])
            .is_err()
        );
    }

    #[test]
    fn decimal_thresholds_match_exact_hundredths_without_float_rounding() {
        for numerator in 0..=100_usize {
            let rate = format!("{}.{:02}", numerator / 100, numerator % 100)
                .parse()
                .unwrap();
            for samples in (1..=300).chain([usize::MAX]) {
                let expected = (numerator as u128 * samples as u128).div_ceil(100);
                assert_eq!(rate_minimum(rate, samples) as u128, expected);
            }
        }
    }
}
