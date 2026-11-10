//! Compatibility argument names; output-producing options are reserved, not implemented.
use std::{collections::HashSet, path::PathBuf};

use crate::{
    filter::ProcessingConfig,
    genotype::{PartialMissingPolicy, Resolution},
};
use clap::{ArgAction, Parser};
use thiserror::Error;

#[derive(Debug, Parser)]
#[command(
    version,
    disable_version_flag = true,
    about = "VCF/BCF validation and scientific inspection (Milestone 2; no alignment conversion)"
)]
pub struct Cli {
    #[arg(short, long)]
    pub input: PathBuf,
    /// Stream and validate input; print diagnostic counts only, never an alignment.
    #[arg(long)]
    pub inspect: bool,
    #[arg(short = 'm', long, default_value = "4", value_parser = positive_usize)]
    pub min_samples_locus: usize,
    /// Exact sample name; repeat for multiple outgroups (no delimiter splitting).
    #[arg(short, long)]
    pub outgroup: Vec<String>,
    /// Reserved conversion option (not supported in inspection mode).
    #[arg(long)]
    pub output_folder: Option<PathBuf>,
    /// Reserved conversion option.
    #[arg(long)]
    pub output_prefix: Option<String>,
    /// Reserved conversion option.
    #[arg(short = 'p', long)]
    pub phylip_disable: bool,
    /// Reserved conversion option.
    #[arg(short = 'f', long)]
    pub fasta: bool,
    /// Reserved conversion option.
    #[arg(short = 'n', long)]
    pub nexus: bool,
    /// Reserved conversion option.
    #[arg(short = 'b', long)]
    pub nexus_binary: bool,
    /// Resolve nucleotide calls by sampling allele copies; requires --seed.
    #[arg(short = 'r', long = "resolve-IUPAC")]
    pub resolve_iupac: bool,
    /// Explicit reproducible seed for --resolve-IUPAC.
    #[arg(long)]
    pub seed: Option<u64>,
    /// Treatment of partially missing nucleotide genotypes (binary is always strict).
    #[arg(long, value_enum, default_value = "missing")]
    pub partial_missing: PartialMissingPolicy,
    /// Reserved conversion option.
    #[arg(short = 'w', long)]
    pub write_used_sites: bool,
    #[arg(short = 'v', long, action = ArgAction::Version)]
    pub version: (),
}

#[derive(Debug, Error)]
pub enum ConfigError {
    #[error("alignment conversion is not implemented; use --inspect for validation only")]
    ConversionUnsupported,
    #[error("conversion/output flags cannot be used with --inspect; no outputs are implemented")]
    OutputUnsupported,
    #[error("minimum samples must be positive")]
    InvalidMinimum,
    #[error("outgroup is not an exact input sample name: {0:?}")]
    UnknownOutgroup(String),
    #[error("outgroup specified more than once: {0:?}")]
    DuplicateOutgroup(String),
    #[error("input has no samples")]
    NoSamples,
    #[error("--resolve-IUPAC requires --seed; --seed requires --resolve-IUPAC")]
    InvalidSeed,
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

impl Cli {
    pub fn validate_mode(&self) -> Result<(), ConfigError> {
        if !self.inspect {
            return Err(ConfigError::ConversionUnsupported);
        }
        if self.output_folder.is_some()
            || self.output_prefix.is_some()
            || self.phylip_disable
            || self.fasta
            || self.nexus
            || self.nexus_binary
            || self.write_used_sites
        {
            return Err(ConfigError::OutputUnsupported);
        }
        if self.min_samples_locus == 0 {
            return Err(ConfigError::InvalidMinimum);
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
        let minimum_samples = self.min_samples_locus.min(samples.len());
        let warnings = if minimum_samples != self.min_samples_locus {
            vec![format!(
                "minimum samples capped from {} to {minimum_samples} (input sample count)",
                self.min_samples_locus
            )]
        } else {
            Vec::new()
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
