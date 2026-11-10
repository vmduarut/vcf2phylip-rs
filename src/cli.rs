//! Compatibility arguments and validation for inspection and conversion.
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
    about = "Convert VCF/BCF to in-memory FASTA and relaxed PHYLIP alignments"
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
    /// Output directory (defaults to the current directory).
    #[arg(long)]
    pub output_folder: Option<PathBuf>,
    /// Output basename (defaults to input basename without known format suffixes).
    #[arg(long)]
    pub output_prefix: Option<String>,
    /// Disable the default relaxed PHYLIP output.
    #[arg(short = 'p', long)]
    pub phylip_disable: bool,
    /// Also write FASTA (80 columns per sequence line).
    #[arg(short = 'f', long)]
    pub fasta: bool,
    /// Reserved for Milestone 4 (not yet supported).
    #[arg(short = 'n', long)]
    pub nexus: bool,
    /// Reserved for Milestone 4 (not yet supported).
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
    /// Write retained nucleotide site coordinates and called sample counts.
    #[arg(short = 'w', long)]
    pub write_used_sites: bool,
    #[arg(short = 'v', long, action = ArgAction::Version)]
    pub version: (),
}

#[derive(Debug, Error)]
pub enum ConfigError {
    #[error("NEXUS and binary NEXUS output are not implemented (Milestone 4)")]
    NexusUnsupported,
    #[error("no alignment output selected; use --fasta with --phylip-disable")]
    NoOutputs,
    #[error("conversion/output flags cannot be used with --inspect")]
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
        if self.nexus || self.nexus_binary {
            return Err(ConfigError::NexusUnsupported);
        }
        if self.inspect
            && (self.output_folder.is_some()
                || self.output_prefix.is_some()
                || self.phylip_disable
                || self.fasta
                || self.nexus
                || self.nexus_binary
                || self.write_used_sites)
        {
            return Err(ConfigError::OutputUnsupported);
        }
        if !self.inspect && self.phylip_disable && !self.fasta {
            return Err(ConfigError::NoOutputs);
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
