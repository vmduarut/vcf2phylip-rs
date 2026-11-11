//! Small CLI-only stderr reporter. Inspection's stdout is deliberately separate.
use std::{fmt::Display, time::Instant};

use vcf2phylip_rs::{cli::Cli, conversion::ConversionReport, filter::Stats, reader::Format};

/// Escape terminal controls in untrusted diagnostic text, without changing names
/// used for filesystem access or scientific output.
pub fn safe(message: impl Display) -> String {
    let mut result = String::new();
    for c in message.to_string().chars() {
        if c.is_control() {
            result.extend(c.escape_default());
        } else {
            result.push(c);
        }
    }
    result
}

/// Keep clap's own multiline layout while escaping values inserted into it.
pub fn clap_error(mut error: clap::Error) -> String {
    use clap::error::{ContextKind, ContextValue};

    let values: Vec<_> = error
        .context()
        .filter(|(kind, _)| *kind != ContextKind::Usage)
        .map(|(kind, value)| {
            let value = match value {
                ContextValue::String(s) => ContextValue::String(safe(s)),
                ContextValue::Strings(items) => {
                    ContextValue::Strings(items.iter().map(safe).collect())
                }
                ContextValue::StyledStr(s) => ContextValue::StyledStr(safe(s).into()),
                ContextValue::StyledStrs(items) => {
                    ContextValue::StyledStrs(items.iter().map(|s| safe(s).into()).collect())
                }
                _ => value.clone(),
            };
            (kind, value)
        })
        .collect();
    for (kind, value) in values {
        error.insert(kind, value);
    }
    error.to_string()
}

pub fn error_chain(error: &anyhow::Error) -> String {
    use std::fmt::Write;

    let mut chain = error.chain();
    let mut result = format!("Error: {}", safe(chain.next().unwrap()));
    let causes: Vec<_> = chain.collect();
    if !causes.is_empty() {
        result.push_str("\n\nCaused by:\n");
        for (index, cause) in causes.iter().enumerate() {
            if causes.len() == 1 {
                let _ = write!(result, "    {}", safe(cause));
            } else {
                let _ = writeln!(result, "    {index}: {}", safe(cause));
            }
        }
        if causes.len() > 1 {
            result.pop();
        }
    }
    result
}

pub struct Reporter {
    quiet: bool,
    verbose: bool,
    started: Instant,
}

impl Reporter {
    pub fn new(quiet: bool, verbose: bool) -> Self {
        Self {
            quiet,
            verbose,
            started: Instant::now(),
        }
    }

    pub fn info(&self, message: impl Display) {
        if !self.quiet {
            eprintln!("{}", safe(message));
        }
    }

    pub fn start(&self, cli: &Cli) {
        self.info(format_args!(
            "{} file '{}':",
            if cli.inspect {
                "Inspecting"
            } else {
                "Converting"
            },
            cli.input.display()
        ));
    }

    pub fn warnings(&self, warnings: &[String]) {
        for warning in warnings {
            eprintln!("warning: {}", safe(warning));
        }
    }

    pub fn header(
        &self,
        samples: usize,
        minimum_samples: usize,
        format: Format,
        compressed: bool,
        cli: &Cli,
    ) {
        self.info(format_args!("Number of samples in input: {}", samples));
        self.info(format_args!(
            "Configured absolute minimum called samples: {}",
            cli.min_samples_locus
        ));
        self.info(format_args!(
            "Minimum call rate: {} ({:.2}%)",
            cli.min_call_rate,
            cli.min_call_rate * 100.0
        ));
        self.info(format_args!(
            "Minimum called samples per locus: {} (effective)",
            minimum_samples
        ));
        self.info(format_args!("Input format: {format:?}"));
        self.info(format_args!("Input gzip/BGZF compressed: {}", compressed));
        if self.verbose {
            self.info("Configuration:");
            self.info(format_args!("Storage: {:?}", cli.storage));
            if cli.force {
                self.info("Replace existing output files: true");
            }
            self.info(format_args!(
                "Disk cell-buffer budget: {} bytes (disk only)",
                cli.memory_limit.unwrap_or(64 * 1024 * 1024)
            ));
            self.info(format_args!(
                "Partial missing policy: {:?}",
                cli.partial_missing
            ));
            self.info(format_args!("Random seed: {:?}", cli.seed));
            self.info(format_args!("Number of outgroups: {}", cli.outgroup.len()));
        }
    }

    fn stats(&self, stats: &Stats) {
        self.info(format_args!("Total records processed: {}", stats.records));
        self.info(format_args!("Records excluded as non-SNP: {} (includes absent ALT, noncanonical and duplicate alleles)", stats.rejected_non_snp));
        self.info(format_args!(
            "SNP loci excluded for insufficient called samples: {}",
            stats.nucleotide_below_minimum()
        ));
        self.info(format_args!(
            "SNP loci retained: {}",
            stats.nucleotide_accepted
        ));
        if self.verbose {
            self.info(format_args!(
                "Nucleotide missing sample calls at SNP-eligible records: {} (not rejected loci).",
                stats.nucleotide_missing_calls
            ));
        }
    }

    pub fn converted(&self, report: &ConversionReport, cli: &Cli) {
        self.warnings(&report.warnings);
        self.start(cli);
        self.header(
            report.samples,
            report.minimum_samples,
            report.input_format,
            report.input_compressed,
            cli,
        );
        self.stats(&report.stats);
        self.info(format_args!(
            "Matrix payload: {} bytes ({:?}, not peak RSS)",
            report.matrix_bytes, report.storage
        ));
        for path in &report.outputs {
            let label = match path.extension().and_then(|ext| ext.to_str()) {
                Some("phy") => "PHYLIP matrix",
                Some("fasta") => "FASTA matrix",
                Some("tsv") => "Used sites",
                _ => "Output",
            };
            self.info(format_args!("{label} saved to {}", path.display()));
        }
        self.finished("Conversion");
    }

    pub fn inspected(&self, stats: &Stats) {
        self.stats(stats);
        self.finished("Inspection");
    }

    fn finished(&self, mode: &str) {
        self.info(format_args!(
            "{mode} elapsed time: {:.3}s",
            self.started.elapsed().as_secs_f64()
        ));
        self.info("Done.");
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn error_chain_escapes_each_message_and_keeps_causes_separate() {
        let error = anyhow::anyhow!("inner\x1b\nforged")
            .context("middle\rmist")
            .context("outer\nline");
        assert_eq!(
            super::error_chain(&error),
            "Error: outer\\nline\n\nCaused by:\n    0: middle\\rmist\n    1: inner\\u{1b}\\nforged"
        );
    }
}
