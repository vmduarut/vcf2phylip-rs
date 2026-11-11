//! Small CLI-only stderr reporter. Inspection's stdout is deliberately separate.
use std::{fmt::Display, time::Instant};

use vcf2phylip_rs::{
    cli::{Cli, ValidatedConfig},
    conversion::ConversionReport,
    filter::Stats,
    reader::VariantReader,
};

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
            eprintln!("{message}");
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

    pub fn header(&self, reader: &VariantReader, config: &ValidatedConfig, cli: &Cli) {
        for warning in &config.warnings {
            eprintln!("warning: {warning}");
        }
        self.info(format_args!(
            "Number of samples in input: {}",
            reader.header().sample_names().len()
        ));
        self.info(format_args!(
            "Minimum called samples per locus: {}",
            config.minimum_samples
        ));
        self.info(format_args!("Input format: {:?}", reader.format()));
        self.info(format_args!(
            "Input gzip/BGZF compressed: {}",
            reader.is_compressed()
        ));
        if self.verbose {
            self.info("Configuration:");
            self.info(format_args!("Storage: {:?}", cli.storage));
            self.info(format_args!(
                "Disk cell-buffer budget: {} bytes (disk only)",
                cli.memory_limit.unwrap_or(64 * 1024 * 1024)
            ));
            self.info(format_args!(
                "Partial missing policy: {:?}",
                cli.partial_missing
            ));
            self.info(format_args!("Random seed: {:?}", cli.seed));
            self.info(format_args!(
                "Number of outgroups: {}",
                config.outgroup_indices.len()
            ));
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

    pub fn converted(&self, report: &ConversionReport) {
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
