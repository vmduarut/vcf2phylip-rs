use anyhow::Result;
use clap::Parser;
use vcf2phylip_rs::filter::Stats;
use vcf2phylip_rs::{cli::Cli, reader::VariantReader};

mod logging;

/// Binary-only options leave the library's public configuration shape unchanged.
#[derive(Debug, Parser)]
#[command(
    version,
    disable_version_flag = true,
    about = "Convert VCF/BCF to FASTA and relaxed PHYLIP alignments"
)]
struct Args {
    #[command(flatten)]
    cli: Cli,
    /// Suppress informational stderr messages (not warnings or errors).
    #[arg(long, conflicts_with = "verbose")]
    quiet: bool,
    /// Include extra configuration and call diagnostics on stderr.
    #[arg(long, conflicts_with = "quiet")]
    verbose: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let cli = args.cli;
    cli.validate_mode()?;
    let log = logging::Reporter::new(args.quiet, args.verbose);
    log.start(&cli);
    if !cli.inspect {
        let report = vcf2phylip_rs::conversion::convert_with_observer(&cli, |reader, config| {
            log.header(reader, config, &cli);
        })?;
        log.converted(&report);
        return Ok(());
    }
    let mut reader = VariantReader::open(&cli.input)?;
    let samples: Vec<_> = reader.header().sample_names().iter().cloned().collect();
    let config = cli.validate_samples(&samples)?;
    log.header(&reader, &config, &cli);
    let mut missing_gt = 0_u64;
    let mut stats = Stats::default();
    while let Some(variant) = reader.next_variant()? {
        missing_gt += variant.genotypes.iter().filter(|gt| gt.is_none()).count() as u64;
        stats.observe(
            &config
                .processing
                .process(&variant, reader.records_read() - 1)?,
        );
    }
    println!("Inspection only; scientific eligibility assessed, no alignment outputs created.");
    println!(
        "Format: {:?}; gzip/BGZF: {}",
        reader.format(),
        reader.is_compressed()
    );
    println!(
        "Samples: {}; records: {}; absent GT calls (unknown ploidy): {missing_gt}",
        samples.len(),
        reader.records_read()
    );
    println!(
        "Minimum called samples: {} (applied independently to nucleotide and binary loci)",
        config.minimum_samples
    );
    println!(
        "Accepted nucleotide loci: {}; accepted binary loci: {}; binary eligible: {}; rejected non-SNP: {}",
        stats.nucleotide_accepted,
        stats.binary_accepted,
        stats.binary_eligible,
        stats.rejected_non_snp
    );
    println!(
        "Nucleotide missing calls: {}; binary missing calls: {}; binary unsupported-ploidy calls: {}",
        stats.nucleotide_missing_calls,
        stats.binary_missing_calls,
        stats.binary_unsupported_ploidy_calls
    );
    println!(
        "Nucleotide partial policy: {:?}; random seed: {:?}",
        cli.partial_missing, cli.seed
    );
    log.inspected(&stats);
    Ok(())
}
