use anyhow::Result;
use clap::Parser;
use vcf2phylip_rs::processing::Stats;
use vcf2phylip_rs::{
    cli::{Cli, DEFAULT_CALL_RATE},
    reader::VariantReader,
};

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

fn main() {
    let args = match Args::try_parse() {
        Ok(args) => args,
        Err(error) if error.use_stderr() => {
            eprint!("{}", logging::clap_error(error));
            std::process::exit(2);
        }
        Err(error) => error.exit(),
    };
    if let Err(error) = run(args) {
        eprintln!("{}", logging::error_chain(&error));
        std::process::exit(1);
    }
}

fn run(args: Args) -> Result<()> {
    let cli = args.cli;
    cli.validate_mode()?;
    let log = logging::Reporter::new(args.quiet, args.verbose);
    if !cli.inspect {
        let report = vcf2phylip_rs::conversion::convert(&cli)?;
        log.converted(&report, &cli);
        return Ok(());
    }
    let mut reader = VariantReader::open(&cli.input)?;
    let samples: Vec<_> = reader.header().sample_names().iter().cloned().collect();
    let config = cli.validate_samples(&samples)?;
    log.warnings(&config.warnings);
    let mut absent_gt_calls = 0_u64;
    let mut stats = Stats::default();
    while let Some(variant) = reader.next_variant()? {
        absent_gt_calls += variant.genotypes.iter().filter(|gt| gt.is_none()).count() as u64;
        stats.observe(
            &config
                .processor
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
        "Samples: {}; records: {}; absent GT calls (unknown ploidy): {absent_gt_calls}",
        samples.len(),
        reader.records_read()
    );
    if let Some(minimum) = cli.min_samples_locus {
        println!("Configured absolute minimum called samples: {minimum}");
    } else {
        let rate = cli.min_call_rate.unwrap_or(DEFAULT_CALL_RATE);
        println!("Minimum call rate: {rate} ({:.2}%)", rate * 100.0);
    }
    println!(
        "Minimum called samples: {} (effective; applied independently to nucleotide and binary loci)",
        config.minimum_samples
    );
    println!(
        "Accepted nucleotide loci: {}; accepted binary loci: {}; binary eligible: {}; rejected non-SNP: {}",
        stats.nucleotide_sites_accepted,
        stats.binary_sites_accepted,
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
    log.start(&cli);
    log.header(
        samples.len(),
        config.minimum_samples,
        reader.format(),
        reader.is_compressed(),
        &cli,
    );
    log.inspected(&stats);
    Ok(())
}
