use anyhow::Result;
use clap::Parser;
use vcf2phylip_rs::filter::Stats;
use vcf2phylip_rs::{cli::Cli, reader::VariantReader};

fn main() -> Result<()> {
    let cli = Cli::parse();
    cli.validate_mode()?;
    let mut reader = VariantReader::open(&cli.input)?;
    let samples: Vec<_> = reader.header().sample_names().iter().cloned().collect();
    let config = cli.validate_samples(&samples)?;
    for warning in config.warnings {
        eprintln!("warning: {warning}");
    }
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
    Ok(())
}
