use anyhow::Result;
use clap::Parser;
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
    while let Some(variant) = reader.next_variant()? {
        missing_gt += variant.genotypes.iter().filter(|gt| gt.is_none()).count() as u64;
    }
    println!("Inspection only; no alignment conversion or filtering performed.");
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
        "Validated minimum samples: {} (not applied)",
        config.minimum_samples
    );
    Ok(())
}
