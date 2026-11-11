//! Reusable conversion with streamed provenance and selectable matrix storage.
use crate::{
    cli::Cli,
    filter::{ProcessingConfig, Stats},
    matrix::{Alignment, Matrix},
    output::{self, Staged},
    reader::{Format, VariantReader},
    storage::{Storage, StorageMode},
    writers,
};
use anyhow::{Context, Result, ensure};
use std::{
    io::{BufWriter, Write},
    path::PathBuf,
};

#[derive(Debug)]
pub struct ConversionReport {
    pub stats: Stats,
    pub samples: usize,
    pub minimum_samples: usize,
    pub input_format: Format,
    pub input_compressed: bool,
    pub warnings: Vec<String>,
    pub matrix_bytes: usize,
    pub outputs: Vec<PathBuf>,
    pub storage: StorageMode,
}

/// Build in original sample order, retaining only accepted nucleotide loci.
/// Optional site reporting is streamed, including repeated coordinates in input order.
pub fn collect(
    reader: &mut VariantReader,
    config: &ProcessingConfig,
    sites: Option<&mut dyn Write>,
) -> Result<(Matrix, Stats)> {
    let matrix = Storage::Memory(Matrix::new(
        reader.header().sample_names().iter().cloned().collect(),
    )?);
    let (matrix, stats) = collect_into(reader, config, sites, matrix)?;
    let Storage::Memory(matrix) = matrix else {
        unreachable!()
    };
    Ok((matrix, stats))
}

fn collect_into(
    reader: &mut VariantReader,
    config: &ProcessingConfig,
    sites: Option<&mut dyn Write>,
    mut matrix: Storage,
) -> Result<(Storage, Stats)> {
    let mut sites = sites.map(BufWriter::new);
    if let Some(writer) = &mut sites {
        writeln!(writer, "#CHROM\tPOS\tNUM_SAMPLES")?;
    }
    let mut stats = Stats::default();
    while let Some(variant) = reader.next_variant()? {
        let locus = config.process(&variant, reader.records_read() - 1)?;
        stats.observe(&locus);
        matrix.append(&locus)?;
        if locus.nucleotide_accepted
            && let Some(writer) = &mut sites
        {
            let chrom = variant.record.reference_sequence_name();
            ensure!(
                !chrom.chars().any(char::is_control),
                "contig cannot be represented in site report: {chrom:?}"
            );
            writeln!(
                writer,
                "{}\t{}\t{}",
                chrom,
                usize::from(variant.record.variant_start().context("missing position")?),
                locus.nucleotide_called
            )?;
        }
    }
    ensure!(
        matrix.sites() > 0,
        "no nucleotide loci retained; check SNP eligibility, missing calls, --min-samples-locus and --min-call-rate"
    );
    ensure!(
        u64::try_from(matrix.sites())? == stats.nucleotide_accepted,
        "matrix/statistics length mismatch"
    );
    if let Some(writer) = &mut sites {
        writer.flush()?;
    }
    matrix.finish()?;
    Ok((matrix, stats))
}

/// Convert with opt-in output replacement. No final output is published before EOF
/// validation and all writes/flushes succeed. Finalization is per-file atomic.
/// Never prints; warnings and input metadata are returned on success without
/// reopening the input. On failure only the error is returned, not warnings.
pub fn convert(cli: &Cli) -> Result<ConversionReport> {
    cli.validate_mode()?;
    ensure!(!cli.inspect, "convert cannot be used with --inspect");
    let prefix = output::prefix(&cli.input, cli.output_prefix.as_deref())?;
    let mut reader = VariantReader::open(&cli.input)?;
    let samples: Vec<_> = reader.header().sample_names().iter().cloned().collect();
    let config = cli.validate_samples(&samples)?;
    // Validate every identifier before any output directory or temporary file.
    Matrix::new(samples.clone())?;
    let matrix = Storage::new(
        samples.clone(),
        cli.storage,
        cli.memory_limit.unwrap_or(64 * 1024 * 1024),
        cli.temp_dir.as_deref(),
    )?;
    let mut order = config.outgroup_indices.clone();
    order.extend((0..samples.len()).filter(|i| !config.outgroup_indices.contains(i)));
    let folder = cli
        .output_folder
        .clone()
        .unwrap_or_else(|| PathBuf::from("."));
    let base = format!("{prefix}.min{}", config.minimum_samples);
    let mut targets = Vec::new();
    if !cli.phylip_disable {
        targets.push(folder.join(format!("{base}.phy")));
    }
    if cli.fasta {
        targets.push(folder.join(format!("{base}.fasta")));
    }
    if cli.write_used_sites {
        targets.push(folder.join(format!("{base}.used_sites.tsv")));
    }
    output::preflight(&targets, &cli.input, cli.force)?;
    std::fs::create_dir_all(&folder)
        .with_context(|| format!("creating output directory {}", folder.display()))?;
    let mut staged = targets
        .into_iter()
        .map(Staged::new)
        .collect::<Result<Vec<_>>>()?;
    let sites = if cli.write_used_sites {
        Some(staged.last_mut().unwrap().file.as_file_mut() as &mut dyn Write)
    } else {
        None
    };
    let (matrix, stats) = collect_into(&mut reader, &config.processing, sites, matrix)?;
    let mut index = 0;
    if !cli.phylip_disable {
        writers::phylip::write(staged[index].file.as_file_mut(), &matrix, &order)?;
        index += 1;
    }
    if cli.fasta {
        writers::fasta::write(staged[index].file.as_file_mut(), &matrix, &order)?;
    }
    // Finish every staged write before publishing any output. Never roll back by
    // deleting paths: in force mode they may have held preexisting user files.
    for file in &mut staged {
        file.file.flush()?;
    }
    let mut outputs = Vec::new();
    for file in staged {
        outputs.push(file.persist(&cli.input, cli.force)?);
    }
    Ok(ConversionReport {
        stats,
        samples: samples.len(),
        minimum_samples: config.minimum_samples,
        input_format: reader.format(),
        input_compressed: reader.is_compressed(),
        warnings: config.warnings,
        matrix_bytes: matrix
            .sites()
            .checked_mul(samples.len())
            .context("matrix size overflow")?,
        storage: cli.storage,
        outputs,
    })
}
