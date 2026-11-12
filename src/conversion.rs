//! Reusable conversion with streamed provenance and selectable matrix storage.
use crate::{
    cli::Cli,
    matrix::{Alignment, Matrix},
    output::{self, StagedOutput},
    processing::{Stats, VariantProcessor},
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
    pub sample_count: usize,
    pub minimum_samples: usize,
    pub input_format: Format,
    pub input_compressed: bool,
    pub warnings: Vec<String>,
    pub matrix_payload_bytes: usize,
    pub outputs: Vec<PathBuf>,
    pub storage: StorageMode,
}

/// Process and filter variants, accumulating accepted sites and statistics.
fn process_variants(
    reader: &mut VariantReader,
    processor: &VariantProcessor,
    site_report_writer: Option<&mut dyn Write>,
    mut matrix: Storage,
) -> Result<(Storage, Stats)> {
    let mut site_report_writer = site_report_writer.map(BufWriter::new);
    if let Some(writer) = &mut site_report_writer {
        writeln!(writer, "#CHROM\tPOS\tNUM_SAMPLES")?;
    }
    let mut stats = Stats::default();
    while let Some(variant) = reader.next_variant()? {
        let locus = processor.process(&variant, reader.records_read() - 1)?;
        stats.observe(&locus);
        matrix.append(&locus)?;
        if locus.nucleotide_accepted
            && let Some(writer) = &mut site_report_writer
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
                locus.nucleotide_called_samples
            )?;
        }
    }
    ensure!(
        matrix.sites() > 0,
        "no nucleotide loci retained; check SNP eligibility, missing calls, and the selected --min-samples-locus or --min-call-rate threshold"
    );
    ensure!(
        u64::try_from(matrix.sites())? == stats.nucleotide_sites_accepted,
        "matrix/statistics length mismatch"
    );
    if let Some(writer) = &mut site_report_writer {
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
    let mut sample_order = config.outgroup_indices.clone();
    sample_order.extend((0..samples.len()).filter(|i| !config.outgroup_indices.contains(i)));
    let folder = &cli.output_folder;
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
    std::fs::create_dir_all(folder)
        .with_context(|| format!("creating output directory {}", folder.display()))?;
    let mut staged = targets
        .into_iter()
        .map(StagedOutput::new)
        .collect::<Result<Vec<_>>>()?;
    let site_report_writer = if cli.write_used_sites {
        Some(staged.last_mut().unwrap().file.as_file_mut() as &mut dyn Write)
    } else {
        None
    };
    let (matrix, stats) =
        process_variants(&mut reader, &config.processor, site_report_writer, matrix)?;
    let mut staged_index = 0;
    if !cli.phylip_disable {
        writers::phylip::write(
            staged[staged_index].file.as_file_mut(),
            &matrix,
            &sample_order,
        )?;
        staged_index += 1;
    }
    if cli.fasta {
        writers::fasta::write(
            staged[staged_index].file.as_file_mut(),
            &matrix,
            &sample_order,
        )?;
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
        sample_count: samples.len(),
        minimum_samples: config.minimum_samples,
        input_format: reader.format(),
        input_compressed: reader.is_compressed(),
        warnings: config.warnings,
        matrix_payload_bytes: matrix
            .sites()
            .checked_mul(samples.len())
            .context("matrix size overflow")?,
        storage: cli.storage,
        outputs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genotype::{NucleotideResolution, PartialMissingPolicy};
    use std::io::Cursor;

    fn reader(records: &str) -> VariantReader {
        let text = format!(
            "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tlong:sample,Ω\n{records}"
        );
        VariantReader::from_reader(Cursor::new(text.into_bytes())).unwrap()
    }

    fn processor() -> VariantProcessor {
        VariantProcessor::new(
            1,
            None,
            PartialMissingPolicy::Missing,
            NucleotideResolution::Iupac,
        )
        .unwrap()
    }

    fn memory() -> Storage {
        Storage::Memory(Matrix::new(vec!["long:sample,Ω".into()]).unwrap())
    }

    struct Failing;
    impl Write for Failing {
        fn write(&mut self, _: &[u8]) -> std::io::Result<usize> {
            Err(std::io::Error::other("injected failure"))
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Err(std::io::Error::other("injected flush failure"))
        }
    }

    #[test]
    fn library_writers_wrap_validate_orders_and_propagate_io_errors() {
        let records = "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n".repeat(81);
        let (matrix, stats) =
            process_variants(&mut reader(&records), &processor(), None, memory()).unwrap();
        assert_eq!(stats.nucleotide_sites_accepted, 81);
        let mut bytes = Vec::new();
        writers::fasta::write(&mut bytes, &matrix, &[0]).unwrap();
        assert_eq!(
            String::from_utf8(bytes).unwrap(),
            format!(">long:sample,Ω\n{}\nM\n", "M".repeat(80))
        );
        assert!(writers::phylip::write(Vec::new(), &matrix, &[]).is_err());
        assert!(writers::fasta::write(Vec::new(), &matrix, &[1]).is_err());
        assert!(writers::fasta::write(Failing, &matrix, &[0]).is_err());
        assert!(writers::phylip::write(Failing, &matrix, &[0]).is_err());
    }

    #[test]
    fn memory_processing_streams_only_retained_sites_in_input_order() {
        let records = "chr1\t9\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\nchr1\t2\t.\tA\tAT\t.\tPASS\t.\tGT\t0/0\nchr1\t1\t.\tA\tC,G,T\t.\tPASS\t.\tGT\t0/1/2/3\nchr1\t1\t.\tG\tT\t.\tPASS\t.\tGT\t./.\nchr1\t1\t.\tG\tT\t.\tPASS\t.\tGT\t0/1\n";
        let mut sites = Vec::new();
        let (matrix, stats) = process_variants(
            &mut reader(records),
            &processor(),
            Some(&mut sites),
            memory(),
        )
        .unwrap();
        let Storage::Memory(matrix) = matrix else {
            panic!("expected memory storage")
        };
        assert_eq!(matrix.samples(), &["long:sample,Ω"]);
        assert_eq!(matrix.sequences(), &[b"MNK".to_vec()]);
        assert_eq!(matrix.sites() as u64, stats.nucleotide_sites_accepted);
        assert_eq!(stats.records_processed, 5);
        assert_eq!(stats.rejected_non_snp, 1);
        assert_eq!(stats.nucleotide_sites_accepted, 3);
        assert_eq!(stats.nucleotide_missing_calls, 1);
        assert_eq!(
            sites,
            b"#CHROM\tPOS\tNUM_SAMPLES\nchr1\t9\t1\nchr1\t1\t1\nchr1\t1\t1\n"
        );
        assert!(
            process_variants(
                &mut reader(records),
                &processor(),
                Some(&mut Failing),
                memory()
            )
            .is_err()
        );
    }

    #[test]
    fn memory_processing_rejects_invalid_and_zero_retained_input() {
        for records in [
            "",
            "chr1\t1\t.\tA\tAT\t.\tPASS\t.\tGT\t0/0\n",
            "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t./.\n",
        ] {
            let error = process_variants(&mut reader(records), &processor(), None, memory())
                .err()
                .unwrap();
            assert!(error.to_string().contains("no nucleotide loci retained"));
        }
        let records = "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\nbad late record\n";
        assert!(process_variants(&mut reader(records), &processor(), None, memory()).is_err());
    }

    #[test]
    fn memory_processing_rejects_matrix_statistics_mismatch() {
        let processor = processor();
        let mut input = reader("chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n");
        let locus = processor
            .process(&input.next_variant().unwrap().unwrap(), 0)
            .unwrap();
        let mut matrix = memory();
        matrix.append(&locus).unwrap();
        let error = process_variants(&mut input, &processor, None, matrix)
            .err()
            .unwrap();
        assert_eq!(error.to_string(), "matrix/statistics length mismatch");
    }
}
