use std::{
    io::{Cursor, Write},
    process::Command,
};

use clap::Parser;
use flate2::{Compression, write::GzEncoder};
use noodles::{bcf, bgzf, vcf};
use vcf2phylip_rs::{
    cli::Cli,
    reader::{Format, VariantReader},
};

const META: &str = "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n";

fn data(samples: &str, records: &str) -> String {
    format!("{META}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\n{records}")
}

fn reader(text: &str) -> anyhow::Result<VariantReader> {
    VariantReader::from_reader(Cursor::new(text.as_bytes().to_vec()))
}

fn one(gt: &str) -> String {
    data(
        "./sample A",
        &format!("chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t{gt}\n"),
    )
}

fn gzip(bytes: &[u8]) -> Vec<u8> {
    let mut writer = GzEncoder::new(Vec::new(), Compression::default());
    writer.write_all(bytes).unwrap();
    writer.finish().unwrap()
}

fn bgzf(bytes: &[u8]) -> Vec<u8> {
    let mut writer = bgzf::io::Writer::new(Vec::new());
    writer.write_all(bytes).unwrap();
    writer.finish().unwrap()
}

fn bcf_bytes(text: &str) -> Vec<u8> {
    let mut input = vcf::io::Reader::new(Cursor::new(text));
    let header = input.read_header().unwrap();
    let mut writer = bcf::io::Writer::new(Vec::new());
    writer.write_header(&header).unwrap();
    for record in input.record_bufs(&header) {
        use vcf::variant::io::Write;
        writer
            .write_variant_record(&header, &record.unwrap())
            .unwrap();
    }
    writer.try_finish().unwrap();
    writer.into_inner().finish().unwrap()
}

/// Independently encoded BCF 2.2, not produced by noodles' BCF writer.
/// Layout follows VCFv4.3.pdf, section "BCF version 2.2": little-endian
/// shared/individual lengths, fixed shared fields, typed strings, then FORMAT
/// key/type and sample-major GT vectors. Int8 vector-end is 0x81 (not missing
/// allele 0); GT alleles are (VCF_index + 1) << 1 with phase bit zero here.
/// Each sample occupies exactly three bytes; padding follows the ENTIRE GT,
/// never individual alleles. Explicit IDX entries fix PASS=0, GT=1, chr1=0.
fn independent_mixed_ploidy_bcf() -> Vec<u8> {
    let header = "##fileformat=VCFv4.3\n##FILTER=<ID=PASS,Description=\"All filters passed\",IDX=0>\n##contig=<ID=chr1,length=1000,IDX=0>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\",IDX=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tH\tD\tP\tM\n\0";
    let mut bcf = b"BCF\x02\x02".to_vec();
    bcf.extend(u32::try_from(header.len()).unwrap().to_le_bytes());
    bcf.extend(header.as_bytes());
    for (pos, genotypes) in [
        // H=1, D=0/1, P=0/1/1, M=./0.
        (
            8_i32,
            [[4, 0x81, 0x81], [2, 4, 0x81], [2, 4, 4], [0, 2, 0x81]],
        ),
        // H=0, D=./., P=0/./1, M=1/1.
        (
            1_i32,
            [[2, 0x81, 0x81], [0, 0, 0x81], [2, 0, 4], [4, 4, 0x81]],
        ),
    ] {
        let mut shared = Vec::new();
        shared.extend(0_i32.to_le_bytes()); // CHROM dictionary index
        shared.extend(pos.to_le_bytes()); // zero-based POS
        shared.extend(1_i32.to_le_bytes()); // REF length
        shared.extend(0x7f80_0001_u32.to_le_bytes()); // missing QUAL
        shared.extend(0_u16.to_le_bytes()); // n_info
        shared.extend(2_u16.to_le_bytes()); // n_allele
        shared.extend(((1_u32 << 24) | 4).to_le_bytes()); // n_fmt=1, n_sample=4
        shared.extend([0x07, 0x17, b'A', 0x17, b'C', 0x11, 0]); // empty ID, REF, ALT, PASS
        let mut individual = vec![0x11, 1, 0x31]; // typed GT key=1; three int8 values/sample
        individual.extend(genotypes.into_iter().flatten());
        bcf.extend(u32::try_from(shared.len()).unwrap().to_le_bytes());
        bcf.extend(u32::try_from(individual.len()).unwrap().to_le_bytes());
        bcf.extend(shared);
        bcf.extend(individual);
    }
    bcf
}

#[test]
fn independent_bcf_mixed_ploidy_parsing_counts_and_conversion() {
    use vcf2phylip_rs::{conversion::convert, genotype::PartialMissingPolicy};
    let vcf = data(
        "H\tD\tP\tM",
        "chr1\t9\t.\tA\tC\t.\tPASS\t.\tGT\t1\t0/1\t0/1/1\t./0\nchr1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t0\t./.\t0/./1\t1/1\n",
    );
    let raw = independent_mixed_ploidy_bcf();
    let expected = [
        vec![
            vec![Some(1)],
            vec![Some(0), Some(1)],
            vec![Some(0), Some(1), Some(1)],
            vec![None, Some(0)],
        ],
        vec![
            vec![Some(0)],
            vec![None, None],
            vec![Some(0), None, Some(1)],
            vec![Some(1), Some(1)],
        ],
    ];
    // Assert actual ploidy and missing-allele positions, not merely symbols:
    // vector-end padding must not become missing calls or truncate other samples.
    for bytes in [vcf.as_bytes().to_vec(), raw.clone(), bgzf(&raw), gzip(&raw)] {
        let mut input = VariantReader::from_reader(Cursor::new(bytes)).unwrap();
        for (ordinal, expected) in expected.iter().enumerate() {
            let variant = input.next_variant().unwrap().unwrap();
            assert_eq!(
                usize::from(variant.record.variant_start().unwrap()),
                [9, 2][ordinal]
            );
            let actual: Vec<Vec<_>> = variant
                .genotypes
                .iter()
                .map(|gt| {
                    gt.as_ref()
                        .unwrap()
                        .as_ref()
                        .iter()
                        .map(|a| a.position())
                        .collect()
                })
                .collect();
            assert_eq!(&actual, expected);
        }
        assert!(input.next_variant().unwrap().is_none());
    }
    for (policy, called, missing, sequences) in [
        (
            PartialMissingPolicy::Missing,
            [3, 2],
            3,
            ["CA", "MN", "MN", "NC"],
        ),
        (
            PartialMissingPolicy::CalledAlleles,
            [4, 3],
            1,
            ["CA", "MN", "MM", "AC"],
        ),
    ] {
        let expected_phy = format!(
            "4 2\nH  {}\nD  {}\nP  {}\nM  {}\n",
            sequences[0], sequences[1], sequences[2], sequences[3]
        );
        let expected_fasta = format!(
            ">H\n{}\n>D\n{}\n>P\n{}\n>M\n{}\n",
            sequences[0], sequences[1], sequences[2], sequences[3]
        );
        for bytes in [vcf.as_bytes().to_vec(), raw.clone(), bgzf(&raw), gzip(&raw)] {
            let dir = tempfile::tempdir().unwrap();
            let input = dir.path().join("input");
            std::fs::write(&input, bytes).unwrap();
            let mut cli = Cli::try_parse_from([
                "app",
                "-i",
                input.to_str().unwrap(),
                "--output-folder",
                dir.path().to_str().unwrap(),
                "-m",
                "2",
                "-f",
                "-w",
            ])
            .unwrap();
            cli.partial_missing = policy;
            let report = convert(&cli).unwrap();
            let stats = report.stats;
            assert_eq!(stats.records_processed, 2);
            assert_eq!(stats.nucleotide_sites_accepted, 2);
            assert_eq!(stats.nucleotide_missing_calls, missing);
            assert_eq!(stats.binary_eligible, 2);
            assert_eq!(stats.binary_sites_accepted, 0);
            assert_eq!(stats.binary_missing_calls, 2);
            assert_eq!(stats.binary_unsupported_ploidy_calls, 4);
            assert_eq!(
                std::fs::read_to_string(dir.path().join("input.min2.used_sites.tsv")).unwrap(),
                format!(
                    "#CHROM\tPOS\tNUM_SAMPLES\nchr1\t9\t{}\nchr1\t2\t{}\n",
                    called[0], called[1]
                )
            );
            let phy = std::fs::read(dir.path().join("input.min2.phy")).unwrap();
            let fasta = std::fs::read(dir.path().join("input.min2.fasta")).unwrap();
            assert_eq!(phy, expected_phy.as_bytes());
            assert_eq!(fasta, expected_fasta.as_bytes());
        }
    }
}

#[test]
fn public_writers_reject_zero_sites_before_writing() {
    use vcf2phylip_rs::{matrix::Matrix, writers};
    let matrix = Matrix::new(vec!["S".into()]).unwrap();
    let mut output = Vec::new();
    assert!(
        writers::fasta::write(&mut output, &matrix, &[0])
            .unwrap_err()
            .to_string()
            .contains("zero-site")
    );
    assert!(output.is_empty());
    assert!(
        writers::phylip::write(&mut output, &matrix, &[0])
            .unwrap_err()
            .to_string()
            .contains("zero-site")
    );
    assert!(output.is_empty());
}

#[test]
fn genotype_indices_missingness_ploidy_and_phase() {
    use vcf::variant::record::samples::series::value::genotype::Phasing;
    for (gt, expected) in [
        ("0", vec![Some(0)]),
        ("0/1", vec![Some(0), Some(1)]),
        ("0|1|1", vec![Some(0), Some(1), Some(1)]),
        ("./.", vec![None, None]),
        ("0/.", vec![Some(0), None]),
        ("./0", vec![None, Some(0)]),
    ] {
        for bytes in [one(gt).into_bytes(), bcf_bytes(&one(gt))] {
            let mut input = VariantReader::from_reader(Cursor::new(bytes)).unwrap();
            let variant = input.next_variant().unwrap().unwrap();
            let alleles = variant.genotypes[0].as_ref().unwrap().as_ref();
            assert_eq!(
                alleles.iter().map(|a| a.position()).collect::<Vec<_>>(),
                expected
            );
            if gt.contains('|') {
                assert!(alleles.iter().all(|a| a.phasing() == Phasing::Phased));
            }
            assert!(input.next_variant().unwrap().is_none());
            assert_eq!(input.records_read(), 1);
        }
    }
}

#[test]
fn multidigit_indices_and_symbolic_alleles_are_not_rewritten() {
    let text = data(
        "S",
        "chr1\t1\t.\tA\tC,G,T,C,G,T,C,G,T,<NON_REF>\t.\tPASS\t.\tGT\t10/0\n",
    );
    let variant = reader(&text).unwrap().next_variant().unwrap().unwrap();
    assert_eq!(
        variant.genotypes[0].as_ref().unwrap().as_ref()[0].position(),
        Some(10)
    );
    assert_eq!(variant.record.alternate_bases().as_ref()[9], "<NON_REF>");
}

#[test]
fn absent_gt_unknown_ploidy_and_missing_allele_are_distinct() {
    let text = data(
        "A\tB\tC",
        "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t.\t./.\t0/.\nchr1\t2\t.\tA\tC\t.\tPASS\t.\tDP\t2\t.\t4\n",
    );
    let mut input = reader(&text).unwrap();
    let first = input.next_variant().unwrap().unwrap();
    assert!(first.genotypes[0].is_none());
    assert!(first.genotypes[1].is_some());
    let second = input.next_variant().unwrap().unwrap();
    assert!(second.genotypes.iter().all(Option::is_none));
}

#[test]
fn formats_detected_by_content_and_fully_streamed() {
    let text = one("0|1");
    let split = text.len() / 2;
    let mut multi = gzip(&text.as_bytes()[..split]);
    multi.extend(gzip(&text.as_bytes()[split..]));
    for (bytes, format, compressed) in [
        (text.as_bytes().to_vec(), Format::Vcf, false),
        (gzip(text.as_bytes()), Format::Vcf, true),
        (multi, Format::Vcf, true),
        (bgzf(text.as_bytes()), Format::Vcf, true),
        (bcf_bytes(&text), Format::Bcf, true),
    ] {
        let mut input = VariantReader::from_reader(Cursor::new(bytes)).unwrap();
        assert_eq!(input.format(), format);
        assert_eq!(input.is_compressed(), compressed);
        assert_eq!(&input.header().sample_names()[0], "./sample A");
        assert!(input.header().contigs().contains_key("chr1"));
        assert!(input.next_variant().unwrap().is_some());
        assert!(input.next_variant().unwrap().is_none());
        assert_eq!(input.records_read(), 1);
    }
}

#[test]
fn header_only_is_valid_including_bcf() {
    let text = data("S", "");
    for bytes in [
        text.as_bytes().to_vec(),
        gzip(text.as_bytes()),
        bgzf(text.as_bytes()),
        bcf_bytes(&text),
    ] {
        let mut input = VariantReader::from_reader(Cursor::new(bytes)).unwrap();
        assert!(input.next_variant().unwrap().is_none());
        assert_eq!(input.records_read(), 0);
    }
}

#[test]
fn malformed_headers_rejected() {
    for text in [
        String::new(),
        "not VCF".into(),
        format!("{META}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"),
        data("S\tS", ""),
        data("", ""),
        data("S", "").replace(
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
            "",
        ),
        data("S", "").replace("ID=GT,Number=1", "ID=GT,Number=2"),
        "##fileformat=VCFv4.3\n".into(),
    ] {
        assert!(reader(&text).is_err(), "accepted {text:?}");
    }
}

#[test]
fn malformed_records_fail_strictly_with_context() {
    for gt in ["2/0", "-1/0", "x", "0//1", "0/", "0:", ""] {
        let mut input = reader(&one(gt)).unwrap();
        let error = format!("{:#}", input.next_variant().unwrap_err());
        assert!(error.contains("record 1"), "{error}");
        assert_eq!(input.records_read(), 0);
        assert!(input.next_variant().unwrap().is_none());
    }
    let good = one("0/1");
    for bad in [
        good.replace("GT\t0/1", "GT\t0/1\t0/0"),
        good.replace("GT\t0/1", "GT"),
        good.replace("GT\t0/1", "GT:GT\t0/1:0/0"),
        good.replace("GT\t0/1", "DP:GT\t4:0/1"),
        good.replace("chr1\t1\t", "chr1\t0\t"),
        good.replace("chr1\t1\t", "chr1\tx\t"),
    ] {
        assert!(
            reader(&bad).unwrap().next_variant().is_err(),
            "accepted {bad:?}"
        );
    }
}

#[test]
fn empty_format_keys_are_rejected() {
    for format in [":GT", "GT:", "GT::DP"] {
        let text = data(
            "S",
            &format!("chr1\t1\t.\tA\tC\t.\tPASS\t.\t{format}\t0/1\n"),
        );
        let error = format!("{:#}", reader(&text).unwrap().next_variant().unwrap_err());
        assert!(error.contains("empty FORMAT key"), "{format}: {error}");
    }
}

#[test]
fn compressed_truncation_and_checksum_errors_rejected() {
    let text = one("0/1");
    for mut bytes in [
        gzip(text.as_bytes()),
        bgzf(text.as_bytes()),
        bcf_bytes(&text),
    ] {
        // Removing beyond the BGZF EOF block ensures a genuinely incomplete data member.
        let remove = if bytes.ends_with(&[0; 8]) { 35 } else { 5 };
        bytes.truncate(bytes.len() - remove);
        assert!(drain(bytes).is_err());
    }
    let mut corrupt = gzip(text.as_bytes());
    let i = corrupt.len() - 8;
    corrupt[i] ^= 1;
    assert!(drain(corrupt).is_err());
}

fn drain(bytes: Vec<u8>) -> anyhow::Result<()> {
    let mut input = VariantReader::from_reader(Cursor::new(bytes))?;
    while input.next_variant()?.is_some() {}
    Ok(())
}

#[test]
fn truncated_bcf_record_prefix_is_not_eof() {
    use std::io::Read;
    let encoded = bcf_bytes(&data("S", ""));
    let mut raw = Vec::new();
    flate2::read::MultiGzDecoder::new(&encoded[..])
        .read_to_end(&mut raw)
        .unwrap();
    for prefix in [vec![1], vec![1, 0], vec![1, 0, 0]] {
        let mut bytes = raw.clone();
        bytes.extend(prefix);
        assert!(drain(bytes).is_err());
    }
}

#[test]
fn bcf_rejects_out_of_range_gt_and_sample_count_mismatch() {
    use std::io::Read;
    assert!(drain(bcf_bytes(&one("2/0"))).is_err());
    let bytes = bcf_bytes(&one("0/1"));
    let mut raw = Vec::new();
    flate2::read::MultiGzDecoder::new(&bytes[..])
        .read_to_end(&mut raw)
        .unwrap();
    let header_len = u32::from_le_bytes(raw[5..9].try_into().unwrap()) as usize;
    let header = String::from_utf8(raw[9..9 + header_len].to_vec()).unwrap();
    let header = header.replace("./sample A\n", "./sample A\tother\n");
    let mut modified = raw[..5].to_vec();
    modified.extend((header.len() as u32).to_le_bytes());
    modified.extend(header.as_bytes());
    modified.extend(&raw[9 + header_len..]);
    let mut input = VariantReader::from_reader(Cursor::new(modified)).unwrap();
    assert_eq!(input.header().sample_names().len(), 2);
    assert!(input.next_variant().is_err());
}

#[test]
fn bcf_missing_gt_at_locus_is_not_reference() {
    let text = data("S", "chr1\t1\t.\tA\tC\t.\tPASS\t.\tDP\t8\n");
    let mut input = VariantReader::from_reader(Cursor::new(bcf_bytes(&text))).unwrap();
    let variant = input.next_variant().unwrap().unwrap();
    assert!(variant.genotypes[0].is_none());
}

#[test]
fn config_validates_names_thresholds_and_reserved_flags() {
    let samples = vec!["./A".into(), "B,C".into()];
    let cli = Cli::try_parse_from([
        "app",
        "-i",
        "x",
        "--inspect",
        "-m",
        "4",
        "-o",
        "B,C",
        "-o",
        "./A",
    ])
    .unwrap();
    let config = cli.validate_samples(&samples).unwrap();
    assert_eq!(config.minimum_samples, 2);
    assert_eq!(config.outgroup_indices, [1, 0]);
    assert_eq!(config.warnings.len(), 1);
    for args in [vec!["-o", "A"], vec!["-o", "./A", "-o", "./A"]] {
        let cli =
            Cli::try_parse_from(["app", "-i", "x", "--inspect"].into_iter().chain(args)).unwrap();
        assert!(cli.validate_samples(&samples).is_err());
    }
    for n in ["0", "-1", "no", "184467440737095516160"] {
        assert!(Cli::try_parse_from(["app", "-i", "x", "--inspect", "-m", n]).is_err());
    }
    for flag in ["-p", "-f", "-r", "-w", "--resolve-IUPAC"] {
        let cli = Cli::try_parse_from(["app", "-i", "x", "--inspect", flag]).unwrap();
        assert!(cli.validate_mode().is_err());
    }
}

#[test]
fn actual_cli_inspection_never_creates_output() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.wrong-extension");
    std::fs::write(&input, gzip(one("0/1").as_bytes())).unwrap();
    let bin = env!("CARGO_BIN_EXE_vcf2phylip-rs");
    let run = |args: &[&str]| {
        Command::new(bin)
            .arg("-i")
            .arg(&input)
            .args(args)
            .current_dir(dir.path())
            .output()
            .unwrap()
    };
    let output = run(&["--inspect", "-m", "4"]);
    assert!(output.status.success(), "{:?}", output);
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("Inspection only"));
    assert!(stdout.contains("records: 1"));
    assert!(String::from_utf8(output.stderr).unwrap().contains("capped"));
    assert!(!run(&["--inspect", "--fasta"]).status.success());
    let seeded = run(&[
        "--inspect",
        "-r",
        "--seed",
        "42",
        "--partial-missing",
        "called-alleles",
    ]);
    assert!(seeded.status.success());
    let text = String::from_utf8(seeded.stdout).unwrap();
    assert!(text.contains("Accepted nucleotide loci: 1; accepted binary loci: 1"));
    assert!(text.contains("random seed: Some(42)"));
    assert!(!run(&["--inspect", "--seed", "42"]).status.success());
    assert!(
        !run(&["--inspect", "--output-folder", "must-not-exist"])
            .status
            .success()
    );
    assert!(!dir.path().join("must-not-exist").exists());
    assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 1);
    for flag in ["--help", "--version", "-v"] {
        assert!(
            Command::new(bin)
                .arg(flag)
                .output()
                .unwrap()
                .status
                .success()
        );
    }
}

fn conversion_fixture() -> String {
    // Uniform within-record ploidy avoids noodles' test-fixture BCF encoder
    // inserting vector-end padding after every allele of a shorter genotype.
    data("long_sample_name:1\tB,C\tΩ", "chr1\t9\t.\ta\tc\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0\nchr2\t2\t.\tA\tAT\t.\tPASS\t.\tGT\t0/0\t0/1\t0/0\nchr1\t1\t.\tA\tC,G,T\t.\tPASS\t.\tGT\t0/1/2\t0/./.\t3/3/3\nchr1\t1\t.\tG\tT\t.\tPASS\t.\tGT\t1/1/1\t0/1/1\t././.\nchr2\t5\t.\tG\t<DEL>\t.\tPASS\t.\tGT\t0\t0\t0\n")
        .replace("##contig=<ID=chr1,length=1000>\n", "##contig=<ID=chr1,length=1000>\n##contig=<ID=chr2,length=1000>\n")
}

#[test]
fn disk_memory_identical_across_tiles_transports_and_policies() {
    use vcf2phylip_rs::{conversion::convert, storage::StorageMode};
    let records: String = (0..163)
        .map(|i| {
            format!(
                "chr1\t{}\t.\tA\tC\t.\tPASS\t.\tGT\t{}\n",
                i + 1,
                ["0|1\t1/1\t./.", "0/1/1\t0/./.\t1|1|0"][i % 2]
            )
        })
        .collect();
    let text = data("A\tB\tC", &records);
    for bytes in [
        text.as_bytes().to_vec(),
        gzip(text.as_bytes()),
        bgzf(text.as_bytes()),
        bcf_bytes(&text),
    ] {
        for random in [false, true] {
            let dir = tempfile::tempdir().unwrap();
            let spool = tempfile::tempdir().unwrap();
            let input = dir.path().join("input");
            std::fs::write(&input, &bytes).unwrap();
            let mut expected = None;
            for budget in [
                None,
                Some(3),
                Some(6),
                Some(21),
                Some(240),
                Some(489),
                Some(600),
            ] {
                let folder = dir.path().join(format!("out{budget:?}"));
                let mut args = vec![
                    "app",
                    "-i",
                    input.to_str().unwrap(),
                    "--output-folder",
                    folder.to_str().unwrap(),
                    "--output-prefix",
                    "test",
                    "-m",
                    "1",
                    "-f",
                    "-w",
                    "-o",
                    "C",
                    "--partial-missing",
                    "called-alleles",
                ];
                if random {
                    args.extend(["-r", "--seed", "42"]);
                }
                let mut cli = Cli::try_parse_from(args).unwrap();
                if let Some(budget) = budget {
                    cli.storage = StorageMode::Disk;
                    cli.memory_limit = Some(budget);
                    cli.temp_dir = Some(spool.path().into());
                }
                let report = convert(&cli).unwrap();
                let outputs = report
                    .outputs
                    .iter()
                    .map(|p| std::fs::read(p).unwrap())
                    .collect::<Vec<_>>();
                let actual = (
                    outputs,
                    format!("{:?}", report.stats),
                    report.matrix_payload_bytes,
                );
                if let Some(expected) = &expected {
                    assert_eq!(expected, &actual);
                } else {
                    expected = Some(actual);
                }
                assert_eq!(std::fs::read_dir(spool.path()).unwrap().count(), 0);
                assert!(
                    convert(&cli)
                        .unwrap_err()
                        .to_string()
                        .contains("existing output")
                );
                assert_eq!(std::fs::read_dir(spool.path()).unwrap().count(), 0);
            }
        }
    }
}

#[test]
fn disk_errors_cleanup_and_budget_validation() {
    use vcf2phylip_rs::{conversion::convert, storage::StorageMode};
    for text in [
        data("S", ""),
        data("S", "chr1\t1\t.\tA\tAT\t.\tPASS\t.\tGT\t0\n"),
        format!(
            "{}bad late record\n",
            data("S", "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0\n")
        ),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let spool = tempfile::tempdir().unwrap();
        let input = dir.path().join("input");
        std::fs::write(&input, text).unwrap();
        let folder = dir.path().join("out");
        let mut cli = Cli::try_parse_from([
            "app",
            "-i",
            input.to_str().unwrap(),
            "--output-folder",
            folder.to_str().unwrap(),
            "-f",
            "-w",
        ])
        .unwrap();
        cli.storage = StorageMode::Disk;
        cli.memory_limit = Some(1);
        cli.temp_dir = Some(spool.path().into());
        assert!(convert(&cli).is_err());
        assert_eq!(std::fs::read_dir(spool.path()).unwrap().count(), 0);
        assert_eq!(std::fs::read_dir(&folder).unwrap().count(), 0);
        cli.temp_dir = Some(spool.path().join("nonexistent"));
        assert!(
            convert(&cli)
                .unwrap_err()
                .to_string()
                .contains("creating matrix spool")
        );
    }
    use vcf2phylip_rs::storage::DiskMatrix;
    let dir = tempfile::tempdir().unwrap();
    assert!(DiskMatrix::new(vec!["A".into(), "B".into()], 1, Some(dir.path())).is_err());
    assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    for value in ["0", "-1", "1GB", "1.5MiB", "18446744073709551616GiB"] {
        assert!(
            Cli::try_parse_from([
                "app",
                "-i",
                "x",
                "--storage",
                "disk",
                "--memory-limit",
                value
            ])
            .is_err()
        );
    }
    for value in ["1", "1B", "1KiB", "2MiB", "1GiB"] {
        assert!(
            Cli::try_parse_from([
                "app",
                "-i",
                "x",
                "--storage",
                "disk",
                "--memory-limit",
                value
            ])
            .is_ok()
        );
    }
    for flags in [
        vec!["--memory-limit", "1KiB"],
        vec!["--temp-dir", "/unused"],
        vec!["--inspect", "--storage", "disk"],
    ] {
        let cli = Cli::try_parse_from(["app", "-i", "unused"].into_iter().chain(flags)).unwrap();
        assert!(cli.validate_mode().is_err());
    }
}

#[test]
fn conversion_transport_equivalence_and_hand_computed_outputs() {
    use vcf2phylip_rs::conversion::convert;
    let text = conversion_fixture();
    for bytes in [
        text.as_bytes().to_vec(),
        gzip(text.as_bytes()),
        bgzf(text.as_bytes()),
        bcf_bytes(&text),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("example.vcf.middle.VCF.GZ");
        std::fs::write(&input, bytes).unwrap();
        let cli = Cli::try_parse_from([
            "app",
            "-i",
            input.to_str().unwrap(),
            "--output-folder",
            dir.path().to_str().unwrap(),
            "-m",
            "2",
            "-f",
            "-w",
            "-o",
            "B,C",
        ])
        .unwrap();
        let report = convert(&cli).unwrap();
        assert_eq!(report.stats.records_processed, 5);
        assert_eq!(report.stats.rejected_non_snp, 2);
        assert_eq!(report.stats.nucleotide_sites_accepted, 3);
        assert_eq!(report.matrix_payload_bytes, 9);
        let base = dir.path().join("example.vcf.middle.min2");
        let read =
            |ext: &str| std::fs::read_to_string(format!("{}.{ext}", base.display())).unwrap();
        assert_eq!(
            read("phy"),
            "3 3\nB,C  CNK\nlong_sample_name:1  MVT\nΩ  ATN\n"
        );
        assert_eq!(
            read("fasta"),
            ">B,C\nCNK\n>long_sample_name:1\nMVT\n>Ω\nATN\n"
        );
        assert_eq!(
            read("used_sites.tsv"),
            "#CHROM\tPOS\tNUM_SAMPLES\nchr1\t9\t3\nchr1\t1\t2\nchr1\t1\t2\n"
        );
        assert!(
            convert(&cli)
                .unwrap_err()
                .to_string()
                .contains("existing output")
        );
        assert_eq!(
            read("phy"),
            "3 3\nB,C  CNK\nlong_sample_name:1  MVT\nΩ  ATN\n"
        );
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 4);
    }
}

#[test]
fn conversion_seed_order_selection_and_partial_policy() {
    use vcf2phylip_rs::conversion::convert;
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.vcf");
    std::fs::write(&input, conversion_fixture()).unwrap();
    let mut baseline = None;
    for (prefix, flags) in [
        ("both", vec!["-f"]),
        ("fasta", vec!["-p", "-f", "-o", "Ω", "-o", "B,C"]),
    ] {
        let cli = Cli::try_parse_from(
            [
                "app",
                "-i",
                input.to_str().unwrap(),
                "--output-folder",
                dir.path().to_str().unwrap(),
                "--output-prefix",
                prefix,
                "-m",
                "1",
                "-r",
                "--seed",
                "42",
            ]
            .into_iter()
            .chain(flags),
        )
        .unwrap();
        convert(&cli).unwrap();
        let fasta =
            std::fs::read_to_string(dir.path().join(format!("{prefix}.min1.fasta"))).unwrap();
        let lines: Vec<_> = fasta.lines().collect();
        let by_name: std::collections::BTreeMap<_, _> = lines
            .chunks(2)
            .map(|c| (c[0].to_owned(), c[1].to_owned()))
            .collect();
        if let Some(expected) = &baseline {
            assert_eq!(&by_name, expected);
        } else {
            baseline = Some(by_name);
        }
        if prefix == "fasta" {
            assert!(fasta.starts_with(">Ω\n"));
            assert!(!dir.path().join("fasta.min1.phy").exists());
        }
    }
    let cli = Cli::try_parse_from([
        "app",
        "-i",
        input.to_str().unwrap(),
        "--output-folder",
        dir.path().to_str().unwrap(),
        "--output-prefix",
        "partial",
        "-m",
        "3",
        "--partial-missing",
        "called-alleles",
    ])
    .unwrap();
    let report = convert(&cli).unwrap();
    assert_eq!(report.stats.nucleotide_sites_accepted, 2);
    assert_eq!(
        std::fs::read_to_string(dir.path().join("partial.min3.phy")).unwrap(),
        "3 2\nlong_sample_name:1  MV\nB,C  CA\nΩ  AT\n"
    );
}

#[test]
fn conversion_failures_leave_no_final_outputs() {
    use vcf2phylip_rs::conversion::convert;
    for (text, flags) in [
        (data("S", ""), vec![]),
        (data("S", "chr1\t1\t.\tA\tAT\t.\tPASS\t.\tGT\t0\n"), vec![]),
        (
            format!(
                "{}bad late record\n",
                data("S", "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0\n")
            ),
            vec![],
        ),
        (one("0/1"), vec![]), // whitespace identifier: inspect permits it, conversion does not
        (data("S", ""), vec!["-p"]),
        (data("S", ""), vec!["-o", "unknown"]),
        (data("S", ""), vec!["--output-prefix", "../escape"]),
        (data("S", ""), vec!["--output-prefix", "/absolute"]),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.vcf");
        std::fs::write(&input, text).unwrap();
        let folder = dir.path().join("out");
        let cli = Cli::try_parse_from(
            [
                "app",
                "-i",
                input.to_str().unwrap(),
                "--output-folder",
                folder.to_str().unwrap(),
                "-f",
                "-w",
            ]
            .into_iter()
            .chain(flags.clone()),
        )
        .unwrap();
        assert!(convert(&cli).is_err(), "{flags:?}");
        if folder.exists() {
            assert_eq!(std::fs::read_dir(folder).unwrap().count(), 0);
        }
    }
}

#[test]
fn conversion_no_clobber_links_collisions_and_bad_directory() {
    use vcf2phylip_rs::conversion::convert;
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.vcf");
    let text = data("S", "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n");
    std::fs::write(&input, &text).unwrap();
    let cli = || {
        Cli::try_parse_from([
            "app",
            "-i",
            input.to_str().unwrap(),
            "--output-folder",
            dir.path().to_str().unwrap(),
            "-f",
        ])
        .unwrap()
    };
    let target = dir.path().join("in.min1.phy");
    std::fs::hard_link(&input, &target).unwrap();
    assert!(convert(&cli()).is_err());
    assert_eq!(std::fs::read_to_string(&input).unwrap(), text);
    std::fs::remove_file(&target).unwrap();
    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(dir.path().join("missing"), &target).unwrap();
        assert!(convert(&cli()).is_err());
        assert!(
            std::fs::symlink_metadata(&target)
                .unwrap()
                .file_type()
                .is_symlink()
        );
        std::fs::remove_file(&target).unwrap();
    }
    std::fs::create_dir(&target).unwrap();
    assert!(convert(&cli()).is_err());
    let mut config = cli();
    config.output_folder = input.join("impossible");
    assert!(convert(&config).is_err());
    config.input = dir.path().join("absent");
    assert!(convert(&config).is_err());
    assert!(!dir.path().join("in.min1.fasta").exists());

    // The input itself can have an alignment-like name; explicit prefixes must
    // not allow the converter to replace it even though detection is by bytes.
    let collision = dir.path().join("collision.min1.phy");
    std::fs::write(&collision, &text).unwrap();
    let mut config = cli();
    config.input = collision.clone();
    config.output_prefix = Some("collision".into());
    assert!(convert(&config).is_err());
    assert_eq!(std::fs::read_to_string(collision).unwrap(), text);
}

#[test]
fn actual_cli_default_conversion_and_unsupported_mode_preflight() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.VCF");
    std::fs::write(&input, data("S", "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n")).unwrap();
    let run = |flags: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_vcf2phylip-rs"))
            .arg("-i")
            .arg(&input)
            .args(flags)
            .current_dir(dir.path())
            .output()
            .unwrap()
    };
    assert!(
        !run(&["-n", "--output-folder", "never-created"])
            .status
            .success()
    );
    assert!(!dir.path().join("never-created").exists());
    let output = run(&[]);
    assert!(output.status.success(), "{output:?}");
    assert_eq!(
        std::fs::read_to_string(dir.path().join("in.min1.phy")).unwrap(),
        "1 1\nS  M\n"
    );
    let output = run(&[]);
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    assert!(
        String::from_utf8(output.stderr)
            .unwrap()
            .contains("existing output")
    );
}

#[test]
fn reader_to_processing_transport_equivalence_and_original_record_order() {
    use vcf2phylip_rs::{
        genotype::{NucleotideResolution, PartialMissingPolicy},
        processing::{Stats, VariantProcessor},
    };
    let text = data("S\tT", "chr1\t9\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t1/1\nchr2\t2\t.\tA\tAT\t.\tPASS\t.\tGT\t0/0\t0/1\nchr1\t1\t.\tA\tC,G,T\t.\tPASS\t.\tGT\t0/1/2/3\t0/.\nchr2\t1\t.\tG\tT\t.\tPASS\t.\tDP\t4\t2\n")
        .replace("##contig=<ID=chr1,length=1000>\n", "##contig=<ID=chr1,length=1000>\n##contig=<ID=chr2,length=1000>\n");
    let processor = VariantProcessor::new(
        2,
        Some(1),
        PartialMissingPolicy::CalledAlleles,
        NucleotideResolution::Random { seed: 9 },
    )
    .unwrap();
    let mut baseline = None;
    for bytes in [
        text.as_bytes().to_vec(),
        gzip(text.as_bytes()),
        bgzf(text.as_bytes()),
        bcf_bytes(&text),
    ] {
        let mut input = VariantReader::from_reader(Cursor::new(bytes)).unwrap();
        let mut results = Vec::new();
        let mut locations = Vec::new();
        let mut stats = Stats::default();
        while let Some(v) = input.next_variant().unwrap() {
            locations.push((
                v.record.reference_sequence_name().to_string(),
                usize::from(v.record.variant_start().unwrap()),
            ));
            let locus = processor.process(&v, input.records_read() - 1).unwrap();
            stats.observe(&locus);
            results.push(locus);
        }
        assert_eq!(
            locations,
            [
                ("chr1".into(), 9),
                ("chr2".into(), 2),
                ("chr1".into(), 1),
                ("chr2".into(), 1)
            ]
        );
        assert_eq!(stats.records_processed, 4);
        assert_eq!(stats.rejected_non_snp, 1);
        assert_eq!(stats.nucleotide_sites_accepted, 2);
        assert_eq!(stats.binary_eligible, 2);
        assert_eq!(stats.binary_sites_accepted, 1);
        if let Some(expected) = &baseline {
            assert_eq!(&results, expected);
        } else {
            baseline = Some(results);
        }
    }
}

#[test]
fn multidigit_out_of_range_is_error_in_vcf_and_bcf() {
    let text = data(
        "S",
        "chr1\t1\t.\tA\tC,G,T,C,G,T,C,G,T,<NON_REF>\t.\tPASS\t.\tGT\t11/0\n",
    );
    for bytes in [text.as_bytes().to_vec(), bcf_bytes(&text)] {
        let mut input = VariantReader::from_reader(Cursor::new(bytes)).unwrap();
        assert!(input.next_variant().is_err());
    }
}
