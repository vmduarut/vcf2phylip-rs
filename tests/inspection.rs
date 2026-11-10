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
    let cli =
        Cli::try_parse_from(["app", "-i", "x", "--inspect", "-o", "B,C", "-o", "./A"]).unwrap();
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
    for flag in ["-p", "-f", "-n", "-b", "-r", "-w", "--resolve-IUPAC"] {
        let cli = Cli::try_parse_from(["app", "-i", "x", "--inspect", flag]).unwrap();
        assert!(cli.validate_mode().is_err());
    }
}

#[test]
fn actual_cli_never_fakes_conversion_or_creates_output() {
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
    let output = run(&["--inspect"]);
    assert!(output.status.success(), "{:?}", output);
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("Inspection only"));
    assert!(stdout.contains("records: 1"));
    assert!(String::from_utf8(output.stderr).unwrap().contains("capped"));
    assert!(!run(&[]).status.success());
    assert!(!run(&["--inspect", "--fasta"]).status.success());
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
