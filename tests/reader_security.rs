use std::io::{Cursor, Write};

use flate2::{Compression, write::GzEncoder};
use vcf2phylip_rs::reader::{BcfLimits, VariantReader};

const META: &str = "##fileformat=VCFv4.3\n##FILTER=<ID=PASS,Description=\"PASS\",IDX=0>\n##contig=<ID=chr1,length=1000,IDX=0>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\",IDX=1>\n";
const COLS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n";

fn header(meta: &str) -> Vec<u8> {
    let text = format!("{meta}{COLS}\0");
    let mut data = b"BCF\x02\x02".to_vec();
    data.extend_from_slice(&(text.len() as u32).to_le_bytes());
    data.extend_from_slice(text.as_bytes());
    data
}

fn bcf(meta: &str, individual: &[u8], formats: u8) -> Vec<u8> {
    let mut shared = Vec::new();
    for n in [0u32, 0, 1, 0x7f800001] {
        shared.extend_from_slice(&n.to_le_bytes());
    }
    shared.extend_from_slice(&0u16.to_le_bytes());
    shared.extend_from_slice(&2u16.to_le_bytes());
    shared.extend_from_slice(&((u32::from(formats) << 24) | 1).to_le_bytes());
    shared.extend_from_slice(&[7, 0x17, b'A', 0x17, b'C', 0x11, 0]);
    let mut data = header(meta);
    data.extend_from_slice(&(shared.len() as u32).to_le_bytes());
    data.extend_from_slice(&(individual.len() as u32).to_le_bytes());
    data.extend_from_slice(&shared);
    data.extend_from_slice(individual);
    data
}

fn gzip(data: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
    encoder.write_all(data).unwrap();
    encoder.finish().unwrap()
}

fn error(data: Vec<u8>) -> String {
    match VariantReader::from_reader(Cursor::new(data)) {
        Err(e) => format!("{e:#}"),
        Ok(mut reader) => {
            let e = reader.next_variant().unwrap_err();
            assert!(reader.next_variant().unwrap().is_none());
            format!("{e:#}")
        }
    }
}

#[test]
fn rejects_header_indices_before_sparse_allocation_or_overflow() {
    for idx in ["100000000", "18446744073709551615"] {
        let data = header(&META.replace("IDX=1", &format!("IDX={idx}")));
        for data in [data.clone(), gzip(&data)] {
            let message = error(data);
            assert!(message.contains("dictionary"), "{message}");
        }
    }
}

#[test]
fn rejects_large_lengths_before_reading_or_allocating_payload() {
    let mut data = b"BCF\x02\x02".to_vec();
    data.extend_from_slice(&u32::MAX.to_le_bytes());
    assert!(error(data).contains("header byte limit"));
    for (shared, individual) in [(1u32 << 30, 1), (31, 1u32 << 30), (u32::MAX, u32::MAX)] {
        let mut data = header(META);
        data.extend_from_slice(&shared.to_le_bytes());
        data.extend_from_slice(&individual.to_le_bytes());
        assert!(error(data).contains("record byte limit"));
    }
}

#[test]
fn rejects_gt_and_auxiliary_panic_payloads_plain_and_gzip() {
    let dp = format!("{META}##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\",IDX=2>\n");
    let character =
        format!("{META}##FORMAT=<ID=X,Number=1,Type=Character,Description=\"X\",IDX=2>\n");
    let cases = [
        bcf(META, &[0x11, 1, 0x22, 2, 0, 4, 0], 1), // valid but unsupported Int16 GT
        bcf(META, &[0x11, 1, 0x15, 0, 0, 0x80, 0x3f], 1),
        bcf(META, &[0x11, 1, 0], 1),
        bcf(META, &[0x11, 1, 0x11, 128], 1),
        bcf(&dp, &[0x11, 1, 0x21, 2, 4, 0x11, 2, 0], 2),
        bcf(&dp, &[0x11, 1, 0x21, 2, 4, 0x11, 2, 0x17, b'1'], 2),
        bcf(&dp, &[0x11, 1, 0x21, 2, 4, 0x11, 2, 0x11, 129], 2),
        bcf(&character, &[0x11, 1, 0x21, 2, 4, 0x11, 2, 0x17, 0], 2),
        // Recursive extended length must not reach noodles' recursive descriptor reader.
        bcf(META, &[0x11, 1, 0xf1, 0xf1, 0xf1], 1),
        bcf(META, &[0x11, 1, 0xf1, 0x13, 0xff, 0xff, 0xff, 0x7f], 1),
    ];
    for data in cases {
        for data in [data.clone(), gzip(&data)] {
            let message = error(data);
            assert!(message.contains("BCF"), "{message}");
        }
    }
}

#[test]
fn bounded_configuration_and_valid_common_records() {
    let dp = format!("{META}##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\",IDX=2>\n");
    let data = bcf(&dp, &[0x11, 1, 0x21, 2, 4, 0x11, 2, 0x12, 42, 0], 2);
    for data in [data.clone(), gzip(&data)] {
        let mut reader = VariantReader::from_reader(Cursor::new(data)).unwrap();
        let record = reader.next_variant().unwrap().unwrap();
        assert_eq!(record.genotypes[0].as_ref().unwrap().as_ref().len(), 2);
        assert!(reader.next_variant().unwrap().is_none());
    }
    let limits = BcfLimits {
        sample_format_cells: 1,
        ..BcfLimits::default()
    };
    let mut reader = VariantReader::from_reader_with_bcf_limits(Cursor::new(data), limits).unwrap();
    assert!(format!("{:#}", reader.next_variant().unwrap_err()).contains("cell limit"));
}

#[test]
fn truncated_bcf_never_succeeds() {
    let data = bcf(META, &[0x11, 1, 0x21, 2, 4], 1);
    let header_len = header(META).len();
    for n in (header_len + 1)..data.len() {
        assert!(error(data[..n].to_vec()).contains("BCF"));
    }
}

#[test]
fn single_byte_record_mutations_return_without_panicking() {
    let base = bcf(META, &[0x11, 1, 0x21, 2, 4], 1);
    // Small explicit limits make this deterministic corpus safe even if mutations
    // turn fixed fields into large sample counts, lengths, or typed descriptors.
    let limits = BcfLimits {
        header_bytes: 4096,
        record_bytes: 4096,
        dictionary_entries: 256,
        sample_format_cells: 256,
        typed_elements: 4096,
    };
    for offset in header(META).len()..base.len() {
        for byte in 0..=255 {
            let mut data = base.clone();
            data[offset] = byte;
            let mut reader =
                VariantReader::from_reader_with_bcf_limits(Cursor::new(data), limits).unwrap();
            while let Ok(Some(_)) = reader.next_variant() {}
        }
    }
}

#[test]
fn reader_identifier_contexts_escape_controls() {
    let data = format!("{META}{COLS}\x1b[2JFORGED\tbad\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n");
    let message = error(data.into_bytes());
    assert!(message.contains("\\u{1b}"));
    assert!(!message.contains('\x1b'));
    let data = format!(
        "{META}{}chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1:1\n",
        COLS.replace("\tS\n", "\t\x1b[2JFORGED\n")
    );
    let message = error(data.into_bytes());
    assert!(message.contains("\\u{1b}"));
    assert!(!message.contains('\x1b'));
}
