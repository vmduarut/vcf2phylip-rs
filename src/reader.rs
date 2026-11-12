//! One-record-at-a-time adapters. No SNP filtering or nucleotide conversion occurs here.
use std::{
    collections::HashSet,
    fs::File,
    io::{BufRead, BufReader, Cursor},
    path::Path,
};

use anyhow::{Context, Result, bail, ensure};
use flate2::bufread::MultiGzDecoder;
use noodles::{bcf, vcf};
use vcf::variant::record_buf::samples::sample::{Value, value::Genotype};

mod bounded_bcf;
pub use bounded_bcf::BcfLimits;

type Input = Box<dyn BufRead>;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Format {
    Vcf,
    Bcf,
}

enum Adapter {
    Vcf(Input),
    Bcf(Input),
}

/// Missing GT (absent FORMAT key or `.`) is None: ploidy is unknown.
/// Present GT retains integer allele indices, per-allele missingness and phasing.
#[derive(Debug)]
pub struct Variant {
    pub record: vcf::variant::RecordBuf,
    pub genotypes: Vec<Option<Genotype>>,
}

pub struct VariantReader {
    adapter: Adapter,
    header: vcf::Header,
    format: Format,
    compressed: bool,
    records_read: u64,
    line: String,
    finished: bool,
    bcf_limits: BcfLimits,
}

impl VariantReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
        Self::from_reader(BufReader::new(file))
            .with_context(|| format!("reading {}", path.display()))
    }

    /// Detect by bytes, not filename. Supports ordinary/multimember gzip and BGZF.
    pub fn from_reader(reader: impl BufRead + 'static) -> Result<Self> {
        Self::from_reader_with_bcf_limits(reader, BcfLimits::default())
    }

    /// Override the BCF safety limits for trusted, unusually large datasets.
    pub fn from_reader_with_bcf_limits(
        reader: impl BufRead + 'static,
        bcf_limits: BcfLimits,
    ) -> Result<Self> {
        let mut input: Input = Box::new(reader);
        // BufRead does not promise any minimum fill size, so sniff using read_exact
        // and chain the consumed bytes back into the stream.
        let (prefix, restored) = sniff(input, 2)?;
        input = restored;
        let compressed = prefix == [0x1f, 0x8b];
        if compressed {
            input = Box::new(BufReader::new(MultiGzDecoder::new(input)));
        }
        let (prefix, restored) = sniff(input, 5).context("detecting decompressed format")?;
        input = restored;
        let (adapter, header, format) = if prefix.starts_with(b"BCF") {
            ensure!(
                prefix == b"BCF\x02\x02",
                "unsupported BCF version (expected 2.2)"
            );
            let header = bounded_bcf::read_header(&mut input, bcf_limits).context("BCF header")?;
            (Adapter::Bcf(input), header, Format::Bcf)
        } else {
            ensure!(
                prefix.starts_with(b"##fil"),
                "missing VCF fileformat header or unrecognized input format"
            );
            let mut reader = vcf::io::Reader::new(input);
            let header = reader.read_header().context("VCF header")?;
            (Adapter::Vcf(reader.into_inner()), header, Format::Vcf)
        };
        validate_header(&header)?;
        Ok(Self {
            adapter,
            header,
            format,
            compressed,
            records_read: 0,
            line: String::new(),
            finished: false,
            bcf_limits,
        })
    }

    pub fn header(&self) -> &vcf::Header {
        &self.header
    }
    pub fn format(&self) -> Format {
        self.format
    }
    pub fn is_compressed(&self) -> bool {
        self.compressed
    }
    pub fn records_read(&self) -> u64 {
        self.records_read
    }

    /// An error terminates this stream. EOF is only success after fully draining input.
    pub fn next_variant(&mut self) -> Result<Option<Variant>> {
        if self.finished {
            return Ok(None);
        }
        let result = self.read_next_variant().with_context(|| {
            format!(
                "record {} ({} records successfully read)",
                self.records_read + 1,
                self.records_read
            )
        });
        match &result {
            Ok(Some(_)) => self.records_read += 1,
            _ => self.finished = true,
        }
        result
    }

    fn read_next_variant(&mut self) -> Result<Option<Variant>> {
        let mut record = vcf::variant::RecordBuf::default();
        match &mut self.adapter {
            Adapter::Vcf(input) => {
                self.line.clear();
                if input.read_line(&mut self.line)? == 0 {
                    return Ok(None);
                }
                let line = self.line.strip_suffix('\n').unwrap_or(&self.line);
                let line = line.strip_suffix('\r').unwrap_or(line);
                let fields: Vec<_> = line.split('\t').collect();
                ensure!(
                    fields.len() == 9 + self.header.sample_names().len(),
                    "VCF field count: expected {}, got {}",
                    9 + self.header.sample_names().len(),
                    fields.len()
                );
                ensure!(fields.iter().all(|s| !s.is_empty()), "empty VCF field");
                let keys: Vec<_> = fields[8].split(':').collect();
                ensure!(keys.iter().all(|k| !k.is_empty()), "empty FORMAT key");
                let mut seen = HashSet::new();
                ensure!(keys.iter().all(|k| seen.insert(*k)), "duplicate FORMAT key");
                if let Some(i) = keys.iter().position(|k| *k == "GT") {
                    ensure!(i == 0, "GT must be the first FORMAT key");
                }
                for (i, sample) in fields[9..].iter().enumerate() {
                    let values: Vec<_> = sample.split(':').collect();
                    ensure!(
                        values.len() <= keys.len() && values.iter().all(|v| !v.is_empty()),
                        "sample {:?}: invalid FORMAT value count or empty value",
                        self.header.sample_names()[i]
                    );
                }
                vcf::io::Reader::new(Cursor::new(line.as_bytes()))
                    .read_record_buf(&self.header, &mut record)
                    .with_context(|| format!("VCF {:?}:{:?}", fields[0], fields[1]))?;
            }
            Adapter::Bcf(input) => {
                if input.fill_buf()?.is_empty() {
                    return Ok(None);
                }
                let frame = bounded_bcf::read_frame(input, &self.header, self.bcf_limits)?;
                let mut decoder = bcf::io::Reader::from(Cursor::new(frame));
                ensure!(
                    decoder.read_record_buf(&self.header, &mut record)? != 0,
                    "truncated BCF record length"
                );
            }
        }
        let context = format!(
            "{:?}:{:?}",
            record.reference_sequence_name(),
            record.variant_start()
        );
        validate_record(record, self.header.sample_names().len())
            .with_context(|| context)
            .map(Some)
    }
}

fn sniff(mut input: Input, n: usize) -> Result<(Vec<u8>, Input)> {
    use std::io::Read;
    let mut prefix = Vec::new();
    input.by_ref().take(n as u64).read_to_end(&mut prefix)?;
    let restored = BufReader::new(Cursor::new(prefix.clone()).chain(input));
    Ok((prefix, Box::new(restored)))
}

fn validate_header(header: &vcf::Header) -> Result<()> {
    use vcf::header::record::value::map::format::{Number, Type};
    ensure!(!header.sample_names().is_empty(), "header has zero samples");
    let mut seen = HashSet::new();
    for name in header.sample_names() {
        ensure!(!name.is_empty(), "empty sample name");
        ensure!(seen.insert(name), "duplicate sample name: {name:?}");
    }
    let gt = header
        .formats()
        .get("GT")
        .context("header must define FORMAT/GT (Number=1, Type=String)")?;
    ensure!(
        gt.number() == Number::Count(1) && gt.ty() == Type::String,
        "FORMAT/GT must have Number=1, Type=String"
    );
    Ok(())
}

fn validate_record(record: vcf::variant::RecordBuf, sample_count: usize) -> Result<Variant> {
    ensure!(record.variant_start().is_some(), "missing variant position");
    ensure!(
        !record.reference_sequence_name().is_empty(),
        "empty contig name"
    );
    ensure!(
        !record.reference_bases().is_empty(),
        "empty reference allele"
    );
    ensure!(
        record.samples().values().count() == sample_count,
        "record sample count differs from header"
    );
    let allele_count = record.alternate_bases().as_ref().len() + 1;
    let genotypes = record
        .samples()
        .values()
        .enumerate()
        .map(|(i, sample)| match sample.get("GT").flatten() {
            None => Ok(None),
            Some(Value::Genotype(gt)) => {
                ensure!(!gt.as_ref().is_empty(), "sample {}: zero-ploidy GT", i + 1);
                for allele in gt.as_ref() {
                    if let Some(index) = allele.position() {
                        ensure!(
                            index < allele_count,
                            "sample {}: GT allele index {index} out of range ({} alleles)",
                            i + 1,
                            allele_count
                        );
                    }
                }
                Ok(Some(gt.clone()))
            }
            Some(_) => bail!("sample {}: GT was not parsed as a genotype", i + 1),
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(Variant { record, genotypes })
}
