//! Preflight for the pinned noodles BCF decoder. Keep in sync with decoder upgrades.
use std::io::{BufRead, Read};

use anyhow::{Context, Result, bail, ensure};
use noodles::vcf;

/// Limits apply to decompressed BCF, independently of alignment storage limits.
/// These bound framing and decoder amplification, not total process memory.
#[derive(Clone, Copy, Debug)]
pub struct BcfLimits {
    pub header_bytes: usize,
    pub record_bytes: usize,
    pub dictionary_entries: usize,
    pub sample_format_cells: usize,
    pub typed_elements: usize,
}

impl Default for BcfLimits {
    fn default() -> Self {
        Self {
            header_bytes: 16 * 1024 * 1024,
            record_bytes: 64 * 1024 * 1024,
            dictionary_entries: 1_000_000,
            sample_format_cells: 1_000_000,
            typed_elements: 4_000_000,
        }
    }
}

fn word(input: &mut impl Read) -> Result<usize> {
    let mut bytes = [0; 4];
    input.read_exact(&mut bytes)?;
    Ok(usize::try_from(u32::from_le_bytes(bytes))?)
}

// Read incrementally: a truncated declaration does not cause an upfront allocation.
fn bounded_bytes(input: &mut impl Read, len: usize) -> Result<Vec<u8>> {
    let mut data = Vec::new();
    let mut chunk = [0; 8192];
    while data.len() < len {
        let n = (len - data.len()).min(chunk.len());
        input.read_exact(&mut chunk[..n])?;
        data.try_reserve(n).context("reserving BCF buffer")?;
        data.extend_from_slice(&chunk[..n]);
    }
    Ok(data)
}

pub(super) fn read_header(input: &mut impl Read, limits: BcfLimits) -> Result<vcf::Header> {
    use vcf::header::{StringMaps, parser::Entry};

    let mut magic = [0; 5];
    input.read_exact(&mut magic)?;
    ensure!(magic == *b"BCF\x02\x02", "unsupported BCF version");
    let len = word(input)?;
    ensure!(len <= limits.header_bytes, "BCF header byte limit exceeded");
    let data = bounded_bytes(input, len)?;
    let end = data
        .iter()
        .position(|&b| b == 0)
        .context("BCF header missing NUL terminator")?;
    let mut parser = vcf::header::Parser::default();
    let mut maps = StringMaps::default();
    let mut string_slots = 1usize; // implicit PASS
    let mut contig_slots = 0usize;
    for line in data[..end].split(|&b| b == b'\n').filter(|s| !s.is_empty()) {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        let entry = parser.parse_partial(line)?;
        let slot = match &entry {
            Entry::Contig(_, m) => Some((&mut contig_slots, m.idx())),
            Entry::Filter(_, m) => Some((&mut string_slots, m.idx())),
            Entry::Format(_, m) => Some((&mut string_slots, m.idx())),
            Entry::Info(_, m) => Some((&mut string_slots, m.idx())),
            _ => None,
        };
        if let Some((slots, idx)) = slot {
            // Check before insert_entry can resize its sparse vector (or overflow i + 1).
            let next = idx
                .unwrap_or(*slots)
                .checked_add(1)
                .context("BCF dictionary index overflow")?;
            *slots = (*slots).max(next);
            ensure!(
                *slots <= limits.dictionary_entries,
                "BCF dictionary entry limit exceeded"
            );
        }
        maps.insert_entry(&entry)?;
    }
    let mut header = parser.finish()?;
    *header.string_maps_mut() = maps;
    Ok(header)
}

pub(super) fn read_frame(
    input: &mut impl BufRead,
    header: &vcf::Header,
    limits: BcfLimits,
) -> Result<Vec<u8>> {
    let shared = word(input).context("BCF shared length")?;
    let individual = word(input).context("BCF individual length")?;
    let len = shared
        .checked_add(individual)
        .context("BCF record length overflow")?;
    ensure!(len <= limits.record_bytes, "BCF record byte limit exceeded");
    ensure!(shared >= 24, "BCF shared block shorter than fixed fields");
    let data = bounded_bytes(input, len).context("truncated BCF record")?;
    validate(&data[..shared], &data[shared..], header, limits)?;
    let mut frame = Vec::new();
    frame.try_reserve_exact(len.checked_add(8).context("BCF frame length overflow")?)?;
    frame.extend_from_slice(&(shared as u32).to_le_bytes());
    frame.extend_from_slice(&(individual as u32).to_le_bytes());
    frame.extend_from_slice(&data);
    Ok(frame)
}

fn take<'a>(src: &mut &'a [u8], n: usize) -> Result<&'a [u8]> {
    ensure!(n <= src.len(), "truncated BCF typed value");
    let (value, rest) = src.split_at(n);
    *src = rest;
    Ok(value)
}

fn integer(ty: u8, bytes: &[u8]) -> Result<i64> {
    Ok(match ty {
        1 => i64::from(bytes[0] as i8),
        2 => i64::from(i16::from_le_bytes(bytes.try_into()?)),
        3 => i64::from(i32::from_le_bytes(bytes.try_into()?)),
        _ => bail!("BCF expected integer"),
    })
}

fn width(ty: u8) -> Result<usize> {
    match ty {
        0 => Ok(0),
        1 | 7 => Ok(1),
        2 => Ok(2),
        3 | 5 => Ok(4),
        _ => bail!("invalid BCF type {ty}"),
    }
}

fn descriptor(src: &mut &[u8]) -> Result<(u8, usize)> {
    let b = take(src, 1)?[0];
    let ty = b & 15;
    width(ty)?;
    let mut len = usize::from(b >> 4);
    if len == 15 {
        // Extended lengths must be scalar integers, never recursively typed arrays.
        let b = take(src, 1)?[0];
        ensure!(
            b >> 4 == 1 && matches!(b & 15, 1..=3),
            "invalid BCF extended length"
        );
        let t = b & 15;
        len = usize::try_from(integer(t, take(src, width(t)?)?)?)?;
    }
    ensure!(ty != 0 || len == 0, "invalid BCF null length");
    Ok((ty, len))
}

struct Typed<'a> {
    ty: u8,
    len: usize,
    data: &'a [u8],
}

fn typed<'a>(src: &mut &'a [u8], samples: usize, budget: &mut usize) -> Result<Typed<'a>> {
    let (ty, len) = descriptor(src)?;
    let count = len
        .checked_mul(samples)
        .context("BCF typed length overflow")?;
    *budget = budget
        .checked_sub(count)
        .context("BCF typed element limit exceeded")?;
    let size = count
        .checked_mul(width(ty)?)
        .context("BCF typed byte length overflow")?;
    let data = take(src, size)?;
    Ok(Typed { ty, len, data })
}

fn key<'a>(src: &mut &[u8], header: &'a vcf::Header, budget: &mut usize) -> Result<&'a str> {
    let value = typed(src, 1, budget)?;
    ensure!(
        value.len == 1 && matches!(value.ty, 1..=3),
        "invalid BCF dictionary key"
    );
    let index = usize::try_from(integer(value.ty, value.data)?)?;
    header
        .string_maps()
        .strings()
        .get_index(index)
        .context("invalid BCF dictionary key index")
}

fn numeric(value: &Typed<'_>, allow_end: bool) -> Result<()> {
    if !matches!(value.ty, 1..=3 | 5) {
        return Ok(());
    }
    for bytes in value.data.chunks_exact(width(value.ty)?) {
        let (reserved, end) = if value.ty == 5 {
            let bits = u32::from_le_bytes(bytes.try_into()?);
            (
                (0x7f800003..=0x7f800007).contains(&bits),
                bits == 0x7f800002,
            )
        } else {
            let n = integer(value.ty, bytes)?;
            let min = match value.ty {
                1 => i8::MIN as i64,
                2 => i16::MIN as i64,
                _ => i32::MIN as i64,
            };
            ((min + 2..=min + 7).contains(&n), n == min + 1)
        };
        ensure!(
            !reserved && (allow_end || !end),
            "invalid BCF numeric sentinel"
        );
    }
    Ok(())
}

fn validate(
    shared: &[u8],
    mut individual: &[u8],
    header: &vcf::Header,
    limits: BcfLimits,
) -> Result<()> {
    use vcf::header::record::value::map::format::{Number, Type};
    let infos = usize::from(u16::from_le_bytes(shared[16..18].try_into()?));
    let alleles = usize::from(u16::from_le_bytes(shared[18..20].try_into()?));
    let counts = u32::from_le_bytes(shared[20..24].try_into()?);
    let formats = (counts >> 24) as usize;
    let samples = (counts & 0xffffff) as usize;
    ensure!(
        samples == header.sample_names().len(),
        "BCF record sample count differs from header"
    );
    ensure!(
        samples
            .checked_mul(formats.max(1))
            .context("BCF sample count overflow")?
            <= limits.sample_format_cells,
        "BCF sample/FORMAT cell limit exceeded"
    );
    let mut budget = limits.typed_elements;
    let mut site = &shared[24..];
    for _ in 0..1 + alleles {
        // ID and each allele
        typed(&mut site, 1, &mut budget)?;
    }
    let filters = typed(&mut site, 1, &mut budget)?;
    numeric(&filters, false)?;
    for _ in 0..infos {
        key(&mut site, header, &mut budget)?;
        let value = typed(&mut site, 1, &mut budget)?;
        numeric(&value, false)?;
    }
    ensure!(site.is_empty(), "trailing BCF shared data");
    let mut seen = std::collections::HashSet::new();
    for i in 0..formats {
        let name = key(&mut individual, header, &mut budget)?;
        ensure!(seen.insert(name), "duplicate BCF FORMAT key");
        let value = typed(&mut individual, samples, &mut budget)?;
        ensure!(
            value.len > 0 && value.ty != 0,
            "BCF FORMAT {name:?}: null or empty value type"
        );
        if name == "GT" {
            ensure!(i == 0, "GT must be the first FORMAT key");
            ensure!(
                value.ty == 1,
                "unsupported BCF GT encoding: only Int8 is supported (use VCF for Int16/Int32 GT)"
            );
            // Negative non-EOV GT codes can overflow noodles' i8 allele arithmetic.
            ensure!(
                value.data.iter().all(|&n| n <= 127 || n == 129),
                "invalid BCF GT code"
            );
        } else {
            let definition = header
                .formats()
                .get(name)
                .context("undefined BCF FORMAT key")?;
            ensure!(
                definition.number() != Number::Count(0),
                "BCF FORMAT {name:?}: invalid Number=0"
            );
            let compatible = match definition.ty() {
                Type::Integer => matches!(value.ty, 1..=3),
                Type::Float => value.ty == 5,
                Type::Character | Type::String => value.ty == 7,
            };
            ensure!(compatible, "BCF FORMAT {name:?}: type mismatch");
            numeric(
                &value,
                definition.number() != Number::Count(1) || value.len != 1,
            )?;
            if definition.ty() == Type::Character {
                for bytes in value.data.chunks_exact(value.len) {
                    let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
                    let text = std::str::from_utf8(&bytes[..end])?;
                    ensure!(
                        text.split(',').all(|s| !s.is_empty()),
                        "BCF FORMAT {name:?}: empty character"
                    );
                }
            }
        }
    }
    // The pinned decoder ignores bytes after the declared FORMAT fields. Retain
    // that compatibility (its writer can emit excess mixed-ploidy GT padding).
    // These bytes remain frame-bounded and are never interpreted by the decoder.
    Ok(())
}
