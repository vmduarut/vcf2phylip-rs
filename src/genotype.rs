//! Genotype interpretation. Missing calls and a called four-base `N` stay distinct.
use anyhow::{Context, Result, ensure};
use noodles::vcf::variant::record_buf::samples::sample::value::Genotype;
use rand_chacha::{
    ChaCha20Rng,
    rand_core::{Rng, SeedableRng},
};

use crate::iupac::BaseMask;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, clap::ValueEnum)]
pub enum PartialMissingPolicy {
    #[default]
    Missing,
    CalledAlleles,
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum Resolution {
    #[default]
    Iupac,
    Random {
        seed: u64,
    },
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NucleotideCall {
    Missing,
    Called(BaseMask),
}

impl NucleotideCall {
    pub fn symbol(self) -> u8 {
        match self {
            Self::Missing => b'N',
            Self::Called(mask) => mask.symbol(),
        }
    }

    pub fn is_called(self) -> bool {
        matches!(self, Self::Called(_))
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BinaryCall {
    Called(u8),
    Missing,
    UnsupportedPloidy,
}

impl BinaryCall {
    pub fn symbol(self) -> u8 {
        match self {
            Self::Called(n) => b'0' + n,
            _ => b'?',
        }
    }
}

/// `bases` contains one single-bit canonical mask per REF/ALT index.
/// Keys must be the original zero-based record ordinal and sample index.
pub fn nucleotide_call(
    gt: Option<&Genotype>,
    bases: &[u8],
    policy: PartialMissingPolicy,
    resolution: Resolution,
    locus: u64,
    sample: usize,
) -> Result<NucleotideCall> {
    ensure!(
        bases.iter().all(|b| [1, 2, 4, 8].contains(b)),
        "invalid canonical base map"
    );
    let Some(gt) = gt else {
        return Ok(NucleotideCall::Missing);
    };
    let alleles = gt.as_ref();
    ensure!(!alleles.is_empty(), "zero-ploidy genotype");
    let mut copies = Vec::with_capacity(alleles.len());
    for allele in alleles {
        if let Some(index) = allele.position() {
            ensure!(index < bases.len(), "GT allele index {index} out of range");
            copies.push(bases[index]);
        }
    }
    if copies.is_empty()
        || (copies.len() != alleles.len() && policy == PartialMissingPolicy::Missing)
    {
        return Ok(NucleotideCall::Missing);
    }
    let bits = match resolution {
        Resolution::Iupac => copies.iter().fold(0, |mask, b| mask | b),
        Resolution::Random { seed } => copies[choose_copy(seed, locus, sample, copies.len())?],
    };
    Ok(NucleotideCall::Called(
        BaseMask::new(bits).context("invalid nucleotide mask")?,
    ))
}

/// Biallelic sites only. Partial policy and nucleotide resolution never apply.
pub fn binary_call(gt: Option<&Genotype>) -> Result<BinaryCall> {
    let Some(gt) = gt else {
        return Ok(BinaryCall::Missing);
    };
    let alleles = gt.as_ref();
    ensure!(!alleles.is_empty(), "zero-ploidy genotype");
    for allele in alleles {
        ensure!(
            allele.position().is_none_or(|i| i < 2),
            "non-biallelic genotype index"
        );
    }
    if alleles.len() != 2 {
        return Ok(BinaryCall::UnsupportedPloidy);
    }
    match (alleles[0].position(), alleles[1].position()) {
        (Some(a), Some(b)) => Ok(BinaryCall::Called((a + b) as u8)),
        _ => Ok(BinaryCall::Missing),
    }
}

// Stable v1 scheme: ChaCha20 key = seed LE || locus LE || 16 zero bytes;
// stream = original sample index as u64; counter starts at zero for each call.
// Draw u64 words, reject [0, 2^64 mod n), then reduce modulo n. This samples
// copies (not unique alleles), without modulo bias or shared mutable RNG state.
fn choose_copy(seed: u64, locus: u64, sample: usize, count: usize) -> Result<usize> {
    let sample = u64::try_from(sample)?;
    let n = u64::try_from(count)?;
    ensure!(n > 0, "cannot choose from zero copies");
    let mut key = [0; 32];
    key[..8].copy_from_slice(&seed.to_le_bytes());
    key[8..16].copy_from_slice(&locus.to_le_bytes());
    let mut rng = ChaCha20Rng::from_seed(key);
    rng.set_stream(sample);
    let threshold = n.wrapping_neg() % n;
    loop {
        let word = rng.next_u64();
        if word >= threshold {
            return Ok(usize::try_from(word % n)?);
        }
    }
}
