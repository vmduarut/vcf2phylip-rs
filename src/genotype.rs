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
pub enum NucleotideResolution {
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

/// `allele_base_masks` contains one single-bit canonical mask per REF/ALT index.
/// Keys must be the original zero-based record ordinal and sample index.
pub fn nucleotide_call(
    gt: Option<&Genotype>,
    allele_base_masks: &[u8],
    policy: PartialMissingPolicy,
    resolution: NucleotideResolution,
    record_ordinal: u64,
    sample: usize,
) -> Result<NucleotideCall> {
    ensure!(
        allele_base_masks.iter().all(|b| [1, 2, 4, 8].contains(b)),
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
            ensure!(
                index < allele_base_masks.len(),
                "GT allele index {index} out of range"
            );
            copies.push(allele_base_masks[index]);
        }
    }
    if copies.is_empty()
        || (copies.len() != alleles.len() && policy == PartialMissingPolicy::Missing)
    {
        return Ok(NucleotideCall::Missing);
    }
    let bits = match resolution {
        NucleotideResolution::Iupac => copies.iter().fold(0, |mask, b| mask | b),
        NucleotideResolution::Random { seed } => {
            copies[choose_copy(seed, record_ordinal, sample, copies.len())?]
        }
    };
    Ok(NucleotideCall::Called(
        BaseMask::new(bits).context("invalid nucleotide mask")?,
    ))
}

// Stable v1 scheme: ChaCha20 key = seed LE || record ordinal LE || 16 zero bytes;
// stream = original sample index as u64; counter starts at zero for each call.
// Draw u64 words, reject [0, 2^64 mod n), then reduce modulo n. This samples
// copies (not unique alleles), without modulo bias or shared mutable RNG state.
fn choose_copy(seed: u64, record_ordinal: u64, sample: usize, count: usize) -> Result<usize> {
    let sample = u64::try_from(sample)?;
    let n = u64::try_from(count)?;
    ensure!(n > 0, "cannot choose from zero copies");
    let mut key = [0; 32];
    key[..8].copy_from_slice(&seed.to_le_bytes());
    key[8..16].copy_from_slice(&record_ordinal.to_le_bytes());
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
