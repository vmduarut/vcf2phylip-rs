//! One growable sequence per sample, stored in original header order.
use anyhow::{Result, ensure};

use crate::filter::ProcessedLocus;

pub struct Matrix {
    samples: Vec<String>,
    sequences: Vec<Vec<u8>>,
    sites: usize,
}

impl Matrix {
    pub fn new(samples: Vec<String>) -> Result<Self> {
        ensure!(!samples.is_empty(), "input has no samples");
        let mut seen = std::collections::HashSet::new();
        for name in &samples {
            ensure!(
                !name.is_empty() && !name.chars().any(|c| c.is_whitespace() || c.is_control()),
                "sample name cannot be represented as an unambiguous FASTA/relaxed PHYLIP identifier: {name:?}"
            );
            ensure!(seen.insert(name), "duplicate sample name: {name:?}");
        }
        Ok(Self {
            sequences: vec![Vec::new(); samples.len()],
            samples,
            sites: 0,
        })
    }

    pub fn append(&mut self, locus: &ProcessedLocus) -> Result<()> {
        if !locus.nucleotide_accepted {
            return Ok(());
        }
        ensure!(
            locus.nucleotide.len() == self.samples.len(),
            "matrix/sample dimension mismatch"
        );
        let sites = self
            .sites
            .checked_add(1)
            .ok_or_else(|| anyhow::anyhow!("matrix length overflow"))?;
        sites
            .checked_mul(self.samples.len())
            .ok_or_else(|| anyhow::anyhow!("matrix size overflow"))?;
        for sequence in &mut self.sequences {
            sequence.try_reserve(1)?;
        }
        for (sequence, call) in self.sequences.iter_mut().zip(&locus.nucleotide) {
            sequence.push(call.symbol());
        }
        self.sites = sites;
        Ok(())
    }

    pub fn samples(&self) -> &[String] {
        &self.samples
    }
    pub fn sequences(&self) -> &[Vec<u8>] {
        &self.sequences
    }
    pub fn sites(&self) -> usize {
        self.sites
    }
    pub fn payload_bytes(&self) -> usize {
        self.sites * self.samples.len()
    }

    pub(crate) fn validate_order(&self, order: &[usize]) -> Result<()> {
        let mut seen = vec![false; self.samples.len()];
        ensure!(
            order.len() == seen.len(),
            "output sample order dimension mismatch"
        );
        for &i in order {
            ensure!(
                i < seen.len() && !seen[i],
                "output order is not a permutation"
            );
            seen[i] = true;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reject_ambiguous_identifiers_without_remapping() {
        for name in ["", "a b", "a\tb", "a\rb", "a\nb", "a\0b", "a\u{2003}b"] {
            assert!(Matrix::new(vec![name.into()]).is_err(), "{name:?}");
        }
        assert!(Matrix::new(vec!["same".into(), "same".into()]).is_err());
        let name = "./long:sample,Ω";
        let matrix = Matrix::new(vec![name.into()]).unwrap();
        assert_eq!(matrix.samples()[0].as_bytes(), name.as_bytes());
    }
}
