//! Per-record scientific eligibility and independent nucleotide/binary thresholds.
use anyhow::{Result, ensure};

use crate::{
    genotype::{self, BinaryCall, NucleotideCall, PartialMissingPolicy, Resolution},
    iupac::canonical_base,
    reader::Variant,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Rejection {
    NoAlternate,
    NonCanonicalAllele,
    DuplicateAllele,
}

/// Returns an index-preserving canonical base map, or a biological rejection.
/// Unsupported biology is not a malformed-record error.
pub fn classify_snp(reference: &str, alternates: &[String]) -> Result<Vec<u8>, Rejection> {
    if alternates.is_empty() {
        return Err(Rejection::NoAlternate);
    }
    let mut bases = Vec::with_capacity(alternates.len() + 1);
    let mut seen = 0;
    for allele in std::iter::once(reference).chain(alternates.iter().map(String::as_str)) {
        let base = canonical_base(allele).ok_or(Rejection::NonCanonicalAllele)?;
        if seen & base != 0 {
            return Err(Rejection::DuplicateAllele);
        }
        seen |= base;
        bases.push(base);
    }
    Ok(bases)
}

#[derive(Clone, Debug)]
pub struct ProcessingConfig {
    sample_count: usize,
    minimum_samples: usize,
    partial: PartialMissingPolicy,
    resolution: Resolution,
}

impl ProcessingConfig {
    /// Low-level count-only constructor; `None` retains the legacy minimum of 4.
    /// CLI callers pass the centralized effective count from `validate_samples`,
    /// combining the absolute minimum and call rate. A capped minimum is exposed
    /// via `minimum_samples`; callers should warn if it differs from the request.
    pub fn new(
        sample_count: usize,
        minimum: Option<usize>,
        partial: PartialMissingPolicy,
        resolution: Resolution,
    ) -> Result<Self> {
        let minimum = minimum.unwrap_or(4);
        ensure!(sample_count > 0, "input has no samples");
        ensure!(minimum > 0, "minimum samples must be positive");
        u64::try_from(sample_count)?;
        Ok(Self {
            sample_count,
            minimum_samples: minimum.min(sample_count),
            partial,
            resolution,
        })
    }

    pub fn minimum_samples(&self) -> usize {
        self.minimum_samples
    }

    /// The ordinal includes rejected records. No sorting, splitting or reordering.
    pub fn process(&self, variant: &Variant, ordinal: u64) -> Result<ProcessedLocus> {
        ensure!(
            variant.genotypes.len() == self.sample_count,
            "genotype/sample dimension mismatch"
        );
        // Validate all GT indices even if the biological site will be rejected.
        let allele_count = variant.record.alternate_bases().as_ref().len() + 1;
        for gt in variant.genotypes.iter().flatten() {
            ensure!(!gt.as_ref().is_empty(), "zero-ploidy genotype");
            ensure!(
                gt.as_ref()
                    .iter()
                    .all(|a| a.position().is_none_or(|i| i < allele_count)),
                "GT allele index out of range"
            );
        }
        let bases = match classify_snp(
            variant.record.reference_bases(),
            variant.record.alternate_bases().as_ref(),
        ) {
            Ok(bases) => bases,
            Err(reason) => {
                return Ok(ProcessedLocus {
                    rejection: Some(reason),
                    nucleotide: vec![],
                    binary: None,
                    nucleotide_called: 0,
                    binary_called: 0,
                    nucleotide_accepted: false,
                    binary_accepted: false,
                });
            }
        };
        let nucleotide = variant
            .genotypes
            .iter()
            .enumerate()
            .map(|(i, gt)| {
                genotype::nucleotide_call(
                    gt.as_ref(),
                    &bases,
                    self.partial,
                    self.resolution,
                    ordinal,
                    i,
                )
            })
            .collect::<Result<Vec<_>>>()?;
        let nucleotide_called = nucleotide.iter().filter(|call| call.is_called()).count();
        let binary = (bases.len() == 2)
            .then(|| {
                variant
                    .genotypes
                    .iter()
                    .map(|gt| genotype::binary_call(gt.as_ref()))
                    .collect::<Result<Vec<_>>>()
            })
            .transpose()?;
        let binary_called = binary.as_ref().map_or(0, |calls| {
            calls
                .iter()
                .filter(|c| matches!(c, BinaryCall::Called(_)))
                .count()
        });
        Ok(ProcessedLocus {
            rejection: None,
            nucleotide,
            binary,
            nucleotide_called,
            binary_called,
            nucleotide_accepted: nucleotide_called >= self.minimum_samples,
            binary_accepted: binary_called >= self.minimum_samples,
        })
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ProcessedLocus {
    pub rejection: Option<Rejection>,
    pub nucleotide: Vec<NucleotideCall>,
    /// None means biologically ineligible; '?' calls occur only at eligible sites.
    pub binary: Option<Vec<BinaryCall>>,
    pub nucleotide_called: usize,
    pub binary_called: usize,
    pub nucleotide_accepted: bool,
    pub binary_accepted: bool,
}

#[derive(Debug, Default)]
pub struct Stats {
    pub records: u64,
    pub rejected_non_snp: u64,
    pub nucleotide_accepted: u64,
    pub binary_eligible: u64,
    pub binary_accepted: u64,
    pub nucleotide_missing_calls: u64,
    pub binary_missing_calls: u64,
    pub binary_unsupported_ploidy_calls: u64,
}

impl Stats {
    /// SNP-eligible records rejected by the nucleotide called-sample threshold.
    ///
    /// For loci from `ProcessingConfig::process`, biological rejection returns
    /// before calling samples and always sets nucleotide acceptance to false.
    /// Otherwise acceptance is exactly the called-sample threshold comparison.
    /// Thus rejected, below-minimum and accepted records partition `records`.
    /// Requires counters accumulated from those loci, without manual mutation.
    pub fn nucleotide_below_minimum(&self) -> u64 {
        self.records - self.rejected_non_snp - self.nucleotide_accepted
    }

    /// Call exactly once per processed record, including biological rejections.
    pub fn observe(&mut self, locus: &ProcessedLocus) {
        self.records += 1;
        self.rejected_non_snp += u64::from(locus.rejection.is_some());
        self.nucleotide_accepted += u64::from(locus.nucleotide_accepted);
        self.binary_eligible += u64::from(locus.binary.is_some());
        self.binary_accepted += u64::from(locus.binary_accepted);
        self.nucleotide_missing_calls +=
            locus.nucleotide.iter().filter(|c| !c.is_called()).count() as u64;
        if let Some(calls) = &locus.binary {
            self.binary_missing_calls += calls
                .iter()
                .filter(|c| matches!(c, BinaryCall::Missing))
                .count() as u64;
            self.binary_unsupported_ploidy_calls += calls
                .iter()
                .filter(|c| matches!(c, BinaryCall::UnsupportedPloidy))
                .count() as u64;
        }
    }
}
