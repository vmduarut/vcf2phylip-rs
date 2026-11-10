//! Canonical nucleotide sets, independent of dosage and phase.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct BaseMask(u8);

impl BaseMask {
    /// Only nonempty subsets of A=1, C=2, G=4, T=8 are valid.
    pub fn new(bits: u8) -> Option<Self> {
        (1..=15).contains(&bits).then_some(Self(bits))
    }

    pub fn bits(self) -> u8 {
        self.0
    }

    pub fn symbol(self) -> u8 {
        b"?ACMGRSVTWYHKDBN"[usize::from(self.0)]
    }
}

pub fn canonical_base(allele: &str) -> Option<u8> {
    if allele.len() != 1 {
        return None;
    }
    match allele.as_bytes()[0].to_ascii_uppercase() {
        b'A' => Some(1),
        b'C' => Some(2),
        b'G' => Some(4),
        b'T' => Some(8),
        _ => None,
    }
}
