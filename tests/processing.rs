use std::io::Cursor;

use vcf2phylip_rs::{
    genotype::{
        NucleotideCall, NucleotideResolution, PartialMissingPolicy as Policy, nucleotide_call,
    },
    iupac::BaseMask,
    processing::{Rejection, Stats, VariantProcessor, classify_snp},
    reader::{Variant, VariantReader},
};

fn variant(reference: &str, alt: &str, calls: &[&str]) -> Variant {
    let samples = (0..calls.len())
        .map(|i| format!("S{i}"))
        .collect::<Vec<_>>()
        .join("\t");
    let text = format!(
        "##fileformat=VCFv4.3\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\nchr1\t1\t.\t{reference}\t{alt}\t.\tPASS\t.\tGT\t{}\n",
        calls.join("\t")
    );
    VariantReader::from_reader(Cursor::new(text.into_bytes()))
        .unwrap()
        .next_variant()
        .unwrap()
        .unwrap()
}

fn processor(
    samples: usize,
    minimum: usize,
    policy: Policy,
    resolution: NucleotideResolution,
) -> VariantProcessor {
    VariantProcessor::new(samples, Some(minimum), policy, resolution).unwrap()
}

#[test]
fn all_fifteen_masks_and_missing_are_independently_expected() {
    let expected = [
        (1, 'A'),
        (2, 'C'),
        (3, 'M'),
        (4, 'G'),
        (5, 'R'),
        (6, 'S'),
        (7, 'V'),
        (8, 'T'),
        (9, 'W'),
        (10, 'Y'),
        (11, 'H'),
        (12, 'K'),
        (13, 'D'),
        (14, 'B'),
        (15, 'N'),
    ];
    for (bits, symbol) in expected {
        let mask = BaseMask::new(bits).unwrap();
        assert_eq!(mask.bits(), bits);
        assert_eq!(mask.symbol(), symbol as u8);
        let gt = (0..4)
            .filter(|i| bits & (1 << i) != 0)
            .map(|i| i.to_string())
            .collect::<Vec<_>>()
            .join("/");
        let v = variant("A", "C,G,T", &[&gt]);
        let result = processor(1, 1, Policy::Missing, NucleotideResolution::Iupac)
            .process(&v, 0)
            .unwrap();
        assert_eq!(result.nucleotide, [NucleotideCall::Called(mask)]);
        assert_eq!(result.nucleotide_called_samples, 1);
        assert!(result.nucleotide_accepted);
    }
    assert!(BaseMask::new(0).is_none());
    assert!(BaseMask::new(16).is_none());
    assert_eq!(NucleotideCall::Missing.symbol(), b'N');
    assert_ne!(
        NucleotideCall::Missing,
        NucleotideCall::Called(BaseMask::new(15).unwrap())
    );
}

#[test]
fn phase_dosage_and_positive_ploidies() {
    let v = variant("a", "c", &["0", "0/0", "0|0|0", "0/1", "1|0", "1/0/1/1"]);
    let result = processor(6, 1, Policy::Missing, NucleotideResolution::Iupac)
        .process(&v, 0)
        .unwrap();
    assert_eq!(
        result
            .nucleotide
            .iter()
            .map(|c| c.symbol())
            .collect::<Vec<_>>(),
        b"AAAMMM"
    );
    assert_eq!(result.nucleotide_called_samples, 6);
}

#[test]
fn partial_missing_policy_is_symmetric() {
    let v = variant("A", "C", &["0/.", "./0", "./.", ".", "1/.", "0/1"]);
    for (policy, expected, count) in [
        (Policy::Missing, b"NNNNNM", 1),
        (Policy::CalledAlleles, b"AANNCM", 4),
    ] {
        let result = processor(6, 2, policy, NucleotideResolution::Iupac)
            .process(&v, 0)
            .unwrap();
        let symbols: Vec<_> = result.nucleotide.iter().map(|c| c.symbol()).collect();
        let expected: Vec<_> = expected.iter().copied().filter(|b| *b != b' ').collect();
        assert_eq!(symbols, expected);
        assert_eq!(result.nucleotide_called_samples, count);
        assert_eq!(result.nucleotide_accepted, count >= 2);
    }
}

#[test]
fn missing_calls_and_positive_ploidies_have_nucleotide_counts() {
    let v = variant(
        "G",
        "T",
        &[
            "0/0", "0/1", "1|0", "1/1", "0", "0/1/1", "0/.", ".", "././.",
        ],
    );
    let result = processor(9, 4, Policy::CalledAlleles, NucleotideResolution::Iupac)
        .process(&v, 1)
        .unwrap();
    assert_eq!(
        result
            .nucleotide
            .iter()
            .map(|c| c.symbol())
            .collect::<Vec<_>>(),
        b"GKKTGKGNN"
    );
    assert_eq!(result.nucleotide_called_samples, 7);
    assert!(result.nucleotide_accepted);
    let mut stats = Stats::default();
    stats.observe(&result);
    assert_eq!(stats.records_processed, 1);
    assert_eq!(stats.nucleotide_missing_calls, 2);
    assert_eq!(stats.nucleotide_sites_accepted, 1);
    assert_eq!(stats.nucleotide_below_minimum(), 0);
}

#[test]
fn biological_classification_is_case_insensitive_and_not_parse_failure() {
    let strings = |xs: &[&str]| xs.iter().map(|s| s.to_string()).collect::<Vec<_>>();
    assert_eq!(
        classify_snp("a", &strings(&["c", "G", "t"])).unwrap(),
        [1, 2, 4, 8]
    );
    assert_eq!(classify_snp("A", &[]), Err(Rejection::NoAlternate));
    for (reference, alts) in [
        ("N", vec!["C"]),
        ("AA", vec!["CC"]),
        ("A", vec!["."]),
        ("A", vec!["N"]),
        ("A", vec!["*"]),
        ("A", vec!["<NON_REF>"]),
        ("A", vec!["C", "<DEL>"]),
        ("A", vec!["C", "AT"]),
    ] {
        assert_eq!(
            classify_snp(reference, &strings(&alts)),
            Err(Rejection::NonCanonicalAllele)
        );
        let v = variant(reference, &alts.join(","), &["0/0"]);
        let result = processor(1, 1, Policy::Missing, NucleotideResolution::Iupac)
            .process(&v, 0)
            .unwrap();
        assert!(result.rejection.is_some());
        assert!(!result.nucleotide_accepted);
        assert!(result.nucleotide.is_empty());
    }
    for alts in [vec!["a"], vec!["C", "c"], vec!["C", "A"]] {
        assert_eq!(
            classify_snp("A", &strings(&alts)),
            Err(Rejection::DuplicateAllele)
        );
    }
}

#[test]
fn multiallelic_sites_keep_one_nucleotide_locus() {
    let v = variant("A", "C,G,T", &["0/1/2/3", "1|2", "3"]);
    let result = processor(3, 3, Policy::Missing, NucleotideResolution::Iupac)
        .process(&v, 20)
        .unwrap();
    assert_eq!(
        result
            .nucleotide
            .iter()
            .map(|c| c.symbol())
            .collect::<Vec<_>>(),
        b"NST"
    );
    assert!(result.nucleotide_accepted);
    assert_eq!(result.nucleotide_called_samples, 3);
}

#[test]
fn genotype_validation_precedes_biological_rejection() {
    let mut rejected = variant("A", "AT", &["0/0"]);
    let out_of_range = variant("A", "C,G", &["2/0"]);
    rejected.genotypes[0] = out_of_range.genotypes[0].clone();
    let processor = processor(1, 1, Policy::Missing, NucleotideResolution::Iupac);
    assert!(processor.process(&rejected, 0).is_err());
    rejected.genotypes[0] = Some(Default::default());
    assert!(processor.process(&rejected, 0).is_err());
}

#[test]
fn validated_dimensions_and_thresholds() {
    assert!(VariantProcessor::new(0, None, Policy::Missing, NucleotideResolution::Iupac).is_err());
    assert!(
        VariantProcessor::new(2, Some(0), Policy::Missing, NucleotideResolution::Iupac).is_err()
    );
    assert_eq!(
        VariantProcessor::new(2, None, Policy::Missing, NucleotideResolution::Iupac)
            .unwrap()
            .minimum_samples(),
        2
    );
    assert_eq!(
        processor(7, 99, Policy::Missing, NucleotideResolution::Iupac).minimum_samples(),
        7
    );
    let v = variant("A", "C", &["0"]);
    assert!(
        processor(2, 1, Policy::Missing, NucleotideResolution::Iupac)
            .process(&v, 0)
            .is_err()
    );
    assert!(
        nucleotide_call(
            v.genotypes[0].as_ref(),
            &[],
            Policy::Missing,
            NucleotideResolution::Iupac,
            0,
            0
        )
        .is_err()
    );
}

#[test]
fn random_resolution_is_keyed_copy_weighted_and_preserves_called_counts() {
    let v = variant(
        "A",
        "C,G,T",
        &["0/0/1", "1|2", "0/1/2/3", "3", "0/.", "./1", "."],
    );
    let processor = processor(
        7,
        1,
        Policy::CalledAlleles,
        NucleotideResolution::Random { seed: 42 },
    );
    let mut golden = Vec::new();
    for record_ordinal in 0..4 {
        let result = processor.process(&v, record_ordinal).unwrap();
        assert_eq!(result, processor.process(&v, record_ordinal).unwrap());
        for (i, call) in result.nucleotide.iter().enumerate() {
            assert!(
                [b"AC".as_slice(), b"CG", b"ACGT", b"T", b"A", b"C", b"N"][i]
                    .contains(&call.symbol())
            );
        }
        golden.extend(result.nucleotide.iter().map(|c| c.symbol()));
        // Direct calls and processing use identical keys.
        for (sample, gt) in v.genotypes.iter().enumerate() {
            assert_eq!(
                result.nucleotide[sample],
                nucleotide_call(
                    gt.as_ref(),
                    &[1, 2, 4, 8],
                    Policy::CalledAlleles,
                    NucleotideResolution::Random { seed: 42 },
                    record_ordinal,
                    sample
                )
                .unwrap()
            );
        }
    }
    assert_eq!(
        String::from_utf8(golden).unwrap(),
        "AGTTACNAGGTACNACTTACNAGATACN"
    );
    let biallelic = variant("A", "C", &["0/0/1", "0/1", "1/1", "0/."]);
    let plain = self::processor(4, 1, Policy::Missing, NucleotideResolution::Iupac)
        .process(&biallelic, 0)
        .unwrap();
    let random = self::processor(
        4,
        1,
        Policy::Missing,
        NucleotideResolution::Random { seed: 42 },
    )
    .process(&biallelic, 0)
    .unwrap();
    assert_eq!(plain.nucleotide_accepted, random.nucleotide_accepted);
    assert_eq!(
        plain.nucleotide_called_samples,
        random.nucleotide_called_samples
    );
    assert_eq!(random.nucleotide[3], NucleotideCall::Missing);
    // Three copies preserve 2:1 dosage rather than a uniform choice of unique bases.
    let mut a = 0;
    for record_ordinal in 0..600 {
        a += usize::from(
            processor.process(&v, record_ordinal).unwrap().nucleotide[0].symbol() == b'A',
        );
    }
    assert!((350..450).contains(&a), "A copies: {a}");
}
