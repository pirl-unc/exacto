use super::*;


#[test]
fn dna_variant_read_support_index_returns_matches() {
    let min_read_support_index: DNAVariantReadSupportIndex = DNAVariantReadSupportIndex::new(
        60u32,
        30u32,
        0.5f64,
        0.001f64,
        0.01f64,
        0.02f64,
        1e-6f64
    );
    assert_eq!(min_read_support_index.get_min_read_support((0, 30)), 6);

    let min_read_support_index: DNAVariantReadSupportIndex = DNAVariantReadSupportIndex::new(
        60u32,
        30u32,
        0.5f64,
        0.001f64,
        0.01f64,
        0.02f64,
        1e-6f64
    );
    assert_eq!(min_read_support_index.get_min_read_support((4, 30)), 11);

    let min_read_support_index: DNAVariantReadSupportIndex = DNAVariantReadSupportIndex::new(
        60u32,
        30u32,
        0.5f64,
        0.001f64,
        0.01f64,
        0.02f64,
        1e-6f64
    );
    assert_eq!(min_read_support_index.get_min_read_support((5, 30)), 13);

    let min_read_support_index: DNAVariantReadSupportIndex = DNAVariantReadSupportIndex::new(
        60u32,
        30u32,
        0.5f64,
        0.001f64,
        0.01f64,
        0.02f64,
        1e-6f64
    );
    assert_eq!(min_read_support_index.get_min_read_support((6, 30)), 14);
}

#[test]
fn calculate_min_dna_read_support_returns_matches() {
    let threshold: VariantReadSupportThreshold = calculate_min_dna_read_support(
        28,
        0.5f64,
        1e-3,
        1e-2,
        1e-6f64
    );
    assert_eq!(threshold.min_read_support, 6);

    let threshold: VariantReadSupportThreshold = calculate_min_dna_read_support(
        30,
        0.5f64,
        1e-3,
        1e-2,
        1e-6f64
    );
    assert_eq!(threshold.min_read_support, 6);

    let threshold: VariantReadSupportThreshold = calculate_min_dna_read_support(
        60,
        0.25f64,
        1e-6,
        1e-2,
        1e-6f64
    );
    assert_eq!(threshold.min_read_support, 9);
}
#[test]
fn dna_variant_read_support_index_matches_a_search_over_every_threshold() {
    // The builder starts each depth at the first threshold within max_fpr and stops early. Here
    // every threshold of every depth is scored and the first k of the peak F1 kept, as the index
    // was built before. The threshold is not monotone in the depth: at the default error model it
    // falls by one read at depths of 444 (repeat length 7), 537 (8) and 625 (3, somatic), so a
    // search that started at the threshold of the depth before would miss it.
    for (allele_fraction, mutation_rate) in [(0.5f64, 1e-3f64), (0.25f64, 1e-6f64)] {
        let index: DNAVariantReadSupportIndex = DNAVariantReadSupportIndex::new(700, 30, allele_fraction, mutation_rate, 0.01, 0.03, 1e-6);
        for repeat_length in [0u32, 3, 4, 7, 8, 9] {
            let error_rate: f64 = if repeat_length == 0 {
                0.01
            } else {
                1.0 - (1.0 - 0.03f64).powf(repeat_length as f64 - 1.0)
            };
            for depth in 1..=700u32 {
                let binom_f: Binomial = Binomial::new(allele_fraction, depth as u64).unwrap();
                let binom_e: Binomial = Binomial::new(error_rate, depth as u64).unwrap();
                let mut peak: Option<(u32, f64)> = None;
                for k in 1..=depth as u64 {
                    let recall: f64 = 1.0 - binom_f.cdf(k - 1);
                    let fpr: f64 = 1.0 - binom_e.cdf(k - 1);
                    let precision: f64 = (recall * mutation_rate) / ((recall * mutation_rate) + (fpr * (1.0 - mutation_rate)));
                    let f1: f64 = f1_score(precision, recall);
                    if fpr <= 1e-6 && f1.is_finite() && peak.map_or(true, |(_, best)| f1 > best) {
                        peak = Some((k as u32, f1));
                    }
                }
                assert_eq!(
                    index.get_min_read_support((repeat_length, depth)),
                    peak.map_or(u32::MAX, |(k, _)| k),
                    "allele fraction {allele_fraction}, repeat length {repeat_length}, depth {depth}"
                );
            }
        }
    }
}
