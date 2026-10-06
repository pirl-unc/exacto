// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


use exacto_core::prelude::*;
use statrs::distribution::{Binomial, DiscreteCDF};
use std::collections::HashMap;

use crate::prelude::*;
use crate::filtering::variant_read_support_threshold::VariantReadSupportThreshold;


#[derive(Clone, Debug, PartialEq)]
pub struct DNAVariantReadSupportIndex {
    /// HashMap<number of repeats, HashMap<total depth, minimum read support>>
    index: HashMap<u32, HashMap<ReadDepth, ReadSupport>>
}

impl DNAVariantReadSupportIndex {
    pub fn new(
        max_depth: ReadDepth,
        max_slippage_repeat_length: u32,
        expected_variant_allele_fraction: f64,
        expected_mutation_rate: f64,
        expected_sequencing_error: f64,
        expected_slippage_rate: f64,
        max_fpr: f64
    ) -> Self {
        assert!(max_depth > 0);
        assert!(max_slippage_repeat_length > 1);
        assert!(expected_variant_allele_fraction > 0.0);
        assert!(expected_variant_allele_fraction < 1.0);
        assert!(expected_mutation_rate > 0.0);
        assert!(expected_mutation_rate < 1.0);
        assert!(expected_sequencing_error >= 0.0);
        assert!(expected_sequencing_error < 1.0);
        assert!(expected_slippage_rate >= 0.0);
        assert!(expected_slippage_rate < 1.0);
        assert!(max_fpr > 0.0);
        assert!(max_fpr <= 1.0);

        let mut index: HashMap<u32, HashMap<ReadDepth, ReadSupport>> = HashMap::new();

        // No repeat.
        for depth in 1..=max_depth {
            let threshold: VariantReadSupportThreshold = calculate_min_dna_read_support(
                depth as u64,
                expected_variant_allele_fraction,
                expected_mutation_rate,
                expected_sequencing_error,
                max_fpr
            );
            index.entry(0)
                .or_insert_with(HashMap::new)
                .insert(depth, threshold.min_read_support);
        }

        // Repeats (from 2-bp).
        // No need to index 1-bp homopolymer repeats.
        for repeat_len in 2..=max_slippage_repeat_length {
            let p_slippage: f64 = 1.0 - (1.0 - expected_slippage_rate).powf(repeat_len as f64 - 1.0);
            for depth in 1..=max_depth {
                let threshold: VariantReadSupportThreshold = calculate_min_dna_read_support(
                    depth as u64,
                    expected_variant_allele_fraction,
                    expected_mutation_rate,
                    p_slippage,
                    max_fpr
                );
                index.entry(repeat_len)
                    .or_insert_with(HashMap::new)
                    .insert(depth, threshold.min_read_support);
            }
        }

        Self {
            index: index
        }
    }
}

impl ReadSupportIndex for DNAVariantReadSupportIndex {
    /// (Number of repeats, read depth)
    type Input = (u32, ReadDepth);
    
    fn get_min_read_support(&self, inputs: (u32, ReadDepth)) -> ReadSupport {
        let num_repeats: u32 = inputs.0;
        let depth: ReadDepth = inputs.1;
        *self
            .index
            .get(&num_repeats)
            .and_then(|by_depth| by_depth.get(&depth))
            .unwrap_or_else(|| panic!(
                "ReadSupportIndex has no entry for (num_repeats {num_repeats}, depth {depth})."
            ))
    }
}


fn calculate_min_dna_read_support(
    total_depth: u64,
    expected_variant_allele_fraction: f64,
    expected_mutation_rate: f64,
    expected_sequencing_error: f64,
    max_fpr: f64
) -> VariantReadSupportThreshold {
    assert!(total_depth > 0);
    assert!((0.0..=1.0).contains(&expected_variant_allele_fraction) && expected_variant_allele_fraction > 0.0);
    assert!((0.0..=1.0).contains(&expected_mutation_rate) && expected_mutation_rate > 0.0);
    assert!((0.0..=1.0).contains(&expected_sequencing_error));
    assert!((0.0..=1.0).contains(&max_fpr));

    // Set up binomial testing.
    let binom_f: Binomial = Binomial::new(expected_variant_allele_fraction, total_depth).unwrap();
    let binom_e: Binomial = Binomial::new(expected_sequencing_error, total_depth).unwrap();

    // The first k with FPR(k) <= max_fpr: FPR(low - 1) > max_fpr >= FPR(high).
    let fpr_at = |k: u64| -> f64 { 1.0 - binom_e.cdf(k - 1) };
    if fpr_at(total_depth) > max_fpr {
        return VariantReadSupportThreshold {
            min_read_support: u32::MAX,
            recall: 0.0,
            precision: 0.0,
            fpr: 0.0,
            f1: 0.0
        };
    }
    let (mut low, mut high): (u64, u64) = (1, total_depth);
    while low < high {
        let middle: u64 = (low + high) / 2;
        if fpr_at(middle) <= max_fpr {
            high = middle;
        } else {
            low = middle + 1;
        }
    }

    // Eligible k: FPR(k) <= max_fpr, with a finite F1. The first k of the peak F1 is kept.
    let mut best: Option<VariantReadSupportThreshold> = None;
    for k in high..=total_depth {
        let recall: f64 = 1.0 - binom_f.cdf(k - 1);
        let fpr: f64    = 1.0 - binom_e.cdf(k - 1);
        let precision: f64 = (recall * expected_mutation_rate) / ((recall * expected_mutation_rate) + (fpr * (1.0 - expected_mutation_rate)));
        let f1: f64 = f1_score(precision, recall);
        if fpr <= max_fpr
            && f1.is_finite()
            && best.map_or(true, |best| f1 > best.f1) {
            best = Some(VariantReadSupportThreshold {
                min_read_support: k as ReadSupport,
                fpr,
                recall,
                precision,
                f1
            });
        }
        if best.is_some_and(|best| f1_score(1.0, recall) <= best.f1) {
            break;
        }
    }

    best.unwrap_or(VariantReadSupportThreshold {
        min_read_support: u32::MAX,
        recall: 0.0f64,
        precision: 0.0f64,
        fpr: 0.0f64,
        f1: 0.0f64
    })
}


#[cfg(test)]
#[path = "../../tests/filtering/dna/variant_read_support_index.rs"]
mod tests;
