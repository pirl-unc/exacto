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
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator};
use rayon::ThreadPool;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use statrs::function::beta::ln_beta;
use statrs::function::factorial::ln_binomial;

use crate::prelude::ReadSupportIndex;


#[derive(Clone, Debug, PartialEq)]
pub struct RNAVariantReadSupportIndex {
    /// HashMap<number of repeats, HashMap<total depth, minimum read support>>
    index: HashMap<u32, HashMap<ReadDepth, ReadSupport>>
}

impl RNAVariantReadSupportIndex {
    pub fn new(
        depths: &HashSet<ReadDepth>,
        max_slippage_repeat_length: u32,
        expected_sequencing_error: f64,
        expected_slippage_rate: f64,
        max_fpr: f64,
        num_threads: usize
    ) -> RNAVariantReadSupportIndex {
        assert!(max_slippage_repeat_length > 1);
        assert!(expected_sequencing_error >= 0.0);
        assert!(expected_sequencing_error < 1.0);
        assert!(expected_slippage_rate >= 0.0);
        assert!(expected_slippage_rate < 1.0);
        assert!(max_fpr > 0.0);
        assert!(max_fpr <= 1.0);

        // Vec<(number of repeats, total depth, error rate)>
        let mut cells: Vec<(u32, ReadDepth, f64)> = Vec::new();

        // No repeat
        for c in depths.iter() {
            cells.push((0, *c, expected_sequencing_error));
        }

        // We do not index 1-bp homopolymer repeats.

        // Repeats (from 2-bp)
        for r in 2..=max_slippage_repeat_length {
            let p_slippage: f64 = 1.0 - (1.0 - expected_slippage_rate).powf(r as f64 - 1.0);
            for c in depths.iter() {
                cells.push((r, *c, p_slippage));
            }
        }

        let min_read_supports: Vec<u32> = calculate_min_rna_variant_read_supports(
            &cells.iter().map(|&(_, c, p)| (c as u64, p, p)).collect::<Vec<(u64, f64, f64)>>(),
            max_fpr,
            num_threads
        );

        let mut index: HashMap<u32, HashMap<ReadDepth, ReadSupport>> = HashMap::new();

        for (&(repeat_len, depth, _), min_read_support) in cells.iter().zip(min_read_supports) {
            index.entry(repeat_len)
                .or_insert_with(HashMap::new)
                .insert(depth, min_read_support as ReadSupport);
        }

        Self {
            index: index
        }
    }
}

impl ReadSupportIndex for RNAVariantReadSupportIndex {
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


pub(crate) fn calculate_min_rna_variant_read_supports(
    cells: &[(u64, f64, f64)],
    max_fpr: f64,
    num_threads: usize
) -> Vec<u32> {
    let calculate = |&(total_depth, rho, expected_sequencing_error): &(u64, f64, f64)| -> u32 {
        calculate_min_rna_variant_read_support(
            total_depth,
            rho,
            expected_sequencing_error,
            max_fpr
        ) as u32
    };
    if rayon::current_thread_index().is_some() {
        return cells.iter().map(calculate).collect();
    }
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("Failed to build Rayon thread pool");
    thread_pool.install(|| {
        cells
            .par_iter()
            .with_max_len(1)
            .map(calculate)
            .collect()
    })
}


fn calculate_min_rna_variant_read_support(
    total_depth: u64,
    rho: f64,
    expected_sequencing_error: f64,
    max_fpr: f64
) -> usize {
    assert!((0.0..1.0).contains(&rho));
    assert!((0.0..1.0).contains(&expected_sequencing_error));
    assert!(max_fpr > 0.0 && max_fpr <= 1.0);

    if total_depth < 2 {
        return 2;
    }

    let (alpha, beta) = eps_rho_to_alpha_beta(expected_sequencing_error, rho);

    let mut tail: f64 = 0.0;
    let mut k: u64 = total_depth;
    while k >= 1 {
        let term: f64 = betabinom_ln_pmf(k, total_depth, alpha, beta).exp();
        if term.is_nan() || tail + term >= max_fpr {
            return ((k + 1) as usize).max(2);
        }
        tail += term;
        k -= 1;
    }

    // Every count down to 1 cleared; the hard floor still demands 2.
    2
}


/// Log PMF of the beta-binomial distribution, `ln P(X = k)` for
/// `X ~ BetaBinomial(n, alpha, beta)`:
///
///   P(X = k) = C(n, k) * B(k + alpha, n - k + beta) / B(alpha, beta)
fn betabinom_ln_pmf(k: u64, n: u64, alpha: f64, beta: f64) -> f64 {
    ln_binomial(n, k)
        + ln_beta(k as f64 + alpha, (n - k) as f64 + beta)
        - ln_beta(alpha, beta)
}


fn eps_rho_to_alpha_beta(eps: f64, rho: f64) -> (f64, f64) {
    let s = (1.0 - rho) / rho;
    (eps * s, (1.0 - eps) * s)
}
