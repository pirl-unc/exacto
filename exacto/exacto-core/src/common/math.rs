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


use statrs::distribution::{Discrete, Hypergeometric};
use statrs::function::beta::ln_beta;
use statrs::function::factorial::ln_binomial;


pub fn f1_score(precision: f64, recall: f64) -> f64 {
    let denom = precision + recall;
    if denom == 0.0 {
        0.0
    } else {
        2.0 * (precision * recall) / denom
    }
}


/// Two-sided Fisher exact test p-value for a 2x2 table.
///
/// Table:
///             FWD    REV
/// ALT          a      b
/// REF          c      d
///
/// Returns p-value in [0, 1].
pub fn fisher_exact_test_two_sided(
    a: u64,
    b: u64,
    c: u64,
    d: u64
) -> f64 {
    let r1: u64 = a + b;    // ALT total
    let r2: u64 = c + d;    // REF total
    let c1: u64 = a + c;    // FWD total
    let c2: u64 = b + d;    // REV total
    let n: u64 = r1 + r2;   // total

    // Population: n
    // Population successes: c1 (FWD reads)
    // Draws: r1 (ALT reads)
    let hg: Hypergeometric = Hypergeometric::new(n, c1, r1).expect("Invalid hypergeometric parameters for Fisher test");
    let observed_p: f64 = hg.pmf(a);

    // Feasible range:
    // max(0, r1 - c2) <= x <= min(r1, c1)
    let min_x: u64 = r1.saturating_sub(n - c1);
    let max_x: u64 = r1.min(c1);

    // Two-sided definition: sum probabilities <= observed probability
    let mut p_two_sided: f64 = 0.0;
    for x in min_x..=max_x {
        let px: f64 = hg.pmf(x);
        if px <= observed_p + 1e-12 {
            p_two_sided += px;
        }
    }

    p_two_sided.min(1.0)
}
