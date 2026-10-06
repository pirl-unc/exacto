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
use std::collections::{HashMap, HashSet};

use crate::prelude::*;
use crate::filtering::rna::variant_read_support_index::calculate_min_rna_variant_read_supports;


#[derive(Clone, Debug, PartialEq)]
pub struct RNAJunctionReadSupportIndex {
    index: HashMap<ReadDepth, ReadSupport>
}

impl RNAJunctionReadSupportIndex {
    pub fn new(
        depths: &HashSet<ReadDepth>,
        expected_sequencing_error: f64,
        max_fpr: f64,
        num_threads: usize
    ) -> RNAJunctionReadSupportIndex {
        assert!(expected_sequencing_error >= 0.0);
        assert!(expected_sequencing_error < 1.0);
        assert!(max_fpr > 0.0);
        assert!(max_fpr <= 1.0);

        // At a sequencing error of 0 the null leaves no read to errors, 
        // so every depth takes the floor of 2.
        let depths: Vec<ReadDepth> = depths.iter().copied().collect();
        let min_read_supports: Vec<ReadSupport> = if expected_sequencing_error == 0.0 {
            vec![2; depths.len()]
        } else {
            // Vec<(total depth, overdispersion, error rate)>
            let cells: Vec<(u64, f64, f64)> = depths
                .iter()
                .map(|&depth| (depth as u64, expected_sequencing_error, expected_sequencing_error))
                .collect();
            calculate_min_rna_variant_read_supports(&cells, max_fpr, num_threads)
        };

        Self {
            index: depths.into_iter().zip(min_read_supports).collect()
        }
    }
}

impl ReadSupportIndex for RNAJunctionReadSupportIndex {
    type Input = ReadDepth;

    fn get_min_read_support(&self, input: ReadDepth) -> ReadSupport {
        *self
            .index
            .get(&input)
            .unwrap_or_else(|| panic!(
                "RNAJunctionReadSupportIndex has no entry for depth {input}."
            ))
    }
}


#[cfg(test)]
#[path = "../../tests/filtering/rna/junction_read_support_index.rs"]
mod tests;
