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


use bimap::BiMap;
use exacto_core::prelude::*;

use crate::prelude::*;


pub struct StrandBiasFilter<'a> {
    read_depths: &'a BAMReadDepths,
    chromosome_names_map: &'a BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    alpha: f64
}

impl <'a> StrandBiasFilter<'a> {
    pub fn new(
        read_depths: &'a BAMReadDepths,
        chromosome_names_map: &'a BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
        alpha: f64
    ) -> Self {
        Self {
            read_depths: read_depths,
            chromosome_names_map: chromosome_names_map,
            alpha: alpha
        }
    }
}

impl VariantFilter for StrandBiasFilter<'_> {
    type Input = VariantCall;

    fn passes(&self, variant_call: &VariantCall) -> bool {
        let operation: &GraphOperation = variant_call.get_consensus_graph_operation();
        if matches!(operation.get_variant_type(), VariantType::Breakpoint | VariantType::Translocation) {
            return true;
        }
        let chromosome_1: &str = self.chromosome_names_map.get_by_right(&operation.get_chromosome_1()).unwrap();
        let chromosome_2: &str = self.chromosome_names_map.get_by_right(&operation.get_chromosome_2()).unwrap();
        let variant_records = variant_call.get_variant_records();
        let alt_fwd_1: u64 = variant_records.iter().filter(|record| *record.get_strand_1() == Strand::Forward).count() as u64;
        let alt_fwd_2: u64 = variant_records.iter().filter(|record| *record.get_strand_2() == Strand::Forward).count() as u64;
        let alt_rev_1: u64 = variant_records.len() as u64 - alt_fwd_1;
        let alt_rev_2: u64 = variant_records.len() as u64 - alt_fwd_2;
        let (ref_fwd_1, ref_rev_1): (ReadSupport, ReadSupport) = self.read_depths.get_strands(chromosome_1, operation.get_position_1());
        let (ref_fwd_2, ref_rev_2): (ReadSupport, ReadSupport) = self.read_depths.get_strands(chromosome_2, operation.get_position_2());

        // Test strand bias for positions 1 and 2.
        !has_strand_bias(alt_fwd_1, alt_rev_1, ref_fwd_1 as u64, ref_rev_1 as u64, self.alpha)
            && !has_strand_bias(alt_fwd_2, alt_rev_2, ref_fwd_2 as u64, ref_rev_2 as u64, self.alpha)
    }
}


fn has_strand_bias(
    alt_fwd: u64,
    alt_rev: u64,
    ref_fwd: u64,
    ref_rev: u64,
    alpha: f64
) -> bool {
    if alt_fwd + alt_rev == 0 || ref_fwd + ref_rev == 0 {
        return false;
    }
    fisher_exact_test_two_sided(
        alt_fwd,
        alt_rev,
        ref_fwd,
        ref_rev
    ) < alpha
}


#[cfg(test)]
#[path = "../../tests/filtering/dna/strand_bias_filter.rs"]
mod tests;
