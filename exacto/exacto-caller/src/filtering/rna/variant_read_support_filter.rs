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
use crate::filtering::repeat_context::is_repeat_variant;
use crate::filtering::rna::variant_read_support_index::RNAVariantReadSupportIndex;


pub struct RNAVariantReadSupportFilter<'a> {
    min_depth: ReadSupport,
    min_homopolymer_len: u32,
    min_dinucleotide_context_len: u32,
    max_slippage_repeat_len: u32,
    read_support_index: &'a RNAVariantReadSupportIndex,
    read_depths: &'a BAMReadDepths,
    chromosome_names_map: &'a BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    fasta_map: &'a FastaMap
}

impl <'a> RNAVariantReadSupportFilter<'a> {
    pub fn new(
        min_depth: ReadSupport,
        min_homopolymer_len: u32,
        min_dinucleotide_context_len: u32,
        max_slippage_repeat_len: u32,
        read_support_index: &'a RNAVariantReadSupportIndex,
        read_depths: &'a BAMReadDepths,
        chromosome_names_map: &'a BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
        fasta_map: &'a FastaMap
    ) -> Self {
        Self {
            min_depth: min_depth,
            min_homopolymer_len: min_homopolymer_len,
            min_dinucleotide_context_len: min_dinucleotide_context_len,
            max_slippage_repeat_len: max_slippage_repeat_len,
            read_support_index: read_support_index,
            read_depths: read_depths,
            chromosome_names_map: chromosome_names_map,
            fasta_map: fasta_map
        }
    }
}

impl VariantFilter for RNAVariantReadSupportFilter<'_> {
    type Input = VariantCall;

    fn passes(&self, variant_call: &VariantCall) -> bool {
        // Minimum total depth.
        let total_depth: i32 = variant_call.get_total_depth();
        if total_depth < self.min_depth as i32 {
            return false;
        }

        // Minimum number of reads.
        let operation: &GraphOperation = variant_call.get_consensus_graph_operation();
        let (is_repeat, repeat_length): (bool, u32) = is_repeat_variant(
            operation,
            self.chromosome_names_map,
            self.fasta_map,
            self.min_homopolymer_len,
            self.min_dinucleotide_context_len
        );
        let repeat_length: u32 = if is_repeat {
            repeat_length.min(self.max_slippage_repeat_len)
        } else {
            0
        };
        let chromosome_1: &str = self.chromosome_names_map.get_by_right(&operation.get_chromosome_1()).unwrap();
        let chromosome_2: &str = self.chromosome_names_map.get_by_right(&operation.get_chromosome_2()).unwrap();
        let depth_1: ReadDepth = self.read_depths.get_depth(chromosome_1, operation.get_position_1());
        let depth_2: ReadDepth = self.read_depths.get_depth(chromosome_2, operation.get_position_2());
        let min_read_support_1: ReadSupport = self.read_support_index.get_min_read_support((repeat_length, depth_1));
        let min_read_support_2: ReadSupport = self.read_support_index.get_min_read_support((repeat_length, depth_2));
        if variant_call.get_num_reads() < min_read_support_1
            || variant_call.get_num_reads() < min_read_support_2 {
            return false;
        }

        true
    }
}
