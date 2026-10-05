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
use exacto_caller::prelude::*;
use exacto_core::prelude::{AnalyteType, FastaMap, GeneAnnotator, ReadDepth, ReadSupport};
use std::collections::HashMap;
use std::sync::Arc;

use crate::options::RNAVariantFilteringOptions;
use crate::variant_filtering::known_variants::KnownVariants;


pub(crate) struct VariantJudge<'a, A: GeneAnnotator + Sync> {
    known_variants: &'a KnownVariants<'a>,
    read_support_index: &'a RNAVariantReadSupportIndex,
    template_switch_filter: Option<TemplateSwitchFilter<'a, A>>,
    filtering: &'a RNAVariantFilteringOptions,
    max_ins_norm_edit_distance: f64,
    chromosome_names_map: &'a BiMap<Box<str>, u16>,
    fasta_map: &'a FastaMap
}

impl<'a, A: GeneAnnotator + Sync> VariantJudge<'a, A> {
    pub(crate) fn new(
        known_variants: &'a KnownVariants<'a>,
        read_support_index: &'a RNAVariantReadSupportIndex,
        transcript_models_map: &'a HashMap<usize, Arc<TranscriptModel>>,
        gene_annotator: &'a A,
        chromosome_names_map: &'a BiMap<Box<str>, u16>,
        fasta_map: &'a FastaMap,
        filtering: &'a RNAVariantFilteringOptions,
        max_ins_norm_edit_distance: f64,
        analyte_type: &AnalyteType
    ) -> Self {
        let template_switch_filter: Option<TemplateSwitchFilter<'a, A>> = (*analyte_type != AnalyteType::RNA)
            .then(|| TemplateSwitchFilter::new(
                transcript_models_map,
                gene_annotator,
                chromosome_names_map,
                fasta_map,
                filtering.template_switch_flank,
                filtering.template_switch_junction_min_homology,
                filtering.template_switch_junction_soft_min_homology,
                filtering.template_switch_junction_max_breakpoint_dispersion,
                filtering.template_switch_foldback_min_stem,
                filtering.template_switch_foldback_max_loop_len,
                filtering.template_switch_foldback_max_distance,
                filtering.template_switch_foldback_slack
            ));
        Self { known_variants, read_support_index, template_switch_filter, filtering, max_ins_norm_edit_distance, chromosome_names_map, fasta_map }
    }

    pub(crate) fn judge(&self, variant_call: &VariantCall, (depth_1, depth_2): (u32, u32)) -> Result<Option<u32>, VariantCallFailure> {
        let op: &GraphOperation = variant_call.get_consensus_graph_operation();

        // Gate 0. A call matching no allowed variant is judged no further, so it never becomes a
        // site to phase on. A listed call faces every gate below.
        if !self.known_variants.allows(op, self.filtering.allowed_variant_max_distance, self.max_ins_norm_edit_distance) {
            return Err(VariantCallFailure::NotListed);
        }

        // Gate 1.
        if self.known_variants.confirms(op) {
            return Ok(None);
        }

        // Gates 2 and 3. A sequence-level call is one site, read at one depth; a junction-type call
        // has to clear both breakpoints. The repeat check reads the inserted sequence, so it has to
        // see the synthesized consensus rather than one member's spelling of it.
        let (is_repeat, repeat_length): (bool, u32) = is_repeat_variant(
            op,
            self.chromosome_names_map,
            self.fasta_map,
            self.filtering.min_homopolymer_len,
            self.filtering.min_dinucleotide_context_len
        );
        let repeat_length: u32 = if is_repeat { repeat_length.min(self.filtering.max_slippage_repeat_len) } else { 0 };

        // Gate 2. The batch sets a call's total depth to the deeper of its two sides.
        let total_depth: i32 = variant_call.get_total_depth();
        let min_total_depth: ReadDepth = self.filtering.min_total_depth as ReadDepth;
        if total_depth < min_total_depth as i32 {
            return Err(VariantCallFailure::LowTotalDepth { total_depth: total_depth as ReadDepth, min_total_depth });
        }

        // Gate 3. The reads are counted, not the records: a read with two records in the call
        // supports it once.
        let num_reads: ReadSupport = variant_call.get_num_reads();
        let min_reads: ReadSupport = self.filtering.min_reads as ReadSupport;
        if num_reads < min_reads {
            return Err(VariantCallFailure::TooFewReads { num_reads, min_reads, repeat_len: None });
        }
        let min_reads: ReadSupport = min_reads_at_site(self.read_support_index, repeat_length, (depth_1, depth_2));
        if num_reads < min_reads {
            return Err(VariantCallFailure::TooFewReads { num_reads, min_reads, repeat_len: Some(repeat_length) });
        }

        // Gate 4. Last, so a call that also fails a read support gate is reported as that.
        if let Some(reason) = self.template_switch_filter.as_ref().and_then(|filter| filter.identify_template_switch(variant_call)) {
            return Err(VariantCallFailure::TemplateSwitch(reason));
        }
        Ok(Some(repeat_length))
    }
}


pub(crate) fn min_reads_at_site(
    read_support_index: &RNAVariantReadSupportIndex,
    repeat_length: u32,
    (depth_1, depth_2): (ReadDepth, ReadDepth)
) -> ReadSupport {
    let floor = |depth: ReadDepth| -> ReadSupport {
        if depth == 0 {
            ReadSupport::MAX
        } else {
            read_support_index.get_min_read_support((repeat_length, depth))
        }
    };
    floor(depth_1).max(floor(depth_2))
}


#[cfg(test)]
#[path = "../tests/variant_filtering/variant_judge.rs"]
mod tests;