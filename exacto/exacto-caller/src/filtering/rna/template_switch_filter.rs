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
use std::collections::HashMap;
use std::sync::Arc;

use crate::prelude::*;


pub struct TemplateSwitchFilter<'a, A: GeneAnnotator + Sync> {
    transcript_models_map: &'a HashMap<ReadID, Arc<TranscriptModel>>,
    gene_annotator: &'a A,
    chromosome_names_map: &'a BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    fasta_map: &'a FastaMap,
    flank: u32,
    junction_min_homology: u32,
    junction_soft_min_homology: u32,
    junction_max_breakpoint_dispersion: u32,
    foldback_min_stem: u32,
    foldback_max_loop_len: u32,
    foldback_max_distance: u32,
    foldback_slack: u32
}

impl <'a, A: GeneAnnotator + Sync> TemplateSwitchFilter<'a, A> {
    pub fn new(
        transcript_models_map: &'a HashMap<ReadID, Arc<TranscriptModel>>,
        gene_annotator: &'a A,
        chromosome_names_map: &'a BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
        fasta_map: &'a FastaMap,
        flank: u32,
        junction_min_homology: u32,
        junction_soft_min_homology: u32,
        junction_max_breakpoint_dispersion: u32,
        foldback_min_stem: u32,
        foldback_max_loop_len: u32,
        foldback_max_distance: u32,
        foldback_slack: u32
    ) -> Self {
        Self {
            transcript_models_map: transcript_models_map,
            gene_annotator: gene_annotator,
            chromosome_names_map: chromosome_names_map,
            fasta_map: fasta_map,
            flank: flank,
            junction_min_homology: junction_min_homology,
            junction_soft_min_homology: junction_soft_min_homology,
            junction_max_breakpoint_dispersion: junction_max_breakpoint_dispersion,
            foldback_min_stem: foldback_min_stem,
            foldback_max_loop_len: foldback_max_loop_len,
            foldback_max_distance: foldback_max_distance,
            foldback_slack: foldback_slack
        }
    }

    pub fn identify_template_switch(&self, variant_call: &VariantCall) -> Option<TemplateSwitchReason> {
        let evidence: TemplateSwitchEvidence = characterize_template_switch(
            variant_call,
            self.transcript_models_map,
            self.gene_annotator,
            self.chromosome_names_map,
            self.fasta_map,
            self.flank,
            self.foldback_max_loop_len,
            self.foldback_slack
        );
        let verdict: TemplateSwitchVerdict = classify_template_switch(
            &evidence,
            self.junction_min_homology,
            self.junction_soft_min_homology,
            self.junction_max_breakpoint_dispersion,
            self.foldback_max_distance,
            self.foldback_min_stem
        );
        match verdict {
            TemplateSwitchVerdict::Flagged(reason) if is_template_switch(verdict, variant_call.get_consensus_graph_operation().get_variant_type()) => Some(reason),
            _ => None
        }
    }
}

impl<A: GeneAnnotator + Sync> VariantFilter for TemplateSwitchFilter<'_, A> {
    type Input = VariantCall;

    fn passes(&self, variant_call: &VariantCall) -> bool {
        self.identify_template_switch(variant_call).is_none()
    }
}


#[cfg(test)]
#[path = "../../tests/filtering/rna/template_switch_filter.rs"]
mod tests;