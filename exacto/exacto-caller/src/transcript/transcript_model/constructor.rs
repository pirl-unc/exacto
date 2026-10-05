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


use super::TranscriptModel;

use bimap::BiMap;
use exacto_core::prelude::*;
use std::collections::HashSet;

use crate::prelude::*;


impl TranscriptModel {
    pub fn new(
        alignment_model: AlignmentModel,
        reference_transcript_matches: Vec<ReferenceTranscriptMatch>,
        gene_annotator: &(impl GeneAnnotator + Sync),
        chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
        fasta_map: &FastaMap
    ) -> Self {
        // Make sure reference transcript matches have unique reference gene IDs and transcript IDs.
        let mut reference_gene_ids: HashSet<ReferenceGeneID> = HashSet::new();
        let mut reference_transcript_ids: HashSet<ReferenceTranscriptID> = HashSet::new();
        for reference_transcript_match in reference_transcript_matches.iter() {
            let inserted: bool = reference_gene_ids.insert(
                reference_transcript_match.get_reference_gene_id().into()
            );
            if !inserted {
                panic!("Only 1 reference transcript match per reference gene ID is allowed.");
            }

            let inserted: bool = reference_transcript_ids.insert(
                reference_transcript_match.get_reference_transcript_id().into()
            );
            if !inserted {
                panic!("Found duplicate reference transcript IDs. This is not allowed.");
            }
        }

        let annotation: TranscriptModelAnnotation = identify_transcript_model_annotation(
            &alignment_model,
            &reference_transcript_matches,
            gene_annotator,
            chromosome_names_map,
            fasta_map
        );

        let exons: Vec<TranscriptModelExon> = identify_transcript_model_exons(&alignment_model);

        let splice_junctions: Vec<TranscriptModelSpliceJunction> = identify_transcript_model_splice_junctions(
            &alignment_model
        );
        
        Self {
            alignment_model: alignment_model,
            reference_transcript_matches: reference_transcript_matches,
            annotation: annotation,
            exons: exons,
            splice_junctions: splice_junctions,
            variant_records: Vec::new(),
            nmd_predictions: Vec::new()
        }
    }
}
