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


pub mod constructor;

use exacto_core::prelude::*;
use serde::{Deserialize, Serialize};
use std::hash::{Hash, Hasher};

use crate::prelude::*;


#[derive(Debug, Serialize, Deserialize)]
pub struct TranscriptModel {
    alignment_model: AlignmentModel,
    reference_transcript_matches: Vec<ReferenceTranscriptMatch>,
    annotation: TranscriptModelAnnotation,
    exons: Vec<TranscriptModelExon>,
    splice_junctions: Vec<TranscriptModelSpliceJunction>,

    /// Empty until `set_variant_records` stores the RNA variant caller's records.
    variant_records: Vec<VariantRecord>,

    /// Empty until `set_nmd_predictions` stores the NMD predictor's calls. Cleared when the
    /// splice junctions change, since the verdicts depend on the last junction.
    nmd_predictions: Vec<NonsenseMediatedDecayCall>
}

impl Eq for TranscriptModel {}

impl PartialEq for TranscriptModel {
    fn eq(&self, other: &Self) -> bool {
        if self.alignment_model.get_read_id() == other.alignment_model.get_read_id() {
            true
        } else {
            false
        }
    }
}

impl Hash for TranscriptModel {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.alignment_model.get_read_id().hash(state);
    }
}


impl TranscriptModel {
    pub fn get_alignment_model(&self) -> &AlignmentModel {
        &self.alignment_model
    }

    pub fn get_annotation(&self) -> &TranscriptModelAnnotation {
        &self.annotation
    }

    pub fn get_exons(&self) -> &Vec<TranscriptModelExon> {
        &self.exons
    }

    pub fn get_read_id(&self) -> ReadID {
        self.alignment_model.get_read_id()
    }

    /// The matched reference transcripts, each once, in the order the matching step set: best first.
    pub fn get_reference_transcript_ids(&self) -> Vec<ReferenceTranscriptID> {
        let mut reference_transcript_ids: Vec<ReferenceTranscriptID> = Vec::new();
        for reference in self.reference_transcript_matches.iter() {
            let reference_transcript_id: ReferenceTranscriptID = reference.get_reference_transcript_id().into();
            if !reference_transcript_ids.contains(&reference_transcript_id) {
                reference_transcript_ids.push(reference_transcript_id);
            }
        }
        reference_transcript_ids
    }

    pub fn get_reference_end_position(&self) -> (ReferenceChromosomeID, ReferencePosition) {
        let last_exon_strand: Strand = self.get_exons().last().unwrap().reference_strand.clone();
        if last_exon_strand == Strand::Forward {
            (self.get_exons().last().unwrap().reference_chromosome_id,
             self.get_exons().last().unwrap().reference_end)
        } else {
            (self.get_exons().first().unwrap().reference_chromosome_id,
             self.get_exons().first().unwrap().reference_end)
        }
    }

    pub fn get_reference_start_position(&self) -> (ReferenceChromosomeID, ReferencePosition) {
        let first_exon_strand = self.get_exons().first().unwrap().reference_strand.clone();
        if first_exon_strand == Strand::Forward {
            (self.get_exons().first().unwrap().reference_chromosome_id,
             self.get_exons().first().unwrap().reference_start)
        } else {
            (self.get_exons().last().unwrap().reference_chromosome_id,
             self.get_exons().last().unwrap().reference_start)
        }
    }

    pub fn get_reference_transcript_matches(&self) -> &Vec<ReferenceTranscriptMatch> {
        &self.reference_transcript_matches
    }

    pub fn get_splice_junctions(&self) -> &Vec<TranscriptModelSpliceJunction> {
        &self.splice_junctions
    }

    pub fn get_splice_junctions_key(&self) -> Vec<SpliceJunction> {
        self.splice_junctions
            .iter()
            .map(|junction| {
                SpliceJunction::new(
                    junction.reference_chromosome_id_1,
                    junction.reference_chromosome_id_2,
                    junction.reference_position_1,
                    junction.reference_position_2,
                    junction.reference_strand_1.clone(),
                    junction.reference_strand_2.clone()
                )
            })
            .collect()
    }

    pub fn get_variant_records(&self) -> &Vec<VariantRecord> {
        &self.variant_records
    }

    /// True if this model represents a reference transcript: matched to exactly
    /// one reference transcript with no variant records.
    pub fn is_reference_transcript(&self) -> bool {
        self.reference_transcript_matches.len() == 1 && self.variant_records.is_empty()
    }

    pub fn is_spliced(&self) -> bool {
        if self.get_splice_junctions().is_empty() {
            false
        } else {
            true
        }
    }
    
    pub fn retain_splice_junctions(&mut self, keep: impl FnMut(&TranscriptModelSpliceJunction) -> bool) {
        self.splice_junctions.retain(keep);
        for (i, splice_junction) in self.splice_junctions.iter_mut().enumerate() {
            splice_junction.splice_junction_number = (i + 1) as SpliceJunctionNumber;
        }
        self.nmd_predictions.clear();
    }

    pub fn set_variant_records(&mut self, variant_records: Vec<VariantRecord>) {
        self.variant_records = variant_records;
    }

    pub fn retain_variant_records(&mut self, keep: impl FnMut(&VariantRecord) -> bool) {
        self.variant_records.retain(keep);
    }

    pub fn replace_variant_record(&mut self, original: &VariantRecord, replacement: VariantRecord) {
        assert_eq!(replacement.get_read_id(), self.get_read_id());
        let record = self.variant_records.iter_mut().find(|record| *record == original)
            .expect("variant record must belong to this transcript model");
        *record = replacement;
    }

    pub fn get_nmd_predictions(&self) -> &Vec<NonsenseMediatedDecayCall> {
        &self.nmd_predictions
    }

    pub fn set_nmd_predictions(&mut self, nmd_predictions: Vec<NonsenseMediatedDecayCall>) {
        self.nmd_predictions = nmd_predictions;
    }
}

impl Clone for TranscriptModel {
    fn clone(&self) -> Self {
        TranscriptModel {
            alignment_model: self.alignment_model.clone(),
            reference_transcript_matches: self.reference_transcript_matches.clone(),
            annotation: self.annotation.clone(),
            exons: self.exons.clone(),
            splice_junctions: self.splice_junctions.clone(),
            variant_records: self.variant_records.clone(),
            nmd_predictions: self.nmd_predictions.clone()
        }
    }
}


#[cfg(test)]
#[path = "../../tests/transcript/transcript_model.rs"]
mod tests;