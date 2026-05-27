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


use std::collections::HashSet;
use serde::{Deserialize, Serialize};

use crate::prelude::*;


#[derive(Debug, Serialize, Deserialize)]
pub struct TranscriptModel {
    /// Globally-unique identifier across all TranscriptModels in a run
    id: usize,

    /// The originating assembled transcript's name (read name)
    assembled_transcript_name: Box<str>,

    /// The reference transcript IDs this contextualization is built against.
    /// Canonical form: sorted, deduplicated; equality on this Vec defines
    /// "same reference transcript combination."
    reference_transcript_ids: HashSet<Box<str>>,

    /// Reference transcript matches - multiple to accommodate fusion genes.
    reference_transcript_matches: Vec<ReferenceTranscriptMatch>,

    /// Contextualized alignment structure.
    alignment_structure: AlignmentStructure,

    variant_records: Vec<VariantRecord>
}

impl TranscriptModel {
    pub fn new(
        id: usize,
        assembled_transcript_name: Box<str>,
        reference_transcript_ids: HashSet<Box<str>>,
        reference_transcript_matches: Vec<ReferenceTranscriptMatch>,
        alignment_structure: AlignmentStructure,
        variant_records: Vec<VariantRecord>
    ) -> Self {
        Self {
            id,
            assembled_transcript_name,
            reference_transcript_ids,
            reference_transcript_matches,
            alignment_structure,
            variant_records,
        }
    }

    pub fn get_id(&self) -> usize {
        self.id
    }

    pub fn get_assembled_transcript_name(&self) -> &str {
        &*self.assembled_transcript_name
    }

    pub fn get_reference_transcript_ids(&self) -> &HashSet<Box<str>> {
        &self.reference_transcript_ids
    }

    pub fn get_reference_transcript_matches(&self) -> &Vec<ReferenceTranscriptMatch> {
        &self.reference_transcript_matches
    }

    pub fn get_alignment_structure(&self) -> &AlignmentStructure {
        &self.alignment_structure
    }

    pub fn get_variant_records(&self) -> &Vec<VariantRecord> {
        &self.variant_records
    }

    /// True if this model represents a reference transcript: matched to
    /// exactly one reference transcript with no variant records.
    /// Moved from `AssembledTranscript::is_reference_transcript` — that
    /// method walked per-combination data which now lives here.
    pub fn is_reference_transcript(&self) -> bool {
        self.reference_transcript_ids.len() == 1 && self.variant_records.is_empty()
    }

    /// Used during `TranscriptModelSet` ingestion to re-assign each
    /// per-combination model a globally-unique id (the constructor only
    /// gives it a local 0..N index scoped to its parent AT).
    pub fn set_id(&mut self, id: usize) {
        self.id = id;
    }
}

impl Clone for TranscriptModel {
    fn clone(&self) -> Self {
        TranscriptModel {
            id: self.id,
            assembled_transcript_name: self.assembled_transcript_name.clone(),
            reference_transcript_ids: self.reference_transcript_ids.clone(),
            reference_transcript_matches: self.reference_transcript_matches.clone(),
            alignment_structure: self.alignment_structure.clone(),
            variant_records: self.variant_records.clone(),
        }
    }
}
