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


use std::cmp::PartialEq;
use std::hash::{Hash, Hasher};
use bimap::BiMap;
use serde::{Deserialize, Serialize};

use crate::prelude::VariantCall;
use crate::structs::transcript_model_exon::TranscriptModelExon;
use crate::structs::transcript_model_splice_junction::TranscriptModelSpliceJunction;


#[derive(Debug,Serialize,Deserialize)]
pub struct TranscriptModel {
    pub transcript_id: usize,
    pub reference_transcript_ids: Vec<Box<str>>,
    pub read_ids: Vec<usize>,
    pub exons: Vec<TranscriptModelExon>,
    pub splice_junctions: Vec<TranscriptModelSpliceJunction>,
    pub variant_calls: Vec<VariantCall>
}

impl PartialEq for TranscriptModel {
    fn eq(&self, other: &Self) -> bool {
        self.transcript_id == other.transcript_id &&
            self.reference_transcript_ids == other.reference_transcript_ids &&
            self.read_ids == other.read_ids &&
            self.exons == other.exons &&
            self.splice_junctions == other.splice_junctions &&
            self.variant_calls == other.variant_calls
    }
}

impl Eq for TranscriptModel {}

impl Hash for TranscriptModel {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.transcript_id.hash(state);
        self.reference_transcript_ids.hash(state);
        self.read_ids.hash(state);
        self.exons.hash(state);
        self.splice_junctions.hash(state);
        self.variant_calls.hash(state);
    }
}

impl TranscriptModel {
    pub fn new(
        transcript_id: usize,
        reference_transcript_ids: Vec<Box<str>>,
        read_ids: Vec<usize>
    ) -> Self {
        Self {
            transcript_id: transcript_id,
            reference_transcript_ids: reference_transcript_ids,
            read_ids: read_ids,
            exons: Vec::new(),
            splice_junctions: Vec::new(),
            variant_calls: Vec::new()
        }
    }

    pub fn add_exon(&mut self, exon: TranscriptModelExon) {
        self.exons.push(exon);
    }

    pub fn add_splice_junction(&mut self, splice_junction: TranscriptModelSpliceJunction) {
        self.splice_junctions.push(splice_junction);
    }

    pub fn add_variant_call(&mut self, variant_call: VariantCall) {
        self.variant_calls.push(variant_call);
    }
}

impl Clone for TranscriptModel {
    fn clone(&self) -> Self {
        TranscriptModel {
            transcript_id: self.transcript_id,
            reference_transcript_ids: self.reference_transcript_ids.clone(),
            read_ids: self.read_ids.clone(),
            exons: self.exons.clone(),
            splice_junctions: self.splice_junctions.clone(),
            variant_calls: self.variant_calls.clone()
        }
    }
}
