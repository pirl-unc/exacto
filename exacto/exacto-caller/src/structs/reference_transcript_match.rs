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


use serde::{Deserialize, Serialize};
use std::cmp::PartialEq;
use std::hash::{Hash, Hasher};

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct ReferenceTranscriptMatch {
    pub reference_gene_id: Box<str>,
    pub reference_transcript_id: Box<str>,
    pub num_overlap_bases: usize,
    pub num_transcript_only_bases: usize,
    pub num_reference_only_bases: usize,
    pub scoring_method: ReferenceTranscriptScoringMethod,
    pub score: f32
}

impl PartialEq for ReferenceTranscriptMatch {
    fn eq(&self, other: &Self) -> bool {
        self.reference_gene_id == other.reference_gene_id && 
            self.reference_transcript_id == other.reference_transcript_id &&
            self.num_overlap_bases == other.num_overlap_bases &&
            self.num_transcript_only_bases == other.num_transcript_only_bases &&
            self.num_reference_only_bases == other.num_reference_only_bases &&
            self.scoring_method == other.scoring_method &&
            self.score == other.score
    }
}

impl Eq for ReferenceTranscriptMatch {}

impl Hash for ReferenceTranscriptMatch {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.reference_gene_id.hash(state);
        self.reference_transcript_id.hash(state);
        self.num_overlap_bases.hash(state);
        self.num_transcript_only_bases.hash(state);
        self.num_reference_only_bases.hash(state);
        self.scoring_method.hash(state);
        self.score.to_bits().hash(state);
    }
}

impl ReferenceTranscriptMatch {
    pub fn new(
        reference_gene_id: &str,
        reference_transcript_id: &str,
        num_overlap_bases: usize,
        num_transcript_only_bases: usize,
        num_reference_only_bases: usize,
        scoring_method: ReferenceTranscriptScoringMethod,
        score: f32
    ) -> Self {
        Self {
            reference_gene_id: reference_gene_id.into(),
            reference_transcript_id: reference_transcript_id.into(),
            num_overlap_bases: num_overlap_bases,
            num_transcript_only_bases: num_transcript_only_bases,
            num_reference_only_bases: num_reference_only_bases,
            scoring_method: scoring_method,
            score: score
        }
    }
}

impl Clone for ReferenceTranscriptMatch {
    fn clone(&self) -> Self {
        ReferenceTranscriptMatch {
            reference_gene_id: self.reference_gene_id.clone(),
            reference_transcript_id: self.reference_transcript_id.clone(),
            num_overlap_bases: self.num_overlap_bases,
            num_transcript_only_bases: self.num_transcript_only_bases,
            num_reference_only_bases: self.num_reference_only_bases,
            scoring_method: self.scoring_method.clone(),
            score: self.score
        }
    }
}
