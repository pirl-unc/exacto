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


use exacto_util::prelude::*;
use std::cmp::PartialEq;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use serde::{Deserialize, Serialize};

use crate::prelude::{Strand, VariantCall};
use crate::structs::reference_transcript_match::ReferenceTranscriptMatch;
use crate::structs::transcript_model_exon::TranscriptModelExon;
use crate::structs::transcript_model_splice_junction::TranscriptModelSpliceJunction;


#[derive(Debug,Serialize,Deserialize)]
pub struct TranscriptModel {
    pub transcript_id: usize,
    pub reference_transcript_matches: Vec<ReferenceTranscriptMatch>,
    pub read_ids: Vec<usize>,
    pub exons: Vec<TranscriptModelExon>,
    pub splice_junctions: Vec<TranscriptModelSpliceJunction>,
    pub sequence_variant_calls: Vec<VariantCall>,
    pub splice_variant_calls: HashMap<Box<str>,Vec<VariantCall>>
}

impl PartialEq for TranscriptModel {
    fn eq(&self, other: &Self) -> bool {
        self.transcript_id == other.transcript_id &&
            self.reference_transcript_matches == other.reference_transcript_matches &&
            self.read_ids == other.read_ids &&
            self.exons == other.exons &&
            self.splice_junctions == other.splice_junctions &&
            self.sequence_variant_calls == other.sequence_variant_calls &&
            self.splice_variant_calls == other.splice_variant_calls
    }
}

impl Eq for TranscriptModel {}

impl Hash for TranscriptModel {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.transcript_id.hash(state);
        self.reference_transcript_matches.hash(state);
        self.read_ids.hash(state);
        self.exons.hash(state);
        self.splice_junctions.hash(state);
        self.sequence_variant_calls.hash(state);
        let mut sorted_keys: Vec<_> = self.splice_variant_calls.keys().collect();
        sorted_keys.sort();
        for key in sorted_keys {
            key.hash(state);
            self.splice_variant_calls[key].hash(state);
        }
    }
}

impl TranscriptModel {
    pub fn new(
        transcript_id: usize,
        reference_transcript_matches: Vec<ReferenceTranscriptMatch>,
        read_ids: Vec<usize>
    ) -> Self {
        Self {
            transcript_id: transcript_id,
            reference_transcript_matches: reference_transcript_matches,
            read_ids: read_ids,
            exons: Vec::new(),
            splice_junctions: Vec::new(),
            sequence_variant_calls: Vec::new(),
            splice_variant_calls: HashMap::new()
        }
    }

    pub fn add_exon(&mut self, exon: TranscriptModelExon) {
        self.exons.push(exon);
    }

    pub fn add_splice_junction(&mut self, splice_junction: TranscriptModelSpliceJunction) {
        self.splice_junctions.push(splice_junction);
    }

    pub fn add_sequence_variant_call(&mut self, variant_call: VariantCall) {
        self.sequence_variant_calls.push(variant_call);
    }
    
    pub fn add_splice_variant_call(&mut self, reference_transcript_id: &str, variant_call: VariantCall) {
        self.splice_variant_calls
            .entry(reference_transcript_id.into())
            .or_insert_with(Vec::new)
            .push(variant_call);
    }

    pub fn get_start_position(&self) -> (u16,u32) {
        let first_exon_strand = self.exons.first().unwrap().strand.clone();
        if first_exon_strand == Strand::Forward {
            (self.exons.first().unwrap().chromosome_id, self.exons.first().unwrap().start)
        } else {
            (self.exons.first().unwrap().chromosome_id, self.exons.first().unwrap().end)
        }
    }

    pub fn get_end_position(&self) -> (u16,u32) {
        let last_exon_strand = self.exons.last().unwrap().strand.clone();
        if last_exon_strand == Strand::Forward {
            (self.exons.last().unwrap().chromosome_id, self.exons.first().unwrap().end)
        } else {
            (self.exons.last().unwrap().chromosome_id, self.exons.first().unwrap().start)
        }
    }

    pub fn vectorize_exons(
        &self,
        chromosome_id: u16,
        start: u32,
        end: u32,
        aligned_value: i8,
        unaligned_value: i8
    ) -> Vec<i8> {
        let v_size: usize = (end - start + 1) as usize;
        let mut v: Vec<i8> = vec![unaligned_value; v_size];
        for exon in self.exons.iter() {
            if chromosome_id == exon.chromosome_id {
                match find_overlap((exon.start as isize, exon.end as isize), (start as isize, end as isize)) {
                    Some((x,y)) => {
                        for pos in x..=y {
                            let i = (pos as usize) - (start as usize);
                            v[i] = aligned_value;
                        }
                    }
                    None => {}
                }
            }
        }
        v
    }
}

impl Clone for TranscriptModel {
    fn clone(&self) -> Self {
        TranscriptModel {
            transcript_id: self.transcript_id,
            reference_transcript_matches: self.reference_transcript_matches.clone(),
            read_ids: self.read_ids.clone(),
            exons: self.exons.clone(),
            splice_junctions: self.splice_junctions.clone(),
            sequence_variant_calls: self.sequence_variant_calls.clone(),
            splice_variant_calls: self.splice_variant_calls.clone()
        }
    }
}
