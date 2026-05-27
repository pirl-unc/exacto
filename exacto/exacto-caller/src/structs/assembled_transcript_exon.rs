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


use exacto_core::prelude::*;
use serde::{Deserialize, Serialize};
use std::cmp::PartialEq;
use std::hash::{Hash, Hasher};


pub fn vectorize_exons(
    exons: &Vec<AssembledTranscriptExon>,
    reference_chromosome_id: u16,
    reference_start: u32,
    reference_end: u32,
    aligned_value: i8,
    unaligned_value: i8
) -> Vec<i8> {
    let v_size: usize = (reference_end - reference_start + 1) as usize;
    let mut v: Vec<i8> = vec![unaligned_value; v_size];
    for exon in exons.iter() {
        if reference_chromosome_id == exon.reference_chromosome_id {
            match find_overlap((exon.reference_start as isize, exon.reference_end as isize), (reference_start as isize, reference_end as isize)) {
                Some((x,y)) => {
                    for pos in x..=y {
                        let i = (pos as usize) - (reference_start as usize);
                        v[i] = aligned_value;
                    }
                }
                None => {}
            }
        }
    }
    v
}

#[derive(Debug,Serialize,Deserialize)]
pub struct AssembledTranscriptExon {
    pub reference_chromosome_id: u16,
    pub reference_start: u32,
    pub reference_end: u32,
    pub reference_strand: Strand,
    pub exon_number: u16,
    pub read_start_position: u32,       // FASTX read sequence start position
    pub read_end_position: u32          // FASTX read sequence end position
}

impl PartialEq for AssembledTranscriptExon {
    fn eq(&self, other: &Self) -> bool {
        self.reference_chromosome_id == other.reference_chromosome_id &&
            self.reference_start == other.reference_start &&
            self.reference_end == other.reference_end &&
            self.reference_strand == other.reference_strand &&
            self.exon_number == other.exon_number &&
            self.read_start_position == other.read_start_position &&
            self.read_end_position == other.read_end_position
    }
}

impl Eq for AssembledTranscriptExon {}

impl Hash for AssembledTranscriptExon {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.reference_chromosome_id.hash(state);
        self.reference_start.hash(state);
        self.reference_end.hash(state);
        self.reference_strand.hash(state);
        self.exon_number.hash(state);
        self.read_start_position.hash(state);
        self.read_end_position.hash(state);
    }
}

impl AssembledTranscriptExon {
    pub fn new(
        reference_chromosome_id: u16,
        reference_start: u32,
        reference_end: u32,
        reference_strand: Strand,
        exon_number: u16,
        read_start_position: u32,
        read_end_position: u32
    ) -> Self {
        assert!(read_start_position <= read_end_position, "read_start_position should be less than or equal to read_end_position.");
        assert!(reference_start <= reference_end, "reference_start should be less than reference_end.");
        Self {
            reference_chromosome_id,
            reference_start,
            reference_end,
            reference_strand: reference_strand.clone(),
            exon_number,
            read_start_position,
            read_end_position
        }
    }
}

impl Clone for AssembledTranscriptExon {
    fn clone(&self) -> Self {
        AssembledTranscriptExon {
            reference_chromosome_id: self.reference_chromosome_id,
            reference_start: self.reference_start,
            reference_end: self.reference_end,
            reference_strand: self.reference_strand.clone(),
            exon_number: self.exon_number,
            read_start_position: self.read_start_position,
            read_end_position: self.read_end_position
        }
    }
}
