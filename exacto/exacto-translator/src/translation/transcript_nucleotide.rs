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
use serde::{Deserialize,Serialize};
use std::collections::HashSet;

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct TranscriptNucleotide {
    nucleotide: Nucleotide,
    transcript_read_position: u32,
    transcript_structure_index: Option<u32>,
    rna_variant_id: Option<u32>,
    dna_variant_ids: Option<HashSet<u32>>,
    preceding_event_rna_variant_id: Option<u32>,
    preceding_event_dna_variant_ids: Option<HashSet<u32>>

}

impl TranscriptNucleotide {
    pub fn new(
        nucleotide: Nucleotide,
        transcript_read_position: u32,
        transcript_structure_index: Option<u32>,
        rna_variant_id: Option<u32>,
        dna_variant_ids: Option<HashSet<u32>>,
        preceding_event_rna_variant_id: Option<u32>,
        preceding_event_dna_variant_ids: Option<HashSet<u32>>
    ) -> Self {
        Self {
            nucleotide: nucleotide,
            transcript_read_position: transcript_read_position,
            transcript_structure_index: transcript_structure_index,
            rna_variant_id: rna_variant_id,
            dna_variant_ids: dna_variant_ids,
            preceding_event_rna_variant_id: preceding_event_rna_variant_id,
            preceding_event_dna_variant_ids: preceding_event_dna_variant_ids
        }
    }

    pub fn get_dna_variant_ids(&self) -> &Option<HashSet<u32>> {
        &self.dna_variant_ids
    }

    pub fn get_nucleotide(&self) -> &Nucleotide {
        &self.nucleotide
    }

    pub fn get_preceding_event_dna_variant_ids(&self) -> &Option<HashSet<u32>> {
        &self.preceding_event_dna_variant_ids
    }

    pub fn get_preceding_event_rna_variant_id(&self) -> Option<u32> {
        self.preceding_event_rna_variant_id
    }

    pub fn get_rna_variant_id(&self) -> Option<u32> {
        self.rna_variant_id
    }

    pub fn get_transcript_read_position(&self) -> u32 {
        self.transcript_read_position
    }

    pub fn get_transcript_structure_index(&self) -> Option<u32> {
        self.transcript_structure_index
    }

    pub fn is_variant(&self) -> bool {
        if self.rna_variant_id.is_some() {
            return true;
        }
        if self.preceding_event_rna_variant_id.is_some() {
            return true;
        }
        false
    }
}

impl Clone for TranscriptNucleotide {
    fn clone(&self) -> Self {
        TranscriptNucleotide {
            nucleotide: self.nucleotide.clone(),
            transcript_read_position: self.transcript_read_position,
            transcript_structure_index: self.transcript_structure_index,
            rna_variant_id: self.rna_variant_id.clone(),
            dna_variant_ids: self.dna_variant_ids.clone(),
            preceding_event_rna_variant_id: self.preceding_event_rna_variant_id.clone(),
            preceding_event_dna_variant_ids: self.preceding_event_dna_variant_ids.clone()
        }
    }
}
