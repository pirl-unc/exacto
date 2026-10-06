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


use bimap::{BiMap};
use exacto_caller::prelude::{GraphOperationView, VariantType};
use exacto_core::prelude::*;
use serde::{Serialize, Deserialize};
use std::collections::{HashMap, HashSet};
use std::str::FromStr;

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct AssembledTranscript {
    pub id: Box<str>,
    pub sequence: Box<str>,
    pub read_ids: Vec<Box<str>>,
    pub assembled_transcript_name: Box<str>,
    pub transcript_alignment: AssembledTranscriptAlignment,
    pub rna_variants: BiMap<u32, GraphOperationView>,
    pub dna_variants: BiMap<u32, GraphOperationView>,
    pub integrated_variant_ids: HashMap<u32, HashSet<u32>>, // HashMap<RNA variant ID, HashSet<DNA variant ID>>
    pub dna_variant_read_names: HashMap<u32, Box<str>>,     // HashMap<DNA variant ID, Box<read names>>
    pub rna_variant_read_spans: Vec<(u32, u32, u32)>,       // (read_start, read_end, canonical RNA variant ID)
    pub rna_variant_types: HashMap<u32, VariantType>,       // HashMap<canonical RNA variant ID, its variant_type>
    pub reference_stitched_intervals: Vec<(u32, u32)>,
    pub proteoforms: Vec<Proteoform>
}

impl AssembledTranscript {
    pub fn new(
        id: Box<str>,
        sequence: Box<str>,
        read_ids: Vec<Box<str>>,
        assembled_transcript_name: Box<str>,
        transcript_alignment: AssembledTranscriptAlignment,
        rna_variants: BiMap<u32, GraphOperationView>,
        dna_variants: BiMap<u32, GraphOperationView>,
        integrated_variant_ids: HashMap<u32, HashSet<u32>>,
        dna_variant_read_names: HashMap<u32, Box<str>>,
        rna_variant_read_spans: Vec<(u32, u32, u32)>,
        rna_variant_types: HashMap<u32, VariantType>,
        reference_stitched_intervals: Vec<(u32, u32)>
    ) -> Self {
        Self {
            id,
            sequence,
            read_ids,
            assembled_transcript_name,
            transcript_alignment,
            rna_variants,
            dna_variants,
            integrated_variant_ids,
            dna_variant_read_names,
            rna_variant_read_spans,
            rna_variant_types,
            reference_stitched_intervals,
            proteoforms: Vec::new()
        }
    }

    pub fn new_with_default(
        id: Box<str>,
        sequence: Box<str>
    ) -> Self {
        let assembled_transcript_name: Box<str> = id.clone();
        Self::new(
            id,
            sequence,
            Vec::new(),
            assembled_transcript_name,
            AssembledTranscriptAlignment::default(),
            BiMap::new(),
            BiMap::new(),
            HashMap::new(),
            HashMap::new(),
            Vec::new(),
            HashMap::new(),
            Vec::new()
        )
    }

    pub fn get_dna_variant(&self, id: u32) -> &GraphOperationView {
        self.dna_variants.get_by_left(&id)
            .expect("DNA variant id not registered on Transcript")
    }

    pub fn get_dna_variant_id(&self, gov: &GraphOperationView) -> u32 {
        *self.dna_variants.get_by_right(gov)
            .expect("DNA variant GraphOperationView not registered on Transcript")
    }

    pub fn get_id(&self) -> &str {
        &*self.id
    }

    pub fn get_assembled_transcript_name(&self) -> &str {
        &*self.assembled_transcript_name
    }

    pub fn get_integrated_dna_variant_ids(&self, rna_variant_id: u32) -> &HashSet<u32> {
        self.integrated_variant_ids.get(&rna_variant_id)
            .expect("RNA variant id has no integrated DNA variant entry")
    }

    pub fn get_nucleotide(&self, position: usize) -> AssembledTranscriptNucleotide {
        // Step 1. Base nucleotide from the assembled sequence.
        let sequence: &str = &self.sequence[position..position + 1];
        let nucleotide: Nucleotide = Nucleotide::from_str(sequence).unwrap_or(Nucleotide::N);

        // Step 2. Find the Base item covering this position.
        let ts_item = self.transcript_alignment.items.iter().find(|item| {
            matches!(item.item_type, TranscriptAlignmentItemType::Base { .. })
                && (item.read_start as usize) <= position
                && position <= (item.read_end as usize)
        });
        let transcript_alignment_index: Option<u32> = ts_item.map(|item| item.index);

        // Step 3. RNA variant carried by THIS base.
        let (span_start, span_end): (u32, u32) = ts_item
            .map_or((position as u32, position as u32), |item| (item.read_start, item.read_end));
        let rna_variant_id: Option<u32> = self.try_get_rna_variant_id_for_read_span(span_start, span_end);

        // Step 4. RNA variant from an event immediately preceding this base
        let preceding_event_rna_variant_id = self.transcript_alignment.items.iter()
            .find(|item| {
                matches!(item.item_type, TranscriptAlignmentItemType::Event { .. })
                    && (item.read_end as usize) == position
            })
            .and_then(|item| self.try_get_rna_variant_id(&item.graph_operation_view));

        // Step 5. DNA variants integrated against the resolved RNA variant ids.
        // Empty integrated_variant_ids yields None from .get() - no panic on absent data.
        let dna_variant_ids: Option<HashSet<u32>> = rna_variant_id
            .and_then(|id| self.integrated_variant_ids.get(&id))
            .cloned();
        let preceding_event_dna_variant_ids: Option<HashSet<u32>> =
            preceding_event_rna_variant_id
                .and_then(|id| self.integrated_variant_ids.get(&id))
                .cloned();

        // Step 6. Provenance: was this base pasted in from the reference transcript?
        let is_reference_stitched: bool = self.reference_stitched_intervals.iter()
            .any(|&(start, end)| (position as u32) >= start && (position as u32) < end);

        AssembledTranscriptNucleotide::new(
            nucleotide,
            position as u32,
            transcript_alignment_index,
            rna_variant_id,
            dna_variant_ids,
            preceding_event_rna_variant_id,
            preceding_event_dna_variant_ids,
            is_reference_stitched
        )
    }

    pub fn get_read_ids(&self) -> &Vec<Box<str>> {
        &self.read_ids
    }

    pub fn get_rna_variant(&self, id: u32) -> &GraphOperationView {
        self.rna_variants.get_by_left(&id)
            .expect("RNA variant id not registered on Transcript")
    }

    pub fn get_rna_variant_id(&self, gov: &GraphOperationView) -> u32 {
        *self.rna_variants.get_by_right(gov)
            .expect("RNA variant GraphOperationView not registered on Transcript")
    }

    pub fn try_get_rna_variant_id(&self, gov: &GraphOperationView) -> Option<u32> {
        self.rna_variants.get_by_right(gov).copied()
    }
    
    pub fn try_get_rna_variant_id_for_read_span(&self, read_start: u32, read_end: u32) -> Option<u32> {
        self.rna_variant_read_spans.iter()
            .filter(|(span_start, span_end, _)| *span_start <= read_start && read_end <= *span_end)
            .min_by_key(|(span_start, span_end, _)| span_end - span_start)
            .map(|(_, _, id)| *id)
    }

    pub fn try_get_dna_variant_id(&self, gov: &GraphOperationView) -> Option<u32> {
        self.dna_variants.get_by_right(gov).copied()
    }

    pub fn get_sequence(&self) -> &str {
        &self.sequence
    }

    pub fn get_transcript_alignment(&self) -> &AssembledTranscriptAlignment {
        &self.transcript_alignment
    }
    
    pub fn keeps_reading_frame(&self, id: u32) -> bool {
        match self.rna_variant_types.get(&id) {
            Some(VariantType::SingleNucleotideVariant | VariantType::MultiNucleotideVariant) => true,
            Some(VariantType::Insertion | VariantType::Deletion) => !self.get_rna_variant(id).is_frameshift(),
            _ => false
        }
    }

    pub fn get_amino_acids<'a>(&'a self, proteoform: &'a Proteoform) -> impl Iterator<Item = AminoAcid> + 'a {
        proteoform.sequence.chars().enumerate().map(move |(index, amino_acid)| {
            let codon_start: usize = proteoform.orf_start as usize + 3 * index;
            let nucleotides: [AssembledTranscriptNucleotide; 3] =
                std::array::from_fn(|offset| self.get_nucleotide(codon_start + offset));
            AminoAcid::new(index as u32, amino_acid, nucleotides)
        })
    }

    pub fn translate(
        &mut self,
        translation_strategy: &TranslationStrategy,
        start_codons: &HashSet<&str>
    ) {
        // The complete ORFs the strategy keeps.
        self.proteoforms = identify_open_reading_frames(&*self.sequence, translation_strategy, start_codons)
            .into_iter()
            .enumerate()
            .map(|(id, (sequence, orf_start, orf_end, _))| Proteoform::new(id as u32, orf_start, orf_end, sequence))
            .collect();
    }
}

impl Clone for AssembledTranscript {
    fn clone(&self) -> Self {
        AssembledTranscript {
            id: self.id.clone(),
            sequence: self.sequence.clone(),
            read_ids: self.read_ids.clone(),
            assembled_transcript_name: self.assembled_transcript_name.clone(),
            transcript_alignment: self.transcript_alignment.clone(),
            rna_variants: self.rna_variants.clone(),
            dna_variants: self.dna_variants.clone(),
            integrated_variant_ids: self.integrated_variant_ids.clone(),
            dna_variant_read_names: self.dna_variant_read_names.clone(),
            rna_variant_read_spans: self.rna_variant_read_spans.clone(),
            rna_variant_types: self.rna_variant_types.clone(),
            reference_stitched_intervals: self.reference_stitched_intervals.clone(),
            proteoforms: self.proteoforms.clone()
        }
    }
}
