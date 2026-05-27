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
use exacto_caller::prelude::{AlignmentStructureBaseKind, GraphOperationView};
use exacto_core::prelude::translate as translate_sequence;
use exacto_core::prelude::*;
use serde::{Serialize, Deserialize};
use std::collections::{HashMap, HashSet};
use std::str::FromStr;

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct Transcript {
    /// Globally unique identifier for this Transcript.
    /// "{assembled_transcript_name}:{transcript_model_id}"
    pub id: Box<str>,
    pub sequence: Box<str>,
    pub read_ids: Vec<Box<str>>,
    pub assembled_transcript_name: Box<str>,
    pub transcript_model_id: u32,
    pub transcript_structure: TranscriptStructure,
    pub rna_variants: BiMap<u32, GraphOperationView>,
    pub dna_variants: BiMap<u32, GraphOperationView>,
    pub integrated_variant_ids: HashMap<u32, HashSet<u32>>, // HashMap<RNA variant ID, HashSet<DNA variant ID>>
    pub dna_variant_read_names: HashMap<u32, Box<str>>,     // HashMap<DNA variant ID, Box<read names>>
    pub primary_structures: Vec<PrimaryStructure>
}

impl Transcript {
    pub fn new(
        id: Box<str>,
        sequence: Box<str>,
        read_ids: Vec<Box<str>>,
        assembled_transcript_name: Box<str>,
        transcript_model_id: u32,
        transcript_structure: TranscriptStructure,
        rna_variants: BiMap<u32, GraphOperationView>,
        dna_variants: BiMap<u32, GraphOperationView>,
        integrated_variant_ids: HashMap<u32, HashSet<u32>>,
        dna_variant_read_names: HashMap<u32, Box<str>>
    ) -> Self {
        Self {
            id,
            sequence,
            read_ids,
            assembled_transcript_name,
            transcript_model_id,
            transcript_structure,
            rna_variants,
            dna_variants,
            integrated_variant_ids,
            dna_variant_read_names,
            primary_structures: Vec::new()
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
            0,
            TranscriptStructure::default(),
            BiMap::new(),
            BiMap::new(),
            HashMap::new(),
            HashMap::new()
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

    pub fn get_transcript_model_id(&self) -> u32 {
        self.transcript_model_id
    }

    pub fn get_integrated_dna_variant_ids(&self, rna_variant_id: u32) -> &HashSet<u32> {
        self.integrated_variant_ids.get(&rna_variant_id)
            .expect("RNA variant id has no integrated DNA variant entry")
    }

    pub fn get_nucleotide(&self, position: usize) -> TranscriptNucleotide {
        // Step 1. Base nucleotide from the assembled sequence
        let sequence: &str = &self.sequence[position..position + 1];
        let nucleotide: Nucleotide = Nucleotide::from_str(sequence).unwrap();

        // Step 2. Find the Base item covering this position.
        // read_start/read_end are inclusive on both ends (verified against
        // the structure-record writer: a base row's read_end == read_start +
        // sequence_length - 1, so a single-base mismatch row has
        // read_start == read_end).
        let ts_item = self.transcript_structure.items.iter().find(|item| {
            matches!(item.item_type, TranscriptStructureItemType::Base { .. })
                && (item.read_start as usize) <= position
                && position <= (item.read_end as usize)
        });
        let transcript_structure_index: Option<u32> = ts_item.map(|item| item.index);

        // Step 3. RNA variant carried by THIS base
        // Mismatch/Insertion bases are themselves the variant; Match bases are not.
        let rna_variant_id = ts_item.and_then(|item| match &item.item_type {
            TranscriptStructureItemType::Base { kind, .. } => match kind {
                AlignmentStructureBaseKind::Mismatch
                | AlignmentStructureBaseKind::Insertion => {
                    self.try_get_rna_variant_id(&item.graph_operation_view)
                }
                _ => None,
            },
            TranscriptStructureItemType::Event { .. } => None,
        });

        // Step 4. RNA variant from an event immediately preceding this base
        let preceding_event_rna_variant_id = self.transcript_structure.items.iter()
            .find(|item| {
                matches!(item.item_type, TranscriptStructureItemType::Event { .. })
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

        TranscriptNucleotide::new(
            nucleotide,
            position as u32,
            transcript_structure_index,
            rna_variant_id,
            dna_variant_ids,
            preceding_event_rna_variant_id,
            preceding_event_dna_variant_ids
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

    pub fn try_get_dna_variant_id(&self, gov: &GraphOperationView) -> Option<u32> {
        self.dna_variants.get_by_right(gov).copied()
    }

    pub fn get_sequence(&self) -> &str {
        &self.sequence
    }

    pub fn get_transcript_structure(&self) -> &TranscriptStructure {
        &self.transcript_structure
    }

    pub fn translate(
        &mut self,
        translation_strategy: &TranslationStrategy,
        start_codons: HashSet<&str>
    ) {
        // Step 1. Free self.primary_structures
        self.primary_structures = Vec::new();

        // Step 2. Find ORFs on the RNA sequence
        let peptides: Vec<(Box<str>, u32, u32, u32)> = translate_sequence(&*self.sequence, start_codons);
        if peptides.is_empty() {
            return;
        }

        // Step 3. Strategy picks which ORFs to keep
        let kept: Vec<(Box<str>, u32, u32, u32)> = match translation_strategy {
            TranslationStrategy::LongestORF => peptides
                .into_iter()
                .max_by_key(|(seq, _, _, _)| seq.len())
                .into_iter()
                .collect(),
            TranslationStrategy::AllORFs => peptides
        };

        // Step 4. Materialize a PrimaryStructure per ORF
        self.primary_structures = kept
            .into_iter()
            .enumerate()
            .map(|(id, (_aa_sequence, orf_start, orf_end, _))| {
                self.build_primary_structure(id as u32, orf_start, orf_end)
            })
            .collect();
    }
}

impl Transcript {
    fn build_primary_structure(
        &self,
        id: u32,
        orf_start: u32,
        orf_end: u32
    ) -> PrimaryStructure {
        let amino_acids: Vec<AminoAcid> = (orf_start..=orf_end)
            .step_by(3)
            .enumerate()
            .map(|(aa_index, codon_start)| {
                let nucleotides: Vec<TranscriptNucleotide> = (0..3u32)
                    .map(|i| self.get_nucleotide((codon_start + i) as usize))
                    .collect();
                AminoAcid::new(aa_index as u32, nucleotides)
            })
            .collect();

        PrimaryStructure::new(
            id,
            orf_start,
            orf_end,
            amino_acids
        )
    }
}

impl Clone for Transcript {
    fn clone(&self) -> Self {
        Transcript {
            id: self.id.clone(),
            sequence: self.sequence.clone(),
            read_ids: self.read_ids.clone(),
            assembled_transcript_name: self.assembled_transcript_name.clone(),
            transcript_model_id: self.transcript_model_id,
            transcript_structure: self.transcript_structure.clone(),
            rna_variants: self.rna_variants.clone(),
            dna_variants: self.dna_variants.clone(),
            integrated_variant_ids: self.integrated_variant_ids.clone(),
            dna_variant_read_names: self.dna_variant_read_names.clone(),
            primary_structures: self.primary_structures.clone()
        }
    }
}
