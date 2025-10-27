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


use bimap::BiMap;
use exacto_core::prelude::*;
use itertools::Itertools;
use std::cmp::PartialEq;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use serde::{Deserialize, Serialize};

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct TranscriptModel {
    id: usize,

    /// Initial alignment structure.
    alignment_structure: AlignmentStructure,

    /// Exons of the initial alignment structure.
    exons: Vec<TranscriptModelExon>,

    /// Introns of the initial alignment structure.
    introns: Vec<TranscriptModelIntron>,

    /// A map between reference transcript IDs and contextualized AlignmentStructure object.
    contextualized_alignment_structures_map: HashMap<Vec<Box<str>>, AlignmentStructure>,

    /// A map between reference transcript ID and ReferenceTranscriptMatch object.
    reference_transcript_matches_map: HashMap<Box<str>, ReferenceTranscriptMatch>,

    /// A map between reference transcript IDs and VariantRecord objects.
    variant_records_map: HashMap<Vec<Box<str>>, Vec<VariantRecord>>,

    /// A map between chromosome names and IDs.
    chromosome_names_map: BiMap<Box<str>, u16>
}

impl PartialEq for TranscriptModel {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id &&
            self.exons == other.exons &&
            self.introns == other.introns &&
            self.reference_transcript_matches_map == other.reference_transcript_matches_map
    }
}

impl Eq for TranscriptModel {}

impl Hash for TranscriptModel {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id.hash(state);
        let mut sorted_keys: Vec<_> = self.reference_transcript_matches_map.keys().collect();
        sorted_keys.sort();
        for key in sorted_keys {
            key.hash(state);
            self.reference_transcript_matches_map[key].hash(state);
        }
    }
}

/// API methods
impl TranscriptModel {
    pub fn new(
        id: usize,
        alignment_structure: &AlignmentStructure,
        chromosome_names_map: &BiMap<Box<str>,u16>,
        reference_genome_fasta_file: &str
    ) -> Self {
        Self {
            id,
            alignment_structure: alignment_structure.clone(),
            exons: alignment_structure.identify_exons(),
            introns: alignment_structure.identify_introns(&chromosome_names_map, reference_genome_fasta_file),
            contextualized_alignment_structures_map: HashMap::new(),
            reference_transcript_matches_map: HashMap::new(),
            variant_records_map: HashMap::new(),
            chromosome_names_map: chromosome_names_map.clone()
        }
    }

    pub fn get_contextualized_structures_map(&self) -> &HashMap<Vec<Box<str>>, AlignmentStructure> {
        &self.contextualized_alignment_structures_map
    }

    pub fn get_chromosome_names_map(&self) -> &BiMap<Box<str>, u16> {
        &self.chromosome_names_map
    }

    pub fn get_exons(&self) -> &Vec<TranscriptModelExon> {
        &self.exons
    }

    pub fn get_id(&self) -> usize {
        self.id
    }

    pub fn get_introns(&self) -> &Vec<TranscriptModelIntron> {
        &self.introns
    }

    pub fn get_read_id(&self) -> usize {
        self.alignment_structure.get_read_id()
    }

    pub fn get_reference_end_position(&self) -> (u16, usize) {
        let last_exon_strand = self.exons.last().unwrap().reference_strand.clone();
        if last_exon_strand == Strand::Forward {
            (self.exons.last().unwrap().reference_chromosome_id, self.exons.first().unwrap().reference_end)
        } else {
            (self.exons.first().unwrap().reference_chromosome_id, self.exons.first().unwrap().reference_end)
        }
    }

    pub fn get_reference_transcript_matches_map(&self) -> &HashMap<Box<str>, ReferenceTranscriptMatch> {
        &self.reference_transcript_matches_map
    }

    pub fn get_reference_start_position(&self) -> (u16, usize) {
        let first_exon_strand = self.exons.first().unwrap().reference_strand.clone();
        if first_exon_strand == Strand::Forward {
            (self.exons.first().unwrap().reference_chromosome_id, self.exons.first().unwrap().reference_start)
        } else {
            (self.exons.last().unwrap().reference_chromosome_id, self.exons.last().unwrap().reference_start)
        }
    }

    pub fn get_variant_records_map(&self) -> &HashMap<Vec<Box<str>>, Vec<VariantRecord>> {
        &self.variant_records_map
    }

    pub fn identify_variants(
        &mut self,
        reference_transcript_matches: &Vec<ReferenceTranscriptMatch>,
        gene_annotator: &(impl GeneAnnotator + Sync),
        reference_genome_fasta_file: &str,
        min_mapping_quality: usize,
        min_base_quality: u8
    ) -> &HashMap<Vec<Box<str>>, Vec<VariantRecord>> {
        self.variant_records_map.clear();

        // Step 1. Get a map of reference transcript IDs, matches, and sequences
        let mut reference_transcript_ids_map: HashMap<Box<str>, Vec<Box<str>>> = HashMap::new();
        let mut reference_transcript_match_map: HashMap<Box<str>, &ReferenceTranscriptMatch> = HashMap::new();
        let mut reference_transcript_sequences_map: HashMap<Box<str>, ReferenceTranscriptSequence> = HashMap::new();
        for reference_transcript_match in reference_transcript_matches.iter() {
            reference_transcript_ids_map
                .entry(reference_transcript_match.reference_gene_id.clone())
                .or_insert(Vec::new())
                .push(reference_transcript_match.reference_transcript_id.clone());
            reference_transcript_match_map.insert(
                reference_transcript_match.reference_transcript_id.clone(),
                reference_transcript_match
            );
            self.reference_transcript_matches_map.insert(
                reference_transcript_match.reference_transcript_id.clone(),
                reference_transcript_match.clone()
            );
            let reference_transcript: &Transcript = gene_annotator
                .get_transcript(&*reference_transcript_match.reference_transcript_id)
                .unwrap();
            let reference_transcript_sequence: ReferenceTranscriptSequence = ReferenceTranscriptSequence::from_reference_transcript(
                reference_transcript,
                &self.chromosome_names_map,
                reference_genome_fasta_file
            );
            reference_transcript_sequences_map.insert(
                reference_transcript_match.reference_transcript_id.clone(),
                reference_transcript_sequence
            );
        }

        // Step 2. Identify RNA variants based on contextualized alignment structures
        if reference_transcript_ids_map.is_empty() {
            let mut alignment_structure = self.alignment_structure.clone();

            alignment_structure.contextualize(
                &vec![],
                gene_annotator,
                &self.chromosome_names_map
            );

            let variant_records: Vec<VariantRecord> = alignment_structure.identify_variant_records(
                min_mapping_quality,
                min_base_quality,
                AnalyteType::RNA
            );

            self.contextualized_alignment_structures_map.insert(
                vec![],
                alignment_structure
            );

            self.variant_records_map.insert(
                vec![],
                variant_records
            );
        } else {
            // Get all combinations of reference transcripts (1 per gene ID)
            let mut keys: Vec<&str> = reference_transcript_ids_map.keys().map(|k| k.as_ref()).collect();
            keys.sort_unstable();
            let reference_transcript_ids_combinations: Vec<Vec<(Box<str>, Box<str>)>>= keys.iter()
                .map(|k| reference_transcript_ids_map[*k].iter())
                .multi_cartesian_product()
                .map(|choice| {
                    keys.iter()
                        .cloned()
                        .zip(choice.into_iter())
                        .map(|(k, v)| (k.into(), v.clone()))
                        .collect::<Vec<(Box<str>, Box<str>)>>()
                })
                .collect();

            for reference_transcript_ids_combination in reference_transcript_ids_combinations.iter() {
                let mut reference_transcript_sequences: Vec<&ReferenceTranscriptSequence> = Vec::new();
                let mut reference_transcript_ids: Vec<Box<str>> = Vec::new();
                for (reference_gene_id, reference_transcript_id) in reference_transcript_ids_combination.iter() {
                    reference_transcript_sequences.push(reference_transcript_sequences_map.get(reference_transcript_id).unwrap());
                    reference_transcript_ids.push(reference_transcript_id.clone());
                }
                let mut alignment_structure = self.alignment_structure.clone();
                alignment_structure.contextualize(
                    &reference_transcript_sequences,
                    gene_annotator,
                    &self.chromosome_names_map
                );
                let variant_records: Vec<VariantRecord> = alignment_structure.identify_variant_records(
                    min_mapping_quality,
                    min_base_quality, 
                    AnalyteType::RNA
                );
                self.contextualized_alignment_structures_map.insert(
                    reference_transcript_ids.clone(),
                    alignment_structure
                );
                self.variant_records_map.insert(
                    reference_transcript_ids.clone(),
                    variant_records
                );
            }
        }

        &self.variant_records_map
    }

    pub fn is_reference_transcript(&self) -> bool {
        for (reference_transcript_ids, alignment_structure) in self.contextualized_alignment_structures_map.iter() {
            let variant_records: &Vec<VariantRecord> = self.variant_records_map.get(reference_transcript_ids).unwrap();
            if reference_transcript_ids.len() == 1 && variant_records.is_empty() {
                return true;
            }
        }
        false
    }
    
    pub fn set_transcript_id(&mut self, transcript_id: usize) {
        self.id = transcript_id;
    }
}

impl Clone for TranscriptModel {
    fn clone(&self) -> Self {
        TranscriptModel {
            id: self.id,
            alignment_structure: self.alignment_structure.clone(),
            exons: self.exons.clone(),
            introns: self.introns.clone(),
            contextualized_alignment_structures_map: self.contextualized_alignment_structures_map.clone(),
            reference_transcript_matches_map: self.reference_transcript_matches_map.clone(),
            variant_records_map: self.variant_records_map.clone(),
            chromosome_names_map: self.chromosome_names_map.clone()
        }
    }
}
