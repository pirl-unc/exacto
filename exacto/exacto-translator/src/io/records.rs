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


use serde::{Deserialize, Serialize, Serializer};
use std::collections::HashSet;

use exacto_core::prelude::LIST_SEPARATOR;

use crate::common::error::TranslatorError;


fn serialize_optional_id_set<S>(value: &Option<HashSet<u32>>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer
{
    let rendered: String = match value {
        Some(ids) => {
            let mut sorted: Vec<u32> = ids.iter().copied().collect();
            sorted.sort_unstable();
            sorted.iter().map(|id| id.to_string()).collect::<Vec<String>>().join(LIST_SEPARATOR)
        }
        None => String::new()
    };
    serializer.serialize_str(&rendered)
}


#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct AssembledTranscriptSupportRecord {
    pub assembled_transcript_name: Box<str>,
    pub sequence: Box<str>,
    pub read_names: Box<str>
}


impl AssembledTranscriptSupportRecord {
    pub fn from_consensus(
        cluster_id: usize,
        consensus_sequence: Box<str>,
        read_names: Box<str>
    ) -> Self {
        Self {
            assembled_transcript_name: cluster_id.to_string().into_boxed_str(),
            sequence: consensus_sequence,
            read_names
        }
    }
}


#[derive(Debug, Clone, PartialEq)]
pub struct ReferenceStitchedSpans {
    pub assembled_transcript_name: Box<str>,
    pub intervals: Vec<(u32, u32)>,
    pub stitched_length: u32
}

impl ReferenceStitchedSpans {
    pub fn from_stitch_points(
        cluster_id: usize,
        five_prime_stitch_point: u32,
        three_prime_stitch_point: u32,
        stitched_length: u32
    ) -> Result<Self, TranslatorError> {
        if five_prime_stitch_point > three_prime_stitch_point || three_prime_stitch_point > stitched_length {
            return Err(TranslatorError::InvalidStitchPoints {
                transcript: cluster_id.to_string().into_boxed_str(),
                five_prime: five_prime_stitch_point,
                three_prime: three_prime_stitch_point,
                stitched_length
            });
        }
        let mut intervals: Vec<(u32, u32)> = Vec::new();
        if five_prime_stitch_point > 0 {
            intervals.push((0, five_prime_stitch_point));
        }
        if three_prime_stitch_point < stitched_length {
            intervals.push((three_prime_stitch_point, stitched_length));
        }
        Ok(Self {
            assembled_transcript_name: cluster_id.to_string().into_boxed_str(),
            intervals,
            stitched_length
        })
    }
}


#[derive(Debug, Serialize)]
pub struct NucleotideRecord {
    pub proteoform_id: u32,
    pub assembled_transcript_name: Box<str>,
    pub amino_acid_index: u32,
    pub amino_acid: Box<str>,
    pub codon_index: u8,
    pub nucleotide: Box<str>,
    pub is_amino_acid_variant: bool,
    pub is_nucleotide_variant: bool,
    pub assembled_transcript_read_position: u32,
    pub assembled_transcript_alignment_index: Option<u32>,
    pub assembled_transcript_variant_id: Option<u32>,
    pub assembled_transcript_variant: Option<Box<str>>,
    #[serde(serialize_with = "serialize_optional_id_set")]
    pub dna_variant_ids: Option<HashSet<u32>>,
    pub dna_variant: Option<Box<str>>,
    pub preceding_event_assembled_transcript_variant_id: Option<u32>,
    pub preceding_event_assembled_transcript_variant: Option<Box<str>>,
    #[serde(serialize_with = "serialize_optional_id_set")]
    pub preceding_event_dna_variant_ids: Option<HashSet<u32>>,
    pub preceding_event_dna_variant: Option<Box<str>>,
    /// Appended last so pre-existing TSVs keep their column positions.
    pub is_reference_stitched: bool
}


#[derive(Debug, Serialize)]
pub struct ProteoformRecord {
    pub proteoform_id: usize,
    pub amino_acid_sequence: String,
    pub amino_acid_sequence_length: usize,
    pub num_mutant_amino_acids: u32,
    pub assembled_transcript_name: Box<str>,
    pub assembled_transcript_sequence: Box<str>,
    pub assembled_transcript_sequence_length: usize,
    pub orf_start: u32,
    pub orf_end: u32,
    pub reference_gene_name: String,
    pub reference_transcript_id: String,
    pub mutant_amino_acid_intervals: String,
    pub assembled_transcript_variant_ids: String,
    pub assembled_transcript_variants: String,
    pub dna_variant_ids: String,
    pub dna_variants: String,
    pub assembled_transcript_read_names: Box<str>,
    pub num_assembled_transcript_read_names: usize,
    pub dna_variant_read_names: Box<str>,
    pub assembled_transcript_read_position: Box<str>,
    /// Appended last so pre-existing TSVs keep their column positions.
    /// Same `start:end` inclusive spelling as `mutant_amino_acid_intervals`.
    pub reference_stitched_amino_acid_intervals: String,
    pub num_reference_stitched_amino_acids: u32
}