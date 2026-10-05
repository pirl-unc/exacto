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
use exacto_core::prelude as exacto_core;
use ::exacto_core::prelude::{ReferenceChromosomeID, ReferenceChromosomeName, ReferenceGeneID, ReferencePosition, ReferenceTranscriptID};
use std::str::FromStr;
use serde::{Deserialize, Serialize};

use crate::prelude::*;


#[derive(Debug, Serialize, Deserialize)]
pub struct ReferenceTranscriptSequence {
    reference_gene_id: ReferenceGeneID,
    reference_transcript_id: ReferenceTranscriptID,
    bases: Vec<ReferenceBase>
}

impl ReferenceTranscriptSequence {
    pub fn new(reference_gene_id: &str, reference_transcript_id: &str) -> Self {
        Self {
            reference_gene_id: reference_gene_id.into(),
            reference_transcript_id: reference_transcript_id.into(),
            bases: Vec::new()
        }
    }

    pub fn get_base(&self, i: usize) -> &ReferenceBase {
        self.bases.get(i).unwrap()
    }
    
    pub fn get_bases(&self) -> &Vec<ReferenceBase> {
        &self.bases
    }
    
    pub fn get_chromosome_id(&self) -> ReferenceChromosomeID {
        self.bases[0].reference_chromosome_id
    }
    
    pub fn get_gene_id(&self) -> &str {
        &self.reference_gene_id
    }
    
    pub fn get_length(&self) -> usize {
        self.bases.len()
    }

    pub fn get_reference_base_position(&self, reference_position: ReferencePosition) -> Option<usize> {
        let descending: bool = self.bases.len() > 1 && self.bases[0].reference_position > self.bases[1].reference_position;
        self.bases
            .binary_search_by(|base| {
                let ordering = base.reference_position.cmp(&reference_position);
                if descending { ordering.reverse() } else { ordering }
            })
            .ok()
    }

    pub fn get_sequence(&self) -> Box<str> {
        let mut s: String = String::new();
        for reference_base in self.bases.iter() {
            s.push_str(reference_base.reference_nucleotide.as_str());
        }
        s.into()
    }
    
    pub fn get_introns(&self) -> Vec<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)> {
        let mut introns: Vec<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)> = Vec::new();
        for i in 0..(self.bases.len() - 1) {
            let curr_base: &ReferenceBase = self.get_base(i);
            let next_base: &ReferenceBase = self.get_base(i + 1);
            if curr_base.reference_position.abs_diff(next_base.reference_position) > 1 {
                if curr_base.reference_strand == exacto_core::Strand::Forward {
                    introns.push(
                        (
                            curr_base.reference_chromosome_id,
                            curr_base.reference_position + 1,
                            next_base.reference_position - 1
                        )
                    );
                } else {
                    introns.push(
                        (
                            curr_base.reference_chromosome_id,
                            next_base.reference_position + 1,
                            curr_base.reference_position - 1
                        )
                    );
                }
            }
        }
        introns.sort_by_key(|&(_, start, _)| start);
        introns
    }

    pub fn get_strand(&self) -> &exacto_core::Strand {
        &self.bases[0].reference_strand
    }
    
    pub fn get_transcript_id(&self) -> &str {
        &*self.reference_transcript_id
    }
    
    pub fn get_transcript_end(&self) -> ReferencePosition {
        let mut end: ReferencePosition = 0;
        for base in self.bases.iter() {
            if base.reference_position > end {
                end = base.reference_position
            }
        }
        end
    }

    pub fn get_transcript_start(&self) -> ReferencePosition {
        let mut start: ReferencePosition = u32::MAX;
        for base in self.bases.iter() {
            if base.reference_position < start {
                start = base.reference_position
            }
        }
        start
    }

    pub fn push(&mut self, reference_base: ReferenceBase) {
        if self.get_length() > 0 {
            assert_eq!(
                reference_base.reference_strand,
                self.get_base(0).reference_strand,
                "Reference strand mismatch."
            );
            assert_eq!(
                reference_base.reference_chromosome_id,
                self.get_base(0).reference_chromosome_id,
                "Reference chromosome mismatch."
            )
        }
        self.bases.push(reference_base);
    }
}

impl ReferenceTranscriptSequence {
    pub fn from_reference_transcript(
        transcript: &exacto_core::Transcript,
        chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
        fasta_map: &exacto_core::FastaMap
    ) -> ReferenceTranscriptSequence {
        let mut reference_transcript_sequence: ReferenceTranscriptSequence = ReferenceTranscriptSequence::new(&*transcript.gene_id, &*transcript.transcript_id);
        let reference_chromosome_name: &str = &*transcript.chromosome;
        let reference_chromosome_id: ReferenceChromosomeID = *chromosome_names_map.get_by_left(&*transcript.chromosome).unwrap();
        let is_forward: bool = transcript.strand == exacto_core::Strand::Forward;
        for exon in transcript.get_sorted_exons() {
            let exon_sequence: Box<str> = fasta_map.get_sequence(
                reference_chromosome_name,
                exon.start as usize,
                exon.end as usize
            ).to_string().into_boxed_str();
            let exon_bytes: &[u8] = exon_sequence.as_bytes();
            assert_eq!(
                exon_bytes.len(),
                (exon.end - exon.start + 1) as usize,
                "Reference sequence length does not match the span of exon {}.",
                exon.exon_id
            );

            // Forward transcripts read the exon 5' -> 3'; reverse transcripts walk it
            // backwards and complement each base.
            let mut positions: Vec<ReferencePosition> = (exon.start..=exon.end).collect();
            if !is_forward {
                positions.reverse();
            }
            for position in positions {
                let offset: usize = (position - exon.start) as usize;
                let sequence: &str = std::str::from_utf8(&exon_bytes[offset..offset + 1])
                    .expect("Failed to convert sequence to UTF-8");
                let mut nucleotide: exacto_core::Nucleotide = exacto_core::Nucleotide::from_str(sequence)
                    .unwrap_or_else(|_| if sequence.chars().next().is_some_and(|c| c.is_ascii_lowercase()) {
                        exacto_core::Nucleotide::n
                    } else {
                        exacto_core::Nucleotide::N
                    });
                if !is_forward {
                    nucleotide = nucleotide.complement();
                }
                let reference_base: ReferenceBase = ReferenceBase::new(
                    reference_chromosome_id,
                    position,
                    nucleotide,
                    exon.strand.clone(),
                    Some(exon.gene_id.clone()),
                    Some(exon.transcript_id.clone()),
                    Some(exon.exon_id.clone())
                );
                reference_transcript_sequence.push(reference_base);
            }
        }
        reference_transcript_sequence
    }
}

impl Clone for ReferenceTranscriptSequence {
    fn clone(&self) -> Self {
        ReferenceTranscriptSequence {
            reference_gene_id: self.reference_gene_id.clone(),
            reference_transcript_id: self.reference_transcript_id.clone(),
            bases: self.bases.clone()
        }
    }
}
