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
use noodles_core::{Region, Position};
use noodles_fasta::io::indexed_reader::Builder;
use std::str::FromStr;

use crate::prelude::*;


#[derive(Debug)]
pub struct ReferenceTranscriptSequence {
    reference_gene_id: Box<str>,
    reference_transcript_id: Box<str>,
    bases: Vec<ReferenceBase>
}

/// API methods
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
    
    pub fn get_chromosome_id(&self) -> u16 {
        self.bases[0].reference_chromosome_id
    }
    
    pub fn get_gene_id(&self) -> &str {
        &self.reference_gene_id
    }
    
    pub fn get_length(&self) -> usize {
        self.bases.len()
    }

    pub fn get_sequence(&self) -> Box<str> {
        let mut s: String = String::new();
        for reference_base in self.bases.iter() {
            s.push_str(reference_base.reference_nucleotide.as_str());
        }
        s.into()
    }

    /// Get introns (reference position, reference position).
    ///
    /// # Returns
    /// A vector of `(u16, usize, usize)` tuples representing intron positions (chromosome ID, start, end).
    pub fn get_introns(&self) -> Vec<(u16, u32, u32)> {
        let mut introns: Vec<(u16, u32, u32)> = Vec::new();
        for i in 0..(self.bases.len() - 1) {
            let curr_base: &ReferenceBase = self.get_base(i);
            let next_base: &ReferenceBase = self.get_base(i + 1);
            if curr_base.reference_position.abs_diff(next_base.reference_position) > 1 {
                if curr_base.reference_strand == Strand::Forward {
                    introns.push((curr_base.reference_chromosome_id, curr_base.reference_position + 1, next_base.reference_position - 1));
                } else {
                    introns.push((curr_base.reference_chromosome_id, next_base.reference_position + 1, curr_base.reference_position - 1));
                }
            }
        }
        introns.sort_by_key(|&(_, start, _)| start);
        introns
    }

    pub fn get_strand(&self) -> &Strand {
        &self.bases[0].reference_strand
    }
    
    pub fn get_transcript_id(&self) -> &str {
        &*self.reference_transcript_id
    }
    
    pub fn get_transcript_end(&self) -> u32 {
        let mut end: u32 = 0;
        for base in self.bases.iter() {
            if base.reference_position > end {
                end = base.reference_position
            }
        }
        end
    }

    pub fn get_transcript_start(&self) -> u32 {
        let mut start: u32 = u32::MAX;
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

// Static functions
impl ReferenceTranscriptSequence {
    pub fn from_reference_transcript(
        transcript: &Transcript,
        chromosome_names_map: &BiMap<Box<str>, u16>,
        reference_genome_fasta_file: &str
    ) -> ReferenceTranscriptSequence {
        let mut fasta_reader = Builder::default()
            .build_from_path(reference_genome_fasta_file)
            .unwrap();
        let mut reference_transcript_sequence: ReferenceTranscriptSequence = ReferenceTranscriptSequence::new(&*transcript.gene_id, &*transcript.transcript_id);
        let reference_chromosome_name: &str = &*transcript.chromosome;
        let reference_chromosome_id: u16 = *chromosome_names_map.get_by_left(&*transcript.chromosome).unwrap();
        for exon in transcript.get_sorted_exons() {
            if transcript.strand == Strand::Forward {
                for i in exon.start..=exon.end {
                    let position_start: Position = Position::try_from(i as usize).unwrap();
                    let position_end: Position = Position::try_from(i as usize).unwrap();
                    let region: Region = Region::new(reference_chromosome_name, position_start..=position_end);
                    let ref_record = fasta_reader.query(&region).unwrap();
                    let ref_sequence_bytes: &[u8] = ref_record.sequence().as_ref();
                    let sequence: &str = std::str::from_utf8(ref_sequence_bytes).expect("Failed to convert sequence to UTF-8");
                    let nucleotide: Nucleotide = Nucleotide::from_str(sequence).unwrap();
                    let reference_base: ReferenceBase = ReferenceBase::new(
                        reference_chromosome_id,
                        i,
                        nucleotide,
                        exon.strand.clone(),
                        Some(exon.gene_id.clone()),
                        Some(exon.transcript_id.clone()),
                        Some(exon.exon_id.clone())
                    );
                    reference_transcript_sequence.push(reference_base);
                }
            } else {
                for i in (exon.start..=exon.end).rev() {
                    let position_start: Position = Position::try_from(i as usize).unwrap();
                    let position_end: Position = Position::try_from(i as usize).unwrap();
                    let region: Region = Region::new(reference_chromosome_name, position_start..=position_end);
                    let ref_record = fasta_reader.query(&region).unwrap();
                    let ref_sequence_bytes: &[u8] = ref_record.sequence().as_ref();
                    let sequence: &str = std::str::from_utf8(ref_sequence_bytes).expect("Failed to convert sequence to UTF-8");
                    let mut nucleotide: Nucleotide = Nucleotide::from_str(sequence).unwrap();
                    nucleotide = nucleotide.complement();
                    let reference_base: ReferenceBase = ReferenceBase::new(
                        reference_chromosome_id,
                        i,
                        nucleotide,
                        exon.strand.clone(),
                        Some(exon.gene_id.clone()),
                        Some(exon.transcript_id.clone()),
                        Some(exon.exon_id.clone())
                    );
                    reference_transcript_sequence.push(reference_base);
                }
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
