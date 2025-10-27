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


use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use crate::prelude::*;


#[derive(Debug, Serialize, Deserialize)]
pub struct Transcript {
    pub gene_id: Box<str>,
    pub transcript_id: Box<str>,
    pub source: Box<str>,
    pub chromosome: Box<str>,
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
    pub level: u8,
    pub transcript_name: Box<str>,
    pub transcript_type: Box<str>,
    pub transcript_support_level: Box<str>,
    pub exons: HashMap<Box<str>, Exon>,      // key = exon ID
    pub utrs: HashMap<Box<str>, UTR>,        // key = exon ID
    pub start_codons: Vec<StartCodon>,
    pub stop_codons: Vec<StopCodon>,
    pub coding_sequences: Vec<CDS>
}

impl PartialEq for Transcript {
    fn eq(&self, other: &Self) -> bool {
        self.gene_id == other.gene_id &&
            self.transcript_id == other.transcript_id &&
            self.source == other.source &&
            self.chromosome == other.chromosome &&
            self.start == other.start &&
            self.end == other.end &&
            self.strand == other.strand &&
            self.level == other.level &&
            self.transcript_name == other.transcript_name &&
            self.transcript_type == other.transcript_type &&
            self.transcript_support_level == other.transcript_support_level
    }
}

impl Transcript {
    pub fn new(
        gene_id: &str,
        transcript_id: &str,
        source: &str,
        chromosome: &str,
        start: usize,
        end: usize,
        strand: Strand,
        level: u8,
        transcript_name: &str,
        transcript_type: &str,
        transcript_support_level: &str
    ) -> Self {
        Self {
            gene_id: gene_id.to_string().into_boxed_str(),
            transcript_id: transcript_id.to_string().into_boxed_str(),
            source: source.to_string().into_boxed_str(),
            chromosome: chromosome.to_string().into_boxed_str(),
            start: start,
            end: end,
            strand: strand.clone(),
            level: level,
            transcript_name: transcript_name.to_string().into_boxed_str(),
            transcript_type: transcript_type.to_string().into_boxed_str(),
            transcript_support_level: transcript_support_level.to_string().into_boxed_str(),
            exons: HashMap::new(),
            utrs: HashMap::new(),
            start_codons: Vec::new(),
            stop_codons: Vec::new(),
            coding_sequences: Vec::new()
        }
    }

    pub fn add_cds(&mut self, cds: CDS) {
        self.coding_sequences.push(cds);
    }

    pub fn add_exon(&mut self, exon: Exon) {
        self.exons.insert(exon.exon_id.to_string().into_boxed_str(), exon);
    }

    pub fn add_start_codon(&mut self, start_codon: StartCodon) {
        self.start_codons.push(start_codon);
    }

    pub fn add_stop_codon(&mut self, stop_codon: StopCodon) {
        self.stop_codons.push(stop_codon);
    }

    pub fn add_utr(&mut self, utr: UTR) {
        self.utrs.insert(utr.exon_id.to_string().into_boxed_str(), utr);
    }

    pub fn get_exon(&self, exon_id: &str) -> Option<&Exon> {
        self.exons.get(exon_id)
    }
    
    pub fn get_exon_by_number(&self, exon_number: u32) -> Option<&Exon> {
        for exon in self.exons.values() {
            if exon.exon_number == exon_number {
                return Some(exon);
            }
        }
        None
    }

    pub fn get_sorted_exons(&self) -> Vec<&Exon> {
        let mut exons: Vec<&Exon> = self.exons.values().collect();
        exons.sort_by_key(|exon| exon.exon_number);
        exons
    }

    pub fn get_introns(&self) -> Vec<Intron> {
        let mut introns: Vec<Intron> = Vec::new();
        let mut intron_number: u32 = 1;
        if self.strand == Strand::Forward {
            let mut intron_start: usize = 0;
            for exon in self.get_sorted_exons() {
                if exon.exon_number == 1 {
                    intron_start = exon.end + 1;
                } else {
                    let intron_end: usize = exon.start - 1;
                    let intron: Intron = Intron::new(
                        &*self.gene_id,
                        &*self.transcript_id,
                        &*self.source,
                        &*self.chromosome,
                        intron_start,
                        intron_end,
                        self.strand.clone(),
                        intron_number
                    );
                    introns.push(intron);
                    intron_start = exon.end + 1;
                    intron_number += 1;
                }
            }
        } else {
            let mut intron_end: usize = 0;
            for exon in self.get_sorted_exons() {
                if exon.exon_number == 1 {
                    intron_end = exon.start - 1;
                } else {
                    let intron_start: usize = exon.end + 1;
                    let intron: Intron = Intron::new(
                        &*self.gene_id,
                        &*self.transcript_id,
                        &*self.source,
                        &*self.chromosome,
                        intron_start,
                        intron_end,
                        self.strand.clone(),
                        intron_number
                    );
                    introns.push(intron);
                    intron_end = exon.start - 1;
                    intron_number += 1;
                }
            }
        }

        introns
    }

    pub fn get_exon_ids(&self) -> Vec<Box<str>> {
        self.exons.keys().cloned().collect()
    }

    pub fn get_size(&self) -> usize {
        self.end - self.start + 1
    }

    pub fn get_5prime_utr(&self) -> &UTR {
        for utr in self.utrs.values() {
            if utr.exon_number == 1 {
                return utr;
            }
        }
        panic!("Could not find 5prime UTR for {}", self.transcript_id);
    }

    pub fn get_3prime_utr(&self) -> &UTR {
        for utr in self.utrs.values() {
            if utr.exon_number != 1 {
                return utr;
            }
        }
        panic!("Could not find 3prime UTR for {}", self.transcript_id);
    }

    /// Determines the genic context of a given genomic position with respect to the transcript.
    ///
    /// This method classifies the position as `Exonic`, `Intronic`, or `Intergenic`, based on
    /// its location relative to the transcript's exons. If the position is exonic or intronic,
    /// it also returns the corresponding exon(s) involved.
    ///
    /// # Arguments
    ///
    /// * `position` - A `u32` representing the genomic position to query.
    ///
    /// # Returns
    ///
    /// A tuple consisting of:
    /// - `GenicRegion`: An enum indicating the genic context (`Exonic`, `Intronic`, `Intergenic`).
    /// - `Option<(Exon, Exon)>`: 
    ///     - For `Exonic`, returns `Some((exon, exon))` where both elements are the same.
    ///     - For `Intronic`, returns `Some((exon_5p, exon_3p))` indicating the flanking exons, 
    ///       ordered according to the transcript's strand (5' to 3').
    ///     - For `Intergenic`, returns `None`.
    pub fn locate_position(&self, position: usize) -> (GenicRegion, Option<(Exon, Exon)>) {
        let mut sorted_exons: Vec<&Exon> = self.exons.values().collect();
        sorted_exons.sort_by_key(|e| e.exon_number);
        for i in 0..sorted_exons.len() {
            let curr_exon = sorted_exons[i];
            
            // Check if position is within exon bounds
            if position >= curr_exon.start && position <= curr_exon.end {
                return (GenicRegion::Exonic, Some((curr_exon.clone(), curr_exon.clone())));
            }

            // Check if position is in an intron
            if i > 0 {
                let prev_exon: &Exon = sorted_exons[i - 1];
                
                if self.strand == Strand::Forward {
                    if position > prev_exon.end && position < curr_exon.start {
                        return (GenicRegion::Intronic, Some((prev_exon.clone(), curr_exon.clone())));
                    }
                } else {
                    if position < prev_exon.start && position > curr_exon.end {
                        return (GenicRegion::Intronic, Some((prev_exon.clone(), curr_exon.clone())));
                    }
                }
            }
        }

        (GenicRegion::Intergenic, None)
    }
    
    pub fn vectorize_exons(
        &self,
        chromosome: Box<str>,
        start: usize,
        end: usize,
        aligned_value: i8,
        unaligned_value: i8
    ) -> Vec<i8> {
        let v_size: usize = end - start + 1;
        let mut v: Vec<i8> = vec![unaligned_value; v_size];
        for exon in self.exons.values() {
            if chromosome == exon.chromosome {
                match find_overlap((exon.start as isize, exon.end as isize), (start as isize, end as isize)) {
                    Some((x,y)) => {
                        for pos in x..=y {
                            let i = (pos as usize) - start;
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

impl Clone for Transcript {
    fn clone(&self) -> Self {
        Transcript {
            gene_id: self.gene_id.to_string().into_boxed_str(),
            transcript_id: self.transcript_id.to_string().into_boxed_str(),
            source: self.source.to_string().into_boxed_str(),
            chromosome: self.chromosome.to_string().into_boxed_str(),
            start: self.start,
            end: self.end,
            strand: self.strand.clone(),
            level: self.level,
            transcript_name: self.transcript_name.to_string().into_boxed_str(),
            transcript_type: self.transcript_type.to_string().into_boxed_str(),
            transcript_support_level: self.transcript_support_level.to_string().into_boxed_str(),
            exons: self.exons.clone(),
            utrs: self.utrs.clone(),
            start_codons: self.start_codons.clone(),
            stop_codons: self.stop_codons.clone(),
            coding_sequences: self.coding_sequences.clone()
        }
    }
}
