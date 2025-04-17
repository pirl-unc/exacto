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


use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use serde::{Serialize, Deserialize};

use crate::common::constants::Strands;
use crate::prelude::find_overlap;
use crate::structs::gene_annotation::cds::CDS;
use crate::structs::gene_annotation::exon::Exon;
use crate::structs::gene_annotation::intron::Intron;
use crate::structs::gene_annotation::utr::UTR;
use crate::structs::gene_annotation::start_codon::StartCodon;
use crate::structs::gene_annotation::stop_codon::StopCodon;


#[derive(Debug, Serialize, Deserialize)]
pub struct Transcript {
    pub gene_id: Box<str>,
    pub transcript_id: Box<str>,
    pub source: Box<str>,
    pub chromosome: Box<str>,
    pub start: u32,
    pub end: u32,
    pub strand: Strands,
    pub level: u16,
    pub transcript_name: Box<str>,
    pub transcript_type: Box<str>,
    pub transcript_support_level: Box<str>,
    pub exons: HashMap<Box<str>,Exon>,      // key = exon ID
    pub utrs: HashMap<Box<str>,UTR>,        // key = exon ID
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

impl Eq for Transcript {}

impl Hash for Transcript {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.gene_id.hash(state);
        self.transcript_id.hash(state);
        self.source.hash(state);
        self.chromosome.hash(state);
        self.start.hash(state);
        self.end.hash(state);
        self.strand.hash(state);
        self.level.hash(state);
        self.transcript_name.hash(state);
        self.transcript_type.hash(state);
        self.transcript_support_level.hash(state);
    }
}

impl Transcript {
    pub fn new(
        gene_id: &str,
        transcript_id: &str,
        source: &str,
        chromosome: &str,
        start: u32,
        end: u32,
        strand: Strands,
        level: u16,
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

    pub fn get_sorted_exons(&self) -> Vec<&Exon> {
        let mut exons: Vec<&Exon> = self.exons.values().collect();
        exons.sort_by_key(|exon| exon.exon_number);
        exons
    }

    pub fn get_introns(&self) -> Vec<Intron> {
        let mut introns: Vec<Intron> = Vec::new();
        let mut intron_number: u16 = 1;
        if self.strand == Strands::Forward {
            let mut intron_start: u32 = 0;
            for exon in self.get_sorted_exons() {
                if exon.exon_number == 1 {
                    intron_start = exon.end + 1;
                } else {
                    let intron_end: u32 = exon.start - 1;
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
            let mut intron_end: u32 = 0;
            for exon in self.get_sorted_exons() {
                if exon.exon_number == 1 {
                    intron_end = exon.start - 1;
                } else {
                    let intron_start: u32 = exon.end + 1;
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

    pub fn get_size(&self) -> u32 {
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

    pub fn vectorize_exons(
        &self,
        chromosome: Box<str>,
        start: u32,
        end: u32,
        aligned_value: i8,
        unaligned_value: i8
    ) -> Vec<i8> {
        let v_size: usize = (end - start + 1) as usize;
        let mut v: Vec<i8> = vec![unaligned_value; v_size];
        for exon in self.exons.values() {
            if chromosome == exon.chromosome {
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
