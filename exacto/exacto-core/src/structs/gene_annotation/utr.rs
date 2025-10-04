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

use crate::common::constants::Strand;


#[derive(Debug, Serialize, Deserialize)]
pub struct UTR {
    pub gene_id: Box<str>,
    pub transcript_id: Box<str>,
    pub exon_id: Box<str>,
    pub source: Box<str>,
    pub chromosome: Box<str>,
    pub start: u32,
    pub end: u32,
    pub strand: Strand,
    pub level: u8,
    pub exon_number: u16
}

impl UTR {
    pub fn new(
        gene_id: &str,
        transcript_id: &str,
        exon_id: &str,
        source: &str,
        chromosome: &str,
        start: u32,
        end: u32,
        strand: Strand,
        level: u8,
        exon_number: u16
    ) -> Self {
        Self {
            gene_id: gene_id.to_string().into_boxed_str(),
            transcript_id: transcript_id.to_string().into_boxed_str(),
            exon_id: exon_id.to_string().into_boxed_str(),
            source: source.to_string().into_boxed_str(),
            chromosome: chromosome.to_string().into_boxed_str(),
            start: start,
            end: end,
            strand: strand.clone(),
            level: level,
            exon_number: exon_number
        }
    }

    pub fn get_size(&self) -> u32 {
        self.end - self.start + 1
    }
}

impl Clone for UTR {
    fn clone(&self) -> Self {
        UTR {
            gene_id: self.gene_id.to_string().into_boxed_str(),
            transcript_id: self.transcript_id.to_string().into_boxed_str(),
            exon_id: self.exon_id.to_string().into_boxed_str(),
            source: self.source.to_string().into_boxed_str(),
            chromosome: self.chromosome.to_string().into_boxed_str(),
            start: self.start,
            end: self.end,
            strand: self.strand.clone(),
            level: self.level,
            exon_number: self.exon_number
        }
    }
}
