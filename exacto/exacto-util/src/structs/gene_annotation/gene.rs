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

use crate::common::constants::Strands;
use crate::structs::gene_annotation::transcript::Transcript;


#[derive(Debug, Serialize, Deserialize)]
pub struct Gene {
    pub gene_id: Box<str>,
    pub source: Box<str>,
    pub assembly: Box<str>,
    pub chromosome: Box<str>,
    pub start: u32,
    pub end: u32,
    pub strand: Strands,
    pub gene_name: Box<str>,
    pub level: u16,
    pub gene_type: Box<str>,
    pub transcripts: HashMap<Box<str>,Transcript>
}

impl Gene {
    pub fn new(
        gene_id: &str,
        source: &str,
        assembly: &str,
        chromosome: &str,
        start: u32,
        end: u32,
        strand: Strands,
        gene_name: &str,
        level: u16,
        gene_type: &str
    ) -> Self {
        Self {
            gene_id: gene_id.to_string().into_boxed_str(),
            source: source.to_string().into_boxed_str(),
            assembly: assembly.to_string().into_boxed_str(),
            chromosome: chromosome.to_string().into_boxed_str(),
            start: start,
            end: end,
            strand: strand.clone(),
            gene_name: gene_name.to_string().into_boxed_str(),
            level: level,
            gene_type: gene_type.to_string().into_boxed_str(),
            transcripts: HashMap::new()
        }
    }

    pub fn add_transcript(&mut self, transcript: Transcript) {
        self.transcripts.insert(transcript.transcript_id.clone(), transcript);
    }

    pub fn get_size(&self) -> u32 {
        self.end - self.start + 1
    }

    pub fn get_transcript(&self, transcript_id: &str) -> Option<&Transcript> {
        self.transcripts.get(transcript_id)
    }

    pub fn get_mut_transcript(&mut self, transcript_id: &str) -> Option<&mut Transcript> {
        self.transcripts.get_mut(transcript_id)
    }

    pub fn get_transcript_ids(&self) -> Vec<Box<str>> {
        self.transcripts.keys().cloned().collect()
    }
}

impl Clone for Gene {
    fn clone(&self) -> Self {
        Gene {
            gene_id: self.gene_id.to_string().into_boxed_str(),
            source: self.source.to_string().into_boxed_str(),
            assembly: self.assembly.to_string().into_boxed_str(),
            chromosome: self.chromosome.to_string().into_boxed_str(),
            start: self.start,
            end: self.end,
            strand: self.strand.clone(),
            gene_name: self.gene_name.to_string().into_boxed_str(),
            level: self.level,
            gene_type: self.gene_type.to_string().into_boxed_str(),
            transcripts: self.transcripts.clone()
        }
    }
}
