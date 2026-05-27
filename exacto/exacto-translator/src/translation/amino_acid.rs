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


use exacto_core::prelude::*;
use serde::{Deserialize,Serialize};

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct AminoAcid {
    pub index: u32,
    pub amino_acid: Box<str>,
    pub nucleotides: Vec<TranscriptNucleotide>
}

impl AminoAcid {
    pub fn new(
        index: u32,
        nucleotides: Vec<TranscriptNucleotide>,
    ) -> Self {
        assert!(nucleotides.len() == 3);
        let codon: String = nucleotides.iter()
            .map(|n| n.get_nucleotide().as_str())
            .collect::<String>()
            .to_uppercase()
            .replace('T', "U");
        let amino_acid: &str = CODON_TABLE[codon.as_str()];
        Self {
            index: index,
            amino_acid: amino_acid.into(),
            nucleotides: nucleotides
        }
    }

    pub fn get_amino_acid(&self) -> &str {
        &*self.amino_acid
    }

    pub fn get_index(&self) -> u32 {
        self.index
    }

    pub fn get_nucleotides(&self) -> &Vec<TranscriptNucleotide> {
        &self.nucleotides
    }

    pub fn is_variant(&self) -> bool {
        for nucleotide in self.nucleotides.iter() {
            if nucleotide.is_variant() {
                return true;
            }
        }
        false
    }
}

impl Clone for AminoAcid {
    fn clone(&self) -> Self {
        AminoAcid {
            index: self.index,
            amino_acid: self.amino_acid.clone(),
            nucleotides: self.nucleotides.clone()
        }
    }
}
