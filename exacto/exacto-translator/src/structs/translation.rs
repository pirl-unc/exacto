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


use serde::{Deserialize, Serialize};
use crate::structs::rna::RNA;
use crate::structs::peptide::Peptide;


#[derive(Debug,Serialize,Deserialize)]
pub struct Translation {
    pub rna: RNA,
    pub peptides: Vec<Peptide>
}

impl Translation {
    pub fn new(rna: RNA) -> Self {
        Self {
            rna: rna,
            peptides: Vec::new()
        }
    }

    pub fn add_peptide(&mut self, peptide: Peptide) {
        self.peptides.push(peptide);
    }

    pub fn get_longest_orf_peptide(&self) -> &Peptide {
        self.peptides
            .iter()
            .max_by_key(|p| p.get_length())
            .unwrap()
    }

    pub fn get_peptides_count(&self) -> usize {
        self.peptides.len()
    }
}

impl Clone for Translation {
    fn clone(&self) -> Self {
        Translation {
            rna: self.rna.clone(),
            peptides: self.peptides.clone()
        }
    }
}
