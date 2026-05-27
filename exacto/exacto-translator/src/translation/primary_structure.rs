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

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct PrimaryStructure {
    pub id: u32,
    pub orf_start: u32,
    pub orf_end: u32,
    pub amino_acids: Vec<AminoAcid>
}

impl PrimaryStructure {
    pub fn new(
        id: u32,
        orf_start: u32,
        orf_end: u32,
        amino_acids: Vec<AminoAcid>
    ) -> Self {
        assert!(orf_start < orf_end);
        Self {
            id: id,
            orf_start: orf_start,
            orf_end: orf_end,
            amino_acids: amino_acids
        }
    }

    pub fn get_id(&self) -> u32 {
        self.id
    }

    pub fn get_orf_start(&self) -> u32 {
        self.orf_start
    }

    pub fn get_orf_end(&self) -> u32 {
        self.orf_end
    }

    pub fn get_amino_acids(&self) -> &Vec<AminoAcid> {
        &self.amino_acids
    }

    pub fn get_length(&self) -> usize {
        self.amino_acids.len()
    }
}

impl Clone for PrimaryStructure {
    fn clone(&self) -> Self {
        PrimaryStructure {
            id: self.id,
            orf_start: self.orf_start,
            orf_end: self.orf_end,
            amino_acids: self.amino_acids.clone()
        }
    }
}
