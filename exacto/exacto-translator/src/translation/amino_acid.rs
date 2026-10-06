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


use crate::prelude::*;


#[derive(Clone,Debug)]
pub struct AminoAcid {
    pub index: u32,
    pub amino_acid: char,
    pub nucleotides: [AssembledTranscriptNucleotide; 3]
}

impl AminoAcid {
    pub fn new(
        index: u32,
        amino_acid: char,
        nucleotides: [AssembledTranscriptNucleotide; 3]
    ) -> Self {
        Self {
            index: index,
            amino_acid: amino_acid,
            nucleotides: nucleotides
        }
    }

    pub fn get_amino_acid(&self) -> char {
        self.amino_acid
    }

    pub fn get_index(&self) -> u32 {
        self.index
    }

    pub fn get_nucleotides(&self) -> &[AssembledTranscriptNucleotide; 3] {
        &self.nucleotides
    }

    pub fn is_variant(&self) -> bool {
        self.nucleotides.iter().any(|nucleotide| nucleotide.is_variant())
    }

    pub fn is_reference_stitched(&self) -> bool {
        self.nucleotides.iter().any(|nucleotide| nucleotide.is_reference_stitched())
    }
}
