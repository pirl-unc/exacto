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


#[derive(Debug,Serialize,Deserialize)]
pub struct MutantPeptide {
    pub id: Box<str>,
    pub peptide_sequence: Box<str>,
    pub peptide_sequence_read_names: Vec<Box<str>>,
    pub k: usize,
    pub kmer_mask: Vec<bool>,
    pub rna_variants: Vec<Box<str>>,
    pub rna_variants_read_names: Vec<Box<str>>,
    pub dna_variants: Vec<Box<str>>,
    pub dna_variants_read_names: Vec<Box<str>>
}

impl MutantPeptide {
    pub fn new(
        id: Box<str>,
        peptide_sequence: Box<str>,
        peptide_sequence_read_names: Vec<Box<str>>,
        k: usize,
        kmer_mask: Vec<bool>,
        rna_variants: Vec<Box<str>>,
        rna_variants_read_names: Vec<Box<str>>,
        dna_variants: Vec<Box<str>>,
        dna_variants_read_names: Vec<Box<str>>
    ) -> Self {
        Self {
            id,
            peptide_sequence,
            peptide_sequence_read_names,
            k,
            kmer_mask,
            rna_variants,
            rna_variants_read_names,
            dna_variants,
            dna_variants_read_names
        }
    }

    pub fn get_kmer_mask_string(&self) -> String {
        let bitmask: String = self.kmer_mask
            .iter()
            .map(|b| if *b { '1'.to_string() } else { '0'.to_string() })
            .collect();
        bitmask
    }

    pub fn get_peptide_sequence_read_names_string(&self) -> String {
        self.peptide_sequence_read_names.join(",").to_string()
    }

    pub fn get_rna_variants_string(&self) -> String {
        self.rna_variants.join(",").to_string()
    }

    pub fn get_rna_variants_read_names_string(&self) -> String {
        self.rna_variants_read_names.join(",").to_string()
    }

    pub fn get_dna_variants_string(&self) -> String {
        self.dna_variants.join(",").to_string()
    }

    pub fn get_dna_variants_read_names_string(&self) -> String {
        self.dna_variants_read_names.join(",").to_string()
    }

    pub fn to_tsv_string(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            self.id,
            self.peptide_sequence,
            self.get_peptide_sequence_read_names_string(),
            self.k,
            self.get_kmer_mask_string(),
            self.get_rna_variants_string(),
            self.get_rna_variants_read_names_string(),
            self.get_dna_variants_string(),
            self.get_dna_variants_read_names_string()
        )
    }
}
