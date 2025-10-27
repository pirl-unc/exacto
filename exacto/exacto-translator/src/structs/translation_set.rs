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


use polars::prelude::*;
use serde::{Serialize, Deserialize};
use std::collections::HashMap;

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct TranslationSet {
    pub translations: Vec<Translation>
}

impl TranslationSet {
    pub fn new() -> Self {
        Self {
            translations: Vec::new()
        }
    }

    pub fn add_translation(&mut self, translation: Translation) {
        self.translations.push(translation);
    }

    pub fn get_count(&self) -> usize {
        self.translations.len()
    }

    pub fn get_peptides_by_strategy(&self, strategy: TranslationStrategy) -> Vec<Translation> {
        let mut translations: Vec<Translation> = Vec::new();
        match strategy {
            TranslationStrategy::AllORFs => {
                for translated_peptides in self.translations.iter() {
                    translations.push(translated_peptides.clone());
                }
            },
            TranslationStrategy::LongestORF => {
                for translation in self.translations.iter() {
                    let mut translation_: Translation = Translation::new(translation.rna.clone());
                    translation_.add_peptide(translation.get_longest_orf_peptide().clone());
                    translations.push(translation_);
                }
            },
            _ => {
                panic!("Unsupported strategy: {}", strategy.as_str());
            }
        }

        translations
    }

    pub fn to_dataframe(&self, strategy: TranslationStrategy) -> DataFrame {
        let mut peptide_id: usize = 1;
        let mut peptide_ids_map: HashMap<&str, usize> = HashMap::new();

        let mut peptide_ids: Vec<u64> = Vec::new();
        let mut peptide_sequences: Vec<String> = Vec::new();
        let mut rna_ids: Vec<String> = Vec::new();
        let mut rna_sequences: Vec<String> = Vec::new();
        let mut orf_start_positions: Vec<u64> = Vec::new();
        let mut orf_end_positions: Vec<u64> = Vec::new();

        for translation in self.get_peptides_by_strategy(strategy).iter() {
            for peptide in translation.peptides.iter() {
                rna_ids.push(translation.rna.id.to_string());
                rna_sequences.push(translation.rna.sequence.to_string());
                peptide_sequences.push(peptide.sequence.to_string());
                orf_start_positions.push(peptide.orf_start as u64);
                orf_end_positions.push(peptide.orf_end as u64);

                let curr_peptide_id: usize;
                if peptide_ids_map.contains_key(&*peptide.sequence) {
                    curr_peptide_id = *peptide_ids_map.get(&*peptide.sequence).unwrap();
                } else {
                    curr_peptide_id = peptide_id;
                    peptide_ids_map.insert(peptide.sequence.as_ref(), peptide_id);
                    peptide_id += 1;
                }

                peptide_ids.push(curr_peptide_id as u64);
            }
        }

        // Convert Vec<usize> to UInt32Chunked and then to Series
        let peptide_ids_ca: UInt64Chunked = ChunkedArray::from_vec("peptide_id".into(), peptide_ids);
        let orf_start_positions_ca: UInt64Chunked = ChunkedArray::from_vec("orf_start".into(), orf_start_positions);
        let orf_end_positions_ca: UInt64Chunked = ChunkedArray::from_vec("orf_end".into(), orf_end_positions);

        // Create the DataFrame with Series
        let df = DataFrame::new(vec![
            Column::from(peptide_ids_ca.into_series()),
            Column::from(Series::new("peptide_sequence".into(), peptide_sequences)),
            Column::from(Series::new("rna_id".into(), rna_ids)),
            Column::from(Series::new("rna_sequence".into(), rna_sequences)),
            Column::from(orf_start_positions_ca.into_series()),
            Column::from(orf_end_positions_ca.into_series())
        ]).unwrap();

        df
    }
}

impl Clone for TranslationSet {
    fn clone(&self) -> Self {
        TranslationSet {
            translations: self.translations.clone()
        }
    }
}
