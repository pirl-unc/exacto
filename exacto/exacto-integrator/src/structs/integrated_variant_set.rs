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
use std::collections::HashSet;
use std::fs::File;
use polars::prelude::*;

use crate::prelude::*;


#[derive(Debug, Serialize, Deserialize)]
pub struct IntegratedVariantSet {
    pub integrated_variants: HashSet<IntegratedVariant>

}

impl IntegratedVariantSet {
    pub fn new() -> Self {
        Self {
            integrated_variants: HashSet::new()
        }
    }

    pub fn add_integrated_variant(&mut self, integrated_variant: IntegratedVariant) {
        self.integrated_variants.insert(integrated_variant);
    }
    
    pub fn get_size(&self) -> usize {
        self.integrated_variants.len()
    }

    pub fn to_dataframe(&self) -> DataFrame {
        let mut transcript_model_ids: Vec<u64> = Vec::new();
        let mut reference_transcript_ids: Vec<String> = Vec::new();
        let mut rna_variant_call_ids: Vec<u64> = Vec::new();
        let mut dna_variant_call_ids: Vec<u64> = Vec::new();
        let mut distances: Vec<u64> = Vec::new();
        let mut rna_variant_positions: Vec<String> = Vec::new();
        let mut dna_variant_positions: Vec<String> = Vec::new();

        for integrated_variant in self.integrated_variants.iter() {
            for (dna_variant_call_id, integrated_variant_distance) in integrated_variant.dna_variant_call_ids.iter() {
                let reference_transcript_ids_string: String = integrated_variant.reference_transcript_ids.iter()
                    .map(|s| s.as_ref()) // &Box<str> -> &str
                    .collect::<Vec<&str>>()
                    .join(",");
                transcript_model_ids.push(integrated_variant.transcript_model_id as u64);
                reference_transcript_ids.push(reference_transcript_ids_string);
                rna_variant_call_ids.push(integrated_variant.rna_variant_call_id as u64);
                dna_variant_call_ids.push(*dna_variant_call_id as u64);
                distances.push(integrated_variant_distance.distance as u64);
                rna_variant_positions.push(integrated_variant_distance.rna_variant_position_used.as_str().to_string());
                dna_variant_positions.push(integrated_variant_distance.dna_variant_position_used.as_str().to_string());
            }
        }

        DataFrame::new(vec![
            Column::from(Series::new("transcript_model_id".into(), transcript_model_ids)),
            Column::from(Series::new("reference_transcript_ids".into(), reference_transcript_ids)),
            Column::from(Series::new("rna_variant_call_id".into(), rna_variant_call_ids)),
            Column::from(Series::new("dna_variant_call_id".into(), dna_variant_call_ids)),
            Column::from(Series::new("distance".into(), distances)),
            Column::from(Series::new("rna_variant_position".into(), rna_variant_positions)),
            Column::from(Series::new("dna_variant_position".into(), dna_variant_positions))
        ]).unwrap()
    }

    pub fn to_tsv_file(&self, output_file: &str) {
        let mut df = self.to_dataframe();
        let mut file = File::create(output_file).unwrap();
        CsvWriter::new(&mut file)
            .include_header(true)
            .with_separator(b'\t')
            .finish(&mut df)
            .unwrap();
    }
}

impl Clone for IntegratedVariantSet {
    fn clone(&self) -> Self {
        IntegratedVariantSet {
            integrated_variants: self.integrated_variants.clone()
        }
    }
}
