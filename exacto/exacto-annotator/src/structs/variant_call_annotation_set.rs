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


use exacto_caller::prelude::*;
use polars::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::str::FromStr;

use crate::prelude::*;


#[derive(Debug, Serialize, Deserialize)]
pub struct VariantCallAnnotationSet {
    /// A map between variant call ID and its associated variant call annotation
    pub annotations: HashMap<usize, VariantCallAnnotation>,

    /// Nested structure for position indexing:
    /// - Outer HashMap: Maps chromosome names to their position index
    /// - Inner BTreeMap: Maps positions to variant call IDs
    position_index: HashMap<Box<str>, BTreeMap<u32, HashSet<usize>>>,

    /// A map between reference transcript ID and its associated variant call IDs
    transcript_index: HashMap<Box<str>, HashSet<usize>>
}

impl PartialEq for VariantCallAnnotationSet {
    fn eq(&self, other: &Self) -> bool {
        if self.annotations == other.annotations {
            true
        } else {
            false
        }
    }
}

impl VariantCallAnnotationSet {
    pub fn new() -> Self {
        Self {
            annotations: HashMap::new(),
            transcript_index: HashMap::new(),
            position_index: HashMap::new()
        }
    }

    pub fn add_annotation(&mut self, variant_annotation: VariantCallAnnotation) {
        // Add to transcript index
        for reference_transcript_ids in variant_annotation.position_1_annotation.transcript_ids.values() {
            for reference_transcript_id in reference_transcript_ids {
                self.transcript_index
                    .entry(reference_transcript_id.clone())
                    .or_insert(HashSet::new())
                    .insert(variant_annotation.variant_call_id);
            }
        }
        for reference_transcript_ids in variant_annotation.position_2_annotation.transcript_ids.values() {
            for reference_transcript_id in reference_transcript_ids {
                self.transcript_index
                    .entry(reference_transcript_id.clone())
                    .or_insert(HashSet::new())
                    .insert(variant_annotation.variant_call_id);
            }
        }

        // Add to position index
        self.position_index
            .entry(variant_annotation.chromosome_1.clone())
            .or_insert_with(BTreeMap::new)
            .entry(variant_annotation.position_1)
            .or_insert_with(HashSet::new)
            .insert(variant_annotation.variant_call_id);
        self.position_index
            .entry(variant_annotation.chromosome_2.clone())
            .or_insert_with(BTreeMap::new)
            .entry(variant_annotation.position_2)
            .or_insert_with(HashSet::new)
            .insert(variant_annotation.variant_call_id);

        // Add to self.annotations
        self.annotations.insert(variant_annotation.variant_call_id, variant_annotation);
    }

    pub fn get(&self, variant_call_id: usize) -> &VariantCallAnnotation {
        self.annotations.get(&variant_call_id).unwrap()
    }
    
    pub fn get_by_range(&self, chromosome: &str, start: u32, end: u32) -> Vec<&VariantCallAnnotation> {
        let mut result_ids = HashSet::new();

        if let Some(position_map) = self.position_index.get(chromosome) {
            for (_pos, variant_call_ids) in position_map.range(start..=end) {
                result_ids.extend(variant_call_ids.iter().cloned());
            }
        }

        result_ids
            .into_iter()
            .map(|id| self.get(id))
            .collect()
    }
    
    pub fn get_by_reference_transcript(&self, reference_transcript_id: &str) -> Option<Vec<&VariantCallAnnotation>> {
        let mut variant_call_annotations: Vec<&VariantCallAnnotation> = Vec::new();
        if self.transcript_index.contains_key(reference_transcript_id) == false {
            None
        } else {
            for variant_call_id in self.transcript_index.get(reference_transcript_id).unwrap() {
                variant_call_annotations.push(self.annotations.get(variant_call_id).unwrap());
            }
            Some(variant_call_annotations)
        }
    }
    
    pub fn read_tsv_file(tsv_file: &str) -> VariantCallAnnotationSet {
        let parse_options = CsvParseOptions::default()
            .with_separator(b'\t');
        let df = CsvReadOptions::default()
            .with_parse_options(parse_options)
            .with_has_header(true)
            .try_into_reader_with_file_path(Some(tsv_file.into()))
            .unwrap()
            .finish()
            .unwrap();

        let col_variant_call_id = df.column("variant_call_id").unwrap().i64().unwrap();
        let col_chromosome_1 = df.column("chromosome_1").unwrap().str().unwrap();
        let col_position_1 = df.column("position_1").unwrap().i64().unwrap();
        let col_chromosome_2 = df.column("chromosome_2").unwrap().str().unwrap();
        let col_position_2 = df.column("position_2").unwrap().i64().unwrap();
        let col_variant_type = df.column("variant_type").unwrap().str().unwrap();
        let col_variant_sequence = df.column("sequence").unwrap().str().unwrap();
        let col_position_1_genic_region = df.column("position_1_genic_region").unwrap().str().unwrap();
        let col_position_1_annotation = df.column("position_1_annotation").unwrap().str().unwrap();
        let col_position_2_genic_region = df.column("position_2_genic_region").unwrap().str().unwrap();
        let col_position_2_annotation = df.column("position_2_annotation").unwrap().str().unwrap();

        let mut variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::new();
        for i in 0 ..df.height() {
            let variant_call_id: usize = col_variant_call_id.get(i).unwrap() as usize;
            let chromosome_1: Box<str> = col_chromosome_1.get(i).unwrap().into();
            let position_1: u32 = col_position_1.get(i).unwrap() as u32;
            let chromosome_2: Box<str> = col_chromosome_2.get(i).unwrap().into();
            let position_2: u32 = col_position_2.get(i).unwrap() as u32;
            let variant_type: VariantType = VariantType::from_str(col_variant_type.get(i).unwrap()).unwrap();
            let variant_sequence: Box<str> = col_variant_sequence.get(i).unwrap().into();
            let position_1_annotation: PositionAnnotation = PositionAnnotation::from_string(col_position_1_annotation.get(i).unwrap());
            let position_2_annotation: PositionAnnotation = PositionAnnotation::from_string(col_position_2_annotation.get(i).unwrap());
            
            let variant_call_annotation: VariantCallAnnotation = VariantCallAnnotation::new(
                variant_call_id,
                &*chromosome_1,
                position_1,
                &*chromosome_2,
                position_2,
                variant_type,
                &*variant_sequence,
                position_1_annotation,
                position_2_annotation
            );

            variant_call_annotation_set.add_annotation(variant_call_annotation);
        }

        variant_call_annotation_set
    }

    pub fn to_dataframe(&self) -> DataFrame {
        let mut variant_call_id_values: Vec<u64> = Vec::new();
        let mut chromosome_1_values: Vec<&str> = Vec::new();
        let mut position_1_values: Vec<u32> = Vec::new();
        let mut chromosome_2_values: Vec<&str> = Vec::new();
        let mut position_2_values: Vec<u32> = Vec::new();
        let mut variant_type_values: Vec<&str> = Vec::new();
        let mut variant_sequence_values: Vec<&str> = Vec::new();
        let mut position_1_annotation_genic_region_values: Vec<String> = Vec::new();
        let mut position_1_annotation_values: Vec<String> = Vec::new();
        let mut position_2_annotation_genic_region_values: Vec<String> = Vec::new();
        let mut position_2_annotation_values: Vec<String> = Vec::new();
        
        for variant_call_annotation in self.annotations.values() {
            variant_call_id_values.push(variant_call_annotation.variant_call_id as u64);
            chromosome_1_values.push(&*variant_call_annotation.chromosome_1);
            position_1_values.push(variant_call_annotation.position_1);
            chromosome_2_values.push(&*variant_call_annotation.chromosome_2);
            position_2_values.push(variant_call_annotation.position_2);
            variant_type_values.push(variant_call_annotation.variant_type.as_str());
            variant_sequence_values.push(&*variant_call_annotation.variant_sequence);
            position_1_annotation_genic_region_values.push(variant_call_annotation.position_1_annotation.genic_region.as_str().to_string());
            position_1_annotation_values.push(variant_call_annotation.position_1_annotation.to_string());
            position_2_annotation_genic_region_values.push(variant_call_annotation.position_2_annotation.genic_region.as_str().to_string());
            position_2_annotation_values.push(variant_call_annotation.position_2_annotation.to_string());
        }

        DataFrame::new(vec![
            Column::from(Series::new("variant_call_id".into(), variant_call_id_values)),
            Column::from(Series::new("chromosome_1".into(), chromosome_1_values)),
            Column::from(Series::new("position_1".into(), position_1_values)),
            Column::from(Series::new("chromosome_2".into(), chromosome_2_values)),
            Column::from(Series::new("position_2".into(), position_2_values)),
            Column::from(Series::new("variant_type".into(), variant_type_values)),
            Column::from(Series::new("sequence".into(), variant_sequence_values)),
            Column::from(Series::new("position_1_genic_region".into(), position_1_annotation_genic_region_values)),
            Column::from(Series::new("position_1_annotation".into(), position_1_annotation_values)),
            Column::from(Series::new("position_2_genic_region".into(), position_2_annotation_genic_region_values)),
            Column::from(Series::new("position_2_annotation".into(), position_2_annotation_values))
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

impl Clone for VariantCallAnnotationSet {
    fn clone(&self) -> Self {
        VariantCallAnnotationSet {
            annotations: self.annotations.clone(),
            position_index: self.position_index.clone(),
            transcript_index: self.transcript_index.clone()
        }
    }
}
