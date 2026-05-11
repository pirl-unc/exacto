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
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct PrimaryStructure {
    pub records: Vec<PrimaryStructureRecord>
}

impl PrimaryStructure {
    pub fn new() -> Self {
        Self {
            records: Vec::new()
        }
    }

    pub fn add_record(&mut self, record: PrimaryStructureRecord) {
        self.records.push(record);
    }

    pub fn to_record(&self) -> HashMap<Box<str>, Vec<AnyValue>> {
        let mut records_map: HashMap<Box<str>, Vec<AnyValue>> = HashMap::new();
        for record in self.records.iter() {
            let reference_transcript_ids: String = record.get_reference_transcript_ids()
                .iter().map(|s| s.as_ref())
                .collect::<Vec<&str>>()
                .join(",");
            let rna_variant_call_ids: String = record.get_rna_variant_call_ids()
                .iter()
                .map(|id| id.to_string())  // convert each usize to a String
                .collect::<Vec<String>>()  // collect into Vec<String>
                .join(",");
            let dna_variant_call_ids: String = record.get_dna_variant_call_ids()
                .iter()
                .map(|id| id.to_string())
                .collect::<Vec<String>>()
                .join(",");
            let codon_rna_variant_call_ids: String = record.get_codon_rna_variant_call_ids()
                .iter()
                .map(|id| id.to_string())
                .collect::<Vec<String>>()
                .join(",");
            let codon_dna_variant_call_ids: String = record.get_codon_dna_variant_call_ids()
                .iter()
                .map(|id| id.to_string())
                .collect::<Vec<String>>()
                .join(",");
            let frameshift_rna_variant_call_ids: String = record.get_frameshift_rna_variant_call_ids()
                .iter()
                .map(|id| id.to_string())
                .collect::<Vec<String>>()
                .join(",");
            let frameshift_dna_variant_call_ids: String = record.get_frameshift_dna_variant_call_ids()
                .iter()
                .map(|id| id.to_string())
                .collect::<Vec<String>>()
                .join(",");
            records_map.entry("primary_structure_index".into()).or_insert_with(Vec::new).push(AnyValue::UInt64(record.get_index() as u64));
            records_map.entry("type".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_record_type().as_str().into()));
            records_map.entry("amino_acid".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_amino_acid().as_ref().map_or("".into(),|k| k.to_string().into())));
            records_map.entry("codon_index".into()).or_insert_with(Vec::new).push(record.get_codon_index().as_ref().map_or(AnyValue::Int32(-1), |k| AnyValue::Int32(*k as i32)));
            records_map.entry("nucleotide".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_nucleotide().as_ref().map_or("".into(),|k| k.as_str().into())));
            records_map.entry("transcript_model_id".into()).or_insert_with(Vec::new).push(AnyValue::UInt64(record.get_transcript_model_id() as u64));
            records_map.entry("reference_transcript_ids".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(reference_transcript_ids.into()));
            records_map.entry("transcript_structure_index".into()).or_insert_with(Vec::new).push(AnyValue::UInt64(record.get_transcript_structure_index() as u64));
            records_map.entry("read_start".into()).or_insert_with(Vec::new).push(AnyValue::UInt64(record.get_read_start() as u64));
            records_map.entry("read_end".into()).or_insert_with(Vec::new).push(AnyValue::UInt64(record.get_read_end() as u64));
            records_map.entry("net_variant_nucleotides_count".into()).or_insert_with(Vec::new).push(AnyValue::Int64(record.get_net_variant_nucleotides_count() as i64));
            records_map.entry("frameshift_state".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_frameshift_state().as_ref().map_or("".into(),|k| k.as_str().into())));
            records_map.entry("rna_variant_call_ids".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(rna_variant_call_ids.into()));
            records_map.entry("dna_variant_call_ids".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(dna_variant_call_ids.into()));
            records_map.entry("codon_rna_variant_call_ids".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(codon_rna_variant_call_ids.into()));
            records_map.entry("codon_dna_variant_call_ids".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(codon_dna_variant_call_ids.into()));
            records_map.entry("frameshift_rna_variant_call_ids".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(frameshift_rna_variant_call_ids.into()));
            records_map.entry("frameshift_dna_variant_call_ids".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(frameshift_dna_variant_call_ids.into()));
            records_map.entry("amino_acid_change".into()).or_insert_with(Vec::new).push(AnyValue::StringOwned(record.get_amino_acid_change().as_ref().map_or("".into(),|k| k.as_str().into())));
        }
        records_map
    }
}

impl Clone for PrimaryStructure {
    fn clone(&self) -> Self {
        PrimaryStructure {
            records: self.records.clone()
        }
    }
}
