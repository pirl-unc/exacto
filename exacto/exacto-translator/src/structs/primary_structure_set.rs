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


use std::collections::HashMap;
use std::fs::File;
use polars::prelude::*;
use serde::{Deserialize,Serialize};
use exacto_core::prelude::write_fasta_file;
use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct PrimaryStructureSet {
    pub primary_structures: Vec<PrimaryStructure>
}

impl PrimaryStructureSet {
    pub fn new() -> Self {
        Self {
            primary_structures: Vec::new()
        }
    }

    pub fn add(&mut self, primary_structure: PrimaryStructure) {
        self.primary_structures.push(primary_structure);
    }

    pub fn to_dataframe(&self) -> DataFrame {
        let mut merged_map: HashMap<Box<str>, Vec<AnyValue>> = HashMap::new();
        merged_map.insert("peptide_id".into(), Vec::new());
        let mut peptide_id: u64 = 1;
        for primary_structure in self.primary_structures.iter() {
            let primary_structure_records: HashMap<Box<str>, Vec<AnyValue>> = primary_structure.to_record();
            let mut length: usize = 0;
            for (key, values) in primary_structure_records.iter() {
                merged_map
                    .entry(key.clone())
                    .or_default()
                    .extend(values.clone());
                length = values.len();
            }
            let peptide_ids: Vec<AnyValue> = vec![AnyValue::StringOwned(peptide_id.to_string().into()); length];
            merged_map
                .entry("peptide_id".into())
                .or_default()
                .extend(peptide_ids);
            peptide_id += 1;
        }
        
        DataFrame::new(vec![
            Column::from(Series::new("peptide_id".into(), merged_map.get("peptide_id").unwrap())),
            Column::from(Series::new("primary_structure_index".into(), merged_map.get("primary_structure_index").unwrap())),
            Column::from(Series::new("type".into(), merged_map.get("type").unwrap())),
            Column::from(Series::new("amino_acid".into(), merged_map.get("amino_acid").unwrap())),
            Column::from(Series::new("codon_index".into(), merged_map.get("codon_index").unwrap())),
            Column::from(Series::new("nucleotide".into(), merged_map.get("nucleotide").unwrap())),
            Column::from(Series::new("transcript_model_id".into(), merged_map.get("transcript_model_id").unwrap())),
            Column::from(Series::new("reference_transcript_ids".into(), merged_map.get("reference_transcript_ids").unwrap())),
            Column::from(Series::new("transcript_structure_index".into(), merged_map.get("transcript_structure_index").unwrap())),
            Column::from(Series::new("read_start".into(), merged_map.get("read_start").unwrap())),
            Column::from(Series::new("read_end".into(), merged_map.get("read_end").unwrap())),
            Column::from(Series::new("net_variant_nucleotides_count".into(), merged_map.get("net_variant_nucleotides_count").unwrap())),
            Column::from(Series::new("frameshift_state".into(), merged_map.get("frameshift_state").unwrap())),
            Column::from(Series::new("rna_variant_call_ids".into(), merged_map.get("rna_variant_call_ids").unwrap())),
            Column::from(Series::new("dna_variant_call_ids".into(), merged_map.get("dna_variant_call_ids").unwrap()))
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

    pub fn to_fasta_file(&self, output_file: &str) {
        let mut peptide_id: u64 = 1;
        let mut sequences: Vec<(Box<str>, Box<str>)> = Vec::new();
        for primary_structure in self.primary_structures.iter() {
            let mut sequence: String = String::new();
            for record in primary_structure.records.iter() {
                if *record.get_record_type() == PrimaryStructureRecordType::Base {
                    if record.get_codon_index().unwrap() == 0 {
                        sequence.push_str(record.get_amino_acid().as_ref().unwrap());
                    }
                }
            }
            sequences.push((peptide_id.to_string().into(), sequence.into()));
            peptide_id += 1;
        }
        write_fasta_file(&sequences, output_file);
    }
}

impl Clone for PrimaryStructureSet {
    fn clone(&self) -> Self {
        PrimaryStructureSet {
            primary_structures: self.primary_structures.clone()
        }
    }
}
