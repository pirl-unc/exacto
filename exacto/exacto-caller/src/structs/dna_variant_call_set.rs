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


use bimap::BiMap;
use polars::prelude::*;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use serde::{Serialize, Deserialize};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct DNAVariantCallSet {
    pub variant_calls: HashMap<usize, VariantCall>,
    
    /// A bidirectional map between read names and read IDs.
    ///
    /// - Left: `Box<str>` - human-readable read name (e.g., from FASTQ or BAM)
    /// - Right: `usize` - internal numeric read ID
    pub read_names_map: BiMap<Box<str>, usize>,

    /// A bidirectional map between chromosome names and internal chromosome IDs.
    ///
    /// - Left: `Box<str>` - chromosome name (e.g., "chr1", "chrX")
    /// - Right: `u16` - numeric chromosome ID used internally
    pub chromosome_names_map: BiMap<Box<str>, u16>,

    /// Nested structure for position indexing:
    /// - Outer HashMap: Maps chromosome IDs to their position index
    /// - Inner BTreeMap: Maps positions to variant call IDs
    position_index: HashMap<u16, BTreeMap<u32, HashSet<usize>>>
}

impl DNAVariantCallSet {
    pub fn new() -> Self {
        Self {
            variant_calls: HashMap::new(),
            position_index: HashMap::new(),
            read_names_map: BiMap::new(),
            chromosome_names_map: BiMap::new()
        }
    }

    pub fn add_variant_call(&mut self, variant_call: VariantCall) {
        if self.variant_calls.contains_key(&variant_call.id) == false {
            let consensus_record = variant_call.get_consensus_record().0;
            let chromosome_1_id: u16 = consensus_record.get_chromosome_1();
            let chromosome_2_id: u16 = consensus_record.get_chromosome_2();
            let position_1: u32 = consensus_record.graph_operation.get_position_1();
            let position_2: u32 = consensus_record.graph_operation.get_position_2();

            // Index the variant by its first breakpoint position
            self.position_index
                .entry(chromosome_1_id)
                .or_insert_with(BTreeMap::new)
                .entry(position_1)
                .or_insert_with(HashSet::new)
                .insert(variant_call.id);

            // Index the variant by its second breakpoint position
            self.position_index
                .entry(chromosome_2_id)
                .or_insert_with(BTreeMap::new)
                .entry(position_2)
                .or_insert_with(HashSet::new)
                .insert(variant_call.id);

            // Add to self.variant_calls
            self.variant_calls.insert(variant_call.id, variant_call);
        }
    }

    pub fn get_size(&self) -> usize {
        self.variant_calls.len()
    }

    pub fn get_variant_calls(&self) -> Vec<&VariantCall> {
        self.variant_calls.values().collect()
    }

    pub fn get_variant_calls_by_range(&self, chromosome_id: u16, start: u32, end: u32) -> Vec<&VariantCall> {
        let mut result_ids = HashSet::new();

        if let Some(position_map) = self.position_index.get(&chromosome_id) {
            for (_pos, variant_call_ids) in position_map.range(start..=end) {
                result_ids.extend(variant_call_ids.iter().cloned());
            }
        }

        result_ids.iter()
            .filter_map(|id| self.variant_calls.get(id))
            .collect()
    }

    pub fn get_variant_records(&self) -> Vec<&VariantRecord> {
        let mut variant_records: Vec<&VariantRecord> = Vec::new();
        for variant_call in self.variant_calls.values() {
            for variant_record in variant_call.variant_records.iter() {
                variant_records.push(variant_record);
            }
        }
        variant_records
    }

    pub fn load_chromosome_names(&mut self, chromosome_names_map: BiMap<Box<str>,u16>) {
        self.chromosome_names_map = chromosome_names_map;
    }

    pub fn load_read_names(&mut self, read_names_map: BiMap<Box<str>,usize>) {
        self.read_names_map = read_names_map;
    }

    pub fn remove_variant_call(&mut self, variant_call: &VariantCall) {
        self.variant_calls.remove(&variant_call.id);
    }

    pub fn to_dataframe(&self, num_threads: usize) -> DataFrame {
        assert!(
            !self.chromosome_names_map.is_empty(),
            "self.chromosome_names_map is empty."
        );
        assert!(
            !self.read_names_map.is_empty(),
            "self.read_names_map is empty."
        );

        let thread_pool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();

        let mut variant_calls: Vec<&VariantCall> = Vec::new();
        for variant_call in self.variant_calls.values() {
            variant_calls.push(variant_call);
        }
        let chunk_size = (variant_calls.len() + num_threads - 1) / num_threads;
        let rows: Vec<_> = thread_pool.install(|| {
            variant_calls
                .par_chunks(chunk_size)
                .flat_map_iter(|chunk| {
                    chunk.iter().map(|variant_call| {
                        let (consensus_record, consensus_read_names) = variant_call.get_named_consensus_record(&self.read_names_map);
                        let chromosome_1: &str = &*self
                            .chromosome_names_map
                            .get_by_right(&consensus_record.get_chromosome_1())
                            .unwrap();
                        let chromosome_2: &str = &*self
                            .chromosome_names_map
                            .get_by_right(&consensus_record.get_chromosome_2())
                            .unwrap();
                        let read_names: Vec<&str> = variant_call
                            .get_read_ids()
                            .iter()
                            .map(|read_id| &**self.read_names_map.get_by_right(read_id).unwrap())
                            .collect();
                        (
                            variant_call.id as u64,
                            chromosome_1,
                            consensus_record.graph_operation.get_position_1(),
                            consensus_record.graph_operation.get_strand_1().as_str(),
                            consensus_record.graph_operation.get_operation_type_1().as_str(),
                            chromosome_2,
                            consensus_record.graph_operation.get_position_2(),
                            consensus_record.graph_operation.get_strand_2().as_str(),
                            consensus_record.graph_operation.get_operation_type_2().as_str(),
                            consensus_record.get_variant_size() as i64,
                            consensus_record.get_variant_type().as_str().to_string(),
                            consensus_record.get_standardized_sequence(),
                            consensus_read_names.join(","),
                            consensus_read_names.len() as u32,
                            read_names.join(","),
                            read_names.len() as u32
                        )
                    })
                })
                .collect()
        });

        DataFrame::new(vec![
            Column::from(Series::new("variant_call_id".into(), rows.iter().map(|r| r.0).collect::<Vec<_>>())),
            Column::from(Series::new("chromosome_1".into(), rows.iter().map(|r| r.1).collect::<Vec<_>>())),
            Column::from(Series::new("position_1".into(), rows.iter().map(|r| r.2).collect::<Vec<_>>())),
            Column::from(Series::new("strand_1".into(), rows.iter().map(|r| r.3).collect::<Vec<_>>())),
            Column::from(Series::new("operation_1".into(), rows.iter().map(|r| r.4).collect::<Vec<_>>())),
            Column::from(Series::new("chromosome_2".into(), rows.iter().map(|r| r.5).collect::<Vec<_>>())),
            Column::from(Series::new("position_2".into(), rows.iter().map(|r| r.6).collect::<Vec<_>>())),
            Column::from(Series::new("strand_2".into(), rows.iter().map(|r| r.7).collect::<Vec<_>>())),
            Column::from(Series::new("operation_2".into(), rows.iter().map(|r| r.8).collect::<Vec<_>>())),
            Column::from(Series::new("variant_size".into(), rows.iter().map(|r| r.9).collect::<Vec<_>>())),
            Column::from(Series::new("variant_type".into(), rows.iter().map(|r| r.10.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("variant_sequence".into(), rows.iter().map(|r| r.11.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("consensus_read_names".into(), rows.iter().map(|r| r.12.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("consensus_read_names_count".into(), rows.iter().map(|r| r.13).collect::<Vec<_>>())),
            Column::from(Series::new("read_names".into(), rows.iter().map(|r| r.14.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("read_names_count".into(), rows.iter().map(|r| r.15).collect::<Vec<_>>()))
        ])
        .unwrap()
    }

    pub fn to_tsv_file(&self, output_file: &str, num_threads: usize) {
        let mut df = self.to_dataframe(num_threads);
        let mut file = File::create(output_file).unwrap();
        CsvWriter::new(&mut file)
            .include_header(true)
            .with_separator(b'\t')
            .finish(&mut df)
            .unwrap();
    }
}

impl Clone for DNAVariantCallSet {
    fn clone(&self) -> Self {
        DNAVariantCallSet {
            variant_calls: self.variant_calls.clone(),
            position_index: self.position_index.clone(),
            read_names_map: self.read_names_map.clone(),
            chromosome_names_map: self.chromosome_names_map.clone()
        }
    }
}
