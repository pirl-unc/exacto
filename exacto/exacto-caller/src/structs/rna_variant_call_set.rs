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
use exacto_core::prelude::*;
use polars::prelude::*;
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};
use serde::{Serialize, Deserialize};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::str::FromStr;

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct RNAVariantCallSet {
    /// A map between a variant call ID and its associated `VariantCall` object.
    pub variant_calls: HashMap<usize, VariantCall>,

    /// A map between reference transcript IDs and variant call IDs.
    /// Nested structure for transcript indexing:
    /// - Outer HashMap: Maps transcript model IDs to their reference transcript IDs
    /// - Inner HashMap: Maps reference transcript IDs to variant call IDs
    pub variant_calls_index: HashMap<usize, HashMap<Vec<Box<str>>, HashSet<usize>>>,

    /// A bidirectional map between read names and read IDs.
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

impl RNAVariantCallSet {
    pub fn new() -> Self {
        Self {
            variant_calls: HashMap::new(),
            variant_calls_index: HashMap::new(),
            read_names_map: BiMap::new(),
            chromosome_names_map: BiMap::new(),
            position_index: HashMap::new()
        }
    }

    pub fn add_variant_call(&mut self, transcript_model_id: usize, reference_transcript_ids: Vec<Box<str>>, variant_call: VariantCall) {
        // Sort reference transcript IDs
        let mut reference_transcript_ids_: Vec<Box<str>> = reference_transcript_ids.clone();
        reference_transcript_ids_.sort();

        if self.variant_calls.contains_key(&variant_call.id) {
            panic!("VariantCallSet already contains variant call with id: {}", variant_call.id);
        }

        // Index the variant call
        self.variant_calls_index
            .entry(transcript_model_id)
            .or_insert_with(HashMap::new)
            .entry(reference_transcript_ids_.clone())
            .or_insert_with(HashSet::new)
            .insert(variant_call.id);

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

    pub fn get_size(&self) -> usize {
        self.variant_calls.len()
    }

    pub fn get_transcript_model_ids(&self) -> Vec<usize> {
        self.variant_calls_index.keys().cloned().collect()
    }

    pub fn get_variant_call(&self, variant_call_id: usize) -> &VariantCall {
        self.variant_calls.get(&variant_call_id).unwrap()
    }

    pub fn get_variant_calls_index(&self) -> &HashMap<usize, HashMap<Vec<Box<str>>, HashSet<usize>>> {
        &self.variant_calls_index
    }

    pub fn get_variant_calls_for_transcript_model_id(&self, transcript_model_id: usize) -> HashMap<&Vec<Box<str>>, Vec<&VariantCall>> {
        let mut variant_calls_map: HashMap<&Vec<Box<str>>, Vec<&VariantCall>> = HashMap::new();
        if self.variant_calls_index.contains_key(&transcript_model_id) {
            for (reference_transcript_ids, variant_call_ids) in self.variant_calls_index.get(&transcript_model_id).unwrap().iter() {
                for variant_call_id in variant_call_ids.iter() {
                    let variant_call: &VariantCall = self.get_variant_call(*variant_call_id);
                    variant_calls_map
                        .entry(reference_transcript_ids)
                        .or_insert_with(Vec::new)
                        .push(variant_call);
                }
            }
        }
        variant_calls_map
    }

    pub fn get_variant_calls_by_range(&self, chromosome_id: u16, start: u32, end: u32) -> Vec<&VariantCall> {
        let mut result_ids = HashSet::new();

        if let Some(position_map) = self.position_index.get(&chromosome_id) {
            for (_, variant_call_ids) in position_map.range(start..=end) {
                result_ids.extend(variant_call_ids.iter().cloned());
            }
        }

        result_ids.iter()
            .filter_map(|id| self.variant_calls.get(id))
            .collect()
    }

    pub fn load_chromosome_names(&mut self, chromosome_names_map: BiMap<Box<str>,u16>) {
        self.chromosome_names_map = chromosome_names_map;
    }

    pub fn load_read_names(&mut self, read_names_map: BiMap<Box<str>,usize>) {
        self.read_names_map = read_names_map;
    }
    
    pub fn read_tsv_file(tsv_file: &str) -> Self {
        let parse_options = CsvParseOptions::default()
            .with_separator(b'\t');
        let df: DataFrame = CsvReadOptions::default()
            .with_parse_options(parse_options)
            .with_has_header(true)
            .try_into_reader_with_file_path(Some(tsv_file.into()))
            .unwrap()
            .finish()
            .unwrap();

        let variant_call_id_col = df.column("variant_call_id").unwrap().i64().unwrap();
        let transcript_model_id_col = df.column("transcript_model_id").unwrap().i64().unwrap();
        let reference_transcript_ids_col = df.column("reference_transcript_ids").unwrap().str().unwrap();
        let chromosome_1_col = df.column("chromosome_1").unwrap().str().unwrap();
        let position_1_col = df.column("position_1").unwrap().i64().unwrap();
        let strand_1_col = df.column("strand_1").unwrap().str().unwrap();
        let operation_1_col = df.column("operation_1").unwrap().str().unwrap();
        let chromosome_2_col = df.column("chromosome_2").unwrap().str().unwrap();
        let position_2_col = df.column("position_2").unwrap().i64().unwrap();
        let strand_2_col = df.column("strand_2").unwrap().str().unwrap();
        let operation_2_col = df.column("operation_2").unwrap().str().unwrap();
        let variant_type_col = df.column("variant_type").unwrap().str().unwrap();
        let variant_sequence_col = df.column("sequence").unwrap().str().unwrap();
        let consensus_read_names_col = df
            .column("consensus_read_names")
            .unwrap()
            .cast(&DataType::String)
            .unwrap();
        let consensus_read_names_col = consensus_read_names_col.str().unwrap();
        let read_start_col = df.column("read_start").unwrap().i64().unwrap();
        let read_end_col = df.column("read_end").unwrap().i64().unwrap();

        let mut variant_call_set_rna: RNAVariantCallSet = RNAVariantCallSet::new();

        // Chromosome names
        let mut chromosomes: HashSet<Box<str>> = HashSet::new();
        for i in 0 ..df.height() {
            chromosomes.insert(chromosome_1_col.get(i).unwrap().into());
            chromosomes.insert(chromosome_2_col.get(i).unwrap().into());
        }
        let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
        for (i,chromosome) in chromosomes.iter().enumerate() {
            chromosome_names_map.insert(chromosome.clone(), i as u16);
        }
        
        // Read names
        let mut read_names: HashSet<Box<str>> = HashSet::new();
        for i in 0 ..df.height() {
            let consensus_read_names: Vec<Box<str>> = consensus_read_names_col.get(i).unwrap().split(",").map(|s| s.into()).collect();
            for consensus_read_name in consensus_read_names {
                read_names.insert(consensus_read_name);
            }
        }
        let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
        for (i,read_name) in read_names.iter().enumerate() {
            read_names_map.insert(read_name.clone(), i);
        }

        for i in 0 ..df.height() {
            let variant_call_id: usize = variant_call_id_col.get(i).unwrap() as usize;
            let transcript_model_id: usize = transcript_model_id_col.get(i).unwrap() as usize;
            let mut reference_transcript_ids: Vec<Box<str>> = reference_transcript_ids_col.get(i).unwrap().split(",").map(|s| s.into()).collect();
            reference_transcript_ids.sort();
            let chromosome_1: Box<str> = chromosome_1_col.get(i).unwrap().into();
            let position_1: u32 = position_1_col.get(i).unwrap() as u32;
            let strand_1: Strand = Strand::from_str(strand_1_col.get(i).unwrap().into()).unwrap();
            let operation_1: GraphOperationType = GraphOperationType::from_str(operation_1_col.get(i).unwrap().into()).unwrap();
            let chromosome_2: Box<str> = chromosome_2_col.get(i).unwrap().into();
            let position_2: u32 = position_2_col.get(i).unwrap() as u32;
            let strand_2: Strand = Strand::from_str(strand_2_col.get(i).unwrap().into()).unwrap();
            let operation_2: GraphOperationType = GraphOperationType::from_str(operation_2_col.get(i).unwrap().into()).unwrap();
            let variant_type: VariantType = VariantType::from_str(variant_type_col.get(i).unwrap()).unwrap();
            let variant_sequence: Box<str> = variant_sequence_col.get(i).unwrap().into();
            let consensus_read_names: Box<str> = consensus_read_names_col.get(i).unwrap().into();
            let read_start: u32 = read_start_col.get(i).unwrap() as u32;
            let read_end: u32 = read_end_col.get(i).unwrap() as u32;

            let sequence_operation: GraphOperation = GraphOperation::new(
                *chromosome_names_map.get_by_left(&*chromosome_1).unwrap(),
                position_1,
                strand_1,
                operation_1,
                *chromosome_names_map.get_by_left(&*chromosome_2).unwrap(),
                position_2,
                strand_2,
                operation_2,
                variant_sequence,
                variant_type
            );

            let mut variant_call: VariantCall = VariantCall::new(variant_call_id);

            for consensus_read_name in consensus_read_names.split(",") {
                let read_id: usize = *read_names_map.get_by_left(consensus_read_name).unwrap();
                let variant_record: VariantRecord = VariantRecord::new(
                    read_id,
                    read_start,
                    read_end,
                    sequence_operation.clone()
                );
                variant_call.add_variant_record(variant_record);
            }

            variant_call_set_rna.add_variant_call(transcript_model_id, reference_transcript_ids, variant_call);
        }

        variant_call_set_rna.load_chromosome_names(chromosome_names_map);
        variant_call_set_rna.load_read_names(read_names_map);

        variant_call_set_rna
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

        let thread_pool: ThreadPool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();
        
        let transcript_model_ids: Vec<usize> = self.variant_calls_index.keys().cloned().collect();

        let rows: Vec<_> = thread_pool.install(|| {
            transcript_model_ids
                .into_par_iter()
                .flat_map_iter(|transcript_model_id| {
                    let mut rows = Vec::new();
                    for (reference_transcript_ids, variant_call_ids) in self.variant_calls_index.get(&transcript_model_id).unwrap().iter() {
                        for variant_call_id in variant_call_ids.iter() {
                            let variant_call: &VariantCall = self.variant_calls.get(variant_call_id).unwrap();
                            let reference_transcript_ids_string: String = reference_transcript_ids.iter()
                                .map(|s| s.as_ref()) // &Box<str> -> &str
                                .collect::<Vec<&str>>()
                                .join(",");
                            assert!(variant_call.variant_records.len() == 1, "variant_call.variant_records.len() != 1");
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
                            rows.push(
                                (
                                    variant_call.id as u64,
                                    transcript_model_id.clone() as u64,
                                    reference_transcript_ids_string,
                                    chromosome_1,
                                    consensus_record.graph_operation.get_position_1() as u64,
                                    consensus_record.graph_operation.get_strand_1().as_str(),
                                    consensus_record.graph_operation.get_operation_type_1().as_str(),
                                    chromosome_2,
                                    consensus_record.graph_operation.get_position_2() as u64,
                                    consensus_record.graph_operation.get_strand_2().as_str(),
                                    consensus_record.graph_operation.get_operation_type_2().as_str(),
                                    consensus_record.get_variant_size() as i64,
                                    consensus_record.get_variant_type().as_str().to_string(),
                                    &*consensus_record.graph_operation.get_sequence(),
                                    consensus_read_names.join(","),
                                    consensus_read_names.len() as u64,
                                    read_names.join(","),
                                    read_names.len() as u64,
                                    consensus_record.read_position_1 as u64,
                                    consensus_record.read_position_2 as u64
                                )
                            );
                        }
                    }

                    rows
                })
                .collect()
        });

        if rows.is_empty() {
            return DataFrame::new(vec![
                Column::from(Series::new("variant_call_id".into(), Vec::<u64>::new())),
                Column::from(Series::new("transcript_model_id".into(), Vec::<u64>::new())),
                Column::from(Series::new("reference_transcript_ids".into(), Vec::<String>::new())),
                Column::from(Series::new("chromosome_1".into(), Vec::<&str>::new())),
                Column::from(Series::new("position_1".into(), Vec::<u64>::new())),
                Column::from(Series::new("strand_1".into(), Vec::<&str>::new())),
                Column::from(Series::new("operation_1".into(), Vec::<&str>::new())),
                Column::from(Series::new("chromosome_2".into(), Vec::<&str>::new())),
                Column::from(Series::new("position_2".into(), Vec::<u64>::new())),
                Column::from(Series::new("strand_2".into(), Vec::<&str>::new())),
                Column::from(Series::new("operation_2".into(), Vec::<&str>::new())),
                Column::from(Series::new("variant_size".into(), Vec::<i64>::new())),
                Column::from(Series::new("variant_type".into(), Vec::<String>::new())),
                Column::from(Series::new("sequence".into(), Vec::<&str>::new())),
                Column::from(Series::new("consensus_read_names".into(), Vec::<String>::new())),
                Column::from(Series::new("consensus_read_names_count".into(), Vec::<u64>::new())),
                Column::from(Series::new("read_names".into(), Vec::<String>::new())),
                Column::from(Series::new("read_names_count".into(), Vec::<u64>::new())),
                Column::from(Series::new("read_start".into(), Vec::<u64>::new())),
                Column::from(Series::new("read_end".into(), Vec::<u64>::new()))
            ]).unwrap();
        }

        let df: DataFrame = DataFrame::new(vec![
            Column::from(Series::new("variant_call_id".into(), rows.iter().map(|r| r.0).collect::<Vec<_>>())),
            Column::from(Series::new("transcript_model_id".into(), rows.iter().map(|r| r.1).collect::<Vec<_>>())),
            Column::from(Series::new("reference_transcript_ids".into(), rows.iter().map(|r| r.2.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("chromosome_1".into(), rows.iter().map(|r| r.3).collect::<Vec<_>>())),
            Column::from(Series::new("position_1".into(), rows.iter().map(|r| r.4).collect::<Vec<_>>())),
            Column::from(Series::new("strand_1".into(), rows.iter().map(|r| r.5).collect::<Vec<_>>())),
            Column::from(Series::new("operation_1".into(), rows.iter().map(|r| r.6).collect::<Vec<_>>())),
            Column::from(Series::new("chromosome_2".into(), rows.iter().map(|r| r.7).collect::<Vec<_>>())),
            Column::from(Series::new("position_2".into(), rows.iter().map(|r| r.8).collect::<Vec<_>>())),
            Column::from(Series::new("strand_2".into(), rows.iter().map(|r| r.9).collect::<Vec<_>>())),
            Column::from(Series::new("operation_2".into(), rows.iter().map(|r| r.10).collect::<Vec<_>>())),
            Column::from(Series::new("variant_size".into(), rows.iter().map(|r| r.11).collect::<Vec<_>>())),
            Column::from(Series::new("variant_type".into(), rows.iter().map(|r| r.12.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("sequence".into(), rows.iter().map(|r| r.13).collect::<Vec<_>>())),
            Column::from(Series::new("consensus_read_names".into(), rows.iter().map(|r| r.14.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("consensus_read_names_count".into(), rows.iter().map(|r| r.15).collect::<Vec<_>>())),
            Column::from(Series::new("read_names".into(), rows.iter().map(|r| r.16.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("read_names_count".into(), rows.iter().map(|r| r.17).collect::<Vec<_>>())),
            Column::from(Series::new("read_start".into(), rows.iter().map(|r| r.18.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("read_end".into(), rows.iter().map(|r| r.19).collect::<Vec<_>>()))
        ]).unwrap();

        df.sort(["variant_call_id"], SortMultipleOptions::default()).unwrap()
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

impl Clone for RNAVariantCallSet {
    fn clone(&self) -> Self {
        RNAVariantCallSet {
            variant_calls: self.variant_calls.clone(),
            variant_calls_index: self.variant_calls_index.clone(),
            read_names_map: self.read_names_map.clone(),
            chromosome_names_map: self.chromosome_names_map.clone(),
            position_index: self.position_index.clone()
        }
    }
}
