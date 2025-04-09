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
use flate2::write::GzEncoder;
use flate2::Compression;
use polars::prelude::*;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use serde::{Serialize, Deserialize};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Mutex;

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct VariantCallSet {
    pub variant_calls: HashSet<VariantCall>,

    // left     =   read name
    // right    =   read ID
    pub read_names_map: BiMap<Box<str>,usize>,

    // left     =   chromosome name
    // right    =   chromosome ID
    pub chromosome_names_map: BiMap<Box<str>,u16>
}

impl VariantCallSet {
    pub fn new() -> Self {
        Self {
            variant_calls: HashSet::new(),
            read_names_map: BiMap::new(),
            chromosome_names_map: BiMap::new()
        }
    }

    pub fn add_variant_call(&mut self, variant_call: VariantCall) {
        self.variant_calls.insert(variant_call);
    }

    pub fn get_size(&self) -> usize {
        self.variant_calls.len()
    }

    pub fn get_variant_calls(&self) -> &HashSet<VariantCall> {
        &self.variant_calls
    }

    pub fn get_variant_records(&self) -> Vec<&VariantRecord> {
        let mut variant_records: Vec<&VariantRecord> = Vec::new();
        for variant_call in self.variant_calls.iter() {
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
        self.variant_calls.remove(variant_call);
    }

    pub fn to_dataframe(self, num_threads: usize) -> DataFrame {
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

        let variant_calls: Vec<VariantCall> = self.variant_calls.into_iter().collect();
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
                            chromosome_1,
                            consensus_record.sequence_operation.position_1,
                            consensus_record.sequence_operation.strand_1.as_str(),
                            consensus_record.sequence_operation.operation_1.as_str(),
                            chromosome_2,
                            consensus_record.sequence_operation.position_2,
                            consensus_record.sequence_operation.strand_2.as_str(),
                            consensus_record.sequence_operation.operation_2.as_str(),
                            consensus_record.get_variant_size() as i64,
                            consensus_record.get_variant_type().as_str().to_string(),
                            &*consensus_record.sequence_operation.sequence,
                            consensus_read_names.join(","),
                            consensus_read_names.len() as u32,
                            read_names.join(","),
                            read_names.len() as u32,
                        )
                    })
                })
                .collect()
        });

        let mut variant_call_ids: Vec<i64> = Vec::new();
        for (i,variant_call) in variant_calls.iter().enumerate() {
            variant_call_ids.push((i+1) as i64);
        }

        DataFrame::new(vec![
            Column::from(Series::new("variant_id".into(), variant_call_ids)),
            Column::from(Series::new("chromosome_1".into(), rows.iter().map(|r| r.0).collect::<Vec<_>>())),
            Column::from(Series::new("position_1".into(), rows.iter().map(|r| r.1).collect::<Vec<_>>())),
            Column::from(Series::new("strand_1".into(), rows.iter().map(|r| r.2).collect::<Vec<_>>())),
            Column::from(Series::new("operation_1".into(), rows.iter().map(|r| r.3).collect::<Vec<_>>())),
            Column::from(Series::new("chromosome_2".into(), rows.iter().map(|r| r.4).collect::<Vec<_>>())),
            Column::from(Series::new("position_2".into(), rows.iter().map(|r| r.5).collect::<Vec<_>>())),
            Column::from(Series::new("strand_2".into(), rows.iter().map(|r| r.6).collect::<Vec<_>>())),
            Column::from(Series::new("operation_2".into(), rows.iter().map(|r| r.7).collect::<Vec<_>>())),
            Column::from(Series::new("variant_size".into(), rows.iter().map(|r| r.8).collect::<Vec<_>>())),
            Column::from(Series::new("variant_type".into(), rows.iter().map(|r| r.9.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("variant_sequence".into(), rows.iter().map(|r| r.10).collect::<Vec<_>>())),
            Column::from(Series::new("consensus_read_names".into(), rows.iter().map(|r| r.11.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("consensus_read_names_count".into(), rows.iter().map(|r| r.12).collect::<Vec<_>>())),
            Column::from(Series::new("read_names".into(), rows.iter().map(|r| r.13.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("read_names_count".into(), rows.iter().map(|r| r.14).collect::<Vec<_>>()))
        ])
        .unwrap()
    }

    pub fn to_tsv(&self, file: &str, buffer_size: usize, num_threads: usize, gzip: bool) {
        let file = File::create(file).unwrap();

        let thread_pool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();

        let header: String = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            "variant_id",
            "chromosome_1",
            "position_1",
            "strand_1",
            "orientation_1",
            "chromosome_2",
            "position_2",
            "strand_2",
            "orientation_2",
            "variant_size",
            "variant_type",
            "variant_sequence",
            "consensus_read_names",
            "consensus_read_names_count",
            "read_names",
            "read_names_count"
        );

        let variant_calls: Vec<&VariantCall> = self.variant_calls.iter().collect();

        if gzip {
            let buf_writer = BufWriter::new(file);
            let writer = Arc::new(Mutex::new(GzEncoder::new(buf_writer, Compression::default())));
            writer.lock().unwrap().write_all(header.as_bytes()).unwrap();
            thread_pool.scope(|s| {
                for (chunk_idx, chunk) in variant_calls.chunks(buffer_size).enumerate() {
                    let writer = writer.clone();
                    s.spawn(move |_| {
                        let mut local_buffer = String::new();
                        for (i, variant_call) in chunk.iter().enumerate() {
                            let variant_call_id = format!("{}\t", chunk_idx * buffer_size + i + 1);
                            let row = variant_call.to_tsv_string(&self.chromosome_names_map, &self.read_names_map);
                            local_buffer.push_str(&(variant_call_id + &row));
                        }
                        let mut writer_guard = writer.lock().unwrap();
                        writer_guard.write_all(local_buffer.as_bytes()).unwrap();
                    });
                }
            });
            writer.lock().unwrap().flush().unwrap();
        } else {
            let writer = Arc::new(Mutex::new(BufWriter::new(file)));
            writer.lock().unwrap().write_all(header.as_bytes()).unwrap();
            thread_pool.scope(|s| {
                for (chunk_idx, chunk) in variant_calls.chunks(buffer_size).enumerate() {
                    let writer = writer.clone();
                    s.spawn(move |_| {
                        let mut local_buffer = String::new();
                        for (i, variant_call) in chunk.iter().enumerate() {
                            let variant_call_id = format!("{}\t", chunk_idx * buffer_size + i + 1);
                            let row = variant_call.to_tsv_string(&self.chromosome_names_map, &self.read_names_map);
                            local_buffer.push_str(&(variant_call_id + &row));
                        }
                        let mut writer_guard = writer.lock().unwrap();
                        writer_guard.write_all(local_buffer.as_bytes()).unwrap();
                    });
                }
            });
            writer.lock().unwrap().flush().unwrap();
        }
    }
}

impl Clone for VariantCallSet {
    fn clone(&self) -> Self {
        VariantCallSet {
            variant_calls: self.variant_calls.clone(),
            read_names_map: self.read_names_map.clone(),
            chromosome_names_map: self.chromosome_names_map.clone()
        }
    }
}
