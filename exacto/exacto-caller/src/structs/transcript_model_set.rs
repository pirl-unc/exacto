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


use serde::{Deserialize, Serialize};
use bimap::BiMap;
use flate2::Compression;
use flate2::write::GzEncoder;
use polars::prelude::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Mutex;
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSlice;
use rayon::ThreadPoolBuilder;
use crate::prelude::VariantCall;
use crate::structs::transcript_model::TranscriptModel;


#[derive(Debug,Serialize,Deserialize)]
pub struct TranscriptModelSet {
    pub transcript_models: HashSet<TranscriptModel>,

    // left     =   read name
    // right    =   read ID
    pub read_names_map: BiMap<Box<str>,usize>,

    // left     =   chromosome name
    // right    =   chromosome ID
    pub chromosome_names_map: BiMap<Box<str>,u16>
}

impl TranscriptModelSet {
    pub fn new() -> Self {
        Self {
            transcript_models: HashSet::new(),
            read_names_map: BiMap::new(),
            chromosome_names_map: BiMap::new()
        }
    }

    pub fn add_transcript_model(&mut self, transcript_model: TranscriptModel) {
        self.transcript_models.insert(transcript_model);
    }

    pub fn load_chromosome_names(&mut self, chromosome_names_map: BiMap<Box<str>,u16>) {
        self.chromosome_names_map = chromosome_names_map;
    }

    pub fn load_read_names(&mut self, read_names_map: BiMap<Box<str>,usize>) {
        self.read_names_map = read_names_map;
    }

    pub fn to_dataframes(self, num_threads: usize) -> (DataFrame,DataFrame,DataFrame) {
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

        // DataFrame of exons
        let mut transcript_id_values: Vec<String> = Vec::new();
        let mut chromosome_values: Vec<String> = Vec::new();
        let mut start_values: Vec<u32> = Vec::new();
        let mut end_values: Vec<u32> = Vec::new();
        let mut exon_number_values: Vec<u32> = Vec::new();
        let mut strand_values: Vec<String> = Vec::new();
        let mut read_names_values: Vec<String> = Vec::new();
        let mut read_names_count_values: Vec<u32> = Vec::new();
        for transcript_model in self.transcript_models.iter() {
            let read_names: Vec<String> = transcript_model
                .read_ids
                .iter()
                .map(|read_id| self.read_names_map.get_by_right(read_id).unwrap().to_string())
                .collect();
            for exon in transcript_model.exons.iter() {
                let chromosome: Box<str> = self.chromosome_names_map.get_by_right(&exon.chromosome_id).unwrap().clone();
                transcript_id_values.push(transcript_model.transcript_id.to_string());
                chromosome_values.push(chromosome.to_string());
                start_values.push(exon.start);
                end_values.push(exon.end);
                exon_number_values.push(exon.number as u32);
                strand_values.push(exon.strand.as_str().to_string());
                read_names_values.push(read_names.join(","));
                read_names_count_values.push(read_names.len() as u32);
            }
        }
        let df_exons: DataFrame = DataFrame::new(vec![
            Column::from(Series::new("transcript_id".into(), transcript_id_values)),
            Column::from(Series::new("chromosome".into(), chromosome_values)),
            Column::from(Series::new("start".into(), start_values)),
            Column::from(Series::new("end".into(), end_values)),
            Column::from(Series::new("exon_number".into(), exon_number_values)),
            Column::from(Series::new("strand_values".into(), strand_values)),
            Column::from(Series::new("read_names_values".into(), read_names_values)),
            Column::from(Series::new("read_names_count_values".into(), read_names_count_values))
        ]).unwrap();

        // DataFrame of splice junctions
        let mut transcript_id_values: Vec<String> = Vec::new();
        let mut chromosome_values: Vec<String> = Vec::new();
        let mut start_values: Vec<u32> = Vec::new();
        let mut end_values: Vec<u32> = Vec::new();
        let mut splice_junction_number_values: Vec<u32> = Vec::new();
        let mut strand_values: Vec<String> = Vec::new();
        let mut read_names_values: Vec<String> = Vec::new();
        let mut read_names_count_values: Vec<u32> = Vec::new();
        for transcript_model in self.transcript_models.iter() {
            let read_names: Vec<String> = transcript_model
                .read_ids
                .iter()
                .map(|read_id| self.read_names_map.get_by_right(read_id).unwrap().to_string())
                .collect();
            for splice_junction in transcript_model.splice_junctions.iter() {
                let chromosome: Box<str> = self.chromosome_names_map.get_by_right(&splice_junction.chromosome_id).unwrap().clone();
                transcript_id_values.push(transcript_model.transcript_id.to_string());
                chromosome_values.push(chromosome.to_string());
                start_values.push(splice_junction.start);
                end_values.push(splice_junction.end);
                splice_junction_number_values.push(splice_junction.number as u32);
                strand_values.push(splice_junction.strand.as_str().to_string());
                read_names_values.push(read_names.join(","));
                read_names_count_values.push(read_names.len() as u32);
            }
        }
        let df_splice_junctions: DataFrame = DataFrame::new(vec![
            Column::from(Series::new("transcript_id".into(), transcript_id_values)),
            Column::from(Series::new("chromosome".into(), chromosome_values)),
            Column::from(Series::new("start".into(), start_values)),
            Column::from(Series::new("end".into(), end_values)),
            Column::from(Series::new("splice_junction_number".into(), splice_junction_number_values)),
            Column::from(Series::new("strand_values".into(), strand_values)),
            Column::from(Series::new("read_names_values".into(), read_names_values)),
            Column::from(Series::new("read_names_count_values".into(), read_names_count_values))
        ]).unwrap();

        // DataFrame of variant calls
        let mut variant_calls: Vec<&VariantCall> = Vec::new();
        for transcript_model in self.transcript_models.iter() {
            for variant_call in transcript_model.variant_calls.iter() {
                variant_calls.push(variant_call);
            }
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

        let df_variant_calls: DataFrame = DataFrame::new(vec![
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
        ]).unwrap();

        (df_exons, df_splice_junctions, df_variant_calls)
    }

    pub fn to_tsv_files(
        &self,
        exon_tsv_file: &str,
        splice_junction_tsv_file: &str,
        variant_calls_tsv_file: &str,
        gzip: bool
    ) {
        // Step 1. Exon TSV file
        let exon_file = File::create(exon_tsv_file).unwrap();
        let header: String = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            "transcript_id",
            "chromosome",
            "start",
            "end",
            "exon_number",
            "strand_values",
            "read_names_values",
            "read_names_count_values"
        );
        if gzip {
            let buf_writer = BufWriter::new(exon_file);
            let mut writer = GzEncoder::new(buf_writer, Compression::default());
            writer.write_all(header.as_bytes()).unwrap();
            for transcript_model in self.transcript_models.iter() {
                for exon in transcript_model.exons.iter() {
                    let chromosome: &str = self.chromosome_names_map.get_by_right(&exon.chromosome_id).unwrap();
                    let read_names: Vec<String> = transcript_model
                        .read_ids
                        .iter()
                        .map(|read_id| self.read_names_map.get_by_right(read_id).unwrap().to_string())
                        .collect();
                    let row = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        transcript_model.transcript_id,
                        chromosome,
                        exon.start,
                        exon.end,
                        exon.number,
                        exon.strand.as_str(),
                        read_names.join(","),
                        read_names.len() as u32
                    );
                    writer.write_all((&row).as_bytes()).unwrap();
                }
            }
            writer.flush().unwrap();
        } else {
            let mut writer = BufWriter::new(exon_file);
            writer.write_all(header.as_bytes()).unwrap();
            for transcript_model in self.transcript_models.iter() {
                for exon in transcript_model.exons.iter() {
                    let chromosome: &str = self.chromosome_names_map.get_by_right(&exon.chromosome_id).unwrap();
                    let read_names: Vec<String> = transcript_model
                        .read_ids
                        .iter()
                        .map(|read_id| self.read_names_map.get_by_right(read_id).unwrap().to_string())
                        .collect();
                    let row = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        transcript_model.transcript_id,
                        chromosome,
                        exon.start,
                        exon.end,
                        exon.number,
                        exon.strand.as_str(),
                        read_names.join(","),
                        read_names.len() as u32
                    );
                    writer.write_all((&row).as_bytes()).unwrap();
                }
            }
            writer.flush().unwrap();
        }

        // Step 2. Splice junction TSV file
        let sj_file = File::create(splice_junction_tsv_file).unwrap();
        let header: String = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            "transcript_id",
            "chromosome",
            "start",
            "end",
            "splice_junction_number",
            "strand_values",
            "read_names_values",
            "read_names_count_values"
        );
        if gzip {
            let buf_writer = BufWriter::new(sj_file);
            let mut writer = GzEncoder::new(buf_writer, Compression::default());
            writer.write_all(header.as_bytes()).unwrap();
            for transcript_model in self.transcript_models.iter() {
                for splice_junction in transcript_model.splice_junctions.iter() {
                    let chromosome: &str = self.chromosome_names_map.get_by_right(&splice_junction.chromosome_id).unwrap();
                    let read_names: Vec<String> = transcript_model
                        .read_ids
                        .iter()
                        .map(|read_id| self.read_names_map.get_by_right(read_id).unwrap().to_string())
                        .collect();
                    let row = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        transcript_model.transcript_id,
                        chromosome,
                        splice_junction.start,
                        splice_junction.end,
                        splice_junction.number,
                        splice_junction.strand.as_str(),
                        read_names.join(","),
                        read_names.len() as u32
                    );
                    writer.write_all((&row).as_bytes()).unwrap();
                }
            }
            writer.flush().unwrap();
        } else {
            let mut writer = BufWriter::new(sj_file);
            writer.write_all(header.as_bytes()).unwrap();
            for transcript_model in self.transcript_models.iter() {
                for splice_junction in transcript_model.splice_junctions.iter() {
                    let chromosome: &str = self.chromosome_names_map.get_by_right(&splice_junction.chromosome_id).unwrap();
                    let read_names: Vec<String> = transcript_model
                        .read_ids
                        .iter()
                        .map(|read_id| self.read_names_map.get_by_right(read_id).unwrap().to_string())
                        .collect();
                    let row = format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                        transcript_model.transcript_id,
                        chromosome,
                        splice_junction.start,
                        splice_junction.end,
                        splice_junction.number,
                        splice_junction.strand.as_str(),
                        read_names.join(","),
                        read_names.len() as u32
                    );
                    writer.write_all((&row).as_bytes()).unwrap();
                }
            }
            writer.flush().unwrap();
        }

        // Step 3. Variant calls TSV file
        let variants_file = File::create(variant_calls_tsv_file).unwrap();
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
        let mut variant_calls: Vec<&VariantCall> = Vec::new();
        for transcript_model in self.transcript_models.iter() {
            variant_calls.extend(transcript_model.variant_calls.iter());
        }
        if gzip {
            let buf_writer = BufWriter::new(variants_file);
            let mut writer = GzEncoder::new(buf_writer, Compression::default());
            writer.write_all(header.as_bytes()).unwrap();
            for (i,variant_call) in variant_calls.iter().enumerate() {
                let variant_call_id = format!("{}\t", i + 1);
                let row = variant_call.to_tsv_string(
                    &self.chromosome_names_map,
                    &self.read_names_map
                );
                writer.write_all((variant_call_id + &row).as_bytes()).unwrap();
            }
            writer.flush().unwrap();
        } else {
            let mut writer = BufWriter::new(variants_file);
            writer.write_all(header.as_bytes()).unwrap();
            for (i,variant_call) in variant_calls.iter().enumerate() {
                let variant_call_id = format!("{}\t", i + 1);
                let row = variant_call.to_tsv_string(
                    &self.chromosome_names_map,
                    &self.read_names_map
                );
                writer.write_all((variant_call_id + &row).as_bytes()).unwrap();
            }
            writer.flush().unwrap();
        }
    }
}
