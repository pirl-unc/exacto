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
use polars::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSlice;
use rayon::ThreadPoolBuilder;
use exacto_util::prelude::Transcript;
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

    pub fn get_exons_dataframe(&self) -> DataFrame {
        assert!(
            !self.chromosome_names_map.is_empty(),
            "self.chromosome_names_map is empty."
        );
        assert!(
            !self.read_names_map.is_empty(),
            "self.read_names_map is empty."
        );

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

        df_exons
    }

    fn get_excluded_read_ids(&self) -> HashSet<usize> {
        let read_ids: Vec<usize> = self.read_names_map.right_values().cloned().collect();
        let included_read_ids: HashSet<usize> = self.get_included_read_ids();
        let mut excluded_read_ids: HashSet<usize> = HashSet::new();
        for read_id in read_ids.iter() {
            if included_read_ids.contains(read_id) == false {
                excluded_read_ids.insert(*read_id);
            }
        }
        excluded_read_ids
    }

    fn get_included_read_ids(&self) -> HashSet<usize> {
        let mut included_read_ids: HashSet<usize> = HashSet::new();
        for transcript_model in self.transcript_models.iter() {
            for read_id in transcript_model.read_ids.iter() {
                included_read_ids.insert(*read_id);
            }
        }
        included_read_ids
    }

    pub fn get_matched_reference_transcripts_dataframe(&self) -> DataFrame {
        let mut transcript_id_values: Vec<String> = Vec::new();
        let mut reference_gene_id_values: Vec<String> = Vec::new();
        let mut reference_transcript_id_values: Vec<String> = Vec::new();
        let mut reference_transcript_chromosome_values: Vec<String> = Vec::new();
        let mut reference_transcript_start_values: Vec<u32> = Vec::new();
        let mut reference_transcript_end_values: Vec<u32> = Vec::new();
        let mut reference_transcript_strand_values: Vec<String> = Vec::new();
        let mut reference_transcript_num_exons_values: Vec<u32> = Vec::new();
        let mut num_overlap_bases_values: Vec<u32> = Vec::new();
        let mut num_transcript_only_bases_values: Vec<u32> = Vec::new();
        let mut num_reference_transcript_only_bases_values: Vec<u32> = Vec::new();
        let mut score_values: Vec<f32> = Vec::new();
        let mut scoring_method_values: Vec<String> = Vec::new();
        for transcript_model in self.transcript_models.iter() {
            if transcript_model.reference_transcript_matches.is_empty() {
                transcript_id_values.push(transcript_model.transcript_id.to_string());
                reference_gene_id_values.push("".to_string());
                reference_transcript_id_values.push("".to_string());
                reference_transcript_chromosome_values.push("".to_string());
                reference_transcript_start_values.push(0);
                reference_transcript_end_values.push(0);
                reference_transcript_strand_values.push("".to_string());
                reference_transcript_num_exons_values.push(0);
                num_overlap_bases_values.push(0);
                num_transcript_only_bases_values.push(0);
                num_reference_transcript_only_bases_values.push(0);
                score_values.push(0.0);
                scoring_method_values.push("".to_string());
            } else {
                for reference_transcript_match in transcript_model.reference_transcript_matches.iter() {
                    transcript_id_values.push(transcript_model.transcript_id.to_string());
                    reference_gene_id_values.push(reference_transcript_match.reference_transcript.gene_id.to_string());
                    reference_transcript_id_values.push(reference_transcript_match.reference_transcript.transcript_id.to_string());
                    reference_transcript_chromosome_values.push(reference_transcript_match.reference_transcript.chromosome.to_string());
                    reference_transcript_start_values.push(reference_transcript_match.reference_transcript.start);
                    reference_transcript_end_values.push(reference_transcript_match.reference_transcript.end);
                    reference_transcript_strand_values.push(reference_transcript_match.reference_transcript.strand.as_str().to_string());
                    reference_transcript_num_exons_values.push(reference_transcript_match.reference_transcript.exons.values().len() as u32);
                    num_overlap_bases_values.push(reference_transcript_match.num_overlap_bases);
                    num_transcript_only_bases_values.push(reference_transcript_match.num_transcript_only_bases);
                    num_reference_transcript_only_bases_values.push(reference_transcript_match.num_reference_only_bases);
                    score_values.push(reference_transcript_match.score);
                    scoring_method_values.push(reference_transcript_match.scoring_method.as_str().to_string());
                }
            }
        }

        let df_reference_transcript_matches: DataFrame = DataFrame::new(vec![
            Column::from(Series::new("transcript_id".into(), transcript_id_values)),
            Column::from(Series::new("reference_gene_id".into(), reference_gene_id_values)),
            Column::from(Series::new("reference_transcript_id".into(), reference_transcript_id_values)),
            Column::from(Series::new("reference_transcript_chromosome".into(), reference_transcript_chromosome_values)),
            Column::from(Series::new("reference_transcript_start".into(), reference_transcript_start_values)),
            Column::from(Series::new("reference_transcript_end".into(), reference_transcript_end_values)),
            Column::from(Series::new("reference_transcript_strand".into(), reference_transcript_strand_values)),
            Column::from(Series::new("reference_transcript_num_exons".into(), reference_transcript_num_exons_values)),
            Column::from(Series::new("num_overlap_bases".into(), num_overlap_bases_values)),
            Column::from(Series::new("num_transcript_only_bases".into(), num_transcript_only_bases_values)),
            Column::from(Series::new("num_reference_transcript_only_bases".into(), num_reference_transcript_only_bases_values)),
            Column::from(Series::new("score".into(), score_values)),
            Column::from(Series::new("scoring_method".into(), scoring_method_values))
        ]).unwrap();

        df_reference_transcript_matches
    }

    pub fn get_read_filter_status_dataframe(&self) -> DataFrame {
        let included_read_ids: HashSet<usize> = self.get_included_read_ids();
        let excluded_read_ids: HashSet<usize> = self.get_excluded_read_ids();

        let mut read_name_values: Vec<String> = Vec::new();
        let mut excluded_values: Vec<bool> = Vec::new();
        for read_id in included_read_ids.iter() {
            let read_name: Box<str> = self.read_names_map.get_by_right(read_id).unwrap().clone();
            read_name_values.push(read_name.to_string());
            excluded_values.push(false);
        }
        for read_id in excluded_read_ids.iter() {
            let read_name: Box<str> = self.read_names_map.get_by_right(read_id).unwrap().clone();
            read_name_values.push(read_name.to_string());
            excluded_values.push(true);
        }

        let df_read_filter_status: DataFrame = DataFrame::new(vec![
            Column::from(Series::new("read_name".into(), read_name_values)),
            Column::from(Series::new("excluded".into(), excluded_values))
        ]).unwrap();

        df_read_filter_status
    }

    pub fn get_read_names_dataframe(&self) -> DataFrame {
        assert!(
            !self.read_names_map.is_empty(),
            "self.read_names_map is empty."
        );

        // DataFrame of read names
        let mut transcript_id_values: Vec<String> = Vec::new();
        let mut read_name_values: Vec<String> = Vec::new();
        for transcript_model in self.transcript_models.iter() {
            for read_id in transcript_model.read_ids.iter() {
                let read_name: Box<str> = self.read_names_map.get_by_right(read_id).unwrap().clone();
                transcript_id_values.push(transcript_model.transcript_id.to_string());
                read_name_values.push(read_name.to_string());
            }
        }

        let df_read_names: DataFrame = DataFrame::new(vec![
            Column::from(Series::new("transcript_id".into(), transcript_id_values)),
            Column::from(Series::new("read_name".into(), read_name_values))
        ]).unwrap();

        df_read_names
    }

    pub fn get_splice_junctions_dataframe(&self) -> DataFrame {
        assert!(
            !self.chromosome_names_map.is_empty(),
            "self.chromosome_names_map is empty."
        );
        assert!(
            !self.read_names_map.is_empty(),
            "self.read_names_map is empty."
        );

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

        df_splice_junctions
    }

    pub fn get_transcripts_dataframe(&self) -> DataFrame {
        assert!(
            !self.chromosome_names_map.is_empty(),
            "self.chromosome_names_map is empty."
        );
        assert!(
            !self.read_names_map.is_empty(),
            "self.read_names_map is empty."
        );

        let mut transcript_id_values: Vec<String> = Vec::new();
        let mut start_chromosome_values: Vec<String> = Vec::new();
        let mut start_values: Vec<u32> = Vec::new();
        let mut end_chromosome_values: Vec<String> = Vec::new();
        let mut end_values: Vec<u32> = Vec::new();
        let mut num_exons_values: Vec<String> = Vec::new();
        let mut num_splice_junction_values: Vec<String> = Vec::new();

        for transcript_model in self.transcript_models.iter() {
            let (start_chromosome_id,start) = transcript_model.get_start_position();
            let (end_chromosome_id,end) = transcript_model.get_end_position();
            let start_chromosome_name = self.chromosome_names_map.get_by_right(&start_chromosome_id).unwrap().clone();
            let end_chromosome_name = self.chromosome_names_map.get_by_right(&end_chromosome_id).unwrap().clone();
            let num_exons: usize = transcript_model.exons.len();
            let num_splice_junctions: usize = transcript_model.splice_junctions.len();
            transcript_id_values.push(transcript_model.transcript_id.to_string());
            start_chromosome_values.push(start_chromosome_name.to_string());
            start_values.push(start);
            end_chromosome_values.push(end_chromosome_name.to_string());
            end_values.push(end);
            num_exons_values.push(num_exons.to_string());
            num_splice_junction_values.push(num_splice_junctions.to_string());
        }

        let df_transcripts: DataFrame = DataFrame::new(vec![
            Column::from(Series::new("transcript_id".into(), transcript_id_values)),
            Column::from(Series::new("start_chromosome".into(), start_chromosome_values)),
            Column::from(Series::new("start".into(), start_values)),
            Column::from(Series::new("end_chromosome".into(), end_chromosome_values)),
            Column::from(Series::new("end".into(), end_values)),
            Column::from(Series::new("num_exons".into(), num_exons_values)),
            Column::from(Series::new("num_splice_junctions".into(), num_splice_junction_values))
        ]).unwrap();

        df_transcripts
    }

    pub fn get_variant_calls_dataframe(&self, num_threads: usize) -> DataFrame {
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
    
        // Transcript ID, reference gene ID, reference transcript ID, VariantCall
        let mut variant_calls: Vec<(usize, Box<str>, Box<str>, &VariantCall)> = Vec::new();

        // Get sequence and splice variant calls
        for transcript_model in self.transcript_models.iter() {
            // Sequence variant calls
            variant_calls.extend(
                transcript_model
                    .sequence_variant_calls
                    .iter()
                    .map(|variant_call| (transcript_model.transcript_id, "".into(), "".into(), variant_call)),
            );
            
            // Splice variant calls
            let reference_transcripts_map: HashMap<Box<str>, &Transcript> = transcript_model
                .reference_transcript_matches
                .iter()
                .map(|m| (m.reference_transcript.transcript_id.clone(), &m.reference_transcript))
                .collect();
            for (reference_transcript_id, variant_call_) in &transcript_model.splice_variant_calls {
                if let Some(reference_transcript) = reference_transcripts_map.get(reference_transcript_id) {
                    variant_calls.extend(variant_call_.iter().map(|variant_call| (transcript_model.transcript_id, reference_transcript.gene_id.clone(), reference_transcript.transcript_id.clone(), variant_call)));
                } else {
                    // Fusion gene
                    variant_calls.extend(variant_call_.iter().map(|variant_call| (transcript_model.transcript_id, reference_transcript_id.clone(), reference_transcript_id.clone(), variant_call)));
                }
            }
        }
        
        let chunk_size = (variant_calls.len() + num_threads - 1) / num_threads;
        let rows: Vec<_> = thread_pool.install(|| {
            variant_calls
                .par_chunks(chunk_size)
                .flat_map_iter(|chunk| {
                    chunk.iter().map(|(transcript_id, reference_gene_id, reference_transcript_id, variant_call)| {
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
                            *transcript_id as u32,
                            reference_gene_id.to_string(),
                            reference_transcript_id.to_string(),
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
            Column::from(Series::new("transcript_id".into(), rows.iter().map(|r| r.0).collect::<Vec<_>>())),
            Column::from(Series::new("reference_gene_id".into(), rows.iter().map(|r| r.1.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("reference_transcript_id".into(), rows.iter().map(|r| r.2.clone()).collect::<Vec<_>>())),
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
            Column::from(Series::new("variant_sequence".into(), rows.iter().map(|r| r.13).collect::<Vec<_>>())),
            Column::from(Series::new("consensus_read_names".into(), rows.iter().map(|r| r.14.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("consensus_read_names_count".into(), rows.iter().map(|r| r.15).collect::<Vec<_>>())),
            Column::from(Series::new("read_names".into(), rows.iter().map(|r| r.16.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("read_names_count".into(), rows.iter().map(|r| r.17).collect::<Vec<_>>()))
        ]).unwrap();

        df_variant_calls
    }

    pub fn load_chromosome_names(&mut self, chromosome_names_map: BiMap<Box<str>,u16>) {
        self.chromosome_names_map = chromosome_names_map;
    }

    pub fn load_read_names(&mut self, read_names_map: BiMap<Box<str>,usize>) {
        self.read_names_map = read_names_map;
    }

    pub fn to_tsv_files(
        &self,
        output_dir: &str,
        prefix: &str
    ) {
        // Step 1. Define output TSV file paths
        let output_dir = if output_dir.ends_with('/') {
            output_dir.to_string()
        } else {
            format!("{}/", output_dir)
        };

        let prefix = if prefix.is_empty() {
            "sample".to_string()
        } else {
            prefix.to_string()
        };

        let make_path = |name: &str| format!("{}{}_exacto_{}.tsv", output_dir, prefix, name);

        let exons_tsv_file: String = make_path("exons");
        let matched_reference_transcripts_tsv_file: String = make_path("reference_transcript_matches");
        let read_filter_status_tsv_file: String = make_path("read_filter_status");
        let read_names_tsv_file: String = make_path("transcripts_read_support");
        let splice_junctions_tsv_file: String = make_path("splice_junctions");
        let transcripts_tsv_file: String = make_path("transcripts");
        let variant_calls_tsv_file: String = make_path("variant_calls");

        // Step 2. Get all DataFrames to output
        let mut df_exons = self.get_exons_dataframe();
        let mut df_read_filter_status: DataFrame = self.get_read_filter_status_dataframe();
        let mut df_read_names: DataFrame = self.get_read_names_dataframe();
        let mut df_matched_reference_transcripts: DataFrame = self.get_matched_reference_transcripts_dataframe();
        let mut df_splice_junctions: DataFrame = self.get_splice_junctions_dataframe();
        let mut df_transcripts: DataFrame = self.get_transcripts_dataframe();
        let mut df_variant_calls: DataFrame = self.get_variant_calls_dataframe(1);

        // Step 3. Write to TSV files
        fn write_df_to_tsv(df: &mut DataFrame, path: &str) {
            let mut file = File::create(path).unwrap();
            CsvWriter::new(&mut file)
                .include_header(true)
                .with_separator(b'\t')
                .finish(df)
                .unwrap();
        }
        write_df_to_tsv(&mut df_exons, exons_tsv_file.as_str());
        write_df_to_tsv(&mut df_read_filter_status, read_filter_status_tsv_file.as_str());
        write_df_to_tsv(&mut df_read_names, read_names_tsv_file.as_str());
        write_df_to_tsv(&mut df_matched_reference_transcripts, matched_reference_transcripts_tsv_file.as_str());
        write_df_to_tsv(&mut df_splice_junctions, splice_junctions_tsv_file.as_str());
        write_df_to_tsv(&mut df_transcripts, transcripts_tsv_file.as_str());
        write_df_to_tsv(&mut df_variant_calls, variant_calls_tsv_file.as_str());
    }
}
