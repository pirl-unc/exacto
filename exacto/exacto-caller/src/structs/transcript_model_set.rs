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
use rayon::iter::ParallelIterator;
use rayon::prelude::ParallelSlice;
use rayon::ThreadPoolBuilder;
use serde::{Deserialize, Serialize};
use std::collections::{HashSet};
use std::fs::File;

use crate::prelude::*;


#[derive(Debug,Serialize,Deserialize)]
pub struct TranscriptModelSet {
    pub transcript_models: HashSet<TranscriptModel>,

    /// A map between read names and read IDs.
    pub read_names_map: BiMap<Box<str>, usize>,

    /// A map between chromosome names and chromosome IDs.
    pub chromosome_names_map: BiMap<Box<str>, u16>
}

/// API methods
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
            let read_name: String = self.read_names_map
                .get_by_right(&transcript_model.get_read_id())
                .unwrap()
                .to_string();
            for exon in transcript_model.get_exons().iter() {
                let chromosome: Box<str> = self.chromosome_names_map.get_by_right(&exon.reference_chromosome_id).unwrap().clone();
                transcript_id_values.push(transcript_model.get_id().to_string());
                chromosome_values.push(chromosome.to_string());
                start_values.push(exon.reference_start);
                end_values.push(exon.reference_end);
                exon_number_values.push(exon.exon_number as u32);
                strand_values.push(exon.reference_strand.as_str().to_string());
                read_names_values.push(read_name.clone());
                read_names_count_values.push(1);
            }
        }

        DataFrame::new(vec![
            Column::from(Series::new("transcript_id".into(), transcript_id_values)),
            Column::from(Series::new("chromosome".into(), chromosome_values)),
            Column::from(Series::new("start".into(), start_values)),
            Column::from(Series::new("end".into(), end_values)),
            Column::from(Series::new("exon_number".into(), exon_number_values)),
            Column::from(Series::new("strand_values".into(), strand_values)),
            Column::from(Series::new("read_names_values".into(), read_names_values)),
            Column::from(Series::new("read_names_count_values".into(), read_names_count_values))
        ]).unwrap()
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
            included_read_ids.insert(transcript_model.get_read_id());
        }
        included_read_ids
    }

    pub fn get_reference_transcript_matches_dataframe(&self) -> DataFrame {
        let mut transcript_id_values: Vec<String> = Vec::new();
        let mut reference_transcript_id_values: Vec<String> = Vec::new();
        let mut num_overlap_bases_values: Vec<u32> = Vec::new();
        let mut num_transcript_only_bases_values: Vec<u32> = Vec::new();
        let mut num_reference_transcript_only_bases_values: Vec<u32> = Vec::new();
        let mut score_values: Vec<f32> = Vec::new();
        let mut scoring_method_values: Vec<String> = Vec::new();
        for transcript_model in self.transcript_models.iter() {
            if transcript_model.get_reference_transcript_matches_map().is_empty() {
                transcript_id_values.push(transcript_model.get_id().to_string());
                reference_transcript_id_values.push("".to_string());
                num_overlap_bases_values.push(0);
                num_transcript_only_bases_values.push(0);
                num_reference_transcript_only_bases_values.push(0);
                score_values.push(0.0);
                scoring_method_values.push("".to_string());
            } else {
                for (reference_transcript_id, reference_transcript_match) in transcript_model.get_reference_transcript_matches_map().iter() {
                    transcript_id_values.push(transcript_model.get_id().to_string());
                    reference_transcript_id_values.push(reference_transcript_match.reference_transcript_id.to_string());
                    num_overlap_bases_values.push(reference_transcript_match.num_overlap_bases);
                    num_transcript_only_bases_values.push(reference_transcript_match.num_transcript_only_bases);
                    num_reference_transcript_only_bases_values.push(reference_transcript_match.num_reference_only_bases);
                    score_values.push(reference_transcript_match.score);
                    scoring_method_values.push(reference_transcript_match.scoring_method.as_str().to_string());
                }
            }
        }

        DataFrame::new(vec![
            Column::from(Series::new("transcript_id".into(), transcript_id_values)),
            Column::from(Series::new("reference_transcript_id".into(), reference_transcript_id_values)),
            Column::from(Series::new("num_overlap_bases".into(), num_overlap_bases_values)),
            Column::from(Series::new("num_transcript_only_bases".into(), num_transcript_only_bases_values)),
            Column::from(Series::new("num_reference_transcript_only_bases".into(), num_reference_transcript_only_bases_values)),
            Column::from(Series::new("score".into(), score_values)),
            Column::from(Series::new("scoring_method".into(), scoring_method_values))
        ]).unwrap()
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

        DataFrame::new(vec![
            Column::from(Series::new("read_name".into(), read_name_values)),
            Column::from(Series::new("excluded".into(), excluded_values))
        ]).unwrap()
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
            let read_name: Box<str> = self.read_names_map.get_by_right(&transcript_model.get_read_id()).unwrap().clone();
            transcript_id_values.push(transcript_model.get_id().to_string());
            read_name_values.push(read_name.to_string());
        }

        DataFrame::new(vec![
            Column::from(Series::new("transcript_id".into(), transcript_id_values)),
            Column::from(Series::new("read_name".into(), read_name_values))
        ]).unwrap()
    }

    pub fn get_size(&self) -> usize {
        self.transcript_models.len()
    }

    pub fn get_introns_dataframe(&self) -> DataFrame {
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
        let mut intron_number_values: Vec<u32> = Vec::new();
        let mut strand_values: Vec<String> = Vec::new();
        let mut read_names_values: Vec<String> = Vec::new();
        let mut read_names_count_values: Vec<u32> = Vec::new();
        for transcript_model in self.transcript_models.iter() {
            let read_name: String = self.read_names_map
                .get_by_right(&transcript_model.get_read_id())
                .unwrap()
                .to_string();
            for intron in transcript_model.get_introns().iter() {
                let chromosome: Box<str> = self.chromosome_names_map.get_by_right(&intron.reference_chromosome_id).unwrap().clone();
                transcript_id_values.push(transcript_model.get_id().to_string());
                chromosome_values.push(chromosome.to_string());
                start_values.push(intron.reference_start);
                end_values.push(intron.reference_end);
                intron_number_values.push(intron.intron_number as u32);
                strand_values.push(intron.reference_strand.as_str().to_string());
                read_names_values.push(read_name.clone());
                read_names_count_values.push(1);
            }
        }

        DataFrame::new(vec![
            Column::from(Series::new("transcript_id".into(), transcript_id_values)),
            Column::from(Series::new("chromosome".into(), chromosome_values)),
            Column::from(Series::new("start".into(), start_values)),
            Column::from(Series::new("end".into(), end_values)),
            Column::from(Series::new("intron_number".into(), intron_number_values)),
            Column::from(Series::new("strand_values".into(), strand_values)),
            Column::from(Series::new("read_names_values".into(), read_names_values)),
            Column::from(Series::new("read_names_count_values".into(), read_names_count_values))
        ]).unwrap()
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
        let mut num_intron_values: Vec<String> = Vec::new();

        for transcript_model in self.transcript_models.iter() {
            let (start_chromosome_id,start) = transcript_model.get_reference_start_position();
            let (end_chromosome_id,end) = transcript_model.get_reference_end_position();
            let start_chromosome_name = self.chromosome_names_map.get_by_right(&start_chromosome_id).unwrap().clone();
            let end_chromosome_name = self.chromosome_names_map.get_by_right(&end_chromosome_id).unwrap().clone();
            let num_exons: usize = transcript_model.get_exons().len();
            let num_introns: usize = transcript_model.get_introns().len();
            transcript_id_values.push(transcript_model.get_id().to_string());
            start_chromosome_values.push(start_chromosome_name.to_string());
            start_values.push(start);
            end_chromosome_values.push(end_chromosome_name.to_string());
            end_values.push(end);
            num_exons_values.push(num_exons.to_string());
            num_intron_values.push(num_introns.to_string());
        }

        DataFrame::new(vec![
            Column::from(Series::new("transcript_id".into(), transcript_id_values)),
            Column::from(Series::new("start_chromosome".into(), start_chromosome_values)),
            Column::from(Series::new("start".into(), start_values)),
            Column::from(Series::new("end_chromosome".into(), end_chromosome_values)),
            Column::from(Series::new("end".into(), end_values)),
            Column::from(Series::new("num_exons".into(), num_exons_values)),
            Column::from(Series::new("num_introns".into(), num_intron_values))
        ]).unwrap()
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

        // Transcript ID, reference transcript ID, VariantCall
        let mut variant_calls_list: Vec<(usize, Vec<Box<str>>, VariantCall)> = Vec::new();

        // Get variant calls
        let mut variant_call_id: usize = 1;
        for transcript_model in self.transcript_models.iter() {
            for (reference_transcript_ids, variant_records) in transcript_model.get_variant_records_map().iter() {
                for variant_record in variant_records.iter() {
                    let mut variant_call: VariantCall = VariantCall::new(variant_call_id);
                    variant_call.add_variant_record(variant_record.clone());
                    variant_call_id += 1;
                    variant_calls_list.push(
                        (
                            transcript_model.get_id(),
                            reference_transcript_ids.clone(),
                            variant_call
                        )
                    );
                }
            }
        }

        let chunk_size = (variant_calls_list.len() + num_threads - 1) / num_threads;
        let rows: Vec<_> = thread_pool.install(|| {
            variant_calls_list
                .par_chunks(chunk_size)
                .flat_map_iter(|chunk| {
                    chunk.iter().map(|(transcript_id, reference_transcript_ids, variant_call)| {
                        let reference_transcript_ids_string: String = reference_transcript_ids.iter()
                            .map(|s| s.as_ref()) // &Box<str> -> &str
                            .collect::<Vec<&str>>()
                            .join(",");
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
                            variant_call.id.to_string(),
                            *transcript_id as u32,
                            reference_transcript_ids_string,
                            chromosome_1,
                            consensus_record.get_graph_operation().get_position_1() as u32,
                            consensus_record.get_graph_operation().get_strand_1().as_str(),
                            consensus_record.get_graph_operation().get_operation_type_1().as_str(),
                            chromosome_2,
                            consensus_record.get_graph_operation().get_position_2(),
                            consensus_record.get_graph_operation().get_strand_2().as_str(),
                            consensus_record.get_graph_operation().get_operation_type_2().as_str(),
                            consensus_record.get_variant_size() as i64,
                            consensus_record.get_variant_type().as_str().to_string(),
                            consensus_record.get_graph_operation().get_standardized_sequence(),
                            consensus_read_names.join(","),
                            consensus_read_names.len() as u32,
                            read_names.join(","),
                            read_names.len() as u32,
                        )
                    })
                })
                .collect()
        });

        DataFrame::new(vec![
            Column::from(Series::new("variant_call_id".into(), rows.iter().map(|r| r.0.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("transcript_id".into(), rows.iter().map(|r| r.1).collect::<Vec<_>>())),
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
            Column::from(Series::new("variant_sequence".into(), rows.iter().map(|r| r.13.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("consensus_read_names".into(), rows.iter().map(|r| r.14.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("consensus_read_names_count".into(), rows.iter().map(|r| r.15).collect::<Vec<_>>())),
            Column::from(Series::new("read_names".into(), rows.iter().map(|r| r.16.clone()).collect::<Vec<_>>())),
            Column::from(Series::new("read_names_count".into(), rows.iter().map(|r| r.17).collect::<Vec<_>>()))
        ]).unwrap()
    }

    pub fn get_variant_call_set(&self) -> RNAVariantCallSet {
        assert!(
            !self.chromosome_names_map.is_empty(),
            "self.chromosome_names_map is empty."
        );

        assert!(
            !self.read_names_map.is_empty(),
            "self.read_names_map is empty."
        );

        let mut rna_variant_call_set: RNAVariantCallSet = RNAVariantCallSet::new();
        let mut variant_call_id: usize = 1;
        for transcript_model in self.transcript_models.iter() {
            for (reference_transcript_ids, variant_records) in transcript_model.get_variant_records_map().iter() {
                for variant_record in variant_records.iter() {
                    let mut variant_call: VariantCall = VariantCall::new(variant_call_id);
                    variant_call.add_variant_record(variant_record.clone());
                    variant_call_id += 1;
                    rna_variant_call_set.add_variant_call(
                        transcript_model.get_id(),
                        reference_transcript_ids.clone(),
                        variant_call
                    );
                }
            }
        }

        rna_variant_call_set.load_chromosome_names(self.chromosome_names_map.clone());
        rna_variant_call_set.load_read_names(self.read_names_map.clone());

        rna_variant_call_set
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
        let introns_tsv_file: String = make_path("introns");
        let transcripts_tsv_file: String = make_path("transcripts");
        let variant_calls_tsv_file: String = make_path("rna_variant_calls");

        // Step 2. Get all DataFrames to output
        let mut df_exons = self.get_exons_dataframe();
        let mut df_read_filter_status: DataFrame = self.get_read_filter_status_dataframe();
        let mut df_read_names: DataFrame = self.get_read_names_dataframe();
        let mut df_matched_reference_transcripts: DataFrame = self.get_reference_transcript_matches_dataframe();
        let mut df_introns: DataFrame = self.get_introns_dataframe();
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
        write_df_to_tsv(&mut df_introns, introns_tsv_file.as_str());
        write_df_to_tsv(&mut df_transcripts, transcripts_tsv_file.as_str());
        write_df_to_tsv(&mut df_variant_calls, variant_calls_tsv_file.as_str());
    }
}
