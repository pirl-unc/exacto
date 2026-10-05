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


extern crate exacto;
extern crate polars;
extern crate pyo3;
extern crate pyo3_polars;

use exacto::cluster::prelude as cluster;
use exacto::consensus::prelude as consensus;
use exacto::core::prelude as core;
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use std::collections::{HashMap,HashSet};
use std::path::{Path, PathBuf};

use crate::functions::io_error;


#[pyfunction]
pub fn determine_rna_consensus(
    py: Python,
    tsv_file: String,
    fastq_file: String,
    output_dir: String,
    output_prefix: String,
    match_score: i32,
    mismatch_score: i32,
    gap_open_score: i32,
    gap_extend_score: i32,
    max_reads_per_cluster: usize,
    seed: usize,
    orientation_kmer_size: usize,
    num_threads: usize,
    output_type: String
) -> PyResult<PyDataFrame> {
    let _ = py;
    let cluster_records: Vec<cluster::RNAReadClusterIDRecord> = cluster::load_rna_read_cluster_id_records(
        tsv_file.as_str()
    );
    let mut sequence_clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::new();
    for record in cluster_records.iter() {
        sequence_clusters
            .entry(record.cluster_id)
            .or_insert_with(HashSet::new)
            .insert(record.read_name.clone());
    }

    let options: consensus::IdentifyConsensusSequencesOptions = consensus::IdentifyConsensusSequencesOptions {
        poa_match_score: match_score,
        poa_mismatch_score: mismatch_score,
        poa_gap_open_score: gap_open_score,
        poa_gap_extend_score: gap_extend_score,
        max_reads_per_cluster: max_reads_per_cluster,
        seed: seed,
        orientation_kmer_size: orientation_kmer_size
    };
    let consensus_set: consensus::ConsensusSequenceSet = consensus::identify_consensus_sequences(
        &sequence_clusters,
        fastq_file.as_str(),
        &options,
        num_threads
    );

    match output_type.as_str() {
        "dataframe" => {
            let df_consensus: DataFrame = consensus::consensus_sequence_records_to_dataframe(
                consensus::build_consensus_sequence_records(&consensus_set)
            );
            Ok(PyDataFrame(df_consensus))
        }
        "file" => {
            let consensus_tsv_file: PathBuf = if output_prefix.is_empty() {
                Path::new(&output_dir).join("exacto_rna_consensus.tsv")
            } else {
                Path::new(&output_dir).join(format!("{}_exacto_rna_consensus.tsv", output_prefix))
            };
            let consensus_fasta_file: PathBuf = if output_prefix.is_empty() {
                Path::new(&output_dir).join("exacto_rna_consensus.fasta")
            } else {
                Path::new(&output_dir).join(format!("{}_exacto_rna_consensus.fasta", output_prefix))
            };

            core::write_tsv_file(
                consensus::build_consensus_sequence_records(&consensus_set),
                &consensus_tsv_file
            ).map_err(io_error)?;

            // (name, sequence) pairs: name = cluster ID, sequence = the consensus sequence.
            let fasta_sequences: Vec<(Box<str>, Box<str>)> = consensus_set.sequences
                .iter()
                .map(|sequence| (
                    sequence.cluster_id.to_string().into_boxed_str(),
                    sequence.consensus_sequence.clone()
                ))
                .collect();

            core::write_fasta_file(
                &fasta_sequences,
                consensus_fasta_file.to_str().unwrap()
            );

            Ok(PyDataFrame(DataFrame::new(vec![]).unwrap()))
        }
        other => {
            let error_message = format!("Unsupported value for output_type: {}", other);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message))
        }
    }
}
