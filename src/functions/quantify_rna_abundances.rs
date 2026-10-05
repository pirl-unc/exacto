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

use exacto::caller::prelude as caller;
use exacto::cluster::prelude as cluster;
use exacto::core::prelude as core;
use exacto::quantifier::pipeline::options::QuantifyRNAAbundancesOptions;
use exacto::quantifier::prelude as quantifier;
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use std::collections::{HashMap,HashSet};
use std::path::{Path, PathBuf};

use crate::functions::io_error;


#[pyfunction]
pub fn quantify_rna_abundances(
    py: Python,
    clusters_tsv_file: String,
    cluster_reference_matches_tsv_file: String,
    cluster_splice_junctions_tsv_file: String,
    cluster_variants_tsv_file: String,
    consensus_rna_reference_matches_tsv_file: String,
    output_dir: String,
    output_prefix: String,
    pseudo_count: f64,
    max_iter: usize,
    tol: f64,
    output_type: String
) -> PyResult<(PyDataFrame, PyDataFrame)> {
    let _ = py;
    // Load the cluster and read IDs
    let cluster_id_records: Vec<cluster::RNAReadClusterIDRecord> = cluster::load_rna_read_cluster_id_records(
        clusters_tsv_file.as_str()
    );
    let mut clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::new();
    for record in cluster_id_records.iter() {
        clusters
            .entry(record.cluster_id)
            .or_insert_with(HashSet::new)
            .insert(record.read_name.clone());
    }

    // Quantify RNA abundances.
    // Every field is a parameter, so the Python defaults are the only defaults.
    let options: QuantifyRNAAbundancesOptions = QuantifyRNAAbundancesOptions {
        pseudo_count: pseudo_count,
        max_iter: max_iter,
        tol: tol
    };
    let cluster_quantification_set: quantifier::ClusterQuantificationSet = quantifier::quantify_rna_abundances(
        &clusters,
        &options
    );

    // Load RNAReadClusterSet
    let mut rna_read_cluster_set: cluster::RNAReadClusterSet = cluster::load_rna_read_cluster_set(
        clusters_tsv_file.as_str(),
        cluster_reference_matches_tsv_file.as_str(),
        cluster_splice_junctions_tsv_file.as_str(),
        cluster_variants_tsv_file.as_str()
    );

    // Load the consensus RNA reference transcript match records
    let consensus_rna_reference_match_records: Vec<caller::AssembledTranscriptReferenceTranscriptMatchRecord> = caller::load_assembled_transcript_reference_transcript_match_records(
        consensus_rna_reference_matches_tsv_file.as_str()
    );

    let mut consensus_reference_gene_transcript_ids: HashMap<usize, HashSet<(Box<str>, Box<str>)>> = HashMap::new();
    for record in consensus_rna_reference_match_records.iter() {
        let cluster_id: usize = record.assembled_transcript_name.parse::<usize>()
            .unwrap_or_else(|_| panic!(
                "Assembled transcript name is not a cluster id in {}: {}",
                consensus_rna_reference_matches_tsv_file,
                record.assembled_transcript_name
            ));
        consensus_reference_gene_transcript_ids
            .entry(cluster_id)
            .or_default()
            .insert((
                record.reference_gene_id.clone(),
                record.reference_transcript_id.clone()
            ));
    }

    // The consensus sequence's reference matches replace the raw reads' per cluster.
    for (cluster_id, reference_gene_transcript_ids) in consensus_reference_gene_transcript_ids {
        rna_read_cluster_set.set_reference_gene_transcript_ids(cluster_id, reference_gene_transcript_ids);
    }

    match output_type.as_str() {
        "dataframe" => {
            let df_abundances: DataFrame = quantifier::quantification_cluster_records_to_dataframe(
                quantifier::build_quantification_cluster_records(
                    &cluster_quantification_set,
                    &rna_read_cluster_set
                )
            );
            let df_reference_transcripts: DataFrame = quantifier::quantification_reference_transcript_records_to_dataframe(
                quantifier::build_quantification_reference_transcript_records(
                    &cluster_quantification_set,
                    &rna_read_cluster_set
                )
            );

            Ok((PyDataFrame(df_abundances),
                PyDataFrame(df_reference_transcripts)))
        }
        "file" => {
            let abundances_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_rna_abundances.tsv", output_prefix));
            let abundances_reference_transcripts_tsv_file: PathBuf = Path::new(&output_dir)
                .join(format!("{}_exacto_rna_abundances_reference_transcripts.tsv", output_prefix));

            core::write_tsv_file(
                quantifier::build_quantification_cluster_records(
                    &cluster_quantification_set,
                    &rna_read_cluster_set
                ),
                &abundances_tsv_file
            ).map_err(io_error)?;
            core::write_tsv_file(
                quantifier::build_quantification_reference_transcript_records(
                    &cluster_quantification_set,
                    &rna_read_cluster_set
                ),
                &abundances_reference_transcripts_tsv_file
            ).map_err(io_error)?;

            Ok((PyDataFrame(DataFrame::new(vec![]).unwrap()),
                PyDataFrame(DataFrame::new(vec![]).unwrap())))
        }
        other => {
            let error_message = format!("Unsupported value for output_type: {}", other);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message))
        }
    }
}
