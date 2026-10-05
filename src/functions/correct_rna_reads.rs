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
extern crate pyo3;

use exacto::cluster::prelude as cluster;
use exacto::core::prelude as core;
use exacto::correction::prelude as correction;
use pyo3::prelude::*;
use std::path::{Path, PathBuf};


/// Error-correct RNA reads against the cluster each was assigned to.
///
/// The corrected reads are streamed to
/// `{output_dir}/{output_prefix}_exacto_rna_corrected_reads.fastq.gz` (BGZF) as correction
/// runs; the path of that file is returned.
#[pyfunction]
pub fn correct_rna_reads(
    py: Python,
    bam_file: String,
    reference_gene_annotation_file: String,
    reference_gene_annotation_source: String,
    reference_gene_annotation_assembly: String,
    reference_gene_annotation_version: String,
    clusters_tsv_file: String,
    cluster_reference_transcripts_tsv_file: String,
    cluster_splice_junctions_tsv_file: String,
    cluster_variants_tsv_file: String,
    output_dir: String,
    output_prefix: String,
    max_correctable_event_len: usize,
    corrected_base_quality: u8,
    trim_transcript_ends: bool,
    num_threads: usize,
    chunk_size: usize
) -> PyResult<String> {
    let _ = py;
    // Clusters come in as the four TSVs `cluster-rna-reads` emits rather than being recomputed,
    // the same way `quantify-rna-abundances` takes them.
    let cluster_set: cluster::RNAReadClusterSet = cluster::load_rna_read_cluster_set(
        clusters_tsv_file.as_str(),
        cluster_reference_transcripts_tsv_file.as_str(),
        cluster_splice_junctions_tsv_file.as_str(),
        cluster_variants_tsv_file.as_str()
    );

    // The corrected FASTQ is written as correction runs, so it needs a home.
    if output_dir.is_empty() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "correct-rna-reads needs --output-dir: the corrected FASTQ is always written to a file."
        ));
    }
    // The handoff to `determine-rna-consensus`, which takes a FASTQ and resolves sequences by
    // read name. BGZF, not plain gzip.
    let corrected_reads_fastq_file: PathBuf = Path::new(&output_dir)
        .join(format!("{}_exacto_rna_corrected_reads.fastq.gz", output_prefix));

    // Only transcript-end trimming reads the annotation, so it is loaded only for trimming.
    let gene_annotator: Option<core::Gencode> = if !trim_transcript_ends {
        None
    } else if reference_gene_annotation_file.is_empty() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "--trim-transcript-ends needs --reference-gene-annotation-file."
        ));
    } else if reference_gene_annotation_source.as_str() == "gencode" {
        Some(core::Gencode::new(
            reference_gene_annotation_file.as_str(),
            reference_gene_annotation_assembly.as_str(),
            reference_gene_annotation_version.as_str(),
            None,
            None,
            None,
            None
        ))
    } else {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Unsupported annotation source: {}", reference_gene_annotation_source
        )));
    };

    // Every field is a parameter, so the Python defaults are the only defaults.
    let options: correction::CorrectRNAReadsOptions = correction::CorrectRNAReadsOptions {
        max_correctable_event_len: max_correctable_event_len,
        corrected_base_quality: corrected_base_quality,
        trim_transcript_ends: trim_transcript_ends
    };

    correction::correct_rna_reads(
        bam_file.as_str(),
        corrected_reads_fastq_file.to_str().unwrap(),
        &cluster_set,
        gene_annotator.as_ref().map(|gene_annotator| gene_annotator as &(dyn core::GeneAnnotator + Sync)),
        &options,
        num_threads,
        chunk_size
    ).map_err(|error| PyErr::new::<pyo3::exceptions::PyValueError, _>(
        format!("{}; the clusters are those of {}", error, clusters_tsv_file)
    ))?;

    Ok(corrected_reads_fastq_file.to_string_lossy().into_owned())
}
