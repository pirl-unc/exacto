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

use exacto::core::prelude as core;
use exacto::stitcher::io::builders::build_stitched_transcript_records;
use exacto::stitcher::io::dataframes::stitched_transcript_records_to_dataframe;
use exacto::stitcher::io::records::StitchedTranscriptRecord;
use exacto::stitcher::prelude as stitcher;
use polars::prelude::DataFrame;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::functions::io_error;


/// Stitch reference transcript sequence onto the degraded ends of aligned RNA sequences.
///
/// Every mapped read in the BAM file is one transcript, identified by its read name. The
/// stitched sequences' coordinates deliberately diverge from any `call-rna-transcript-vars`
/// output computed on the same BAM: realign the stitched FASTA and re-run
/// `call-rna-transcript-vars` before translating.
#[pyfunction]
pub fn stitch_reference_transcripts(
    _py: Python,
    bam_file: String,
    reference_genome_fasta_file: String,
    reference_gene_annotation_file: String,
    reference_gene_annotation_source: String,
    reference_gene_annotation_assembly: String,
    reference_gene_annotation_version: String,
    min_mapping_quality: u16,
    min_num_splice_junction_matches: usize,
    pas_search_size: usize,
    pas_hexamers: Vec<String>,
    pas_start_offset_range: (usize, usize),
    polya_tail_min_adenosine_fraction: f64,
    polya_window_size: usize,
    polya_min_consecutive_adenosine: usize,
    num_threads: usize,
    output_dir: String,
    output_prefix: String,
    output_type: String
) -> PyResult<PyDataFrame> {
    let gene_annotator = if reference_gene_annotation_source.as_str() == "gencode" {
        core::Gencode::new(
            reference_gene_annotation_file.as_str(),
            reference_gene_annotation_assembly.as_str(),
            reference_gene_annotation_version.as_str(),
            None,
            None,
            None,
            None
        )
    } else {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Unsupported annotation source: {}", reference_gene_annotation_source
        )));
    };

    let (pas_start_offset_min, pas_start_offset_max): (usize, usize) = pas_start_offset_range;
    if pas_start_offset_min > pas_start_offset_max {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Unsupported value for pas_start_offset_range: ({}, {}). Expected (min, max) with min <= max.",
            pas_start_offset_min, pas_start_offset_max
        )));
    }

    // Every field is a parameter, so the Python defaults are the only defaults.
    let options: stitcher::StitchReferenceTranscriptsOptions = stitcher::StitchReferenceTranscriptsOptions {
        min_mapping_quality: min_mapping_quality,
        min_num_splice_junction_matches: min_num_splice_junction_matches,
        pas_search_size: pas_search_size,
        pas_hexamers: pas_hexamers.into_iter().map(|hexamer| hexamer.into_boxed_str()).collect(),
        pas_start_offset_range: pas_start_offset_min..=pas_start_offset_max,
        polya_tail_min_adenosine_fraction: polya_tail_min_adenosine_fraction,
        polya_window_size: polya_window_size,
        polya_min_consecutive_adenosine: polya_min_consecutive_adenosine
    };

    let stitched_transcript_set: stitcher::StitchedTranscriptSet = stitcher::stitch_reference_transcripts(
        bam_file.as_str(),
        reference_genome_fasta_file.as_str(),
        &gene_annotator,
        &options,
        num_threads
    );

    // A stitched transcript carries its read ID, not its read name. Read IDs are assigned in
    // BAM order, so indexing the BAM the way the pipeline does recovers the same map.
    let (_, read_names_map) = core::index_bam_records(bam_file.as_str(), true, num_threads);
    let stitched_transcript_records: Vec<StitchedTranscriptRecord> = build_stitched_transcript_records(
        &stitched_transcript_set,
        &read_names_map
    ).collect();

    match output_type.as_str() {
        "dataframe" => {
            let df_stitched: DataFrame = stitched_transcript_records_to_dataframe(stitched_transcript_records);
            Ok(PyDataFrame(df_stitched))
        }
        "file" => {
            let stitched_tsv_file: PathBuf = if output_prefix.is_empty() {
                Path::new(&output_dir).join("exacto_stitched_transcripts.tsv")
            } else {
                Path::new(&output_dir).join(format!("{}_exacto_stitched_transcripts.tsv", output_prefix))
            };
            let stitched_fasta_file: PathBuf = if output_prefix.is_empty() {
                Path::new(&output_dir).join("exacto_stitched_transcripts.fasta")
            } else {
                Path::new(&output_dir).join(format!("{}_exacto_stitched_transcripts.fasta", output_prefix))
            };

            core::write_tsv_file(stitched_transcript_records.iter(), &stitched_tsv_file).map_err(io_error)?;
            write_stitched_fasta_file(&stitched_transcript_records, &stitched_fasta_file)?;

            Ok(PyDataFrame(DataFrame::empty()))
        }
        other => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Unsupported value for output_type: {:?}. Expected \"file\" or \"dataframe\".",
            other
        )))
    }
}


/// Write one FASTA record per stitched transcript, named by its read name. This is the file
/// to realign for the second `call-rna-transcript-vars` pass.
fn write_stitched_fasta_file(
    stitched_transcript_records: &Vec<StitchedTranscriptRecord>,
    path: &Path
) -> PyResult<()> {
    let file: File = File::create(path).map_err(io_error)?;
    let mut writer = BufWriter::new(file);
    for record in stitched_transcript_records.iter() {
        // An empty FASTA record is unindexable and unalignable, so skip it here; the row
        // stays in the TSV.
        if record.stitched_sequence.is_empty() {
            continue;
        }
        writeln!(writer, ">{}", record.read_name).map_err(io_error)?;
        writeln!(writer, "{}", record.stitched_sequence).map_err(io_error)?;
    }
    writer.flush().map_err(io_error)?;
    core::create_fai_file(path.to_str().unwrap()).map_err(io_error)?;
    Ok(())
}
