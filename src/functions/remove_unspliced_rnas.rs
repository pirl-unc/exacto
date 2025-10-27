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

use exacto::qc::prelude as qc;
use exacto::core::prelude as core;
use pyo3::prelude::*;
use std::collections::HashSet;


#[pyfunction]
pub fn remove_unspliced_rnas(
    py: Python,
    bam_file: String,
    bam_bai_file: String,
    fasta_file: String,
    reference_gene_annotation_file: String,
    reference_gene_annotation_source: String,
    reference_gene_annotation_assembly: String,
    reference_gene_annotation_version: String,
    gene_types: Vec<String>,
    gene_levels: Vec<u8>,
    transcript_types: Vec<String>,
    transcript_levels: Vec<u8>,
    output_bam_file: String,
    output_bam_bai_file: String,
    output_fasta_file: String,
    num_threads: usize,
    min_mapping_quality: usize
) -> PyResult<()> {
    let gene_annotator = if reference_gene_annotation_source.as_str() == "gencode" {
        let gene_types_: Option<HashSet<&str>> = (!gene_types.is_empty()).then(|| gene_types.iter().map(String::as_str).collect());
        let gene_levels_: Option<HashSet<u8>> = (!gene_levels.is_empty()).then(|| gene_levels.iter().copied().collect());
        let transcript_types_: Option<HashSet<&str>> = (!transcript_types.is_empty()).then(|| transcript_types.iter().map(String::as_str).collect());
        let transcript_levels_: Option<HashSet<u8>> = (!transcript_levels.is_empty()).then(|| transcript_levels.iter().copied().collect());
        core::Gencode::new(
            reference_gene_annotation_file.as_str(),
            reference_gene_annotation_assembly.as_str(),
            reference_gene_annotation_version.as_str(),
            gene_types_,
            gene_levels_,
            transcript_types_,
            transcript_levels_
        )
    } else {
        panic!("Unsupported annotation source: {}", reference_gene_annotation_source);
    };
    qc::remove_unspliced_rnas(
        bam_file.as_str(),
        bam_bai_file.as_str(),
        fasta_file.as_str(),
        &gene_annotator,
        output_bam_file.as_str(),
        output_bam_bai_file.as_str(),
        output_fasta_file.as_str(),
        num_threads,
        min_mapping_quality
    );
    Ok(())
}
