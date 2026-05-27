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
use exacto::core::prelude as core;
use exacto::integrator::prelude as integrator;
use exacto::translator::prelude as translator;
use polars::prelude::DataFrame;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::str::FromStr;
use std::path::{Path, PathBuf};


/// Translate the assembled-transcript pipeline outputs (the five
/// post-caller/integrator TSVs) into primary structures.
///
/// `output_type` selects the return shape:
/// - `"file"`     → write TSV and (if `output_fasta_file` is non-empty) FASTA;
///                  return an empty `PyDataFrame`.
/// - `"dataframe"` → build a polars DataFrame in memory and return it;
///                  no files are written.
#[pyfunction]
pub fn translate_structures(
    _py: Python,
    rna_assembly_support_tsv_file: String,
    transcript_model_structures_tsv_file: String,
    rna_variants_tsv_file: String,
    dna_variants_tsv_file: String,
    integrated_variants_tsv_file: String,
    strategy: String,
    start_codons: Vec<String>,
    output_dir: String,
    output_prefix: String,
    num_threads: usize,
    output_type: String,
) -> PyResult<(PyDataFrame, PyDataFrame)> {
    // Step 1. Load every record stream.
    let assembled_transcript_support_records: Vec<translator::AssembledTranscriptSupportRecord> =
        translator::load_assembled_transcript_support_records(&rna_assembly_support_tsv_file);
    let transcript_model_structure_records: Vec<caller::TranscriptModelStructureRecord> =
        caller::load_transcript_model_structure_records(&transcript_model_structures_tsv_file);
    let rna_variant_records: Vec<caller::RNAVariantRecord> =
        caller::load_rna_variant_records(&rna_variants_tsv_file);
    let dna_variant_records: Vec<caller::DNAVariantRecord> =
        caller::load_dna_variant_records(&dna_variants_tsv_file);
    let integrated_variant_records: Vec<integrator::IntegratedVariantRecord> =
        integrator::load_integrated_variant_records(&integrated_variants_tsv_file);

    // Step 2. Translate. `translate_structures` builds the TranscriptSet
    // from the five record streams and runs translation.
    let translation_strategy: translator::TranslationStrategy =
        translator::TranslationStrategy::from_str(&strategy).unwrap();
    let start_codons_set: HashSet<&str> = start_codons.iter().map(|s| s.as_str()).collect();
    let transcript_set: translator::TranscriptSet = translator::translate_structures(
        &assembled_transcript_support_records,
        &transcript_model_structure_records,
        &rna_variant_records,
        &dna_variant_records,
        &integrated_variant_records,
        translation_strategy,
        start_codons_set,
        num_threads,
    );

    // Step 3. Branch on the requested output shape.
    match output_type.as_str() {
        "file" => {
            let primary_structures_tsv_file: PathBuf = if output_prefix.is_empty() {
                Path::new(&output_dir).join("exacto_primary_structures.tsv")
            } else {
                Path::new(&output_dir).join(format!("{}_exacto_primary_structures.tsv", output_prefix))
            };

            let nucleotide_tsv_file: PathBuf = if output_prefix.is_empty() {
                Path::new(&output_dir).join("exacto_primary_structure_nucleotides.tsv")
            } else {
                Path::new(&output_dir).join(format!("{}_exacto_primary_structure_nucleotides.tsv", output_prefix))
            };

            let fasta_file: PathBuf = if output_prefix.is_empty() {
                Path::new(&output_dir).join("exacto_primary_structures.fasta")
            } else {
                Path::new(&output_dir).join(format!("{}_exacto_primary_structures.fasta", output_prefix))
            };

            core::write_tsv_file(
                translator::build_primary_structure_records(&transcript_set),
                &primary_structures_tsv_file
            );

            core::write_tsv_file(
                translator::build_nucleotide_records(&transcript_set),
                &nucleotide_tsv_file
            );

            write_protein_fasta(&transcript_set, &fasta_file)?;

            Ok((PyDataFrame(DataFrame::empty()), PyDataFrame(DataFrame::empty())))
        }
        "dataframe" => {
            let ps_records: Vec<translator::PrimaryStructureRecord> =
                translator::build_primary_structure_records(&transcript_set).collect();
            let nuc_records: Vec<translator::NucleotideRecord> =
                translator::build_nucleotide_records(&transcript_set).collect();

            let df_ps: DataFrame = translator::primary_structure_records_to_dataframe(ps_records);
            let df_ps_nucleotides: DataFrame = translator::nucleotide_records_to_dataframe(nuc_records);

            Ok((PyDataFrame(df_ps), PyDataFrame(df_ps_nucleotides)))
        }
        other => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!(
                "Unsupported value for output_type: {:?}. Expected \"file\" or \"dataframe\".",
                other
            ),
        )),
    }
}


/// Write one FASTA record per `PrimaryStructure` across all transcripts.
/// Header format: `>{transcript_id}|orf_{orf_start}-{orf_end}`.
fn write_protein_fasta(
    transcript_set: &translator::TranscriptSet,
    path: &Path,
) -> PyResult<()> {
    let file: File = File::create(path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    let mut writer = BufWriter::new(file);
    for transcript in transcript_set.iter() {
        for primary_structure in transcript.primary_structures.iter() {
            writeln!(
                writer,
                ">{}|orf_{}-{}",
                transcript.get_id(),
                primary_structure.get_orf_start(),
                primary_structure.get_orf_end(),
            )
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
            for amino_acid in primary_structure.amino_acids.iter() {
                write!(writer, "{}", amino_acid.get_amino_acid())
                    .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
            }
            writeln!(writer)
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        }
    }
    writer.flush()
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    Ok(())
}
