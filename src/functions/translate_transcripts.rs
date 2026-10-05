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
use exacto::consensus::prelude as consensus;
use exacto::core::prelude as core;
use exacto::integrator::prelude as integrator;
use exacto::stitcher::io::loaders::load_stitched_transcript_records;
use exacto::stitcher::io::records::StitchedTranscriptRecord;
use exacto::translator::prelude as translator;
use polars::prelude::DataFrame;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::str::FromStr;
use std::path::{Path, PathBuf};

use super::translator_error;


/// Translate the assembled-transcript pipeline outputs (the five
/// post-caller/integrator TSVs) into proteoforms.
///
/// Transcript sequences and their supporting reads come from exactly one of two files, since
/// there are two ways to arrive at an assembled transcript:
/// - `assembled_transcript_support_tsv_file` — transcripts assembled outside exacto
///   (RATTLE, RNA-Bloom2), named by that pipeline.
/// - `rna_consensus_tsv_file` — `determine-rna-consensus`' own output, one consensus per read
///   cluster, named by cluster id.
///
/// `stitched_transcripts_tsv_file` is optional: the `stitch-reference-rnas` table of the same
/// transcripts. The other TSVs must then come from the second `call-rna-transcript-vars` pass,
/// run on the realigned stitched FASTA, because their coordinates index the stitched sequences.
///
/// `output_type` selects the return shape:
/// - `"file"`     → write TSV and (if `output_fasta_file` is non-empty) FASTA;
///                  return an empty `PyDataFrame`.
/// - `"dataframe"` → build a polars DataFrame in memory and return it;
///                  no files are written.
#[pyfunction]
pub fn translate_transcripts(
    _py: Python,
    assembled_transcript_support_tsv_file: String,
    rna_consensus_tsv_file: String,
    stitched_transcripts_tsv_file: String,
    assembled_transcript_model_alignments_tsv_file: String,
    assembled_transcript_variants_tsv_file: String,
    dna_variants_tsv_file: String,
    integrated_variants_tsv_file: String,
    strategy: String,
    start_codons: Vec<String>,
    output_dir: String,
    output_prefix: String,
    num_threads: usize,
    output_type: String,
) -> PyResult<(PyDataFrame, PyDataFrame)> {
    // Step 1. Load every record stream. Transcript support comes from whichever of the two
    // sources was supplied; requiring exactly one keeps the transcript sequences and their
    // read names from being drawn from two files that disagree.
    let assembled_transcript_support_records: Vec<translator::AssembledTranscriptSupportRecord> =
        match (
            assembled_transcript_support_tsv_file.is_empty(),
            rna_consensus_tsv_file.is_empty()
        ) {
            (false, true) => translator::load_assembled_transcript_support_records(
                &assembled_transcript_support_tsv_file
            ).map_err(translator_error)?,
            (true, false) => consensus::load_consensus_sequence_records(&rna_consensus_tsv_file)
                .into_iter()
                .map(|record| translator::AssembledTranscriptSupportRecord::from_consensus(
                    record.cluster_id,
                    record.consensus_sequence,
                    record.read_names
                ))
                .collect(),
            (true, true) => {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    "No transcript sequences supplied: pass exactly one of \
                     assembled_transcript_support_tsv_file (externally assembled transcripts) \
                     or rna_consensus_tsv_file (output of determine-rna-consensus)."
                ));
            }
            (false, false) => {
                return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                    "Both assembled_transcript_support_tsv_file and rna_consensus_tsv_file \
                     were supplied; pass exactly one, so the transcript sequences and their \
                     read names come from a single source."
                ));
            }
        };
    // Reference-stitched transcripts, optional. `stitch-reference-rnas` writes one row per
    // transcript, keyed by the name its sequence was aligned under, with the stitched sequence
    // and the terminal spans pasted in from the reference transcript. The second-pass tables
    // index the STITCHED sequence, so a stitched row replaces its support row's sequence. A
    // support row without a stitched row is dropped: that transcript never entered the stitched
    // FASTA, so it cannot be in the second-pass tables; if it is, the files come from different
    // runs and `build_transcript_set` fails on the unmatched name.
    let mut reference_stitched_spans: Vec<translator::ReferenceStitchedSpans> = Vec::new();
    let assembled_transcript_support_records: Vec<translator::AssembledTranscriptSupportRecord> =
        if stitched_transcripts_tsv_file.is_empty() {
            assembled_transcript_support_records
        } else {
            let mut stitched_transcript_records: HashMap<Box<str>, StitchedTranscriptRecord> = HashMap::new();
            for record in load_stitched_transcript_records(&stitched_transcripts_tsv_file) {
                if let Some(duplicate) = stitched_transcript_records.insert(record.read_name.clone(), record) {
                    return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                        "Duplicate row for transcript {} in {} — the stitcher writes exactly one \
                         row per transcript; the file is concatenated or corrupt.",
                        duplicate.read_name, stitched_transcripts_tsv_file
                    )));
                }
            }

            let mut stitched_support_records: Vec<translator::AssembledTranscriptSupportRecord> = Vec::new();
            for mut support_record in assembled_transcript_support_records {
                let record: &StitchedTranscriptRecord = match stitched_transcript_records
                    .get(&support_record.assembled_transcript_name) {
                    Some(record) => record,
                    None => continue
                };
                if record.stitched_sequence.len() != record.stitched_length as usize
                    || record.five_prime_stitch_end > record.three_prime_stitch_start
                    || record.three_prime_stitch_end != record.stitched_length {
                    return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                        "Transcript {}'s stitch coordinates disagree with its stitched sequence in {} — \
                         the file is edited or corrupt.",
                        record.read_name, stitched_transcripts_tsv_file
                    )));
                }

                // Half-open `[start, end)` spans of the stitched sequence; an unstitched end
                // has an empty span, which is dropped.
                let mut intervals: Vec<(u32, u32)> = Vec::new();
                if record.five_prime_stitch_start < record.five_prime_stitch_end {
                    intervals.push((record.five_prime_stitch_start, record.five_prime_stitch_end));
                }
                if record.three_prime_stitch_start < record.three_prime_stitch_end {
                    intervals.push((record.three_prime_stitch_start, record.three_prime_stitch_end));
                }
                reference_stitched_spans.push(translator::ReferenceStitchedSpans {
                    assembled_transcript_name: record.read_name.clone(),
                    intervals: intervals,
                    stitched_length: record.stitched_length
                });

                support_record.sequence = record.stitched_sequence.clone();
                stitched_support_records.push(support_record);
            }
            stitched_support_records
        };
    let assembled_transcript_model_alignments_records: Vec<caller::AssembledTranscriptModelAlignmentRecord> =
        caller::load_assembled_transcript_model_alignment_records(&assembled_transcript_model_alignments_tsv_file);
    let assembled_transcript_variant_records: Vec<caller::AssembledTranscriptVariantRecord> =
        caller::load_assembled_transcript_variant_records(&assembled_transcript_variants_tsv_file);
    let dna_variant_records: Vec<caller::DNAVariantRecord> =
        caller::load_dna_variant_records(&dna_variants_tsv_file);
    let integrated_variant_records: Vec<integrator::IntegratedVariantRecord> =
        integrator::load_integrated_variant_records(&integrated_variants_tsv_file);

    // Step 2. Translate. `translate_transcripts` builds the TranscriptSet
    // from the five record streams and runs translation.
    let translation_strategy: translator::TranslationStrategy =
        translator::TranslationStrategy::from_str(&strategy).map_err(|_| PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("Unsupported value for strategy: {:?}.", strategy)
        ))?;
    let start_codons_set: HashSet<&str> = start_codons.iter().map(|s| s.as_str()).collect();
    let transcript_set: translator::AssembledTranscriptSet = translator::translate_transcripts(
        &assembled_transcript_support_records,
        &assembled_transcript_model_alignments_records,
        &assembled_transcript_variant_records,
        &dna_variant_records,
        &integrated_variant_records,
        &reference_stitched_spans,
        translation_strategy,
        &start_codons_set,
        num_threads,
    ).map_err(translator_error)?;

    // Step 3. Branch on the requested output shape.
    match output_type.as_str() {
        "file" => {
            let proteoforms_tsv_file: PathBuf = if output_prefix.is_empty() {
                Path::new(&output_dir).join("exacto_proteoforms.tsv")
            } else {
                Path::new(&output_dir).join(format!("{}_exacto_proteoforms.tsv", output_prefix))
            };

            let proteoform_nucleotides_tsv_file: PathBuf = if output_prefix.is_empty() {
                Path::new(&output_dir).join("exacto_proteoform_nucleotides.tsv")
            } else {
                Path::new(&output_dir).join(format!("{}_exacto_proteoform_nucleotides.tsv", output_prefix))
            };

            let fasta_file: PathBuf = if output_prefix.is_empty() {
                Path::new(&output_dir).join("exacto_proteoforms.fasta")
            } else {
                Path::new(&output_dir).join(format!("{}_exacto_proteoforms.fasta", output_prefix))
            };

            core::write_tsv_file(
                translator::build_proteoform_records(&transcript_set),
                &proteoforms_tsv_file
            ).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Failed to write proteoforms TSV: {}", e)
            ))?;

            core::write_tsv_file(
                translator::build_nucleotide_records(&transcript_set),
                &proteoform_nucleotides_tsv_file
            ).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("Failed to write nucleotides TSV: {}", e)
            ))?;

            write_proteoforms_fasta_file(&transcript_set, &fasta_file)?;

            Ok((PyDataFrame(DataFrame::empty()), PyDataFrame(DataFrame::empty())))
        }
        "dataframe" => {
            // The records stream into the column builders; collecting them first would hold
            // every nucleotide record at once.
            let df_proteoforms: DataFrame = translator::proteoform_records_to_dataframe(
                translator::build_proteoform_records(&transcript_set)
            );
            let df_proteoform_nucleotides: DataFrame = translator::nucleotide_records_to_dataframe(
                translator::build_nucleotide_records(&transcript_set)
            );

            Ok((PyDataFrame(df_proteoforms), PyDataFrame(df_proteoform_nucleotides)))
        }
        other => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!(
                "Unsupported value for output_type: {:?}. Expected \"file\" or \"dataframe\".",
                other
            ),
        )),
    }
}


/// Write one FASTA record per `Proteoform` across all transcripts.
/// Header format: `>{transcript_id}|orf_{orf_start}-{orf_end}`.
fn write_proteoforms_fasta_file(
    transcript_set: &translator::AssembledTranscriptSet,
    path: &Path,
) -> PyResult<()> {
    let file: File = File::create(path)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    let mut writer = BufWriter::new(file);
    for transcript in transcript_set.iter() {
        for proteoform in transcript.proteoforms.iter() {
            writeln!(
                writer,
                ">{}|orf_{}-{}",
                transcript.get_id(),
                proteoform.get_orf_start(),
                proteoform.get_orf_end(),
            )
            .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
            writeln!(writer, "{}", proteoform.get_sequence())
                .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
        }
    }
    writer.flush()
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;
    Ok(())
}
