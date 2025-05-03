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


extern crate chrono;
extern crate env_logger;
extern crate exacto;
extern crate exitcode;
extern crate flate2;
extern crate log;
extern crate noodles_bgzf;
extern crate noodles_fasta;
extern crate noodles_fastq;
extern crate polars;
extern crate pyo3_polars;
extern crate pyo3;
extern crate sysinfo;
extern crate tempfile;

use chrono::Local;
use env_logger::{Builder, Env};
use exacto::caller::prelude as caller;
use exacto::translator::prelude as translator;
use exacto::util::prelude as util;
use flate2::read::GzDecoder;
use log::{info, LevelFilter};
use polars::prelude::*;
use polars::io::ipc::IpcWriter;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use pyo3_polars::PyDataFrame;
use noodles_bgzf as bgzf;
use noodles_fasta::{self as fasta, record::{Definition, Sequence}};
use noodles_fastq as fastq;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufReader,Cursor,Read,Write};
use std::path::Path;
use sysinfo::System;
use tempfile::NamedTempFile;


pub fn capture_memory_usage(message: &str) {
    let mut sys = System::new_all();
    sys.refresh_all();
    let pid = sysinfo::get_current_pid().unwrap();
    if let Some(process) = sys.process(pid) {
        let memory_usage = process.memory();
        let memory_usage_gb = memory_usage as f64 / (1024.0 * 1024.0 * 1024.0);
        println!("{}: {:.2} GB", message, memory_usage_gb);
    } else {
        println!("Could not get process memory usage");
    }
}


/// This function identifies DNA variants in a long-read WGS BAM file.
#[pyfunction]
fn identify_dna_variants(
    py: Python,
    bam_file: String,
    bam_bai_file: String,
    output_tsv_file: String,
    gzip: bool,
    min_reads: usize,
    min_mapping_quality: usize,
    min_average_base_quality: f32,
    min_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    max_intrachromosomal_distance_tau: usize,
    max_intrachromosomal_distance: usize,
    max_interchromosomal_distance: usize,
    num_threads: usize,
    chromosomes: Vec<String>,
    temp_dir: String,
    output_type: String
) -> PyResult<PyDataFrame> {
    let variant_call_set: caller::VariantCallSet = caller::identify_dna_variants(
        bam_file.as_str(),
        bam_bai_file.as_str(),
        min_reads,
        min_mapping_quality,
        min_average_base_quality,
        min_size_proportion,
        max_ins_norm_edit_distance,
        max_intrachromosomal_distance_tau as u32,
        max_intrachromosomal_distance as u32,
        max_interchromosomal_distance as u32,
        num_threads,
        chromosomes.iter().map(|s| s.as_str()).collect(),
        temp_dir.as_str()
    );
    capture_memory_usage("Successfully ran Exacto DNA variant calling (1)");
    match output_type.as_str() {
        "dataframe" => {
            Ok(PyDataFrame(variant_call_set.to_dataframe(num_threads)))
        }
        "file" => {
            variant_call_set.to_tsv(
                output_tsv_file.as_str(),
                100_000,
                num_threads,
                gzip,
            );
            capture_memory_usage("Successfully wrote to TSV file");
            Ok(PyDataFrame(DataFrame::new(vec![]).unwrap()))
        }
        other => {
            let error_message = format!("Unsupported value for output_type: {}", other);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message))
        }
    }
}


/// This function identifies case-specific DNA variants in a long-read WGS BAM file.
#[pyfunction]
fn identify_case_specific_dna_variants(
    py: Python,
    case_bam_file: String,
    case_bam_bai_file: String,
    control_bam_files: Vec<String>,
    control_bam_bai_files: Vec<String>,
    output_tsv_file: String,
    gzip: bool,
    min_reads: usize,
    min_mapping_quality: usize,
    min_average_base_quality: f32,
    min_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    max_intrachromosomal_distance_tau: usize,
    max_intrachromosomal_distance: usize,
    max_interchromosomal_distance: usize,
    apply_infinite_sites_assumption: bool,
    num_threads: usize,
    chromosomes: Vec<String>,
    temp_dir: String,
    output_type: String
) -> PyResult<PyDataFrame> {
    let variant_call_set: caller::VariantCallSet = caller::identify_case_specific_dna_variants(
        case_bam_file.as_str(),
        case_bam_bai_file.as_str(),
        control_bam_files.iter().map(|s| s.as_str()).collect(),
        control_bam_bai_files.iter().map(|s| s.as_str()).collect(),
        min_reads,
        min_mapping_quality,
        min_average_base_quality,
        min_size_proportion,
        max_ins_norm_edit_distance,
        max_intrachromosomal_distance_tau as u32,
        max_intrachromosomal_distance as u32,
        max_interchromosomal_distance as u32,
        apply_infinite_sites_assumption,
        num_threads,
        chromosomes.iter().map(|s| s.as_str()).collect(),
        temp_dir.as_str()
    );
    capture_memory_usage("Successfully ran Exacto DNA variant calling (1)");
    match output_type.as_str() {
        "dataframe" => {
            Ok(PyDataFrame(variant_call_set.to_dataframe(num_threads)))
        }
        "file" => {
            variant_call_set.to_tsv(
                output_tsv_file.as_str(),
                100_000,
                num_threads,
                gzip,
            );
            capture_memory_usage("Successfully wrote to TSV file");
            Ok(PyDataFrame(DataFrame::new(vec![]).unwrap()))
        }
        other => {
            let error_message = format!("Unsupported value for output_type: {}", other);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message))
        }
    }
}


/// This function identifies RNA variants in a long-read RNA BAM file.
#[pyfunction]
fn identify_rna_variants(
    py: Python,
    bam_file: String,
    bam_bai_file: String,
    reference_genome_fasta_file: String,
    gene_annotation_file: String,
    gene_annotation_source: String,
    output_dir: String,
    output_prefix: String,
    reference_transcript_scoring_method: String,
    reference_transcript_selection_strategy: String,
    reference_transcript_top_k: usize,
    reference_transcript_threshold: f32,
    min_mapping_quality: usize,
    min_average_base_quality: f32,
    num_threads: usize,
    temp_dir: String,
    output_type: String
) -> PyResult<(PyDataFrame,PyDataFrame,PyDataFrame,PyDataFrame,PyDataFrame,PyDataFrame,PyDataFrame)> {
    let gene_annotator = if gene_annotation_source.as_str() == "gencode" {
        util::Gencode::new(gene_annotation_file.as_str(), "hg38")
    } else {
        panic!("Unsupported annotation source: {}", gene_annotation_source);
    };

    let scoring_method: caller::ReferenceTranscriptScoringMethod = reference_transcript_scoring_method.as_str().parse().unwrap();
    let selection_strategy: caller::ReferenceTranscriptSelectionStrategy = reference_transcript_selection_strategy.as_str().parse().unwrap();

    let transcript_model_set: caller::TranscriptModelSet = caller::identify_variant_transcripts(
        bam_file.as_str(),
        bam_bai_file.as_str(),
        reference_genome_fasta_file.as_str(),
        &gene_annotator,
        scoring_method,
        selection_strategy,
        reference_transcript_top_k,
        reference_transcript_threshold,
        min_mapping_quality,
        min_average_base_quality,
        num_threads
    );

    match output_type.as_str() {
        "dataframe" => {
            let df_exons: DataFrame = transcript_model_set.get_exons_dataframe();
            let df_read_filter_status: DataFrame = transcript_model_set.get_read_filter_status_dataframe();
            let df_read_names: DataFrame = transcript_model_set.get_read_names_dataframe();
            let df_reference_matches: DataFrame = transcript_model_set.get_matched_reference_transcripts_dataframe();
            let df_splice_junctions: DataFrame = transcript_model_set.get_splice_junctions_dataframe();
            let df_transcripts: DataFrame = transcript_model_set.get_transcripts_dataframe();
            let df_variant_calls: DataFrame = transcript_model_set.get_variant_calls_dataframe(num_threads);
            Ok((PyDataFrame(df_exons),
                PyDataFrame(df_read_filter_status),
                PyDataFrame(df_read_names),
                PyDataFrame(df_reference_matches),
                PyDataFrame(df_splice_junctions),
                PyDataFrame(df_transcripts),
                PyDataFrame(df_variant_calls)))
        }
        "file" => {
            transcript_model_set.to_tsv_files(
                output_dir.as_str(),
                output_prefix.as_str()
            );
            Ok((PyDataFrame(DataFrame::new(vec![]).unwrap()),
                PyDataFrame(DataFrame::new(vec![]).unwrap()),
                PyDataFrame(DataFrame::new(vec![]).unwrap()),
                PyDataFrame(DataFrame::new(vec![]).unwrap()),
                PyDataFrame(DataFrame::new(vec![]).unwrap()),
                PyDataFrame(DataFrame::new(vec![]).unwrap()),
                PyDataFrame(DataFrame::new(vec![]).unwrap())))
        }
        other => {
            let error_message = format!("Unsupported value for output_type: {}", other);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message))
        }
    }
}


// /// This function identifies peptides variants.
// #[pyfunction]
// fn identify_peptide_variants(
//     py: Python,
//     fasta_file: String,
//     rna_bam_file: String,
//     rna_bam_bai_file: String,
//     reference_fasta_file: String,
//     translations_tsv_file: String,
//     rna_variants_tsv_file: String,
//     dna_variants_tsv_file: String,
//     exclude_bed_file: String,
//     min_reads: usize,
//     k: usize,
//     num_threads: usize,
//     dna_variant_padding: usize,
//     output_tsv_file: String,
//     output_fasta_file: String,
//     gzip: bool,
//     output_type: String
// ) -> PyResult<PyDataFrame> {
//     let mutant_peptides_set: caller::MutantPeptidesSet = caller::identify_peptide_variants(
//         fasta_file.as_str(),
//         rna_bam_file.as_str(),
//         rna_bam_bai_file.as_str(),
//         reference_fasta_file.as_str(),
//         translations_tsv_file.as_str(),
//         rna_variants_tsv_file.as_str(),
//         dna_variants_tsv_file.as_str(),
//         exclude_bed_file.as_str(),
//         min_reads,
//         k,
//         num_threads,
//         dna_variant_padding
//     );
//     capture_memory_usage("Successfully ran Exacto peptide variant calling (1)");
//     match output_type.as_str() {
//         "dataframe" => {
//             Ok(PyDataFrame(mutant_peptides_set.to_dataframe(num_threads)))
//         }
//         "file" => {
//             // FASTA file
//             if gzip {
//                 let output_fasta_file_ = File::create(output_fasta_file).unwrap();
//                 let bgzf_writer = bgzf::Writer::new(output_fasta_file_);
//                 let mut fasta_writer = fasta::Writer::new(bgzf_writer);
//                 for mutant_peptide in mutant_peptides_set.mutant_peptides.iter() {
//                     let definition = Definition::new(mutant_peptide.id.to_string(), None);
//                     let sequence = Sequence::from(mutant_peptide.peptide_sequence.as_bytes().to_vec());
//                     let record = fasta::Record::new(definition, sequence);
//                     fasta_writer.write_record(&record).unwrap();
//                 }
//             } else {
//                 let output_fasta_file_ = File::create(output_fasta_file).unwrap();
//                 let mut fasta_writer = fasta::Writer::new(output_fasta_file_);
//                 for mutant_peptide in mutant_peptides_set.mutant_peptides.iter() {
//                     let definition = Definition::new(mutant_peptide.id.to_string(), None);
//                     let sequence = Sequence::from(mutant_peptide.peptide_sequence.as_bytes().to_vec());
//                     let record = fasta::Record::new(definition, sequence);
//                     fasta_writer.write_record(&record).unwrap();
//                 }
//             }
//
//             // TSV file
//             mutant_peptides_set.to_tsv(
//                 output_tsv_file.as_str(),
//                 100_000,
//                 num_threads,
//                 gzip
//             );
//             capture_memory_usage("Successfully wrote to TSV file");
//             Ok(PyDataFrame(DataFrame::new(vec![]).unwrap()))
//         }
//         other => {
//             let error_message = format!("Unsupported value for output_type: {}", other);
//             Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message))
//         }
//     }
// }


/// This function translates a RNA sequence.
#[pyfunction]
fn translate_rna_sequence(
    py: Python,
    rna_sequence: String,
    strategy: String
) -> PyResult<(String,usize,usize)> {
    let rnas = vec![translator::RNA::new("sequence".into(), rna_sequence.into_boxed_str())];
    let mut translation_set: translator::TranslationSet = translator::translate(rnas, 1);
    if strategy == translator::TranslationStrategies::LONGEST_ORF {
        if translation_set.get_count() > 0 {
            let peptide: &translator::Peptide = translation_set.translations[0].get_longest_orf_peptide();
            Ok((peptide.sequence.to_string(), peptide.orf_start, peptide.orf_end))
        } else {
            Ok(("".to_string(),0,0))
        }
    } else {
        panic!("Unsupported translation strategy: {}", strategy);
    }
}


/// This function translates RNA sequences in a long-read RNA-seq FASTA file to peptides.
#[pyfunction]
fn translate_rna_fasta_file(
    py: Python,
    fasta_file: String,
    strategy: String,
    num_threads: usize,
    temp_dir: String
) -> PyResult<String> {
    // Step 1. Read the FASTA file
    let gzipped = util::is_gzipped(fasta_file.as_str());
    let file = File::open(fasta_file.as_str()).expect("Unable to open FASTA file");
    let reader: Box<dyn Read> = if gzipped {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let buffered_reader = BufReader::new(reader);
    let mut fasta_reader = fasta::Reader::new(buffered_reader);
    let mut rnas: Vec<translator::RNA> = Vec::new();
    for result in fasta_reader.records() {
        match result {
            Ok(record) => {
                // Convert &[u8] name to String (strict UTF-8)
                let name_bytes = record.name();
                let name = String::from_utf8(name_bytes.to_vec())
                    .expect("Invalid UTF-8 in FASTA record name");

                // Convert sequence to &[u8] then to String (strict UTF-8)
                let sequence_bytes = record.sequence().as_ref(); // &[u8]
                let sequence = String::from_utf8(sequence_bytes.to_vec())
                    .expect("Invalid UTF-8 in FASTA sequence");

                let rna = translator::RNA::new(name.into_boxed_str(), sequence.into_boxed_str());
                rnas.push(rna);
            }
            Err(e) => {
                panic!("Error reading record: {}", e);
            }
        }
    }

    // Step 2. Translate the RNA sequences
    let mut translation_set: translator::TranslationSet = translator::translate(
        rnas,
        num_threads
    );

    // Step 3. Write the translated peptides set to a temporary file
    let dir: String = if temp_dir.is_empty() {
        env::var("TMPDIR").unwrap_or_else(|_| env::temp_dir().to_string_lossy().to_string())
    } else {
        temp_dir.to_string()
    };
    let dir_path = Path::new(&dir);
    if !dir_path.exists() {
        panic!("Directory does not exist: {}", dir);
    }
    let temp_file = NamedTempFile::new_in(dir_path).unwrap();
    let persisted_file_path = temp_file.into_temp_path().keep().unwrap();
    let file = File::create(&persisted_file_path).unwrap();
    let mut writer = IpcWriter::new(file);
    writer.finish(&mut translation_set.to_dataframe(strategy.as_str())).unwrap();
    drop(writer);
    translation_set.translations.clear();
    translation_set.translations.shrink_to_fit();
    drop(translation_set);
    Ok(persisted_file_path.to_str().unwrap().to_string())
}


/// This function translates RNA sequences in a long-read RNA-seq FASTQ file to peptides.
#[pyfunction]
fn translate_rna_fastq_file(
    py: Python,
    fastq_file: String,
    strategy: String,
    num_threads: usize,
    temp_dir: String
) -> PyResult<String> {
    // Step 1. Read the FASTQ file
    let gzipped = util::is_gzipped(fastq_file.as_str());
    let file = File::open(fastq_file.as_str()).expect("Unable to open FASTQ file");
    let reader: Box<dyn Read> = if gzipped {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let buffered_reader = BufReader::new(reader);
    let mut fastq_reader = fastq::Reader::new(buffered_reader);
    let mut rnas: Vec<translator::RNA> = Vec::new();
    for result in fastq_reader.records() {
        match result {
            Ok(record) => {
                let sequence_result = String::from_utf8(record.sequence().to_vec());
                let sequence: String = match sequence_result {
                    Ok(seq) => seq,
                    Err(e) => {
                        panic!("Error converting sequence to UTF-8: {}", e);
                    }
                };
                let rna: translator::RNA = translator::RNA::new(
                    record.name().to_string().into_boxed_str(),
                    sequence.into_boxed_str(),
                );
                rnas.push(rna);
            }
            Err(e) => {
                panic!("Error reading record: {}", e);
            }
        }
    }

    // Step 2. Translate the RNA sequences
    let mut translation_set: translator::TranslationSet = translator::translate(
        rnas,
        num_threads
    );

    // Step 3. Write the translated peptides set to a temporary file
    let dir: String = if temp_dir.is_empty() {
        env::var("TMPDIR").unwrap_or_else(|_| env::temp_dir().to_string_lossy().to_string())
    } else {
        temp_dir.to_string()
    };
    let dir_path = Path::new(&dir);
    if !dir_path.exists() {
        panic!("Directory does not exist: {}", dir);
    }
    let temp_file = NamedTempFile::new_in(dir_path).unwrap();
    let persisted_file_path = temp_file.into_temp_path().keep().unwrap();
    let file = File::create(&persisted_file_path).unwrap();
    let mut writer = IpcWriter::new(file);
    writer.finish(&mut translation_set.to_dataframe(strategy.as_str())).unwrap();
    drop(writer);
    translation_set.translations.clear();
    translation_set.translations.shrink_to_fit();
    drop(translation_set);
    Ok(persisted_file_path.to_str().unwrap().to_string())
}


#[pymodule]
fn exactolibrs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(identify_dna_variants, m)?)?;
    m.add_function(wrap_pyfunction!(identify_case_specific_dna_variants, m)?)?;
    m.add_function(wrap_pyfunction!(identify_rna_variants, m)?)?;
//     m.add_function(wrap_pyfunction!(identify_peptide_variants, m)?)?;
    m.add_function(wrap_pyfunction!(translate_rna_sequence, m)?)?;
    m.add_function(wrap_pyfunction!(translate_rna_fasta_file, m)?)?;
    m.add_function(wrap_pyfunction!(translate_rna_fastq_file, m)?)?;
    Ok(())
}
