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
extern crate noodles_fasta;
extern crate polars;
extern crate pyo3;
extern crate rayon;

use exacto::core::prelude as core;
use exacto::graph::prelude as graph;
use noodles_fasta::{self as fasta, record::{Definition, Record, Sequence}};
use noodles_fasta::io::Writer;
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use rayon::prelude::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::BufWriter;


/// This function builds a genome variation graph.
#[pyfunction]
pub fn build_genome_variation_graph(
    py: Python,
    df_variants: PyDataFrame,
    fasta_file: String,
    output_fasta_file: String,
    sequence_prefix: String,
    remove_unknown_bases: bool,
    only_variant_sequences: bool,
    graph_type: String,
    num_threads: usize,
    output_type: String,
    verbose: bool
) -> PyResult<PyObject> {
    core::init_logging(verbose);

    let mut df_variants_: DataFrame = df_variants.into();

    let graph_type: graph::VarGraphTypes = graph_type.as_str().parse().unwrap();

    // Step 1. Build a variation graph with the variants
    let vargraphs: Vec<graph::VarGraph> = graph::build_genome_variation_graph(
        &fasta_file,
        &df_variants_,
        graph_type,
        num_threads
    );

    // Step 2. Find paths
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let mut paths_list: Vec<HashSet<graph::VarGraphPath>> = thread_pool.install(|| {
        vargraphs
            .par_iter()
            .map(|vargraph| {
                vargraph.find_genome_paths(&vargraph.get_variant_node_ids().into_iter().collect(), &HashSet::new())
            })
            .collect()
    });

    match output_type.as_str() {
        "dataframe" => {
            let mut rows: Vec<(String, String, bool)> = Vec::new();
            let mut idx: usize = 1;

            // Variant sequences
            for paths in paths_list.iter() {
                for path in paths.iter() {
                    let sequence = path.get_sequence();
                    if remove_unknown_bases {
                        for sub_sequence in sequence.split(|c| c == 'N' || c == 'n').filter(|s| !s.is_empty()) {
                            let name = format!("{sequence_prefix}_{idx}");
                            rows.push((name.to_string(), sub_sequence.to_string(), true));
                            idx += 1;
                        }
                    } else {
                        let name: String = format!("{sequence_prefix}_{idx}");
                        rows.push((name.to_string(), sequence.to_string(), true));
                        idx += 1;
                    }
                }
            }

            if only_variant_sequences == false {
                // Reference chromosomes without variants
                let mut included_chromosomes: HashSet<Box<str>> = HashSet::new();
                let col_chromosome_1 = df_variants_.column("chromosome_1").unwrap().str().unwrap();
                let col_chromosome_2 = df_variants_.column("chromosome_2").unwrap().str().unwrap();
                for i in 0..df_variants_.height() {
                    let chromosome_1: Box<str> = col_chromosome_1.get(i).unwrap().into();
                    let chromosome_2: Box<str> = col_chromosome_2.get(i).unwrap().into();
                    included_chromosomes.insert(chromosome_1);
                    included_chromosomes.insert(chromosome_2);
                }
                let fasta_sequence_ids: Vec<(Box<str>, u32)> = core::get_fasta_sequence_ids(fasta_file.as_str());
                for (sequence_id, length) in fasta_sequence_ids.iter() {
                    if included_chromosomes.contains(sequence_id) == false {
                        let sequence: Box<str> = core::get_fasta_sequence(&*sequence_id, 1, *length, fasta_file.as_str());
                        if remove_unknown_bases {
                            for sub_sequence in sequence
                                .split(|c| c == 'N' || c == 'n')
                                .filter(|s| !s.is_empty()) {
                                let name: String = format!("{sequence_prefix}_{idx}");
                                rows.push((name.to_string(), sub_sequence.to_string(), false));
                                idx += 1;
                            }
                        } else {
                            let name: String = format!("{sequence_prefix}_{idx}");
                            rows.push((name.to_string(), sequence.to_string(), false));
                            idx += 1;
                        }
                    }
                }
            }

            let df = DataFrame::new(vec![
                Column::from(Series::new("id".into(), rows.iter().map(|r| r.0.clone()).collect::<Vec<_>>())),
                Column::from(Series::new("sequence".into(), rows.iter().map(|r| r.1.clone()).collect::<Vec<_>>())),
                Column::from(Series::new("is_variant".into(), rows.iter().map(|r| r.2).collect::<Vec<_>>()))
            ]).unwrap();

            Ok(PyDataFrame(df).into_py(py))
        },
        "file" => {
            // Write the variant sequences to FASTA file
            let f = File::create(output_fasta_file).unwrap();
            let mut writer = Writer::new(BufWriter::new(f));
            let mut idx: usize = 1;
            for paths in paths_list.iter() {
                for path in paths.iter() {
                    let sequence = path.get_sequence();
                    if remove_unknown_bases {
                        for sub_sequence in sequence.split(|c| c == 'N' || c == 'n').filter(|s| !s.is_empty()) {
                            let name = format!("{sequence_prefix}_{idx}");
                            let def = Definition::new(name, None);
                            let sequence_ = Sequence::from(sub_sequence.as_bytes().to_vec());
                            let record = Record::new(def, sequence_);
                            writer.write_record(&record).unwrap();
                            idx += 1;
                        }
                    } else {
                        let name: String = format!("{sequence_prefix}_{idx}");
                        let def = Definition::new(name, None);
                        let sequence_ = Sequence::from(sequence.as_bytes().to_vec());
                        let record = Record::new(def, sequence_);
                        writer.write_record(&record).unwrap();
                        idx += 1;
                    }
                }
            }

            // Write reference chromosomes without variants
            let mut included_chromosomes: HashSet<Box<str>> = HashSet::new();
            let col_chromosome_1 = df_variants_.column("chromosome_1").unwrap().str().unwrap();
            let col_chromosome_2 = df_variants_.column("chromosome_2").unwrap().str().unwrap();
            for i in 0..df_variants_.height() {
                let chromosome_1: Box<str> = col_chromosome_1.get(i).unwrap().into();
                let chromosome_2: Box<str> = col_chromosome_2.get(i).unwrap().into();
                included_chromosomes.insert(chromosome_1);
                included_chromosomes.insert(chromosome_2);
            }
            let fasta_sequence_ids: Vec<(Box<str>, u32)> = core::get_fasta_sequence_ids(fasta_file.as_str());
            for (sequence_id, length) in fasta_sequence_ids.iter() {
                if included_chromosomes.contains(sequence_id) == false {
                    let sequence: Box<str> = core::get_fasta_sequence(&*sequence_id, 1, *length, fasta_file.as_str());
                    if remove_unknown_bases {
                        for sub_sequence in sequence
                            .split(|c| c == 'N' || c == 'n')
                            .filter(|s| !s.is_empty()) {
                            let name: String = format!("{sequence_prefix}_{idx}");
                            let def = Definition::new(name, None);
                            let sequence_ = Sequence::from(sub_sequence.as_bytes().to_vec());
                            let record = Record::new(def, sequence_);
                            writer.write_record(&record).unwrap();
                            idx += 1;
                        }
                    } else {
                        let name: String = format!("{sequence_prefix}_{idx}");
                        let def = Definition::new(name, None);
                        let sequence_ = Sequence::from(sequence.as_bytes().to_vec());
                        let record = Record::new(def, sequence_);
                        writer.write_record(&record).unwrap();
                        idx += 1;
                    }
                }
            }

            Ok(py.None().into_py(py))
        },
        "vector" => {
            // Return: Vec<(id, is_variant, sequence)>
            let mut rows: Vec<(String, bool, String)> = Vec::new();
            let mut idx: usize = 1;

            // Variant sequences
            for paths in paths_list.iter() {
                for path in paths.iter() {
                    let sequence = path.get_sequence();
                    if remove_unknown_bases {
                        for sub_sequence in sequence
                            .split(|c| c == 'N' || c == 'n')
                            .filter(|s| !s.is_empty())
                        {
                            let name = format!("{sequence_prefix}_{idx}");
                            rows.push((name, true, sub_sequence.to_string()));
                            idx += 1;
                        }
                    } else {
                        let name = format!("{sequence_prefix}_{idx}");
                        rows.push((name, true, sequence.to_string()));
                        idx += 1;
                    }
                }
            }

            if only_variant_sequences == false {
                // Reference chromosomes without variants
                let mut included_chromosomes: HashSet<Box<str>> = HashSet::new();
                let col_chromosome_1 = df_variants_.column("chromosome_1").unwrap().str().unwrap();
                let col_chromosome_2 = df_variants_.column("chromosome_2").unwrap().str().unwrap();
                for i in 0..df_variants_.height() {
                    let chromosome_1: Box<str> = col_chromosome_1.get(i).unwrap().into();
                    let chromosome_2: Box<str> = col_chromosome_2.get(i).unwrap().into();
                    included_chromosomes.insert(chromosome_1);
                    included_chromosomes.insert(chromosome_2);
                }

                let fasta_sequence_ids: Vec<(Box<str>, u32)> =
                    core::get_fasta_sequence_ids(fasta_file.as_str());

                for (sequence_id, length) in fasta_sequence_ids.iter() {
                    if !included_chromosomes.contains(sequence_id) {
                        let sequence: Box<str> = core::get_fasta_sequence(
                            &*sequence_id,
                            1,
                            *length,
                            fasta_file.as_str(),
                        );

                        if remove_unknown_bases {
                            for sub_sequence in sequence
                                .split(|c| c == 'N' || c == 'n')
                                .filter(|s| !s.is_empty())
                            {
                                let name = format!("{sequence_prefix}_{idx}");
                                rows.push((name, false, sub_sequence.to_string()));
                                idx += 1;
                            }
                        } else {
                            let name = format!("{sequence_prefix}_{idx}");
                            rows.push((name, false, sequence.to_string()));
                            idx += 1;
                        }
                    }
                }
            }

            Ok(rows.into_py(py))
        },
        other => {
            let error_message = format!("Unsupported value for output_type: {}", other);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message))
        }
    }
}