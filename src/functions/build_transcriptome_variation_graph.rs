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
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufWriter;
use std::sync::Mutex;


/// This function builds a genome variation graph.
#[pyfunction]
pub fn build_transcriptome_variation_graph(
    py: Python,
    df_transcript_structures: PyDataFrame,
    fasta_file: String,
    output_fasta_file: String,
    graph_type: String,
    num_threads: usize,
    batch_size: usize,
    output_type: String,
    verbose: bool
) -> PyResult<PyObject> {
    core::init_logging(verbose);

    let df_transcript_structures_: DataFrame = df_transcript_structures.into();

    let graph_type: graph::VarGraphTypes = graph_type.as_str().parse().unwrap();

    // Step 1. Prepare once: chromosome lengths, partition groups, thread pool
    let mut chromosome_lengths_map: HashMap<Box<str>, u32> = HashMap::new();
    for (chromosome, length) in core::get_fasta_sequence_ids(&fasta_file).iter() {
        chromosome_lengths_map.insert(chromosome.clone(), *length);
    }
    let groups: Vec<DataFrame> = df_transcript_structures_
        .partition_by(["transcript_model_id", "reference_transcript_ids"], true)
        .unwrap();
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let effective_batch_size = if batch_size == 0 { groups.len() } else { batch_size };

    // Step 2. Process in batches — each thread builds graph + finds path + drops graph
    match output_type.as_str() {
        "file" => {
            // Stream to disk — zero accumulation (mutex protects the writer)
            let writer = Mutex::new(
                Writer::new(BufWriter::new(File::create(&output_fasta_file).unwrap()))
            );
            for chunk in groups.chunks(effective_batch_size) {
                thread_pool.install(|| {
                    chunk.par_iter().for_each(|df_group| {
                        let (transcript_model_id, reference_transcript_ids, vargraph) =
                            graph::build_single_transcriptome_variation_graph(
                                &fasta_file, df_group, &chromosome_lengths_map, graph_type.clone()
                            );
                        let path = vargraph.find_transcript_path(
                            transcript_model_id,
                            &*reference_transcript_ids,
                            &vargraph.get_variant_node_ids().into_iter().collect()
                        );
                        let sequence = path.get_sequence();
                        let name: String = format!("{transcript_model_id}");
                        let def = Definition::new(name, None);
                        let sequence_ = Sequence::from(sequence.as_bytes().to_vec());
                        let record = Record::new(def, sequence_);
                        writer.lock().unwrap().write_record(&record).unwrap();
                        // lock released, graph + path dropped
                    });
                });
            }
            Ok(py.None().into_py(py))
        },
        output_kind @ ("dataframe" | "vector") => {
            // Accumulate lightweight (id, sequence) pairs only
            let mut all_rows: Vec<(usize, String)> = Vec::new();
            for chunk in groups.chunks(effective_batch_size) {
                let batch_results: Vec<(usize, String)> = thread_pool.install(|| {
                    chunk
                        .par_iter()
                        .map(|df_group| {
                            let (transcript_model_id, reference_transcript_ids, vargraph) =
                                graph::build_single_transcriptome_variation_graph(
                                    &fasta_file, df_group, &chromosome_lengths_map, graph_type.clone()
                                );
                            let path = vargraph.find_transcript_path(
                                transcript_model_id,
                                &*reference_transcript_ids,
                                &vargraph.get_variant_node_ids().into_iter().collect()
                            );
                            (transcript_model_id, path.get_sequence().to_string())
                            // graph + path dropped
                        })
                        .collect()
                });
                all_rows.extend(batch_results);
            }
            if output_kind == "dataframe" {
                let df = DataFrame::new(vec![
                    Column::from(Series::new("transcript_model_id".into(), all_rows.iter().map(|r| r.0 as u64).collect::<Vec<_>>())),
                    Column::from(Series::new("sequence".into(), all_rows.iter().map(|r| r.1.clone()).collect::<Vec<_>>())),
                ]).unwrap();
                Ok(PyDataFrame(df).into_py(py))
            } else {
                Ok(all_rows.into_py(py))
            }
        },
        other => {
            let error_message = format!("Unsupported value for output_type: {}", other);
            Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message))
        }
    }
}
