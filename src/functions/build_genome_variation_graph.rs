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

use exacto::core::prelude as core;
use exacto::graph::prelude as graph;
use noodles_fasta::record::{Definition, Record, Sequence};
use noodles_fasta::io::Writer;
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use std::fs::File;
use std::io::{BufWriter, Write};

use super::{graph_error, io_error};


/// This function builds a genome variation graph and writes its sequences: the variant sequences,
/// then (unless `only_variant_sequences`) the reference sequence of every contig. Records are named
/// `<sequence_prefix>_<n>` in that order; with `remove_unknown_bases` each run of bases between
/// N or n is a record of its own.
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

    let df_variants_: DataFrame = df_variants.into();

    let graph_type: graph::VarGraphTypes = graph_type.as_str().parse().map_err(|_| {
        PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Unsupported value for graph_type: {}", graph_type))
    })?;
    let mut writer = match output_type.as_str() {
        "file" => Some(Writer::new(BufWriter::new(File::create(&output_fasta_file).map_err(io_error)?))),
        "dataframe" | "vector" => None,
        other => {
            let error_message = format!("Unsupported value for output_type: {}", other);
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(error_message));
        }
    };

    // (id, sequence, is_variant) of each record, when they are returned rather than written
    let mut rows: Vec<(String, String, bool)> = Vec::new();
    let mut idx: usize = 1;
    graph::find_genome_variation_graph_sequences(
        &fasta_file,
        &df_variants_,
        graph_type,
        !only_variant_sequences,
        num_threads,
        |sequence, is_variant| {
            let pieces: Vec<&str> = if remove_unknown_bases {
                sequence.split(|c| c == 'N' || c == 'n').filter(|s| !s.is_empty()).collect()
            } else {
                vec![&*sequence]
            };
            for piece in pieces {
                let name: String = format!("{sequence_prefix}_{idx}");
                idx += 1;
                match writer.as_mut() {
                    Some(writer) => {
                        let record = Record::new(Definition::new(name, None), Sequence::from(piece.as_bytes().to_vec()));
                        writer.write_record(&record).map_err(|error| graph::GraphError::File {
                            file: output_fasta_file.as_str().into(),
                            reason: error.to_string().into()
                        })?;
                    },
                    None => rows.push((name, piece.to_string(), is_variant))
                }
            }
            Ok(())
        }
    ).map_err(graph_error)?;

    match writer {
        Some(mut writer) => {
            writer.get_mut().flush().map_err(io_error)?;
            Ok(py.None().into_py(py))
        },
        None if output_type == "dataframe" => {
            let df = DataFrame::new(vec![
                Column::from(Series::new("id".into(), rows.iter().map(|r| r.0.clone()).collect::<Vec<_>>())),
                Column::from(Series::new("sequence".into(), rows.iter().map(|r| r.1.clone()).collect::<Vec<_>>())),
                Column::from(Series::new("is_variant".into(), rows.iter().map(|r| r.2).collect::<Vec<_>>()))
            ]).unwrap();
            Ok(PyDataFrame(df).into_py(py))
        },
        None => {
            // Return: Vec<(id, is_variant, sequence)>
            let rows: Vec<(String, bool, String)> = rows.into_iter().map(|(id, sequence, is_variant)| (id, is_variant, sequence)).collect();
            Ok(rows.into_py(py))
        }
    }
}
