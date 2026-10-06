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


/// This function builds a transcriptome variation graph from the transcript model alignments table
/// of call-rna-transcript-vars and writes one sequence per transcript model, named by its
/// `assembled_transcript_name`, in the order of the table. A model whose rows do not join into one
/// path is logged and skipped.
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

    // (assembled_transcript_name, sequence), when they are returned rather than written
    let mut rows: Vec<(String, String)> = Vec::new();
    graph::find_transcriptome_variation_graph_sequences(
        &fasta_file,
        &df_transcript_structures_,
        graph_type,
        num_threads,
        batch_size,
        |name, sequence| {
            match writer.as_mut() {
                Some(writer) => {
                    let record = Record::new(
                        Definition::new(name.to_string(), None),
                        Sequence::from(sequence.as_bytes().to_vec())
                    );
                    writer.write_record(&record).map_err(|error| graph::GraphError::File {
                        file: output_fasta_file.as_str().into(),
                        reason: error.to_string().into()
                    })?;
                },
                None => rows.push((name.to_string(), sequence.to_string()))
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
                Column::from(Series::new("assembled_transcript_name".into(), rows.iter().map(|r| r.0.clone()).collect::<Vec<_>>())),
                Column::from(Series::new("sequence".into(), rows.iter().map(|r| r.1.clone()).collect::<Vec<_>>())),
            ]).unwrap();
            Ok(PyDataFrame(df).into_py(py))
        },
        None => Ok(rows.into_py(py))
    }
}
