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


extern crate pyo3;

use pyo3::prelude::*;

mod functions;

use functions::annotate_variant_calls::*;
use functions::identify_case_specific_dna_variants::*;
use functions::identify_dna_variants::*;
use functions::identify_rna_variants::*;
use functions::integrate_dna_rna_variants::*;
use functions::remove_unspliced_rnas::*;
use functions::translate_fasta_file::*;
use functions::translate_fastq_file::*;
use functions::translate_sequence::*;
use functions::translate_structures::*;


// extern crate chrono;
// extern crate env_logger;
// extern crate exacto;
// extern crate exitcode;
// extern crate flate2;
// extern crate log;
// extern crate noodles_bgzf;
// extern crate noodles_fasta;
// extern crate noodles_fastq;
// extern crate polars;
// extern crate pyo3_polars;
// extern crate sysinfo;
// extern crate tempfile;
//
// use chrono::Local;
// use env_logger::{Builder, Env};
// use exacto::annotator::prelude as annotator;
// use exacto::caller::prelude as caller;
// use exacto::core::prelude as core;
// use exacto::graph::prelude as graph;
// use exacto::integrator::prelude as integrator;
// use exacto::translator::prelude as translator;
// use flate2::read::GzDecoder;
// use log::{info, LevelFilter};
// use polars::prelude::*;
// use polars::io::ipc::IpcWriter;
// use pyo3::types::{PyDict, PyList};
// use pyo3_polars::PyDataFrame;
// use noodles_bgzf as bgzf;
// use noodles_fasta::{self as fasta, record::{Definition, Record, Sequence}};
// use noodles_fasta::io::Writer;
// use noodles_fastq as fastq;
// use std::collections::HashMap;
// use std::env;
// use std::fs::File;
// use std::io::{BufReader,BufWriter,Cursor,Read,Write};
// use std::path::Path;
// use sysinfo::System;
// use tempfile::NamedTempFile;



#[pymodule]
fn exactolibrs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(annotate_variant_calls, m)?)?;
    m.add_function(wrap_pyfunction!(identify_case_specific_dna_variants, m)?)?;
    m.add_function(wrap_pyfunction!(identify_dna_variants, m)?)?;
    m.add_function(wrap_pyfunction!(identify_rna_variants, m)?)?;
    m.add_function(wrap_pyfunction!(integrate_dna_rna_variants, m)?)?;
    m.add_function(wrap_pyfunction!(remove_unspliced_rnas, m)?)?;
    m.add_function(wrap_pyfunction!(translate_fasta_file, m)?)?;
    m.add_function(wrap_pyfunction!(translate_fastq_file, m)?)?;
    m.add_function(wrap_pyfunction!(translate_sequence, m)?)?;
    m.add_function(wrap_pyfunction!(translate_structures, m)?)?;

//     m.add_function(wrap_pyfunction!(build_variation_graph, m)?)?;
//     m.add_function(wrap_pyfunction!(integrate_dna_rna_variants, m)?)?;
//     m.add_function(wrap_pyfunction!(identify_peptide_variants, m)?)?;
    Ok(())
}
