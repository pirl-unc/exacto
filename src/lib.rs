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
use functions::build_genome_variation_graph::*;
use functions::build_transcriptome_variation_graph::*;
use functions::identify_somatic_dna_variants::*;
use functions::identify_germline_dna_variants::*;
use functions::identify_rna_variants::*;
use functions::integrate_dna_rna_variants::*;
use functions::remove_unspliced_rnas::*;
use functions::translate_fasta_file::*;
use functions::translate_fastq_file::*;
use functions::translate_sequence::*;
use functions::translate_structures::*;


#[pymodule]
fn exactolibrs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(annotate_variant_calls, m)?)?;
    m.add_function(wrap_pyfunction!(build_genome_variation_graph, m)?)?;
    m.add_function(wrap_pyfunction!(build_transcriptome_variation_graph, m)?)?;
    m.add_function(wrap_pyfunction!(identify_somatic_dna_variants, m)?)?;
    m.add_function(wrap_pyfunction!(identify_germline_dna_variants, m)?)?;
    m.add_function(wrap_pyfunction!(identify_rna_variants, m)?)?;
    m.add_function(wrap_pyfunction!(integrate_dna_rna_variants, m)?)?;
    m.add_function(wrap_pyfunction!(remove_unspliced_rnas, m)?)?;
    m.add_function(wrap_pyfunction!(translate_fasta_file, m)?)?;
    m.add_function(wrap_pyfunction!(translate_fastq_file, m)?)?;
    m.add_function(wrap_pyfunction!(translate_sequence, m)?)?;
    m.add_function(wrap_pyfunction!(translate_structures, m)?)?;
    Ok(())
}
