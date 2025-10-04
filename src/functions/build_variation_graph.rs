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


extern crate polars;
extern crate pyo3;

use pyo3::prelude::*;


/// This function builds a variation graph.
#[pyfunction]
fn build_variation_graph(
    py: Python,
    variants_tsv_file: String,
    fasta_file: String,
    output_fasta_file: String
) {
    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    let df_variants: DataFrame = CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .try_into_reader_with_file_path(Some(variants_tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();

    let vargraph: graph::VarGraph = graph::build_variation_graph(&fasta_file, &df_variants);
    let paths: Vec<graph::VarGraphPath> = vargraph.get_linearized_contigs(vargraph.get_variant_node_ids());

    let f = File::create(output_fasta_file).unwrap();
    let mut writer = Writer::new(BufWriter::new(f));

    for (i, path) in paths.iter().enumerate() {
        let name: String = format!("variant_sequence_{}", i + 1);
        let def = Definition::new(name, None);
        let sequence = Sequence::from(path.get_sequence().as_bytes().to_vec());
        let record = Record::new(def, sequence);
        writer.write_record(&record).unwrap();
    }
}