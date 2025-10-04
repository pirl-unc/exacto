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


use polars::prelude::*;
use std::str::FromStr;

use crate::prelude::*;


pub fn build_variation_graph(
    fasta_file: &str,
    df_variants: &DataFrame
) -> VarGraph {
    // Step 1. Get DataFrame columns
    let variant_id_col_casted = df_variants.column("variant_id").unwrap().cast(&DataType::String).unwrap();
    let variant_id_col = variant_id_col_casted.str().unwrap();
    let chromosome_1_col = df_variants.column("chromosome_1").unwrap().str().unwrap();
    let position_1_col = df_variants.column("position_1").unwrap().i64().unwrap();
    let operation_1_col = df_variants.column("operation_1").unwrap().str().unwrap();
    let strand_1_col = df_variants.column("strand_1").unwrap().str().unwrap();
    let chromosome_2_col = df_variants.column("chromosome_2").unwrap().str().unwrap();
    let position_2_col = df_variants.column("position_2").unwrap().i64().unwrap();
    let operation_2_col = df_variants.column("operation_2").unwrap().str().unwrap();
    let strand_2_col = df_variants.column("strand_2").unwrap().str().unwrap();
    let sequence_col = df_variants.column("sequence").unwrap().str().unwrap();

    // Step 2. Build the variation graph
    let mut vargraph: VarGraph = VarGraph::from_fasta_file(fasta_file);
    for i in 0..df_variants.height() {
        let variant_id: &str = variant_id_col.get(i).unwrap();
        let chromosome_1: &str = chromosome_1_col.get(i).unwrap();
        let position_1: usize = position_1_col.get(i).unwrap() as usize;
        let orientation_1: VarGraphOrientations = VarGraphOrientations::from_str(operation_1_col.get(i).unwrap()).unwrap();
        let strand_1: VarGraphStrands = VarGraphStrands::from_str(strand_1_col.get(i).unwrap()).unwrap();
        let chromosome_2: &str = chromosome_2_col.get(i).unwrap();
        let position_2: usize = position_2_col.get(i).unwrap() as usize;
        let orientation_2: VarGraphOrientations = VarGraphOrientations::from_str(operation_2_col.get(i).unwrap()).unwrap();
        let strand_2: VarGraphStrands = VarGraphStrands::from_str(strand_2_col.get(i).unwrap()).unwrap();
        let sequence: &str = sequence_col.get(i).unwrap_or("");

        vargraph.add_variant(
            chromosome_1,
            position_1,
            orientation_1,
            strand_1,
            chromosome_2,
            position_2,
            orientation_2,
            strand_2,
            sequence
        );
    }

    vargraph
}
