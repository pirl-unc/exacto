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

use crate::io::records::*;


pub fn integrated_variant_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = IntegratedVariantRecord>
{
    let mut assembled_transcript_name: Vec<String> = Vec::new();
    let mut transcript_model_id: Vec<u32> = Vec::new();
    let mut reference_gene_name: Vec<String> = Vec::new();
    let mut reference_transcript_id: Vec<String> = Vec::new();
    let mut rna_variant_id: Vec<u32> = Vec::new();
    let mut dna_variant_id: Vec<u32> = Vec::new();
    let mut distance: Vec<u32> = Vec::new();
    let mut rna_variant_position: Vec<String> = Vec::new();
    let mut dna_variant_position: Vec<String> = Vec::new();

    for r in records {
        assembled_transcript_name.push(r.assembled_transcript_name.to_string());
        transcript_model_id.push(r.transcript_model_id);
        reference_gene_name.push(r.reference_gene_name.to_string());
        reference_transcript_id.push(r.reference_transcript_id.to_string());
        rna_variant_id.push(r.rna_variant_id);
        dna_variant_id.push(r.dna_variant_id);
        distance.push(r.distance);
        rna_variant_position.push(r.rna_variant_position.to_string());
        dna_variant_position.push(r.dna_variant_position.to_string());
    }

    DataFrame::new(vec![
        Column::from(Series::new("assembled_transcript_name".into(), assembled_transcript_name)),
        Column::from(Series::new("transcript_model_id".into(), transcript_model_id)),
        Column::from(Series::new("reference_gene_name".into(), reference_gene_name)),
        Column::from(Series::new("reference_transcript_id".into(), reference_transcript_id)),
        Column::from(Series::new("rna_variant_id".into(), rna_variant_id)),
        Column::from(Series::new("dna_variant_id".into(), dna_variant_id)),
        Column::from(Series::new("distance".into(), distance)),
        Column::from(Series::new("rna_variant_position".into(), rna_variant_position)),
        Column::from(Series::new("dna_variant_position".into(), dna_variant_position))
    ])
    .unwrap()
}
