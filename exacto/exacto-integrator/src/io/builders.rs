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


use exacto_core::prelude::join_ids;

use crate::prelude::*;


pub fn build_integrated_variant_records<'a>(
    integrated_variants: &'a Vec<IntegratedVariant>
) -> impl Iterator<Item = IntegratedVariantRecord> + 'a {
    let mut integrated_variant_records: Vec<IntegratedVariantRecord> = Vec::new();
    for integrated_variant in integrated_variants {
        for (dna_variant_id, integrated_variant_distance) in integrated_variant.dna_variant_ids.iter() {
            let reference_transcript_ids_string: String = integrated_variant.reference_transcript_ids.iter()
                .map(|s| s.as_ref()) // &Box<str> -> &str
                .collect::<Vec<&str>>()
                .join(",");
            integrated_variant_records.push(
                IntegratedVariantRecord {
                    assembled_transcript_name: integrated_variant.assembled_transcript_name.clone(),
                    transcript_model_id: integrated_variant.transcript_model_id,
                    reference_gene_name: join_ids(integrated_variant.reference_gene_names.clone()).into(),
                    reference_transcript_id: join_ids(integrated_variant.reference_transcript_ids.clone()).into(),
                    rna_variant_id: integrated_variant.rna_variant_id,
                    dna_variant_id: *dna_variant_id,
                    distance: integrated_variant_distance.distance,
                    rna_variant_position: integrated_variant_distance.rna_variant_position_used.as_str().into(),
                    dna_variant_position: integrated_variant_distance.dna_variant_position_used.as_str().into()
                }
            )
        }
    }
    
    integrated_variant_records.into_iter()
}
