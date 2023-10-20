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


extern crate serde;
use serde::{Serialize, Deserialize};


#[derive(Debug, Serialize, Deserialize)]
pub struct VariantAnnotation {
    region: String,
    source: String,
    source_version: String,
    gene_id: String,
    gene_stable_id: String,
    gene_version: String,
    gene_name: String,
    gene_type: String,
    gene_strand: String,
    species: String
}

impl VariantAnnotation {
    fn new(
        region: String,
        source: String,
        source_version: String,
        gene_id: String,
        gene_stable_id: String,
        gene_version: String,
        gene_name: String,
        gene_type: String,
        gene_strand: String,
        species: String) -> Self {
        VariantAnnotation {
            region: region,
            source: source,
            source_version: source_version,
            gene_id: gene_id,
            gene_stable_id: gene_stable_id,
            gene_version: gene_version,
            gene_name: gene_name,
            gene_type: gene_type,
            gene_strand: gene_strand,
            species: species
        }
    }
}

impl Clone for VariantAnnotation {
    fn clone(&self) -> Self {
        VariantAnnotation {
            region: self.region.clone(),
            source: self.source.clone(),
            source_version: self.source_version.clone(),
            gene_id: self.gene_id.clone(),
            gene_stable_id: self.gene_stable_id.clone(),
            gene_version: self.gene_version.clone(),
            gene_name: self.gene_name.clone(),
            gene_type: self.gene_type.clone(),
            gene_strand: self.gene_strand.clone(),
            species: self.species.clone()
        }
    }
}