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
use std::collections::HashMap;
use crate::variant_annotation::VariantAnnotation;


#[derive(Debug, Serialize, Deserialize)]
pub struct VariantCall {
    pub id: String,
    pub source_id: String,
    pub sample_id: String,
    pub phase_block_id: String,
    pub clone_id: String,
    pub nucleic_acid: String,
    pub variant_calling_method: String,
    pub sequencing_platform: String,
    pub chromosome_1: String,
    pub position_1: isize,
    pub chromosome_2: String,
    pub position_2: isize,
    pub reference_allele: String,
    pub alternate_allele: String,
    pub filter: String,
    pub quality_score: f64,
    pub precise: bool,
    pub variant_type: String,
    pub variant_subtype: String,
    pub variant_size: isize,
    pub total_read_count: isize,
    pub reference_allele_read_count: isize,
    pub alternate_allele_read_count: isize,
    pub alternate_allele_fraction: f64,
    pub alternate_allele_read_ids: Vec<String>,
    pub variant_sequences: Vec<String>,
    pub tool_attributes: HashMap<String, String>,
    pub position_1_annotations: Vec<VariantAnnotation>,
    pub position_2_annotations: Vec<VariantAnnotation>,
    pub tags: Vec<String>
}

impl VariantCall {
    pub fn new(
        id: String,
        source_id: String,
        sample_id: String,
        phase_block_id: String,
        clone_id: String,
        nucleic_acid: String,
        variant_calling_method: String,
        sequencing_platform: String,
        chromosome_1: String,
        position_1: isize,
        chromosome_2: String,
        position_2: isize,
        reference_allele: String,
        alternate_allele: String,
        filter: String,
        quality_score: f64,
        precise: bool,
        variant_type: String,
        variant_subtype: String,
        variant_size: isize,
        total_read_count: isize,
        reference_allele_read_count: isize,
        alternate_allele_read_count: isize,
        alternate_allele_fraction: f64) -> Self {
        Self {
            id: id,
            source_id: source_id,
            sample_id: sample_id,
            phase_block_id: phase_block_id,
            clone_id: clone_id,
            nucleic_acid: nucleic_acid,
            variant_calling_method: variant_calling_method,
            sequencing_platform: sequencing_platform,
            chromosome_1: chromosome_1,
            position_1: position_1,
            chromosome_2: chromosome_2,
            position_2: position_2,
            reference_allele: reference_allele,
            alternate_allele: alternate_allele,
            filter: filter,
            quality_score: quality_score,
            precise: precise,
            variant_type: variant_type,
            variant_subtype: variant_subtype,
            variant_size: variant_size,
            total_read_count: total_read_count,
            reference_allele_read_count: reference_allele_read_count,
            alternate_allele_read_count: alternate_allele_read_count,
            alternate_allele_fraction: alternate_allele_fraction,
            variant_sequences: Vec::new(),
            alternate_allele_read_ids: Vec::new(),
            tool_attributes: HashMap::new(),
            position_1_annotations: Vec::new(),
            position_2_annotations: Vec::new(),
            tags: Vec::new()
        }
    }

    pub fn add_alternate_allele_read_id(&mut self, alternate_allele_read_id: String) {
        self.alternate_allele_read_ids.push(alternate_allele_read_id);
    }

    pub fn add_position_1_annotation(&mut self, variant_annotation: VariantAnnotation) {
        self.position_1_annotations.push(variant_annotation);
    }

    pub fn add_position_2_annotation(&mut self, variant_annotation: VariantAnnotation) {
        self.position_2_annotations.push(variant_annotation);
    }

    pub fn add_tag(&mut self, tag: String) {
        self.tags.push(tag);
    }

    pub fn add_tool_attribute(&mut self, key: String, value: String) {
        self.tool_attributes.insert(key, value);
    }

    pub fn add_variant_sequence(&mut self, variant_sequence: String) {
        self.variant_sequences.push(variant_sequence);
    }
}

impl Clone for VariantCall {
    fn clone(&self) -> Self {
        let mut tool_attributes: HashMap<String, String> = HashMap::new();
        for (key, value) in self.tool_attributes.iter() {
            tool_attributes.insert(key.clone(), value.clone());
        }
        let mut position_1_annotations: Vec<VariantAnnotation> = Vec::new();
        for variant_annotation in &self.position_1_annotations {
            position_1_annotations.push(variant_annotation.clone());
        }
        let mut position_2_annotations: Vec<VariantAnnotation> = Vec::new();
        for variant_annotation in &self.position_2_annotations {
            position_2_annotations.push(variant_annotation.clone());
        }
        VariantCall {
            id: self.id.clone(),
            source_id: self.source_id.clone(),
            sample_id: self.sample_id.clone(),
            phase_block_id: self.phase_block_id.clone(),
            clone_id: self.clone_id.clone(),
            nucleic_acid: self.nucleic_acid.clone(),
            variant_calling_method: self.variant_calling_method.clone(),
            sequencing_platform: self.sequencing_platform.clone(),
            chromosome_1: self.chromosome_1.clone(),
            position_1: self.position_1,
            chromosome_2: self.chromosome_2.clone(),
            position_2: self.position_2,
            reference_allele: self.reference_allele.clone(),
            alternate_allele: self.alternate_allele.clone(),
            filter: self.filter.clone(),
            quality_score: self.quality_score,
            precise: self.precise,
            variant_type: self.variant_type.clone(),
            variant_subtype: self.variant_subtype.clone(),
            variant_size: self.variant_size,
            total_read_count: self.total_read_count,
            reference_allele_read_count: self.reference_allele_read_count,
            alternate_allele_read_count: self.alternate_allele_read_count,
            alternate_allele_fraction: self.alternate_allele_fraction,
            alternate_allele_read_ids: self.alternate_allele_read_ids.clone(),
            variant_sequences: self.variant_sequences.clone(),
            tool_attributes: tool_attributes,
            position_1_annotations: position_1_annotations,
            position_2_annotations: position_2_annotations,
            tags: self.tags.clone()
        }
    }
}