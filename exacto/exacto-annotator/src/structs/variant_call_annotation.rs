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


use exacto_caller::prelude::VariantType;
use serde::{Serialize, Deserialize};

use crate::prelude::*;


#[derive(Debug, Serialize, Deserialize)]
pub struct VariantCallAnnotation {
    pub variant_call_id: usize,
    pub chromosome_1: Box<str>,
    pub position_1: u32,
    pub chromosome_2: Box<str>,
    pub position_2: u32,
    pub variant_type: VariantType,
    pub variant_sequence: Box<str>,
    pub position_1_annotation: PositionAnnotation,
    pub position_2_annotation: PositionAnnotation
}

impl VariantCallAnnotation {
    pub fn new(
        variant_call_id: usize,
        chromosome_1: &str,
        position_1: u32,
        chromosome_2: &str,
        position_2: u32,
        variant_type: VariantType,
        variant_sequence: &str,
        position_1_annotation: PositionAnnotation,
        position_2_annotation: PositionAnnotation
    ) -> Self {
        VariantCallAnnotation {
            variant_call_id,
            chromosome_1: chromosome_1.into(),
            position_1: position_1,
            chromosome_2: chromosome_2.into(),
            position_2: position_2,
            variant_type: variant_type,
            variant_sequence: variant_sequence.into(),
            position_1_annotation,
            position_2_annotation
        }
    }
}

impl Clone for VariantCallAnnotation {
    fn clone(&self) -> Self {
        VariantCallAnnotation {
            variant_call_id: self.variant_call_id,
            chromosome_1: self.chromosome_1.clone(),
            position_1: self.position_1,
            chromosome_2: self.chromosome_2.clone(),
            position_2: self.position_2,
            variant_type: self.variant_type.clone(),
            variant_sequence: self.variant_sequence.clone(),
            position_1_annotation: self.position_1_annotation.clone(),
            position_2_annotation: self.position_2_annotation.clone()
        }
    }
}
