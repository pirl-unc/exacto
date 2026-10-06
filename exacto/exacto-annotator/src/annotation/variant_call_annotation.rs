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


use exacto_caller::prelude::*;
use serde::{Serialize, Deserialize};
use exacto_core::prelude::*;

use crate::prelude::*;


#[derive(Debug, Serialize, Deserialize)]
pub struct VariantCallAnnotation {
    pub id: VariantID,
    pub chromosome_1: ReferenceChromosomeName,
    pub position_1: ReferencePosition,
    pub chromosome_2: ReferenceChromosomeName,
    pub position_2: ReferencePosition,
    pub variant_type: VariantType,
    pub variant_sequence: Box<str>,
    pub position_1_annotation: PositionAnnotation,
    pub position_2_annotation: PositionAnnotation,

    /// Annotation of `position_1 + 1`. The first altered base of an SNV, an MNV or a deletion.
    pub position_1_plus_1_annotation: PositionAnnotation,

    /// Annotation of `position_2 - 1`. The last altered base of an SNV, an MNV or a deletion.
    pub position_2_minus_1_annotation: PositionAnnotation
}

impl PartialEq for VariantCallAnnotation {
    fn eq(&self, other: &Self) -> bool {
        if self.id == other.id &&
            self.chromosome_1 == other.chromosome_1 &&
            self.position_1 == other.position_1 &&
            self.chromosome_2 == other.chromosome_2 &&
            self.position_2 == other.position_2 &&
            self.variant_type == other.variant_type &&
            self.variant_sequence == other.variant_sequence &&
            self.position_1_annotation == other.position_1_annotation &&
            self.position_1_plus_1_annotation == other.position_1_plus_1_annotation &&
            self.position_2_minus_1_annotation == other.position_2_minus_1_annotation &&
            self.position_2_annotation == other.position_2_annotation {
                true
        } else {
            false
        }
    }
}

impl VariantCallAnnotation {
    pub fn new(
        id: VariantID,
        chromosome_1: ReferenceChromosomeName,
        position_1: ReferencePosition,
        chromosome_2: ReferenceChromosomeName,
        position_2: ReferencePosition,
        variant_type: VariantType,
        variant_sequence: Box<str>,
        position_1_annotation: PositionAnnotation,
        position_1_plus_1_annotation: PositionAnnotation,
        position_2_minus_1_annotation: PositionAnnotation,
        position_2_annotation: PositionAnnotation
    ) -> Self {
        VariantCallAnnotation {
            id: id,
            chromosome_1: chromosome_1,
            position_1: position_1,
            chromosome_2: chromosome_2,
            position_2: position_2,
            variant_type: variant_type,
            variant_sequence: variant_sequence,
            position_1_annotation,
            position_1_plus_1_annotation,
            position_2_minus_1_annotation,
            position_2_annotation
        }
    }
}

impl Clone for VariantCallAnnotation {
    fn clone(&self) -> Self {
        VariantCallAnnotation {
            id: self.id,
            chromosome_1: self.chromosome_1.clone(),
            position_1: self.position_1,
            chromosome_2: self.chromosome_2.clone(),
            position_2: self.position_2,
            variant_type: self.variant_type.clone(),
            variant_sequence: self.variant_sequence.clone(),
            position_1_annotation: self.position_1_annotation.clone(),
            position_1_plus_1_annotation: self.position_1_plus_1_annotation.clone(),
            position_2_minus_1_annotation: self.position_2_minus_1_annotation.clone(),
            position_2_annotation: self.position_2_annotation.clone()
        }
    }
}


#[cfg(test)]
#[path = "../tests/annotation/variant_call_annotation.rs"]
mod tests;