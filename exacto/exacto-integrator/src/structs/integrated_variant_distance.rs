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


use std::hash::{Hash, Hasher};
use serde::{Serialize, Deserialize};

use crate::prelude::*;


#[derive(Debug, Eq, PartialEq, Serialize, Deserialize)]
pub struct IntegratedVariantDistance {
    pub distance: usize,
    pub rna_variant_position_used: VariantPosition,
    pub dna_variant_position_used: VariantPosition
}

impl Hash for IntegratedVariantDistance {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.distance.hash(state);
        self.rna_variant_position_used.as_str().hash(state);
        self.dna_variant_position_used.as_str().hash(state);
    }
}

impl IntegratedVariantDistance {
    pub fn new(distance: usize, rna_variant_position_used: VariantPosition, dna_variant_position_used: VariantPosition) -> Self {
        Self {
            distance: distance,
            rna_variant_position_used: rna_variant_position_used,
            dna_variant_position_used: dna_variant_position_used
        }
    }
}

impl Clone for IntegratedVariantDistance {
    fn clone(&self) -> Self {
        IntegratedVariantDistance {
            distance: self.distance,
            rna_variant_position_used: self.rna_variant_position_used.clone(),
            dna_variant_position_used: self.dna_variant_position_used.clone()
        }
    }
}
