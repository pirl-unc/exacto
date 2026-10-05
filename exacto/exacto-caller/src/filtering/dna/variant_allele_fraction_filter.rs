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


use crate::prelude::*;


pub struct VariantAlleleFractionFilter {
    pub min_alt_allele_fraction: f64
}

impl VariantAlleleFractionFilter {
    pub fn new(min_alt_allele_fraction: f64) -> Self {
        Self {
            min_alt_allele_fraction: min_alt_allele_fraction
        }
    }
}

impl VariantFilter for VariantAlleleFractionFilter {
    type Input = VariantCall;

    fn passes(&self, variant_call: &VariantCall) -> bool {
        if variant_call.get_alternate_allele_fraction() < self.min_alt_allele_fraction {
            false
        } else {
            true
        }
    }
}
