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


use exacto_core::prelude::*;

use crate::prelude::*;


pub struct RNAJunctionReadSupportFilter<'a, C: SplicingEventReadCounts> {
    read_counts: &'a C,
    read_support_index: &'a dyn ReadSupportIndex<Input = ReadDepth>
}

impl <'a, C: SplicingEventReadCounts> RNAJunctionReadSupportFilter<'a, C> {
    pub fn new(
        read_counts: &'a C,
        read_support_index: &'a dyn ReadSupportIndex<Input = ReadDepth>
    ) -> Self {
        Self {
            read_counts: read_counts,
            read_support_index: read_support_index
        }
    }
}

impl <C: SplicingEventReadCounts> ReadFilter for RNAJunctionReadSupportFilter<'_, C> {
    type Input = C::Input;

    fn passes(&self, input: &C::Input) -> bool {
        self.read_counts
            .get_read_counts(input)
            .into_iter()
            .all(|(num_reads, depth)| num_reads >= self.read_support_index.get_min_read_support(depth))
    }
}


#[cfg(test)]
#[path = "../../tests/filtering/rna/junction_read_support_filter.rs"]
mod tests;
