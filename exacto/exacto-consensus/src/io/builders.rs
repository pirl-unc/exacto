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


use exacto_core::prelude::LIST_SEPARATOR;

use crate::prelude::*;


pub fn build_consensus_sequence_records<'a>(
    consensus_sequence_set: &'a ConsensusSequenceSet
) -> impl Iterator<Item = ConsensusSequenceRecord> + 'a {
    assert!(
        !consensus_sequence_set.sequences.is_empty(),
        "consensus_sequence_set.sequences is empty."
    );
    consensus_sequence_set.sequences.iter().map(move |consensus| {
        // Sorted, not just collected: `read_names` is a `HashSet`, so emitting it in iteration
        // order would make consecutive runs over identical input produce different files.
        let mut read_names: Vec<&str> = consensus.read_names.iter().map(|name| name.as_ref()).collect();
        read_names.sort_unstable();
        ConsensusSequenceRecord {
            cluster_id: consensus.cluster_id,
            consensus_sequence: consensus.consensus_sequence.clone(),
            num_reads: consensus.read_names.len(),
            read_names: read_names.join(LIST_SEPARATOR).into_boxed_str()
        }
    })
}


#[cfg(test)]
#[path = "../tests/io/builders.rs"]
mod tests;