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


use exacto_caller::prelude::{SpliceJunction, SplicingEventReadCounts};
use exacto_core::prelude::{ReadDepth, ReadSupport};
use std::collections::{BTreeSet, HashMap, HashSet};


pub struct SpliceTransitionReadCounts {
    /// HashMap<(junction, next junction), (number of supporting reads, depth)>
    transitions: HashMap<(SpliceJunction, SpliceJunction), (ReadSupport, ReadDepth)>
}

impl SpliceTransitionReadCounts {
    pub fn new(
        chains: &HashMap<Vec<SpliceJunction>, HashSet<usize>>,
        annotated_transitions: &HashSet<(SpliceJunction, SpliceJunction)>,
        junctions_to_genes: &HashMap<(u16, u32, u32), BTreeSet<Box<str>>>
    ) -> Self {
        // Step 1. The reads spelling each transition, and the reads leaving and entering each
        // junction.
        let mut support: HashMap<(SpliceJunction, SpliceJunction), HashSet<usize>> = HashMap::new();
        let mut outgoing: HashMap<SpliceJunction, HashSet<usize>> = HashMap::new();
        let mut incoming: HashMap<SpliceJunction, HashSet<usize>> = HashMap::new();
        for (chain, reads) in chains.iter() {
            for pair in chain.windows(2) {
                support
                    .entry((pair[0].clone(), pair[1].clone()))
                    .or_default()
                    .extend(reads);
                outgoing.entry(pair[0].clone()).or_default().extend(reads);
                incoming.entry(pair[1].clone()).or_default().extend(reads);
            }
        }

        // Step 2. Keep the transitions to test. If both junctions have gene assignments, do not
        // pool unrelated genes as controls.
        let transitions: HashMap<(SpliceJunction, SpliceJunction), (ReadSupport, ReadDepth)> = support
            .into_iter()
            .filter(|(transition, _)| {
                if annotated_transitions.contains(transition) {
                    return false;
                }
                match (
                    junctions_to_genes.get(&transition.0.intron_span_key()),
                    junctions_to_genes.get(&transition.1.intron_span_key())
                ) {
                    (Some(genes_1), Some(genes_2)) => !genes_1.is_disjoint(genes_2),
                    _ => true
                }
            })
            .map(|(transition, reads)| {
                let depth: ReadDepth = outgoing[&transition.0].union(&incoming[&transition.1]).count() as ReadDepth;
                (transition, (reads.len() as ReadSupport, depth))
            })
            .collect();

        Self {
            transitions: transitions
        }
    }

    pub fn get_depths(&self) -> HashSet<ReadDepth> {
        self.transitions.values().map(|&(_, depth)| depth).collect()
    }
}

impl SplicingEventReadCounts for SpliceTransitionReadCounts {
    type Input = Vec<SpliceJunction>;

    fn get_read_counts(&self, chain: &Vec<SpliceJunction>) -> Vec<(ReadSupport, ReadDepth)> {
        chain
            .windows(2)
            .filter_map(|pair| self.transitions.get(&(pair[0].clone(), pair[1].clone())).copied())
            .collect()
    }
}


#[cfg(test)]
#[path = "../tests/clustering/splice_transition.rs"]
mod tests;
