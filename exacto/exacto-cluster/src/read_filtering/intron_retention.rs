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


use bimap::BiMap;
use exacto_caller::prelude::*;
use exacto_core::log_info;
use exacto_core::prelude::{FastaMap, Strand};
use std::collections::HashSet;

use crate::reference::junction_orientation::is_aligned_against_transcript;


pub(crate) fn identify_nascent_rna_read_ids_per_strand(
    summaries: &mut [RNAReadCharacterizationSummary],
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap,
    sequencing_error: f64,
    max_fpr: f64,
    num_threads: usize
) -> HashSet<usize> {
    // 0 forward, 1 no strand, 2 reverse
    let strand_rank = |summary: &RNAReadCharacterizationSummary| -> u8 {
        let Some(junction) = summary.splice_junctions.first() else {
            return 1;
        };
        let is_against: bool = is_aligned_against_transcript(&summary.splice_junctions, chromosome_names_map, fasta_map);
        match (&junction.strand_1, is_against) {
            (Strand::Forward, false) | (Strand::Reverse, true) => 0,
            (Strand::Reverse, false) | (Strand::Forward, true) => 2,
            _ => 1
        }
    };
    summaries.sort_by_cached_key(|summary| (strand_rank(summary), summary.read_id));
    let first_junction_less: usize = summaries.partition_point(|summary| strand_rank(summary) < 1);
    let first_reverse: usize = summaries.partition_point(|summary| strand_rank(summary) < 2);

    let mut nascent_read_ids: HashSet<usize> = HashSet::new();
    for strand_summaries in [&summaries[..first_reverse], &summaries[first_junction_less..]] {
        log_info!("Building minimum RNA splice junction read support index.");
        let read_counts: RNAJunctionReadCounts = RNAJunctionReadCounts::new(strand_summaries, num_threads);
        let read_support_index: RNAJunctionReadSupportIndex = RNAJunctionReadSupportIndex::new(
            &read_counts.get_depths(),
            sequencing_error,
            max_fpr,
            num_threads
        );
        log_info!("Identifying reads with unsupported intron retention.");
        let read_filter: RNAJunctionReadSupportFilter<RNAJunctionReadCounts> = RNAJunctionReadSupportFilter::new(
            &read_counts,
            &read_support_index
        );
        nascent_read_ids.extend(
            strand_summaries
                .iter()
                .filter(|summary| !read_filter.passes(summary))
                .map(|summary| summary.read_id)
        );
    }
    log_info!("Removed {} read(s) with unsupported intron retention.", nascent_read_ids.len());
    nascent_read_ids
}


#[cfg(test)]
#[path = "../tests/read_filtering/intron_retention.rs"]
mod tests;
