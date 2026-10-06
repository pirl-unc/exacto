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
use exacto_core::prelude::*;
use noodles_bgzf::VirtualPosition;
use std::collections::HashMap;
use std::ops::RangeInclusive;

use crate::modeling::rna_read_modeling::{map_rna_reads, model_rna_read};
use crate::modeling::transcript_model_topology::classify_topology;
use crate::prelude::*;
use crate::stitching::end_matching::*;
use crate::stitching::end_resolution::*;


pub(crate) fn stitch_rna_reads(
    bam_file: &str,
    read_names_map: &BiMap<Box<str>, usize>,
    record_positions_map: &HashMap<usize, Vec<VirtualPosition>>,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap,
    gene_annotator: &(impl GeneAnnotator + Sync),

    // Read filtering
    min_mapping_quality: u16,

    // Reference transcript match selection
    min_num_splice_junction_matches: usize,

    // Polyadenylation hexamer signal
    pas_search_size: usize,
    pas_hexamers: &[&str],
    pas_start_offset_range: &RangeInclusive<usize>,

    // Polyadenylation tail
    polya_tail_min_adenosine_fraction: f64,
    polya_window_size: usize,
    polya_min_consecutive_a: usize,
    num_threads: usize
) -> Vec<StitchedTranscript> {
    map_rna_reads(bam_file, read_names_map, record_positions_map, num_threads, |read_id, records| {
        // The read as sequenced, for a read that passes through.
        let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);

        // Modeling runs the caller on the read's alignment records, and a read it cannot
        // handle passes through. Stitching runs outside the guard: a panic there is a bug or
        // a bad option, and fails the call.
        let transcript_model = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| model_rna_read(
            read_id,
            records,
            gene_annotator,
            chromosome_names_map,
            fasta_map,
            min_mapping_quality
        )));
        match transcript_model {
            Ok(Some(transcript_model)) => stitch_transcript_model(
                &transcript_model,
                chromosome_names_map,
                gene_annotator,
                fasta_map,
                min_num_splice_junction_matches,
                pas_search_size,
                pas_hexamers,
                pas_start_offset_range,
                polya_tail_min_adenosine_fraction,
                polya_window_size,
                polya_min_consecutive_a
            ),
            Ok(None) => build_unmodeled_stitched_transcript(read_id, read_sequence, "mapping_quality"),
            Err(_) => {
                eprintln!(
                    "Failed to model the read name {}; it passes through unstitched.",
                    read_names_map.get_by_right(&read_id).unwrap()
                );
                build_unmodeled_stitched_transcript(read_id, read_sequence, "modeling_error")
            }
        }
    })
}


/// A read without a transcript model, as sequenced.
fn build_unmodeled_stitched_transcript(
    read_id: usize,
    read_sequence: Box<str>,
    unmodeled_reason: &'static str
) -> StitchedTranscript {
    StitchedTranscript {
        read_id: read_id,
        original_length: read_sequence.len() as u32,
        stitched_length: read_sequence.len() as u32,
        stitched_sequence: read_sequence,
        five_prime_stitch_end_index: None,
        three_prime_stitch_start_index: None,
        five_prime_end_resolution: None,
        three_prime_end_resolution: None,
        unmodeled_reason: Some(unmodeled_reason)
    }
}


fn build_stitched_transcript(
    transcript_model: &TranscriptModel,
    five_prime_end_resolution: &Option<EndResolution>,
    three_prime_end_resolution: &Option<EndResolution>
) -> StitchedTranscript {
    let bases: &[AlignmentModelBase] = transcript_model
        .get_alignment_model()
        .get_bases();

    // Step 1. Decide which read bases are retained.
    // A resolved end keeps the read from its join position;
    // an unresolved end keeps the read as sequenced.
    let first_read_position: usize = match five_prime_end_resolution {
        Some(resolution) => resolution.read_join_position as usize,
        None => 0
    };
    let last_read_position: usize = match three_prime_end_resolution {
        Some(resolution) => resolution.read_join_position as usize,
        None => bases.len() - 1
    };

    // A 3' roll-back moves its join inward, so on a short read it could pass the 5' join.
    assert!(
        first_read_position <= last_read_position,
        "read {}: the 5' join position ({}) passes the 3' join position ({})",
        transcript_model.get_read_id(), first_read_position, last_read_position
    );

    // Step 2. Assemble the sequence: 5' stitch, retained read bases, 3' stitch.
    // Each boundary is read off the sequence's length at the moment it is crossed.
    // Nucleotides are ASCII, so a byte offset is a base offset.
    let mut stitched_sequence: String = String::new();

    let mut five_prime_stitch_end_index: Option<usize> = None;
    if let Some(resolution) = five_prime_end_resolution {
        stitched_sequence.push_str(&resolution.stitch_sequence);
        five_prime_stitch_end_index = Some(stitched_sequence.len());
    }

    for base in bases[first_read_position..=last_read_position].iter() {
        stitched_sequence.push_str(base.get_nucleotide().as_str());
    }

    let mut three_prime_stitch_start_index: Option<usize> = None;
    if let Some(resolution) = three_prime_end_resolution {
        three_prime_stitch_start_index = Some(stitched_sequence.len());
        stitched_sequence.push_str(&resolution.stitch_sequence);
    }

    StitchedTranscript {
        read_id: transcript_model.get_read_id(),
        original_length: bases.len() as u32,
        stitched_length: stitched_sequence.len() as u32,
        stitched_sequence: stitched_sequence.into_boxed_str(),
        five_prime_stitch_end_index: five_prime_stitch_end_index,
        three_prime_stitch_start_index: three_prime_stitch_start_index,
        five_prime_end_resolution: five_prime_end_resolution.clone(),
        three_prime_end_resolution: three_prime_end_resolution.clone(),
        unmodeled_reason: None
    }
}


fn stitch_transcript_model(
    transcript_model: &TranscriptModel,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    gene_annotator: &(impl GeneAnnotator + Sync),
    fasta_map: &FastaMap,

    // Reference transcript match selection
    min_num_splice_junction_matches: usize,

    // Polyadenylation hexamer signal
    pas_search_size: usize,
    pas_hexamers: &[&str],
    pas_start_offset_range: &RangeInclusive<usize>,

    // Polyadenylation tail
    polya_tail_min_adenosine_fraction: f64,
    polya_window_size: usize,
    polya_min_consecutive_a: usize
) -> StitchedTranscript {
    // Step 1. Classify the transcript model topology.
    let topology: TranscriptModelTopology = classify_topology(transcript_model);
    if topology == TranscriptModelTopology::BackSplicing {
        return build_stitched_transcript(
            transcript_model,
            &None,
            &None
        );
    }

    // Step 2. Identify the 5' end reference transcript matches.
    let five_prime_end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        chromosome_names_map,
        gene_annotator,
        &TranscriptTerminus::FivePrime,
        min_num_splice_junction_matches
    );

    // Step 3. Identify the 3' end reference transcript matches.
    let three_prime_end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        chromosome_names_map,
        gene_annotator,
        &TranscriptTerminus::ThreePrime,
        min_num_splice_junction_matches
    );

    // Step 4. Resolve the 5' end.
    let five_prime_end_resolution: Option<EndResolution> = resolve_five_prime_end(
        transcript_model,
        &five_prime_end_matches,
        gene_annotator,
        chromosome_names_map,
        fasta_map
    );

    // Step 5. Resolve the 3' end.
    let three_prime_end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &three_prime_end_matches,
        gene_annotator,
        chromosome_names_map,
        fasta_map,
        pas_search_size,
        pas_hexamers,
        pas_start_offset_range,
        polya_tail_min_adenosine_fraction,
        polya_window_size,
        polya_min_consecutive_a
    );

    // Step 6. Build a stitched transcript.
    build_stitched_transcript(
        transcript_model,
        &five_prime_end_resolution,
        &three_prime_end_resolution
    )
}


#[cfg(test)]
#[path = "../tests/stitching/transcript_model_stitching.rs"]
mod tests;