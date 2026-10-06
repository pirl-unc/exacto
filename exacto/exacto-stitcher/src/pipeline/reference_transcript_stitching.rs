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
use exacto_core::log_info;
use exacto_core::prelude::*;
use noodles_bam as bam;
use noodles_bgzf::VirtualPosition;
use std::collections::HashMap;

use crate::prelude::*;
use crate::pipeline::options::StitchReferenceTranscriptsOptions;
use crate::stitching::transcript_model_stitching::stitch_rna_reads;


pub fn stitch_reference_transcripts(
    bam_file: &str,
    reference_genome_fasta_file: &str,
    gene_annotator: &(impl GeneAnnotator + Sync),
    options: &StitchReferenceTranscriptsOptions,
    num_threads: usize
) -> StitchedTranscriptSet {
    // Step 1. Index the BAM file.
    log_info!("Fetching BAM index.");
    let (record_positions_map, read_names_map):
        (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(
        bam_file,
        true,
        num_threads
    );

    // Step 2. Check for the cs tag. Every read is modeled from it, so without it every read
    // would pass through unstitched; fail instead.
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    reader.read_header().unwrap();
    if let Some(record) = reader
        .records()
        .map(|result| result.unwrap())
        .find(|record| !record.flags().is_unmapped()) {
        assert!(
            has_tag(&record, "cs"),
            "{} has no cs tag on its first mapped record ({}). Exacto needs an aligner that emits one (e.g. minimap2 --cs).",
            bam_file,
            record.name().map(|name| name.to_string()).unwrap_or_default()
        );
    }

    // Step 3. Create a map of chromosome names and IDs.
    log_info!("Creating a map of chromosome names and IDs.");
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);

    // Step 4. Create a map of the reference genome FASTA.
    log_info!("Creating a map of the reference genome FASTA file.");
    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);

    // Step 5. Model and stitch every read, in read ID order so that the output keeps the
    // BAM's order.
    log_info!("Modeling and stitching RNA reads.");
    let pas_hexamers: Vec<&str> = options.pas_hexamers
        .iter()
        .map(|hexamer| &**hexamer)
        .collect();
    let stitched_transcripts: Vec<StitchedTranscript> = stitch_rna_reads(
        bam_file,
        &read_names_map,
        &record_positions_map,
        &chromosome_names_map,
        &fasta_map,
        gene_annotator,
        options.min_mapping_quality,
        options.min_num_splice_junction_matches,
        options.pas_search_size,
        &pas_hexamers,
        &options.pas_start_offset_range,
        options.polya_tail_min_adenosine_fraction,
        options.polya_window_size,
        options.polya_min_consecutive_adenosine,
        num_threads
    );

    let mut stitched_transcript_set: StitchedTranscriptSet = StitchedTranscriptSet::new();
    for stitched_transcript in stitched_transcripts {
        stitched_transcript_set.add(stitched_transcript);
    }

    stitched_transcript_set
}


#[cfg(test)]
#[path = "../tests/pipeline/transcript_stitching.rs"]
mod tests;