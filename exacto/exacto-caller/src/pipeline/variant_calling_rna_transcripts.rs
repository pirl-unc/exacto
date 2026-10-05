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
use exacto_core::prelude::*;
use exacto_core::log_info;
use noodles_bgzf::VirtualPosition;
use std::collections::{HashMap, HashSet};
use tempfile::TempPath;

use crate::prelude::*;


pub fn identify_rna_transcript_variants(
    bam_file: &str,
    reference_genome_fasta_file: &str,
    gene_annotator: &(impl GeneAnnotator + Sync),
    options: &IdentifyRNATranscriptVariantsOptions,
    num_threads: usize,
    temp_dir: &str
) -> TranscriptModelSet {
    // Step 1. Check the inputs.
    check_alignment_inputs(
        bam_file,
        None,
        reference_genome_fasta_file
    );

    // Step 2. Get chromosome IDs and names.
    log_info!("Getting chromosomes.");
    let chromosome_names_map: BiMap<ReferenceChromosomeName, ReferenceChromosomeID> = create_chromosome_names_map(bam_file);

    // Step 3. Index BAM records by read names.
    log_info!("Indexing BAM records by read names.");
    let (record_positions_map, read_names_map): (HashMap<ReadID, Vec<VirtualPosition>>, BiMap<ReadName, ReadID>) = index_bam_records(
        bam_file,
        true,
        num_threads
    );

    // Step 4. Load the reference genome FASTA file.
    log_info!("Loading the reference genome FASTA file.");
    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);

    // Step 5. Characterize RNA reads.
    log_info!("Characterizing RNA reads.");
    let read_ids: Vec<ReadID> = read_names_map.right_values().copied().collect();
    let start_codons: HashSet<&str> = options.annotation.start_codons
        .iter()
        .map(|codon| codon.as_str())
        .collect();
    let nmd_predictor: Option<NonsenseMediatedDecayPredictor> = options.annotation.predict_nonsense_mediated_decay
        .then(|| NonsenseMediatedDecayPredictor {
            translation_strategy: &options.annotation.translation_strategy,
            start_codons: &start_codons,
            distance_threshold: options.annotation.nmd_distance_threshold
        });
    let (temp_file, summaries, failed_read_ids): (TempPath, Vec<RNAReadCharacterizationSummary>, Vec<ReadID>) = characterize_rna_reads(
        bam_file,
        &read_ids,
        &read_names_map,
        &record_positions_map,
        &chromosome_names_map,
        &fasta_map,
        gene_annotator,
        num_threads,
        options.calling.min_mapping_quality,
        options.calling.min_terminal_soft_clip_ins_len,
        options.calling.bkpt_rescue,
        options.calling.bkpt_rescue_min_ins_len,
        options.calling.bkpt_rescue_realignment_gap_open_score,
        options.calling.bkpt_rescue_realignment_gap_extend_score,
        options.calling.bkpt_rescue_realignment_k,
        options.calling.bkpt_rescue_realignment_band_width,
        options.calling.bkpt_rescue_realignment_min_score_fraction,
        options.calling.bkpt_rescue_realignment_min_query_coverage,
        options.calling.bkpt_rescue_realignment_min_placed_fraction,
        options.calling.bkpt_rescue_realignment_max_pieces,
        nmd_predictor.as_ref(),
        options.calling.chunk_size,
        temp_dir
    );

    // A read that fails is left out and named on stderr. When most reads fail, the input is at
    // fault rather than the reads, and a set built from the rest would look complete.
    if failed_read_ids.len() * 2 > read_ids.len() {
        panic!(
            "{} of the {} reads of {bam_file} could not be characterized, {} among them.",
            failed_read_ids.len(),
            read_ids.len(),
            read_names_map.get_by_right(&failed_read_ids[0]).unwrap()
        );
    }

    // Step 6. Load transcript models.
    log_info!("Loading transcript models.");
    let mut transcript_model_set: TranscriptModelSet = TranscriptModelSet::new();
    let mut transcript_models: Vec<TranscriptModel> = load_transcript_models(&temp_file, &summaries);
    transcript_models.sort_unstable_by_key(|transcript_model| transcript_model.get_read_id());
    for transcript_model in transcript_models {
        transcript_model_set.add(transcript_model);
    }
    transcript_model_set.load_read_names(read_names_map);
    transcript_model_set.load_chromosome_names(chromosome_names_map);

    transcript_model_set
}


#[cfg(test)]
#[path = "../tests/pipeline/variant_calling_rna_transcripts.rs"]
mod tests;