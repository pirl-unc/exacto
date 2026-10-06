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
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicUsize, Ordering};
use noodles_bam as bam;
use noodles_bgzf::VirtualPosition;
use rayon::prelude::*;
use rayon::ThreadPool;
use serde::{Deserialize, Serialize};
use tempfile::{NamedTempFile, TempPath};

use crate::prelude::*;


pub fn characterize_rna_reads(
    bam_file: &str,
    read_ids: &Vec<ReadID>,
    read_names_map: &BiMap<ReadName, ReadID>,
    record_positions_map: &HashMap<ReadID, Vec<VirtualPosition>>,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    fasta_map: &FastaMap,
    gene_annotator: &(impl GeneAnnotator + Sync),
    num_threads: usize,
    min_mapping_quality: u16,
    min_terminal_soft_clip_ins_len: u32,

    // Rescue breakpoint
    bkpt_rescue: bool,
    bkpt_rescue_min_ins_len: u32,
    bkpt_rescue_gap_open: i32,
    bkpt_rescue_gap_extend: i32,
    bkpt_rescue_k: u32,
    bkpt_rescue_band_width: u32,
    bkpt_rescue_min_score_fraction: f64,
    bkpt_rescue_min_query_coverage: f64,
    bkpt_rescue_min_placed_fraction: f64,
    bkpt_rescue_max_pieces: u32,

    // Nonsense mediated decay
    nmd_predictor: Option<&NonsenseMediatedDecayPredictor>,

    chunk_size: usize,
    temp_dir: &str
) -> (TempPath, Vec<RNAReadCharacterizationSummary>, Vec<ReadID>) {
    let num_reads: usize = record_positions_map.len();
    let failed_read_ids: Mutex<Vec<ReadID>> = Mutex::new(Vec::new());

    // Characterize each read's junctions and variants.
    let num_characterized: AtomicUsize = AtomicUsize::new(0);
    let mut writer: BufWriter<NamedTempFile> = BufWriter::new(NamedTempFile::new_in(temp_dir).unwrap());
    let mut offset: u64 = 0;
    let mut summaries: Vec<RNAReadCharacterizationSummary> = Vec::with_capacity(num_reads);
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    for chunk in read_ids.chunks(chunk_size) {
        let results: Vec<TranscriptModel> = thread_pool.install(|| {
            chunk
                .par_iter()
                .map_init(
                    // One reader per worker thread - reused across every read ID it handles.
                    || {
                        let mut reader = bam::io::reader::Builder::default()
                            .build_from_path(bam_file)
                            .unwrap();
                        reader.read_header().unwrap(); // consume header before seeking
                        reader
                    },
                    |reader, &read_id| -> Option<TranscriptModel> {
                        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            // Get the BAM records.
                            let records: Vec<bam::Record> = fetch_bam_records_for_read_id(
                                reader,
                                read_id,
                                &record_positions_map
                            );

                            characterize_rna_read(
                                read_id,
                                records,
                                gene_annotator,
                                &chromosome_names_map,
                                &fasta_map,
                                min_mapping_quality,
                                min_terminal_soft_clip_ins_len,
                                bkpt_rescue,
                                bkpt_rescue_min_ins_len,
                                bkpt_rescue_gap_open,
                                bkpt_rescue_gap_extend,
                                bkpt_rescue_k,
                                bkpt_rescue_band_width,
                                bkpt_rescue_min_score_fraction,
                                bkpt_rescue_min_query_coverage,
                                bkpt_rescue_min_placed_fraction,
                                bkpt_rescue_max_pieces,
                                nmd_predictor
                            )
                        }));

                        let result: Option<TranscriptModel> = match result {
                            Ok(result) => result,
                            Err(_) => {
                                let read_name: ReadName = read_names_map.get_by_right(&read_id).unwrap().clone();
                                eprintln!("Failed to characterize junctions and variants for the read name {}", read_name);
                                failed_read_ids.lock().unwrap().push(read_id);
                                None
                            }
                        };

                        // Log the progress.
                        let num_done: usize = num_characterized.fetch_add(1, Ordering::Relaxed) + 1;
                        if (num_done * 10) / num_reads != ((num_done - 1) * 10) / num_reads {
                            log_info!(
                                "Characterized {}% of reads ({}/{})",
                                (num_done * 100) / num_reads,
                                num_done,
                                num_reads
                            );
                        }

                        result
                    }
                )
                .flatten()
                .collect()
        });

        // Append each model and record its address and junctions.
        for transcript_model in results.iter() {
            let bytes: Vec<u8> = bincode::serialize(transcript_model)
                .expect("Failed to serialize transcript model.");
            writer.write_all(&bytes).unwrap();
            let (start_chromosome, reference_start): (ReferenceChromosomeID, ReferencePosition) = transcript_model.get_reference_start_position();
            let (end_chromosome, reference_end): (ReferenceChromosomeID, ReferencePosition) = transcript_model.get_reference_end_position();

            // The read ends at a position when the exon there holds the first or the last base
            // of the read. On the reverse strand the read runs toward lower positions.
            let first_exon: &TranscriptModelExon = transcript_model.get_exons().first().unwrap();
            let last_exon: &TranscriptModelExon = transcript_model.get_exons().last().unwrap();
            let starts_at_first_exon: bool = first_exon.read_start_position == 0;
            let ends_at_last_exon: bool = last_exon.read_end_position + 1 == transcript_model.get_alignment_model().num_bases();
            summaries.push(RNAReadCharacterizationSummary {
                chromosome: (start_chromosome == end_chromosome).then_some(start_chromosome),
                reference_start,
                reference_end,
                ends_at_reference_start: match (&first_exon.reference_strand, &last_exon.reference_strand) {
                    (Strand::Forward, _) => starts_at_first_exon,
                    (_, Strand::Forward) => false,
                    _ => ends_at_last_exon
                },
                ends_at_reference_end: match (&first_exon.reference_strand, &last_exon.reference_strand) {
                    (_, Strand::Forward) => ends_at_last_exon,
                    (Strand::Forward, _) => false,
                    _ => starts_at_first_exon
                },
                read_id: transcript_model.get_read_id(),
                splice_junctions: transcript_model.get_splice_junctions_key(),
                exons: transcript_model
                    .get_exons()
                    .iter()
                    .map(|exon| (exon.reference_chromosome_id, exon.reference_start, exon.reference_end))
                    .collect(),
                offset,
                length: bytes.len() as u64
            });
            offset += bytes.len() as u64;
        }
    }

    let transcript_models_file: TempPath = writer.into_inner().unwrap().into_temp_path();
    let mut failed_read_ids: Vec<ReadID> = failed_read_ids.into_inner().unwrap();
    failed_read_ids.sort_unstable();
    if !failed_read_ids.is_empty() {
        log_info!("{} of {} reads could not be characterized.", failed_read_ids.len(), read_ids.len());
    }
    (transcript_models_file, summaries, failed_read_ids)
}


pub fn load_transcript_models<'a>(
    file: &TempPath,
    summaries: impl IntoIterator<Item = &'a RNAReadCharacterizationSummary>
) -> Vec<TranscriptModel> {
    let mut addresses: Vec<(u64, u64)> = summaries
        .into_iter()
        .map(|summary| (summary.offset, summary.length))
        .collect();
    addresses.sort_unstable();
    let mut reader: BufReader<File> = BufReader::new(File::open(file).unwrap());
    let mut models: Vec<TranscriptModel> = Vec::with_capacity(addresses.len());
    for (offset, length) in addresses.into_iter() {
        reader.seek(SeekFrom::Start(offset)).unwrap();
        models.push(
            bincode::deserialize_from(reader.by_ref().take(length))
                .expect("Failed to deserialize transcript model.")
        );
    }
    models
}


/// Slide an insertion/junction boundary as far as it will go toward the far end of the intron,
/// and return that boundary's new position.
///
/// Normalize junction
pub fn normalise_junction_boundary(
    position_1: ReferencePosition,
    position_2: ReferencePosition,
    insertion_sequence: &str,
    boundary: JunctionBoundary,
    chromosome: &str,
    fasta_map: &FastaMap
) -> ReferencePosition {
    let start: ReferencePosition = match boundary {
        JunctionBoundary::Donor => position_1,
        JunctionBoundary::Acceptor => position_2
    };
    if insertion_sequence.is_empty() || position_1 == position_2 {
        return start;
    }

    // A junction's `position_1` is the first intron base in transcript order and `position_2` the
    // last, so their order is the strand: `position_1 < position_2` forward, the reverse otherwise.
    let forward: bool = position_1 < position_2;
    let transcript_step: i64 = if forward {
        1
    } else {
        -1
    };

    // Both boundaries slide toward the far end of the intron, which for the acceptor means against
    // transcript order. `limit` is that far end: reaching it would dissolve the junction.
    let (limit, step): (ReferencePosition, i64) = match boundary {
        JunctionBoundary::Donor => (position_2, transcript_step),
        JunctionBoundary::Acceptor => (position_1, -transcript_step)
    };

    // The insertion bases that can cross the boundary are the ones nearest it: the leading bases for
    // a donor, the trailing bases for an acceptor.
    let insertion_bases: Vec<char> = match boundary {
        JunctionBoundary::Donor => insertion_sequence.chars().collect(),
        JunctionBoundary::Acceptor => insertion_sequence.chars().rev().collect()
    };

    let mut position: i64 = start as i64;
    for base in insertion_bases.into_iter() {
        if position == limit as i64 {
            break;
        }
        let Some(reference_base) = transcript_strand_base(fasta_map, chromosome, position, forward) else {
            break;
        };
        if reference_base != base.to_ascii_uppercase() {
            break;
        }
        position += step;
    }

    position as ReferencePosition
}


fn characterize_rna_read(
    read_id: ReadID,
    records: Vec<bam::Record>,
    gene_annotator: &(impl GeneAnnotator + Sync),
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    fasta_map: &FastaMap,
    min_mapping_quality: u16,
    min_terminal_soft_clip_ins_len: u32,

    // Breakpoint rescue
    bkpt_rescue: bool,
    bkpt_rescue_min_ins_len: u32,
    bkpt_rescue_gap_open: i32,
    bkpt_rescue_gap_extend: i32,
    bkpt_rescue_k: u32,
    bkpt_rescue_band_width: u32,
    bkpt_rescue_min_score_fraction: f64,
    bkpt_rescue_min_query_coverage: f64,
    bkpt_rescue_min_placed_fraction: f64,
    bkpt_rescue_max_pieces: u32,

    // Nonsense mediated decay
    nmd_predictor: Option<&NonsenseMediatedDecayPredictor>
) -> Option<TranscriptModel> {
    // Step 1. Get the FASTX read sequence.
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);

    // Step 2. Get the base quality scores. AlignmentModel gives Q60 to a read stored without them.
    let base_quality_scores: Vec<BaseQuality> = get_bam_fastx_base_quality_scores(&records);

    // Step 3. Get references to the BAM records.
    let bam_records: Vec<Arc<bam::Record>> = records
        .into_iter()
        .map(Arc::new)
        .collect();

    // Step 4. Construct an instance of AlignmentModel.
    let mut alignment_model: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &base_quality_scores,
        &bam_records
    );

    // Step 5. Filter based on mapping quality.
    let max_mapping_quality: MappingQuality = alignment_model.get_records().iter()
        .map(|r| get_alignment_mapping_quality(&r.record))
        .max()
        .unwrap_or(0);
    if min_mapping_quality > max_mapping_quality {
        return None;
    }

    // Step 6. Identify reference transcript matches.
    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &identify_transcript_model_exons(&alignment_model),
        &identify_transcript_model_splice_junctions(&alignment_model),
        gene_annotator,
        chromosome_names_map
    );

    // Step 7. Construct an instance of TranscriptModel.
    let mut transcript_model: TranscriptModel = TranscriptModel::new(
        alignment_model,
        reference_transcript_matches,
        gene_annotator,
        chromosome_names_map,
        fasta_map
    );

    // Step 8. Identify RNA variants.
    let breakpoint_rescue: Option<RNABreakpointRescue> = bkpt_rescue.then(|| RNABreakpointRescue {
        gene_annotator,
        chromosome_names_map,
        fasta_map,
        min_insertion_length: bkpt_rescue_min_ins_len,
        gap_open: bkpt_rescue_gap_open,
        gap_extend: bkpt_rescue_gap_extend,
        k: bkpt_rescue_k,
        band_width: bkpt_rescue_band_width,
        min_score_fraction: bkpt_rescue_min_score_fraction,
        min_query_coverage: bkpt_rescue_min_query_coverage,
        min_placed_fraction: bkpt_rescue_min_placed_fraction,
        max_pieces: bkpt_rescue_max_pieces
    });
    let variant_records: Vec<VariantRecord> = RNAVariantRecordCaller::new(
        min_mapping_quality,
        0,
        min_terminal_soft_clip_ins_len,
        breakpoint_rescue
    ).call(&transcript_model);
    transcript_model.set_variant_records(variant_records);

    // Step 9. Predict nonsense mediated decay.
    if let Some(nmd_predictor) = nmd_predictor {
        transcript_model.set_nmd_predictions(nmd_predictor.predict(&transcript_model));
    }

    Some(transcript_model)
}


fn transcript_strand_base(
    fasta_map: &FastaMap,
    chromosome: &str,
    position: i64,
    forward: bool
) -> Option<char> {
    if position < 1 || position > fasta_map.get_length(chromosome) as i64 {
        return None;
    }
    let base: char = fasta_map
        .get_sequence(chromosome, position as usize, position as usize)
        .chars()
        .next()?
        .to_ascii_uppercase();
    Some(if forward {
        base
    } else {
        match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            other => other
        }
    })
}


#[cfg(test)]
#[path = "../tests/transcript/rna_read_characterization.rs"]
mod tests;