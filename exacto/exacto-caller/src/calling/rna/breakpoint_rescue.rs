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
use std::ops::Range;

use crate::prelude::*;
use crate::reference::reference_transcript_alignment::{place_on_reference_transcript, ReferenceTranscriptPlacement};


pub(crate) fn retype_insertion_on_reference_transcripts(
    variant_record: &VariantRecord,
    read_length: u32,
    reference_transcript_sequences: &Vec<ReferenceTranscriptSequence>,
    min_insertion_length: u32,
    gap_open: i32,
    gap_extend: i32,
    k: u32,
    band_width: u32,
    min_score_fraction: f64,
    min_query_coverage: f64,
    min_placed_fraction: f64,
    max_pieces: u32
) -> Option<Vec<VariantRecord>> {
    assert_eq!(
        *variant_record.get_variant_type(),
        VariantType::Insertion
    );

    // The inserted bases, in read order.
    let query: &str = variant_record.get_sequence();
    let query_length: u32 = query.len() as u32;

    if query_length < min_insertion_length {
        return None;
    }

    // A terminal soft clip is an insertion on the first or the last bases of the read. The read
    // reaches the reference on the clip's anchor side only, so the walk below takes no hop from or
    // to the other side.
    let is_leading_clip: bool = variant_record.get_read_position_1() == 0;
    let is_trailing_clip: bool = variant_record.get_read_position_2() + 1 == read_length;

    for rts in reference_transcript_sequences {
        // Both insertion breakpoints must lie on this transcript.
        // These positions are retained as the original left and right flanks when
        // constructing breakpoint adjacencies below.
        let (Some(a), Some(b)) = (
            rts.get_reference_base_position(variant_record.get_position_1()),
            rts.get_reference_base_position(variant_record.get_position_2()),
        ) else {
            continue;
        };

        let (t_left, t_right): (ReferenceTranscriptPosition, ReferenceTranscriptPosition) = (a.min(b) as ReferenceTranscriptPosition, a.max(b) as ReferenceTranscriptPosition);

        // Align against the entire spliced reference-transcript sequence.
        let window: Range<ReferenceTranscriptPosition> = 0..rts.get_sequence().len() as u32;

        // Greedy decomposition: find the best placement, mask it in the query, and repeat.
        let mut masked: Vec<u8> = query.as_bytes().to_vec();
        let mut pieces: Vec<ReferenceTranscriptPlacement> = Vec::new();

        while max_pieces > pieces.len() as u32  {
            let remaining: &str = std::str::from_utf8(&masked).unwrap();
            let Some(piece) = place_on_reference_transcript(
                remaining,
                rts,
                window.clone(),
                gap_open,
                gap_extend,
                k,
                band_width,
                min_score_fraction,
                min_query_coverage
            ) else {
                break;
            };

            // A reverse hit reports coordinates on the reverse-complemented
            // query. Convert them back to coordinates on the original query.
            let query_range: Range<u32> = if piece.is_forward {
                piece.query_range.clone()
            } else {
                (query_length - piece.query_range.end)..(query_length - piece.query_range.start)
            };

            // Prevent the next alignment from rediscovering this piece.
            masked[query_range.start as usize..query_range.end as usize].fill(b'N');

            pieces.push(ReferenceTranscriptPlacement {
                query_range: query_range,
                reference_range: piece.reference_range,
                is_forward: piece.is_forward,
                score: piece.score
            });
        }

        // Require the placed pieces to account for enough of the insertion.
        let placed: u32 = pieces
            .iter()
            .map(|piece| piece.query_range.len() as u32)
            .sum();

        if (placed as f64) < min_placed_fraction * query_length as f64 {
            continue;
        }

        // Put pieces in insertion/read order rather than transcript order.
        pieces.sort_by_key(|piece| piece.query_range.start);

        // Walk, in read order:
        // left insertion flank -> placed pieces -> right insertion flank
        let strand: Strand = rts.get_strand().clone();

        let flipped: Strand = match strand {
            Strand::Forward => Strand::Reverse,
            _ => Strand::Forward
        };

        let (first_flank, last_flank, read_strand): (ReferenceTranscriptPosition, ReferenceTranscriptPosition, Strand) =
            if *variant_record.get_strand_1() == strand {
                (t_left, t_right, strand.clone())
            } else {
                (t_right, t_left, flipped.clone())
            };

        let enter = |piece: &ReferenceTranscriptPlacement| -> (ReferenceTranscriptPosition, Strand) {
            if piece.is_forward {
                (
                    piece.reference_range.start,
                    strand.clone()
                )
            } else {
                (
                    piece.reference_range.end - 1,
                    flipped.clone()
                )
            }
        };

        let leave = |piece: &ReferenceTranscriptPlacement| -> (ReferenceTranscriptPosition, Strand) {
            if piece.is_forward {
                (
                    piece.reference_range.end - 1,
                    strand.clone(),
                )
            } else {
                (
                    piece.reference_range.start,
                    flipped.clone(),
                )
            }
        };

        let mut adjacencies: Vec<(
            (ReferenceTranscriptPosition, Strand),
            (ReferenceTranscriptPosition, Strand),
        )> = Vec::new();

        // A leading clip has no flank before it in the read, a trailing clip none after it.
        let mut previous: Option<(ReferenceTranscriptPosition, Strand)> = (!is_leading_clip).then(|| (
            first_flank,
            read_strand.clone(),
        ));

        for piece in &pieces {
            if let Some(previous) = previous {
                adjacencies.push((
                    previous,
                    enter(piece)
                ));
            }

            previous = Some(leave(piece));
        }

        if let Some(previous) = previous.filter(|_| !is_trailing_clip) {
            adjacencies.push((
                previous,
                (last_flank, read_strand.clone())
            ));
        }

        // Ordinary continuation through the transcript is not a junction: forward on
        // the transcript's strand, backward on the flipped one. Convert every other
        // adjacency into a breakpoint.
        let mut records: Vec<VariantRecord> = Vec::new();

        for ((leave_offset, leave_strand), (enter_offset, enter_strand)) in adjacencies {
            let is_continuation: bool = leave_strand == enter_strand
                && if leave_strand == strand {
                    enter_offset == leave_offset + 1
                } else {
                    enter_offset + 1 == leave_offset
                };

            if is_continuation {
                continue;
            }

            let leaving_base: &ReferenceBase = rts.get_base(leave_offset as usize);
            let chromosome: ReferenceChromosomeID = leaving_base.reference_chromosome_id;
            let leave_position: ReferencePosition = leaving_base.reference_position;

            let entering_base: &ReferenceBase = rts.get_base(enter_offset as usize);
            let enter_position: ReferencePosition = entering_base.reference_position;

            records.push(VariantRecord::new(
                variant_record.get_read_id(),
                variant_record.get_read_position_1(),
                variant_record.get_read_position_2(),
                GraphOperation::new(
                    chromosome,
                    leave_position,
                    leave_strand.clone(),
                    GraphOperationType::for_breakpoint(
                        &leave_strand,
                        false,
                    ),
                    chromosome,
                    enter_position,
                    enter_strand.clone(),
                    GraphOperationType::for_breakpoint(
                        &enter_strand,
                        true,
                    ),
                    "".into(),
                    VariantType::Breakpoint
                ),
            ));
        }

        if !records.is_empty() {
            return Some(records);
        }
    }

    None
}


pub struct RNABreakpointRescue<'a> {
    pub gene_annotator: &'a (dyn GeneAnnotator + Sync),
    pub chromosome_names_map: &'a BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    pub fasta_map: &'a FastaMap,
    pub min_insertion_length: u32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub k: u32,
    pub band_width: u32,
    pub min_score_fraction: f64,
    pub min_query_coverage: f64,
    pub min_placed_fraction: f64,
    pub max_pieces: u32
}

impl RNABreakpointRescue<'_> {
    pub fn retype_insertions(
        &self,
        transcript_model: &TranscriptModel,
        variant_records: Vec<VariantRecord>
    ) -> Vec<VariantRecord> {
        if transcript_model.get_reference_transcript_matches().is_empty() {
            return variant_records;
        }

        // No insertion to retype: the reference transcript sequences are not needed.
        let has_large_insertion: bool = variant_records.iter().any(|record| {
            record.get_variant_type() == &VariantType::Insertion
                && self.min_insertion_length <= record.get_sequence().len() as u32
        });
        if !has_large_insertion {
            return variant_records;
        }

        // Get reference transcript sequences and indices.
        let mut reference_transcript_sequences: Vec<ReferenceTranscriptSequence> = Vec::new();
        for reference_transcript_match in transcript_model.get_reference_transcript_matches().iter() {
            let reference_transcript: &Transcript = self.gene_annotator.get_transcript(
                reference_transcript_match.get_reference_transcript_id()
            ).unwrap();
            let sequence: ReferenceTranscriptSequence = ReferenceTranscriptSequence::from_reference_transcript(
                reference_transcript,
                self.chromosome_names_map,
                self.fasta_map
            );
            reference_transcript_sequences.push(sequence);
        }
        
        // Retype insertions.
        let mut records: Vec<VariantRecord> = Vec::with_capacity(variant_records.len());
        for record in variant_records {
            let is_large_insertion: bool = record.get_variant_type() == &VariantType::Insertion 
                && self.min_insertion_length <= record.get_sequence().len() as u32;
            
            if !is_large_insertion {
                records.push(record);
                continue;
            }
            
            // Retype insertion.
            let results: Option<Vec<VariantRecord>> = retype_insertion_on_reference_transcripts(
                &record,
                transcript_model.get_alignment_model().num_bases(),
                &reference_transcript_sequences,
                self.min_insertion_length,
                self.gap_open,
                self.gap_extend,
                self.k,
                self.band_width,
                self.min_score_fraction,
                self.min_query_coverage,
                self.min_placed_fraction,
                self.max_pieces
            );
            
            match results {
                Some(retyped) => records.extend(retyped),
                None => records.push(record)
            }
        }
        records
    }
}


#[cfg(test)]
#[path = "../../tests/calling/rna/breakpoint_rescue.rs"]
mod tests;