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
use std::collections::HashSet;
use std::ops::{Range, RangeInclusive};


pub(crate) fn identify_transcript_end_trims(
    model: &AlignmentModel,
    blocks: &[TranscriptModelExon],
    chromosome_names_map: &BiMap<Box<str>, u16>,
    gene_annotator: &dyn GeneAnnotator,
    protected: &HashSet<&GraphOperation>
) -> Vec<Range<u32>> {
    terminal_blocks(blocks).into_iter().filter_map(|terminal| {
        // A gap in reference coordinates alone is insufficient: require an actual
        // splice event. Same-chromosome split alignments must not license trimming.
        if !has_anchoring_splice(model, &terminal) {
            return None;
        }

        let block: &TranscriptModelExon = terminal.block;
        let chromosome_name: &Box<str> = chromosome_names_map.get_by_right(&block.reference_chromosome_id)?;

        // Local: only a transcript overlapping the anchoring intron can carry it.
        let transcript_ids: Vec<Box<str>> = gene_annotator
            .get_transcript_ids_overlapping_region(chromosome_name, terminal.intron.0, terminal.intron.1);

        let transcripts: Vec<&Transcript> = transcript_ids
            .iter()
            .filter_map(|transcript_id| gene_annotator.get_transcript(transcript_id))
            .collect();

        let boundary: u32 = resolve_cut_boundary(&terminal, &transcripts)?;

        // Preserve the archived site protection, including flanking anchors. Test
        // intervals directly rather than expanding a large deletion into single sites.
        let removed: RangeInclusive<u32> = if terminal.rightward {
            boundary + 1..=block.reference_end
        } else {
            block.reference_start..=boundary - 1
        };

        let contains_protected_variant: bool = protected.iter().any(|operation| {
            overlaps_protected_operation(
                operation,
                block.reference_chromosome_id,
                &removed
            )
        });

        if contains_protected_variant {
            return None;
        }

        let placements: Vec<(u32, u32)> = (block.read_start_position..=block.read_end_position)
            .filter_map(|read_position| {
                let (_, reference_position, _) = model.get_base(read_position).get_placement().get_coordinate()?;
                Some((read_position, reference_position))
            })
            .collect();

        let (read_position, length): (u32, usize) = cut_from_placements(
            &placements,
            boundary,
            terminal.rightward,
            terminal.end,
            model.num_bases()
        )?;

        Some(read_position..read_position + length as u32)
    }).collect()
}


#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ReadEnd {
    Head,
    Tail
}


#[derive(Clone, Debug)]
struct TerminalBlock<'a> {
    end: ReadEnd,
    block: &'a TranscriptModelExon,
    /// The intron between this block and its neighbour, 1-based inclusive `(first, last)`.
    intron: (u32, u32),
    /// Whether the read leaves the block toward higher reference coordinates: a tail on the
    /// forward strand or a head on the reverse strand.
    rightward: bool
}


#[derive(Clone, Debug, PartialEq, Eq)]
enum EndVerdict {
    /// The exon beyond the intron is the transcript's own terminus: the read may run past it.
    TranscriptEnd,
    /// The exon's outer boundary. The read overhangs it when its block extends beyond.
    Boundary(u32)
}


fn cut_from_placements(
    placements: &[(u32, u32)],
    boundary: u32,
    rightward: bool,
    end: ReadEnd,
    read_length: u32
) -> Option<(u32, usize)> {
    let beyond = |&&(_, reference_position): &&(u32, u32)| -> bool {
        if rightward {
            reference_position > boundary
        } else {
            reference_position < boundary
        }
    };

    match end {
        ReadEnd::Tail => {
            let first: u32 = placements.iter().find(beyond)?.0;
            Some((first, (read_length - first) as usize))
        },
        ReadEnd::Head => {
            let last: u32 = placements.iter().rev().find(beyond)?.0;
            Some((0, last as usize + 1))
        }
    }
}


fn has_anchoring_splice(model: &AlignmentModel, terminal: &TerminalBlock) -> bool {
    let block: &TranscriptModelExon = terminal.block;

    model.get_events().values().any(|event| {
        if *event.get_kind() != AlignmentModelEventKind::Splicing {
            return false;
        }

        let previous: u32 = event.get_prev_read_position();
        let next: u32 = event.get_next_read_position();
        let block_read_range: RangeInclusive<u32> = block.read_start_position..=block.read_end_position;

        if !block_read_range.contains(&previous) && !block_read_range.contains(&next) {
            return false;
        }

        let (Some((chr_a, pos_a, strand_a)), Some((chr_b, pos_b, strand_b))) = (
            model.get_base(previous).get_placement().get_coordinate(),
            model.get_base(next).get_placement().get_coordinate()
        ) else {
            return false;
        };

        chr_a == block.reference_chromosome_id
            && chr_b == chr_a
            && strand_a == &block.reference_strand
            && strand_b == strand_a
            && pos_a.min(pos_b).checked_add(1) == Some(terminal.intron.0)
            && pos_a.max(pos_b).checked_sub(1) == Some(terminal.intron.1)
    })
}


/// Judge an end only against transcripts sharing its strand and anchoring intron.
fn judge_end(terminal: &TerminalBlock, transcript: &Transcript) -> Option<EndVerdict> {
    if transcript.strand != terminal.block.reference_strand {
        return None;
    }

    let intron: Intron = transcript
        .get_introns()
        .into_iter()
        .find(|intron| (intron.start, intron.end) == terminal.intron)?;

    // Intron `k` lies between exon `k` and exon `k + 1` in transcript order, and GENCODE numbers
    // exons 5' to 3' on both strands. The tail sits on the intron's 3' side, the head on its 5'.
    let exon_number: u16 = match terminal.end {
        ReadEnd::Tail => intron.intron_number + 1,
        ReadEnd::Head => intron.intron_number
    };

    let last_exon_number: u16 = transcript.get_sorted_exons().last()?.exon_number;

    let is_transcript_end: bool = match terminal.end {
        ReadEnd::Tail => exon_number == last_exon_number,
        ReadEnd::Head => exon_number == 1
    };

    if is_transcript_end {
        return Some(EndVerdict::TranscriptEnd);
    }

    let exon: &Exon = transcript.get_exon_by_number(exon_number)?;

    Some(EndVerdict::Boundary(if terminal.rightward { exon.end } else { exon.start }))
}


fn overlaps_protected_operation(
    operation: &GraphOperation,
    chromosome: u16,
    removed: &RangeInclusive<u32>
) -> bool {
    let chromosome_1: u16 = operation.get_chromosome_1();
    let chromosome_2: u16 = operation.get_chromosome_2();
    let position_1: u32 = operation.get_position_1();
    let position_2: u32 = operation.get_position_2();

    if (chromosome_1 == chromosome && removed.contains(&position_1))
        || (chromosome_2 == chromosome && removed.contains(&position_2)) {
        return true;
    }

    chromosome_1 == chromosome && chromosome_2 == chromosome
        && matches!(operation.get_variant_type(),
            VariantType::SingleNucleotideVariant
            | VariantType::MultiNucleotideVariant
            | VariantType::Insertion
            | VariantType::Deletion)
        && position_1.min(position_2) <= *removed.end()
        && position_1.max(position_2) >= *removed.start()
}


fn resolve_cut_boundary(terminal: &TerminalBlock, transcripts: &[&Transcript]) -> Option<u32> {
    let mut boundaries: Vec<u32> = Vec::new();

    for verdict in transcripts.iter().filter_map(|transcript| judge_end(terminal, transcript)) {
        match verdict {
            EndVerdict::TranscriptEnd => return None,
            EndVerdict::Boundary(boundary) => boundaries.push(boundary)
        }
    }

    let boundary: u32 = *if terminal.rightward { boundaries.iter().max() } else { boundaries.iter().min() }?;

    let overhangs: bool = if terminal.rightward {
        terminal.block.reference_end > boundary
    } else {
        terminal.block.reference_start < boundary
    };

    overhangs.then_some(boundary)
}


fn terminal_blocks(blocks: &[TranscriptModelExon]) -> Vec<TerminalBlock<'_>> {
    if blocks.len() < 2 {
        return Vec::new();
    }
    
    let pairs: [(ReadEnd, &TranscriptModelExon, &TranscriptModelExon); 2] = [
        (ReadEnd::Head, &blocks[0], &blocks[1]),
        (ReadEnd::Tail, &blocks[blocks.len() - 1], &blocks[blocks.len() - 2])
    ];
    
    pairs.into_iter().filter_map(|(end, block, neighbour)| {
        if block.reference_chromosome_id != neighbour.reference_chromosome_id
            || block.reference_strand != neighbour.reference_strand {
            return None;
        }
        
        let rightward: bool = match block.reference_strand {
            Strand::Forward => end == ReadEnd::Tail,
            Strand::Reverse => end == ReadEnd::Head,
            _ => return None
        };
        
        // A back-splice reverses this order and fails the gap test.
        let (lower, higher): (&TranscriptModelExon, &TranscriptModelExon) = if rightward { 
            (neighbour, block) 
        } else { 
            (block, neighbour) 
        };
        
        if lower.reference_end + 1 >= higher.reference_start {
            return None;
        }
        
        Some(TerminalBlock {
            end: end, 
            block: block, 
            intron: (lower.reference_end + 1, higher.reference_start - 1), 
            rightward: rightward
        })
    }).collect()
}


#[cfg(test)]
#[path = "../tests/correction/transcript_end_trimming.rs"]
mod tests;