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
use bio::alignment::pairwise::banded::Aligner;
use bio::alignment::Alignment;
use exacto_core::prelude::*;
use std::collections::HashMap;

use crate::prelude::*;


pub(crate) fn retype_insertion_as_breakpoints(
    variant_record: &VariantRecord,
    query_positions: &HashMap<ReferenceChromosomeID, Vec<ReferencePosition>>,
    fasta_map: &FastaMap,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    min_query_length: u32,
    max_query_length: u32,
    search_distance: u32,
    gap_open: i32,
    gap_extend: i32,
    k: u32,
    band_width: u32,
    min_score_fraction: f64,
    min_query_coverage: f64,
    min_span_proportion: f64
) -> Option<Vec<VariantRecord>> {
    assert_eq!(*variant_record.get_variant_type(), VariantType::Insertion);
    assert_eq!(variant_record.get_chromosome_1(), variant_record.get_chromosome_2());
    assert!(!query_positions.is_empty());
    assert!(min_query_length > 0);
    assert!(max_query_length >= min_query_length);
    assert!(k > 0);
    assert!((0.0..=1.0).contains(&min_score_fraction));
    assert!((0.0..=1.0).contains(&min_query_coverage));

    // Step 1. Check if the insertion sequence is long enough to attempt to align.
    let sequence: String = variant_record.get_standardized_sequence();
    if min_query_length > sequence.len() as u32 {
        return None;
    }
    let anchor_chromosome_id: ReferenceChromosomeID = variant_record.get_chromosome_1();
    let Some(chromosome) = chromosome_names_map.get_by_right(&anchor_chromosome_id) else {
        return None;
    };

    // Step 2. Normalize the insertion locus.
    let (anchor_position_1, anchor_position_2, sequence): (ReferencePosition, ReferencePosition, String) =
        normalize_insertion_locus(
            sequence,
            variant_record.get_position_1(),
            variant_record.get_position_2(),
            chromosome,
            fasta_map
        );
    let sequence_length: u32 = sequence.len() as u32;
    let query_length: u32 = sequence_length.min(max_query_length);

    // Step 3. Build reference windows.
    let reference_windows: Vec<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)> = build_reference_windows(
        query_positions,
        search_distance
    );
    if reference_windows.is_empty() {
        return None;
    }

    // Step 4. Place the allele.
    let whole_placement: Option<QueryPlacement> = if sequence_length <= max_query_length {
        place_query(
            sequence.as_str(),
            &reference_windows,
            fasta_map,
            chromosome_names_map,
            gap_open,
            gap_extend,
            k,
            band_width,
            min_score_fraction,
            min_query_coverage
        )
    } else {
        None
    };
    let (entry, exit): (WalkBreakend, WalkBreakend) = match whole_placement {
        Some(placement) => (segment_entry(&placement), segment_exit(&placement)),
        None => {
            let front_query: String = sequence[..query_length as usize].to_ascii_uppercase();
            let back_query: String =
                sequence[sequence_length as usize - query_length as usize..].to_ascii_uppercase();
            let front: QueryPlacement = place_query(
                front_query.as_str(),
                &reference_windows,
                fasta_map,
                chromosome_names_map,
                gap_open,
                gap_extend,
                k,
                band_width,
                min_score_fraction,
                min_query_coverage
            )?;
            let back: QueryPlacement = place_query(
                back_query.as_str(),
                &reference_windows,
                fasta_map,
                chromosome_names_map,
                gap_open,
                gap_extend,
                k,
                band_width,
                min_score_fraction,
                min_query_coverage
            )?;
            match merge_end_placements(&front, &back, sequence_length, min_span_proportion) {
                Some(merged) => (segment_entry(&merged), segment_exit(&merged)),
                None => (segment_entry(&front), segment_exit(&back))
            }
        }
    };

    // Step 5. Read the two candidate adjacencies off the walk and emit the novel,
    // representable ones.
    let anchor_exit: WalkBreakend = WalkBreakend {
        chromosome_id: anchor_chromosome_id,
        position: anchor_position_1,
        strand: variant_record.get_strand_1().clone(),
        operation: variant_record.get_operation_1().clone()
    };
    let anchor_entry: WalkBreakend = WalkBreakend {
        chromosome_id: anchor_chromosome_id,
        position: anchor_position_2,
        strand: variant_record.get_strand_2().clone(),
        operation: variant_record.get_operation_2().clone()
    };

    let mut retyped_variant_records: Vec<VariantRecord> = Vec::new();
    if !is_reference_continuation(&anchor_exit, &entry) {
        let record: VariantRecord = build_adjacency_record(
            variant_record,
            &anchor_exit,
            &entry,
            variant_record.get_read_position_1(),
            variant_record.get_read_position_1() + 1
        );
        if !is_deletion_shaped(record.get_graph_operation()) {
            retyped_variant_records.push(record);
        }
    }
    if !is_reference_continuation(&exit, &anchor_entry) {
        let record: VariantRecord = build_adjacency_record(
            variant_record,
            &exit,
            &anchor_entry,
            variant_record.get_read_position_2() - 1,
            variant_record.get_read_position_2()
        );
        if !is_deletion_shaped(record.get_graph_operation()) {
            retyped_variant_records.push(record);
        }
    }

    // Step 6. Compose a residual duplication cycle.
    if retyped_variant_records.len() == 2 {
        let front_operation: &GraphOperation = retyped_variant_records[0].get_graph_operation();
        let back_operation: &GraphOperation = retyped_variant_records[1].get_graph_operation();
        let same_chromosome: bool =
            front_operation.get_chromosome_1() == front_operation.get_chromosome_2()
                && front_operation.get_chromosome_1() == back_operation.get_chromosome_1()
                && front_operation.get_chromosome_1() == back_operation.get_chromosome_2();
        let cycle_shaped: bool =
            *front_operation.get_operation_type_1() == GraphOperationType::Upstream
                && *front_operation.get_operation_type_2() == GraphOperationType::Downstream
                && *back_operation.get_operation_type_1() == GraphOperationType::Upstream
                && *back_operation.get_operation_type_2() == GraphOperationType::Downstream;
        let meets_at_anchors: bool =
            front_operation.get_position_2() == anchor_position_1
                && back_operation.get_position_1() == anchor_position_2
                && anchor_position_1.checked_add(1) == Some(anchor_position_2);
        if same_chromosome && cycle_shaped && meets_at_anchors {
            let outer_span: u32 =
                back_operation.get_position_2() - front_operation.get_position_1() + 1;
            let shorter: u32 = outer_span.min(sequence_length);
            let longer: u32 = outer_span.max(sequence_length);
            if (shorter as f64) >= min_span_proportion * longer as f64 {
                let merged_operation: GraphOperation = GraphOperation::new(
                    front_operation.get_chromosome_1(),
                    front_operation.get_position_1(),
                    front_operation.get_strand_1().clone(),
                    front_operation.get_operation_type_1().clone(),
                    back_operation.get_chromosome_2(),
                    back_operation.get_position_2(),
                    back_operation.get_strand_2().clone(),
                    back_operation.get_operation_type_2().clone(),
                    "".into(),
                    VariantType::Breakpoint
                );
                return Some(vec![VariantRecord::new(
                    variant_record.get_read_id(),
                    variant_record.get_read_position_1(),
                    variant_record.get_read_position_2(),
                    merged_operation
                )]);
            }
        }
    }

    if retyped_variant_records.is_empty() {
        None
    } else {
        Some(retyped_variant_records)
    }
}


#[derive(Clone, Debug)]
struct QueryPlacement {
    chromosome_id: ReferenceChromosomeID,
    low: ReferencePosition,
    high: ReferencePosition,
    strand: Strand,
    score: i32,
    aligned_query_length: usize
}


#[derive(Clone, Debug)]
struct WalkBreakend {
    chromosome_id: ReferenceChromosomeID,
    position: ReferencePosition,
    strand: Strand,
    operation: GraphOperationType
}


fn build_adjacency_record(
    variant_record: &VariantRecord,
    exit: &WalkBreakend,
    entry: &WalkBreakend,
    read_position_1: ReadPosition,
    read_position_2: ReadPosition
) -> VariantRecord {
    let variant_type: VariantType = if exit.chromosome_id == entry.chromosome_id {
        VariantType::Breakpoint
    } else {
        VariantType::Translocation
    };
    VariantRecord::new(
        variant_record.get_read_id(),
        read_position_1,
        read_position_2,
        GraphOperation::new(
            exit.chromosome_id,
            exit.position,
            exit.strand.clone(),
            exit.operation.clone(),
            entry.chromosome_id,
            entry.position,
            entry.strand.clone(),
            entry.operation.clone(),
            "".into(),
            variant_type
        )
    )
}


fn build_reference_windows(
    query_positions: &HashMap<ReferenceChromosomeID, Vec<ReferencePosition>>,
    search_distance: u32
) -> Vec<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)> {
    let mut windows: Vec<(ReferenceChromosomeID, ReferencePosition, ReferencePosition)> = Vec::new();
    for &chromosome_id in query_positions.keys() {
        let ranges: Vec<(isize, isize)> = query_positions[&chromosome_id]
            .iter()
            .map(|position| {
                let start: ReferencePosition = position.saturating_sub(search_distance).max(1);
                let end: ReferencePosition = position.saturating_add(search_distance).max(start);
                (start as isize, end as isize)
            })
            .collect();
        windows.extend(
            merge_regions(ranges)
                .into_iter()
                .map(|(start, end)| (chromosome_id, start as ReferencePosition, end as ReferencePosition))
        );
    }
    windows
}


fn is_reference_continuation(exit: &WalkBreakend, entry: &WalkBreakend) -> bool {
    exit.chromosome_id == entry.chromosome_id
        && exit.operation == GraphOperationType::Downstream
        && entry.operation == GraphOperationType::Upstream
        && exit.position.checked_add(1) == Some(entry.position)
}


fn is_deletion_shaped(operation: &GraphOperation) -> bool {
    operation.get_chromosome_1() == operation.get_chromosome_2()
        && *operation.get_operation_type_1() == GraphOperationType::Downstream
        && *operation.get_operation_type_2() == GraphOperationType::Upstream
}


fn merge_end_placements(
    front: &QueryPlacement,
    back: &QueryPlacement,
    sequence_length: u32,
    min_span_proportion: f64
) -> Option<QueryPlacement> {
    if front.chromosome_id != back.chromosome_id || front.strand != back.strand {
        return None;
    }
    let low: ReferencePosition = front.low.min(back.low);
    let high: ReferencePosition = front.high.max(back.high);
    let span: u32 = high - low + 1;
    let shorter: u32 = span.min(sequence_length);
    let longer: u32 = span.max(sequence_length);
    if (shorter as f64) < min_span_proportion * longer as f64 {
        return None;
    }
    Some(QueryPlacement {
        chromosome_id: front.chromosome_id,
        low,
        high,
        strand: front.strand.clone(),
        score: front.score + back.score,
        aligned_query_length: front.aligned_query_length + back.aligned_query_length
    })
}


fn normalize_insertion_locus(
    sequence: String,
    anchor_position_1: ReferencePosition,
    anchor_position_2: ReferencePosition,
    chromosome: &str,
    fasta_map: &FastaMap
) -> (ReferencePosition, ReferencePosition, String) {
    let sequence_length: usize = sequence.len();
    let flank_start: usize = anchor_position_2 as usize;
    let flank_end: usize = flank_start + sequence_length - 1;
    let Some(flank) = fasta_map.try_get_sequence(chromosome, flank_start, flank_end) else {
        return (anchor_position_1, anchor_position_2, sequence);
    };
    let flank: String = flank.to_ascii_uppercase();

    let sequence_bytes: &[u8] = sequence.as_bytes();
    let flank_bytes: &[u8] = flank.as_bytes();
    let mut shift: usize = 0;
    while shift < sequence_length
        && shift < flank_bytes.len()
        && sequence_bytes[shift] == flank_bytes[shift] {
        shift += 1;
    }
    if shift == 0 {
        return (anchor_position_1, anchor_position_2, sequence);
    }

    let rotated: String = format!("{}{}", &sequence[shift..], &sequence[..shift]);
    (
        anchor_position_1 + shift as u32,
        anchor_position_2 + shift as u32,
        rotated
    )
}


fn place_query(
    query: &str,
    reference_windows: &[(ReferenceChromosomeID, ReferencePosition, ReferencePosition)],
    fasta_map: &FastaMap,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    gap_open: i32,
    gap_extend: i32,
    k: u32,
    band_width: u32,
    min_score_fraction: f64,
    min_query_coverage: f64
) -> Option<QueryPlacement> {
    let query_reversed: String = reverse_complement(query).to_string();

    let mut best: Option<QueryPlacement> = None;
    let mut best_is_ambiguous: bool = false;
    let score_fn = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

    for reference_window in reference_windows {
        let Some(chromosome) = chromosome_names_map.get_by_right(&reference_window.0)
        else {
            continue;
        };
        let Some(reference_sequence) = fasta_map.try_get_sequence(
            chromosome,
            reference_window.1 as usize,
            reference_window.2 as usize
        ) else {
            continue;
        };
        if k > reference_sequence.len() as u32 {
            continue;
        }
        let reference_sequence: String = reference_sequence.to_ascii_uppercase();

        let mut aligner: Aligner<_> = Aligner::new(gap_open, gap_extend, score_fn, k as usize, band_width as usize);
        let forward: Alignment = aligner.semiglobal(query.as_bytes(), reference_sequence.as_bytes());
        let reverse: Alignment = aligner.semiglobal(query_reversed.as_bytes(), reference_sequence.as_bytes());

        for (alignment, is_forward) in [(forward, true), (reverse, false)] {
            let aligned_query_length: usize = alignment.xend.saturating_sub(alignment.xstart);
            if aligned_query_length == 0
                || (aligned_query_length as f64) < min_query_coverage * query.len() as f64
                || (alignment.score as f64) < min_score_fraction * aligned_query_length as f64 {
                continue;
            }
            let Some(high_offset) = alignment.yend.checked_sub(1) else {
                continue;
            };
            let (Ok(low_offset), Ok(high_offset)) =
                (u32::try_from(alignment.ystart), u32::try_from(high_offset))
            else {
                continue;
            };
            let (Some(low), Some(high)) = (
                reference_window.1.checked_add(low_offset),
                reference_window.1.checked_add(high_offset)
            ) else {
                continue;
            };

            let candidate = QueryPlacement {
                chromosome_id: reference_window.0,
                low,
                high,
                strand: if is_forward { Strand::Forward } else { Strand::Reverse },
                score: alignment.score,
                aligned_query_length
            };

            // The query is aligned end to end, so the score alone ranks the placements, and the
            // verdict does not depend on the order of the windows.
            match &best {
                Some(current) if candidate.score < current.score => {},
                Some(current) if candidate.score == current.score => {
                    if !same_placement(&candidate, current) {
                        best_is_ambiguous = true;
                    }
                },
                _ => {
                    best = Some(candidate);
                    best_is_ambiguous = false;
                }
            }
        }
    }

    if best_is_ambiguous {
        None
    } else {
        best
    }
}


fn same_placement(a: &QueryPlacement, b: &QueryPlacement) -> bool {
    a.chromosome_id == b.chromosome_id
        && a.low == b.low
        && a.high == b.high
        && a.strand == b.strand
}


fn segment_entry(placement: &QueryPlacement) -> WalkBreakend {
    let operation: GraphOperationType =
        GraphOperationType::for_breakpoint(&placement.strand, true);
    let position: ReferencePosition = match operation {
        GraphOperationType::Upstream => placement.low,
        GraphOperationType::Downstream => placement.high,
        _ => unreachable!("for_breakend must return Upstream or Downstream")
    };
    WalkBreakend {
        chromosome_id: placement.chromosome_id,
        position,
        strand: placement.strand.clone(),
        operation
    }
}


fn segment_exit(placement: &QueryPlacement) -> WalkBreakend {
    let operation: GraphOperationType =
        GraphOperationType::for_breakpoint(&placement.strand, false);
    let position: ReferencePosition = match operation {
        GraphOperationType::Upstream => placement.low,
        GraphOperationType::Downstream => placement.high,
        _ => unreachable!("for_breakend must return Upstream or Downstream")
    };
    WalkBreakend {
        chromosome_id: placement.chromosome_id,
        position,
        strand: placement.strand.clone(),
        operation
    }
}


#[cfg(test)]
#[path = "../../tests/calling/dna/breakpoint_rescue.rs"]
mod tests;