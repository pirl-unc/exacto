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


use exacto_core::prelude::*;
use std::collections::BTreeSet;

use crate::prelude::*;


pub fn identify_transcript_model_exons(alignment_model: &AlignmentModel) -> Vec<TranscriptModelExon> {
    // Step 1. Cluster adjacent bases by reference position.
    let is_exonic = |base: &AlignmentModelBase| matches!(
        base.get_kind(),
        AlignmentModelBaseKind::Match |
        AlignmentModelBaseKind::Mismatch |
        AlignmentModelBaseKind::Insertion
    );
    let mut uf: UnionFind = UnionFind::new();
    for i in 0..alignment_model.num_bases() {
        let curr_base: &AlignmentModelBase = alignment_model.get_base(i);
        if is_exonic(curr_base) == false {
            continue;
        }
        uf.union(i, i);

        if i == 0 {
            continue;
        }

        let prev_base: &AlignmentModelBase = alignment_model.get_base(i - 1);
        if is_exonic(prev_base) == false {
            continue;
        }

        // The event carrying this boundary is not always keyed on the two adjacent read
        // positions. An event next to an inserted run is stated over the aligned bases
        // either side of the run, so its key straddles the run and is found from either end.
        let event: Option<&AlignmentModelEvent> = alignment_model
            .get_event(i - 1, i)
            .or_else(|| alignment_model.get_event_at(i - 1))
            .or_else(|| alignment_model.get_event_at(i))
            .filter(|event| {
                event.get_prev_read_position() <= i - 1 && event.get_next_read_position() >= i
            });

        // The two bases continue one exon when they are on the same chromosome and strand
        // and at most one reference position apart.
        let continues_reference: bool = match (
            prev_base.get_placement().get_coordinate(),
            curr_base.get_placement().get_coordinate()
        ) {
            (
                Some((prev_chromosome_id, prev_position, prev_strand)),
                Some((curr_chromosome_id, curr_position, curr_strand))
            ) => prev_chromosome_id == curr_chromosome_id
                && prev_strand == curr_strand
                && prev_position.abs_diff(curr_position) <= 1,
            _ => false
        };

        let should_union: bool = match event {
            // Inserted bases are placed at the reference base to their left, so beside an
            // inserted run the event lies where the placement jumps. Where it does not
            // jump, the inserted bases continue the exon of the base they are placed at.
            Some(event) if (event.get_prev_read_position(), event.get_next_read_position()) != (i - 1, i)
                && continues_reference => true,
            // An event separates the two bases. Only a deletion keeps them inside
            // one exon, and only when both sit on the same chromosome and strand.
            Some(event) => event.get_kind() == &AlignmentModelEventKind::Deletion
                && match (prev_base.get_placement(), curr_base.get_placement()) {
                (
                    AlignmentModelBasePlacement::Placed {
                        chromosome_id: prev_chromosome_id,
                        strand: prev_strand,
                        ..
                    },
                    AlignmentModelBasePlacement::Placed {
                        chromosome_id: curr_chromosome_id,
                        strand: curr_strand,
                        ..
                    }
                ) => prev_chromosome_id == curr_chromosome_id
                    && prev_strand == curr_strand,
                (AlignmentModelBasePlacement::Unplaced, AlignmentModelBasePlacement::Unplaced) => true,
                _ => false
            },
            None => continues_reference
        };

        if should_union {
            uf.union(i - 1, i);
        }
    }

    // Step 2. Identify exonic boundaries.
    let mut exons: Vec<TranscriptModelExon> = Vec::new();
    let clusters: Vec<BTreeSet<ReadPosition>> = uf.get_clusters();
    for cluster in clusters.iter() {
        // Sort the read positions.
        let mut read_positions: Vec<ReadPosition> = cluster.iter().map(|&pos| pos).collect();
        read_positions.sort();

        let mut bases: Vec<&AlignmentModelBase> = Vec::new();
        for read_position in read_positions.iter() {
            bases.push(alignment_model.get_base(*read_position));
        }

        // A cluster is homogeneous in placement: Step 1 unions two placed bases
        // or two unplaced ones, never one of each. So the first base settles it
        // for the whole cluster.
        let Some((reference_chromosome_id, _, reference_strand)) =
            bases.first().unwrap().get_placement().get_coordinate()
        else {
            continue;
        };

        let reference_strand: Strand = reference_strand.clone();

        // The bounds come from the matched and mismatched bases. An inserted base is placed
        // at the reference base to its left, which after an intron is the intron's last
        // base; a cluster of inserted bases alone keeps their placement.
        let is_aligned = |base: &AlignmentModelBase| matches!(
            base.get_kind(),
            AlignmentModelBaseKind::Match | AlignmentModelBaseKind::Mismatch
        );
        let bounding_bases: Vec<&AlignmentModelBase> = if bases.iter().any(|base| is_aligned(base)) {
            bases.iter().copied().filter(|base| is_aligned(base)).collect()
        } else {
            bases.clone()
        };
        let reference_positions: Vec<ReferencePosition> = bounding_bases.iter()
            .filter_map(|base| base.get_placement().get_coordinate())
            .map(|(_, reference_position, _)| reference_position)
            .collect();
        let reference_start: ReferencePosition = reference_positions.iter().copied().min().unwrap();
        let reference_end: ReferencePosition = reference_positions.iter().copied().max().unwrap();
        let read_start_position: ReadPosition = bases.first().unwrap().get_read_position();
        let read_end_position: ReadPosition = bases.last().unwrap().get_read_position();

        // Use 0 as a placeholder - exon_number assigned after sorting.
        let exon: TranscriptModelExon = TranscriptModelExon::new(
            reference_chromosome_id,
            reference_start,
            reference_end,
            reference_strand,
            0,
            read_start_position,
            read_end_position
        );

        exons.push(exon);
    }

    // Step 3. Sort by read order, then assign exon numbers.
    exons.sort_by_key(|e| e.read_start_position);
    for (i, exon) in exons.iter_mut().enumerate() {
        exon.exon_number = (i + 1) as ExonNumber;
    }

    exons
}