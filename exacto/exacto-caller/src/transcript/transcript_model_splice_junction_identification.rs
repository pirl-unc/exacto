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

use crate::prelude::*;


pub fn identify_transcript_model_splice_junctions(
    alignment_model: &AlignmentModel
) -> Vec<TranscriptModelSpliceJunction> {
    let mut splice_junctions: Vec<TranscriptModelSpliceJunction> = Vec::new();
    for ((read_position_1, read_position_2), event) in alignment_model.get_events() {
        match event.get_kind() {
            AlignmentModelEventKind::Splicing => {
                let base_1: &AlignmentModelBase = alignment_model.get_base(*read_position_1);
                let base_2: &AlignmentModelBase = alignment_model.get_base(*read_position_2);

                // A splicing event is only ever recorded between two bases that an
                // alignment record placed, so in a well-formed model both flanks
                // are `Placed` and neither arm below is taken.
                let AlignmentModelBasePlacement::Placed {
                    chromosome_id: chromosome_id_1,
                    position: position_1,
                    strand: strand_1,
                    ..
                } = base_1.get_placement() else {
                    continue;
                };
                let AlignmentModelBasePlacement::Placed {
                    chromosome_id: chromosome_id_2,
                    position: position_2,
                    strand: strand_2,
                    ..
                } = base_2.get_placement() else {
                    continue;
                };

                let reference_chromosome_id_1: ReferenceChromosomeID = *chromosome_id_1;
                let reference_chromosome_id_2: ReferenceChromosomeID = *chromosome_id_2;
                let reference_strand_1: Strand = strand_1.clone();
                let reference_strand_2: Strand = strand_2.clone();

                // Determine splice junctions positions 1 and 2.
                let reference_position_1: ReferencePosition = if reference_strand_1 == Strand::Forward {
                    *position_1 + 1
                } else {
                    *position_1 - 1
                };
                let reference_position_2: ReferencePosition = if reference_strand_2 == Strand::Forward {
                    *position_2 - 1
                } else {
                    *position_2 + 1
                };

                // Use 0 as a placeholder - splice_juncion_number assigned after sorting.
                let splice_junction: TranscriptModelSpliceJunction = TranscriptModelSpliceJunction::new(
                    reference_chromosome_id_1,
                    reference_chromosome_id_2,
                    reference_position_1,
                    reference_position_2,
                    reference_strand_1,
                    reference_strand_2,
                    0,
                    *read_position_1,
                    *read_position_2
                );

                splice_junctions.push(splice_junction);
            },
            _ => {
                // Do nothing.
            }
        }
    }

    // Step 2. Sort by read order, then assign splice junction numbers.
    splice_junctions.sort_by_key(|e| e.read_position_1);
    for (i, splice_junction) in splice_junctions.iter_mut().enumerate() {
        splice_junction.splice_junction_number = (i + 1) as SpliceJunctionNumber;
    }

    splice_junctions
}
