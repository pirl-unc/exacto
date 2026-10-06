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


use exacto_caller::prelude::*;
use exacto_core::prelude::*;
use std::collections::HashMap;

use crate::correction::rna_read_correction_plan::RNAReadCorrectionPlan;


pub(crate) struct ReferenceColumn {
    pub chromosome: u16,
    pub position: u32,
    pub strand: Strand,

    /// Some for an original aligned base; None for a missing reference base.
    pub original: Option<usize>,

    /// Original read boundary where a missing base would be emitted.
    pub boundary: usize,

    /// Nucleotide in read orientation.
    pub reference: u8
}


pub(crate) fn identify_reference_columns(
    model: &AlignmentModel,
    plan: &mut RNAReadCorrectionPlan,
    trimmed: &[bool],
    max_event_len: usize
) -> (Vec<ReferenceColumn>, Vec<(u16, u32, u32)>) {
    let mut restored: HashMap<usize, Vec<ReferenceColumn>> = HashMap::new();
    let mut long_deletions: Vec<(u16, u32, u32)> = Vec::new();
    for record in identify_alignment_model_records(model, None) {
        let start = record.get_start() as usize;
        let end = record.get_end() as usize;
        match record.get_kind() {
            AlignmentModelKind::Base(AlignmentModelBaseKind::Insertion) => {
                // Junction sequence is not an ordinary insertion to erase.
                let in_junction = model.get_events().values().any(|event| {
                    *event.get_kind() != AlignmentModelEventKind::Deletion
                        && start < event.get_next_read_position() as usize
                        && end >= event.get_prev_read_position() as usize + 1
                });
                if end - start + 1 <= max_event_len && !in_junction {
                    plan.skip(start..end + 1);
                }
            },
            AlignmentModelKind::Event(AlignmentModelEventKind::Deletion) => {
                // Do not restore across a discarded transcript end.
                if trimmed[start] || trimmed[end] {
                    continue;
                }
                let left: u32 = record.get_position_1();
                let right: u32 = record.get_position_2();
                let length: usize = right
                    .checked_sub(left)
                    .and_then(|span| span.checked_sub(1))
                    .unwrap() as usize;
                if length == 0 {
                    continue;
                }
                if length > max_event_len {
                    long_deletions.push((record.get_chromosome_1(), left, right));
                    continue;
                }
                let forward: bool = *record.get_strand_1() == Strand::Forward;
                let flank: usize = if forward { start } else { end };
                let deleted: &[Nucleotide] = model
                    .get_base(flank as u32)
                    .get_deleted_reference_bases();

                // Fill only what the cs tag spelled out.
                if deleted.len() != length {
                    continue;
                }

                // Distinguish I-D from D-I, including reverse-strand alignments.
                let preceding_anchor: u32 = if forward { left } else { right - 1 };
                let mut boundary: usize = start + 1;
                while boundary < end {
                    let base: &AlignmentModelBase = model.get_base(boundary as u32);
                    let precedes_gap: bool =
                        *base.get_kind() == AlignmentModelBaseKind::Insertion
                            && base.get_placement().get_coordinate().is_some_and(
                            |(chromosome, position, strand)| {
                                chromosome == record.get_chromosome_1()
                                    && position == preceding_anchor
                                    && strand == record.get_strand_1()
                            }
                        );
                    if !precedes_gap {
                        break;
                    }
                    boundary += 1;
                }

                let columns = deleted
                    .iter()
                    .enumerate()
                    .map(|(offset, nucleotide)| ReferenceColumn {
                        chromosome: record.get_chromosome_1(),
                        position: if forward {
                            left + 1 + offset as u32
                        } else {
                            right - 1 - offset as u32
                        },
                        strand: record.get_strand_1().clone(),
                        original: None,
                        boundary,
                        reference: nucleotide.as_str().as_bytes()[0],
                    })
                    .collect();

                restored.insert(boundary, columns);
            },
            _ => {}
        }
    }

    let mut columns: Vec<ReferenceColumn> = Vec::new();

    for (index, base) in model.get_bases().iter().enumerate() {
        if let Some(missing) = restored.remove(&index) {
            columns.extend(missing);
        }

        if trimmed[index]
            || !matches!(
                base.get_kind(),
                AlignmentModelBaseKind::Match | AlignmentModelBaseKind::Mismatch
            ) {
            continue;
        }

        let Some((chromosome, position, strand)) =
            base.get_placement().get_coordinate()
        else {
            continue;
        };

        if !matches!(strand, Strand::Forward | Strand::Reverse) {
            continue;
        }

        let reference: &Nucleotide =
            if *base.get_kind() == AlignmentModelBaseKind::Mismatch
                && max_event_len > 0 {
                base.get_reference_nucleotide().unwrap()
            } else {
                base.get_nucleotide()
            };

        columns.push(ReferenceColumn {
            chromosome,
            position,
            strand: strand.clone(),
            original: Some(index),
            boundary: index,
            reference: reference.as_str().as_bytes()[0]
        });
    }

    (columns, long_deletions)
}