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
use rayon::iter::ParallelIterator;
use std::collections::HashSet;
use std::ops::Range;

use crate::common::enums::BaseAction;
use crate::correction::reference_column::{identify_reference_columns, ReferenceColumn};
use crate::correction::rna_read_correction_plan::*;
use crate::correction::transcript_end_trimming::identify_transcript_end_trims;
use crate::correction::trusted_alleles::TrustedAlleles;
use crate::prelude::CorrectedRNARead;


pub(crate) fn correct_rna_read(
    cluster_id: usize,
    read_name: &str,
    alignment_model: &AlignmentModel,
    cluster_graph_operations: &HashSet<&GraphOperation>,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    gene_annotator: Option<&(dyn GeneAnnotator + Sync)>,
    max_correctable_event_len: usize,
    corrected_base_quality: u8
) -> CorrectedRNARead {
    // Step 1. Identify the trimmed ranges of the AlignmentModel, if an annotation is given.
    let transcript_trims: Vec<Range<u32>> = match gene_annotator {
        Some(gene_annotator) => identify_transcript_end_trims(
            alignment_model,
            &identify_transcript_model_exons(alignment_model),
            chromosome_names_map,
            gene_annotator,
            cluster_graph_operations
        ),
        None => Vec::new()
    };

    // Step 2. Plan RNA read correction.
    let plan: RNAReadCorrectionPlan = plan_rna_read_correction(
        alignment_model,
        cluster_graph_operations,
        &transcript_trims,
        max_correctable_event_len,
        corrected_base_quality,
    );

    // Step 3. Emit the corrected sequence and base quality scores.
    let (sequence, corrected_qualities) = plan.emit(alignment_model);

    let corrected_sequence: Box<str> = String::from_utf8(sequence)
        .map_err(|error| error.to_string())
        .unwrap()
        .into_boxed_str();

    CorrectedRNARead::new(
        cluster_id,
        read_name.into(),
        corrected_sequence,
        corrected_qualities
    )
}


fn insertion_calls(
    model: &AlignmentModel,
    range: Range<usize>,
    sequence: &str,
    quality: u8
) -> Vec<BaseCall> {
    let observed: &[AlignmentModelBase] = &model.get_bases()[range];

    let same_allele: bool = observed.len() == sequence.len()
        && observed.iter().zip(sequence.bytes()).all(|(base, nucleotide)| {
        observed_call(base)
            .nucleotide
            .eq_ignore_ascii_case(&nucleotide)
    });

    if same_allele {
        observed.iter().map(observed_call).collect()
    } else {
        sequence
            .bytes()
            .map(|nucleotide| BaseCall {
                nucleotide,
                quality,
            })
            .collect()
    }
}


fn plan_rna_read_correction(
    model: &AlignmentModel,
    operations: &HashSet<&GraphOperation>,
    transcript_trims: &[Range<u32>],
    max_event_len: usize,
    quality: u8
) -> RNAReadCorrectionPlan {
    let mut plan: RNAReadCorrectionPlan = RNAReadCorrectionPlan::new(model);

    // Step 1. Get trusted alleles.
    let mut trusted: TrustedAlleles = TrustedAlleles::new(operations);

    // Step 2. Create a trim mask.
    // The trim mask distinguishes intentional trimming from other Skip actions.
    // Every original base still has its own emission row.
    let mut trimmed: Vec<bool> = vec![false; model.get_bases().len()];
    for range in transcript_trims {
        let range = range.start as usize..range.end as usize;
        trimmed[range.clone()].fill(true);
        plan.skip(range);
    }

    // Step 3. Plan terminal soft clips: replace a called one, keep the others.
    let clip_anchors: Vec<(u16, u32)> = plan_terminal_clips(
        model,
        &mut plan,
        &trimmed,
        &trusted,
        quality
    );

    // Step 4. Identify reference-coordinate columns for RNA read correction, in read order.
    // One column per reference position the read covers, including bases the read lacks
    // because of a deletion no longer than `max_event_len`.
    let (reference_columns, long_deletions): (Vec<ReferenceColumn>, Vec<(u16, u32, u32)>) = identify_reference_columns(
        model,
        &mut plan,
        &trimmed,
        max_event_len
    );

    // The read's bases between a column and the one before it, when the two are consecutive
    // reference positions: (left anchor, first read position, end read position). Either column
    // may be a restored deletion base.
    let gap = |index: usize| -> Option<(u32, usize, usize)> {
        let previous: &ReferenceColumn = &reference_columns[index.checked_sub(1)?];
        let column: &ReferenceColumn = &reference_columns[index];
        let consecutive: bool =
            previous.chromosome == column.chromosome
                && previous.strand == column.strand
                && match column.strand {
                Strand::Forward => {
                    previous.position.checked_add(1) == Some(column.position)
                }
                Strand::Reverse => {
                    column.position.checked_add(1) == Some(previous.position)
                }
                _ => false
            };
        let start: usize = previous
            .original
            .map_or(previous.boundary, |position| position + 1);
        consecutive.then_some((previous.position.min(column.position), start, column.boundary))
    };
    let ordinary_insertion = |start: usize, end: usize| -> bool {
        let in_junction: bool =
            model.get_events().values().any(|event| {
                *event.get_kind() != AlignmentModelEventKind::Deletion
                    && (event.get_prev_read_position() as usize) < end
                    && (event.get_next_read_position() as usize)
                    >= start
            });

        start <= end
            && !in_junction
            && model.get_bases()[start..end].iter().all(|base| {
            *base.get_kind() == AlignmentModelBaseKind::Insertion
        })
    };

    // Step 5. Aligners place a long insertion or deletion at different copies of its repeat in
    // different reads.
    let mut own_insertions: Vec<(u16, u32)> = clip_anchors;
    for index in 1..reference_columns.len() {
        if let Some((left, start, end)) = gap(index) {
            if end.saturating_sub(start) > max_event_len && ordinary_insertion(start, end) {
                own_insertions.push((reference_columns[index].chromosome, left));
            }
        }
    }
    let own_deletions: Vec<(u16, u32)> = long_deletions
        .iter()
        .map(|&(chromosome, left, _)| (chromosome, left))
        .collect();
    let holds = |own: &[(u16, u32)], chromosome: u16, anchor: u32, length: usize| -> bool {
        own.iter().any(|&(chr, position)| {
            chr == chromosome && position != anchor && position.abs_diff(anchor) as usize <= length
        })
    };
    trusted.insertions.retain(|&(chromosome, anchor), sequence| {
        !holds(&own_insertions, chromosome, anchor, sequence.len())
    });
    trusted.deletions.retain(|&(chromosome, left, right)| {
        !holds(&own_deletions, chromosome, left, (right - left - 1) as usize)
    });

    // Remove observed insertion sequence inside a chosen deletion.
    for (index, base) in model.get_bases().iter().enumerate() {
        if *base.get_kind() != AlignmentModelBaseKind::Insertion {
            continue;
        }

        if let Some((chromosome, anchor, _)) = base.get_placement().get_coordinate() {
            let inside_deletion = trusted.deletions.iter().any(|&(chr, left, right)| {
                chr == chromosome && left <= anchor && anchor < right - 1
            });
            if inside_deletion {
                plan.bases[index].action = BaseAction::Skip;
            }
        }
    }

    for (index, column) in reference_columns.iter().enumerate() {
        // A confident insertion belongs between consecutive reference columns.
        if let Some((left, start, end)) = gap(index) {
            if let Some(sequence) = trusted.insertions.get(&(column.chromosome, left)) {
                if ordinary_insertion(start, end) {
                    let oriented: Box<str> = if column.strand == Strand::Reverse {
                        reverse_complement(sequence)
                    } else {
                        sequence.clone()
                    };
                    let calls: Vec<BaseCall> = insertion_calls(model, start..end, &oriented, quality);
                    plan.skip(start..end);
                    plan.insert_before(column.boundary, calls);
                }
            }
        }

        // Apply a confident deletion even when this read contains a match here.
        if trusted.deletes(column.chromosome, column.position) {
            if let Some(original) = column.original {
                plan.bases[original].action = BaseAction::Skip;
            }

            // A virtual reference column is simply not inserted.
            continue;
        }

        // Choose the confident substitution, otherwise the reference nucleotide.
        let nucleotide = if let Some(&allele) = trusted
            .substitutions
            .get(&(column.chromosome, column.position)) {
            if column.strand == Strand::Reverse {
                reverse_complement(
                    std::str::from_utf8(&[allele]).expect("ASCII allele"),
                )
                    .as_bytes()[0]
            } else {
                allele
            }
        } else {
            column.reference
        };

        if let Some(original) = column.original {
            let observed: BaseCall = observed_call(model.get_base(original as u32));
            plan.bases[original].action = if observed.nucleotide.eq_ignore_ascii_case(&nucleotide) {
                BaseAction::Keep
            } else {
                BaseAction::Replace(BaseCall {
                    nucleotide,
                    quality,
                })
            };
        } else {
            plan.insert_before(
                column.boundary,
                [BaseCall {
                    nucleotide,
                    quality,
                }],
            );
        }
    }

    plan
}


fn plan_terminal_clips(
    model: &AlignmentModel,
    plan: &mut RNAReadCorrectionPlan,
    trimmed: &[bool],
    trusted: &TrustedAlleles,
    quality: u8
) -> Vec<(u16, u32)> {
    let mut kept: Vec<(u16, u32)> = Vec::new();
    for range in model.terminal_soft_clip_ranges() {
        let leading: bool = range.start == 0;

        let anchor: &AlignmentModelBase = model.get_base(if leading {
            range.end
        } else {
            range.start - 1
        });

        let Some((chromosome, position, strand)) =
            anchor.get_placement().get_coordinate()
        else {
            continue;
        };

        if !matches!(strand, Strand::Forward | Strand::Reverse) {
            continue;
        }

        let left: u32 = match GraphOperationType::for_breakpoint(strand, leading) {
            GraphOperationType::Upstream => {
                let Some(left) = position.checked_sub(1) else {
                    continue;
                };
                left
            },
            _ => position
        };

        let range: Range<usize> = range.start as usize..range.end as usize;

        if trimmed[range.clone()].iter().any(|&trimmed| trimmed) {
            continue;
        }

        // A clip is bases the aligner did not place, not a claim the variant list rules on, so
        // an uncalled clip is kept as sequenced. A clip the cluster called as an insertion at
        // its anchor is written in the cluster's spelling.
        if let Some(sequence) = trusted.insertions.get(&(chromosome, left)) {
            let oriented = if *strand == Strand::Reverse {
                reverse_complement(sequence)
            } else {
                sequence.clone()
            };

            let calls = insertion_calls(model, range.clone(), &oriented, quality);
            plan.skip(range.clone());
            plan.insert_before(range.start, calls);
        } else {
            kept.push((chromosome, left));
        }
    }
    kept
}


#[cfg(test)]
#[path = "../tests/correction/rna_read_correction.rs"]
mod tests;