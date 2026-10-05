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


fn base_context<'a>(
    annotation: Option<&'a TranscriptModelAnnotation>, 
    read_position: ReadPosition
) -> Option<&'a AlignmentModelBaseContext> {
    annotation?.get_base_context(read_position)
}

fn base_annotation<'a>(
    annotation: Option<&'a TranscriptModelAnnotation>, 
    read_position: ReadPosition
) -> Option<&'a TranscriptModelBaseAnnotation> {
    annotation?.get_base(read_position)
}

fn event_annotation<'a>(
    annotation: Option<&'a TranscriptModelAnnotation>, 
    read_position_1: ReadPosition,
    read_position_2: ReadPosition
) -> Option<&'a TranscriptModelEventAnnotation> {
    annotation?.get_event(read_position_1, read_position_2)
}


pub fn identify_alignment_model_records(
    model: &AlignmentModel,
    annotation: Option<&TranscriptModelAnnotation>
) -> Vec<AlignmentModelRecord> {
    // Step 1. Cluster Match, Mismatch, Insertion bases.
    let num_bases: u32 = model.num_bases();
    let mut uf_bases: UnionFind = UnionFind::new();
    for i in 0..num_bases {
        if matches!(model.get_base(i).get_kind(),
            AlignmentModelBaseKind::Match |
            AlignmentModelBaseKind::Mismatch |
            AlignmentModelBaseKind::Insertion) {
            uf_bases.union(i, i);
        }
    }
    for i in 0..num_bases {
        if i > 0 {
            let prev_base: &AlignmentModelBase = model.get_base(i - 1);
            let curr_base: &AlignmentModelBase = model.get_base(i);
            let is_prev_base_aligned: bool = matches!(prev_base.get_kind(),
                AlignmentModelBaseKind::Match |
                AlignmentModelBaseKind::Mismatch |
                AlignmentModelBaseKind::Insertion
            );
            let is_curr_base_aligned: bool = matches!(curr_base.get_kind(),
                AlignmentModelBaseKind::Match |
                AlignmentModelBaseKind::Mismatch |
                AlignmentModelBaseKind::Insertion
            );
            if is_prev_base_aligned && is_curr_base_aligned {
                // Both bases must be placed on the reference for the coordinate
                // comparison below to mean anything.
                let Some((prev_chromosome_id, prev_position, prev_strand)) = prev_base
                    .get_placement()
                    .get_coordinate()
                else {
                    continue;
                };
                let Some((curr_chromosome_id, curr_position, curr_strand)) = curr_base
                    .get_placement()
                    .get_coordinate()
                else {
                    continue;
                };
                
                // Absent for every base of a DNA read.
                let prev_context: Option<&AlignmentModelBaseContext> = base_context(
                    annotation,
                    prev_base.get_read_position()
                );
                let curr_context: Option<&AlignmentModelBaseContext> = base_context(
                    annotation,
                    curr_base.get_read_position()
                );

                match (prev_context, curr_context) {
                    (Some(prev_context), Some(curr_context)) => {
                        // Both bases are annotated: the contexts must agree too,
                        // so an exon/intron boundary splits the cluster.
                        if *prev_base.get_kind() == *curr_base.get_kind() &&
                            *prev_context == *curr_context &&
                            prev_chromosome_id == curr_chromosome_id &&
                            prev_position.abs_diff(curr_position) <= 1 &&
                            prev_strand == curr_strand {
                            uf_bases.union(i - 1, i);
                        }
                    },
                    _ => {
                        // At least one of the two bases carries no context:
                        // compare kind and coordinates only.
                        if *prev_base.get_kind() == *curr_base.get_kind() &&
                            prev_chromosome_id == curr_chromosome_id &&
                            prev_position.abs_diff(curr_position) <= 1 &&
                            prev_strand == curr_strand {
                            uf_bases.union(i - 1, i);
                        }
                    }
                }
            }
        }
    }

    // Step 2. Record bases.
    let mut records: Vec<AlignmentModelRecord> = Vec::new();
    for cluster in uf_bases.get_clusters() {
        let mut read_positions: Vec<ReadPosition> = cluster.into_iter().collect();
        read_positions.sort();
        let bases: Vec<&AlignmentModelBase> = read_positions
            .iter()
            .map(|&i| model.get_base(i))
            .collect();
        let first_base: &AlignmentModelBase = bases.first().unwrap();
        let last_base: &AlignmentModelBase = bases.last().unwrap();
        let mut sequence: String = String::new();
        let mut base_quality_scores: Vec<BaseQuality> = Vec::new();
        for base in bases.iter() {
            assert_eq!(
                matches!(base.get_kind(),
                    AlignmentModelBaseKind::Match |
                    AlignmentModelBaseKind::Mismatch |
                    AlignmentModelBaseKind::Insertion
                ),
                true
            );
            sequence.push_str(base.get_nucleotide().as_str());
            base_quality_scores.push(base.get_base_quality());
        }

        // Orientation is read off the first base of the cluster. A cluster of one
        // unplaced base has no reference coordinates and so nothing to record.
        let Some((_, _, first_base_strand)) = first_base
            .get_placement()
            .get_coordinate()
        else {
            continue;
        };
        let (base_1, base_2) = if *first_base_strand == Strand::Forward {
            (first_base, last_base)
        } else {
            (last_base, first_base)
        };

        // Step 1 only ever merges two placed bases, so every base of a multi-base
        // cluster is placed, and a single-base cluster was just checked above.
        let AlignmentModelBasePlacement::Placed {
            chromosome_id: chromosome_id_1,
            position: position_1,
            strand: strand_1,
            mapping_quality: mapping_quality_1
        } = base_1.get_placement() else {
            continue;
        };
        let AlignmentModelBasePlacement::Placed {
            chromosome_id: chromosome_id_2,
            position: position_2,
            strand: strand_2,
            mapping_quality: mapping_quality_2
        } = base_2.get_placement() else {
            continue;
        };

        // Gene, transcript and exon IDs come from the annotation overlay, keyed
        // by read position, instead of from three fields on the base itself.
        // All are None for DNA.
        let base_annotation_1: Option<&TranscriptModelBaseAnnotation> = base_annotation(
            annotation,
            base_1.get_read_position()
        );
        let base_annotation_2: Option<&TranscriptModelBaseAnnotation> = base_annotation(
            annotation,
            base_2.get_read_position()
        );
        let gene_id_1: Option<ReferenceGeneID> = base_annotation_1
            .and_then(|annotation| annotation.reference_gene_id.clone());
        let transcript_id_1: Option<ReferenceTranscriptID> = base_annotation_1
            .and_then(|annotation| annotation.reference_transcript_id.clone());
        let exon_id_1: Option<ReferenceExonID> = base_annotation_1
            .and_then(|annotation| annotation.reference_exon_id.clone());
        let gene_id_2: Option<ReferenceGeneID> = base_annotation_2
            .and_then(|annotation| annotation.reference_gene_id.clone());
        let transcript_id_2: Option<ReferenceTranscriptID> = base_annotation_2
            .and_then(|annotation| annotation.reference_transcript_id.clone());
        let exon_id_2: Option<ReferenceExonID> = base_annotation_2
            .and_then(|annotation| annotation.reference_exon_id.clone());

        // The record's context is the first base's context.
        let context: Option<AlignmentModelContext> = base_context(annotation, first_base.get_read_position())
            .cloned()
            .map(AlignmentModelContext::Base);

        match first_base.get_kind() {
            AlignmentModelBaseKind::Match => {
                let record: AlignmentModelRecord = AlignmentModelRecord::new(
                    first_base.get_read_position(),
                    last_base.get_read_position(),
                    sequence.as_str(),
                    base_quality_scores,
                    AlignmentModelRecordType::Base,
                    AlignmentModelKind::Base(first_base.get_kind().clone()),
                    context,
                    *chromosome_id_1,
                    *position_1,
                    GraphOperationType::Include,
                    strand_1.clone(),
                    *mapping_quality_1,
                    *chromosome_id_2,
                    *position_2,
                    GraphOperationType::Include,
                    strand_2.clone(),
                    *mapping_quality_2,
                    gene_id_1,
                    transcript_id_1,
                    exon_id_1,
                    gene_id_2,
                    transcript_id_2,
                    exon_id_2,
                    None
                );
                records.push(record);
            },
            AlignmentModelBaseKind::Mismatch => {
                let record: AlignmentModelRecord = AlignmentModelRecord::new(
                    first_base.get_read_position(),
                    last_base.get_read_position(),
                    sequence.as_str(),
                    base_quality_scores,
                    AlignmentModelRecordType::Base,
                    AlignmentModelKind::Base(first_base.get_kind().clone()),
                    context,
                    *chromosome_id_1,
                    *position_1 - 1,
                    GraphOperationType::Downstream,
                    strand_1.clone(),
                    *mapping_quality_1,
                    *chromosome_id_2,
                    *position_2 + 1,
                    GraphOperationType::Upstream,
                    strand_2.clone(),
                    *mapping_quality_2,
                    gene_id_1,
                    transcript_id_1,
                    exon_id_1,
                    gene_id_2,
                    transcript_id_2,
                    exon_id_2,
                    None
                );
                records.push(record);
            },
            AlignmentModelBaseKind::Insertion => {
                let record: AlignmentModelRecord = AlignmentModelRecord::new(
                    first_base.get_read_position(),
                    last_base.get_read_position(),
                    sequence.as_str(),
                    base_quality_scores,
                    AlignmentModelRecordType::Base,
                    AlignmentModelKind::Base(first_base.get_kind().clone()),
                    context,
                    *chromosome_id_1,
                    *position_1,
                    GraphOperationType::Downstream,
                    strand_1.clone(),
                    *mapping_quality_1,
                    *chromosome_id_2,
                    *position_2 + 1,
                    GraphOperationType::Upstream,
                    strand_2.clone(),
                    *mapping_quality_2,
                    gene_id_1,
                    transcript_id_1,
                    exon_id_1,
                    gene_id_2,
                    transcript_id_2,
                    exon_id_2,
                    None
                );
                records.push(record);
            },
            _ => {
                // Do nothing.
            }
        }
    }

    // Step 3. Record events.
    for ((read_position_1, read_position_2), event) in model.get_events().iter() {
        let prev_base: &AlignmentModelBase = model.get_base(event.get_prev_read_position());
        let next_base: &AlignmentModelBase = model.get_base(event.get_next_read_position());

        // Get the sequence between the two read positions.
        let mut sequence: String = "".to_string();
        let mut base_quality_scores: Vec<BaseQuality> = Vec::new();
        if read_position_1 < read_position_2 {
            for k in read_position_1 + 1..=read_position_2 - 1 {
                let base: &AlignmentModelBase = model.get_base(k);
                sequence.push_str(base.get_nucleotide().as_str());
                base_quality_scores.push(base.get_base_quality());
            }
        }

        // Both flanks must be placed: the record anchors the event at reference
        // coordinates on either side, and the ordering below is decided by them.
        // An event with an unplaced flank cannot be anchored, so it is skipped.
        let Some((prev_chromosome_id, prev_position, _)) = prev_base.get_placement().get_coordinate() else {
            continue;
        };
        let Some((next_chromosome_id, next_position, _)) = next_base.get_placement().get_coordinate() else {
            continue;
        };

        let (
            base_1,
            base_2,
            operation_1,
            operation_2
        ) = if prev_chromosome_id == next_chromosome_id {
            if prev_position < next_position {
                (
                    prev_base,
                    next_base,
                    event.get_prev_graph_operation_type(),
                    event.get_next_graph_operation_type()
                )
            } else {
                (
                    next_base,
                    prev_base,
                    event.get_next_graph_operation_type(),
                    event.get_prev_graph_operation_type()
                )
            }
        } else {
            (
                prev_base,
                next_base,
                event.get_prev_graph_operation_type(),
                event.get_next_graph_operation_type()
            )
        };

        // `base_1` and `base_2` are `prev_base` and `next_base` in one order or
        // the other, and both were placed above, so both destructures succeed.
        let AlignmentModelBasePlacement::Placed {
            chromosome_id: chromosome_id_1,
            position: position_1,
            strand: strand_1,
            mapping_quality: mapping_quality_1
        } = base_1.get_placement() else {
            continue;
        };
        let AlignmentModelBasePlacement::Placed {
            chromosome_id: chromosome_id_2,
            position: position_2,
            strand: strand_2,
            mapping_quality: mapping_quality_2
        } = base_2.get_placement() else {
            continue;
        };
        
        let event_annotation: Option<&TranscriptModelEventAnnotation> = event_annotation(
            annotation, 
            *read_position_1, 
            *read_position_2
        );
        let context: Option<AlignmentModelContext> = event_annotation
            .and_then(|annotation| annotation.get_context().clone())
            .map(AlignmentModelContext::Event);
        
        let skipped: Option<Vec<Vec<ReferenceBase>>> = Some(
            event_annotation
                .map(|annotation| annotation.get_skipped_reference_bases())
                .unwrap_or_default()
        );

        // Gene, transcript and exon IDs for either flank, from the same overlay.
        let base_annotation_1: Option<&TranscriptModelBaseAnnotation> = base_annotation(annotation, base_1.get_read_position());
        let base_annotation_2: Option<&TranscriptModelBaseAnnotation> = base_annotation(annotation, base_2.get_read_position());

        let record: AlignmentModelRecord = AlignmentModelRecord::new(
            prev_base.get_read_position(),
            next_base.get_read_position(),
            sequence.as_str(),
            base_quality_scores,
            AlignmentModelRecordType::Event,
            AlignmentModelKind::Event(event.get_kind().clone()),
            context,
            *chromosome_id_1,
            *position_1,
            operation_1.clone(),
            strand_1.clone(),
            *mapping_quality_1,
            *chromosome_id_2,
            *position_2,
            operation_2.clone(),
            strand_2.clone(),
            *mapping_quality_2,
            base_annotation_1.and_then(|annotation| annotation.reference_gene_id.clone()),
            base_annotation_1.and_then(|annotation| annotation.reference_transcript_id.clone()),
            base_annotation_1.and_then(|annotation| annotation.reference_exon_id.clone()),
            base_annotation_2.and_then(|annotation| annotation.reference_gene_id.clone()),
            base_annotation_2.and_then(|annotation| annotation.reference_transcript_id.clone()),
            base_annotation_2.and_then(|annotation| annotation.reference_exon_id.clone()),
            skipped
        );
        records.push(record);
    }

    // Step 4. Sort the records.
    records.sort_by(|a, b| {
        a.get_start().cmp(&b.get_start()).then(a.get_end().cmp(&b.get_end()))
    });

    records
}


/// Mismatch and insertion records as variant records. A mismatch run is cut at every base
/// below `min_base_quality`, and each stretch that remains keeps its own positions; an insertion
/// is kept whole or dropped.
///
/// This function identifies the following variant types:
/// - Single-nucleotide variant
/// - Multi-nucleotide variant
/// - Insertion
pub(crate) fn identify_base_variant_records<'a>(
    read_id: ReadID,
    records: impl Iterator<Item = &'a AlignmentModelRecord>,
    min_mapping_quality: MappingQuality,
    min_base_quality: BaseQuality
) -> Vec<VariantRecord> {
    let mut variant_records: Vec<VariantRecord> = Vec::new();
    for curr_record in records {
        if *curr_record.get_record_type() == AlignmentModelRecordType::Event {
            continue;
        }

        // Check the mapping quality scores.
        if curr_record.get_mapping_quality_1() < min_mapping_quality
            || curr_record.get_mapping_quality_2() < min_mapping_quality {
            continue;
        }

        match curr_record.get_kind() {
            AlignmentModelKind::Base(AlignmentModelBaseKind::Mismatch) => {
                // A base below the minimum base quality is removed and the run is cut
                // there. Each stretch that remains is its own record and keeps the
                // positions of its own bases.
                let bases: Vec<char> = curr_record.get_sequence().chars().collect();
                let base_quality_scores: &Vec<BaseQuality> = curr_record.get_base_quality_scores();
                let mut i: usize = 0;
                while i < bases.len() {
                    if base_quality_scores[i] < min_base_quality {
                        i += 1;
                        continue;
                    }
                    let mut j: usize = i;
                    while j + 1 < bases.len() && base_quality_scores[j + 1] >= min_base_quality {
                        j += 1;
                    }

                    let sequence: String = bases[i..=j].iter().collect();
                    let (position_1, position_2): (ReferencePosition, ReferencePosition) = if *curr_record.get_strand_1() == Strand::Forward {
                        (
                            curr_record.get_position_1() + i as u32,
                            curr_record.get_position_1() + j as u32 + 2
                        )
                    } else {
                        (
                            curr_record.get_position_2() - j as u32 - 2,
                            curr_record.get_position_2() - i as u32
                        )
                    };

                    let variant_type = if sequence.len() == 1 {
                        VariantType::SingleNucleotideVariant
                    } else {
                        VariantType::MultiNucleotideVariant
                    };

                    let graph_operation: GraphOperation = GraphOperation::new(
                        curr_record.get_chromosome_1(),
                        position_1,
                        curr_record.get_strand_1().clone(),
                        curr_record.get_operation_1().clone(),
                        curr_record.get_chromosome_2(),
                        position_2,
                        curr_record.get_strand_2().clone(),
                        curr_record.get_operation_2().clone(),
                        sequence.into(),
                        variant_type
                    );

                    variant_records.push(
                        VariantRecord::new(
                            read_id,
                            curr_record.get_start() + i as u32,
                            curr_record.get_start() + j as u32,
                            graph_operation
                        )
                    );

                    i = j + 1;
                }
            },
            AlignmentModelKind::Base(AlignmentModelBaseKind::Insertion) => {
                if !passes_base_quality(curr_record.get_base_quality_scores(), min_base_quality) {
                    continue;
                }
                let sequence: String = curr_record.get_sequence().to_string();

                let graph_operation: GraphOperation = GraphOperation::new(
                    curr_record.get_chromosome_1(),
                    curr_record.get_position_1(),
                    curr_record.get_strand_1().clone(),
                    curr_record.get_operation_1().clone(),
                    curr_record.get_chromosome_2(),
                    curr_record.get_position_2(),
                    curr_record.get_strand_2().clone(),
                    curr_record.get_operation_2().clone(),
                    sequence.into(),
                    VariantType::Insertion
                );

                variant_records.push(
                    VariantRecord::new(
                        read_id,
                        curr_record.get_start(),
                        curr_record.get_end(),
                        graph_operation
                    )
                );
            },
            _ => {
                // Do nothing.
            }
        }
    }

    variant_records
}

/// Identifies event variant records.
///
/// This function identifies the following variant types:
/// - Breakpoint
/// - Circular RNA
/// - Deletion
/// - Exon truncation
/// - Fusion gene
/// - Noncanonical splicing
/// - Translocation
pub(crate) fn identify_event_variant_records(
    read_id: ReadID,
    records: &[AlignmentModelRecord],
    min_mapping_quality: MappingQuality,
    min_base_quality: BaseQuality
) -> Vec<VariantRecord> {
    let mut variant_records: Vec<VariantRecord> = Vec::new();
    for record in records.iter() {
        if *record.get_record_type() == AlignmentModelRecordType::Base {
            continue;
        }

        // Check the mapping quality scores.
        if record.get_mapping_quality_1() < min_mapping_quality
            || record.get_mapping_quality_2() < min_mapping_quality {
            continue;
        }

        match record.get_kind() {
            AlignmentModelKind::Event(AlignmentModelEventKind::Breakpoint) => {
                // The read bases between the two sides, kept whole or not at all.
                let has_passing_sequence: bool = passes_base_quality(record.get_base_quality_scores(), min_base_quality);

                // Get variant type.
                let mut variant_type: VariantType = VariantType::Breakpoint;
                if record.get_context().is_some() {
                    match record.get_context().as_ref().unwrap() {
                        AlignmentModelContext::Event(AlignmentModelEventContext::BackSplicing) => {
                            variant_type = VariantType::CircularRNA;
                        },
                        AlignmentModelContext::Event(AlignmentModelEventContext::FusionGene) => {
                            variant_type = VariantType::FusionGene;
                        },
                        AlignmentModelContext::Event(AlignmentModelEventContext::NonCanonicalSplicing) => {
                            if record.get_chromosome_1() == record.get_chromosome_2() {
                                variant_type = VariantType::Breakpoint;
                            } else {
                                variant_type = VariantType::Translocation;
                            }
                        },
                        _ => {
                            // Do nothing.
                        }
                    }
                } else {
                    if record.get_chromosome_1() == record.get_chromosome_2() {
                        variant_type = VariantType::Breakpoint;
                    } else {
                        variant_type = VariantType::Translocation;
                    }
                }

                let graph_operation: GraphOperation = GraphOperation::new(
                    record.get_chromosome_1(),
                    record.get_position_1(),
                    record.get_strand_1().clone(),
                    record.get_operation_1().clone(),
                    record.get_chromosome_2(),
                    record.get_position_2(),
                    record.get_strand_2().clone(),
                    record.get_operation_2().clone(),
                    record.get_sequence().clone(),
                    variant_type
                );

                if has_passing_sequence {
                    variant_records.push(
                        VariantRecord::new(
                            read_id,
                            record.get_start(),
                            record.get_end(),
                            graph_operation
                        )
                    );
                }

                if record.get_skipped().is_some() {
                    for reference_bases in record.get_skipped().as_ref().unwrap().iter() {
                        let reference_chromosome_id: ReferenceChromosomeID = reference_bases.first().unwrap().reference_chromosome_id;
                        let reference_position_1: ReferencePosition = reference_bases.first().unwrap().reference_position;
                        let reference_position_2: ReferencePosition = reference_bases.last().unwrap().reference_position;
                        let reference_strand: &Strand = &reference_bases.first().unwrap().reference_strand;
                        let sequence: String = if reference_bases.first().unwrap().reference_strand == Strand::Forward {
                            reference_bases
                                .iter()
                                .map(|x| x.reference_nucleotide.as_str())
                                .collect()
                        } else {
                            reference_bases
                                .iter()
                                .rev()
                                .map(|x| x.reference_nucleotide.as_str())
                                .collect()
                        };

                        let graph_operation: GraphOperation = GraphOperation::new(
                            reference_chromosome_id,
                            reference_position_1,
                            reference_strand.clone(),
                            GraphOperationType::Skip,
                            reference_chromosome_id,
                            reference_position_2,
                            reference_strand.clone(),
                            GraphOperationType::Skip,
                            sequence.into(),
                            VariantType::ExonTruncation
                        );

                        variant_records.push(
                            VariantRecord::new(
                                read_id,
                                record.get_start(),
                                record.get_end(),
                                graph_operation
                            )
                        );
                    }
                }
            },
            AlignmentModelKind::Event(AlignmentModelEventKind::Boundary) => {
                for reference_bases in record.get_skipped().as_ref().unwrap().iter() {
                    let reference_chromosome_id: ReferenceChromosomeID = reference_bases.first().unwrap().reference_chromosome_id;
                    let reference_position_1: ReferencePosition = reference_bases.first().unwrap().reference_position;
                    let reference_position_2: ReferencePosition = reference_bases.last().unwrap().reference_position;
                    let reference_strand: &Strand = &reference_bases.first().unwrap().reference_strand;
                    let sequence: String = if reference_bases.first().unwrap().reference_strand == Strand::Forward {
                        reference_bases
                            .iter()
                            .map(|x| x.reference_nucleotide.as_str())
                            .collect()
                    } else {
                        reference_bases
                            .iter()
                            .rev()
                            .map(|x| x.reference_nucleotide.as_str())
                            .collect()
                    };

                    let graph_operation: GraphOperation = GraphOperation::new(
                        reference_chromosome_id,
                        reference_position_1,
                        reference_strand.clone(),
                        GraphOperationType::Skip,
                        reference_chromosome_id,
                        reference_position_2,
                        reference_strand.clone(),
                        GraphOperationType::Skip,
                        sequence.into(),
                        VariantType::ExonTruncation
                    );

                    variant_records.push(
                        VariantRecord::new(
                            read_id,
                            record.get_start(),
                            record.get_end(),
                            graph_operation
                        )
                    );
                }
            },
            AlignmentModelKind::Event(AlignmentModelEventKind::Deletion) => {
                let graph_operation: GraphOperation = GraphOperation::new(
                    record.get_chromosome_1(),
                    record.get_position_1(),
                    record.get_strand_1().clone(),
                    record.get_operation_1().clone(),
                    record.get_chromosome_2(),
                    record.get_position_2(),
                    record.get_strand_2().clone(),
                    record.get_operation_2().clone(),
                    "".into(),
                    VariantType::Deletion
                );

                variant_records.push(
                    VariantRecord::new(
                        read_id,
                        record.get_start(),
                        record.get_end(),
                        graph_operation
                    )
                );
            },
            AlignmentModelKind::Event(AlignmentModelEventKind::Splicing) => {
                if record.get_context().is_some() {
                    // The read bases between the event's two flanks: empty for a clean
                    // junction, and the untemplated insert when the aligner spelled one as an
                    // `I` op against the junction. Unlike the splicing record below, a fusion
                    // breakend is defined partly by those bases, so they belong on it, kept
                    // whole or not at all.
                    if record.get_context().as_ref().unwrap() == &AlignmentModelContext::Event(AlignmentModelEventContext::FusionGene)
                        && passes_base_quality(record.get_base_quality_scores(), min_base_quality) {
                        let graph_operation: GraphOperation = GraphOperation::new(
                            record.get_chromosome_1(),
                            record.get_position_1(),
                            record.get_strand_1().clone(),
                            record.get_operation_1().clone(),
                            record.get_chromosome_2(),
                            record.get_position_2(),
                            record.get_strand_2().clone(),
                            record.get_operation_2().clone(),
                            record.get_sequence().clone(),
                            VariantType::FusionGene
                        );

                        variant_records.push(
                            VariantRecord::new(
                                read_id,
                                record.get_start(),
                                record.get_end(),
                                graph_operation
                            )
                        );
                    }

                    // Identify non-canonical junction.
                    if record.get_context().as_ref().unwrap() == &AlignmentModelContext::Event(AlignmentModelEventContext::NonCanonicalSplicing) {
                        let graph_operation: GraphOperation = GraphOperation::new(
                            record.get_chromosome_1(),
                            record.get_position_1(),
                            record.get_strand_1().clone(),
                            record.get_operation_1().clone(),
                            record.get_chromosome_2(),
                            record.get_position_2(),
                            record.get_strand_2().clone(),
                            record.get_operation_2().clone(),
                            "".into(),
                            VariantType::NonCanonicalSplicing
                        );

                        variant_records.push(
                            VariantRecord::new(
                                read_id,
                                record.get_start(),
                                record.get_end(),
                                graph_operation
                            )
                        );
                    }
                }

                if record.get_skipped().is_some() {
                    for reference_bases in record.get_skipped().as_ref().unwrap().iter() {
                        let reference_chromosome_id: ReferenceChromosomeID = reference_bases.first().unwrap().reference_chromosome_id;
                        let reference_position_1: ReferencePosition = reference_bases.first().unwrap().reference_position;
                        let reference_position_2: ReferencePosition = reference_bases.last().unwrap().reference_position;
                        let reference_strand: &Strand = &reference_bases.first().unwrap().reference_strand;
                        let sequence: String = if reference_bases.first().unwrap().reference_strand == Strand::Forward {
                            reference_bases
                                .iter()
                                .map(|x| x.reference_nucleotide.as_str())
                                .collect()
                        } else {
                            reference_bases
                                .iter()
                                .rev()
                                .map(|x| x.reference_nucleotide.as_str())
                                .collect()
                        };

                        let graph_operation: GraphOperation = GraphOperation::new(
                            reference_chromosome_id,
                            reference_position_1,
                            reference_strand.clone(),
                            GraphOperationType::Skip,
                            reference_chromosome_id,
                            reference_position_2,
                            reference_strand.clone(),
                            GraphOperationType::Skip,
                            sequence.into(),
                            VariantType::ExonTruncation
                        );

                        variant_records.push(
                            VariantRecord::new(
                                read_id,
                                record.get_start(),
                                record.get_end(),
                                graph_operation
                            )
                        );
                    }
                }
            },
            _ => {
                continue;
            }
        };
    }

    variant_records
}


pub(crate) fn identify_terminal_soft_clip_variant_records(
    model: &AlignmentModel,
    min_mapping_quality: MappingQuality,
    min_base_quality: BaseQuality,
    min_length: u32,
    graph_operation: impl Fn(ReferenceChromosomeID, ReferencePosition, Strand, GraphOperationType, Box<str>) -> GraphOperation
) -> Vec<VariantRecord> {
    let mut variant_records: Vec<VariantRecord> = Vec::new();
    for range in model.terminal_soft_clip_ranges() {
        if range.end - range.start < min_length.max(1) {
            continue;
        }
        let leading: bool = range.start == 0;
        let anchor: &AlignmentModelBase = model.get_base(if leading { range.end } else { range.start - 1 });
        let AlignmentModelBasePlacement::Placed { chromosome_id, position, strand, mapping_quality } = anchor.get_placement() else {
            continue;
        };
        if *mapping_quality < min_mapping_quality {
            continue;
        }
        let bases: &[AlignmentModelBase] = &model.get_bases()[range.start as usize..range.end as usize];
        let base_quality_scores: Vec<BaseQuality> = bases.iter().map(AlignmentModelBase::get_base_quality).collect();
        if !passes_base_quality(&base_quality_scores, min_base_quality) {
            continue;
        }
        let sequence: String = bases.iter().map(|base| base.get_nucleotide().as_str()).collect();
        let operation: GraphOperationType = GraphOperationType::for_breakpoint(strand, leading);
        variant_records.push(VariantRecord::new(
            model.get_read_id(),
            range.start,
            range.end - 1,
            graph_operation(*chromosome_id, *position, strand.clone(), operation, sequence.into())
        ));
    }
    variant_records
}


pub(crate) fn sort_variant_records(variant_records: &mut [VariantRecord]) {
    variant_records.sort_by(|a, b| {
        a.get_chromosome_1()
            .cmp(&b.get_chromosome_1())
            .then(a.get_position_1().cmp(&b.get_position_1()))
    });
}


pub(crate) fn passes_base_quality(base_quality_scores: &[BaseQuality], min_base_quality: BaseQuality) -> bool {
    2 * base_quality_scores.iter().filter(|&&quality| quality >= min_base_quality).count() >= base_quality_scores.len()
}
