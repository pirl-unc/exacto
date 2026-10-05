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
use std::collections::HashSet;

use crate::prelude::*;
use crate::calling::alignment_model_record_calling::*;


pub struct RNAVariantRecordCaller<'a> {
    min_mapping_quality: MappingQuality,
    min_base_quality: BaseQuality,
    min_terminal_soft_clip_ins_len: u32,
    breakpoint_rescue: Option<RNABreakpointRescue<'a>>
}

impl<'a> RNAVariantRecordCaller<'a> {
    pub fn new(
        min_mapping_quality: MappingQuality,
        min_base_quality: BaseQuality,
        min_terminal_soft_clip_ins_len: u32,
        breakpoint_rescue: Option<RNABreakpointRescue<'a>>
    ) -> Self {
        RNAVariantRecordCaller { 
            min_mapping_quality,
            min_base_quality,
            min_terminal_soft_clip_ins_len,
            breakpoint_rescue
        }
    }
}

impl VariantRecordCaller for RNAVariantRecordCaller<'_> {
    type Input = TranscriptModel;

    fn call(&self, transcript_model: &TranscriptModel) -> Vec<VariantRecord> {
        let model: &AlignmentModel = transcript_model.get_alignment_model();
        let annotation: &TranscriptModelAnnotation = transcript_model.get_annotation();
        let read_id: ReadID = model.get_read_id();
        let records: Vec<AlignmentModelRecord> = identify_alignment_model_records(model, Some(annotation));

        // Step 1. Mismatches and insertions.
        let num_bases: u32 = model.num_bases();
        let mut variant_records: Vec<VariantRecord> = identify_base_variant_records(
            read_id,
            records.iter().filter(|record| {
                *record.get_kind() != AlignmentModelKind::Base(AlignmentModelBaseKind::Mismatch)
                    || (record.get_start() != 0 && record.get_end() + 1 != num_bases)
            }),
            self.min_mapping_quality,
            self.min_base_quality
        );

        // Step 2. Cryptic exons, intron retentions and UTR extensions.
        let context_variant_records: Vec<VariantRecord> = identify_context_variant_records(
            read_id,
            annotation,
            &records,
            &variant_records
        );
        variant_records.extend(context_variant_records);

        // Step 3. Terminal soft clips, each an insertion beside its anchor.
        // A clip of at most three bases holding CC or GG is a template-switching artifact and is left out.
        let mut soft_clip_variant_records: Vec<VariantRecord> = identify_terminal_soft_clip_variant_records(
            model,
            self.min_mapping_quality,
            self.min_base_quality,
            self.min_terminal_soft_clip_ins_len,
            |chromosome_id, position, strand, operation, sequence| {
                let (position_1, position_2): (ReferencePosition, ReferencePosition) = match operation {
                    GraphOperationType::Upstream => (position.saturating_sub(1), position),
                    _ => (position, position + 1)
                };
                GraphOperation::new(
                    chromosome_id,
                    position_1,
                    strand.clone(),
                    GraphOperationType::Downstream,
                    chromosome_id,
                    position_2,
                    strand,
                    GraphOperationType::Upstream,
                    sequence,
                    VariantType::Insertion
                )
            }
        );
        soft_clip_variant_records.retain(|vr| {
            let sequence: String = vr.get_sequence().to_uppercase();
            sequence.len() > 3 || !(sequence.contains("CC") || sequence.contains("GG"))
        });
        variant_records.extend(soft_clip_variant_records);

        // Step 4. Deletions, splicing events, fusions and circular RNAs.
        variant_records.extend(identify_event_variant_records(
            read_id,
            &records,
            self.min_mapping_quality,
            self.min_base_quality
        ));

        sort_variant_records(&mut variant_records);

        // Step 5. Retype large insertions as breakpoints.
        match &self.breakpoint_rescue {
            Some(breakpoint_rescue) => breakpoint_rescue.retype_insertions(transcript_model, variant_records),
            None => variant_records
        }
    }
}


/// Runs of records off the exons (intronic or intergenic bases, and the deletions between
/// them) as variant records. A run that is already a mismatch or insertion record is left out.
///
/// This function identifies the following variant types:
/// - Cryptic exon
/// - Intron retention
/// - UTR extension
fn identify_context_variant_records(
    read_id: ReadID,
    annotation: &TranscriptModelAnnotation,
    records: &[AlignmentModelRecord],
    base_variant_records: &[VariantRecord]
) -> Vec<VariantRecord> {
    let mut variant_records: Vec<VariantRecord> = Vec::new();

    // Step 1. Get a set of all included read bases.
    let mut included_read_bases: HashSet<(ReadPosition, ReadPosition)> = HashSet::new();
    for variant_record in base_variant_records {
        included_read_bases.insert(
            (variant_record.get_read_position_1(), variant_record.get_read_position_2())
        );
    }

    // Step 2. Join the records off the exons.
    let mut uf: UnionFind = UnionFind::new();
    for i in 0..records.len() {
        let record: &AlignmentModelRecord = records.get(i).unwrap();
        if *record.get_record_type() == AlignmentModelRecordType::Base &&
            *record.get_context().as_ref().unwrap() != AlignmentModelContext::Base(AlignmentModelBaseContext::Exonic) {
            uf.union(i as u32, i as u32);
        }
    }
    for i in 1..records.len() {
        let prev_record: &AlignmentModelRecord = records.get(i - 1).unwrap();
        let curr_record: &AlignmentModelRecord = records.get(i).unwrap();
        match (prev_record.get_record_type(), curr_record.get_record_type()) {
            (AlignmentModelRecordType::Event, AlignmentModelRecordType::Event) => {
                continue;
            },
            (AlignmentModelRecordType::Event, AlignmentModelRecordType::Base) => {
                if *curr_record.get_context().as_ref().unwrap() != AlignmentModelContext::Base(AlignmentModelBaseContext::Exonic)
                    && *prev_record.get_kind() == AlignmentModelKind::Event(AlignmentModelEventKind::Deletion)
                    && prev_record.get_chromosome_1() == curr_record.get_chromosome_1()
                    && prev_record.get_chromosome_2() == curr_record.get_chromosome_2()
                    && prev_record.get_strand_1() == curr_record.get_strand_1()
                    && prev_record.get_strand_2() == curr_record.get_strand_2()
                    && (prev_record.get_position_1().abs_diff(curr_record.get_position_2()) <= 1
                        || prev_record.get_position_2().abs_diff(curr_record.get_position_1()) <= 1) {
                    // Contexts of the two bases flanking the event record.
                    let base_context_1: Option<&AlignmentModelBaseContext> = annotation.get_base_context(prev_record.get_start());
                    let base_context_2: Option<&AlignmentModelBaseContext> = annotation.get_base_context(prev_record.get_end());
                    if let (Some(base_context_1), Some(base_context_2)) = (base_context_1, base_context_2) {
                        if *base_context_1 != AlignmentModelBaseContext::Exonic
                            && *base_context_1 == *base_context_2
                            && *base_context_1 == *curr_record.get_context().as_ref().unwrap().as_base().unwrap() {
                            uf.union(i as u32 - 1, i as u32);
                        }
                    }
                }
            },
            (AlignmentModelRecordType::Base, AlignmentModelRecordType::Event) => {
                if *prev_record.get_context().as_ref().unwrap() != AlignmentModelContext::Base(AlignmentModelBaseContext::Exonic)
                    && *curr_record.get_kind() == AlignmentModelKind::Event(AlignmentModelEventKind::Deletion)
                    && prev_record.get_chromosome_1() == curr_record.get_chromosome_1()
                    && prev_record.get_chromosome_2() == curr_record.get_chromosome_2()
                    && prev_record.get_strand_1() == curr_record.get_strand_1()
                    && prev_record.get_strand_2() == curr_record.get_strand_2()
                    && (prev_record.get_position_1().abs_diff(curr_record.get_position_2()) <= 1
                        || prev_record.get_position_2().abs_diff(curr_record.get_position_1()) <= 1) {
                    // Same substitution as the branch above, on the other side:
                    // here the *current* record is the event, so its flanking
                    // bases are the ones whose context is compared.
                    let base_context_1: Option<&AlignmentModelBaseContext> = annotation.get_base_context(curr_record.get_start());
                    let base_context_2: Option<&AlignmentModelBaseContext> = annotation.get_base_context(curr_record.get_end());
                    if let (Some(base_context_1), Some(base_context_2)) = (base_context_1, base_context_2) {
                        if *base_context_1 != AlignmentModelBaseContext::Exonic
                            && *base_context_1 == *base_context_2
                            && *base_context_1 == *prev_record.get_context().as_ref().unwrap().as_base().unwrap() {
                            uf.union(i as u32 - 1, i as u32);
                        }
                    }
                }
            }
            (AlignmentModelRecordType::Base, AlignmentModelRecordType::Base) => {
                if *prev_record.get_context().as_ref().unwrap() != AlignmentModelContext::Base(AlignmentModelBaseContext::Exonic)
                    && *curr_record.get_context().as_ref().unwrap() != AlignmentModelContext::Base(AlignmentModelBaseContext::Exonic)
                    && *prev_record.get_context().as_ref().unwrap().as_base().unwrap() == *curr_record.get_context().as_ref().unwrap().as_base().unwrap()
                    && prev_record.get_chromosome_1() == curr_record.get_chromosome_1()
                    && prev_record.get_chromosome_2() == curr_record.get_chromosome_2()
                    && prev_record.get_strand_1() == curr_record.get_strand_1()
                    && prev_record.get_strand_2() == curr_record.get_strand_2()
                    && (prev_record.get_position_1().abs_diff(curr_record.get_position_2()) <= 1
                        || prev_record.get_position_2().abs_diff(curr_record.get_position_1()) <= 1) {
                    uf.union(i as u32 - 1, i as u32);
                }
            }
        }
    }
    for cluster in uf.get_clusters().iter() {
        // Sort the record positions.
        let mut record_indices: Vec<u32> = cluster.iter().map(|&pos| pos).collect();
        record_indices.sort();

        let first_record: &AlignmentModelRecord = records.get(*record_indices.first().unwrap() as usize).unwrap();
        let last_record: &AlignmentModelRecord = records.get(*record_indices.last().unwrap() as usize).unwrap();

        let prev_record: Option<&AlignmentModelRecord> = if *record_indices.first().unwrap() > 0u32 {
            records.get(record_indices[0] as usize - 1)
        } else {
            None
        };

        let next_record: Option<&AlignmentModelRecord> = if *record_indices.last().unwrap() < records.len() as u32 - 1 {
            records.get(*record_indices.last().unwrap() as usize + 1)
        } else {
            None
        };

        let mut variant_type: VariantType = VariantType::CrypticExon;
        if prev_record.is_some() {
            let prev_record_: &AlignmentModelRecord = prev_record.unwrap();
            if *prev_record_.get_record_type() == AlignmentModelRecordType::Base
                && *prev_record_.get_context().as_ref().unwrap() == AlignmentModelContext::Base(AlignmentModelBaseContext::Exonic)
                && prev_record_.get_chromosome_1() == first_record.get_chromosome_1()
                && prev_record_.get_chromosome_2() == first_record.get_chromosome_2()
                && prev_record_.get_strand_1() == first_record.get_strand_1()
                && prev_record_.get_strand_2() == first_record.get_strand_2()
                && (prev_record_.get_position_1().abs_diff(first_record.get_position_2()) <= 1
                    || prev_record_.get_position_2().abs_diff(first_record.get_position_1()) <= 1) {
                if *first_record.get_context().as_ref().unwrap() == AlignmentModelContext::Base(AlignmentModelBaseContext::Intronic) {
                    variant_type = VariantType::IntronRetention;
                }
                if *first_record.get_context().as_ref().unwrap() == AlignmentModelContext::Base(AlignmentModelBaseContext::Intergenic) {
                    variant_type = VariantType::UTRExtension;
                }
            }
        }
        if next_record.is_some() {
            let next_record_: &AlignmentModelRecord = next_record.unwrap();
            if *next_record_.get_record_type() == AlignmentModelRecordType::Base
                && *next_record_.get_context().as_ref().unwrap() == AlignmentModelContext::Base(AlignmentModelBaseContext::Exonic)
                && next_record_.get_chromosome_1() == last_record.get_chromosome_1()
                && next_record_.get_chromosome_2() == last_record.get_chromosome_2()
                && next_record_.get_strand_1() == last_record.get_strand_1()
                && next_record_.get_strand_2() == last_record.get_strand_2()
                && (next_record_.get_position_1().abs_diff(last_record.get_position_2()) <= 1
                    || next_record_.get_position_2().abs_diff(last_record.get_position_1()) <= 1) {
                if *last_record.get_context().as_ref().unwrap() == AlignmentModelContext::Base(AlignmentModelBaseContext::Intronic) {
                    variant_type = VariantType::IntronRetention;
                }
                if *last_record.get_context().as_ref().unwrap() == AlignmentModelContext::Base(AlignmentModelBaseContext::Intergenic) {
                    variant_type = VariantType::UTRExtension;
                }
            }
        }

        assert_eq!(
            first_record.get_chromosome_1(),
            last_record.get_chromosome_1(),
            "The first ({}) and last ({}) records must have the same chromosome 1.",
            first_record.get_chromosome_1(),
            last_record.get_chromosome_1()
        );
        assert_eq!(
            first_record.get_chromosome_2(),
            last_record.get_chromosome_2(),
            "The first ({}) and last ({}) records must have the same chromosome 2.",
            first_record.get_chromosome_2(),
            last_record.get_chromosome_2()
        );
        assert_eq!(
            first_record.get_strand_1(),
            last_record.get_strand_1(),
            "The first ({}) and last ({}) records must have the same strand 1.",
            first_record.get_strand_1().as_str(),
            last_record.get_strand_1().as_str()
        );
        assert_eq!(
            first_record.get_strand_2(),
            last_record.get_strand_2(),
            "The first ({}) and last ({}) records must have the same strand 2.",
            first_record.get_strand_2().as_str(),
            last_record.get_strand_2().as_str()
        );

        let graph_operation: GraphOperation = if *first_record.get_strand_1() == Strand::Forward {
            GraphOperation::new(
                first_record.get_chromosome_1(),
                first_record.get_position_1(),
                first_record.get_strand_1().clone(),
                GraphOperationType::Include,
                last_record.get_chromosome_1(),
                last_record.get_position_2(),
                last_record.get_strand_1().clone(),
                GraphOperationType::Include,
                "".into(),
                variant_type
            )
        } else {
            GraphOperation::new(
                first_record.get_chromosome_1(),
                first_record.get_position_2(),
                first_record.get_strand_1().clone(),
                GraphOperationType::Include,
                last_record.get_chromosome_1(),
                last_record.get_position_1(),
                last_record.get_strand_1().clone(),
                GraphOperationType::Include,
                "".into(),
                variant_type
            )
        };

        if included_read_bases.contains(&(first_record.get_start(), last_record.get_end())) == false {
            variant_records.push(
                VariantRecord::new(
                    read_id,
                    first_record.get_start(),
                    last_record.get_end(),
                    graph_operation
                )
            );
        }
    }

    variant_records
}
