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
use exacto_core::prelude::{reverse_complement, Strand, LIST_SEPARATOR};
use exacto_integrator::prelude::*;
use std::collections::{HashMap, HashSet};
use std::str::FromStr;

use crate::prelude::*;


pub fn build_transcript_set(
    assembled_transcript_support_records: &Vec<AssembledTranscriptSupportRecord>,
    assembled_transcript_model_alignment_records: &Vec<AssembledTranscriptModelAlignmentRecord>,
    assembled_transcript_variant_records: &Vec<AssembledTranscriptVariantRecord>,
    dna_variant_records: &Vec<DNAVariantRecord>,
    integrated_variant_records: &Vec<IntegratedVariantRecord>,
    reference_stitched_spans: &Vec<ReferenceStitchedSpans>
) -> Result<AssembledTranscriptSet, TranslatorError> {
    // Step 1. Index records
    let assembled_transcript_support_index: HashMap<&str, &AssembledTranscriptSupportRecord> =
        assembled_transcript_support_records
            .iter()
            .map(|r| (r.assembled_transcript_name.as_ref(), r))
            .collect();
    let reference_stitched_spans_index: HashMap<&str, &ReferenceStitchedSpans> =
        reference_stitched_spans
            .iter()
            .map(|r| (r.assembled_transcript_name.as_ref(), r))
            .collect();
    let mut alignment_records_index: HashMap<Box<str>, Vec<&AssembledTranscriptModelAlignmentRecord>> = HashMap::new();
    for r in assembled_transcript_model_alignment_records {
        alignment_records_index
            .entry(r.assembled_transcript_name.clone())
            .or_default()
            .push(r);
    }
    let mut rna_variant_records_index: HashMap<Box<str>, Vec<&AssembledTranscriptVariantRecord>> = HashMap::new();
    for r in assembled_transcript_variant_records {
        rna_variant_records_index
            .entry(r.assembled_transcript_name.clone())
            .or_default()
            .push(r);
    }
    let mut integrated_variant_records_index: HashMap<Box<str>, Vec<&IntegratedVariantRecord>> = HashMap::new();
    for r in integrated_variant_records {
        integrated_variant_records_index
            .entry(r.assembled_transcript_name.clone())
            .or_default()
            .push(r);
    }
    let dna_variant_records_by_id: HashMap<u32, &DNAVariantRecord> = dna_variant_records.iter().map(|r| (r.variant_id, r)).collect();

    // Step 2. Build one Transcript per assembled_transcript_name. A read yields
    // exactly one transcript model, so the assembled transcript name is the
    // whole join key.
    let mut keys: Vec<Box<str>> = alignment_records_index.keys().cloned().collect();
    keys.sort_unstable();

    // Every aligned transcript must have a support row. Falling back to concatenated
    // alignment-row sequences (which is what this did) hides the failure that actually
    // happens in practice: the support file and the alignments file disagree about naming —
    // bare cluster ids from `determine-rna-consensus` against `cid_*`/`gid_*` from an external
    // assembler — and then *every* transcript silently gets a stitched-together sequence and
    // an empty read list, which translates into plausible, wrong proteoforms. Surplus support
    // rows are fine and deliberately not checked: a cluster whose consensus never aligned has
    // no alignment rows to answer for.
    let unmatched: Vec<&str> = keys
        .iter()
        .map(|key| key.as_ref())
        .filter(|key| !assembled_transcript_support_index.contains_key(key))
        .collect();
    if !unmatched.is_empty() {
        const NUM_NAMES_SHOWN: usize = 10;
        let shown: Vec<&str> = unmatched.iter().take(NUM_NAMES_SHOWN).copied().collect();
        return Err(TranslatorError::UnmatchedTranscripts {
            num_unmatched: unmatched.len(),
            num_transcripts: keys.len(),
            names: format!(
                "{}{}",
                shown.join(", "),
                if unmatched.len() > NUM_NAMES_SHOWN { ", ..." } else { "" }
            ).into_boxed_str()
        });
    }
    let mut transcripts: Vec<AssembledTranscript> = Vec::with_capacity(keys.len());
    for (transcript_index, assembled_transcript_name) in keys.into_iter().enumerate() {
        let alignment_rows: &Vec<&AssembledTranscriptModelAlignmentRecord> = &alignment_records_index[&assembled_transcript_name];
        // Guaranteed present: the unmatched check above panics otherwise.
        let support: &AssembledTranscriptSupportRecord =
            assembled_transcript_support_index[assembled_transcript_name.as_ref()];
        let (sequence, read_ids): (Box<str>, Vec<Box<str>>) =
            (support.sequence.clone(), split_read_names(&support.read_names));
        // Read coordinates and codons index the sequence by byte.
        if !sequence.is_ascii() {
            return Err(TranslatorError::NonAsciiSequence { transcript: assembled_transcript_name.clone() });
        }

        // TranscriptAlignment
        let (reference_gene_names, reference_transcript_ids): (Vec<Box<str>>, Vec<Box<str>>) =
            alignment_rows
                .first()
                .map(|r| {
                    let gene_names: Vec<Box<str>> = r.reference_gene_name
                        .split(LIST_SEPARATOR).filter(|s| !s.is_empty()).map(|s| s.into()).collect();
                    let transcript_ids: Vec<Box<str>> = r.reference_transcript_id
                        .split(LIST_SEPARATOR).filter(|s| !s.is_empty()).map(|s| s.into()).collect();
                    (gene_names, transcript_ids)
                })
                .unwrap_or_default();

        let mut transcript_alignment: AssembledTranscriptAlignment = AssembledTranscriptAlignment::new(
            transcript_index,
            reference_gene_names,
            reference_transcript_ids
        );

        // Every Base row carries the bases it aligned, so a sequence from another run (the same
        // cluster ids, different consensus) shows up as a row that does not spell it.
        for row in alignment_rows.iter() {
            let item: AssembledTranscriptAlignmentRecord = build_transcript_alignment_item(row)?;
            if matches!(item.item_type, TranscriptAlignmentItemType::Base { .. })
                && !spells_sequence(&sequence, row, item.graph_operation_view.get_strand_1()) {
                return Err(TranslatorError::SequenceMismatch {
                    transcript: assembled_transcript_name.clone(),
                    index: row.index,
                    read_start: row.read_start,
                    read_end: row.read_end
                });
            }
            transcript_alignment.add_item(item);
        }

        // BiMap<rna_variant_id, GraphOperationView>
        //
        // RNA variant records are denormalized by reference_transcript_id —
        // the same biological variant (same GOV) may appear multiple times
        // with distinct variant_ids, one per matching reference transcript.
        // BiMap enforces unique values, so a naive insert loop would evict
        // every prior (id, GOV) pair and leave only the last id, breaking
        // the integration lookup whenever it references a non-last id.
        //
        // Resolve this by dedup-on-GOV with "first id wins" as canonical,
        // and build `rna_variant_id_remap` so the integration loop below
        // can translate any non-canonical id back to its canonical id.
        let mut rna_variants_bimap: BiMap<u32, GraphOperationView> = BiMap::new();
        let mut rna_variant_id_remap: HashMap<u32, u32> = HashMap::new();
        // (read_start, read_end, canonical rna_variant_id) for every canonical
        // RNA variant. The per-nucleotide lookup matches a alignment base to its
        // variant by read-coordinate containment (see
        // `Transcript::try_get_rna_variant_id_for_read_span`), which is the only
        // reliable key for cryptic-exon / UTR calls — their GOVs do not match
        // the matched-sequence alignment bases.
        let mut rna_variant_read_spans: Vec<(u32, u32, u32)> = Vec::new();
        // The variant type decides whether a variant keeps the reading frame
        // (`AssembledTranscript::keeps_reading_frame`): a junction and a deletion
        // have the same descriptor shape.
        let mut rna_variant_types: HashMap<u32, VariantType> = HashMap::new();
        if let Some(rna_rows) = rna_variant_records_index.get(&assembled_transcript_name) {
            for r in rna_rows {
                let gov: GraphOperationView = graph_operation_view_from_rna_variant_record(r)?;
                let canonical_id: u32 = match rna_variants_bimap.get_by_right(&gov).copied() {
                    Some(existing) => existing,
                    None => {
                        rna_variants_bimap.insert(r.variant_id, gov);
                        rna_variant_read_spans.push((r.read_start, r.read_end, r.variant_id));
                        rna_variant_types.insert(
                            r.variant_id,
                            parse(&assembled_transcript_name, "variant_type", &r.variant_type)?
                        );
                        r.variant_id
                    }
                };
                rna_variant_id_remap.insert(r.variant_id, canonical_id);
            }
        }

        // Collect DNA ids. Key `integrated_variant_ids` on the canonical
        // rna_variant_id so the per-nucleotide lookup (which receives the
        // canonical id from `rna_variants_bimap`) finds its DNA partners.
        let mut integrated_variant_ids: HashMap<u32, HashSet<u32>> = HashMap::new(); // HashMap<canonical_rna_variant_id, HashSet<dna_variant_id>>
        let mut dna_variant_ids_needed: HashSet<u32> = HashSet::new();
        if let Some(records) = integrated_variant_records_index.get(&assembled_transcript_name) {
            for record in records {
                dna_variant_ids_needed.insert(record.dna_variant_id);
                let canonical_rna_id: u32 = rna_variant_id_remap
                    .get(&record.rna_variant_id)
                    .copied()
                    .unwrap_or(record.rna_variant_id);
                integrated_variant_ids
                    .entry(canonical_rna_id)
                    .or_default()
                    .insert(record.dna_variant_id);
            }
        }

        // DNA variants BiMap + read_names map - only the ids this transcript
        // integrates against. read_names are denormalized so downstream
        // proteoform record builders don't need the global DNA stream.
        //
        // Apply the same dedup-on-GOV / canonical-id pattern as for RNA, in
        // case two DNA variant records share a GOV (rare in practice but
        // would silently drop integration linkages if not handled).
        let mut dna_variants_bimap: BiMap<u32, GraphOperationView> = BiMap::new();
        let mut dna_variant_read_names: HashMap<u32, Box<str>> = HashMap::new();
        let mut dna_variant_id_remap: HashMap<u32, u32> = HashMap::new();
        //
        // Sorted so that which of two equal descriptors becomes canonical does not
        // depend on hash order.
        let mut dna_variant_ids_needed: Vec<u32> = dna_variant_ids_needed.into_iter().collect();
        dna_variant_ids_needed.sort_unstable();
        for dna_id in dna_variant_ids_needed {
            let dna_record: &DNAVariantRecord = dna_variant_records_by_id.get(&dna_id)
                .ok_or_else(|| TranslatorError::UnknownDnaVariant {
                    transcript: assembled_transcript_name.clone(),
                    dna_variant_id: dna_id
                })?;
            let gov: GraphOperationView = graph_operation_view_from_dna_variant_record(dna_record, &assembled_transcript_name)?;
            let canonical_id: u32 = match dna_variants_bimap.get_by_right(&gov).copied() {
                Some(existing) => existing,
                None => {
                    dna_variants_bimap.insert(dna_id, gov);
                    dna_variant_read_names.insert(dna_id, dna_record.read_names.clone());
                    dna_id
                }
            };
            dna_variant_id_remap.insert(dna_id, canonical_id);
        }

        // Rewrite `integrated_variant_ids` so each DNA id is canonical too.
        if !dna_variant_id_remap.is_empty() {
            for dna_ids in integrated_variant_ids.values_mut() {
                let canonicalized: HashSet<u32> = dna_ids.iter()
                    .map(|id| dna_variant_id_remap.get(id).copied().unwrap_or(*id))
                    .collect();
                *dna_ids = canonicalized;
            }
        }

        // Reference-stitched provenance spans for this transcript, if supplied. The
        // stitched_length cross-check makes a stale annotations file (from a different
        // stitch run) abort instead of silently mislabeling provenance.
        let reference_stitched_intervals: Vec<(u32, u32)> =
            match reference_stitched_spans_index.get(&*assembled_transcript_name) {
                Some(spans) if spans.stitched_length as usize != sequence.len() => {
                    return Err(TranslatorError::StitchedLengthMismatch {
                        transcript: assembled_transcript_name.clone(),
                        sequence_length: sequence.len(),
                        stitched_length: spans.stitched_length
                    });
                }
                Some(spans) => spans.intervals.clone(),
                None => Vec::new()
            };

        // The assembled transcript name is the join key, so it doubles as the
        // Transcript.id.
        let transcript_id: Box<str> = assembled_transcript_name.clone();

        transcripts.push(AssembledTranscript::new(
            transcript_id,
            sequence,
            read_ids,
            assembled_transcript_name,
            transcript_alignment,
            rna_variants_bimap,
            dna_variants_bimap,
            integrated_variant_ids,
            dna_variant_read_names,
            rna_variant_read_spans,
            rna_variant_types,
            reference_stitched_intervals
        ));
    }

    Ok(AssembledTranscriptSet::new(transcripts))
}


/// Whether a Base row spells `sequence` at its read span `[read_start, read_end]`. The row
/// carries its bases in forward (reference) orientation, so a minus-strand row is
/// reverse-complemented first. U reads as T and case is ignored.
fn spells_sequence(sequence: &str, row: &AssembledTranscriptModelAlignmentRecord, strand: &Strand) -> bool {
    let span: &str = match sequence.get(row.read_start as usize..=row.read_end as usize) {
        Some(span) => span,
        None => return false
    };
    let is_nucleotide = |base: char| "ACGTUNRYKMBVDHSW".contains(base.to_ascii_uppercase());
    let bases: Box<str> = match strand {
        Strand::Reverse if row.sequence.chars().all(is_nucleotide) => reverse_complement(&row.sequence),
        Strand::Reverse => return false,
        _ => row.sequence.clone()
    };
    let normalize = |base: u8| match base.to_ascii_uppercase() { b'U' => b'T', other => other };
    span.len() == bases.len() && span.bytes().zip(bases.bytes()).all(|(a, b)| normalize(a) == normalize(b))
}


/// Parses one field of an input row, naming the transcript, the field and the value when it fails.
fn parse<T: FromStr>(transcript: &str, field: &'static str, value: &str) -> Result<T, TranslatorError> {
    T::from_str(value).map_err(|_| TranslatorError::UnknownValue {
        transcript: transcript.into(),
        field,
        value: value.into()
    })
}


pub(crate) fn build_transcript_alignment_item(
    record: &AssembledTranscriptModelAlignmentRecord
) -> Result<AssembledTranscriptAlignmentRecord, TranslatorError> {
    let transcript: &str = &record.assembled_transcript_name;
    let num_cycles: Option<u16> = if is_cycle_creating_descriptor(
        &record.chromosome_1,
        record.position_1,
        &record.operation_1,
        &record.chromosome_2,
        record.position_2,
        &record.operation_2
    ) {
        Some(1)
    } else {
        None
    };

    let gov: GraphOperationView = graph_operation_view_from_descriptor(
        transcript,
        &record.chromosome_1,
        record.position_1,
        &record.operation_1,
        &record.strand_1,
        &record.chromosome_2,
        record.position_2,
        &record.operation_2,
        &record.strand_2,
        &record.sequence,
        num_cycles
    )?;

    let record_type: AlignmentModelRecordType = parse(transcript, "type", &record.record_type)?;
    let item_type = match record_type {
        AlignmentModelRecordType::Base => {
            let kind = parse(transcript, "kind", &record.kind)?;
            let context = parse(transcript, "context", &record.context)?;
            TranscriptAlignmentItemType::Base { kind, context }
        }
        AlignmentModelRecordType::Event => {
            let kind = parse(transcript, "kind", &record.kind)?;
            // An unclassified event is serialized with an empty context; map it
            // back to `None`. A non-empty but unrecognized value is a genuine
            // schema error, so surface it clearly.
            let context: Option<AlignmentModelEventContext> = if record.context.trim().is_empty() {
                None
            } else {
                Some(parse(transcript, "context", &record.context)?)
            };
            TranscriptAlignmentItemType::Event { kind, context }
        }
    };

    let annotation = TranscriptAlignmentAnnotation {
        position_1_annotation: Annotation {
            gene_id: optional_id(&record.reference_gene_id_1),
            transcript_id: optional_id(&record.reference_transcript_id_1),
            exon_id: optional_id(&record.reference_exon_id_1)
        },
        position_2_annotation: Annotation {
            gene_id: optional_id(&record.reference_gene_id_2),
            transcript_id: optional_id(&record.reference_transcript_id_2),
            exon_id: optional_id(&record.reference_exon_id_2)
        },
    };

    Ok(AssembledTranscriptAlignmentRecord::new(
        record.index,
        record.read_start,
        record.read_end,
        item_type,
        gov,
        annotation
    ))
}


fn graph_operation_view_from_rna_variant_record(
    record: &AssembledTranscriptVariantRecord
) -> Result<GraphOperationView, TranslatorError> {
    let num_cycles: Option<u16> = if is_cycle_creating_descriptor(
        &record.chromosome_1,
        record.position_1,
        &record.operation_1,
        &record.chromosome_2,
        record.position_2,
        &record.operation_2
    ) {
        Some(1)
    } else {
        None
    };

    graph_operation_view_from_descriptor(
        &record.assembled_transcript_name,
        &record.chromosome_1,
        record.position_1,
        &record.operation_1,
        &record.strand_1,
        &record.chromosome_2,
        record.position_2,
        &record.operation_2,
        &record.strand_2,
        &record.sequence,
        num_cycles
    )
}


/// `transcript` is the assembled transcript the DNA variant is integrated with, named in an error.
pub(crate) fn graph_operation_view_from_dna_variant_record(
    record: &DNAVariantRecord,
    transcript: &str
) -> Result<GraphOperationView, TranslatorError> {
    // A cycle-creating operation (Upstream -> Downstream on the same chromosome
    // with position_1 < position_2, e.g. a tandem duplication) requires a cycle
    // count. DNA variant records do not carry one, so default such operations to
    // a single cycle. This keeps the cycle-count handling in the DNA load path
    // and leaves the `GraphOperationView` invariant untouched.
    let num_cycles: Option<u16> = if is_cycle_creating_descriptor(
        &record.chromosome_1,
        record.position_1,
        &record.operation_1,
        &record.chromosome_2,
        record.position_2,
        &record.operation_2
    ) {
        Some(1)
    } else {
        None
    };

    graph_operation_view_from_descriptor(
        transcript,
        &record.chromosome_1,
        record.position_1,
        &record.operation_1,
        &record.strand_1,
        &record.chromosome_2,
        record.position_2,
        &record.operation_2,
        &record.strand_2,
        &record.sequence,
        num_cycles
    )
}


fn graph_operation_view_from_descriptor(
    transcript: &str,
    chromosome_1: &str,
    position_1: u32,
    operation_1: &str,
    strand_1: &str,
    chromosome_2: &str,
    position_2: u32,
    operation_2: &str,
    strand_2: &str,
    sequence: &str,
    num_cycles: Option<u16>
) -> Result<GraphOperationView, TranslatorError> {
    let strand_1: Strand = parse(transcript, "strand_1", strand_1)?;
    let strand_2: Strand = parse(transcript, "strand_2", strand_2)?;
    let operation_type_1: GraphOperationType = parse(transcript, "operation_1", operation_1)?;
    let operation_type_2: GraphOperationType = parse(transcript, "operation_2", operation_2)?;
    Ok(GraphOperationView::new(
        chromosome_1,
        position_1,
        operation_type_1,
        strand_1,
        chromosome_2,
        position_2,
        operation_type_2,
        strand_2,
        sequence,
        num_cycles
    ))
}

/// Returns true when a variant descriptor describes a cycle-creating operation:
/// an Upstream -> Downstream pair on the same chromosome with
/// `position_1 < position_2` (e.g. a tandem duplication). This mirrors the
/// invariant enforced by `GraphOperationView::new`, which requires a cycle count
/// for such operations.
fn is_cycle_creating_descriptor(
    chromosome_1: &str,
    position_1: u32,
    operation_1: &str,
    chromosome_2: &str,
    position_2: u32,
    operation_2: &str
) -> bool {
    matches!(GraphOperationType::from_str(operation_1), Ok(GraphOperationType::Upstream))
        && matches!(GraphOperationType::from_str(operation_2), Ok(GraphOperationType::Downstream))
        && chromosome_1 == chromosome_2
        && position_1 < position_2
}


fn optional_id(s: &str) -> Option<Box<str>> {
    if s.is_empty() { None } else { Some(s.into()) }
}


pub fn build_nucleotide_records<'a>(
    transcript_set: &'a AssembledTranscriptSet
) -> impl Iterator<Item =NucleotideRecord> + 'a {
    transcript_set.iter().flat_map(move |transcript| {
        let assembled_transcript_name: Box<str> = transcript.get_assembled_transcript_name().into();

        transcript.proteoforms.iter().flat_map(move |proteoform| {
            let assembled_transcript_name = assembled_transcript_name.clone();
            transcript.get_amino_acids(proteoform).flat_map(move |amino_acid| {
                let assembled_transcript_name = assembled_transcript_name.clone();
                let amino_acid_letter: Box<str> = amino_acid.get_amino_acid().to_string().into_boxed_str();
                let is_amino_acid_variant: bool = amino_acid.is_variant();
                let amino_acid_index: u32 = amino_acid.get_index();
                amino_acid.nucleotides.into_iter().enumerate().map(move |(codon_idx, transcript_nucleotide)| {
                    let rna_variant: Option<Box<str>> = transcript_nucleotide
                        .get_rna_variant_id()
                        .and_then(|id|{
                                Some(transcript.get_rna_variant(id).as_boxed_str())
                        });
                    // Descriptors in ascending id order, the order `dna_variant_ids` is written in;
                    // a HashSet's own order changes from run to run.
                    let dna_variant: Option<Box<str>> = transcript_nucleotide
                        .get_dna_variant_ids()
                        .as_ref()
                        .and_then(|ids| {
                            let mut ids: Vec<u32> = ids.iter().copied().collect();
                            ids.sort_unstable();
                            let descriptors: Vec<Box<str>> = ids.iter()
                                .map(|id| transcript.get_dna_variant(*id).as_boxed_str())
                                .collect();
                            (!descriptors.is_empty())
                                .then(|| descriptors.join(LIST_SEPARATOR).into_boxed_str())
                        });
                    let preceding_event_rna_variant: Option<Box<str>> = transcript_nucleotide
                        .get_preceding_event_rna_variant_id()
                        .and_then(|id|{
                            Some(transcript.get_rna_variant(id).as_boxed_str())
                        });
                    let preceding_event_dna_variant: Option<Box<str>> = transcript_nucleotide
                        .get_preceding_event_dna_variant_ids()
                        .as_ref()
                        .and_then(|ids| {
                            let mut ids: Vec<u32> = ids.iter().copied().collect();
                            ids.sort_unstable();
                            let descriptors: Vec<Box<str>> = ids.iter()
                                .map(|id| transcript.get_dna_variant(*id).as_boxed_str())
                                .collect();
                            (!descriptors.is_empty())
                                .then(|| descriptors.join(LIST_SEPARATOR).into_boxed_str())
                        });

                    NucleotideRecord {
                        proteoform_id: proteoform.get_id(),
                        assembled_transcript_name: assembled_transcript_name.clone(),
                        amino_acid_index,
                        amino_acid: amino_acid_letter.clone(),
                        codon_index: codon_idx as u8,
                        nucleotide: transcript_nucleotide.get_nucleotide().as_str().into(),
                        is_amino_acid_variant,
                        is_nucleotide_variant: transcript_nucleotide.is_variant(),
                        assembled_transcript_read_position: transcript_nucleotide.get_transcript_read_position(),
                        assembled_transcript_alignment_index: transcript_nucleotide.get_transcript_alignment_index(),
                        assembled_transcript_variant_id: transcript_nucleotide.get_rna_variant_id(),
                        assembled_transcript_variant: rna_variant,
                        dna_variant_ids: transcript_nucleotide.get_dna_variant_ids().clone(),
                        dna_variant: dna_variant,
                        preceding_event_assembled_transcript_variant_id: transcript_nucleotide.get_preceding_event_rna_variant_id(),
                        preceding_event_assembled_transcript_variant: preceding_event_rna_variant,
                        preceding_event_dna_variant_ids: transcript_nucleotide.get_preceding_event_dna_variant_ids().clone(),
                        preceding_event_dna_variant: preceding_event_dna_variant,
                        is_reference_stitched: transcript_nucleotide.is_reference_stitched()
                    }
                })
            })
        })
    })
}

pub fn build_proteoform_records<'a>(
    transcript_set: &'a AssembledTranscriptSet
) -> impl Iterator<Item =ProteoformRecord> + 'a {
    transcript_set.iter().flat_map(move |transcript| {
        let assembled_transcript_name: Box<str> = transcript.get_assembled_transcript_name().into();
        let assembled_transcript_sequence: Box<str> = transcript.sequence.clone();
        let assembled_transcript_sequence_length: usize = assembled_transcript_sequence.len();
        let rna_read_names: Box<str> = join_box_strs(transcript.get_read_ids());
        let num_rna_read_names: usize = transcript.get_read_ids().len();
        let reference_gene_names: String = transcript
            .transcript_alignment.reference_gene_names.iter()
            .map(|s| s.as_ref()).collect::<Vec<&str>>().join(LIST_SEPARATOR);
        let reference_transcript_ids: String = transcript
            .transcript_alignment.reference_transcript_ids.iter()
            .map(|s| s.as_ref()).collect::<Vec<&str>>().join(LIST_SEPARATOR);

        transcript.proteoforms.iter().map(move |proteoform| {
            // Walk all amino_acids and nucleotides once to collect everything
            let mut rna_variant_id_set: HashSet<u32> = HashSet::new();
            let mut dna_variant_id_set: HashSet<u32> = HashSet::new();
            let mut variant_read_positions: Vec<u32> = Vec::new();
            let mut mutant_amino_acid_indices: Vec<u32> = Vec::new();
            let mut reference_stitched_amino_acid_indices: Vec<u32> = Vec::new();

            // Earliest amino acid carrying an RNA variant that is not known to
            // keep the reading frame (`AssembledTranscript::keeps_reading_frame`):
            // a frameshifting indel, a junction, a region call or a breakend.
            // Every codon after it may be read in another frame than the
            // reference's, so the whole downstream region is marked mutant, not
            // just the codons the variant lands in (which is all that
            // `is_variant()` flags directly). Where the frame did not in fact
            // change, the downstream residues are wild type, and call-peptide-vars
            // drops their windows against the reference proteome.
            let mut mutant_region_start_index: Option<u32> = None;
            for amino_acid in transcript.get_amino_acids(proteoform) {
                if amino_acid.is_variant() {
                    mutant_amino_acid_indices.push(amino_acid.get_index());
                }
                // Pushed in the walk's ascending index order, so no sort is needed
                // before interval formatting (and no downstream expansion applies:
                // stitched provenance is positional, not frame-dependent).
                if amino_acid.is_reference_stitched() {
                    reference_stitched_amino_acid_indices.push(amino_acid.get_index());
                }
                for nucleotide in amino_acid.nucleotides.iter() {
                    if let Some(id) = nucleotide.get_rna_variant_id() {
                        rna_variant_id_set.insert(id);
                    }
                    if let Some(ids) = nucleotide.get_dna_variant_ids().as_ref() {
                        dna_variant_id_set.extend(ids.iter().copied());
                    }
                    if nucleotide.is_variant() {
                        variant_read_positions.push(nucleotide.get_transcript_read_position());
                    }
                    // The walk is in ascending index order, so the first hit is
                    // the earliest. A variant is carried directly, or by the base
                    // immediately after it via a preceding event.
                    let opens_region: bool = [
                        nucleotide.get_rna_variant_id(),
                        nucleotide.get_preceding_event_rna_variant_id()
                    ].into_iter().flatten().any(|id| !transcript.keeps_reading_frame(id));
                    if opens_region && mutant_region_start_index.is_none() {
                        mutant_region_start_index = Some(amino_acid.get_index());
                    }
                }
            }

            // Expand the mutant set from that codon through the final residue of
            // the ORF.
            if let Some(start) = mutant_region_start_index {
                mutant_amino_acid_indices.extend(start..proteoform.get_length() as u32);
            }
            // Direct-hit indices were pushed in order, but the downstream
            // expansion can overlap and reorder them; sort + dedup so the
            // interval formatter sees a clean ascending sequence and the count
            // is accurate.
            mutant_amino_acid_indices.sort_unstable();
            mutant_amino_acid_indices.dedup();

            // Sort and dedup so all joined strings are deterministic
            let mut rna_ids_sorted: Vec<u32> = rna_variant_id_set.into_iter().collect();
            rna_ids_sorted.sort_unstable();
            let mut dna_ids_sorted: Vec<u32> = dna_variant_id_set.into_iter().collect();
            dna_ids_sorted.sort_unstable();
            variant_read_positions.sort_unstable();
            variant_read_positions.dedup();

            // Look up GOV descriptors for each variant id
            let rna_variants: String = rna_ids_sorted.iter()
                .map(|id| transcript.get_rna_variant(*id).as_boxed_str())
                .collect::<Vec<Box<str>>>()
                .join(LIST_SEPARATOR);
            let dna_variants: String = dna_ids_sorted.iter()
                .map(|id| transcript.get_dna_variant(*id).as_boxed_str())
                .collect::<Vec<Box<str>>>()
                .join(LIST_SEPARATOR);

            // DNA-variant read_names: dedup across all integrated DNA ids
            let mut dna_read_name_set: HashSet<&str> = HashSet::new();
            for id in dna_ids_sorted.iter() {
                if let Some(names) = transcript.dna_variant_read_names.get(id) {
                    // Not `split(',')`: `build_dna_variant_records` joins this field with
                    // `join_ids`, which is `LIST_SEPARATOR`. Splitting on the wrong character
                    // matched nothing, so each variant's whole read list entered the set as
                    // one blob — reads shared by two integrated DNA variants were then
                    // deduplicated at blob granularity and reported twice.
                    for name in split_read_names_borrowed(names) {
                        dna_read_name_set.insert(name);
                    }
                }
            }
            let mut dna_read_names_sorted: Vec<&str> =
                dna_read_name_set.into_iter().collect();
            dna_read_names_sorted.sort_unstable();
            let dna_variant_read_names: Box<str> =
                dna_read_names_sorted.join(LIST_SEPARATOR).into_boxed_str();

            let assembled_transcript_read_position: Box<str> = variant_read_positions
                .iter().map(|p| p.to_string())
                .collect::<Vec<String>>()
                .join(LIST_SEPARATOR)
                .into_boxed_str();

            ProteoformRecord {
                proteoform_id: proteoform.get_id() as usize,
                amino_acid_sequence: proteoform.get_sequence().to_string(),
                amino_acid_sequence_length: proteoform.get_length(),
                num_mutant_amino_acids: mutant_amino_acid_indices.len() as u32,
                assembled_transcript_name: assembled_transcript_name.clone(),
                assembled_transcript_sequence: assembled_transcript_sequence.clone(),
                assembled_transcript_sequence_length,
                orf_start: proteoform.get_orf_start(),
                orf_end: proteoform.get_orf_end(),
                reference_gene_name: reference_gene_names.clone(),
                reference_transcript_id: reference_transcript_ids.clone(),
                mutant_amino_acid_intervals: format_intervals(&mutant_amino_acid_indices),
                assembled_transcript_variant_ids: join_u32_ids(&rna_ids_sorted),
                assembled_transcript_variants: rna_variants,
                dna_variant_ids: join_u32_ids(&dna_ids_sorted),
                dna_variants,
                assembled_transcript_read_names: rna_read_names.clone(),
                num_assembled_transcript_read_names: num_rna_read_names,
                dna_variant_read_names,
                assembled_transcript_read_position: assembled_transcript_read_position,
                reference_stitched_amino_acid_intervals: format_intervals(&reference_stitched_amino_acid_indices),
                num_reference_stitched_amino_acids: reference_stitched_amino_acid_indices.len() as u32
            }
        })
    })
}
