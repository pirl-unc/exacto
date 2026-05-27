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
use exacto_core::prelude::Strand;
use exacto_integrator::prelude::*;
use std::collections::{HashMap, HashSet};
use std::str::FromStr;

use crate::prelude::*;


pub fn build_transcript_set(
    assembled_transcript_support_records: &Vec<AssembledTranscriptSupportRecord>,
    transcript_model_structure_records: &Vec<TranscriptModelStructureRecord>,
    rna_variant_records: &Vec<RNAVariantRecord>,
    dna_variant_records: &Vec<DNAVariantRecord>,
    integrated_variant_records: &Vec<IntegratedVariantRecord>
) -> TranscriptSet {
    // Step 1. Index records
    let assembled_transcript_support_index: HashMap<&str, &AssembledTranscriptSupportRecord> =
        assembled_transcript_support_records
            .iter()
            .map(|r| (r.assembled_transcript_name.as_ref(), r))
            .collect();
    let mut structure_records_index: HashMap<(Box<str>, u32), Vec<&TranscriptModelStructureRecord>> = HashMap::new();
    for r in transcript_model_structure_records {
        structure_records_index
            .entry((r.assembled_transcript_name.clone(), r.transcript_model_id))
            .or_default()
            .push(r);
    }
    let mut rna_variant_records_index: HashMap<(Box<str>, u32), Vec<&RNAVariantRecord>> = HashMap::new();
    for r in rna_variant_records {
        rna_variant_records_index
            .entry((r.assembled_transcript_name.clone(), r.transcript_model_id))
            .or_default()
            .push(r);
    }
    let mut integrated_variant_records_index: HashMap<(Box<str>, u32), Vec<&IntegratedVariantRecord>> = HashMap::new();
    for r in integrated_variant_records {
        integrated_variant_records_index
            .entry((r.assembled_transcript_name.clone(), r.transcript_model_id))
            .or_default()
            .push(r);
    }
    let dna_variant_records_by_id: HashMap<u32, &DNAVariantRecord> = dna_variant_records.iter().map(|r| (r.variant_id, r)).collect();

    // Step 2. Build one Transcript per unique (assembled_transcript_name, transcript_model_id)
    let mut keys: Vec<(Box<str>, u32)> = structure_records_index.keys().cloned().collect();
    keys.sort_unstable();
    let mut transcripts: Vec<Transcript> = Vec::with_capacity(keys.len());
    for (assembled_transcript_name, transcript_model_id) in keys {
        let structure_rows: &Vec<&TranscriptModelStructureRecord> = &structure_records_index[&(assembled_transcript_name.clone(), transcript_model_id)];
        let (sequence, read_ids): (Box<str>, Vec<Box<str>>) =
            match assembled_transcript_support_index.get(assembled_transcript_name.as_ref()) {
                Some(support) => {
                    let read_ids: Vec<Box<str>> = support.read_names
                        .split(',')
                        .filter(|s| !s.is_empty())
                        .map(|s| s.into())
                        .collect();
                    (support.sequence.clone(), read_ids)
                }
                None => {
                    // No assembled-transcript support row; concatenate
                    // structure-row sequences and use an empty read_ids list.
                    let concat: String = structure_rows.iter()
                        .map(|r| r.sequence.as_ref())
                        .collect();
                    (concat.into_boxed_str(), Vec::new())
                }
            };

        // TranscriptStructure
        let (reference_gene_names, reference_transcript_ids): (Vec<Box<str>>, Vec<Box<str>>) =
            structure_rows
                .first()
                .map(|r| {
                    let gene_names: Vec<Box<str>> = r.reference_gene_name
                        .split(',').filter(|s| !s.is_empty()).map(|s| s.into()).collect();
                    let transcript_ids: Vec<Box<str>> = r.reference_transcript_id
                        .split(',').filter(|s| !s.is_empty()).map(|s| s.into()).collect();
                    (gene_names, transcript_ids)
                })
                .unwrap_or_default();

        let mut transcript_structure: TranscriptStructure = TranscriptStructure::new(
            transcript_model_id as usize,
            reference_gene_names,
            reference_transcript_ids
        );

        for row in structure_rows.iter() {
            transcript_structure.add_item(build_transcript_structure_item(row));
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
        if let Some(rna_rows) = rna_variant_records_index.get(&(assembled_transcript_name.clone(), transcript_model_id)) {
            for r in rna_rows {
                let gov = graph_operation_view_from_rna_variant_record(r);
                let canonical_id: u32 = match rna_variants_bimap.get_by_right(&gov).copied() {
                    Some(existing) => existing,
                    None => {
                        rna_variants_bimap.insert(r.variant_id, gov);
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
        if let Some(records) = integrated_variant_records_index.get(&(assembled_transcript_name.clone(), transcript_model_id)) {
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
        // primary-structure record builders don't need the global DNA stream.
        //
        // Apply the same dedup-on-GOV / canonical-id pattern as for RNA, in
        // case two DNA variant records share a GOV (rare in practice but
        // would silently drop integration linkages if not handled).
        let mut dna_variants_bimap: BiMap<u32, GraphOperationView> = BiMap::new();
        let mut dna_variant_read_names: HashMap<u32, Box<str>> = HashMap::new();
        let mut dna_variant_id_remap: HashMap<u32, u32> = HashMap::new();
        for dna_id in dna_variant_ids_needed {
            if let Some(dna_record) = dna_variant_records_by_id.get(&dna_id) {
                let gov: GraphOperationView = graph_operation_view_from_dna_variant_record(dna_record);
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

        // Composite Transcript.id "{at_name}:{tm_id}" — unique per
        // (assembled_transcript_name, transcript_model_id) pair (the join key).
        let transcript_id: Box<str> = format!("{}:{}", assembled_transcript_name, transcript_model_id).into_boxed_str();

        transcripts.push(Transcript::new(
            transcript_id,
            sequence,
            read_ids,
            assembled_transcript_name.clone(),
            transcript_model_id,
            transcript_structure,
            rna_variants_bimap,
            dna_variants_bimap,
            integrated_variant_ids,
            dna_variant_read_names
        ));
    }

    TranscriptSet::new(transcripts)
}


fn build_transcript_structure_item(
    record: &TranscriptModelStructureRecord
) -> TranscriptStructureItem {
    let gov: GraphOperationView = graph_operation_view_from_descriptor(
        &record.chromosome_1,
        record.position_1,
        &record.operation_1,
        &record.strand_1,
        &record.chromosome_2,
        record.position_2,
        &record.operation_2,
        &record.strand_2,
        &record.sequence
    );

    let record_type = AlignmentStructureRecordType::from_str(&record.record_type)
        .unwrap_or_else(|_| panic!(
            "Unknown TranscriptModelStructureRecord.record_type: {:?}",
            record.record_type
        ));
    let item_type = match record_type {
        AlignmentStructureRecordType::Base => {
            let kind = AlignmentStructureBaseKind::from_str(&record.kind).unwrap();
            let context = AlignmentStructureBaseContext::from_str(&record.context).unwrap();
            TranscriptStructureItemType::Base { kind, context }
        }
        AlignmentStructureRecordType::Event => {
            let kind = AlignmentStructureEventKind::from_str(&record.kind).unwrap();
            let context = AlignmentStructureEventContext::from_str(&record.context).unwrap();
            TranscriptStructureItemType::Event { kind, context }
        }
    };

    let annotation = TranscriptStructureAnnotation {
        position_1_annotation: Annotation {
            gene_id: optional_id(&record.gene_id_1),
            transcript_id: optional_id(&record.transcript_id_1),
            exon_id: optional_id(&record.exon_id_1)
        },
        position_2_annotation: Annotation {
            gene_id: optional_id(&record.gene_id_2),
            transcript_id: optional_id(&record.transcript_id_2),
            exon_id: optional_id(&record.exon_id_2)
        },
    };

    TranscriptStructureItem::new(
        record.index,
        record.read_start,
        record.read_end,
        item_type,
        gov,
        annotation
    )
}


fn graph_operation_view_from_rna_variant_record(record: &RNAVariantRecord) -> GraphOperationView {
    graph_operation_view_from_descriptor(
        &record.chromosome_1,
        record.position_1,
        &record.operation_1,
        &record.strand_1,
        &record.chromosome_2,
        record.position_2,
        &record.operation_2,
        &record.strand_2,
        &record.sequence
    )
}


fn graph_operation_view_from_dna_variant_record(record: &DNAVariantRecord) -> GraphOperationView {
    graph_operation_view_from_descriptor(
        &record.chromosome_1,
        record.position_1,
        &record.operation_1,
        &record.strand_1,
        &record.chromosome_2,
        record.position_2,
        &record.operation_2,
        &record.strand_2,
        &record.sequence
    )
}


fn graph_operation_view_from_descriptor(
    chromosome_1: &str,
    position_1: u32,
    operation_1: &str,
    strand_1: &str,
    chromosome_2: &str,
    position_2: u32,
    operation_2: &str,
    strand_2: &str,
    sequence: &str
) -> GraphOperationView {
    let strand_1: Strand = Strand::from_str(strand_1).unwrap();
    let strand_2: Strand = Strand::from_str(strand_2).unwrap();
    let operation_type_1: GraphOperationType = GraphOperationType::from_str(operation_1).unwrap();
    let operation_type_2: GraphOperationType = GraphOperationType::from_str(operation_2).unwrap();
    GraphOperationView::new(
        chromosome_1,
        position_1,
        &operation_type_1,
        &strand_1,
        chromosome_2,
        position_2,
        &operation_type_2,
        &strand_2,
        sequence,
        None
    )
}


fn optional_id(s: &str) -> Option<Box<str>> {
    if s.is_empty() { None } else { Some(s.into()) }
}


pub fn build_nucleotide_records<'a>(
    transcript_set: &'a TranscriptSet
) -> impl Iterator<Item =NucleotideRecord> + 'a {
    transcript_set.iter().flat_map(move |transcript| {
        // Identity fields: pulled directly from the Transcript so the
        // record carries the two underlying keys (assembled-transcript name
        // and transcript-model id) rather than the composite Transcript.id.
        let assembled_transcript_name: Box<str> = transcript.get_assembled_transcript_name().into();
        let transcript_model_id: u32 = transcript.get_transcript_model_id();

        transcript.primary_structures.iter().flat_map(move |primary_structure| {
            let assembled_transcript_name = assembled_transcript_name.clone();
            primary_structure.amino_acids.iter().flat_map(move |amino_acid| {
                let assembled_transcript_name = assembled_transcript_name.clone();
                amino_acid.nucleotides.iter().enumerate().map(move |(codon_idx, transcript_nucleotide)| {
                    let rna_variant: Option<Box<str>> = transcript_nucleotide
                        .get_rna_variant_id()
                        .and_then(|id|{
                                Some(transcript.get_rna_variant(id).as_boxed_str())
                        });
                    let dna_variant: Option<Box<str>> = transcript_nucleotide
                        .get_dna_variant_ids()
                        .as_ref()
                        .and_then(|ids| {
                            let descriptors: Vec<Box<str>> = ids.iter()
                                .map(|id| transcript.get_dna_variant(*id).as_boxed_str())
                                .collect();
                            (!descriptors.is_empty())
                                .then(|| descriptors.join(",").into_boxed_str())
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
                            let descriptors: Vec<Box<str>> = ids.iter()
                                .map(|id| transcript.get_dna_variant(*id).as_boxed_str())
                                .collect();
                            (!descriptors.is_empty())
                                .then(|| descriptors.join(",").into_boxed_str())
                        });

                    NucleotideRecord {
                        primary_structure_id: primary_structure.get_id(),
                        assembled_transcript_name: assembled_transcript_name.clone(),
                        transcript_model_id,
                        amino_acid_index: amino_acid.get_index(),
                        amino_acid: amino_acid.get_amino_acid().into(),
                        codon_index: codon_idx as u8,
                        nucleotide: transcript_nucleotide.get_nucleotide().as_str().into(),
                        is_amino_acid_variant: amino_acid.is_variant(),
                        is_nucleotide_variant: transcript_nucleotide.is_variant(),
                        transcript_read_position: transcript_nucleotide.get_transcript_read_position(),
                        transcript_structure_index: transcript_nucleotide.get_transcript_structure_index(),
                        rna_variant_id: transcript_nucleotide.get_rna_variant_id(),
                        rna_variant: rna_variant,
                        dna_variant_ids: transcript_nucleotide.get_dna_variant_ids().clone(),
                        dna_variant: dna_variant,
                        preceding_event_rna_variant_id: transcript_nucleotide.get_preceding_event_rna_variant_id(),
                        preceding_event_rna_variant: preceding_event_rna_variant,
                        preceding_event_dna_variant_ids: transcript_nucleotide.get_preceding_event_dna_variant_ids().clone(),
                        preceding_event_dna_variant: preceding_event_dna_variant
                    }
                })
            })
        })
    })
}

pub fn build_primary_structure_records<'a>(
    transcript_set: &'a TranscriptSet
) -> impl Iterator<Item =PrimaryStructureRecord> + 'a {
    transcript_set.iter().flat_map(move |transcript| {
        let transcript_model_id: u32 = transcript.get_transcript_model_id();
        let assembled_transcript_name: Box<str> = transcript.get_assembled_transcript_name().into();
        let assembled_transcript_sequence: Box<str> = transcript.sequence.clone();
        let assembled_transcript_sequence_length: usize = assembled_transcript_sequence.len();
        let rna_read_names: Box<str> = join_box_strs(transcript.get_read_ids());
        let num_rna_read_names: usize = transcript.get_read_ids().len();
        let reference_gene_names: String = transcript
            .transcript_structure.reference_gene_names.iter()
            .map(|s| s.as_ref()).collect::<Vec<&str>>().join(",");
        let reference_transcript_ids: String = transcript
            .transcript_structure.reference_transcript_ids.iter()
            .map(|s| s.as_ref()).collect::<Vec<&str>>().join(",");

        transcript.primary_structures.iter().map(move |primary_structure| {
            // Walk all amino_acids and nucleotides once to collect everything
            let mut rna_variant_id_set: HashSet<u32> = HashSet::new();
            let mut dna_variant_id_set: HashSet<u32> = HashSet::new();
            let mut variant_read_positions: Vec<u32> = Vec::new();
            let mut mutant_amino_acid_indices: Vec<u32> = Vec::new();
            let mut amino_acid_sequence: String =
                String::with_capacity(primary_structure.amino_acids.len());

            for amino_acid in primary_structure.amino_acids.iter() {
                amino_acid_sequence.push_str(amino_acid.get_amino_acid());
                if amino_acid.is_variant() {
                    mutant_amino_acid_indices.push(amino_acid.get_index());
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
                }
            }

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
                .join(",");
            let dna_variants: String = dna_ids_sorted.iter()
                .map(|id| transcript.get_dna_variant(*id).as_boxed_str())
                .collect::<Vec<Box<str>>>()
                .join(",");

            // DNA-variant read_names: dedup across all integrated DNA ids
            let mut dna_read_name_set: HashSet<&str> = HashSet::new();
            for id in dna_ids_sorted.iter() {
                if let Some(names) = transcript.dna_variant_read_names.get(id) {
                    for name in names.split(',').filter(|s| !s.is_empty()) {
                        dna_read_name_set.insert(name);
                    }
                }
            }
            let mut dna_read_names_sorted: Vec<&str> =
                dna_read_name_set.into_iter().collect();
            dna_read_names_sorted.sort_unstable();
            let dna_variant_read_names: Box<str> =
                dna_read_names_sorted.join(",").into_boxed_str();

            let transcript_read_position: Box<str> = variant_read_positions
                .iter().map(|p| p.to_string())
                .collect::<Vec<String>>()
                .join(",")
                .into_boxed_str();

            PrimaryStructureRecord {
                primary_structure_id: primary_structure.get_id() as usize,
                amino_acid_sequence,
                amino_acid_sequence_length: primary_structure.get_length(),
                num_mutant_amino_acids: mutant_amino_acid_indices.len() as u32,
                assembled_transcript_name: assembled_transcript_name.clone(),
                assembled_transcript_sequence: assembled_transcript_sequence.clone(),
                assembled_transcript_sequence_length,
                transcript_model_id,
                orf_start: primary_structure.get_orf_start(),
                orf_end: primary_structure.get_orf_end(),
                reference_gene_names: reference_gene_names.clone(),
                reference_transcript_ids: reference_transcript_ids.clone(),
                mutant_amino_acid_intervals: format_intervals(&mutant_amino_acid_indices),
                rna_variant_ids: join_u32_ids(&rna_ids_sorted),
                rna_variants,
                dna_variant_ids: join_u32_ids(&dna_ids_sorted),
                dna_variants,
                rna_read_names: rna_read_names.clone(),
                num_rna_read_names,
                dna_variant_read_names,
                transcript_read_position,
            }
        })
    })
}
