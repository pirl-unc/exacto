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


use exacto_core::prelude::join_ids;
use std::collections::{HashMap, HashSet};

use crate::prelude::*;


pub fn build_assembled_transcript_records<'a>(
    tms: &'a TranscriptModelSet
) -> impl Iterator<Item = AssembledTranscriptRecord> + 'a {
    assert!(
        !tms.chromosome_names_map.is_empty(),
        "tms.chromosome_names_map is empty."
    );

    tms.assembled_transcripts.iter().map(move |at| {
        let (start_chr_id, start) = at.get_reference_start_position();
        let (end_chr_id, end) = at.get_reference_end_position();
        let start_chromosome: Box<str> = tms
            .chromosome_names_map
            .get_by_right(&start_chr_id)
            .unwrap()
            .clone();
        let end_chromosome: Box<str> = tms
            .chromosome_names_map
            .get_by_right(&end_chr_id)
            .unwrap()
            .clone();
        AssembledTranscriptRecord {
            assembled_transcript_name: at.get_name().into(),
            start_chromosome,
            start: start,
            end_chromosome,
            end: end,
            num_exons: at.get_exons().len() as u32,
            num_introns: at.get_introns().len() as u32
        }
    })
}


pub fn build_dna_variant_records<'a>(
    dvcs: &'a DNAVariantCallSet
) -> impl Iterator<Item = DNAVariantRecord> + 'a {
    assert!(
        !dvcs.chromosome_names_map.is_empty(),
        "dvcs.chromosome_names_map is empty."
    );
    assert!(
        !dvcs.read_names_map.is_empty(),
        "dvcs.read_names_map is empty."
    );

    dvcs.variant_calls.iter().map(move |(variant_id, variant_call)| {
        let (consensus_record, consensus_read_names) = variant_call.get_named_consensus_record(&dvcs.read_names_map);
        let chromosome_1: &str = dvcs
            .chromosome_names_map
            .get_by_right(&consensus_record.get_chromosome_1())
            .unwrap();
        let chromosome_2: &str = dvcs
            .chromosome_names_map
            .get_by_right(&consensus_record.get_chromosome_2())
            .unwrap();
        let read_names: Vec<&str> = variant_call
            .get_read_ids()
            .iter()
            .map(|read_id| dvcs.read_names_map.get_by_right(read_id).unwrap().as_ref())
            .collect();

        let num_consensus_read_names: u32 = consensus_read_names.len() as u32;
        let num_read_names: u32 = read_names.len() as u32;

        DNAVariantRecord {
            variant_id: *variant_id as u32,
            chromosome_1: chromosome_1.into(),
            position_1: consensus_record.get_position_1(),
            strand_1: consensus_record.get_strand_1().as_str().into(),
            operation_1: consensus_record.get_operation_1().as_str().into(),
            chromosome_2: chromosome_2.into(),
            position_2: consensus_record.get_position_2(),
            strand_2: consensus_record.get_strand_2().as_str().into(),
            operation_2: consensus_record.get_operation_2().as_str().into(),
            sequence: consensus_record.get_standardized_sequence().into(),
            variant_size: consensus_record.get_variant_size().abs() as u32,
            variant_type: consensus_record.get_variant_type().as_str().into(),
            consensus_read_names: join_ids(consensus_read_names).into(),
            num_consensus_read_names: num_consensus_read_names,
            read_names: join_ids(read_names).into(),
            num_read_names: num_read_names
        }
    })
}

pub fn build_exon_records<'a>(
    tms: &'a TranscriptModelSet
) -> impl Iterator<Item = ExonRecord> + 'a {
    assert!(
        !tms.chromosome_names_map.is_empty(),
        "tms.chromosome_names_map is empty."
    );
    assert!(
        !tms.read_names_map.is_empty(),
        "tms.read_names_map is empty."
    );

    tms.assembled_transcripts.iter().flat_map(move |assembled_transcript| {
        let assembled_transcript_name: Box<str> = assembled_transcript.get_name().into();
        assembled_transcript.get_exons().iter().map(move |exon| {
            let chromosome: Box<str> = tms
                .chromosome_names_map
                .get_by_right(&exon.reference_chromosome_id)
                .unwrap()
                .clone();
            ExonRecord {
                assembled_transcript_name: assembled_transcript_name.clone(),
                chromosome,
                start: exon.reference_start,
                end: exon.reference_end,
                exon_number: exon.exon_number as u32,
                strand: exon.reference_strand.as_str().into()
            }
        })
    })
}

pub fn build_intron_records<'a>(
    tms: &'a TranscriptModelSet
) -> impl Iterator<Item = IntronRecord> + 'a {
    assert!(
        !tms.chromosome_names_map.is_empty(),
        "tms.chromosome_names_map is empty."
    );
    assert!(
        !tms.read_names_map.is_empty(),
        "tms.read_names_map is empty."
    );

    tms.assembled_transcripts.iter().flat_map(move |assembled_transcript| {
        let assembled_transcript_name: Box<str> = assembled_transcript.get_name().into();
        assembled_transcript.get_introns().iter().map(move |intron| {
            let chromosome: Box<str> = tms
                .chromosome_names_map
                .get_by_right(&intron.reference_chromosome_id)
                .unwrap()
                .clone();
            IntronRecord {
                assembled_transcript_name: assembled_transcript_name.clone(),
                chromosome,
                start: intron.reference_start,
                end: intron.reference_end,
                intron_number: intron.intron_number as u32,
                strand: intron.reference_strand.as_str().into()
            }
        })
    })
}


pub fn build_reference_transcript_match_records<'a>(
    tms: &'a TranscriptModelSet
) -> impl Iterator<Item = ReferenceTranscriptMatchRecord> + 'a {
    let mut rows: Vec<ReferenceTranscriptMatchRecord> = Vec::new();

    for tm in tms.transcript_models.iter() {
        let assembled_transcript_name: Box<str> = tm.get_assembled_transcript_name().into();
        let transcript_model_id: u32 = tm.get_id() as u32;
        for reference_transcript_match in tm.get_reference_transcript_matches().iter() {
            rows.push(ReferenceTranscriptMatchRecord {
                assembled_transcript_name: assembled_transcript_name.clone(),
                transcript_model_id,
                reference_gene_name: reference_transcript_match.reference_gene_name.clone(),
                reference_transcript_id: reference_transcript_match.reference_transcript_id.clone(),
                num_overlap_bases: reference_transcript_match.num_overlap_bases,
                num_transcript_only_bases: reference_transcript_match.num_transcript_only_bases,
                num_reference_transcript_only_bases: reference_transcript_match.num_reference_only_bases,
                score: reference_transcript_match.score,
                scoring_method: reference_transcript_match.scoring_method.as_str().into(),
            });
        }
    }

    rows.into_iter()
}


pub fn build_read_filter_status_records<'a>(
    tms: &'a TranscriptModelSet
) -> impl Iterator<Item = ReadFilterStatusRecord> + 'a {
    let included_read_names: HashSet<Box<str>> = tms
        .assembled_transcripts
        .iter()
        .map(|at| at.get_name().into())
        .collect();

    tms.read_names_map
        .iter()
        .map(move |(read_name, _read_id)| ReadFilterStatusRecord {
            read_name: read_name.clone(),
            excluded: !included_read_names.contains(read_name)
        })
}



pub fn build_transcript_model_structure_records<'a>(
    tms: &'a TranscriptModelSet
) -> impl Iterator<Item = TranscriptModelStructureRecord> + 'a {
    assert!(
        !tms.chromosome_names_map.is_empty(),
        "tms.chromosome_names_map is empty."
    );

    tms.transcript_models.iter().flat_map(move |tm| {
        let assembled_transcript_name: Box<str> = tm.get_assembled_transcript_name().into();

        let transcript_id_to_gene_name: HashMap<String, String> = tm
            .get_reference_transcript_matches()
            .iter()
            .map(|m| {
                (
                    m.reference_transcript_id.to_string(),
                    m.reference_gene_name.to_string(),
                )
            })
            .collect();

        let reference_transcript_ids: Vec<&str> = tm
            .get_reference_transcript_ids()
            .iter()
            .map(|s| s.as_ref())
            .collect();
        let reference_transcript_id_cell: Box<str> = join_ids(&reference_transcript_ids).into();
        let reference_gene_name_cell: Box<str> = reference_transcript_ids
            .iter()
            .map(|tid| {
                transcript_id_to_gene_name
                    .get(*tid)
                    .map(|gn| gn.as_str())
                    .unwrap_or("")
            })
            .collect::<Vec<&str>>()
            .join(",")
            .into();

        let chromosome_names_map = &tms.chromosome_names_map;

        let records: Vec<AlignmentStructureRecord> = tm.get_alignment_structure().identify_records();
        records.into_iter().enumerate().map(move |(index, record)| {
            let chromosome_1: Box<str> = chromosome_names_map
                .get_by_right(&record.get_chromosome_1())
                .unwrap()
                .clone();
            let chromosome_2: Box<str> = chromosome_names_map
                .get_by_right(&record.get_chromosome_2())
                .unwrap()
                .clone();
            TranscriptModelStructureRecord {
                assembled_transcript_name: assembled_transcript_name.clone(),
                transcript_model_id: tm.get_id() as u32,
                reference_gene_name: reference_gene_name_cell.clone(),
                reference_transcript_id: reference_transcript_id_cell.clone(),
                index: index as u32,
                read_start: record.get_start(),
                read_end: record.get_end(),
                sequence: record.get_standardized_sequence(),
                record_type: record.get_record_type().as_str().into(),
                kind: record.get_kind().as_str().into(),
                context: record
                    .get_context()
                    .as_ref()
                    .map(|c| c.as_str().into())
                    .unwrap_or_else(|| "".into()),
                chromosome_1,
                position_1: record.get_position_1(),
                operation_1: record.get_operation_1().as_str().into(),
                strand_1: record.get_strand_1().as_str().into(),
                chromosome_2,
                position_2: record.get_position_2(),
                operation_2: record.get_operation_2().as_str().into(),
                strand_2: record.get_strand_2().as_str().into(),
                gene_id_1: record
                    .get_gene_id_1()
                    .as_ref()
                    .map(|k| k.clone())
                    .unwrap_or_else(|| "".into()),
                transcript_id_1: record
                    .get_transcript_id_1()
                    .as_ref()
                    .map(|k| k.clone())
                    .unwrap_or_else(|| "".into()),
                exon_id_1: record
                    .get_exon_id_1()
                    .as_ref()
                    .map(|k| k.clone())
                    .unwrap_or_else(|| "".into()),
                gene_id_2: record
                    .get_gene_id_2()
                    .as_ref()
                    .map(|k| k.clone())
                    .unwrap_or_else(|| "".into()),
                transcript_id_2: record
                    .get_transcript_id_2()
                    .as_ref()
                    .map(|k| k.clone())
                    .unwrap_or_else(|| "".into()),
                exon_id_2: record
                    .get_exon_id_2()
                    .as_ref()
                    .map(|k| k.clone())
                    .unwrap_or_else(|| "".into()),
                skipped: record.get_skipped_string(chromosome_names_map).into(),
            }
        }).collect::<Vec<_>>().into_iter()
    })
}


pub fn build_rna_variant_records<'a>(
    tms: &'a TranscriptModelSet
) -> impl Iterator<Item = RNAVariantRecord> + 'a {
    assert!(
        !tms.chromosome_names_map.is_empty(),
        "tms.chromosome_names_map is empty."
    );
    assert!(
        !tms.read_names_map.is_empty(),
        "tms.read_names_map is empty."
    );

    // gene_name lookup is per-AT (keyed by AT name). Per-tm matches resolve via this.
    let mut gene_name_lookup: HashMap<Box<str>, HashMap<String, String>> = HashMap::new();
    for tm in tms.transcript_models.iter() {
        let inner = gene_name_lookup
            .entry(tm.get_assembled_transcript_name().into())
            .or_insert_with(HashMap::new);
        for m in tm.get_reference_transcript_matches().iter() {
            inner
                .entry(m.reference_transcript_id.to_string())
                .or_insert_with(|| m.reference_gene_name.to_string());
        }
    }

    tms.transcript_models
        .iter()
        .flat_map(move |tm| {
            let assembled_transcript_name: Box<str> = tm.get_assembled_transcript_name().into();
            let reference_transcript_ids: Vec<&str> = tm
                .get_reference_transcript_ids()
                .iter()
                .map(|s| s.as_ref())
                .collect();
            let reference_transcript_id_cell: Box<str> = join_ids(&reference_transcript_ids).into();
            let reference_gene_name_cell: Box<str> = match gene_name_lookup.get(&assembled_transcript_name) {
                Some(inner) => reference_transcript_ids
                    .iter()
                    .map(|tid| inner.get(*tid).map(|gn| gn.as_str()).unwrap_or(""))
                    .collect::<Vec<&str>>()
                    .join(",")
                    .into(),
                None => "".into(),
            };

            tm.get_variant_records().iter().map(move |variant_record| {
                let chromosome_1: Box<str> = tms
                    .chromosome_names_map
                    .get_by_right(&variant_record.get_chromosome_1())
                    .unwrap()
                    .clone();
                let chromosome_2: Box<str> = tms
                    .chromosome_names_map
                    .get_by_right(&variant_record.get_chromosome_2())
                    .unwrap()
                    .clone();
                RNAVariantRecord {
                    variant_id: 0, // overwritten below via scan
                    assembled_transcript_name: assembled_transcript_name.clone(),
                    transcript_model_id: tm.get_id() as u32,
                    reference_gene_name: reference_gene_name_cell.clone(),
                    reference_transcript_id: reference_transcript_id_cell.clone(),
                    chromosome_1,
                    position_1: variant_record.get_graph_operation().get_position_1(),
                    strand_1: variant_record.get_graph_operation().get_strand_1().as_str().into(),
                    operation_1: variant_record.get_graph_operation().get_operation_type_1().as_str().into(),
                    chromosome_2,
                    position_2: variant_record.get_graph_operation().get_position_2(),
                    strand_2: variant_record.get_graph_operation().get_strand_2().as_str().into(),
                    operation_2: variant_record.get_graph_operation().get_operation_type_2().as_str().into(),
                    variant_size: variant_record.get_variant_size().abs() as u32,
                    variant_type: variant_record.get_variant_type().as_str().into(),
                    sequence: variant_record.get_graph_operation().get_standardized_sequence().into(),
                    read_start: variant_record.get_read_position_1(),
                    read_end: variant_record.get_read_position_2()
                }
            })
        })
        .scan(1u32, |next_id, mut record| {
            record.variant_id = *next_id;
            *next_id += 1;
            Some(record)
        })
}
