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


use exacto_core::prelude::{
    get_fasta_sequence, 
    get_fasta_sequence_ids, 
    join_ids, 
    LIST_SEPARATOR
};
use exacto_core::prelude::{ReadID, ReadName, ReadSupport, ReferenceChromosomeName, ReferencePosition, ReferenceTranscriptID, VariantID};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::rc::Rc;

use crate::prelude::*;


fn build_junction_homology_fetcher(
    reference_genome_fasta_file: &str
) -> impl Fn(&str, ReferencePosition, ReferencePosition) -> Option<Box<str>> {
    let contig_lengths: HashMap<ReferenceChromosomeName, u32> =
        get_fasta_sequence_ids(reference_genome_fasta_file)
            .into_iter()
            .collect();
    let fasta_file: String = reference_genome_fasta_file.to_string();
    move |chromosome: &str, start: ReferencePosition, end: ReferencePosition| -> Option<Box<str>> {
        let contig_length: u32 = *contig_lengths.get(chromosome)?;
        let clamped_start: ReferencePosition = start.max(1);
        let clamped_end: ReferencePosition = end.min(contig_length);
        if clamped_start > clamped_end {
            return Some("".into());
        }
        Some(get_fasta_sequence(chromosome, clamped_start, clamped_end, &fasta_file))
    }
}


pub fn build_assembled_transcript_records<'a>(
    tms: &'a TranscriptModelSet
) -> impl Iterator<Item = AssembledTranscriptRecord> + 'a {
    assert!(
        !tms.chromosome_names_map.is_empty(),
        "tms.chromosome_names_map is empty."
    );
    tms.transcript_models.iter().map(move |tm| {
        let read_name: ReadName = tms.read_names_map.get_by_right(&tm.get_read_id()).unwrap().clone();
        let (start_chr_id, start) = tm.get_reference_start_position();
        let (end_chr_id, end) = tm.get_reference_end_position();
        let start_chromosome: ReferenceChromosomeName = tms
            .chromosome_names_map
            .get_by_right(&start_chr_id)
            .unwrap()
            .clone();
        let end_chromosome: ReferenceChromosomeName = tms
            .chromosome_names_map
            .get_by_right(&end_chr_id)
            .unwrap()
            .clone();
        AssembledTranscriptRecord {
            assembled_transcript_name: read_name,
            start_chromosome,
            start: start,
            end_chromosome,
            end: end,
            num_exons: tm.get_exons().len() as u32,
            num_splice_junctions: tm.get_splice_junctions().len() as u32,
            is_variant: !tm.is_reference_transcript()
        }
    })
}


pub fn build_assembled_transcript_exon_records<'a>(
    tms: &'a TranscriptModelSet
) -> impl Iterator<Item =AssembledTranscriptExonRecord> + 'a {
    assert!(
        !tms.chromosome_names_map.is_empty(),
        "tms.chromosome_names_map is empty."
    );
    assert!(
        !tms.read_names_map.is_empty(),
        "tms.read_names_map is empty."
    );
    tms.transcript_models.iter().flat_map(move |tm| {
        let read_name: ReadName = tms.read_names_map.get_by_right(&tm.get_read_id()).unwrap().clone();
        tm.get_exons().iter().map(move |exon| {
            let chromosome: ReferenceChromosomeName = tms
                .chromosome_names_map
                .get_by_right(&exon.reference_chromosome_id)
                .unwrap()
                .clone();
            AssembledTranscriptExonRecord {
                assembled_transcript_name: read_name.clone(),
                chromosome,
                start: exon.reference_start,
                end: exon.reference_end,
                exon_number: exon.exon_number as u32,
                strand: exon.reference_strand.as_str().into()
            }
        })
    })
}


pub fn build_assembled_transcript_splice_junction_records<'a>(
    tms: &'a TranscriptModelSet
) -> impl Iterator<Item =AssembledTranscriptSpliceJunctionRecord> + 'a {
    assert!(
        !tms.chromosome_names_map.is_empty(),
        "tms.chromosome_names_map is empty."
    );
    assert!(
        !tms.read_names_map.is_empty(),
        "tms.read_names_map is empty."
    );
    tms.transcript_models.iter().flat_map(move |tm| {
        let read_name: ReadName = tms.read_names_map.get_by_right(&tm.get_read_id()).unwrap().clone();
        tm.get_splice_junctions().iter().map(move |splice_junction| {
            let chromosome_1: ReferenceChromosomeName = tms
                .chromosome_names_map
                .get_by_right(&splice_junction.reference_chromosome_id_1)
                .unwrap()
                .clone();
            let chromosome_2: ReferenceChromosomeName = tms
                .chromosome_names_map
                .get_by_right(&splice_junction.reference_chromosome_id_2)
                .unwrap()
                .clone();
            AssembledTranscriptSpliceJunctionRecord {
                assembled_transcript_name: read_name.clone(),
                chromosome_1: chromosome_1,
                chromosome_2: chromosome_2,
                position_1: splice_junction.reference_position_1,
                position_2: splice_junction.reference_position_2,
                strand_1: splice_junction.reference_strand_1.as_str().into(),
                strand_2: splice_junction.reference_strand_2.as_str().into(),
                splice_junction_number: splice_junction.splice_junction_number as u32
            }
        })
    })
}


pub fn build_assembled_transcript_reference_transcript_match_records<'a>(
    tms: &'a TranscriptModelSet
) -> impl Iterator<Item =AssembledTranscriptReferenceTranscriptMatchRecord> + 'a {
    let mut rows: Vec<AssembledTranscriptReferenceTranscriptMatchRecord> = Vec::new();
    for tm in tms.transcript_models.iter() {
        let read_name: ReadName = tms.read_names_map.get_by_right(&tm.get_read_id()).unwrap().clone();
        for reference_transcript_match in tm.get_reference_transcript_matches().iter() {
            rows.push(AssembledTranscriptReferenceTranscriptMatchRecord {
                assembled_transcript_name: read_name.clone(),
                reference_gene_id: reference_transcript_match.get_reference_gene_id().into(),
                reference_gene_name: reference_transcript_match.get_reference_gene_name().into(),
                reference_transcript_id: reference_transcript_match.get_reference_transcript_id().into(),
                num_splice_junction_matches: reference_transcript_match.num_splice_junction_matches() as u32,
                num_overlapping_bases: reference_transcript_match.num_overlapping_bases(),
                num_reference_only_bases: reference_transcript_match.num_reference_only_bases(),
                num_query_only_bases: reference_transcript_match.num_query_only_bases(),
                num_terminal_query_only_bases: reference_transcript_match.num_terminal_query_only_bases(),
                num_internal_query_only_bases: reference_transcript_match.num_internal_query_only_bases()
            });
        }
    }

    rows.into_iter()
}


pub fn build_assembled_transcript_nonsense_mediated_decay_records<'a>(
    tms: &'a TranscriptModelSet
) -> impl Iterator<Item = AssembledTranscriptNonsenseMediatedDecayRecord> + 'a {
    tms.transcript_models.iter().flat_map(move |tm| {
        let read_name: ReadName = tms.read_names_map.get_by_right(&tm.get_read_id()).unwrap().clone();
        tm.get_nmd_predictions().iter().map(move |call| {
            let (nmd_predicted, distance_to_last_junction): (bool, Option<u32>) = match call.verdict {
                NonsenseMediatedDecayVerdict::Predicted { distance_to_last_junction } => (true, Some(distance_to_last_junction)),
                NonsenseMediatedDecayVerdict::NotPredicted { distance_to_last_junction } => (false, distance_to_last_junction)
            };
            AssembledTranscriptNonsenseMediatedDecayRecord {
                assembled_transcript_name: read_name.clone(),
                orf_start: call.orf_start,
                orf_end: call.orf_end,
                nmd_predicted,
                distance_to_last_junction
            }
        })
    })
}


pub fn build_assembled_transcript_filter_status_records<'a>(
    tms: &'a TranscriptModelSet
) -> impl Iterator<Item =AssembledTranscriptFilterStatusRecord> + 'a {
    let included_read_names: HashSet<ReadName> = tms
        .transcript_models
        .iter()
        .map(|tm| tms.read_names_map.get_by_right(&tm.get_read_id()).unwrap().clone())
        .collect();
    let mut read_names: Vec<(&ReadID, &ReadName)> = tms.read_names_map
        .iter()
        .map(|(read_name, read_id)| (read_id, read_name))
        .collect();
    read_names.sort_unstable();
    read_names
        .into_iter()
        .map(move |(_read_id, read_name)| AssembledTranscriptFilterStatusRecord {
            read_name: read_name.clone(),
            excluded: !included_read_names.contains(read_name)
        })
}


pub fn build_assembled_transcript_alignment_records<'a>(
    tms: &'a TranscriptModelSet
) -> impl Iterator<Item =AssembledTranscriptModelAlignmentRecord> + 'a {
    assert!(
        !tms.chromosome_names_map.is_empty(),
        "tms.chromosome_names_map is empty."
    );

    tms.transcript_models.iter().flat_map(move |tm| {
        let read_name: ReadName = tms.read_names_map.get_by_right(&tm.get_read_id()).unwrap().clone();

        let transcript_id_to_gene_name: HashMap<String, String> = tm
            .get_reference_transcript_matches()
            .iter()
            .map(|m| {
                (
                    m.get_reference_transcript_id().to_string(),
                    m.get_reference_gene_name().to_string(),
                )
            })
            .collect();

        let reference_transcript_ids: Vec<ReferenceTranscriptID> = tm
            .get_reference_transcript_ids()
            .into_iter()
            .collect();
        let reference_transcript_id_cell: Box<str> = join_ids(&reference_transcript_ids).into();
        let reference_gene_name_cell: Box<str> = reference_transcript_ids
            .iter()
            .map(|tid| {
                transcript_id_to_gene_name
                    .get(&tid.to_string())
                    .map(|gn| gn.as_str())
                    .unwrap_or("")
            })
            .collect::<Vec<&str>>()
            .join(LIST_SEPARATOR)
            .into();

        let chromosome_names_map = &tms.chromosome_names_map;

        let records: Vec<AlignmentModelRecord> = identify_alignment_model_records(tm.get_alignment_model(), Some(tm.get_annotation()));
        records.into_iter().enumerate().map(move |(index, record)| {
            let chromosome_1: ReferenceChromosomeName = chromosome_names_map
                .get_by_right(&record.get_chromosome_1())
                .unwrap()
                .clone();
            let chromosome_2: ReferenceChromosomeName = chromosome_names_map
                .get_by_right(&record.get_chromosome_2())
                .unwrap()
                .clone();
            AssembledTranscriptModelAlignmentRecord {
                assembled_transcript_name: read_name.clone(),
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
                reference_gene_id_1: record
                    .get_gene_id_1()
                    .as_ref()
                    .map(|k| k.clone())
                    .unwrap_or_else(|| "".into()),
                reference_transcript_id_1: record
                    .get_transcript_id_1()
                    .as_ref()
                    .map(|k| k.clone())
                    .unwrap_or_else(|| "".into()),
                reference_exon_id_1: record
                    .get_exon_id_1()
                    .as_ref()
                    .map(|k| k.clone())
                    .unwrap_or_else(|| "".into()),
                reference_gene_id_2: record
                    .get_gene_id_2()
                    .as_ref()
                    .map(|k| k.clone())
                    .unwrap_or_else(|| "".into()),
                reference_transcript_id_2: record
                    .get_transcript_id_2()
                    .as_ref()
                    .map(|k| k.clone())
                    .unwrap_or_else(|| "".into()),
                reference_exon_id_2: record
                    .get_exon_id_2()
                    .as_ref()
                    .map(|k| k.clone())
                    .unwrap_or_else(|| "".into()),
                skipped: record.get_skipped_string(chromosome_names_map).into(),
            }
        }).collect::<Vec<_>>().into_iter()
    })
}


pub fn build_assembled_transcript_variant_records<'a>(
    tms: &'a TranscriptModelSet,
    reference_genome_fasta_file: Option<&str>,
    flank: u32
) -> impl Iterator<Item =AssembledTranscriptVariantRecord> + 'a {
    assert!(
        !tms.chromosome_names_map.is_empty(),
        "tms.chromosome_names_map is empty."
    );
    assert!(
        !tms.read_names_map.is_empty(),
        "tms.read_names_map is empty."
    );

    // `None` FASTA means the homology columns are all null; `Rc` because the fetcher
    // is shared into each per-transcript-model inner closure.
    let fetch: Option<Rc<dyn Fn(&str, ReferencePosition, ReferencePosition) -> Option<Box<str>>>> =
        reference_genome_fasta_file.map(|fasta_file| {
            Rc::new(build_junction_homology_fetcher(fasta_file))
                as Rc<dyn Fn(&str, ReferencePosition, ReferencePosition) -> Option<Box<str>>>
        });

    // gene_name lookup is per-AT (keyed by AT name). Per-tm matches resolve via this.
    let mut gene_name_lookup: HashMap<ReadName, HashMap<String, String>> = HashMap::new();
    for tm in tms.transcript_models.iter() {
        let read_name: ReadName = tms.read_names_map.get_by_right(&tm.get_read_id()).unwrap().clone();
        let inner = gene_name_lookup
            .entry(read_name)
            .or_insert_with(HashMap::new);
        for m in tm.get_reference_transcript_matches().iter() {
            inner
                .entry(m.get_reference_transcript_id().to_string())
                .or_insert_with(|| m.get_reference_gene_name().to_string());
        }
    }

    tms.transcript_models
        .iter()
        .flat_map(move |tm| {
            let fetch: Option<Rc<dyn Fn(&str, ReferencePosition, ReferencePosition) -> Option<Box<str>>>> = fetch.clone();
            let assembled_transcript_name: ReadName = tms.read_names_map
                .get_by_right(&tm.get_read_id())
                .unwrap()
                .clone();
            let reference_transcript_ids: Vec<ReferenceTranscriptID> = tm
                .get_reference_transcript_ids()
                .into_iter()
                .collect();
            let reference_transcript_id_cell: Box<str> = join_ids(&reference_transcript_ids).into();
            let reference_gene_name_cell: Box<str> = match gene_name_lookup.get(&assembled_transcript_name) {
                Some(inner) => reference_transcript_ids
                    .iter()
                    .map(|tid| inner.get(&tid.to_string()).map(|gn| gn.as_str()).unwrap_or(""))
                    .collect::<Vec<&str>>()
                    .join(LIST_SEPARATOR)
                    .into(),
                None => "".into(),
            };

            // ExonTruncation records are deliberately not published.
            tm.get_variant_records()
                .iter()
                .filter(|variant_record| {
                    *variant_record.get_variant_type() != VariantType::ExonTruncation
                })
                .map(move |variant_record| {
                    let chromosome_1: ReferenceChromosomeName = tms
                        .chromosome_names_map
                        .get_by_right(&variant_record.get_chromosome_1())
                        .unwrap()
                        .clone();
                    let chromosome_2: ReferenceChromosomeName = tms
                        .chromosome_names_map
                        .get_by_right(&variant_record.get_chromosome_2())
                        .unwrap()
                        .clone();
                    AssembledTranscriptVariantRecord {
                        variant_id: 0, // overwritten below via scan
                        assembled_transcript_name: assembled_transcript_name.clone(),
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
                        variant_size: u32::try_from(variant_record.get_variant_size()).ok(),
                        variant_type: variant_record.get_variant_type().as_str().into(),
                        sequence: variant_record.get_graph_operation().get_standardized_sequence().into(),
                        read_start: variant_record.get_read_position_1(),
                        read_end: variant_record.get_read_position_2(),
                        origin: "".into()
                    }
                })
        })
        .scan(1u32, |next_id, mut record| {
            record.variant_id = *next_id;
            *next_id += 1;
            Some(record)
        })
}


pub fn annotate_dna_variant_origins<'a, I>(
    records: I,
    dna_variant_records: &'a [DNAVariantRecord],
    buffer: u32
) -> impl Iterator<Item = AssembledTranscriptVariantRecord> + 'a
where
    I: IntoIterator<Item = AssembledTranscriptVariantRecord>,
    I::IntoIter: 'a
{
    // (chromosome_1, position_1) -> DNA variants, in the order of the tables.
    let mut index: BTreeMap<(ReferenceChromosomeName, ReferencePosition), Vec<&'a DNAVariantRecord>> = BTreeMap::new();
    for dna_variant_record in dna_variant_records.iter() {
        index
            .entry((dna_variant_record.chromosome_1.clone(), dna_variant_record.position_1))
            .or_default()
            .push(dna_variant_record);
    }

    records.into_iter().map(move |mut record| {
        let lo: (ReferenceChromosomeName, ReferencePosition) = (record.chromosome_1.clone(), record.position_1.saturating_sub(buffer));
        let hi: (ReferenceChromosomeName, ReferencePosition) = (record.chromosome_1.clone(), record.position_1.saturating_add(buffer));
        let mut origins: Vec<&str> = Vec::new();
        for dna_variant_record in index.range(lo..=hi).flat_map(|(_, dna_variant_records)| dna_variant_records.iter()) {
            if dna_variant_record.chromosome_2 == record.chromosome_2
                && dna_variant_record.position_2.abs_diff(record.position_2) <= buffer
                && dna_variant_record.operation_1 == record.operation_1
                && dna_variant_record.operation_2 == record.operation_2
                && dna_variant_record.variant_type == record.variant_type
                && dna_variant_record.sequence.eq_ignore_ascii_case(&record.sequence)
                && !origins.contains(&&*dna_variant_record.origin)
            {
                origins.push(&dna_variant_record.origin);
            }
        }
        record.origin = origins.join(LIST_SEPARATOR).into();
        record
    })
}


pub fn build_dna_variant_records<'a>(
    dvcs: &'a DNAVariantCallSet
) -> impl Iterator<Item = DNAVariantRecord> + 'a {
    // The maps are read only for a call, and a file without reads has no call and no read names.
    assert!(
        dvcs.variant_calls.is_empty() || !dvcs.chromosome_names_map.is_empty(),
        "dvcs.chromosome_names_map is empty."
    );
    assert!(
        dvcs.variant_calls.is_empty() || !dvcs.read_names_map.is_empty(),
        "dvcs.read_names_map is empty."
    );

    dvcs.get_variant_calls().into_iter().map(move |variant_call| {
        let variant_id: VariantID = variant_call.get_id();
        let consensus_graph_operation: &GraphOperation = variant_call.get_consensus_graph_operation();
        let consensus_read_names: Vec<ReadName> = variant_call.get_consensus_read_names(&dvcs.read_names_map);
        let chromosome_1: &str = dvcs
            .chromosome_names_map
            .get_by_right(&consensus_graph_operation.get_chromosome_1())
            .unwrap();
        let chromosome_2: &str = dvcs
            .chromosome_names_map
            .get_by_right(&consensus_graph_operation.get_chromosome_2())
            .unwrap();
        let read_names: Vec<&str> = variant_call
            .get_read_ids()
            .iter()
            .map(|read_id| dvcs.read_names_map.get_by_right(read_id).unwrap().as_ref())
            .collect();

        let num_consensus_read_names: u32 = consensus_read_names.len() as u32;
        let num_read_names: ReadSupport = read_names.len() as ReadSupport;

        DNAVariantRecord {
            origin: dvcs.origin.as_str().into(),
            variant_id: variant_id as u32,
            chromosome_1: chromosome_1.into(),
            position_1: consensus_graph_operation.get_position_1(),
            strand_1: consensus_graph_operation.get_strand_1().as_str().into(),
            operation_1: consensus_graph_operation.get_operation_type_1().as_str().into(),
            chromosome_2: chromosome_2.into(),
            position_2: consensus_graph_operation.get_position_2(),
            strand_2: consensus_graph_operation.get_strand_2().as_str().into(),
            operation_2: consensus_graph_operation.get_operation_type_2().as_str().into(),
            sequence: consensus_graph_operation.get_standardized_sequence().into(),
            variant_size: i32::try_from(consensus_graph_operation.get_variant_size()).ok().filter(|size| *size >= 0),
            variant_type: consensus_graph_operation.get_variant_type().as_str().into(),
            consensus_read_names: join_ids(consensus_read_names).into(),
            num_consensus_read_names: num_consensus_read_names,
            read_names: join_ids(read_names).into(),
            num_read_names: num_read_names
        }
    })
}


#[cfg(test)]
#[path = "../tests/io/builders.rs"]
mod tests;