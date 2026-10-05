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


use exacto_core::prelude::{ReadPosition, ReadSupport, ReferencePosition};
use polars::prelude::*;

use crate::io::records::*;


pub fn assembled_transcript_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = AssembledTranscriptRecord>
{
    let mut assembled_transcript_name: Vec<String> = Vec::new();
    let mut start_chromosome: Vec<String> = Vec::new();
    let mut start: Vec<ReferencePosition> = Vec::new();
    let mut end_chromosome: Vec<String> = Vec::new();
    let mut end: Vec<ReferencePosition> = Vec::new();
    let mut num_exons: Vec<u32> = Vec::new();
    let mut num_splice_junctions: Vec<u32> = Vec::new();
    let mut is_variant: Vec<bool> = Vec::new();

    for r in records {
        assembled_transcript_name.push(r.assembled_transcript_name.to_string());
        start_chromosome.push(r.start_chromosome.to_string());
        start.push(r.start);
        end_chromosome.push(r.end_chromosome.to_string());
        end.push(r.end);
        num_exons.push(r.num_exons);
        num_splice_junctions.push(r.num_splice_junctions);
        is_variant.push(r.is_variant);
    }

    DataFrame::new(vec![
        Column::from(Series::new("assembled_transcript_name".into(), assembled_transcript_name)),
        Column::from(Series::new("start_chromosome".into(), start_chromosome)),
        Column::from(Series::new("start".into(), start)),
        Column::from(Series::new("end_chromosome".into(), end_chromosome)),
        Column::from(Series::new("end".into(), end)),
        Column::from(Series::new("num_exons".into(), num_exons)),
        Column::from(Series::new("num_splice_junctions".into(), num_splice_junctions)),
        Column::from(Series::new("is_variant".into(), is_variant)),
    ])
    .unwrap()
}


pub fn assembled_transcript_exon_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item =AssembledTranscriptExonRecord>
{
    let mut assembled_transcript_name: Vec<String> = Vec::new();
    let mut chromosome: Vec<String> = Vec::new();
    let mut start: Vec<ReferencePosition> = Vec::new();
    let mut end: Vec<ReferencePosition> = Vec::new();
    let mut exon_number: Vec<u32> = Vec::new();
    let mut strand: Vec<String> = Vec::new();
    for r in records {
        assembled_transcript_name.push(r.assembled_transcript_name.to_string());
        chromosome.push(r.chromosome.to_string());
        start.push(r.start);
        end.push(r.end);
        exon_number.push(r.exon_number);
        strand.push(r.strand.to_string());
    }

    DataFrame::new(vec![
        Column::from(Series::new("assembled_transcript_name".into(), assembled_transcript_name)),
        Column::from(Series::new("chromosome".into(), chromosome)),
        Column::from(Series::new("start".into(), start)),
        Column::from(Series::new("end".into(), end)),
        Column::from(Series::new("exon_number".into(), exon_number)),
        Column::from(Series::new("strand".into(), strand)),
    ])
    .unwrap()
}


pub fn assembled_transcript_splice_junction_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item =AssembledTranscriptSpliceJunctionRecord>
{
    let mut assembled_transcript_name: Vec<String> = Vec::new();
    let mut chromosome_1: Vec<String> = Vec::new();
    let mut chromosome_2: Vec<String> = Vec::new();
    let mut position_1: Vec<ReferencePosition> = Vec::new();
    let mut position_2: Vec<ReferencePosition> = Vec::new();
    let mut strand_1: Vec<String> = Vec::new();
    let mut strand_2: Vec<String> = Vec::new();
    let mut splice_junction_number: Vec<u32> = Vec::new();
    for r in records {
        assembled_transcript_name.push(r.assembled_transcript_name.to_string());
        chromosome_1.push(r.chromosome_1.to_string());
        chromosome_2.push(r.chromosome_2.to_string());
        position_1.push(r.position_1);
        position_2.push(r.position_2);
        strand_1.push(r.strand_1.to_string());
        strand_2.push(r.strand_2.to_string());
        splice_junction_number.push(r.splice_junction_number);
    }

    DataFrame::new(vec![
        Column::from(Series::new("assembled_transcript_name".into(), assembled_transcript_name)),
        Column::from(Series::new("chromosome_1".into(), chromosome_1)),
        Column::from(Series::new("chromosome_2".into(), chromosome_2)),
        Column::from(Series::new("position_1".into(), position_1)),
        Column::from(Series::new("position_2".into(), position_2)),
        Column::from(Series::new("strand_1".into(), strand_1)),
        Column::from(Series::new("strand_2".into(), strand_2)),
        Column::from(Series::new("splice_junction_number".into(), splice_junction_number))
    ])
    .unwrap()
}


pub fn assembled_transcript_reference_transcript_match_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item =AssembledTranscriptReferenceTranscriptMatchRecord>
{
    let mut assembled_transcript_name: Vec<String> = Vec::new();
    let mut reference_gene_id: Vec<String> = Vec::new();
    let mut reference_gene_name: Vec<String> = Vec::new();
    let mut reference_transcript_id: Vec<String> = Vec::new();
    let mut num_splice_junction_matches: Vec<u32> = Vec::new();
    let mut num_overlapping_bases: Vec<u32> = Vec::new();
    let mut num_reference_only_bases: Vec<u32> = Vec::new();
    let mut num_query_only_bases: Vec<u32> = Vec::new();
    let mut num_terminal_query_only_bases: Vec<u32> = Vec::new();
    let mut num_internal_query_only_bases: Vec<u32> = Vec::new();
    for r in records {
        assembled_transcript_name.push(r.assembled_transcript_name.to_string());
        reference_gene_id.push(r.reference_gene_id.to_string());
        reference_gene_name.push(r.reference_gene_name.to_string());
        reference_transcript_id.push(r.reference_transcript_id.to_string());
        num_splice_junction_matches.push(r.num_splice_junction_matches);
        num_overlapping_bases.push(r.num_overlapping_bases);
        num_reference_only_bases.push(r.num_reference_only_bases);
        num_query_only_bases.push(r.num_query_only_bases);
        num_terminal_query_only_bases.push(r.num_terminal_query_only_bases);
        num_internal_query_only_bases.push(r.num_internal_query_only_bases);
    }

    DataFrame::new(vec![
        Column::from(Series::new("assembled_transcript_name".into(), assembled_transcript_name)),
        Column::from(Series::new("reference_gene_id".into(), reference_gene_id)),
        Column::from(Series::new("reference_gene_name".into(), reference_gene_name)),
        Column::from(Series::new("reference_transcript_id".into(), reference_transcript_id)),
        Column::from(Series::new("num_splice_junction_matches".into(), num_splice_junction_matches)),
        Column::from(Series::new("num_overlapping_bases".into(), num_overlapping_bases)),
        Column::from(Series::new("num_reference_only_bases".into(), num_reference_only_bases)),
        Column::from(Series::new("num_query_only_bases".into(), num_query_only_bases)),
        Column::from(Series::new("num_terminal_query_only_bases".into(), num_terminal_query_only_bases)),
        Column::from(Series::new("num_internal_query_only_bases".into(), num_internal_query_only_bases))
    ])
    .unwrap()
}


pub fn assembled_transcript_nonsense_mediated_decay_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = AssembledTranscriptNonsenseMediatedDecayRecord>
{
    let mut assembled_transcript_name: Vec<String> = Vec::new();
    let mut orf_start: Vec<ReadPosition> = Vec::new();
    let mut orf_end: Vec<ReadPosition> = Vec::new();
    let mut nmd_predicted: Vec<bool> = Vec::new();
    let mut distance_to_last_junction: Vec<Option<u32>> = Vec::new();
    for r in records {
        assembled_transcript_name.push(r.assembled_transcript_name.to_string());
        orf_start.push(r.orf_start);
        orf_end.push(r.orf_end);
        nmd_predicted.push(r.nmd_predicted);
        distance_to_last_junction.push(r.distance_to_last_junction);
    }

    DataFrame::new(vec![
        Column::from(Series::new("assembled_transcript_name".into(), assembled_transcript_name)),
        Column::from(Series::new("orf_start".into(), orf_start)),
        Column::from(Series::new("orf_end".into(), orf_end)),
        Column::from(Series::new("nmd_predicted".into(), nmd_predicted)),
        Column::from(Series::new("distance_to_last_junction".into(), distance_to_last_junction))
    ])
    .unwrap()
}


pub fn assembled_transcript_filter_status_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item =AssembledTranscriptFilterStatusRecord>
{
    let mut read_name: Vec<String> = Vec::new();
    let mut excluded: Vec<bool> = Vec::new();
    for r in records {
        read_name.push(r.read_name.to_string());
        excluded.push(r.excluded);
    }

    DataFrame::new(vec![
        Column::from(Series::new("read_name".into(), read_name)),
        Column::from(Series::new("excluded".into(), excluded))
    ])
    .unwrap()
}

pub fn assembled_transcript_model_alignment_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item =AssembledTranscriptModelAlignmentRecord>
{
    let mut assembled_transcript_name: Vec<String> = Vec::new();
    let mut reference_gene_name: Vec<String> = Vec::new();
    let mut reference_transcript_id: Vec<String> = Vec::new();
    let mut index: Vec<u32> = Vec::new();
    let mut read_start: Vec<ReadPosition> = Vec::new();
    let mut read_end: Vec<ReadPosition> = Vec::new();
    let mut sequence: Vec<String> = Vec::new();
    let mut record_type: Vec<String> = Vec::new();
    let mut kind: Vec<String> = Vec::new();
    let mut context: Vec<String> = Vec::new();
    let mut chromosome_1: Vec<String> = Vec::new();
    let mut position_1: Vec<ReferencePosition> = Vec::new();
    let mut operation_1: Vec<String> = Vec::new();
    let mut strand_1: Vec<String> = Vec::new();
    let mut chromosome_2: Vec<String> = Vec::new();
    let mut position_2: Vec<ReferencePosition> = Vec::new();
    let mut operation_2: Vec<String> = Vec::new();
    let mut strand_2: Vec<String> = Vec::new();
    let mut gene_id_1: Vec<String> = Vec::new();
    let mut transcript_id_1: Vec<String> = Vec::new();
    let mut exon_id_1: Vec<String> = Vec::new();
    let mut gene_id_2: Vec<String> = Vec::new();
    let mut transcript_id_2: Vec<String> = Vec::new();
    let mut exon_id_2: Vec<String> = Vec::new();
    let mut skipped: Vec<String> = Vec::new();
    for r in records {
        assembled_transcript_name.push(r.assembled_transcript_name.to_string());
        reference_gene_name.push(r.reference_gene_name.to_string());
        reference_transcript_id.push(r.reference_transcript_id.to_string());
        index.push(r.index);
        read_start.push(r.read_start);
        read_end.push(r.read_end);
        sequence.push(r.sequence.to_string());
        record_type.push(r.record_type.to_string());
        kind.push(r.kind.to_string());
        context.push(r.context.to_string());
        chromosome_1.push(r.chromosome_1.to_string());
        position_1.push(r.position_1);
        operation_1.push(r.operation_1.to_string());
        strand_1.push(r.strand_1.to_string());
        chromosome_2.push(r.chromosome_2.to_string());
        position_2.push(r.position_2);
        operation_2.push(r.operation_2.to_string());
        strand_2.push(r.strand_2.to_string());
        gene_id_1.push(r.reference_gene_id_1.to_string());
        transcript_id_1.push(r.reference_transcript_id_1.to_string());
        exon_id_1.push(r.reference_exon_id_1.to_string());
        gene_id_2.push(r.reference_gene_id_2.to_string());
        transcript_id_2.push(r.reference_transcript_id_2.to_string());
        exon_id_2.push(r.reference_exon_id_2.to_string());
        skipped.push(r.skipped.to_string());
    }

    DataFrame::new(vec![
        Column::from(Series::new("assembled_transcript_name".into(), assembled_transcript_name)),
        Column::from(Series::new("reference_gene_name".into(), reference_gene_name)),
        Column::from(Series::new("reference_transcript_id".into(), reference_transcript_id)),
        Column::from(Series::new("index".into(), index)),
        Column::from(Series::new("read_start".into(), read_start)),
        Column::from(Series::new("read_end".into(), read_end)),
        Column::from(Series::new("sequence".into(), sequence)),
        Column::from(Series::new("type".into(), record_type)),
        Column::from(Series::new("kind".into(), kind)),
        Column::from(Series::new("context".into(), context)),
        Column::from(Series::new("chromosome_1".into(), chromosome_1)),
        Column::from(Series::new("position_1".into(), position_1)),
        Column::from(Series::new("operation_1".into(), operation_1)),
        Column::from(Series::new("strand_1".into(), strand_1)),
        Column::from(Series::new("chromosome_2".into(), chromosome_2)),
        Column::from(Series::new("position_2".into(), position_2)),
        Column::from(Series::new("operation_2".into(), operation_2)),
        Column::from(Series::new("strand_2".into(), strand_2)),
        Column::from(Series::new("reference_gene_id_1".into(), gene_id_1)),
        Column::from(Series::new("reference_transcript_id_1".into(), transcript_id_1)),
        Column::from(Series::new("reference_exon_id_1".into(), exon_id_1)),
        Column::from(Series::new("reference_gene_id_2".into(), gene_id_2)),
        Column::from(Series::new("reference_transcript_id_2".into(), transcript_id_2)),
        Column::from(Series::new("reference_exon_id_2".into(), exon_id_2)),
        Column::from(Series::new("skipped".into(), skipped))
    ])
    .unwrap()
}


pub fn assembled_transcript_variant_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item =AssembledTranscriptVariantRecord>
{
    let mut variant_id: Vec<u32> = Vec::new();
    let mut assembled_transcript_name: Vec<String> = Vec::new();
    let mut reference_gene_name: Vec<String> = Vec::new();
    let mut reference_transcript_id: Vec<String> = Vec::new();
    let mut chromosome_1: Vec<String> = Vec::new();
    let mut position_1: Vec<ReferencePosition> = Vec::new();
    let mut strand_1: Vec<String> = Vec::new();
    let mut operation_1: Vec<String> = Vec::new();
    let mut chromosome_2: Vec<String> = Vec::new();
    let mut position_2: Vec<ReferencePosition> = Vec::new();
    let mut strand_2: Vec<String> = Vec::new();
    let mut operation_2: Vec<String> = Vec::new();
    let mut variant_size: Vec<Option<u32>> = Vec::new();
    let mut variant_type: Vec<String> = Vec::new();
    let mut sequence: Vec<String> = Vec::new();
    let mut read_start: Vec<ReadPosition> = Vec::new();
    let mut read_end: Vec<ReadPosition> = Vec::new();
    let mut origin: Vec<String> = Vec::new();
    for r in records {
        variant_id.push(r.variant_id);
        assembled_transcript_name.push(r.assembled_transcript_name.to_string());
        reference_gene_name.push(r.reference_gene_name.to_string());
        reference_transcript_id.push(r.reference_transcript_id.to_string());
        chromosome_1.push(r.chromosome_1.to_string());
        position_1.push(r.position_1);
        strand_1.push(r.strand_1.to_string());
        operation_1.push(r.operation_1.to_string());
        chromosome_2.push(r.chromosome_2.to_string());
        position_2.push(r.position_2);
        strand_2.push(r.strand_2.to_string());
        operation_2.push(r.operation_2.to_string());
        variant_size.push(r.variant_size);
        variant_type.push(r.variant_type.to_string());
        sequence.push(r.sequence.to_string());
        read_start.push(r.read_start);
        read_end.push(r.read_end);
        origin.push(r.origin.to_string());
    }

    DataFrame::new(vec![
        Column::from(Series::new("variant_id".into(), variant_id)),
        Column::from(Series::new("assembled_transcript_name".into(), assembled_transcript_name)),
        Column::from(Series::new("reference_gene_name".into(), reference_gene_name)),
        Column::from(Series::new("reference_transcript_id".into(), reference_transcript_id)),
        Column::from(Series::new("chromosome_1".into(), chromosome_1)),
        Column::from(Series::new("position_1".into(), position_1)),
        Column::from(Series::new("strand_1".into(), strand_1)),
        Column::from(Series::new("operation_1".into(), operation_1)),
        Column::from(Series::new("chromosome_2".into(), chromosome_2)),
        Column::from(Series::new("position_2".into(), position_2)),
        Column::from(Series::new("strand_2".into(), strand_2)),
        Column::from(Series::new("operation_2".into(), operation_2)),
        Column::from(Series::new("variant_size".into(), variant_size)),
        Column::from(Series::new("variant_type".into(), variant_type)),
        Column::from(Series::new("sequence".into(), sequence)),
        Column::from(Series::new("read_start".into(), read_start)),
        Column::from(Series::new("read_end".into(), read_end)),
        Column::from(Series::new("origin".into(), origin))
    ])
    .unwrap()
}


pub fn dna_variant_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = DNAVariantRecord>
{
    let mut variant_id: Vec<u32> = Vec::new();
    let mut origin: Vec<String> = Vec::new();
    let mut chromosome_1: Vec<String> = Vec::new();
    let mut position_1: Vec<ReferencePosition> = Vec::new();
    let mut strand_1: Vec<String> = Vec::new();
    let mut operation_1: Vec<String> = Vec::new();
    let mut chromosome_2: Vec<String> = Vec::new();
    let mut position_2: Vec<ReferencePosition> = Vec::new();
    let mut strand_2: Vec<String> = Vec::new();
    let mut operation_2: Vec<String> = Vec::new();
    let mut sequence: Vec<String> = Vec::new();
    let mut variant_size: Vec<Option<i32>> = Vec::new();
    let mut variant_type: Vec<String> = Vec::new();
    let mut consensus_read_names: Vec<String> = Vec::new();
    let mut num_consensus_read_names: Vec<u32> = Vec::new();
    let mut read_names: Vec<String> = Vec::new();
    let mut num_read_names: Vec<ReadSupport> = Vec::new();

    for r in records {
        variant_id.push(r.variant_id);
        origin.push(r.origin.to_string());
        chromosome_1.push(r.chromosome_1.to_string());
        position_1.push(r.position_1);
        strand_1.push(r.strand_1.to_string());
        operation_1.push(r.operation_1.to_string());
        chromosome_2.push(r.chromosome_2.to_string());
        position_2.push(r.position_2);
        strand_2.push(r.strand_2.to_string());
        operation_2.push(r.operation_2.to_string());
        sequence.push(r.sequence.to_string());
        variant_size.push(r.variant_size);
        variant_type.push(r.variant_type.to_string());
        consensus_read_names.push(r.consensus_read_names.to_string());
        num_consensus_read_names.push(r.num_consensus_read_names);
        read_names.push(r.read_names.to_string());
        num_read_names.push(r.num_read_names);
    }

    DataFrame::new(vec![
        Column::from(Series::new("variant_id".into(), variant_id)),
        Column::from(Series::new("origin".into(), origin)),
        Column::from(Series::new("chromosome_1".into(), chromosome_1)),
        Column::from(Series::new("position_1".into(), position_1)),
        Column::from(Series::new("strand_1".into(), strand_1)),
        Column::from(Series::new("operation_1".into(), operation_1)),
        Column::from(Series::new("chromosome_2".into(), chromosome_2)),
        Column::from(Series::new("position_2".into(), position_2)),
        Column::from(Series::new("strand_2".into(), strand_2)),
        Column::from(Series::new("operation_2".into(), operation_2)),
        Column::from(Series::new("sequence".into(), sequence)),
        Column::from(Series::new("variant_size".into(), variant_size)),
        Column::from(Series::new("variant_type".into(), variant_type)),
        Column::from(Series::new("consensus_read_names".into(), consensus_read_names)),
        Column::from(Series::new("num_consensus_read_names".into(), num_consensus_read_names)),
        Column::from(Series::new("read_names".into(), read_names)),
        Column::from(Series::new("num_read_names".into(), num_read_names))
    ])
        .unwrap()
}


#[cfg(test)]
#[path = "../tests/io/dataframes.rs"]
mod tests;