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


use polars::prelude::*;
use std::collections::HashSet;

use crate::io::records::*;


fn join_u32_set(ids: &HashSet<u32>) -> String {
    let mut sorted: Vec<u32> = ids.iter().copied().collect();
    sorted.sort_unstable();
    sorted.iter().map(|id| id.to_string()).collect::<Vec<String>>().join(",")
}

pub fn primary_structure_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = PrimaryStructureRecord>
{
    let mut primary_structure_id: Vec<u64> = Vec::new();
    let mut amino_acid_sequence: Vec<String> = Vec::new();
    let mut amino_acid_sequence_length: Vec<u64> = Vec::new();
    let mut num_mutant_amino_acids: Vec<u32> = Vec::new();
    let mut assembled_transcript_name: Vec<String> = Vec::new();
    let mut assembled_transcript_sequence: Vec<String> = Vec::new();
    let mut assembled_transcript_sequence_length: Vec<u64> = Vec::new();
    let mut orf_start: Vec<u32> = Vec::new();
    let mut orf_end: Vec<u32> = Vec::new();
    let mut transcript_model_id: Vec<u32> = Vec::new();
    let mut reference_gene_names: Vec<String> = Vec::new();
    let mut reference_transcript_ids: Vec<String> = Vec::new();
    let mut mutant_amino_acid_intervals: Vec<String> = Vec::new();
    let mut rna_variant_ids: Vec<String> = Vec::new();
    let mut rna_variants: Vec<String> = Vec::new();
    let mut dna_variant_ids: Vec<String> = Vec::new();
    let mut dna_variants: Vec<String> = Vec::new();
    let mut rna_read_names: Vec<String> = Vec::new();
    let mut num_rna_read_names: Vec<u64> = Vec::new();
    let mut dna_variant_read_names: Vec<String> = Vec::new();
    let mut transcript_read_position: Vec<String> = Vec::new();

    for r in records {
        primary_structure_id.push(r.primary_structure_id as u64);
        amino_acid_sequence.push(r.amino_acid_sequence);
        amino_acid_sequence_length.push(r.amino_acid_sequence_length as u64);
        num_mutant_amino_acids.push(r.num_mutant_amino_acids);
        assembled_transcript_name.push(r.assembled_transcript_name.to_string());
        assembled_transcript_sequence.push(r.assembled_transcript_sequence.to_string());
        assembled_transcript_sequence_length.push(r.assembled_transcript_sequence_length as u64);
        orf_start.push(r.orf_start);
        orf_end.push(r.orf_end);
        transcript_model_id.push(r.transcript_model_id);
        reference_gene_names.push(r.reference_gene_names);
        reference_transcript_ids.push(r.reference_transcript_ids);
        mutant_amino_acid_intervals.push(r.mutant_amino_acid_intervals);
        rna_variant_ids.push(r.rna_variant_ids);
        rna_variants.push(r.rna_variants);
        dna_variant_ids.push(r.dna_variant_ids);
        dna_variants.push(r.dna_variants);
        rna_read_names.push(r.rna_read_names.to_string());
        num_rna_read_names.push(r.num_rna_read_names as u64);
        dna_variant_read_names.push(r.dna_variant_read_names.to_string());
        transcript_read_position.push(r.transcript_read_position.to_string());
    }

    DataFrame::new(vec![
        Column::from(Series::new("primary_structure_id".into(), primary_structure_id)),
        Column::from(Series::new("amino_acid_sequence".into(), amino_acid_sequence)),
        Column::from(Series::new("amino_acid_sequence_length".into(), amino_acid_sequence_length)),
        Column::from(Series::new("num_mutant_amino_acids".into(), num_mutant_amino_acids)),
        Column::from(Series::new("assembled_transcript_name".into(), assembled_transcript_name)),
        Column::from(Series::new("assembled_transcript_sequence".into(), assembled_transcript_sequence)),
        Column::from(Series::new("assembled_transcript_sequence_length".into(), assembled_transcript_sequence_length)),
        Column::from(Series::new("transcript_model_id".into(), transcript_model_id)),
        Column::from(Series::new("orf_start".into(), orf_start)),
        Column::from(Series::new("orf_end".into(), orf_end)),
        Column::from(Series::new("reference_gene_names".into(), reference_gene_names)),
        Column::from(Series::new("reference_transcript_ids".into(), reference_transcript_ids)),
        Column::from(Series::new("mutant_amino_acid_intervals".into(), mutant_amino_acid_intervals)),
        Column::from(Series::new("rna_variant_ids".into(), rna_variant_ids)),
        Column::from(Series::new("rna_variants".into(), rna_variants)),
        Column::from(Series::new("dna_variant_ids".into(), dna_variant_ids)),
        Column::from(Series::new("dna_variants".into(), dna_variants)),
        Column::from(Series::new("rna_read_names".into(), rna_read_names)),
        Column::from(Series::new("num_rna_read_names".into(), num_rna_read_names)),
        Column::from(Series::new("dna_variant_read_names".into(), dna_variant_read_names)),
        Column::from(Series::new("transcript_read_position".into(), transcript_read_position)),
    ])
    .unwrap()
}


/// Convert an iterator of `NucleotideRecord` into a polars `DataFrame`
/// whose columns match the record's serde field order (same names, same
/// order as the csv writer emits). Nullable fields use polars' native
/// nullability via `Vec<Option<T>>`; `HashSet<u32>` columns are flattened
/// to comma-joined strings with deterministic (sorted) ordering.
pub fn nucleotide_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = NucleotideRecord>
{
    let mut primary_structure_id: Vec<u32> = Vec::new();
    let mut assembled_transcript_name: Vec<String> = Vec::new();
    let mut transcript_model_id: Vec<u32> = Vec::new();
    let mut amino_acid_index: Vec<u32> = Vec::new();
    let mut amino_acid: Vec<String> = Vec::new();
    // codon_index is u8 on the record (values 0/1/2). Widen to u32 here to
    // keep the index columns in this DataFrame all the same dtype.
    let mut codon_index: Vec<u32> = Vec::new();
    let mut nucleotide: Vec<String> = Vec::new();
    let mut is_amino_acid_variant: Vec<bool> = Vec::new();
    let mut is_nucleotide_variant: Vec<bool> = Vec::new();
    let mut transcript_read_position: Vec<u32> = Vec::new();
    let mut transcript_structure_index: Vec<Option<u32>> = Vec::new();
    let mut rna_variant_id: Vec<Option<u32>> = Vec::new();
    let mut rna_variant: Vec<Option<String>> = Vec::new();
    let mut dna_variant_ids: Vec<Option<String>> = Vec::new();
    let mut dna_variant: Vec<Option<String>> = Vec::new();
    let mut preceding_event_rna_variant_id: Vec<Option<u32>> = Vec::new();
    let mut preceding_event_rna_variant: Vec<Option<String>> = Vec::new();
    let mut preceding_event_dna_variant_ids: Vec<Option<String>> = Vec::new();
    let mut preceding_event_dna_variant: Vec<Option<String>> = Vec::new();

    for r in records {
        primary_structure_id.push(r.primary_structure_id);
        assembled_transcript_name.push(r.assembled_transcript_name.to_string());
        transcript_model_id.push(r.transcript_model_id);
        amino_acid_index.push(r.amino_acid_index);
        amino_acid.push(r.amino_acid.to_string());
        codon_index.push(r.codon_index as u32);
        nucleotide.push(r.nucleotide.to_string());
        is_amino_acid_variant.push(r.is_amino_acid_variant);
        is_nucleotide_variant.push(r.is_nucleotide_variant);
        transcript_read_position.push(r.transcript_read_position);
        transcript_structure_index.push(r.transcript_structure_index);
        rna_variant_id.push(r.rna_variant_id);
        rna_variant.push(r.rna_variant.map(|s| s.to_string()));
        dna_variant_ids.push(r.dna_variant_ids.as_ref().map(join_u32_set));
        dna_variant.push(r.dna_variant.map(|s| s.to_string()));
        preceding_event_rna_variant_id.push(r.preceding_event_rna_variant_id);
        preceding_event_rna_variant.push(r.preceding_event_rna_variant.map(|s| s.to_string()));
        preceding_event_dna_variant_ids.push(r.preceding_event_dna_variant_ids.as_ref().map(join_u32_set));
        preceding_event_dna_variant.push(r.preceding_event_dna_variant.map(|s| s.to_string()));
    }

    DataFrame::new(vec![
        Column::from(Series::new("primary_structure_id".into(), primary_structure_id)),
        Column::from(Series::new("assembled_transcript_name".into(), assembled_transcript_name)),
        Column::from(Series::new("transcript_model_id".into(), transcript_model_id)),
        Column::from(Series::new("amino_acid_index".into(), amino_acid_index)),
        Column::from(Series::new("amino_acid".into(), amino_acid)),
        Column::from(Series::new("codon_index".into(), codon_index)),
        Column::from(Series::new("nucleotide".into(), nucleotide)),
        Column::from(Series::new("is_amino_acid_variant".into(), is_amino_acid_variant)),
        Column::from(Series::new("is_nucleotide_variant".into(), is_nucleotide_variant)),
        Column::from(Series::new("transcript_read_position".into(), transcript_read_position)),
        Column::from(Series::new("transcript_structure_index".into(), transcript_structure_index)),
        Column::from(Series::new("rna_variant_id".into(), rna_variant_id)),
        Column::from(Series::new("rna_variant".into(), rna_variant)),
        Column::from(Series::new("dna_variant_ids".into(), dna_variant_ids)),
        Column::from(Series::new("dna_variant".into(), dna_variant)),
        Column::from(Series::new("preceding_event_rna_variant_id".into(), preceding_event_rna_variant_id)),
        Column::from(Series::new("preceding_event_rna_variant".into(), preceding_event_rna_variant)),
        Column::from(Series::new("preceding_event_dna_variant_ids".into(), preceding_event_dna_variant_ids)),
        Column::from(Series::new("preceding_event_dna_variant".into(), preceding_event_dna_variant)),
    ])
    .unwrap()
}
