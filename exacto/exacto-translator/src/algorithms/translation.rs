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
use polars::prelude::*;
use rayon::prelude::*;
use std::str::FromStr;

use crate::prelude::*;


pub fn translate_rnas(
    rnas: Vec<RNA>,
    num_threads: usize
) -> TranslationSet {
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let translations: Vec<Translation> = thread_pool.install(|| {
        rnas
            .par_iter()
            .filter_map(|rna| rna.translate()) // Filters out `None` and unwraps `Some`.
            .collect()
    });
    let mut translation_set: TranslationSet = TranslationSet::new();
    for translation in translations {
        translation_set.add_translation(translation);
    }
    translation_set
}

pub fn translate_transcript_structures(
    df_transcript_structures: &DataFrame,
    rna_variant_call_set: &RNAVariantCallSet,
    df_integrated_variants: &DataFrame,
    translation_strategy: TranslationStrategy,
    num_threads: usize
) -> PrimaryStructureSet {
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    let partitions: Vec<DataFrame> = df_transcript_structures
        .partition_by(["transcript_model_id", "reference_transcript_ids"], true)
        .unwrap();

    let primary_structures: Vec<PrimaryStructure> = thread_pool.install(|| {
        partitions
            .into_par_iter()
            .map(|df| {
                // Step 1. Sort the DataFrame
                let df_sorted: DataFrame = df.sort(["transcript_structure_index"], SortMultipleOptions::default()).unwrap();

                // Step 2. Fetch the RNA sequence
                let col_sequence = df.column("sequence").unwrap().str().unwrap();
                let mut sequence: String = String::new();
                for i in 0..df_sorted.height() {
                    sequence.push_str(col_sequence.get(i).unwrap());
                }

                // Step 3. Translate the RNA sequence
                let rna: RNA = RNA::new("".into(), sequence.into());
                let translation_result: Option<Translation> = rna.translate();
                if translation_result.is_none() {
                    return Vec::new();
                }
                let translation: Translation = translation_result.unwrap();
                let peptides: Vec<Peptide> = match translation_strategy {
                    TranslationStrategy::LongestORF => {
                        vec![translation.get_longest_orf_peptide().clone()]
                    },
                    TranslationStrategy::AllORFs => {
                        translation.get_peptides().clone()
                    }
                };

                // Step 4. Fetch the DataFrame columns
                let transcript_model_id: usize = df.column("transcript_model_id").unwrap().i64().unwrap().get(0).unwrap() as usize;
                let reference_transcript_ids: Vec<Box<str>> = df.column("reference_transcript_ids").unwrap().str().unwrap().get(0).unwrap().split(",").map(|s| s.into()).collect();
                let reference_transcript_ids_str: Box<str> = df.column("reference_transcript_ids").unwrap().str().unwrap().get(0).unwrap().into();
                let col_type = df.column("type").unwrap().str().unwrap();
                let col_read_start = df.column("read_start").unwrap().i64().unwrap();
                let col_read_end = df.column("read_end").unwrap().i64().unwrap();
                let col_sequence = df.column("sequence").unwrap().str().unwrap();
                let col_tsi = df.column("transcript_structure_index").unwrap().i64().unwrap();
                let col_kind = df.column("kind").unwrap().str().unwrap();
                let col_context = df.column("context").unwrap().str().unwrap();
                let col_skipped = df.column("skipped").unwrap().str().unwrap();

                // Step 5. Fetch RNA variant calls for the current transcript model ID and reference transcript IDs
                let mut rna_variant_calls: Vec<&VariantCall> = Vec::new();
                if rna_variant_call_set.variant_calls_index.contains_key(&(transcript_model_id)) {
                    if rna_variant_call_set.variant_calls_index.get(&(transcript_model_id)).unwrap().contains_key(&reference_transcript_ids) {
                        for variant_call_id in rna_variant_call_set.variant_calls_index.get(&(transcript_model_id)).unwrap().get(&reference_transcript_ids).unwrap() {
                            rna_variant_calls.push(rna_variant_call_set.get_variant_call(*variant_call_id));
                        }
                    }
                }

                // Step 6. Fetch integrated variant rows for the current transcript model ID and reference transcript IDs
                let df_integrated_variants_result: PolarsResult<DataFrame> = df_integrated_variants
                    .clone()
                    .lazy()
                    .filter(col("transcript_model_id").eq(lit(transcript_model_id as u64)))
                    .filter(col("reference_transcript_ids").eq(lit(reference_transcript_ids_str.to_string())))
                    .collect();
                let df_curr_integrated_variants: DataFrame = df_integrated_variants_result.unwrap();

                // Step 7. Identify the underlying primary structure of each peptide
                let mut primary_structures: Vec<PrimaryStructure> = Vec::new();
                for (i, peptide) in peptides.iter().enumerate() {
                    let mut primary_structure: PrimaryStructure = PrimaryStructure::new();
                    let mut primary_structure_index: usize = 0;
                    let mut codon_index: u8 = 0;
                    let mut net_variant_nucleotides_count: i32 = 0;
                    for j in 0..df.height() {
                        let record_type: AlignmentStructureRecordType = AlignmentStructureRecordType::from_str(col_type.get(j).unwrap()).unwrap();
                        let read_start: usize = col_read_start.get(j).unwrap() as usize;
                        let read_end: usize = col_read_end.get(j).unwrap() as usize;
                        let sequence: &str = col_sequence.get(j).unwrap();
                        let transcript_structure_index: usize = col_tsi.get(j).unwrap() as usize;
                        let kind_str: &str = col_kind.get(j).unwrap();
                        let context_str: &str = col_context.get(j).unwrap();
                        let skipped: &str = col_skipped.get(j).unwrap();

                        // Only start recording if the current row is within the open reading frame
                        if overlaps(
                            peptide.orf_start as isize,
                            peptide.orf_end as isize,
                            read_start as isize,
                            read_end as isize) == false {
                            continue;
                        }

                        // Get any RNA and DNA variants matching the current transcript structure row
                        let mut rna_variant_call_ids: Vec<usize> = Vec::new();
                        let mut dna_variant_call_ids: Vec<usize> = Vec::new();
                        for rna_variant_call in rna_variant_calls.iter() {
                            let rna_variant_record: &VariantRecord = rna_variant_call.get_consensus_record().0;
                            if rna_variant_record.get_read_position_1() == read_start && rna_variant_record.get_read_position_2() == read_end {
                                rna_variant_call_ids.push(rna_variant_call.id);
                                // Get DNA variant call IDs
                                let df_curr_integrated_variants_: PolarsResult<DataFrame> = df_curr_integrated_variants
                                    .clone()
                                    .lazy()
                                    .filter(col("rna_variant_call_id").eq(lit(rna_variant_call.id as u64)))
                                    .collect();
                                let df_curr_integrated_variants_matched: DataFrame = df_curr_integrated_variants_.unwrap();
                                if df_curr_integrated_variants_matched.height() > 0 {
                                    let col_dna_variant_call_id = df_curr_integrated_variants_matched.column("dna_variant_call_id").unwrap().i64().unwrap();
                                    for k in 0..df_curr_integrated_variants_matched.height() {
                                        let dna_variant_call_id: usize = col_dna_variant_call_id.get(k).unwrap() as usize;
                                        dna_variant_call_ids.push(dna_variant_call_id);
                                    }
                                }
                            }
                        }

                        // Create primary structure records
                        match record_type {
                            AlignmentStructureRecordType::Base => {
                                let kind: AlignmentStructureBaseKind = AlignmentStructureBaseKind::from_str(kind_str).unwrap();
                                let context: AlignmentStructureBaseContext = AlignmentStructureBaseContext::from_str(context_str).unwrap();

                                let mut read_position: usize = read_start;
                                let mut chars = sequence.char_indices().peekable();
                                while let Some((start, _)) = chars.next() {
                                    if read_position >= peptide.orf_start &&
                                        read_position <= peptide.orf_end {
                                        let end: usize = chars.peek().map(|(k, _)| *k).unwrap_or(sequence.len());
                                        let nucleotide:&str = &sequence[start..end];

                                        match (&kind, &context) {
                                            (AlignmentStructureBaseKind::Insertion, _) => {
                                                net_variant_nucleotides_count += 1;
                                            },
                                            (AlignmentStructureBaseKind::Match, AlignmentStructureBaseContext::Intronic) => {
                                                net_variant_nucleotides_count += 1;
                                            },
                                            (AlignmentStructureBaseKind::Match, AlignmentStructureBaseContext::Intergenic) => {
                                                net_variant_nucleotides_count += 1;
                                            },
                                            (AlignmentStructureBaseKind::Unaligned, _) => {
                                                net_variant_nucleotides_count += 1;
                                            },
                                            (_, _) => {
                                                // Do nothing
                                            }
                                        }

                                        primary_structure.add_record(
                                            PrimaryStructureRecord::new(
                                                primary_structure_index,
                                                PrimaryStructureRecordType::Base,
                                                None,
                                                Some(codon_index),
                                                Some(Nucleotide::from_str(nucleotide).unwrap()),
                                                transcript_model_id,
                                                reference_transcript_ids.clone(),
                                                transcript_structure_index,
                                                read_position,
                                                read_position,
                                                net_variant_nucleotides_count,
                                                None,
                                                rna_variant_call_ids.clone(),
                                                dna_variant_call_ids.clone()
                                            )
                                        );

                                        primary_structure_index += 1;
                                        codon_index += 1;

                                        if codon_index == 3 {
                                            codon_index = 0;
                                        }
                                    }

                                    read_position += 1;
                                }
                            },
                            AlignmentStructureRecordType::Event => {
                                if skipped.is_empty() == false {
                                    for skipped_coordinates in skipped.split(",") {
                                        let skipped_elements: Vec<&str> = skipped_coordinates.split("|").collect();
                                        let chromosome_1: &str = skipped_elements[0].split(":").collect::<Vec<&str>>()[0];
                                        let position_1: usize = skipped_elements[0].split(":").collect::<Vec<&str>>()[1].parse::<usize>().unwrap();
                                        let chromosome_2: &str = skipped_elements[1].split(":").collect::<Vec<&str>>()[0];
                                        let position_2: usize = skipped_elements[1].split(":").collect::<Vec<&str>>()[1].parse::<usize>().unwrap();
                                        let skipped_length: usize = position_1.abs_diff(position_2) + 1;
                                        net_variant_nucleotides_count -= skipped_length as i32;
                                    }
                                }

                                primary_structure.add_record(
                                    PrimaryStructureRecord::new(
                                        primary_structure_index,
                                        PrimaryStructureRecordType::Event,
                                        None,
                                        None,
                                        None,
                                        transcript_model_id,
                                        reference_transcript_ids.clone(),
                                        transcript_structure_index,
                                        read_start,
                                        read_end,
                                        net_variant_nucleotides_count,
                                        None,
                                        rna_variant_call_ids,
                                        dna_variant_call_ids
                                    )
                                );

                                primary_structure_index += 1;
                            }
                        }
                    }

                    // Assign the amino acid and frameshift states
                    let mut base_records: Vec<&mut PrimaryStructureRecord> = primary_structure
                        .records
                        .iter_mut()
                        .filter(|r| *r.get_record_type() == PrimaryStructureRecordType::Base)
                        .collect();
                    for chunk in base_records.chunks_mut(3) {
                        if chunk.len() < 3 {
                            break;
                        }
                        let mut nucleotides: String = String::with_capacity(3);
                        for record in &*chunk {
                            nucleotides.push_str(record.get_nucleotide().as_ref().unwrap().as_str().to_uppercase().as_str());
                        }
                        let frameshift_state: FrameshiftState = if chunk[2].get_net_variant_nucleotides_count() % 3 == 0 {
                            FrameshiftState::InFrame
                        } else {
                            FrameshiftState::FrameShifted
                        };
                        let amino_acid: &str = CODON_TABLE[nucleotides.replace("T", "U").as_str()];
                        for record in chunk {
                            record.set_amino_acid(Some(amino_acid.into()));
                            record.set_frameshift_state(Some(frameshift_state.clone()));
                        }
                    }

                    primary_structures.push(primary_structure);
                }

                primary_structures
            })
            .flatten()
            .collect::<Vec<PrimaryStructure>>()
    });

    let mut primary_structure_set: PrimaryStructureSet = PrimaryStructureSet::new();
    for primary_structure in primary_structures {
        primary_structure_set.add(primary_structure);
    }

    primary_structure_set
}
