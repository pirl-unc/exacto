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
use std::collections::HashMap;
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

    // Step 1. Pre-build integrated variants lookup:
    // (transcript_model_id, reference_transcript_ids_str, rna_variant_call_id) -> Vec<dna_variant_call_id>
    let mut integrated_variants_index: HashMap<(usize, Box<str>, usize), Vec<usize>> = HashMap::new();
    if df_integrated_variants.height() > 0 {
        let iv_col_tmi = df_integrated_variants.column("transcript_model_id").unwrap();
        let iv_col_rti = df_integrated_variants.column("reference_transcript_ids").unwrap().str().unwrap();
        let iv_col_rna = df_integrated_variants.column("rna_variant_call_id").unwrap();
        let iv_col_dna = df_integrated_variants.column("dna_variant_call_id").unwrap();
        for i in 0..df_integrated_variants.height() {
            let tmi: usize = if let Ok(c) = iv_col_tmi.i64() {
                c.get(i).unwrap() as usize
            } else {
                iv_col_tmi.str().unwrap().get(i).unwrap().parse::<usize>().unwrap()
            };
            let rti: Box<str> = iv_col_rti.get(i).unwrap().into();
            let rna_id: usize = if let Ok(c) = iv_col_rna.i64() {
                c.get(i).unwrap() as usize
            } else {
                iv_col_rna.str().unwrap().get(i).unwrap().parse::<usize>().unwrap()
            };
            let dna_id: usize = if let Ok(c) = iv_col_dna.i64() {
                c.get(i).unwrap() as usize
            } else {
                iv_col_dna.str().unwrap().get(i).unwrap().parse::<usize>().unwrap()
            };
            integrated_variants_index
                .entry((tmi, rti, rna_id))
                .or_default()
                .push(dna_id);
        }
    }

    // Step 2. Extract row data from DataFrame into lightweight vectors per partition
    let df_sorted_all: DataFrame = df_transcript_structures
        .sort(["transcript_model_id", "reference_transcript_ids", "index"], SortMultipleOptions::default())
        .unwrap();

    let col_tmi_series = df_sorted_all.column("transcript_model_id").unwrap().cast(&DataType::String).unwrap();
    let col_tmi_str = col_tmi_series.str().unwrap();
    let col_rti_str = df_sorted_all.column("reference_transcript_ids").unwrap().str().unwrap();
    let col_type = df_sorted_all.column("type").unwrap().str().unwrap();
    let col_read_start = df_sorted_all.column("read_start").unwrap().i64().unwrap();
    let col_read_end = df_sorted_all.column("read_end").unwrap().i64().unwrap();
    let col_sequence = df_sorted_all.column("sequence").unwrap().str().unwrap();
    let col_tsi = df_sorted_all.column("index").unwrap().i64().unwrap();
    let col_kind = df_sorted_all.column("kind").unwrap().str().unwrap();
    let col_context = df_sorted_all.column("context").unwrap().str().unwrap();
    let col_skipped = df_sorted_all.column("skipped").unwrap().str().unwrap();
    let col_tmi_val = df_sorted_all.column("transcript_model_id").unwrap();

    struct RowData {
        record_type: Box<str>,
        read_start: u32,
        read_end: u32,
        sequence: Box<str>,
        transcript_structure_index: usize,
        kind: Box<str>,
        context: Box<str>,
        skipped: Box<str>,
    }

    struct PartitionData {
        transcript_model_id: usize,
        reference_transcript_ids: Vec<Box<str>>,
        reference_transcript_ids_str: Box<str>,
        rna_sequence: String,
        rows: Vec<RowData>
    }

    let mut partitions: Vec<PartitionData> = Vec::new();
    if df_sorted_all.height() > 0 {
        let mut group_start: usize = 0;
        for i in 0..=df_sorted_all.height() {
            let is_boundary = i == df_sorted_all.height()
                || col_tmi_str.get(i) != col_tmi_str.get(group_start)
                || col_rti_str.get(i) != col_rti_str.get(group_start);

            if is_boundary && i > group_start {
                let tmi: usize = if let Ok(c) = col_tmi_val.i64() {
                    c.get(group_start).unwrap() as usize
                } else {
                    col_tmi_val.str().unwrap().get(group_start).unwrap().parse::<usize>().unwrap()
                };
                let rti_str: Box<str> = col_rti_str.get(group_start).unwrap().into();
                let ref_ids: Vec<Box<str>> = col_rti_str.get(group_start).unwrap().split(",").map(|s| s.into()).collect();

                let mut rna_sequence = String::new();
                let mut rows: Vec<RowData> = Vec::new();
                for j in group_start..i {
                    let seq = col_sequence.get(j).unwrap_or("");
                    rna_sequence.push_str(seq);
                    rows.push(RowData {
                        record_type: col_type.get(j).unwrap_or("").into(),
                        read_start: col_read_start.get(j).unwrap() as u32,
                        read_end: col_read_end.get(j).unwrap() as u32,
                        sequence: seq.into(),
                        transcript_structure_index: col_tsi.get(j).unwrap() as usize,
                        kind: col_kind.get(j).unwrap_or("").into(),
                        context: col_context.get(j).unwrap_or("").into(),
                        skipped: col_skipped.get(j).unwrap_or("").into(),
                    });
                }

                partitions.push(PartitionData {
                    transcript_model_id: tmi,
                    reference_transcript_ids: ref_ids,
                    reference_transcript_ids_str: rti_str,
                    rna_sequence,
                    rows,
                });
                group_start = i;
            }
        }
    }

    // Step 3. Process partitions in parallel — no DataFrame access in threads
    let primary_structures: Vec<PrimaryStructure> = thread_pool.install(|| {
        partitions
            .into_par_iter()
            .map(|partition| {
                let rna: RNA = RNA::new("".into(), partition.rna_sequence.into());
                let translation_result: Option<Translation> = rna.translate();
                if translation_result.is_none() {
                    return Vec::new();
                }
                let translation: Translation = translation_result.unwrap();
                // Convert orf_start/orf_end from 0-based (translate()) to 1-based (read positions)
                let first_read_start = partition.rows.first().map(|r| r.read_start).unwrap_or(1);
                let peptides: Vec<Peptide> = match translation_strategy {
                    TranslationStrategy::LongestORF => {
                        let p = translation.get_longest_orf_peptide();
                        vec![Peptide::new(p.sequence.clone(), p.orf_start + first_read_start, p.orf_end + first_read_start)]
                    },
                    TranslationStrategy::AllORFs => {
                        translation.get_peptides().iter().map(|p| {
                            Peptide::new(p.sequence.clone(), p.orf_start + first_read_start, p.orf_end + first_read_start)
                        }).collect()
                    }
                };

                // Fetch RNA variant calls
                let mut rna_variant_calls: Vec<&VariantCall> = Vec::new();
                if let Some(gene_map) = rna_variant_call_set.variant_calls_index.get(&partition.transcript_model_id) {
                    if let Some(call_ids) = gene_map.get(&partition.reference_transcript_ids) {
                        for variant_call_id in call_ids {
                            rna_variant_calls.push(rna_variant_call_set.get_variant_call(*variant_call_id));
                        }
                    }
                }

                let mut primary_structures: Vec<PrimaryStructure> = Vec::new();
                for peptide in peptides.iter() {
                    let mut primary_structure: PrimaryStructure = PrimaryStructure::new();
                    let mut primary_structure_index: usize = 0;
                    let mut codon_index: u8 = 0;
                    let mut net_variant_nucleotides_count: i32 = 0;

                    for row in partition.rows.iter() {
                        if row.record_type.is_empty() { continue; }
                        let record_type: AlignmentStructureRecordType = AlignmentStructureRecordType::from_str(&row.record_type).unwrap();
                        let read_start = row.read_start;
                        let read_end = row.read_end;

                        if overlaps(
                            peptide.orf_start as isize,
                            peptide.orf_end as isize,
                            read_start as isize,
                            read_end as isize) == false {
                            continue;
                        }

                        // Get any RNA and DNA variants matching the current row
                        let mut rna_variant_call_ids: Vec<usize> = Vec::new();
                        let mut dna_variant_call_ids: Vec<usize> = Vec::new();
                        for rna_variant_call in rna_variant_calls.iter() {
                            let rna_variant_record: &VariantRecord = rna_variant_call.get_consensus_record().0;
                            if rna_variant_record.get_read_position_1() == read_start && rna_variant_record.get_read_position_2() == read_end {
                                rna_variant_call_ids.push(rna_variant_call.id);
                                let key = (partition.transcript_model_id, partition.reference_transcript_ids_str.clone(), rna_variant_call.id);
                                if let Some(dna_ids) = integrated_variants_index.get(&key) {
                                    dna_variant_call_ids.extend(dna_ids);
                                }
                            }
                        }

                        match record_type {
                            AlignmentStructureRecordType::Base => {
                                let kind: AlignmentStructureBaseKind = AlignmentStructureBaseKind::from_str(&row.kind).unwrap();
                                let context: AlignmentStructureBaseContext = AlignmentStructureBaseContext::from_str(&row.context).unwrap();

                                let mut read_position: u32 = read_start;
                                let mut chars = row.sequence.char_indices().peekable();
                                while let Some((start, _)) = chars.next() {
                                    if read_position >= peptide.orf_start &&
                                        read_position <= peptide.orf_end {
                                        let end: usize = chars.peek().map(|(k, _)| *k).unwrap_or(row.sequence.len());
                                        let nucleotide: &str = &row.sequence[start..end];

                                        let is_mutant_base: bool = match (&kind, &context) {
                                            (AlignmentStructureBaseKind::Insertion, _) => {
                                                net_variant_nucleotides_count += 1;
                                                true
                                            },
                                            (AlignmentStructureBaseKind::Match, AlignmentStructureBaseContext::Intronic) => {
                                                net_variant_nucleotides_count += 1;
                                                true
                                            },
                                            (AlignmentStructureBaseKind::Match, AlignmentStructureBaseContext::Intergenic) => {
                                                net_variant_nucleotides_count += 1;
                                                true
                                            },
                                            (AlignmentStructureBaseKind::Unaligned, _) => {
                                                net_variant_nucleotides_count += 1;
                                                true
                                            },
                                            (AlignmentStructureBaseKind::Mismatch, _) => {
                                                true
                                            },
                                            (_, _) => { false }
                                        };

                                        primary_structure.add_record(
                                            PrimaryStructureRecord::new(
                                                primary_structure_index,
                                                PrimaryStructureRecordType::Base,
                                                None,
                                                Some(codon_index),
                                                Some(Nucleotide::from_str(nucleotide).unwrap()),
                                                partition.transcript_model_id,
                                                partition.reference_transcript_ids.clone(),
                                                row.transcript_structure_index,
                                                read_position,
                                                read_position,
                                                net_variant_nucleotides_count,
                                                None,
                                                rna_variant_call_ids.clone(),
                                                dna_variant_call_ids.clone(),
                                                is_mutant_base
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
                                if row.skipped.is_empty() == false {
                                    for skipped_coordinates in row.skipped.split(",") {
                                        let skipped_elements: Vec<&str> = skipped_coordinates.split("|").collect();
                                        let position_1: usize = skipped_elements[0].split(":").collect::<Vec<&str>>()[1].parse::<usize>().unwrap();
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
                                        partition.transcript_model_id,
                                        partition.reference_transcript_ids.clone(),
                                        row.transcript_structure_index,
                                        read_start,
                                        read_end,
                                        net_variant_nucleotides_count,
                                        None,
                                        rna_variant_call_ids,
                                        dna_variant_call_ids,
                                        false
                                    )
                                );

                                primary_structure_index += 1;
                            }
                        }
                    }

                    // Assign the amino acid, frameshift states, codon variant IDs, and frameshift variant IDs
                    let mut base_records: Vec<&mut PrimaryStructureRecord> = primary_structure
                        .records
                        .iter_mut()
                        .filter(|r| *r.get_record_type() == PrimaryStructureRecordType::Base)
                        .collect();
                    let mut prev_frameshift_state = FrameshiftState::InFrame;
                    let mut active_frameshift_rna_ids: Vec<usize> = Vec::new();
                    let mut active_frameshift_dna_ids: Vec<usize> = Vec::new();
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

                        // Codon-level variant propagation: collect all variant IDs across the 3 bases
                        let mut codon_rna_ids: Vec<usize> = Vec::new();
                        let mut codon_dna_ids: Vec<usize> = Vec::new();
                        for record in &*chunk {
                            codon_rna_ids.extend(record.get_rna_variant_call_ids().iter());
                            codon_dna_ids.extend(record.get_dna_variant_call_ids().iter());
                        }
                        codon_rna_ids.sort_unstable();
                        codon_rna_ids.dedup();
                        codon_dna_ids.sort_unstable();
                        codon_dna_ids.dedup();

                        // Frameshift variant propagation: track causal variants across codons
                        match (&prev_frameshift_state, &frameshift_state) {
                            (FrameshiftState::InFrame, FrameshiftState::FrameShifted) => {
                                // Transition into frameshift: this codon's variants caused it
                                for record in &*chunk {
                                    active_frameshift_rna_ids.extend(record.get_rna_variant_call_ids().iter());
                                    active_frameshift_dna_ids.extend(record.get_dna_variant_call_ids().iter());
                                }
                                active_frameshift_rna_ids.sort_unstable();
                                active_frameshift_rna_ids.dedup();
                                active_frameshift_dna_ids.sort_unstable();
                                active_frameshift_dna_ids.dedup();
                            },
                            (FrameshiftState::FrameShifted, FrameshiftState::FrameShifted) => {
                                // Still frameshifted: accumulate any new variants from this codon
                                let has_new_variants = chunk.iter().any(|r| !r.get_rna_variant_call_ids().is_empty());
                                if has_new_variants {
                                    for record in &*chunk {
                                        active_frameshift_rna_ids.extend(record.get_rna_variant_call_ids().iter());
                                        active_frameshift_dna_ids.extend(record.get_dna_variant_call_ids().iter());
                                    }
                                    active_frameshift_rna_ids.sort_unstable();
                                    active_frameshift_rna_ids.dedup();
                                    active_frameshift_dna_ids.sort_unstable();
                                    active_frameshift_dna_ids.dedup();
                                }
                            },
                            (FrameshiftState::FrameShifted, FrameshiftState::InFrame) => {
                                // Frame restored: clear causal variant tracking
                                active_frameshift_rna_ids.clear();
                                active_frameshift_dna_ids.clear();
                            },
                            _ => {}
                        }
                        prev_frameshift_state = frameshift_state.clone();

                        // Determine amino_acid_change: mutant if any base is mutant or codon is frameshifted
                        let amino_acid_change = if frameshift_state == FrameshiftState::FrameShifted
                            || chunk.iter().any(|r| r.is_mutant_base()) {
                            AminoAcidChange::Mutant
                        } else {
                            AminoAcidChange::Reference
                        };

                        for record in chunk {
                            record.set_amino_acid(Some(amino_acid.into()));
                            record.set_frameshift_state(Some(frameshift_state.clone()));
                            record.set_amino_acid_change(Some(amino_acid_change.clone()));
                            record.set_codon_rna_variant_call_ids(codon_rna_ids.clone());
                            record.set_codon_dna_variant_call_ids(codon_dna_ids.clone());
                            if frameshift_state == FrameshiftState::FrameShifted {
                                record.set_frameshift_rna_variant_call_ids(active_frameshift_rna_ids.clone());
                                record.set_frameshift_dna_variant_call_ids(active_frameshift_dna_ids.clone());
                            }
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
