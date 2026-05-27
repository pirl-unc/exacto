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
use indicatif::{ProgressBar, ProgressStyle};
use polars::prelude::*;
use rayon::prelude::*;
use std::sync::Arc;

use crate::prelude::*;


pub fn integrate_dna_rna_variants(
    dna_variant_records: &Vec<DNAVariantRecord>,
    rna_variant_records: &Vec<RNAVariantRecord>,
    gene_annotator: &(impl GeneAnnotator + Sync),
    max_exon_offset: u16,
    max_transcript_boundary_offset: u32,
    max_intergenic_distance: u32,
    num_threads: usize
) -> Vec<IntegratedVariant> {
    let dna_variant_records_index: DNAVariantIndex = DNAVariantIndex::new(dna_variant_records);
    let rna_variant_records_index: RNAVariantIndex = RNAVariantIndex::new(rna_variant_records);

    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    
    let pb: Arc<ProgressBar> = Arc::new(ProgressBar::new(rna_variant_records_index.len() as u64));

    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-")
    );

    let integrated_variants: Vec<IntegratedVariant> = thread_pool.install(|| {
        rna_variant_records
            .par_iter()
            .map(|rna_variant_record| {
                let reference_gene_names: Vec<Box<str>> = if rna_variant_record.reference_gene_name.is_empty() {
                    Vec::new()
                } else {
                    rna_variant_record
                        .reference_gene_name
                        .split(',')
                        .map(|s| Box::from(s))
                        .collect()
                };
                let reference_transcript_ids: Vec<Box<str>> = if rna_variant_record.reference_transcript_id.is_empty() {
                    Vec::new()
                } else {
                    rna_variant_record
                        .reference_transcript_id
                        .split(',')
                        .map(|s| Box::from(s))
                        .collect()
                };

                let mut integrated_variant: IntegratedVariant = IntegratedVariant::new(
                    rna_variant_record.assembled_transcript_name.clone(),
                    rna_variant_record.transcript_model_id,
                    &reference_gene_names,
                    &reference_transcript_ids,
                    rna_variant_record.variant_id
                );

                if reference_transcript_ids.is_empty() {
                    // Intergenic RNA variant — range-query DNA variants around
                    // each RNA endpoint and accept any whose endpoint falls
                    // within `max_intergenic_distance` bp.
                    let dna_records_1: Vec<&DNAVariantRecord> = dna_variant_records_index.get_by_range(
                        &rna_variant_record.chromosome_1,
                        rna_variant_record.position_1.saturating_sub(max_intergenic_distance),
                        rna_variant_record.position_1 + max_intergenic_distance
                    );
                    let dna_records_2: Vec<&DNAVariantRecord> = dna_variant_records_index.get_by_range(
                        &rna_variant_record.chromosome_2,
                        rna_variant_record.position_2.saturating_sub(max_intergenic_distance),
                        rna_variant_record.position_2 + max_intergenic_distance
                    );
                    for dna_variant_record in dna_records_1.iter().chain(dna_records_2.iter()) {
                        let dna_variant_id: u32 = dna_variant_record.variant_id;
                        if !integrated_variant.dna_variant_ids.contains_key(&dna_variant_id) {
                            let integrated_variant_distance: IntegratedVariantDistance = calculate_distance(
                                rna_variant_record,
                                dna_variant_record
                            );
                            integrated_variant.add_dna_variant_id(dna_variant_id, integrated_variant_distance);
                        }
                    }
                } else {
                    // Intragenic RNA variant — for each annotated reference
                    // transcript, range-query DNA variants that fall within
                    // the RT's genomic span, then apply the exon/intron/
                    // intergenic proximity rules.
                    for reference_transcript_id in reference_transcript_ids.iter() {
                        let reference_transcript: &Transcript = gene_annotator
                            .get_transcript(reference_transcript_id)
                            .unwrap();
                        let dna_variant_records_in_transcript: Vec<&DNAVariantRecord> = dna_variant_records_index.get_by_range(
                            &reference_transcript.chromosome,
                            reference_transcript.start,
                            reference_transcript.end
                        );
                        let (rna_genic_region_1, rna_exons_1) = reference_transcript.locate_position(rna_variant_record.position_1);
                        let (rna_genic_region_2, rna_exons_2) = reference_transcript.locate_position(rna_variant_record.position_2);

                        for dna_variant_record in dna_variant_records_in_transcript.iter() {
                            let (dna_genic_region_1, dna_exons_1) = reference_transcript.locate_position(dna_variant_record.position_1);
                            let (dna_genic_region_2, dna_exons_2) = reference_transcript.locate_position(dna_variant_record.position_2);

                            // RNA position 1 and DNA position 1
                            let integrate_11 = is_dna_variant_near_rna_variant(
                                rna_variant_record.position_1,
                                rna_exons_1.clone(),
                                rna_genic_region_1.clone(),
                                dna_variant_record.position_1,
                                dna_genic_region_1.clone(),
                                dna_exons_1.clone(),
                                reference_transcript,
                                max_exon_offset,
                                max_transcript_boundary_offset,
                                max_intergenic_distance
                            );

                            // RNA position 1 and DNA position 2
                            let integrate_12 = is_dna_variant_near_rna_variant(
                                rna_variant_record.position_1,
                                rna_exons_1.clone(),
                                rna_genic_region_1.clone(),
                                dna_variant_record.position_2,
                                dna_genic_region_2.clone(),
                                dna_exons_2.clone(),
                                reference_transcript,
                                max_exon_offset,
                                max_transcript_boundary_offset,
                                max_intergenic_distance
                            );

                            // RNA position 2 and DNA position 1
                            let integrate_21 = is_dna_variant_near_rna_variant(
                                rna_variant_record.position_2,
                                rna_exons_2.clone(),
                                rna_genic_region_2.clone(),
                                dna_variant_record.position_1,
                                dna_genic_region_1.clone(),
                                dna_exons_1.clone(),
                                reference_transcript,
                                max_exon_offset,
                                max_transcript_boundary_offset,
                                max_intergenic_distance
                            );

                            // RNA position 2 and DNA position 2
                            let integrate_22 = is_dna_variant_near_rna_variant(
                                rna_variant_record.position_2,
                                rna_exons_2.clone(),
                                rna_genic_region_2.clone(),
                                dna_variant_record.position_2,
                                dna_genic_region_2.clone(),
                                dna_exons_2.clone(),
                                reference_transcript,
                                max_exon_offset,
                                max_transcript_boundary_offset,
                                max_intergenic_distance
                            );

                            if integrate_11 || integrate_12 || integrate_21 || integrate_22 {
                                let dna_variant_id: u32 = dna_variant_record.variant_id;
                                if !integrated_variant.dna_variant_ids.contains_key(&dna_variant_id) {
                                    let integrated_variant_distance: IntegratedVariantDistance = calculate_distance(
                                        rna_variant_record,
                                        dna_variant_record
                                    );
                                    integrated_variant.add_dna_variant_id(dna_variant_id, integrated_variant_distance);
                                }
                            }
                        }
                    }
                }

                pb.inc(1);
                integrated_variant
            })
            .filter(|integrated_variant| !integrated_variant.dna_variant_ids.is_empty())
            .collect()
    });

    pb.finish_with_message("Completed annotating variant calls.");

    integrated_variants
}

fn are_exons_proximal(exon_a: &Exon, exon_b: &Exon, max_exon_offset: u16) -> bool {
    exon_a.exon_number.abs_diff(exon_b.exon_number) <= max_exon_offset
}

fn calculate_distance(
    rna_variant_record: &RNAVariantRecord,
    dna_variant_record: &DNAVariantRecord
) -> IntegratedVariantDistance {
    let mut min_distance: u32 = u32::MAX;
    let mut rna_variant_position_used: VariantPosition = VariantPosition::Position1;
    let mut dna_variant_position_used: VariantPosition = VariantPosition::Position1;

    // RNA variant position 1 vs DNA variant position 1
    if rna_variant_record.chromosome_1 == dna_variant_record.chromosome_1 {
        let curr_distance: u32 = rna_variant_record.position_1.abs_diff(dna_variant_record.position_1);
        if curr_distance < min_distance {
            min_distance = curr_distance;
            rna_variant_position_used = VariantPosition::Position1;
            dna_variant_position_used = VariantPosition::Position1;
        }
    }

    // RNA variant position 1 vs DNA variant position 2
    if rna_variant_record.chromosome_1 == dna_variant_record.chromosome_2 {
        let curr_distance: u32 = rna_variant_record.position_1.abs_diff(dna_variant_record.position_2);
        if curr_distance < min_distance {
            min_distance = curr_distance;
            rna_variant_position_used = VariantPosition::Position1;
            dna_variant_position_used = VariantPosition::Position2;
        }
    }

    // RNA variant position 2 vs DNA variant position 1
    if rna_variant_record.chromosome_2 == dna_variant_record.chromosome_1 {
        let curr_distance: u32 = rna_variant_record.position_2.abs_diff(dna_variant_record.position_1);
        if curr_distance < min_distance {
            min_distance = curr_distance;
            rna_variant_position_used = VariantPosition::Position2;
            dna_variant_position_used = VariantPosition::Position1;
        }
    }

    // RNA variant position 2 vs DNA variant position 2
    if rna_variant_record.chromosome_2 == dna_variant_record.chromosome_2 {
        let curr_distance: u32 = rna_variant_record.position_2.abs_diff(dna_variant_record.position_2);
        if curr_distance < min_distance {
            min_distance = curr_distance;
            rna_variant_position_used = VariantPosition::Position2;
            dna_variant_position_used = VariantPosition::Position2;
        }
    }

    assert_ne!(min_distance, u32::MAX);

    IntegratedVariantDistance::new(min_distance, rna_variant_position_used, dna_variant_position_used)
}

fn is_dna_variant_near_rna_variant(
    rna_position: u32,
    rna_exons: Option<(Exon, Exon)>,
    rna_region: GenicRegion,
    dna_position: u32,
    dna_region: GenicRegion,
    dna_exons: Option<(Exon, Exon)>,
    reference_transcript: &Transcript,
    max_exon_offset: u16,
    max_transcript_boundary_offset: u32,
    max_intergenic_distance: u32
) -> bool {
    match (rna_region, dna_region) {
        (GenicRegion::Exonic, GenicRegion::Exonic) => {
            if let (Some((rna_exon, _)), Some((dna_exon, _))) = (rna_exons.as_ref(), dna_exons.as_ref()) {
                return are_exons_proximal(rna_exon, dna_exon, max_exon_offset);
            }
        },
        (GenicRegion::Exonic, GenicRegion::Intronic) => {
            if let (Some((rna_exon, _)), Some((dna_exon_1, dna_exon_2))) = (rna_exons.as_ref(), dna_exons.as_ref()) {
                return are_exons_proximal(rna_exon, dna_exon_1, max_exon_offset) ||
                        are_exons_proximal(rna_exon, dna_exon_2, max_exon_offset);
            }
        },
        (GenicRegion::Exonic, GenicRegion::Intergenic) => {
            if let Some((rna_exon, _)) = rna_exons.as_ref() {
                let exon_count = reference_transcript.get_exon_ids().len();
                if rna_exon.exon_number <= max_exon_offset ||
                    rna_exon.exon_number as usize >= exon_count.saturating_sub(max_exon_offset as usize) {
                    return rna_position.abs_diff(dna_position) <= max_transcript_boundary_offset;
                }
            }
        },
        (GenicRegion::Intronic, GenicRegion::Exonic) => {
            if let (Some((rna_exon_1, rna_exon_2)), Some((dna_exon, _))) = (rna_exons.as_ref(), dna_exons.as_ref()) {
                return are_exons_proximal(rna_exon_1, dna_exon, max_exon_offset) ||
                        are_exons_proximal(rna_exon_2, dna_exon, max_exon_offset);
            }
        },
        (GenicRegion::Intronic, GenicRegion::Intronic) => {
            if let (Some((rna_exon_1, rna_exon_2)), Some((dna_exon_1, dna_exon_2))) = (rna_exons.as_ref(), dna_exons.as_ref()) {
                return are_exons_proximal(rna_exon_1, dna_exon_1, max_exon_offset) ||
                        are_exons_proximal(rna_exon_1, dna_exon_2, max_exon_offset) ||
                        are_exons_proximal(rna_exon_2, dna_exon_1, max_exon_offset) ||
                        are_exons_proximal(rna_exon_2, dna_exon_2, max_exon_offset);
            }
        },
        (GenicRegion::Intronic, GenicRegion::Intergenic) => {
            if let Some((rna_exon_1, rna_exon_2)) = rna_exons.as_ref() {
                let exon_count = reference_transcript.get_exon_ids().len();
                if ((rna_exon_1.exon_number <= max_exon_offset || rna_exon_1.exon_number as usize >= exon_count.saturating_sub(max_exon_offset as usize)) ||
                    (rna_exon_2.exon_number <= max_exon_offset || rna_exon_2.exon_number as usize >= exon_count.saturating_sub(max_exon_offset as usize))) &&
                    rna_position.abs_diff(dna_position) <= max_transcript_boundary_offset {
                    return true;
                } else {
                    return false;
                }
            }
        },
        (GenicRegion::Intergenic, GenicRegion::Exonic) => {
            // Check if the RNA variant is adjacent to the reference transcript
            if rna_position.abs_diff(reference_transcript.start) <= max_transcript_boundary_offset {
                let first_exon: &Exon = reference_transcript.get_exon_by_number(1).unwrap();
                if let Some((dna_exon, _)) = dna_exons.as_ref() {
                    return are_exons_proximal(first_exon, dna_exon, max_exon_offset)
                }
            }
            if rna_position.abs_diff(reference_transcript.end) <= max_transcript_boundary_offset {
                let last_exon: &Exon = reference_transcript.get_exon_by_number(reference_transcript.get_exon_ids().len() as u16).unwrap();
                if let Some((dna_exon, _)) = dna_exons.as_ref() {
                    return are_exons_proximal(last_exon, dna_exon, max_exon_offset)
                }
            }
            return rna_position.abs_diff(dna_position) <= max_intergenic_distance;
        },
        (GenicRegion::Intergenic, GenicRegion::Intronic) => {
            // Check if the RNA variant is adjacent to the reference transcript
            if rna_position.abs_diff(reference_transcript.start) <= max_transcript_boundary_offset {
                let first_exon: &Exon = reference_transcript.get_exon_by_number(1).unwrap();
                if let Some((dna_exon_1, dna_exon_2)) = dna_exons.as_ref() {
                    return are_exons_proximal(first_exon, dna_exon_1, max_exon_offset) ||
                            are_exons_proximal(first_exon, dna_exon_2, max_exon_offset);
                }
            }
            if rna_position.abs_diff(reference_transcript.end) <= max_transcript_boundary_offset {
                let last_exon: &Exon = reference_transcript.get_exon_by_number(reference_transcript.get_exon_ids().len() as u16).unwrap();
                if let Some((dna_exon_1, dna_exon_2)) = dna_exons.as_ref() {
                    return are_exons_proximal(last_exon, dna_exon_1, max_exon_offset) ||
                            are_exons_proximal(last_exon, dna_exon_2, max_exon_offset);
                }
            }
            return rna_position.abs_diff(dna_position) <= max_intergenic_distance;
        },
        (GenicRegion::Intergenic, GenicRegion::Intergenic) => {
            return rna_position.abs_diff(dna_position) <= max_intergenic_distance;
        }
    }

    false
}
