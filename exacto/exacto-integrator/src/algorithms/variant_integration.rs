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
use exacto_annotator::prelude::*;
use exacto_caller::prelude::*;
use exacto_core::prelude::*;
use indicatif::{ProgressBar, ProgressStyle};
use polars::prelude::*;
use rayon::prelude::*;
use std::sync::Arc;

use crate::prelude::*;


pub fn integrate_dna_rna_variants(
    dna_variant_call_annotation_set: &VariantCallAnnotationSet,
    rna_variant_call_set: &RNAVariantCallSet,
    gene_annotator: &(impl GeneAnnotator + Sync),
    max_exon_offset: u32,
    max_transcript_boundary_offset: usize,
    max_intergenic_distance: usize,
    num_threads: usize
) -> IntegratedVariantSet {
    // Step 1. Get RNA variant calls
    let mut rna_variant_calls: Vec<(usize, &Vec<Box<str>>, &VariantCall)> = Vec::new();
    for (transcript_model_id, inner_map) in rna_variant_call_set.variant_calls_index.iter() {
        for (reference_transcript_ids, variant_call_ids) in inner_map.iter() {
            for variant_call_id in variant_call_ids.iter() {
                rna_variant_calls.push((
                    *transcript_model_id,
                    reference_transcript_ids,
                    rna_variant_call_set.get_variant_call(*variant_call_id)
                ));
            }
        }
    }

    // Step 2. Integrate DNA and RNA variants
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let pb = Arc::new(ProgressBar::new(rna_variant_calls.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-")
    );
    let integrated_variants: Vec<IntegratedVariant> = thread_pool.install(|| {
        rna_variant_calls
            .par_iter()
            .map(|(transcript_model_id, reference_transcript_ids, rna_variant_call)| {
                let mut integrated_variant: IntegratedVariant = IntegratedVariant::new(
                    *transcript_model_id,
                    reference_transcript_ids,
                    rna_variant_call.id
                );
                let (rna_variant_record,_) = rna_variant_call.get_consensus_record();
                let chromosome_1: Box<str> = rna_variant_call_set.chromosome_names_map.get_by_right(&rna_variant_record.get_chromosome_1()).unwrap().clone();
                let chromosome_2: Box<str> = rna_variant_call_set.chromosome_names_map.get_by_right(&rna_variant_record.get_chromosome_2()).unwrap().clone();
                if reference_transcript_ids.is_empty() {
                    // Check if this intergenic RNA variant can be integrated with any DNA variants
                    let dna_variant_call_annotations_1: Vec<&VariantCallAnnotation> = dna_variant_call_annotation_set.get_by_range(
                        &*chromosome_1,
                        rna_variant_record.get_position_1().saturating_sub(max_intergenic_distance),
                        rna_variant_record.get_position_1() + max_intergenic_distance
                    );
                    let dna_variant_call_annotations_2: Vec<&VariantCallAnnotation> = dna_variant_call_annotation_set.get_by_range(
                        &*chromosome_2,
                        rna_variant_record.get_position_2().saturating_sub(max_intergenic_distance),
                        rna_variant_record.get_position_2() + max_intergenic_distance
                    );
                    for dna_variant_call_annotation in dna_variant_call_annotations_1.iter().chain(dna_variant_call_annotations_2.iter()) {
                        if integrated_variant.dna_variant_call_ids.contains_key(&dna_variant_call_annotation.variant_call_id) == false {
                            let integrated_variant_distance: IntegratedVariantDistance = calculate_distance(
                                rna_variant_call,
                                dna_variant_call_annotation,
                                &rna_variant_call_set.chromosome_names_map
                            );
                            integrated_variant.add_dna_variant_call_id(dna_variant_call_annotation.variant_call_id, integrated_variant_distance);
                        }
                    }
                } else {
                    // Integrate RNA variants within known reference transcripts with DNA variants within maximum exon offset
                    for reference_transcript_id in reference_transcript_ids.iter() {
                        let reference_transcript: &Transcript = gene_annotator.get_transcript(reference_transcript_id).unwrap();
                        if let Some(dna_variant_call_annotations) = dna_variant_call_annotation_set.get_by_reference_transcript(reference_transcript_id) {
                            let (rna_genic_region_1, rna_exons_1) = reference_transcript.locate_position(rna_variant_record.get_position_1());
                            let (rna_genic_region_2, rna_exons_2) = reference_transcript.locate_position(rna_variant_record.get_position_2());

                            for dna_variant_call_annotation in dna_variant_call_annotations.iter() {
                                let (dna_genic_region_1, dna_exons_1) = reference_transcript.locate_position(dna_variant_call_annotation.position_1);
                                let (dna_genic_region_2, dna_exons_2) = reference_transcript.locate_position(dna_variant_call_annotation.position_2);

                                // RNA position 1 and DNA position 1
                                let integrate_11 = is_dna_variant_near_rna_variant(
                                    rna_variant_record.get_position_1(),
                                    rna_exons_1.clone(),
                                    rna_genic_region_1.clone(),
                                    dna_variant_call_annotation.position_1,
                                    dna_genic_region_1.clone(),
                                    dna_exons_1.clone(),
                                    reference_transcript,
                                    max_exon_offset,
                                    max_transcript_boundary_offset,
                                    max_intergenic_distance
                                );

                                // RNA position 1 and DNA position 2
                                let integrate_12 = is_dna_variant_near_rna_variant(
                                    rna_variant_record.get_position_1(),
                                    rna_exons_1.clone(),
                                    rna_genic_region_1.clone(),
                                    dna_variant_call_annotation.position_2,
                                    dna_genic_region_2.clone(),
                                    dna_exons_2.clone(),
                                    reference_transcript,
                                    max_exon_offset,
                                    max_transcript_boundary_offset,
                                    max_intergenic_distance
                                );

                                // RNA position 2 and DNA position 1
                                let integrate_21 = is_dna_variant_near_rna_variant(
                                    rna_variant_record.get_position_2(),
                                    rna_exons_2.clone(),
                                    rna_genic_region_2.clone(),
                                    dna_variant_call_annotation.position_1,
                                    dna_genic_region_1.clone(),
                                    dna_exons_1.clone(),
                                    reference_transcript,
                                    max_exon_offset,
                                    max_transcript_boundary_offset,
                                    max_intergenic_distance
                                );

                                // RNA position 2 and DNA position 2
                                let integrate_22 = is_dna_variant_near_rna_variant(
                                    rna_variant_record.get_position_2(),
                                    rna_exons_2.clone(),
                                    rna_genic_region_2.clone(),
                                    dna_variant_call_annotation.position_2,
                                    dna_genic_region_2.clone(),
                                    dna_exons_2.clone(),
                                    reference_transcript,
                                    max_exon_offset,
                                    max_transcript_boundary_offset,
                                    max_intergenic_distance
                                );

                                if integrate_11 || integrate_12 || integrate_21 || integrate_22 {
                                    if integrated_variant.dna_variant_call_ids.contains_key(&dna_variant_call_annotation.variant_call_id) == false {
                                        let integrated_variant_distance: IntegratedVariantDistance = calculate_distance(
                                            rna_variant_call,
                                            dna_variant_call_annotation,
                                            &rna_variant_call_set.chromosome_names_map
                                        );
                                        integrated_variant.add_dna_variant_call_id(dna_variant_call_annotation.variant_call_id, integrated_variant_distance);
                                    }
                                }
                            }
                        }
                    }
                }

                pb.inc(1);

                Some(integrated_variant)
            })
            .filter_map(|v| v)
            .collect()
    });
    pb.finish_with_message("Completed annotating variant calls.");

    let mut integrated_variant_set: IntegratedVariantSet = IntegratedVariantSet::new();
    for integrated_variant in integrated_variants {
        if !integrated_variant.dna_variant_call_ids.is_empty() {
            integrated_variant_set.add_integrated_variant(integrated_variant);
        }
    }

    integrated_variant_set
}

fn are_exons_proximal(exon_a: &Exon, exon_b: &Exon, max_exon_offset: u32) -> bool {
    exon_a.exon_number.abs_diff(exon_b.exon_number) <= max_exon_offset
}

fn calculate_distance(
    rna_variant_call: &VariantCall,
    dna_variant_call_annotation: &VariantCallAnnotation,
    rna_chromosome_names_map: &BiMap<Box<str>, u16>
) -> IntegratedVariantDistance {
    let mut min_distance: usize = std::usize::MAX;
    let mut rna_variant_position_used: VariantPosition = VariantPosition::Position1;
    let mut dna_variant_position_used: VariantPosition = VariantPosition::Position1;

    let consensus_rna_variant_record: &VariantRecord = rna_variant_call.get_consensus_record().0;
    let rna_chromosome_1_name: Box<str> = rna_chromosome_names_map.get_by_right(&consensus_rna_variant_record.get_chromosome_1()).unwrap().clone();
    let rna_chromosome_2_name: Box<str> = rna_chromosome_names_map.get_by_right(&consensus_rna_variant_record.get_chromosome_2()).unwrap().clone();

    let dna_chromosome_1_name: Box<str> = dna_variant_call_annotation.chromosome_1.clone();
    let dna_chromosome_2_name: Box<str> = dna_variant_call_annotation.chromosome_2.clone();

    // RNA variant position 1 vs DNA variant position 1
    if rna_chromosome_1_name == dna_chromosome_1_name {
        let curr_distance: usize = consensus_rna_variant_record.get_position_1().abs_diff(dna_variant_call_annotation.position_1);
        if curr_distance < min_distance {
            min_distance = curr_distance;
            rna_variant_position_used = VariantPosition::Position1;
            dna_variant_position_used = VariantPosition::Position1;
        }
    }

    // RNA variant position 1 vs DNA variant position 2
    if rna_chromosome_1_name == dna_chromosome_2_name {
        let curr_distance: usize = consensus_rna_variant_record.get_position_1().abs_diff(dna_variant_call_annotation.position_2);
        if curr_distance < min_distance {
            min_distance = curr_distance;
            rna_variant_position_used = VariantPosition::Position1;
            dna_variant_position_used = VariantPosition::Position2;
        }
    }

    // RNA variant position 2 vs DNA variant position 1
    if rna_chromosome_2_name == dna_chromosome_1_name {
        let curr_distance: usize = consensus_rna_variant_record.get_position_2().abs_diff(dna_variant_call_annotation.position_1);
        if curr_distance < min_distance {
            min_distance = curr_distance;
            rna_variant_position_used = VariantPosition::Position2;
            dna_variant_position_used = VariantPosition::Position1;
        }
    }

    // RNA variant position 2 vs DNA variant position 2
    if rna_chromosome_2_name == dna_chromosome_2_name {
        let curr_distance: usize = consensus_rna_variant_record.get_position_2().abs_diff(dna_variant_call_annotation.position_2);
        if curr_distance < min_distance {
            min_distance = curr_distance;
            rna_variant_position_used = VariantPosition::Position2;
            dna_variant_position_used = VariantPosition::Position2;
        }
    }

    assert_ne!(min_distance, std::usize::MAX);

    IntegratedVariantDistance::new(min_distance, rna_variant_position_used, dna_variant_position_used)
}

fn is_dna_variant_near_rna_variant(
    rna_position: usize,
    rna_exons: Option<(Exon, Exon)>,
    rna_region: GenicRegion,
    dna_position: usize,
    dna_region: GenicRegion,
    dna_exons: Option<(Exon, Exon)>,
    reference_transcript: &Transcript,
    max_exon_offset: u32,
    max_transcript_boundary_offset: usize,
    max_intergenic_distance: usize
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
                    rna_exon.exon_number as usize >= (exon_count - max_exon_offset as usize) {
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
                if ((rna_exon_1.exon_number <= max_exon_offset || rna_exon_1.exon_number as usize >= (exon_count - max_exon_offset as usize)) ||
                    (rna_exon_2.exon_number <= max_exon_offset || rna_exon_2.exon_number as usize >= (exon_count - max_exon_offset as usize))) &&
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
                let last_exon: &Exon = reference_transcript.get_exon_by_number(reference_transcript.get_exon_ids().len() as u32).unwrap();
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
                let last_exon: &Exon = reference_transcript.get_exon_by_number(reference_transcript.get_exon_ids().len() as u32).unwrap();
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
