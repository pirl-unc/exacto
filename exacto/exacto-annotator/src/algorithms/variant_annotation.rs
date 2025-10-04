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
use std::collections::{HashMap, HashSet};
use std::str::FromStr;
use std::sync::Arc;

use crate::prelude::*;


pub fn annotate_position(
    chromosome_name: &str,
    position: u32,
    gene_annotator: &(impl GeneAnnotator + Sync)
) -> PositionAnnotation {
    let gene_ids: Vec<Box<str>> = gene_annotator.get_gene_ids_at_locus(&*chromosome_name, position);
    let mut transcript_ids: HashMap<Box<str>,HashSet<Box<str>>> = HashMap::new();
    let mut exon_ids: HashMap<Box<str>,Box<str>> = HashMap::new();

    if !gene_ids.is_empty() {
        let pos_isize = position as isize;
        for gene_id in gene_ids.iter() {
            if let Some(gene) = gene_annotator.get_gene(&*gene_id) {
                for transcript_id in gene.get_transcript_ids() {
                    if let Some(transcript) = gene_annotator.get_transcript(&*transcript_id) {
                        if overlaps(pos_isize, pos_isize, transcript.start as isize, transcript.end as isize) {
                            transcript_ids
                                .entry(gene_id.clone())
                                .or_insert_with(HashSet::new)
                                .insert(transcript_id.clone());
                        }
                        for exon_id in transcript.get_exon_ids() {
                            if let Some(exon) = gene_annotator.get_exon(&*exon_id) {
                                if overlaps(pos_isize, pos_isize, exon.start as isize, exon.end as isize) {
                                    transcript_ids
                                        .entry(gene_id.clone())
                                        .or_insert_with(HashSet::new)
                                        .insert(transcript_id.clone());
                                    exon_ids.insert(transcript_id.clone(), exon_id.clone());
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    let mut genic_region: GenicRegion = GenicRegion::Intergenic;
    if !transcript_ids.is_empty() {
        genic_region = GenicRegion::Intronic;
    }
    if !exon_ids.is_empty() {
        genic_region = GenicRegion::Exonic;
    }

    let mut position_annotation: PositionAnnotation = PositionAnnotation::new(genic_region);

    for gene_id in gene_ids.iter() {
        position_annotation.add_gene_id(gene_id);
    }
    for (gene_id, transcript_ids) in transcript_ids.iter() {
        for transcript_id in transcript_ids.iter() {
            position_annotation.add_transcript_id(gene_id, transcript_id);
        }
    }
    for (transcript_id, exon_id) in exon_ids.iter() {
        position_annotation.add_exon_id(transcript_id, exon_id);
    }

    position_annotation
}

/// Annotate variant calls.
///
/// # Arguments
///
/// * `df_variant_calls` - A reference to a `DataFrame` with the following columns: 
///     - `variant_call_id`
///     - `chromosome_1` 
///     - `position_1`
///     - `chromosome_2`
///     - `position_2`
///     - `variant_type`
///     - `variant_sequence`
/// * `gene_annotator` - A reference to an object implementing the `GeneAnnotator` trait, which provides
///   methods for querying gene and transcript annotations.
/// * `num_threads` - The number of threads to use for parallel annotation. Useful for improving performance
///   on large datasets.
pub fn annotate_variant_calls(
    df_variant_calls: &DataFrame,
    gene_annotator: &(impl GeneAnnotator + Sync),
    num_threads: usize
) -> VariantCallAnnotationSet {
    let variant_call_id_col = df_variant_calls.column("variant_call_id").unwrap().i64().unwrap();
    let chromosome_1_col = df_variant_calls.column("chromosome_1").unwrap().str().unwrap();
    let position_1_col = df_variant_calls.column("position_1").unwrap().i64().unwrap();
    let chromosome_2_col = df_variant_calls.column("chromosome_2").unwrap().str().unwrap();
    let position_2_col = df_variant_calls.column("position_2").unwrap().i64().unwrap();
    let variant_type_col = df_variant_calls.column("variant_type").unwrap().str().unwrap();
    let variant_sequence_col = df_variant_calls.column("variant_sequence").unwrap().str().unwrap();

    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    
    let pb = Arc::new(ProgressBar::new(df_variant_calls.height() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-")
    );
    
    let variant_call_annotations: Vec<VariantCallAnnotation> = thread_pool.install(|| {
        (0..df_variant_calls.height())
            .into_par_iter()
            .map(|i| {
                let variant_call_id: usize = variant_call_id_col.get(i).unwrap() as usize;
                
                // Position 1
                let position_1: u32 = position_1_col.get(i).unwrap() as u32;
                let chromosome_1: &str = chromosome_1_col.get(i).unwrap();
                let position_annotation_1: PositionAnnotation = annotate_position(
                    chromosome_1,
                    position_1,
                    gene_annotator
                );

                // Position 2
                let position_2: u32 = position_2_col.get(i).unwrap() as u32;
                let chromosome_2: &str = chromosome_2_col.get(i).unwrap();
                let position_annotation_2: PositionAnnotation = annotate_position(
                    chromosome_2,
                    position_2,
                    gene_annotator
                );

                let variant_type: VariantType = VariantType::from_str(variant_type_col.get(i).unwrap()).unwrap();
                let variant_sequence: &str = variant_sequence_col.get(i).unwrap();

                pb.inc(1);

                VariantCallAnnotation::new(
                    variant_call_id,
                    chromosome_1,
                    position_1,
                    chromosome_2,
                    position_2,
                    variant_type,
                    variant_sequence,
                    position_annotation_1,
                    position_annotation_2
                )
            })
            .collect()
    });
    
    pb.finish_with_message("Completed annotating variant calls.");
    
    let mut variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::new();
    for variant_call_annotation in variant_call_annotations {
        variant_call_annotation_set.add_annotation(variant_call_annotation);
    }

    variant_call_annotation_set
}
