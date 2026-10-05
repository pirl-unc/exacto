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
use std::collections::{BTreeSet, HashSet};
use std::str::FromStr;
use std::sync::Arc;

use crate::annotation::variant_annotation::annotate_position;
use crate::prelude::*;


/// Annotate variant calls.
///
/// # Arguments
/// * `df_variant_calls` - A reference to a `DataFrame` with the following columns:
///     - `variant_id`
///     - `chromosome_1`
///     - `position_1`
///     - `chromosome_2`
///     - `position_2`
///     - `variant_type`
///     - `sequence`
/// * `gene_annotator` - A reference to an object implementing the `GeneAnnotator` trait, which provides
///   methods for querying gene and transcript annotations.
/// * `num_threads` - The number of threads to use for parallel annotation. Useful for improving performance
///   on large datasets.
///
/// # Returns
/// * `VariantCallAnnotationSet` object.
pub fn annotate_variant_calls(
    df_variant_calls: &DataFrame,
    gene_annotator: &(impl GeneAnnotator + Sync),
    num_threads: usize
) -> VariantCallAnnotationSet {
    let variant_call_annotations: Vec<VariantCallAnnotation> = annotate_variant_call_rows(
        df_variant_calls,
        gene_annotator,
        num_threads
    );

    let mut variant_call_annotation_set: VariantCallAnnotationSet = VariantCallAnnotationSet::new();
    for variant_call_annotation in variant_call_annotations {
        variant_call_annotation_set.add_annotation(variant_call_annotation);
    }

    variant_call_annotation_set
}

/// Annotate variant calls and append the annotation to the table.
///
/// # Arguments
/// * `df_variant_calls` - A reference to a `DataFrame` with the following columns. Any other
///   column is carried to the result unchanged.
///     - `variant_id`
///     - `chromosome_1`
///     - `position_1`
///     - `chromosome_2`
///     - `position_2`
///     - `variant_type`
///     - `sequence`
/// * `gene_annotator` - A reference to an object implementing the `GeneAnnotator` trait, which provides
///   methods for querying gene and transcript annotations.
/// * `num_threads` - The number of threads to use for parallel annotation. Useful for improving performance
///   on large datasets.
///
/// # Returns
/// * `df_variant_calls` with its columns and rows in their order, and the following columns appended.
///   A column of `df_variant_calls` that has one of these names is replaced.
///     - `position_1_genic_region`
///     - `position_1_annotation`
///     - `position_1_plus_1_genic_region`
///     - `position_1_plus_1_annotation`
///     - `position_2_minus_1_genic_region`
///     - `position_2_minus_1_annotation`
///     - `position_2_genic_region`
///     - `position_2_annotation`
pub fn append_variant_call_annotations(
    df_variant_calls: &DataFrame,
    gene_annotator: &(impl GeneAnnotator + Sync),
    num_threads: usize
) -> DataFrame {
    let variant_call_annotations: Vec<VariantCallAnnotation> = annotate_variant_call_rows(
        df_variant_calls,
        gene_annotator,
        num_threads
    );

    let mut position_1_genic_region_values: Vec<&str> = Vec::new();
    let mut position_1_annotation_values: Vec<String> = Vec::new();
    let mut position_1_plus_1_genic_region_values: Vec<&str> = Vec::new();
    let mut position_1_plus_1_annotation_values: Vec<String> = Vec::new();
    let mut position_2_minus_1_genic_region_values: Vec<&str> = Vec::new();
    let mut position_2_minus_1_annotation_values: Vec<String> = Vec::new();
    let mut position_2_genic_region_values: Vec<&str> = Vec::new();
    let mut position_2_annotation_values: Vec<String> = Vec::new();
    for variant_call_annotation in variant_call_annotations.iter() {
        position_1_genic_region_values.push(variant_call_annotation.position_1_annotation.genic_region.as_str());
        position_1_annotation_values.push(variant_call_annotation.position_1_annotation.to_string());
        position_1_plus_1_genic_region_values.push(variant_call_annotation.position_1_plus_1_annotation.genic_region.as_str());
        position_1_plus_1_annotation_values.push(variant_call_annotation.position_1_plus_1_annotation.to_string());
        position_2_minus_1_genic_region_values.push(variant_call_annotation.position_2_minus_1_annotation.genic_region.as_str());
        position_2_minus_1_annotation_values.push(variant_call_annotation.position_2_minus_1_annotation.to_string());
        position_2_genic_region_values.push(variant_call_annotation.position_2_annotation.genic_region.as_str());
        position_2_annotation_values.push(variant_call_annotation.position_2_annotation.to_string());
    }

    // The annotation columns go to the end of the table in the order of the positions, also when
    // the table was annotated before
    let mut df_annotated: DataFrame = df_variant_calls.drop_many([
        "position_1_genic_region",
        "position_1_annotation",
        "position_1_plus_1_genic_region",
        "position_1_plus_1_annotation",
        "position_2_minus_1_genic_region",
        "position_2_minus_1_annotation",
        "position_2_genic_region",
        "position_2_annotation"
    ]);
    df_annotated.with_column(Series::new("position_1_genic_region".into(), position_1_genic_region_values)).unwrap();
    df_annotated.with_column(Series::new("position_1_annotation".into(), position_1_annotation_values)).unwrap();
    df_annotated.with_column(Series::new("position_1_plus_1_genic_region".into(), position_1_plus_1_genic_region_values)).unwrap();
    df_annotated.with_column(Series::new("position_1_plus_1_annotation".into(), position_1_plus_1_annotation_values)).unwrap();
    df_annotated.with_column(Series::new("position_2_minus_1_genic_region".into(), position_2_minus_1_genic_region_values)).unwrap();
    df_annotated.with_column(Series::new("position_2_minus_1_annotation".into(), position_2_minus_1_annotation_values)).unwrap();
    df_annotated.with_column(Series::new("position_2_genic_region".into(), position_2_genic_region_values)).unwrap();
    df_annotated.with_column(Series::new("position_2_annotation".into(), position_2_annotation_values)).unwrap();

    df_annotated
}


/// Compares the chromosome names of the variant calls with those of the gene annotation.
///
/// # Arguments
/// * `col_chromosome_1` - Column `chromosome_1` of the table of variant calls.
/// * `col_chromosome_2` - Column `chromosome_2` of the table of variant calls.
/// * `gene_annotator` - A reference to an object implementing the `GeneAnnotator` trait. Its
///   chromosome names are the chromosomes of its genes.
///
/// # Panics
/// * If no chromosome name of the variant calls is a chromosome name of the gene annotation.
///
/// # Returns
/// * The chromosome names of the variant calls that the gene annotation does not have, sorted.
fn get_chromosomes_missing_from_gene_annotation(
    col_chromosome_1: &StringChunked,
    col_chromosome_2: &StringChunked,
    gene_annotator: &(impl GeneAnnotator + Sync)
) -> Vec<Box<str>> {
    let mut gene_annotation_chromosomes: BTreeSet<&str> = BTreeSet::new();
    for gene in gene_annotator.get_genes() {
        gene_annotation_chromosomes.insert(&*gene.chromosome);
    }

    let mut chromosomes: BTreeSet<&str> = BTreeSet::new();
    for chromosome in col_chromosome_1.into_no_null_iter() {
        chromosomes.insert(chromosome);
    }
    for chromosome in col_chromosome_2.into_no_null_iter() {
        chromosomes.insert(chromosome);
    }

    let missing_chromosomes: Vec<Box<str>> = chromosomes
        .iter()
        .filter(|chromosome| !gene_annotation_chromosomes.contains(*chromosome))
        .map(|chromosome| (*chromosome).into())
        .collect();

    assert!(
        missing_chromosomes.len() < chromosomes.len(),
        "df_variant_calls and the gene annotation share no chromosome name. df_variant_calls has {}. The gene annotation has {}.",
        chromosomes.iter().copied().collect::<Vec<&str>>().join(", "),
        gene_annotation_chromosomes.iter().copied().collect::<Vec<&str>>().join(", ")
    );

    missing_chromosomes
}

/// Annotate each row of a table of variant calls.
///
/// # Arguments
/// * `df_variant_calls` - A reference to a `DataFrame` with the following columns:
///     - `variant_id`
///     - `chromosome_1`
///     - `position_1`
///     - `chromosome_2`
///     - `position_2`
///     - `variant_type`
///     - `sequence`
/// * `gene_annotator` - A reference to an object implementing the `GeneAnnotator` trait, which provides
///   methods for querying gene and transcript annotations.
/// * `num_threads` - The number of threads to use for parallel annotation. Useful for improving performance
///   on large datasets.
///
/// # Panics
/// * If a chromosome name in `df_variant_calls` does not start with `chr`.
/// * If no chromosome name in `df_variant_calls` is a chromosome name of the gene annotation. A
///   name that the gene annotation does not have, among others that it has, is a warning.
/// * If a `variant_id` is on more than one row of `df_variant_calls`.
///
/// # Returns
/// * One `VariantCallAnnotation` per row of `df_variant_calls`, in the order of the rows. Each holds
///   the annotation of four positions: `position_1`, `position_1 + 1`, `position_2 - 1` and
///   `position_2`. `position_1` and `position_2` are the unchanged bases on each side of the variant,
///   so the two positions between them are the first and the last altered base of an SNV, an MNV
///   or a deletion.
fn annotate_variant_call_rows(
    df_variant_calls: &DataFrame,
    gene_annotator: &(impl GeneAnnotator + Sync),
    num_threads: usize
) -> Vec<VariantCallAnnotation> {
    // A table without rows has no variant call to annotate. It may have no column either.
    if df_variant_calls.height() == 0 {
        return Vec::new();
    }

    // IDs and positions are read as integers and chromosome names as strings, whatever type the
    // table holds them in
    let col_variant_call_id = df_variant_calls.column("variant_id").unwrap().strict_cast(&DataType::Int64).unwrap();
    let col_variant_call_id = col_variant_call_id.i64().unwrap();
    let col_chromosome_1 = df_variant_calls.column("chromosome_1").unwrap().strict_cast(&DataType::String).unwrap();
    let col_chromosome_1 = col_chromosome_1.str().unwrap();
    let col_position_1 = df_variant_calls.column("position_1").unwrap().strict_cast(&DataType::Int64).unwrap();
    let col_position_1 = col_position_1.i64().unwrap();
    let col_chromosome_2 = df_variant_calls.column("chromosome_2").unwrap().strict_cast(&DataType::String).unwrap();
    let col_chromosome_2 = col_chromosome_2.str().unwrap();
    let col_position_2 = df_variant_calls.column("position_2").unwrap().strict_cast(&DataType::Int64).unwrap();
    let col_position_2 = col_position_2.i64().unwrap();
    let col_variant_type = df_variant_calls.column("variant_type").unwrap().str().unwrap();
    let col_variant_sequence = df_variant_calls.column("sequence").unwrap().str().unwrap();

    // Every chromosome name starts with "chr"
    for i in 0..df_variant_calls.height() {
        let chromosome_1: &str = col_chromosome_1.get(i).unwrap();
        assert!(
            chromosome_1.starts_with("chr"),
            "chromosome_1 '{}' on row {} of df_variant_calls does not start with 'chr'.",
            chromosome_1,
            i + 1
        );
        let chromosome_2: &str = col_chromosome_2.get(i).unwrap();
        assert!(
            chromosome_2.starts_with("chr"),
            "chromosome_2 '{}' on row {} of df_variant_calls does not start with 'chr'.",
            chromosome_2,
            i + 1
        );
    }

    // The variant calls and the gene annotation share a chromosome name
    for chromosome in get_chromosomes_missing_from_gene_annotation(col_chromosome_1, col_chromosome_2, gene_annotator) {
        init_logging(true);
        log::warn!(
            "Chromosome '{}' of df_variant_calls is not in the gene annotation. Its positions are annotated as intergenic.",
            chromosome
        );
    }

    // There is exactly one row for each variant ID
    let mut variant_call_ids: HashSet<i64> = HashSet::new();
    for i in 0..df_variant_calls.height() {
        let variant_call_id: i64 = col_variant_call_id.get(i).unwrap();
        assert!(
            variant_call_ids.insert(variant_call_id),
            "variant_id {} is on more than one row of df_variant_calls (row {} is the second).",
            variant_call_id,
            i + 1
        );
    }

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
                let variant_call_id: usize = col_variant_call_id.get(i).unwrap() as usize;
                
                // Position 1 and the base after it
                let position_1: ReferencePosition = col_position_1.get(i).unwrap() as u32;
                let chromosome_1: ReferenceChromosomeName = col_chromosome_1.get(i).unwrap().into();
                let position_annotation_1: PositionAnnotation = annotate_position(
                    &*chromosome_1,
                    position_1,
                    gene_annotator
                );
                let position_annotation_1_plus_1: PositionAnnotation = annotate_position(
                    &*chromosome_1,
                    position_1.saturating_add(1),
                    gene_annotator
                );

                // Position 2 and the base before it
                let position_2: ReferencePosition = col_position_2.get(i).unwrap() as u32;
                let chromosome_2: ReferenceChromosomeName = col_chromosome_2.get(i).unwrap().into();
                let position_annotation_2_minus_1: PositionAnnotation = annotate_position(
                    &*chromosome_2,
                    position_2.saturating_sub(1),
                    gene_annotator
                );
                let position_annotation_2: PositionAnnotation = annotate_position(
                    &*chromosome_2,
                    position_2,
                    gene_annotator
                );

                let variant_type: VariantType = VariantType::from_str(col_variant_type.get(i).unwrap()).unwrap();
                let variant_sequence: Box<str> = col_variant_sequence.get(i).unwrap_or("").into();

                pb.inc(1);

                VariantCallAnnotation::new(
                    variant_call_id,
                    chromosome_1.into(),
                    position_1,
                    chromosome_2.into(),
                    position_2,
                    variant_type,
                    variant_sequence,
                    position_annotation_1,
                    position_annotation_1_plus_1,
                    position_annotation_2_minus_1,
                    position_annotation_2
                )
            })
            .collect()
    });
    
    pb.finish_with_message("Completed annotating variant calls.");

    variant_call_annotations
}


#[cfg(test)]
#[path = "../tests/pipeline/variant_annotation.rs"]
mod tests;