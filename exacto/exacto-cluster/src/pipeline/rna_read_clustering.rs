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
use exacto_caller::prelude::*;
use exacto_core::log_info;
use exacto_core::prelude::*;
use rayon::ThreadPool;
use std::collections::{HashMap, HashSet};
use std::env;
use std::ops::Range;
use std::path::{Path, PathBuf};
use noodles_bgzf::VirtualPosition;
use tempfile::{NamedTempFile, TempPath};
use exacto_qc::prelude::remove_unspliced_rnas;

use crate::clustering::grouping::{group_transcript_models_by_shared_junctions, pack_batches};
use crate::io::records::RNAReadClusterVariantRecord;
use crate::io::writers::write_rna_read_cluster_set;
use crate::options::ClusterRNAReadsOptions;
use crate::pipeline::batch::{cluster_batch, BatchContext};
use crate::read_cluster::rna_read_cluster_set::RNAReadClusterSet;
use crate::read_filtering::intron_retention::identify_nascent_rna_read_ids_per_strand;
use crate::read_filtering::soft_clip::build_intron_boundary_index;
use crate::reference::splice_junction_annotation_index::SpliceJunctionAnnotationIndex;
use crate::variant_filtering::known_variants::KnownVariants;


pub fn cluster_rna_reads(
    bam_file: &str,
    bai_file: &str,
    reference_genome_fasta_file: &str,
    gene_annotator: &(impl GeneAnnotator + Sync),
    dna_variant_records: Option<&Vec<DNAVariantRecord>>,
    allowed_variant_records: Option<&Vec<RNAReadClusterVariantRecord>>,
    options: &ClusterRNAReadsOptions,
    num_threads: usize,
    temp_dir: &str,
    output_dir: &str,
    output_prefix: &str
) -> RNAReadClusterSet {
    // Step 1. Check the inputs.
    check_alignment_inputs(
        bam_file,
        Some(bai_file),
        reference_genome_fasta_file
    );

    // Step 2. Fetch the temp directory.
    let temp_directory: PathBuf = fetch_temp_directory(temp_dir);

    // Step 3. Get chromosome IDs and names.
    log_info!("Getting chromosomes.");
    let chromosome_names_map: BiMap<ReferenceChromosomeName, ReferenceChromosomeID> = create_chromosome_names_map(bam_file);
    
    // Step 4. Put the DNA variant records.
    let known_variants: KnownVariants = KnownVariants::new(
        dna_variant_records,
        allowed_variant_records,
        &chromosome_names_map
    );

    // Step 5. Load the reference genome FASTA file.
    log_info!("Loading the reference genome FASTA file.");
    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);

    // Step 6. Prepare the RNA reads.
    let (read_names_map, models_file, summaries): (BiMap<Box<str>, usize>, TempPath, Vec<RNAReadCharacterizationSummary>) = prepare_reads(
        bam_file,
        bai_file,
        gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        options,
        num_threads,
        &temp_directory
    );

    // Step 6. Index the introns the reads splice. A terminal soft clip that spells the reference
    // across one of them is the end of a degraded read reaching into the neighbouring exon, not
    // an insertion.
    let intron_boundary_index: HashMap<(u16, u32), Vec<(u32, u32)>> = build_intron_boundary_index(&summaries);

    // Step 7. Identify groups of TranscriptModel objects to load together, and pack them into
    // batches.
    log_info!("Identifying groups of TranscriptModel objects to load together.");
    let annotation_index: SpliceJunctionAnnotationIndex = SpliceJunctionAnnotationIndex::new(
        gene_annotator,
        &chromosome_names_map
    );
    let summary_groups: Vec<Vec<&RNAReadCharacterizationSummary>> = group_transcript_models_by_shared_junctions(
        &summaries,
        &annotation_index,
        options.max_locus_gap,
        options.unspliced_bin_size,
        options.junction.min_reads,
        options.filtering.max_slippage_repeat_len,
        options.error_model.sequencing_error,
        options.error_model.slippage_prob,
        options.filtering.max_fpr,
        num_threads
    );
    log_info!("{} summary groups were identified.", summary_groups.len());
    let batches: Vec<Range<usize>> = pack_batches(&summary_groups, options.max_batch_reads);
    log_info!("{} batches to process.", batches.len());

    // Step 8. Build a read support index. Every depth it is asked for counts the reads of one
    // junction cluster: at a call's site, the cluster itself, or one of its cells. No cluster
    // holds more reads than its summary group. A site no read of the cluster holds has no floor
    // a call could clear (see `min_reads_at_site`).
    log_info!("Building minimum RNA read support index.");
    let max_group_size: u32 = summary_groups.iter().map(|group| group.len() as u32).max().unwrap_or(0);
    let min_read_support_index: RNAVariantReadSupportIndex = RNAVariantReadSupportIndex::new(
        &(1..=max_group_size).collect(),
        options.filtering.max_slippage_repeat_len,
        options.error_model.sequencing_error,
        options.error_model.slippage_prob,
        options.filtering.max_fpr,
        num_threads
    );

    // Step 9. Cluster each batch, then identify, genotype, and phase the variants within each
    // junction cluster.
    log_info!("Clustering reads by their splice junctions chains, then identify, genotype, and phase variants within each cluster.");
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let context: BatchContext<_> = BatchContext {
        gene_annotator,
        chromosome_names_map: &chromosome_names_map,
        fasta_map: &fasta_map,
        annotation_index: &annotation_index,
        intron_boundary_index: &intron_boundary_index,
        known_variants: &known_variants,
        read_support_index: &min_read_support_index,
        models_file: &models_file,
        options,
        thread_pool: &thread_pool
    };
    let mut rna_read_cluster_set: RNAReadClusterSet = RNAReadClusterSet::new(
        read_names_map,
        chromosome_names_map.clone()
    );
    let (mut next_cluster_id, mut next_variant_call_id): (usize, usize) = (1, 1);
    let num_batches: usize = batches.len();
    for (i, batch) in batches.into_iter().enumerate() {
        cluster_batch(&context, &summary_groups[batch], &mut next_cluster_id, &mut next_variant_call_id, &mut rna_read_cluster_set);
        if (i + 1) * 10 / num_batches > i * 10 / num_batches {
            log_info!(
                "Processed {} of {} batches ({}%).",
                i + 1,
                num_batches,
                (i + 1) * 100 / num_batches
            );
        }
    }

    // Step 10. Write the tables.
    log_info!("Writing the clustering results to {}.", output_dir);
    write_rna_read_cluster_set(&rna_read_cluster_set, output_dir, output_prefix)
        .unwrap_or_else(|error| panic!("Failed to write the RNA read cluster tables to {}: {}", output_dir, error));

    rna_read_cluster_set
}


fn prepare_reads(
    bam_file: &str,
    bam_bai_file: &str,
    gene_annotator: &(impl GeneAnnotator + Sync),
    chromosome_names_map: &BiMap<Box<str>, u16>,
    fasta_map: &FastaMap,
    options: &ClusterRNAReadsOptions,
    num_threads: usize,
    temp_dir: &PathBuf
) -> (BiMap<Box<str>, usize>, TempPath, Vec<RNAReadCharacterizationSummary>) {
    // Step 1. Remove unspliced RNA reads.
    let unspliced_removed_bam_temp_file: NamedTempFile = NamedTempFile::new_in(temp_dir).unwrap();
    let unspliced_removed_bam_file: &str = unspliced_removed_bam_temp_file.path().to_str().expect("Path is not valid UTF-8.");
    let bam_file: &str = if options.remove_unspliced_rnas {
        log_info!("Removing unspliced RNA reads.");
        let unspliced_removed_bai_file: String = format!("{}.bai", unspliced_removed_bam_file);
        remove_unspliced_rnas(
            bam_file,
            bam_bai_file,
            gene_annotator,
            unspliced_removed_bam_file,
            unspliced_removed_bai_file.as_str(),
            num_threads,
            options.calling.min_mapping_quality,
            true
        );
        unspliced_removed_bam_file
    } else {
        bam_file
    };

    // Step 2. Index the BAM file.
    log_info!("Indexing the BAM file.");
    let (mut record_positions_map, mut read_names_map): (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(
        bam_file,
        true,
        num_threads
    );
    let num_reads: usize = record_positions_map.len();

    // Step 3. Drop reads with too many alignment records.
    log_info!("Dropping reads with more than {} alignment records.", options.calling.max_records);
    record_positions_map.retain(|_read_id, virtual_positions|
        virtual_positions.len() <= options.calling.max_records
    );
    read_names_map.retain(|_read_name, read_id|
        record_positions_map.contains_key(read_id)
    );
    let num_excluded: usize = num_reads - record_positions_map.len();
    log_info!(
        "Excluded {} read(s) with more than {} alignment record(s).",
        num_excluded,
        options.calling.max_records
    );

    // Step 4. Get read IDs and sort their order.
    log_info!("Sorting read IDs.");
    let mut read_ids: Vec<usize> = read_names_map.right_values().copied().collect();
    read_ids.sort_unstable();

    // Step 5. Characterize each read's junctions and variants.
    log_info!("Characterizing RNA reads.");
    let (temp_file, mut summaries, _failed_read_ids): (TempPath, Vec<RNAReadCharacterizationSummary>, Vec<usize>) = characterize_rna_reads(
        bam_file,
        &read_ids,
        &read_names_map,
        &record_positions_map,
        chromosome_names_map,
        fasta_map,
        gene_annotator,
        num_threads,
        options.calling.min_mapping_quality,
        options.calling.min_terminal_soft_clip_ins_len,
        options.calling.bkpt_rescue,
        options.calling.bkpt_rescue_min_ins_len,
        options.calling.bkpt_rescue_realignment_gap_open_score,
        options.calling.bkpt_rescue_realignment_gap_extend_score,
        options.calling.bkpt_rescue_realignment_k,
        options.calling.bkpt_rescue_realignment_band_width,
        options.calling.bkpt_rescue_realignment_min_score_fraction,
        options.calling.bkpt_rescue_realignment_min_query_coverage,
        options.calling.bkpt_rescue_realignment_min_placed_fraction,
        options.calling.bkpt_rescue_realignment_max_pieces,
        None,
        options.calling.chunk_size,
        temp_dir.to_str().unwrap()
    );

    // Step 6. Remove reads with unsupported intron retention.
    let nascent_read_ids: HashSet<usize> = if options.remove_unspliced_rnas {
        identify_nascent_rna_read_ids_per_strand(
            &mut summaries,
            chromosome_names_map,
            fasta_map,
            options.error_model.sequencing_error,
            options.filtering.max_fpr,
            num_threads
        )
    } else {
        HashSet::new()
    };

    // The exon chains have no reader past this point, so they are released with the reads.
    let mut summaries: Vec<RNAReadCharacterizationSummary> = summaries
        .into_iter()
        .filter(|summary| !nascent_read_ids.contains(&summary.read_id))
        .map(|summary| RNAReadCharacterizationSummary { exons: Vec::new(), ..summary })
        .collect();
    summaries.sort_by_key(|summary| summary.read_id);

    (read_names_map, temp_file, summaries)
}


#[cfg(test)]
#[path = "../tests/pipeline/rna_read_clustering.rs"]
mod tests;
