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
use bincode;
use exacto_core::prelude::*;
use exacto_core::log_info;
use noodles_bam as bam;
use noodles_bam::bai;
use noodles_bam::bai::Index;
use noodles_bgzf::VirtualPosition;
use noodles_sam::Header;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;
use tempfile::{NamedTempFile, TempPath};

use crate::prelude::*;
use crate::pipeline::regions::{generate_genomic_regions, get_region_ends};
use crate::pipeline::options::IdentifyGermlineDNAVariantsOptions;
use crate::calling::dna::variant_record_clustering::{cluster_dna_variant_records, count_dna_variant_depths};
use crate::filtering::dna::strand_bias_filter::StrandBiasFilter;
use crate::filtering::dna::variant_allele_fraction_filter::VariantAlleleFractionFilter;
use crate::filtering::dna::variant_read_support_filter::DNAVariantReadSupportFilter;
use crate::filtering::dna::variant_read_support_index::DNAVariantReadSupportIndex;
use crate::variant::dna::variant_call_set::merge_temp_variant_call_sets;


pub fn identify_germline_dna_variants(
    bam_file: &str,
    bai_file: &str,
    reference_genome_fasta_file: &str,
    regions: &Vec<(&str, ReferencePosition, ReferencePosition)>,
    options: &IdentifyGermlineDNAVariantsOptions,
    num_threads: usize,
    temp_dir: &str
) -> DNAVariantCallSet {
    // Step 1. Check the inputs.
    check_alignment_inputs(
        bam_file,
        Some(bai_file),
        reference_genome_fasta_file
    );

    // Step 2. Get chromosome IDs and names.
    log_info!("Getting chromosomes.");
    let chromosome_names_map: BiMap<ReferenceChromosomeName, ReferenceChromosomeID> = create_chromosome_names_map(bam_file);
    let chromosome_lengths: HashMap<ReferenceChromosomeName, u32> = get_chromosome_lengths(bam_file);
    let chromosomes: Vec<ReferenceChromosomeName> = chromosome_lengths
        .keys()
        .cloned()
        .collect();

    // Step 3. Get sequencing depths.
    log_info!("Getting the greatest sequencing depth.");
    let max_depth: ReadDepth = get_bam_max_depth(bam_file, bai_file, num_threads);
    log_info!("\tMax sequencing depth: {max_depth}.");

    // A file without aligned reads has nothing to call.
    if max_depth == 0 {
        let mut variant_call_set: DNAVariantCallSet = DNAVariantCallSet::new(DNAVariantOrigin::Germline);
        variant_call_set.load_chromosome_names(chromosome_names_map);
        return variant_call_set;
    }

    // Step 4. Compute minimum read support index.
    log_info!("Computing minimum read support index.");
    let min_read_support_index: DNAVariantReadSupportIndex = DNAVariantReadSupportIndex::new(
        max_depth,
        options.filtering.max_slippage_repeat_len,
        options.prior.allele_fraction,
        options.prior.mutation_rate,
        options.error_model.sequencing_error,
        options.error_model.slippage_prob,
        options.filtering.max_fpr
    );

    // Step 5. Index BAM records by read names.
    log_info!("Indexing BAM records by read names.");
    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        num_threads
    );

    // Step 6. Generate regions.
    log_info!("Generating genomic regions.");
    let ordered_regions: BTreeMap<ReferenceChromosomeName, Vec<(ReferencePosition, ReferencePosition)>> = generate_genomic_regions(
        regions,
        &chromosomes,
        &chromosome_lengths,
        options.calling.chunk_size
    );

    // Step 7. Load the reference genome FASTA file.
    log_info!("Loading the reference genome FASTA file.");
    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);

    // Step 8. Fetch the temp directory.
    let temp_directory: PathBuf = fetch_temp_directory(temp_dir);

    // Step 9. Identify DNA variants by region.
    log_info!("Identifying DNA variants by region.");
    let temp_files: Vec<TempPath> = identify_germline_dna_variants_by_region(
        bam_file,
        bai_file,
        &record_positions_map,
        &read_names_map,
        &fasta_map,
        &ordered_regions,
        &chromosome_names_map,
        &min_read_support_index,
        options,
        &temp_directory,
        num_threads
    );

    // Step 10. Load VariantCallSet objects and merge them.
    log_info!("Loading temp files and merging them into a variant call set.");
    let variant_call_set: DNAVariantCallSet = merge_temp_variant_call_sets(
        DNAVariantOrigin::Germline,
        &temp_files,
        &read_names_map,
        &chromosome_names_map,
        num_threads
    );

    variant_call_set
}


fn identify_germline_dna_variants_by_region(
    bam_file: &str,
    bai_file: &str,
    record_positions_map: &HashMap<ReadID, Vec<VirtualPosition>>,
    read_names_map: &BiMap<ReadName, ReadID>,
    fasta_map: &FastaMap,
    ordered_regions: &BTreeMap<ReferenceChromosomeName, Vec<(ReferencePosition, ReferencePosition)>>,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    min_read_support_index: &DNAVariantReadSupportIndex,
    options: &IdentifyGermlineDNAVariantsOptions,
    temp_directory: &PathBuf,
    num_threads: usize
) -> Vec<TempPath> {
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let mut variant_call_idx: VariantID = 1;
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bai_file).unwrap();
    let mut temp_files: Vec<TempPath> = Vec::new();
    for (chromosome, curr_regions) in ordered_regions {
        let chromosome_id: ReferenceChromosomeID = chromosome_names_map.get_by_left(chromosome).unwrap().clone();
        let region_ends: Vec<ReferencePosition> = get_region_ends(curr_regions);
        let mut chromosome_variant_records: Vec<VariantRecord> = Vec::new();
        let mut prev_end: ReferencePosition = 0;
        for (chunk_index, (start, end)) in curr_regions.iter().enumerate() {
            // A window owns the calls that start in its chunk, (prev_end, end] (the chromosome's
            // first window, [0, end]). It reads max_clustering_distance past the chunk, so that
            // it sees whole a call that starts in the chunk and ends past it, but fetches no read
            // past the end of the chunk's region.
            let window_end: ReferencePosition = end.saturating_add(options.calling.max_clustering_distance);
            let owned_from: ReferencePosition = if prev_end == 0 { 0 } else { prev_end + 1 };

            // Step 1. Fetch BAM records.
            let records_map: HashMap<ReadID, Vec<bam::Record>> = fetch_bam_records(
                &mut reader,
                &header,
                &index,
                chromosome.clone(),
                *start,
                window_end.min(region_ends[chunk_index]),
                &record_positions_map,
                &read_names_map,
                options.calling.max_records,
                num_threads
            );

            if records_map.is_empty() {
                continue;
            } else {
                log_info!("\tIdentifying variant calls in {}:{}-{}.", chromosome, start, end);
                log_info!("\t\t{} BAM record(s) fetched.", records_map.len());
            }

            // Step 2. Identify variant records in each read.
            let variant_caller: DNAVariantRecordCaller = DNAVariantRecordCaller::new(
                options.calling.min_mapping_quality,
                options.calling.min_base_quality,
                options.calling.min_terminal_soft_clip_ins_len
            );
            let mut variant_records: HashSet<VariantRecord> = thread_pool.install(|| {
                records_map
                    .into_par_iter()
                    .flat_map(|(read_id, mut records)| {
                        // Duplicate and QC-failed records are left out, as the read depths leave
                        // them out; a read whose primary record is one of them is left out whole.
                        records.retain(|record| !record.flags().is_duplicate() && !record.flags().is_qc_fail());
                        if records.iter().all(|record| record.flags().is_supplementary()) {
                            return Vec::new();
                        }

                        // Get the FASTX read sequence.
                        let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);

                        // Get the FASTX base quality scores.
                        let base_quality_scores: Vec<BaseQuality> = get_bam_fastx_base_quality_scores(&records);

                        // Get the BAM records.
                        let bam_records: Vec<Arc<bam::Record>> = records
                            .into_iter()
                            .map(|r| Arc::new(r))
                            .collect();


                        // Construct an instance of AlignmentModel.
                        let alignment_model: AlignmentModel = AlignmentModel::new(
                            read_id,
                            &*read_sequence,
                            &base_quality_scores,
                            &bam_records
                        );

                        // Filter based on mapping quality.
                        let max_mapping_quality = alignment_model.get_records().iter()
                            .map(|r| get_alignment_mapping_quality(&r.record))
                            .max()
                            .unwrap_or(0);
                        if options.calling.min_mapping_quality > max_mapping_quality {
                            return Vec::new();
                        }

                        // Identify DNA variants.
                        variant_caller.call(&alignment_model)
                    })
                    .collect()
            });

            // Filter variant records by chromosome. The windows below are taken by position 1,
            // so a translocation belongs to the chromosome of its side 1 and to no other.
            variant_records = thread_pool.install(|| {
                variant_records
                    .into_par_iter()
                    .filter(|vr| vr.get_chromosome_1() == chromosome_id)
                    .collect()
            });
            log_info!("\t\t{} variant record(s) identified.", variant_records.len());

            // Move into chromosome variant records
            chromosome_variant_records.extend(variant_records.into_iter());

            // Keep sorted by position_1
            chromosome_variant_records.sort_unstable_by_key(|vr| vr.get_position_1());

            // Compute window indices [keep_from, window_end]
            let keep_from: ReferencePosition = prev_end.saturating_sub(options.calling.max_clustering_distance);
            let keep_to: ReferencePosition = window_end;

            // Step 3. Slice variant records for clustering.
            let start_idx: usize = chromosome_variant_records.partition_point(|vr| vr.get_position_1() < keep_from);
            let end_idx: usize = chromosome_variant_records.partition_point(|vr| vr.get_position_1() <= keep_to);
            let variant_records_slice: &[VariantRecord] = &chromosome_variant_records[start_idx..end_idx];
            let variant_records_slice_rc: Vec<Arc<VariantRecord>> = variant_records_slice
                .iter()
                .cloned()
                .map(Arc::new)
                .collect();

            // Step 4. Cluster variant records into the variant calls this window owns: those
            // that start in its chunk.
            let mut variant_calls: Vec<VariantCall> = cluster_dna_variant_records(
                variant_records_slice_rc,
                owned_from..=*end,
                &fasta_map,
                &chromosome_names_map,
                num_threads,
                options.calling.min_size_proportion,
                options.calling.max_ins_norm_edit_distance,
                options.calling.max_clustering_distance,
                options.error_model.sequencing_error,
                options.calling.bkpt_rescue,
                options.calling.bkpt_rescue_min_ins_len,
                options.calling.bkpt_rescue_max_ins_len,
                options.calling.bkpt_rescue_search_distance,
                options.calling.bkpt_rescue_realignment_gap_open_score,
                options.calling.bkpt_rescue_realignment_gap_extend_score,
                options.calling.bkpt_rescue_realignment_k,
                options.calling.bkpt_rescue_realignment_band_width,
                options.calling.bkpt_rescue_realignment_min_score_fraction,
                options.calling.bkpt_rescue_realignment_min_query_coverage,
                options.calling.bkpt_rescue_realignment_min_span_proportion,
                options.calling.poa_match_score,
                options.calling.poa_mismatch_score,
                options.calling.poa_gap_open_score,
                options.calling.poa_gap_extend_score
            );

            // Count depth and strands at the positions the calls read, and set their total depth.
            let read_depths: BAMReadDepths = count_dna_variant_depths(
                &mut variant_calls,
                bam_file,
                bai_file,
                &chromosome_names_map,
                options.calling.read_depth_max_merge_distance
            );

            // Step 5. Filter variant calls based on the following:
            // 1. Read support: total depth, number of reads, and the minimum read support at the
            //    depth of each flank (a flank no read covers fails)
            // 2. Alternate allele fraction
            // 3. Strand bias
            // Total depth comes first: the allele fraction needs a depth above 0.
            let filters: Vec<Box<dyn VariantFilter<Input = VariantCall> + '_>> = vec![
                Box::new(DNAVariantReadSupportFilter::new(
                    options.filtering.min_reads as ReadSupport,
                    options.filtering.min_total_depth as ReadSupport,
                    options.filtering.min_homopolymer_len,
                    options.filtering.min_dinucleotide_context_len,
                    options.filtering.max_slippage_repeat_len,
                    None,
                    &min_read_support_index,
                    &read_depths,
                    &chromosome_names_map,
                    &fasta_map
                )),
                Box::new(VariantAlleleFractionFilter::new(options.filtering.min_alt_allele_fraction)),
                Box::new(StrandBiasFilter::new(&read_depths, &chromosome_names_map, 0.05))
            ];
            variant_calls = thread_pool.install(|| {
                variant_calls
                    .into_par_iter()
                    .filter(|vc| filters.iter().all(|filter| filter.passes(vc)))
                    .collect()
            });
            log_info!("\t\t{} variant call(s) identified.", variant_calls.len());

            // Step 6. Aggregate variant calls into a variant callset.
            let mut variant_call_set = DNAVariantCallSet::new(DNAVariantOrigin::Germline);
            for mut variant_call in variant_calls {
                // Rename the variant call ID
                variant_call.set_id(variant_call_idx);
                variant_call_set.add_variant_call(variant_call);
                variant_call_idx += 1;
            }

            // Store variant_call_set in a temp file.
            let temp_path: TempPath = {
                let temp_file = NamedTempFile::new_in(temp_directory.as_path()).unwrap();
                let mut writer = BufWriter::new(temp_file);
                bincode::serialize_into(&mut writer, &variant_call_set)
                    .expect("Failed to serialize variant_call_set");
                writer.flush().unwrap();
                let temp_file: NamedTempFile = writer.into_inner().unwrap();
                temp_file.into_temp_path()
            };
            temp_files.push(temp_path);

            prev_end = *end;
        }
        chromosome_variant_records.clear();
    }

    temp_files
}


#[cfg(test)]
#[path = "../tests/pipeline/variant_calling_dna_germline.rs"]
mod tests;