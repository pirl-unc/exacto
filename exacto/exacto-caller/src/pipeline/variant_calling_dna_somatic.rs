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
use exacto_core::log_info;
use exacto_core::prelude::*;
use noodles_bam as bam;
use noodles_bam::bai;
use noodles_bam::bai::Index;
use noodles_bgzf::VirtualPosition;
use noodles_sam::Header;
use rayon::iter::IntoParallelIterator;
use rayon::prelude::*;
use rayon::ThreadPool;
use std::cmp::max;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufWriter, Write};
use std::ops::RangeInclusive;
use std::path::PathBuf;
use std::sync::Arc;
use tempfile::{NamedTempFile, TempPath};

use crate::prelude::*;
use crate::pipeline::options::IdentifySomaticDNAVariantsOptions;
use crate::pipeline::regions::{generate_genomic_regions, get_region_ends};
use crate::calling::dna::breakpoint_rescue::retype_insertion_as_breakpoints;
use crate::calling::dna::variant_record_clustering::{build_dna_breakpoint_interval_trees, cluster_dna_variant_records, count_dna_variant_depths};
use crate::filtering::dna::control_subtraction::diff_dna_variant_records;
use crate::filtering::dna::strand_bias_filter::StrandBiasFilter;
use crate::filtering::dna::variant_allele_fraction_filter::VariantAlleleFractionFilter;
use crate::filtering::dna::variant_read_support_filter::DNAVariantReadSupportFilter;
use crate::filtering::dna::variant_read_support_index::DNAVariantReadSupportIndex;
use crate::variant::dna::variant_call_set::{merge_temp_variant_call_sets, DNAVariantCallSet};


pub fn identify_somatic_dna_variants(
    case_bam_file: &str,
    case_bai_file: &str,
    control_bam_files: Vec<&str>,
    control_bai_files: Vec<&str>,
    reference_genome_fasta_file: &str,
    regions: &Vec<(&str, ReferencePosition, ReferencePosition)>,
    options: &IdentifySomaticDNAVariantsOptions,
    num_threads: usize,
    temp_dir: &str
) -> DNAVariantCallSet {
    // Step 1. Check the inputs.
    check_alignment_inputs(
        case_bam_file,
        Some(case_bai_file),
        reference_genome_fasta_file
    );
    assert_eq!(
        control_bam_files.len(),
        control_bai_files.len(),
        "control_bam_files and control_bai_files must have the same number of files."
    );
    for (control_bam_file, control_bam_bai_file) in control_bam_files.iter().zip(control_bai_files.iter()) {
        check_alignment_inputs(
            control_bam_file,
            Some(*control_bam_bai_file),
            reference_genome_fasta_file
        );
    }

    // Step 2. Get chromosome IDs and names.
    let chromosome_names_map: BiMap<ReferenceChromosomeName, ReferenceChromosomeID> = create_chromosome_names_map(case_bam_file);
    let chromosome_lengths: HashMap<ReferenceChromosomeName, u32> = get_chromosome_lengths(case_bam_file);
    let chromosomes: Vec<ReferenceChromosomeName> = chromosome_lengths
        .keys()
        .cloned()
        .collect();

    // Step 3. Make sure chromosome IDs and names are the same in the control BAM files.
    for control_bam_file in control_bam_files.iter() {
        let chromosome_names_map_: BiMap<ReferenceChromosomeName, ReferenceChromosomeID> = create_chromosome_names_map(control_bam_file);
        let chromosome_lengths_: HashMap<ReferenceChromosomeName, u32> = get_chromosome_lengths(control_bam_file);
        for (chromosome_name,chromosome_id) in chromosome_names_map.iter() {
            if chromosome_names_map_.get_by_left(chromosome_name) != Some(chromosome_id) {
                panic!("The contig names of {control_bam_file} differ from those of {case_bam_file}.");
            }
        }
        if chromosome_lengths != chromosome_lengths_ {
            panic!("The contig lengths of {control_bam_file} differ from those of {case_bam_file}.");
        }
    }

    // Step 4. Get sequencing depths.
    log_info!("Getting the greatest sequencing depth.");
    let max_depth: ReadDepth = get_bam_max_depth(case_bam_file, case_bai_file, num_threads);
    log_info!("\tMax sequencing depth: {max_depth}.");

    // A case file without aligned reads has nothing to call.
    if max_depth == 0 {
        let mut variant_call_set: DNAVariantCallSet = DNAVariantCallSet::new(DNAVariantOrigin::Somatic);
        variant_call_set.load_chromosome_names(chromosome_names_map);
        return variant_call_set;
    }

    // Step 5. Compute minimum read support index.
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

    // Step 6. Index BAM records by read names.
    log_info!("Indexing BAM records by read names.");
    let (case_record_positions_map, case_read_names_map) = index_bam_records(
        case_bam_file,
        true,
        num_threads
    );
    log_info!("{} read names/IDs in case BAM file.", case_read_names_map.len());
    let mut control_record_positions_map: HashMap<Box<str>, HashMap<ReadID, Vec<VirtualPosition>>> = HashMap::new();
    let mut control_read_names_map: HashMap<Box<str>, BiMap<ReadName, ReadID>> = HashMap::new();
    for (index,control_bam_file) in control_bam_files.iter().enumerate() {
        let (record_positions_map, read_names_map) = index_bam_records(
            control_bam_file,
            true,
            num_threads
        );
        log_info!("{} read names/IDs in control BAM file {}/{}.", read_names_map.len(), index + 1, control_bam_files.len());
        control_record_positions_map.insert(control_bam_file.to_string().into_boxed_str(), record_positions_map);
        control_read_names_map.insert(control_bam_file.to_string().into_boxed_str(), read_names_map);
    }

    // Step 7. Generate regions.
    log_info!("Generating genomic regions.");
    let ordered_regions: BTreeMap<ReferenceChromosomeName, Vec<(ReferencePosition, ReferencePosition)>> = generate_genomic_regions(
        regions,
        &chromosomes,
        &chromosome_lengths,
        options.calling.chunk_size
    );

    // Step 8. Load the reference genome FASTA file.
    log_info!("Loading the reference genome FASTA file.");
    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);

    // Step 9. Fetch the temp directory.
    let temp_directory: PathBuf = fetch_temp_directory(temp_dir);

    // Step 10. Identify case-specific variant calls.
    log_info!("Identifying case-specific DNA variants.");
    let temp_files: Vec<TempPath> = identify_somatic_dna_variants_by_region(
        case_bam_file,
        case_bai_file,
        &control_bam_files,
        &control_bai_files,
        &case_record_positions_map,
        &case_read_names_map,
        &control_record_positions_map,
        &control_read_names_map,
        &fasta_map,
        &ordered_regions,
        &chromosome_names_map,
        &min_read_support_index,
        options,
        max_depth,
        &temp_directory,
        num_threads
    );

    // Step 11. Load all VariantCallSet objects and merge them.
    log_info!("Loading all temp files and merging them into a variant call set.");
    let variant_call_set: DNAVariantCallSet = merge_temp_variant_call_sets(
        DNAVariantOrigin::Somatic,
        &temp_files,
        &case_read_names_map,
        &chromosome_names_map,
        num_threads
    );

    variant_call_set
}


/// Call the somatic variants of each region chunk, and store each chunk's calls in a temp file.
fn identify_somatic_dna_variants_by_region(
    case_bam_file: &str,
    case_bai_file: &str,
    control_bam_files: &Vec<&str>,
    control_bai_files: &Vec<&str>,
    case_record_positions_map: &HashMap<ReadID, Vec<VirtualPosition>>,
    case_read_names_map: &BiMap<ReadName, ReadID>,
    control_record_positions_map: &HashMap<Box<str>, HashMap<ReadID, Vec<VirtualPosition>>>,
    control_read_names_map: &HashMap<Box<str>, BiMap<ReadName, ReadID>>,
    fasta_map: &FastaMap,
    ordered_regions: &BTreeMap<ReferenceChromosomeName, Vec<(ReferencePosition, ReferencePosition)>>,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    min_read_support_index: &DNAVariantReadSupportIndex,
    options: &IdentifySomaticDNAVariantsOptions,
    max_depth: ReadDepth,
    temp_directory: &PathBuf,
    num_threads: usize
) -> Vec<TempPath> {
    let max_clustering_distance: u32 = options.calling.max_clustering_distance;
    let bin_size: u32 = 10_u32.pow((max_clustering_distance as f32).log10().floor() as u32 + 1);
    let mut case_reader = bam::io::reader::Builder::default()
        .build_from_path(case_bam_file)
        .unwrap();
    let case_header: Header = case_reader.read_header().unwrap();
    let case_index: Index = bai::fs::read(case_bai_file).unwrap();

    // The control files are read in the order given.
    let mut controls: Vec<(&str, bam::io::Reader<_>, Header, Index)> = Vec::new();
    for (control_bam_file, control_bai_file) in control_bam_files.iter().zip(control_bai_files.iter()) {
        let mut reader = bam::io::reader::Builder::default()
            .build_from_path(control_bam_file)
            .unwrap();
        let header: Header = reader.read_header().unwrap();
        let index: Index = bai::fs::read(control_bai_file).unwrap();
        controls.push((*control_bam_file, reader, header, index));
    }

    let mut variant_call_idx: VariantID = 1;
    let mut temp_files: Vec<TempPath> = Vec::new();
    for (chromosome, curr_regions) in ordered_regions {
        let chromosome_id: ReferenceChromosomeID = chromosome_names_map.get_by_left(chromosome).unwrap().clone();
        let region_ends: Vec<ReferencePosition> = get_region_ends(curr_regions);
        // The chromosome's case variant records so far, sorted by position 1.
        let mut chromosome_case_variant_records: Vec<VariantRecord> = Vec::new();
        let mut prev_end: ReferencePosition = 0;
        for (chunk_index, (start, end)) in curr_regions.iter().enumerate() {
            // A window owns the calls that start in its chunk, (prev_end, end] (the chromosome's
            // first window, [0, end]). It reads max_clustering_distance past the chunk, so that
            // it sees whole a call that starts in the chunk and ends past it, but fetches no case
            // read past the end of the chunk's region.
            let window_end: ReferencePosition = end.saturating_add(max_clustering_distance);
            let owned_from: ReferencePosition = if prev_end == 0 { 0 } else { prev_end + 1 };

            // Step 1. Fetch case BAM records.
            let case_records_map: HashMap<ReadID, Vec<bam::Record>> = fetch_bam_records(
                &mut case_reader,
                &case_header,
                &case_index,
                chromosome.clone(),
                *start,
                window_end.min(region_ends[chunk_index]),
                case_record_positions_map,
                case_read_names_map,
                options.calling.max_records,
                num_threads
            );
            if case_records_map.is_empty() {
                continue;
            }
            log_info!("\tProcessing {}:{}-{}.", chromosome, start, end);
            log_info!("\t\t{} BAM record(s) fetched.", case_records_map.len());

            // Step 2. Identify case variant records in each read.
            let case_variant_records: HashSet<VariantRecord> = identify_dna_variant_records(
                case_records_map,
                chromosome_id,
                &options.calling,
                num_threads
            );
            log_info!("\t\t{} case variant record(s) identified.", case_variant_records.len());
            chromosome_case_variant_records.extend(case_variant_records);
            chromosome_case_variant_records.sort_unstable_by_key(|vr| vr.get_position_1());

            // Step 3. Slice the case variant records of the window [keep_from, window_end].
            let keep_from: ReferencePosition = prev_end.saturating_sub(max_clustering_distance);
            let keep_to: ReferencePosition = window_end;
            let start_idx: usize = chromosome_case_variant_records.partition_point(|vr| vr.get_position_1() < keep_from);
            let end_idx: usize = chromosome_case_variant_records.partition_point(|vr| vr.get_position_1() <= keep_to);
            let window_case_variant_records: Vec<Arc<VariantRecord>> = chromosome_case_variant_records[start_idx..end_idx]
                .iter()
                .cloned()
                .map(Arc::new)
                .collect();

            // Step 4. Cluster case variant records into the variant calls this window owns, those
            // that start in its chunk, and filter them.
            let variant_calls: Vec<VariantCall> = call_dna_variants(
                window_case_variant_records,
                owned_from..=*end,
                case_bam_file,
                case_bai_file,
                fasta_map,
                chromosome_names_map,
                min_read_support_index,
                max_depth,
                options,
                num_threads
            );
            log_info!("\t\t{} case variant call(s) identified.", variant_calls.len());
            let mut case_variant_records: Vec<Arc<VariantRecord>> = variant_calls
                .into_iter()
                .flat_map(VariantCall::into_variant_records)
                .map(Arc::new)
                .collect();
            log_info!("\t\t{} case variant record(s) retained after filtering.", case_variant_records.len());

            // Step 5. Filter out case variant records near control variant records.
            for (i, (control_bam_file, control_reader, control_header, control_index)) in controls.iter_mut().enumerate() {
                if case_variant_records.is_empty() {
                    break;
                }
                log_info!("\t\tProcessing control BAM file {}/{}.", i + 1, control_bam_files.len());

                // Fetch control BAM records around the window.
                let control_records_map: HashMap<ReadID, Vec<bam::Record>> = fetch_bam_records(
                    control_reader,
                    control_header,
                    control_index,
                    chromosome.clone(),
                    max(1u32, start.saturating_sub(max_clustering_distance)),
                    window_end.saturating_add(max_clustering_distance),
                    &control_record_positions_map[*control_bam_file],
                    &control_read_names_map[*control_bam_file],
                    options.calling.max_records,
                    num_threads
                );
                log_info!("\t\t\t{} BAM record(s) fetched.", control_records_map.len());

                // Identify control variant records, and give its insertions their breakpoint form as well.
                let control_variant_records: HashSet<VariantRecord> = identify_dna_variant_records(
                    control_records_map,
                    chromosome_id,
                    &options.calling,
                    num_threads
                );
                log_info!("\t\t\t{} control variant record(s) identified.", control_variant_records.len());
                let control_variant_records: Vec<Arc<VariantRecord>> = retype_control_insertions(
                    control_variant_records,
                    &case_variant_records,
                    fasta_map,
                    chromosome_names_map,
                    &options.calling
                );

                // Filter out control variants.
                case_variant_records = diff_dna_variant_records(
                    case_variant_records,
                    control_variant_records,
                    bin_size,
                    num_threads,
                    options.calling.min_size_proportion,
                    options.calling.max_ins_norm_edit_distance,
                    max_clustering_distance,
                    options.error_model.sequencing_error,
                    options.subtraction.max_control_reads.saturating_add(1),
                    options.subtraction.apply_infinite_sites_assumption,
                    false
                );
                log_info!("\t\t\t{} case variant record(s) retained after diffing.", case_variant_records.len());
            }

            // Step 6. Cluster case-specific variant records into variant calls, and filter them.
            // They all belong to calls of this window.
            let variant_calls: Vec<VariantCall> = call_dna_variants(
                case_variant_records,
                0..=ReferencePosition::MAX,
                case_bam_file,
                case_bai_file,
                fasta_map,
                chromosome_names_map,
                min_read_support_index,
                max_depth,
                options,
                num_threads
            );
            log_info!("\t\t\t{} case-specific variant call(s) identified.", variant_calls.len());

            // Step 7. Number the variant calls, and store them in a temp file.
            let mut variant_call_set: DNAVariantCallSet = DNAVariantCallSet::new(DNAVariantOrigin::Somatic);
            for mut variant_call in variant_calls {
                variant_call.set_id(variant_call_idx);
                variant_call_set.add_variant_call(variant_call);
                variant_call_idx += 1;
            }
            let temp_path: TempPath = {
                let temp_file = NamedTempFile::new_in(temp_directory).unwrap();
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
    }

    temp_files
}


fn identify_dna_variant_records(
    records_map: HashMap<ReadID, Vec<bam::Record>>,
    chromosome_id: ReferenceChromosomeID,
    options: &DNAVariantCallingOptions,
    num_threads: usize
) -> HashSet<VariantRecord> {
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    // min_terminal_soft_clip_ins_len 0: DNA clip evidence belongs to breakpoint clustering.
    let variant_caller: DNAVariantRecordCaller = DNAVariantRecordCaller::new(
        options.min_mapping_quality,
        options.min_base_quality,
        0
    );
    thread_pool.install(|| {
        records_map
            .into_par_iter()
            .flat_map(|(read_id, mut records)| {
                // Duplicate and QC-failed records are left out, as the read depths leave them
                // out; a read whose primary record is one of them is left out whole.
                records.retain(|record| !record.flags().is_duplicate() && !record.flags().is_qc_fail());
                if records.iter().all(|record| record.flags().is_supplementary()) {
                    return Vec::new();
                }

                // Get the FASTX read sequence and base quality scores.
                let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
                let base_quality_scores: Vec<BaseQuality> = get_bam_fastx_base_quality_scores(&records);

                // Construct an instance of AlignmentModel.
                let bam_records: Vec<Arc<bam::Record>> = records
                    .into_iter()
                    .map(Arc::new)
                    .collect();
                let alignment_model: AlignmentModel = AlignmentModel::new(
                    read_id,
                    &read_sequence,
                    &base_quality_scores,
                    &bam_records
                );

                // Filter based on mapping quality.
                let max_mapping_quality: MappingQuality = alignment_model.get_records()
                    .iter()
                    .map(|r| get_alignment_mapping_quality(&r.record))
                    .max()
                    .unwrap_or(0);
                if options.min_mapping_quality > max_mapping_quality {
                    return Vec::new();
                }

                // Identify DNA variants.
                variant_caller.call(&alignment_model)
            })
            .filter(|vr| vr.get_chromosome_1() == chromosome_id)
            .collect()
    })
}


fn call_dna_variants(
    variant_records: Vec<Arc<VariantRecord>>,
    owned_positions: RangeInclusive<ReferencePosition>,
    case_bam_file: &str,
    case_bai_file: &str,
    fasta_map: &FastaMap,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    min_read_support_index: &DNAVariantReadSupportIndex,
    max_depth: ReadDepth,
    options: &IdentifySomaticDNAVariantsOptions,
    num_threads: usize
) -> Vec<VariantCall> {
    // Step 1. Cluster variant records into variant calls.
    let mut variant_calls: Vec<VariantCall> = cluster_dna_variant_records(
        variant_records,
        owned_positions,
        fasta_map,
        chromosome_names_map,
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

    // Step 2. Count depth and strands at the positions the calls read, and set their total depth.
    let read_depths: BAMReadDepths = count_dna_variant_depths(
        &mut variant_calls,
        case_bam_file,
        case_bai_file,
        chromosome_names_map,
        options.calling.read_depth_max_merge_distance
    );

    // Step 3. Filter variant calls for the following:
    // 1. Read support: total depth, number of reads, and the minimum read support at the depth of
    //    each flank (a flank no read covers is held to the support at the greatest depth)
    // 2. Alternate allele fraction
    // 3. Strand bias
    // Total depth comes first: the allele fraction needs a depth above 0.
    let filtering: &DNAVariantFilteringOptions = &options.filtering;
    let filters: Vec<Box<dyn VariantFilter<Input = VariantCall> + '_>> = vec![
        Box::new(DNAVariantReadSupportFilter::new(
            filtering.min_reads as ReadSupport,
            filtering.min_total_depth as ReadSupport,
            filtering.min_homopolymer_len,
            filtering.min_dinucleotide_context_len,
            filtering.max_slippage_repeat_len,
            Some(max_depth),
            min_read_support_index,
            &read_depths,
            chromosome_names_map,
            fasta_map
        )),
        Box::new(VariantAlleleFractionFilter::new(filtering.min_alt_allele_fraction)),
        Box::new(StrandBiasFilter::new(&read_depths, chromosome_names_map, 0.05))
    ];
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    thread_pool.install(|| {
        variant_calls
            .into_par_iter()
            .filter(|vc| filters.iter().all(|filter| filter.passes(vc)))
            .collect()
    })
}


fn retype_control_insertions(
    control_variant_records: HashSet<VariantRecord>,
    case_variant_records: &Vec<Arc<VariantRecord>>,
    fasta_map: &FastaMap,
    chromosome_names_map: &BiMap<ReferenceChromosomeName, ReferenceChromosomeID>,
    options: &DNAVariantCallingOptions
) -> Vec<Arc<VariantRecord>> {
    let control_variant_records: Vec<Arc<VariantRecord>> = control_variant_records
        .into_iter()
        .map(Arc::new)
        .collect();
    if options.bkpt_rescue == false {
        return control_variant_records;
    }

    // Create an IntervalTree of the BND and TRA positions of the control and of the case.
    let breakpoint_variant_records: Vec<Arc<VariantRecord>> = control_variant_records
        .iter()
        .chain(case_variant_records.iter())
        .cloned()
        .collect();
    // HashMap<chromosome ID, IntervalTree<(chromosome 1 ID, position 1, chromosome 2 ID, position 2, record index)>>
    let interval_tree_map: HashMap<ReferenceChromosomeID, IntervalTree<(ReferenceChromosomeID, ReferencePosition, ReferenceChromosomeID, ReferencePosition, usize)>> = build_dna_breakpoint_interval_trees(
        &breakpoint_variant_records,
        options.bkpt_rescue_search_distance
    );

    control_variant_records
        .into_iter()
        .flat_map(|vr| -> Vec<Arc<VariantRecord>> {
            if vr.get_variant_type() != &VariantType::Insertion {
                return vec![vr];
            }

            // Get mate 2 positions of BND or TRA whose mate 1 position
            // is near the insertion position.
            let breakpoint_positions: Vec<&(ReferenceChromosomeID, ReferencePosition, ReferenceChromosomeID, ReferencePosition, usize)> = match interval_tree_map.get(&vr.get_chromosome_1()) {
                Some(interval_tree) => interval_tree.overlaps(vr.get_position_1() as isize, vr.get_position_2() as isize),
                None => Vec::new()
            };
            if breakpoint_positions.is_empty() {
                // No BND or TRA exists that is near the insertion.
                return vec![vr];
            }

            // Align INS sequence.
            let mut query_positions: HashMap<ReferenceChromosomeID, Vec<ReferencePosition>> = HashMap::new();
            for (chromosome_1, position_1, chromosome_2, position_2, _) in breakpoint_positions {
                query_positions.entry(*chromosome_1).or_default().push(*position_1);
                query_positions.entry(*chromosome_2).or_default().push(*position_2);
            }
            match retype_insertion_as_breakpoints(
                &vr,
                &query_positions,
                fasta_map,
                chromosome_names_map,
                options.bkpt_rescue_min_ins_len,
                options.bkpt_rescue_max_ins_len,
                options.bkpt_rescue_search_distance,
                options.bkpt_rescue_realignment_gap_open_score,
                options.bkpt_rescue_realignment_gap_extend_score,
                options.bkpt_rescue_realignment_k,
                options.bkpt_rescue_realignment_band_width,
                options.bkpt_rescue_realignment_min_score_fraction,
                options.bkpt_rescue_realignment_min_query_coverage,
                options.bkpt_rescue_realignment_min_span_proportion
            ) {
                Some(retyped_records) => retyped_records
                    .into_iter()
                    .map(Arc::new)
                    .chain(std::iter::once(vr))
                    .collect(),
                None => vec![vr]
            }
        })
        .collect()
}


#[cfg(test)]
#[path = "../tests/pipeline/variant_calling_dna_somatic.rs"]
mod tests;