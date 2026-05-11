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
use noodles_fasta::io::indexed_reader::Builder;
use noodles_sam::Header;
use rayon::prelude::*;
use rayon::ThreadPool;
use std::cmp::max;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::mem;
use std::sync::Arc;
use noodles_bgzf::io::{BufRead, Seek};
use noodles_sam::alignment::Record;
use sysinfo::System;
use tempfile::{NamedTempFile, TempPath};

use crate::prelude::*;


/// Identify germline DNA variants.
///
/// # Parameters:
///
/// * `min_reads` is the minimum number of supporting reads for a variant call.
/// * `min_size_proportion` is the minimum proportion that two variant sizes must share in
/// order to be considered clusterable. The value should be between 0.0 and 1.0.
/// * `max_ins_norm_edit_distance` is the maximum normalized edit distance for two insertions to be
/// considered clusterable. The value should be between 0.0 and 1.0.
/// * `max_intrachromosomal_distance_coefficient` is the `a` term of the formula for determining the `max_breakpoint_distance`:
/// `(mx)/(10m+x)`. The term `x` is the variant size. The bigger of the two variant sizes
/// is used for `x`. Both breakpoint pairs must be within max_breakpoint_distance to be considered
/// clusterable.
/// * `max_interchromosomal_distance` is the maximum distance between two breakpoints for
/// interchromosomal translocations.
/// * `chromosomes` is a vector of chromosomes in which variants should be called.
pub fn identify_germline_dna_variants(
    bam_file: &str,
    bam_bai_file: &str,
    fasta_file: &str,
    regions: &Vec<(&str, u32, u32)>,
    min_reads: usize,
    min_mapping_quality: u16,
    min_base_quality: u8,
    min_total_depth: u32,
    min_alt_allele_fraction: f32,
    min_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    max_intrachromosomal_distance_tau: u32,
    max_intrachromosomal_distance: u32,
    max_interchromosomal_distance: u32,
    max_slippage_repeat_length: u32,
    num_threads: usize,
    chunk_size: u32,
    max_records: usize,
    expected_variant_allele_fraction: f64,
    expected_mutation_rate: f64,
    expected_sequencing_error: f64,
    expected_slippage_probability: f64,
    max_f1_fraction: f64,
    max_fpr: f64,
    temp_dir: &str
) -> DNAVariantCallSet {
    // Step 1. Get all chromosome IDs and names
    log_info!("Getting chromosomes.");
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosomes: Vec<&str> = chromosome_lengths
        .keys()
        .map(|k| k.as_ref())
        .collect();

    // Step 2. Get sequencing depths
    log_info!("Getting sequencing depths and strands.");
    let depths_map: Arc<HashMap<Box<str>, Vec<u32>>> = Arc::new(get_bam_depths_map(bam_file, num_threads));
    let strands_map: Arc<HashMap<Box<str>, Vec<(u32, u32)>>> = Arc::new(get_bam_strands_map(bam_file, num_threads));
    let max_depth: u32 = depths_map
        .values()
        .flat_map(|v| v.iter())
        .copied()
        .max()
        .unwrap();
    log_info!("\tMax sequencing depth: {max_depth}.");

    // Step 3. Compute minimum read support index
    log_info!("Computing minimum read support index.");
    let min_read_support_index: Vec<Vec<u32>> = compute_min_read_support_index(
        max_depth as u64,
        max_slippage_repeat_length,
        expected_variant_allele_fraction,
        expected_mutation_rate,
        expected_sequencing_error,
        expected_slippage_probability,
        max_f1_fraction,
        max_fpr
    );

    // Step 4. Index BAM records by read names
    log_info!("Indexing BAM records by read names.");
    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        num_threads
    );

    // Step 5. Generate regions
    log_info!("Generating genomic regions.");
    let ordered_regions: BTreeMap<Box<str>, Vec<(u32, u32)>> = if regions.is_empty() {
        let regions: HashMap<Box<str>, Vec<(u32, u32)>> = generate_regions(
            &chromosomes,
            &chromosome_lengths,
            chunk_size
        );
        let mut ordered_regions: BTreeMap<Box<str>, Vec<(u32, u32)>> = regions.into_iter().collect();
        for (_, vec) in ordered_regions.iter_mut() {
            vec.sort_by(|a, b| a.0.cmp(&b.0));
        }
        ordered_regions
    } else {
        let regions_split: Vec<(Box<str>, u32, u32)> = split_regions(regions, chunk_size);
        let mut ordered_regions: BTreeMap<Box<str>, Vec<(u32, u32)>> = regions_split
            .iter()
            .map(|(contig, start, end)| {
                (
                    contig.to_string().into_boxed_str(),
                    vec![(*start, *end)]
                )
            })
            .collect();
        for vec in ordered_regions.values_mut() {
            vec.sort_by_key(|(start, _)| *start);
        }
        ordered_regions
    };

    // Step 6. Load the FASTA file
    log_info!("Loading FASTA file.");
    let fasta_map: FastaMap = FastaMap::new(fasta_file);

    // Step 7. Fetch the temp directory
    let dir: String = if temp_dir.is_empty() {
        // Use TMPDIR environment variable or fallback to system temp directory
        env::var("TMPDIR").unwrap_or_else(|_| env::temp_dir().to_string_lossy().to_string())
    } else {
        temp_dir.to_string()
    };
    let dir_path = Path::new(&dir);
    if !dir_path.exists() {
        panic!("Directory does not exist: {}", dir);
    }

    // Step 8. Identify variant calls
    log_info!("Identifying DNA variants.");
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let mut variant_call_idx: usize = 1;
    let padding: u32 = 2 * max(max_interchromosomal_distance, max(max_intrachromosomal_distance, max_intrachromosomal_distance_tau));
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();
    let mut temp_files: Vec<TempPath> = Vec::new();
    for (chromosome, curr_regions) in &ordered_regions {
        let chromosome_id: u16 = chromosome_names_map.get_by_left(chromosome).unwrap().clone();
        let mut chromosome_variant_records: Vec<VariantRecord> = Vec::new();
        let mut prev_end: u32 = 0;
        for (start,end) in curr_regions {
            log_info!("\tIdentifying variant calls in {}:{}-{}.", chromosome, start, end);

            // Fetch BAM records
            let mut records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
                &mut reader,
                &header,
                &index,
                chromosome,
                *start,
                *end,
                &record_positions_map,
                &read_names_map,
                max_records,
                num_threads
            );
            log_info!("\t\t{} BAM record(s) fetched.", records_map.len());
            
            // Identify variant records
            let mut variant_records: HashSet<VariantRecord> = thread_pool.install(|| {
                records_map
                    .into_par_iter()
                    .flat_map(|(read_id, records)| {
                        let read_sequence: Box<str> = get_fastx_read_sequence(&records);
                        let quality_scores: Vec<u8> = get_fastx_base_quality_scores(&records);
                        let arc_records: Vec<Arc<bam::Record>> = records
                            .into_iter()
                            .map(|r| Arc::new(r))
                            .collect();
                        let alignment: Alignment = Alignment::new(
                            read_id,
                            &*read_sequence,
                            &quality_scores,
                            &arc_records
                        );
                        alignment
                            .get_alignment_structure()
                            .identify_variant_records(
                                min_mapping_quality,
                                min_base_quality,
                                AnalyteType::DNA
                            )
                    })
                    .collect()
            });

            // Filter variant records by chromosome (allow inter-chromosomal translocations)
            variant_records = thread_pool.install(|| {
                variant_records
                    .into_par_iter()
                    .filter(|vr| match vr.get_variant_type() {
                        VariantType::Translocation => {
                            vr.get_chromosome_1() == chromosome_id || vr.get_chromosome_2() == chromosome_id
                        },
                        _ => vr.get_chromosome_1() == chromosome_id
                    })
                    .collect()
            });
            log_info!("\t\t{} variant record(s) identified.", variant_records.len());

            // Move into chromosome variant records
            chromosome_variant_records.extend(variant_records.into_iter());

            // Keep sorted by position_1
            chromosome_variant_records.sort_unstable_by_key(|vr| vr.get_position_1());

            // Compute window indices [keep_from, end]
            let keep_from: u32 = prev_end.saturating_sub(padding);
            let keep_to: u32 = *end;

            // Slice variant records for clustering
            let start_idx: usize = chromosome_variant_records
                .binary_search_by_key(&keep_from, |vr| vr.get_position_1())
                .unwrap_or_else(|idx| idx);
            let end_idx: usize = chromosome_variant_records
                .binary_search_by_key(&keep_to, |vr| vr.get_position_1())
                .map(|idx| idx + 1)
                .unwrap_or_else(|idx| idx);
            let slice_for_clustering = &chromosome_variant_records[start_idx..end_idx];

            // Cluster variant records into variant calls
            let mut variant_calls: Vec<VariantCall> = cluster_variant_records(
                slice_for_clustering
                    .iter()
                    .cloned()
                    .map(Arc::new)
                    .collect(),
                &depths_map,
                &chromosome_names_map,
                num_threads,
                min_size_proportion,
                max_ins_norm_edit_distance,
                max_intrachromosomal_distance_tau,
                max_intrachromosomal_distance,
                max_interchromosomal_distance,
                false
            );

            // Filter variant calls for the following:
            // 1. Total depth
            // 2. Alternate allele fraction
            // 3. Variant read count
            // 4. Strand bias
            variant_calls = thread_pool.install(|| {
                variant_calls
                    .into_par_iter()
                    .filter(|vc| {
                        // Total depth
                        if vc.get_total_depth() < min_total_depth as i32 {
                            return false;
                        }

                        // Alternate allele fraction
                        if vc.get_alternate_allele_fraction() < min_alt_allele_fraction {
                            return false;
                        }

                        // Read count
                        let curr_read_support: usize = vc.get_read_ids().len();
                        if curr_read_support < min_reads {
                            return false;
                        }
                        let (vr, consensus_read_ids) = vc.get_consensus_record();
                        let chromosome_1: &Box<str> = chromosome_names_map.get_by_right(&vr.get_chromosome_1()).unwrap();
                        let chromosome_2: &Box<str> = chromosome_names_map.get_by_right(&vr.get_chromosome_2()).unwrap();
                        let position_1_depth: u32 = depths_map[chromosome_1][(vr.get_position_1() - 1) as usize];
                        let position_2_depth: u32 = depths_map[chromosome_2][(vr.get_position_2() - 1) as usize];
                        let (is_repeat, repeat_length_) = is_repeat_indel(vr, &chromosome_names_map, &fasta_map);
                        let repeat_length: usize = repeat_length_.min(max_slippage_repeat_length) as usize;

                        let min_read_support_1: usize = if position_1_depth == 0 {
                            usize::MAX
                        } else {
                            let repeat_idx: usize = if is_repeat { repeat_length } else { 0 };
                            min_read_support_index[repeat_idx][position_1_depth as usize - 1] as usize
                        };
                        let min_read_support_2: usize = if position_2_depth == 0 {
                            usize::MAX
                        } else {
                            let repeat_idx: usize = if is_repeat { repeat_length } else { 0 };
                            min_read_support_index[repeat_idx][position_2_depth as usize - 1] as usize
                        };
                        if curr_read_support < min_read_support_1 || curr_read_support < min_read_support_2 {
                            return false;
                        }

                        // Strand bias — skip for variants whose alt support is intrinsically
                        // strand-asymmetric (breakpoints / translocations). For these, all reads
                        // spanning the breakend share the same strand pair *by construction*, so
                        // a one-sided strand distribution is the variant signature, not an
                        // artifact. The Fisher's-exact check is only meaningful for SNVs/indels.
                        match vr.get_variant_type() {
                            VariantType::Breakpoint | VariantType::Translocation => {}
                            _ => {
                                let mut alt_fwd_1: u64 = 0u64;
                                let mut alt_rev_1: u64 = 0u64;
                                let mut alt_fwd_2: u64 = 0u64;
                                let mut alt_rev_2: u64 = 0u64;
                                for vr in &vc.variant_records {
                                    if *vr.get_strand_1() == Strand::Forward {
                                        alt_fwd_1 += 1;
                                    } else {
                                        alt_rev_1 += 1;
                                    }
                                    if *vr.get_strand_2() == Strand::Forward {
                                        alt_fwd_2 += 1;
                                    } else {
                                        alt_rev_2 += 1;
                                    }
                                }
                                let (ref_fwd_1, ref_rev_1) = strands_map.get(chromosome_1).unwrap().get(vr.get_position_1() as usize - 1).unwrap();
                                let (ref_fwd_2, ref_rev_2) = strands_map.get(chromosome_2).unwrap().get(vr.get_position_2() as usize - 1).unwrap();
                                if has_strand_bias(alt_fwd_1, alt_rev_1, *ref_fwd_1 as u64, *ref_rev_1 as u64, 0.05f64) ||
                                    has_strand_bias(alt_fwd_2, alt_rev_2, *ref_fwd_2 as u64, *ref_rev_2 as u64, 0.05f64) {
                                    return false;
                                }
                            }
                        }

                        true
                    })
                    .collect()
            });
            log_info!("\t\t{} variant call(s) identified.", variant_calls.len());

            // Aggregate variant calls into a variant callset
            let mut variant_call_set = DNAVariantCallSet::new();
            for mut variant_call in variant_calls {
                // Rename the variant call ID
                variant_call.id = variant_call_idx;
                variant_call_set.add_variant_call(variant_call);
                variant_call_idx += 1;
            }

            // Store variant_call_set in a temp file
            let temp_path: TempPath = {
                let temp_file = NamedTempFile::new_in(dir_path).unwrap();
                let mut writer = BufWriter::new(temp_file);
                bincode::serialize_into(&mut writer, &variant_call_set)
                    .expect("Failed to serialize variant_call_set");
                writer.flush().unwrap();
                let temp_file: NamedTempFile = writer.into_inner().unwrap();
                temp_file.into_temp_path()
            };
            temp_files.push(temp_path);

            prev_end = keep_to;
        }
        chromosome_variant_records.clear();
    }

    // Step 8. Load all VariantCallSet objects and merge them
    log_info!("Loading all temp files and merging them into a variant call set.");
    let variant_call_set: DNAVariantCallSet = merge_temp_variant_call_sets(
        &temp_files,
        &read_names_map,
        &chromosome_names_map,
        num_threads
    );

    variant_call_set
}

/// Identify somatic DNA variants.
///
/// # Parameters:
///
/// * `min_reads` is the minimum number of supporting reads for a variant call.
/// * `min_size_proportion` is the minimum proportion that two variant sizes must share in
/// order to be considered clusterable. The value should be between 0.0 and 1.0.
/// * `max_ins_norm_edit_distance` is the maximum normalized edit distance for two insertions to be
/// considered clusterable. The value should be between 0.0 and 1.0.
/// * `max_intrachromosomal_distance_coefficient` is the `a` term of the formula for determining the `max_breakpoint_distance`:
/// `(mx)/(10m+x)`. The term `x` is the variant size. The bigger of the two variant sizes
/// is used for `x`. Both breakpoint pairs must be within max_breakpoint_distance to be considered
/// clusterable.
/// * `max_interchromosomal_distance` is the maximum distance between two breakpoints for
/// interchromosomal translocations.
/// * `chunk_size` is the window of a chromosome to process each time (increasing this value
/// increases the maximum peak memory used).
/// * `chromosomes` is a vector of chromosomes in which variants should be called.
/// * If `apply_infinite_sites_assumption` is true, any `a` variant record that shares a breakpoint
/// (either position_1 or position_2) with any of the `b` variant record will be filtered out.
pub fn identify_somatic_dna_variants(
    case_bam_file: &str,
    case_bam_bai_file: &str,
    control_bam_files: Vec<&str>,
    control_bam_bai_files: Vec<&str>,
    fasta_file: &str,
    regions: &Vec<(&str, u32, u32)>,
    min_reads: usize,
    min_mapping_quality: u16,
    min_base_quality: u8,
    min_total_depth: u32,
    min_alt_allele_fraction: f32,
    min_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    max_intrachromosomal_distance_tau: u32,
    max_intrachromosomal_distance: u32,
    max_interchromosomal_distance: u32,
    max_slippage_repeat_length: u32,
    num_threads: usize,
    chunk_size: u32,
    max_records: usize,
    expected_variant_allele_fraction: f64,
    expected_mutation_rate: f64,
    expected_sequencing_error: f64,
    expected_slippage_probability: f64,
    max_f1_fraction: f64,
    max_fpr: f64,
    apply_infinite_sites_assumption: bool,
    temp_dir: &str
) -> DNAVariantCallSet {
    assert!(control_bam_files.len() == control_bam_bai_files.len());

    // Step 1. Get all chromosome IDs and names
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(case_bam_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(case_bam_file);
    let chromosomes: Vec<&str> = chromosome_lengths
        .keys()
        .map(|k| k.as_ref())
        .collect();

    // Step 2. Make sure the chromosome IDs and names are the same in the control BAM files
    for control_bam_file in control_bam_files.iter() {
        let chromosome_names_map_: BiMap<Box<str>, u16> = create_chromosome_names_map(control_bam_file);
        let chromosome_lengths_: HashMap<Box<str>, u32> = get_chromosome_lengths(control_bam_file);
        for (chromosome_name,chromosome_id) in chromosome_names_map.iter() {
            if chromosome_names_map_.contains_left(chromosome_name) {
                let chromosome_id_: u16 = *chromosome_names_map_.get_by_left(chromosome_name).unwrap();
                if chromosome_id_ != *chromosome_id {
                    panic!("Mismatch in chromosome names map for BAM file: {:?}", control_bam_file);
                }
            } else {
                panic!("Mismatch in chromosome names map for BAM file: {:?}", control_bam_file);
            }
        }
        if chromosome_lengths != chromosome_lengths_ {
            panic!("Mismatch in chromosome lengths for BAM file: {:?}", control_bam_file);
        }
    }

    // Step 3. Get sequencing depths
    log_info!("Getting sequencing depths and strands.");
    let depths_map: Arc<HashMap<Box<str>, Vec<u32>>> = Arc::new(get_bam_depths_map(case_bam_file, num_threads));
    let strands_map: Arc<HashMap<Box<str>, Vec<(u32, u32)>>> = Arc::new(get_bam_strands_map(case_bam_file, num_threads));
    let max_depth: u32 = depths_map
        .values()
        .flat_map(|v| v.iter())
        .copied()
        .max()
        .unwrap();
    log_info!("\tMax sequencing depth: {max_depth}.");

    // Step 4. Compute minimum read support index
    log_info!("Computing minimum read support index.");
    let min_read_support_index: Vec<Vec<u32>> = compute_min_read_support_index(
        max_depth as u64,
        max_slippage_repeat_length,
        expected_variant_allele_fraction,
        expected_mutation_rate,
        expected_sequencing_error,
        expected_slippage_probability,
        max_f1_fraction,
        max_fpr
    );

    // Step 5. Index BAM records by read names
    log_info!("Indexing BAM records by read names.");
    let (case_record_positions_map, case_read_names_map) = index_bam_records(
        case_bam_file,
        num_threads
    );
    log_info!("{} read names/IDs in case BAM file.", case_read_names_map.len());
    let mut control_record_positions_map: HashMap<Box<str>, HashMap<usize, Vec<VirtualPosition>>> = HashMap::new();
    let mut control_read_names_map: HashMap<Box<str>, BiMap<Box<str>, usize>> = HashMap::new();
    for (index,control_bam_file) in control_bam_files.iter().enumerate() {
        let (record_positions_map, read_names_map) = index_bam_records(
            control_bam_file,
            num_threads
        );
        log_info!("{} read names/IDs in control BAM file {}/{}.", read_names_map.len(), index + 1, control_bam_files.len());
        control_record_positions_map.insert(control_bam_file.to_string().into_boxed_str(), record_positions_map);
        control_read_names_map.insert(control_bam_file.to_string().into_boxed_str(), read_names_map);
    }

    // Step 6. Generate regions
    log_info!("Generating genomic regions.");
    let ordered_regions: BTreeMap<Box<str>, Vec<(u32, u32)>> = if regions.is_empty() {
        let regions: HashMap<Box<str>, Vec<(u32, u32)>> = generate_regions(
            &chromosomes,
            &chromosome_lengths,
            chunk_size
        );
        let mut ordered_regions: BTreeMap<Box<str>, Vec<(u32, u32)>> = regions.into_iter().collect();
        for (_, vec) in ordered_regions.iter_mut() {
            vec.sort_by(|a, b| a.0.cmp(&b.0));
        }
        ordered_regions
    } else {
        let regions_split: Vec<(Box<str>, u32, u32)> = split_regions(regions, chunk_size);
        let mut ordered_regions: BTreeMap<Box<str>, Vec<(u32, u32)>> = regions_split
            .iter()
            .map(|(contig, start, end)| {
                (
                    contig.to_string().into_boxed_str(),
                    vec![(*start, *end)]
                )
            })
            .collect();
        for vec in ordered_regions.values_mut() {
            vec.sort_by_key(|(start, _)| *start);
        }
        ordered_regions
    };

    // Step 7. Load the FASTA file
    log_info!("Loading FASTA file.");
    let fasta_map: FastaMap = FastaMap::new(fasta_file);

    // Step 8. Fetch the temp directory
    let dir: String = if temp_dir.is_empty() {
        // Use TMPDIR environment variable or fallback to system temp directory
        env::var("TMPDIR").unwrap_or_else(|_| env::temp_dir().to_string_lossy().to_string())
    } else {
        temp_dir.to_string()
    };
    let dir_path: &Path = Path::new(&dir);
    if !dir_path.exists() {
        panic!("Directory does not exist: {}", dir);
    }

    // Step 9. Identify case-specific variant calls
    log_info!("Identifying case-specific DNA variants.");
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let mut variant_call_idx: usize = 1;
    let padding: u32 = 2 * max(max_interchromosomal_distance, max(max_intrachromosomal_distance, max_intrachromosomal_distance_tau));
    let mut case_reader = bam::io::reader::Builder::default()
        .build_from_path(case_bam_file)
        .unwrap();
    let mut max_distance: u32 = max(max_intrachromosomal_distance_tau, max_intrachromosomal_distance);
    max_distance = max(max_distance, max_interchromosomal_distance);
    let bin_size: u32 = 10_u32.pow((max_distance as f32).log10().floor() as u32 + 1);
    let case_header: Header = case_reader.read_header().unwrap();
    let case_index: Index = bai::fs::read(case_bam_bai_file).unwrap();
    let mut control_readers = HashMap::new();
    let mut control_headers = HashMap::new();
    let mut control_indices = HashMap::new();
    for (i, control_bam_file) in control_bam_files.iter().enumerate() {
        let control_bam_bai_file: &str = control_bam_bai_files[i];
        let mut reader = bam::io::reader::Builder::default()
            .build_from_path(control_bam_file)
            .unwrap();
        let header: Header = reader.read_header().unwrap();
        let index = bai::fs::read(control_bam_bai_file).unwrap();
        control_readers.insert(control_bam_file, reader);
        control_headers.insert(control_bam_file, header);
        control_indices.insert(control_bam_file, index);
    }
    let mut temp_files: Vec<TempPath> = Vec::new();
    for (chromosome, curr_regions) in &ordered_regions {
        let chromosome_id: u16 = chromosome_names_map.get_by_left(chromosome).unwrap().clone();
        let mut chromosome_case_variant_records: Vec<VariantRecord> = Vec::new();
        let mut prev_end: u32 = 0;
        for (start,end) in curr_regions {
            log_info!("\tProcessing {}:{}-{}.", chromosome, start, end);

            // Fetch case BAM records
            let mut case_records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
                &mut case_reader,
                &case_header,
                &case_index,
                chromosome,
                *start,
                *end,
                &case_record_positions_map,
                &case_read_names_map,
                max_records,
                num_threads
            );
            log_info!("\t\t{} BAM record(s) fetched.", case_records_map.len());

            // Identify case variant records
            let mut case_variant_records: HashSet<VariantRecord> = thread_pool.install(|| {
                case_records_map
                    .into_par_iter()
                    .flat_map(|(read_id, records)| {
                        let read_sequence: Box<str> =  get_fastx_read_sequence(&records);
                        let quality_scores: Vec<u8> = get_fastx_base_quality_scores(&records);
                        let arc_records: Vec<Arc<bam::Record>> = records
                            .into_iter()
                            .map(|r| Arc::new(r)) // move, not clone
                            .collect();
                        let alignment: Alignment = Alignment::new(
                            read_id,
                            &*read_sequence,
                            &quality_scores,
                            &arc_records
                        );
                        alignment
                            .get_alignment_structure()
                            .identify_variant_records(
                                min_mapping_quality,
                                min_base_quality,
                                AnalyteType::DNA
                            )
                    })
                    .collect()
            });

            // Filter case variant records by chromosome (allow inter-chromosomal translocations)
            case_variant_records = thread_pool.install(|| {
                case_variant_records
                    .into_par_iter()
                    .filter(|vr| match vr.get_variant_type() {
                        VariantType::Translocation => {
                            vr.get_chromosome_1() == chromosome_id || vr.get_chromosome_2() == chromosome_id
                        },
                        _ => vr.get_chromosome_1() == chromosome_id
                    })
                    .collect()
            });
            log_info!("\t\t{} case variant record(s) identified.", case_variant_records.len());

            // Move into chromosome variant records
            chromosome_case_variant_records.extend(case_variant_records.into_iter());

            // Keep sorted by position_1
            chromosome_case_variant_records.sort_unstable_by_key(|vr| vr.get_position_1());

            // Compute window indices [keep_from, end]
            let keep_from: u32 = prev_end.saturating_sub(padding);
            let keep_to: u32 = *end;

            // Slice case variant records for clustering
            let start_idx: usize = chromosome_case_variant_records
                .binary_search_by_key(&keep_from, |vr| vr.get_position_1())
                .unwrap_or_else(|idx| idx);
            let end_idx: usize = chromosome_case_variant_records
                .binary_search_by_key(&keep_to, |vr| vr.get_position_1())
                .map(|idx| idx + 1)
                .unwrap_or_else(|idx| idx);
            let slice_for_clustering: &[VariantRecord] = &chromosome_case_variant_records[start_idx..end_idx];

            // Cluster case variant records into variant calls
            let mut variant_calls: Vec<VariantCall> = cluster_variant_records(
                slice_for_clustering
                    .iter()
                    .cloned()
                    .map(Arc::new)
                    .collect(),
                &depths_map,
                &chromosome_names_map,
                num_threads,
                min_size_proportion,
                max_ins_norm_edit_distance,
                max_intrachromosomal_distance_tau,
                max_intrachromosomal_distance,
                max_interchromosomal_distance,
                false
            );

            // Filter variant calls for the following:
            // 1. Total depth
            // 2. Alternate allele fraction
            // 3. Variant read count
            // 4. Strand bias
            variant_calls = thread_pool.install(|| {
                variant_calls
                    .into_par_iter()
                    .filter(|vc| {
                        // Total depth
                        if vc.get_total_depth() < min_total_depth as i32 {
                            return false;
                        }

                        // Alternate allele fraction
                        if vc.get_alternate_allele_fraction() < min_alt_allele_fraction {
                            return false;
                        }

                        // Read count
                        let curr_read_support: usize = vc.get_read_ids().len();
                        if curr_read_support < min_reads {
                            return false;
                        }
                        let (vr, consensus_read_ids) = vc.get_consensus_record();
                        let chromosome_1: &Box<str> = chromosome_names_map.get_by_right(&vr.get_chromosome_1()).unwrap();
                        let chromosome_2: &Box<str> = chromosome_names_map.get_by_right(&vr.get_chromosome_2()).unwrap();
                        let position_1_depth: u32 = depths_map[chromosome_1][(vr.get_position_1() - 1) as usize];
                        let position_2_depth: u32 = depths_map[chromosome_2][(vr.get_position_2() - 1) as usize];
                        let (is_repeat, repeat_length_) = is_repeat_indel(vr, &chromosome_names_map, &fasta_map);
                        let repeat_length: usize = repeat_length_.min(max_slippage_repeat_length) as usize;
                        let min_read_support_1: usize = if position_1_depth == 0 {
                            let repeat_idx: usize = if is_repeat { repeat_length } else { 0 };
                            min_read_support_index[repeat_idx][max_depth as usize - 1] as usize
                        } else {
                            let repeat_idx: usize = if is_repeat { repeat_length } else { 0 };
                            min_read_support_index[repeat_idx][position_1_depth as usize - 1] as usize
                        };
                        let min_read_support_2: usize = if position_2_depth == 0 {
                            let repeat_idx: usize = if is_repeat { repeat_length } else { 0 };
                            min_read_support_index[repeat_idx][max_depth as usize - 1] as usize
                        } else {
                            let repeat_idx: usize = if is_repeat { repeat_length } else { 0 };
                            min_read_support_index[repeat_idx][position_2_depth as usize - 1] as usize
                        };
                        if curr_read_support < min_read_support_1 || curr_read_support < min_read_support_2 {
                            return false;
                        }

                        // Strand bias — skip for variants whose alt support is intrinsically
                        // strand-asymmetric (breakpoints / translocations). For these, all reads
                        // spanning the breakend share the same strand pair *by construction*, so
                        // a one-sided strand distribution is the variant signature, not an
                        // artifact. The Fisher's-exact check is only meaningful for SNVs/indels.
                        match vr.get_variant_type() {
                            VariantType::Breakpoint | VariantType::Translocation => {}
                            _ => {
                                let mut alt_fwd_1: u64 = 0u64;
                                let mut alt_rev_1: u64 = 0u64;
                                let mut alt_fwd_2: u64 = 0u64;
                                let mut alt_rev_2: u64 = 0u64;
                                for vr in &vc.variant_records {
                                    if *vr.get_strand_1() == Strand::Forward {
                                        alt_fwd_1 += 1;
                                    } else {
                                        alt_rev_1 += 1;
                                    }
                                    if *vr.get_strand_2() == Strand::Forward {
                                        alt_fwd_2 += 1;
                                    } else {
                                        alt_rev_2 += 1;
                                    }
                                }
                                let (ref_fwd_1, ref_rev_1) = strands_map.get(chromosome_1).unwrap().get(vr.get_position_1() as usize - 1).unwrap();
                                let (ref_fwd_2, ref_rev_2) = strands_map.get(chromosome_2).unwrap().get(vr.get_position_2() as usize - 1).unwrap();
                                if has_strand_bias(alt_fwd_1, alt_rev_1, *ref_fwd_1 as u64, *ref_rev_1 as u64, 0.05f64) ||
                                    has_strand_bias(alt_fwd_2, alt_rev_2, *ref_fwd_2 as u64, *ref_rev_2 as u64, 0.05f64) {
                                    return false;
                                }
                            }
                        }

                        true
                    })
                    .collect()
            });
            log_info!("\t\t{} case variant call(s) identified.", variant_calls.len());

            let mut case_variant_records_filtered: Vec<Arc<VariantRecord>> = Vec::new();
            for variant_call in variant_calls {
                for variant_record in variant_call.variant_records {
                    case_variant_records_filtered.push(Arc::new(variant_record));
                }
            }
            log_info!("\t\t{} case variant record(s) retained after filtering.", case_variant_records_filtered.len());

            // Filter out variant records near control variant records
            for (i, (control_bam_file, mut control_reader)) in control_readers.iter_mut().enumerate() {
                if case_variant_records_filtered.is_empty() == false {
                    log_info!("\t\tProcessing control BAM file {}/{}.", i+1, control_bam_files.len());

                    // Fetch case BAM records
                    let mut control_records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
                        &mut control_reader,
                        control_headers.get(control_bam_file).unwrap(),
                        control_indices.get(control_bam_file).unwrap(),
                        chromosome,
                        max(1u32, (*start).saturating_sub(padding)),
                        (*end).saturating_add(padding),
                        &control_record_positions_map.get(&control_bam_file.to_string().into_boxed_str()).unwrap(),
                        &control_read_names_map.get(&control_bam_file.to_string().into_boxed_str()).unwrap(),
                        max_records,
                        num_threads
                    );
                    log_info!("\t\t\t{} BAM record(s) fetched.", control_records_map.len());

                    // Identify control variant records
                    let mut control_variant_records: HashSet<VariantRecord> = thread_pool.install(|| {
                        control_records_map
                            .into_par_iter()
                            .flat_map(|(read_id, records)| {
                                let read_sequence: Box<str> = get_fastx_read_sequence(&records);
                                let quality_scores: Vec<u8> = get_fastx_base_quality_scores(&records);
                                let arc_records: Vec<Arc<bam::Record>> = records
                                    .into_iter()
                                    .map(|r| Arc::new(r)) // move, not clone
                                    .collect();
                                let alignment: Alignment = Alignment::new(
                                    read_id,
                                    &*read_sequence,
                                    &quality_scores,
                                    &arc_records
                                );
                                alignment
                                    .get_alignment_structure()
                                    .identify_variant_records(
                                        min_mapping_quality,
                                        min_base_quality,
                                        AnalyteType::DNA
                                    )
                            })
                            .collect()
                    });
                    log_info!("\t\t\t{} control variant record(s) identified.", control_variant_records.len());

                    // Filter by chromosome (allow inter-chromosomal translocations)
                    control_variant_records.retain(|vr| match vr.get_variant_type() {
                        VariantType::Translocation => {
                            vr.get_chromosome_1() == chromosome_id || vr.get_chromosome_2() == chromosome_id
                        },
                        _ => vr.get_chromosome_1() == chromosome_id
                    });

                    // Filter out control variants
                    case_variant_records_filtered = diff_variant_records(
                        case_variant_records_filtered,
                        control_variant_records.into_iter().map(Arc::new).collect(),
                        bin_size,
                        num_threads,
                        min_size_proportion,
                        max_ins_norm_edit_distance,
                        max_intrachromosomal_distance_tau,
                        max_intrachromosomal_distance,
                        max_interchromosomal_distance,
                        apply_infinite_sites_assumption,
                        false
                    );
                    log_info!("\t\t\t{} case variant record(s) retained after diffing.", case_variant_records_filtered.len());
                }
            }

            let mut variant_calls: Vec<VariantCall> = cluster_variant_records(
                case_variant_records_filtered,
                &depths_map,
                &chromosome_names_map,
                num_threads,
                min_size_proportion,
                max_ins_norm_edit_distance,
                max_intrachromosomal_distance_tau,
                max_intrachromosomal_distance,
                max_interchromosomal_distance,
                false
            );
            log_info!("\t\t\t{} case-specific variant call(s) identified.", variant_calls.len());

            // Filter variant calls for the following:
            // 1. Total depth
            // 2. Alternate allele fraction
            // 3. Variant read count
            // 4. Strand bias
            variant_calls = thread_pool.install(|| {
                variant_calls
                    .into_par_iter()
                    .filter(|vc| {
                        // Total depth
                        if vc.get_total_depth() < min_total_depth as i32 {
                            return false;
                        }

                        // Alternate allele fraction
                        if vc.get_alternate_allele_fraction() < min_alt_allele_fraction {
                            return false;
                        }

                        // Read count
                        let curr_read_support: usize = vc.get_read_ids().len();
                        if curr_read_support < min_reads {
                            return false;
                        }
                        let (vr, consensus_read_ids) = vc.get_consensus_record();
                        let chromosome_1: &Box<str> = chromosome_names_map.get_by_right(&vr.get_chromosome_1()).unwrap();
                        let chromosome_2: &Box<str> = chromosome_names_map.get_by_right(&vr.get_chromosome_2()).unwrap();
                        let position_1_depth: u32 = depths_map[chromosome_1][(vr.get_position_1() - 1) as usize];
                        let position_2_depth: u32 = depths_map[chromosome_2][(vr.get_position_2() - 1) as usize];
                        let (is_repeat, repeat_length_) = is_repeat_indel(vr, &chromosome_names_map, &fasta_map);
                        let repeat_length: usize = repeat_length_.min(max_slippage_repeat_length) as usize;
                        let min_read_support_1: usize = if position_1_depth == 0 {
                            let repeat_idx: usize = if is_repeat { repeat_length } else { 0 };
                            min_read_support_index[repeat_idx][max_depth as usize - 1] as usize
                        } else {
                            let repeat_idx: usize = if is_repeat { repeat_length } else { 0 };
                            min_read_support_index[repeat_idx][position_1_depth as usize - 1] as usize
                        };
                        let min_read_support_2: usize = if position_2_depth == 0 {
                            let repeat_idx: usize = if is_repeat { repeat_length } else { 0 };
                            min_read_support_index[repeat_idx][max_depth as usize - 1] as usize
                        } else {
                            let repeat_idx: usize = if is_repeat { repeat_length } else { 0 };
                            min_read_support_index[repeat_idx][position_2_depth as usize - 1] as usize
                        };
                        if curr_read_support < min_read_support_1 || curr_read_support < min_read_support_2 {
                            return false;
                        }

                        // Strand bias — skip for variants whose alt support is intrinsically
                        // strand-asymmetric (breakpoints / translocations). For these, all reads
                        // spanning the breakend share the same strand pair *by construction*, so
                        // a one-sided strand distribution is the variant signature, not an
                        // artifact. The Fisher's-exact check is only meaningful for SNVs/indels.
                        match vr.get_variant_type() {
                            VariantType::Breakpoint | VariantType::Translocation => {}
                            _ => {
                                let mut alt_fwd_1: u64 = 0u64;
                                let mut alt_rev_1: u64 = 0u64;
                                let mut alt_fwd_2: u64 = 0u64;
                                let mut alt_rev_2: u64 = 0u64;
                                for vr in &vc.variant_records {
                                    if *vr.get_strand_1() == Strand::Forward {
                                        alt_fwd_1 += 1;
                                    } else {
                                        alt_rev_1 += 1;
                                    }
                                    if *vr.get_strand_2() == Strand::Forward {
                                        alt_fwd_2 += 1;
                                    } else {
                                        alt_rev_2 += 1;
                                    }
                                }
                                let (ref_fwd_1, ref_rev_1) = strands_map.get(chromosome_1).unwrap().get(vr.get_position_1() as usize - 1).unwrap();
                                let (ref_fwd_2, ref_rev_2) = strands_map.get(chromosome_2).unwrap().get(vr.get_position_2() as usize - 1).unwrap();
                                if has_strand_bias(alt_fwd_1, alt_rev_1, *ref_fwd_1 as u64, *ref_rev_1 as u64, 0.05f64) ||
                                    has_strand_bias(alt_fwd_2, alt_rev_2, *ref_fwd_2 as u64, *ref_rev_2 as u64, 0.05f64) {
                                    return false;
                                }
                            }
                        }

                        true
                    })
                    .collect()
            });

            let mut variant_call_set: DNAVariantCallSet = DNAVariantCallSet::new();
            for mut variant_call in variant_calls {
                // Rename the variant call ID
                variant_call.id = variant_call_idx;
                variant_call_set.add_variant_call(variant_call);
                variant_call_idx += 1;
            }

            // Store variant_call_set in a temp file
            let temp_path: TempPath = {
                let temp_file = NamedTempFile::new_in(dir_path).unwrap();
                let mut writer = BufWriter::new(temp_file);
                bincode::serialize_into(&mut writer, &variant_call_set)
                    .expect("Failed to serialize variant_call_set");
                writer.flush().unwrap();
                let temp_file: NamedTempFile = writer.into_inner().unwrap();
                temp_file.into_temp_path()
            };
            temp_files.push(temp_path);

            prev_end = keep_to;

        }
        chromosome_case_variant_records.clear();
    }

    // Step 10. Load all VariantCallSet objects and merge them
    log_info!("Loading all temp files and merging them into a variant call set.");
    let variant_call_set: DNAVariantCallSet = merge_temp_variant_call_sets(
        &temp_files,
        &case_read_names_map,
        &chromosome_names_map,
        num_threads
    );

    variant_call_set
}

fn merge_temp_variant_call_sets(
    temp_files: &[TempPath],
    read_names_map: &BiMap<Box<str>, usize>,
    chromosome_names_map: &BiMap<Box<str>, u16>,
    num_threads: usize
) -> DNAVariantCallSet {
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    let variant_calls: HashSet<VariantCall> = thread_pool.install(|| {
        temp_files
            .par_iter()
            .map(|temp_path| {
                let file = File::open(temp_path).unwrap();
                let mut reader = BufReader::new(file);
                let set: DNAVariantCallSet = bincode::deserialize_from(&mut reader).expect("Failed to deserialize data");
                let mut local_variant_calls: HashSet<VariantCall> = HashSet::with_capacity(set.variant_calls.len());
                for (_, vc) in set.variant_calls {
                    local_variant_calls.insert(vc);
                }
                local_variant_calls
            })
            .reduce(
                || HashSet::new(),
                |mut acc, local| {
                    // Reduce rehashing during merge
                    acc.reserve(local.len());
                    acc.extend(local);
                    acc
                }
            )
    });

    let mut variant_call_set: DNAVariantCallSet = DNAVariantCallSet::new();

    let mut variant_call_id: usize = 1;
    for mut variant_call in variant_calls {
        variant_call.id = variant_call_id;
        variant_call_set.add_variant_call(variant_call);
        variant_call_id += 1;
    }

    variant_call_set.load_read_names(read_names_map.clone());
    variant_call_set.load_chromosome_names(chromosome_names_map.clone());

    variant_call_set
}
