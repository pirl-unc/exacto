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
use indicatif::{ProgressBar, ProgressStyle};
use noodles_bam as bam;
use noodles_fasta::io::indexed_reader::Builder;
use rayon::prelude::*;
use std::cmp::max;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::mem;
use std::sync::Arc;
use sysinfo::System;
use tempfile::NamedTempFile;

use crate::prelude::*;
use crate::log_info;


/// Identify DNA variants.
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
pub fn identify_dna_variants(
    bam_file: &str,
    bam_bai_file: &str,
    min_reads: usize,
    min_mapping_quality: u32,
    min_base_quality: u8,
    min_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    max_intrachromosomal_distance_tau: u32,
    max_intrachromosomal_distance: u32,
    max_interchromosomal_distance: u32,
    num_threads: usize,
    chromosomes: Vec<&str>,
    temp_dir: &str
) -> DNAVariantCallSet {
    // Step 1. Get read IDs map
    log_info!("Converting read names to IDs.");
    let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        bam_file,
        bam_bai_file,
        num_threads
    );

    // Step 2. Get all chromosome IDs and names
    let chromosome_names_map: BiMap<Box<str>,u16> = create_chromosome_names_map(bam_file);
    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);

    // Step 3. Generate regions
    let regions: HashMap<Box<str>,Vec<(usize,usize)>> = generate_regions(
        bam_file,
        &chromosomes,
        *chromosome_lengths.values().max().unwrap()
    );
    let mut ordered_regions: BTreeMap<Box<str>,Vec<(usize,usize)>> = regions.into_iter().collect();
    for (_, vec) in ordered_regions.iter_mut() {
        vec.sort_by(|a, b| a.0.cmp(&b.0));
    }

    // Step 4. Fetch the temp directory
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

    // Step 5. Identify variant calls
    let mut temp_files: Vec<NamedTempFile> = Vec::new();
    let padding: isize = 2 * std::cmp::max(max_interchromosomal_distance, std::cmp::max(max_intrachromosomal_distance,max_intrachromosomal_distance_tau)) as isize;
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    log_info!("Identifying DNA variants.");
    let pb = Arc::new(ProgressBar::new(ordered_regions.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-")
    );
    let mut variant_call_idx: usize = 1;
    for (chromosome,curr_regions) in &ordered_regions {
        let mut variant_records_btree: BTreeMap<usize,HashSet<VariantRecord>> = BTreeMap::new();
        let mut prev_end: isize = 0;
        for (start,end) in curr_regions.iter() {
            log_info!("Identifying variant calls in {}:{}-{}.", chromosome, start, end);

            // Fetch BAM records
            let mut records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
                bam_file,
                bam_bai_file,
                chromosome,
                *start,
                *end,
                &read_names_map,
                num_threads
            );

            // Identify variant records
            log_info!("\tIdentifying variant records.");
            let mut variant_records: HashSet<VariantRecord> = thread_pool.install(|| {
                records_map
                    .par_iter()
                    .map(|(read_id, records)| {
                        let read_sequence: Box<str> =  get_fastx_read_sequence(records.iter().collect::<Vec<_>>().as_slice());
                        let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records.iter().collect::<Vec<_>>().as_slice());
                        let alignment: Alignment = Alignment::new(
                            *read_id,
                            &*read_sequence,
                            &quality_scores,
                            records
                        );
                        let variant_records: Vec<VariantRecord> = alignment
                            .get_alignment_structure()
                            .identify_dna_variant_records(
                            min_mapping_quality,
                            min_base_quality
                        );
                        variant_records
                    })
                    .flatten()
                    .collect()
            });

            records_map.clear();
            records_map.shrink_to_fit();
            drop(records_map);

            // Filter by chromosome (allow inter-chromosomal translocations)
            let chromosome_id: u16 = chromosome_names_map.get_by_left(chromosome).unwrap().clone();
            variant_records.retain(|vr| {
                if vr.get_variant_type().clone() == VariantType::Translocation {
                    if vr.get_chromosome_1() == chromosome_id || vr.get_chromosome_2() == chromosome_id {
                        true
                    } else {
                        false
                    }
                } else {
                    if vr.get_chromosome_1() == chromosome_id {
                        true
                    } else {
                        // intra-chromosomal variants in other chromosomes should not be called
                        false
                    }
                }
            });

            // Add the variant records to the binary search tree
            for variant_record in variant_records.iter() {
                variant_records_btree
                    .entry(variant_record.get_position_1() as usize)
                    .or_insert_with(HashSet::new)
                    .insert(variant_record.clone());
            }

            // Keep VariantRecord objects whose position_1 is equal to or greater than (prev_end - padding)
            prev_end = if prev_end - padding < 0 { 0 } else { prev_end - padding };
            variant_records_btree = variant_records_btree.split_off(&(prev_end as usize));

            // Deduplicate variant records
            variant_records_btree.values()
                .flat_map(|curr_variant_records| curr_variant_records.iter())
                .for_each(|variant_record| {
                    variant_records.insert(variant_record.clone());
                });

            // Update prev_end
            prev_end = *end as isize;

            // Identify variant calls
            log_info!("\tClustering variant records into variant calls.");
            let variant_calls: Vec<VariantCall> = cluster_variant_records(
                variant_records.into_iter().map(Arc::new).collect(),
                num_threads,
                min_size_proportion,
                max_ins_norm_edit_distance,
                max_intrachromosomal_distance_tau,
                max_intrachromosomal_distance,
                max_interchromosomal_distance,
                false
            );

            // Filter variant calls by the minimum read count and then store
            // them into a variant_call_set
            log_info!("\tFiltering variant calls by the minimum read count.");
            let mut variant_call_set: DNAVariantCallSet = DNAVariantCallSet::new();
            for mut variant_call in variant_calls {
                if variant_call.get_consensus_record().1.len() >= min_reads {
                    // Rename the variant call ID
                    variant_call.id = variant_call_idx;
                    variant_call_set.add_variant_call(variant_call);
                    variant_call_idx += 1;
                }
            }

            // Store variant_call_set in a temp file
            log_info!("\tStoring variant calls into a temp file.");
            let mut temp_file = NamedTempFile::new_in(dir_path).unwrap();
            let mut encoded: Vec<u8> = bincode::serialize(&variant_call_set).expect("Failed to serialize data");
            temp_file.write_all(&encoded).unwrap();
            temp_file.flush().unwrap();
            temp_files.push(temp_file);

            encoded.clear();
            encoded.shrink_to_fit();
            drop(encoded);
            variant_call_set.variant_calls.clear();
            variant_call_set.variant_calls.shrink_to_fit();
            drop(variant_call_set);
        }
        variant_records_btree.clear();
        drop(variant_records_btree);
        pb.inc(1);
    }
    pb.finish_with_message("Completed identifying DNA variants.");
    drop(thread_pool);

    // Step 6. Load all VariantCallSet objects and merge them
    log_info!("Loading all temp files and merging them into a variant call set.");
    let mut variant_calls: HashSet<VariantCall> = HashSet::new();
    for temp_file in temp_files.iter() {
        let mut file = File::open(temp_file.path().to_str().unwrap()).unwrap();
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer).unwrap();
        let variant_call_set_: DNAVariantCallSet = bincode::deserialize(&buffer).expect("Failed to deserialize data");
        buffer.clear();
        buffer.shrink_to_fit();
        for (_,variant_call) in variant_call_set_.variant_calls {
            variant_calls.insert(variant_call);
        }
    }
    let mut variant_call_set: DNAVariantCallSet = DNAVariantCallSet::new();
    let mut variant_call_id: usize = 1;
    for mut variant_call in variant_calls {
        variant_call.id = variant_call_id;
        variant_call_set.add_variant_call(variant_call);
        variant_call_id += 1;
    }
    variant_call_set.load_read_names(read_names_map);
    variant_call_set.load_chromosome_names(chromosome_names_map);

    variant_call_set
}

/// Identify case-specific DNA variants.
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
pub fn identify_case_specific_dna_variants(
    case_bam_file: &str,
    case_bam_bai_file: &str,
    control_bam_files: Vec<&str>,
    control_bam_bai_files: Vec<&str>,
    min_reads: usize,
    min_mapping_quality: u32,
    min_base_quality: u8,
    min_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    max_intrachromosomal_distance_tau: u32,
    max_intrachromosomal_distance: u32,
    max_interchromosomal_distance: u32,
    apply_infinite_sites_assumption: bool,
    num_threads: usize,
    chromosomes: Vec<&str>,
    temp_dir: &str
) -> DNAVariantCallSet {
    assert!(control_bam_files.len() == control_bam_bai_files.len());

    // Step 1. Get all chromosome IDs and names
    let chromosome_names_map: BiMap<Box<str>,u16> = create_chromosome_names_map(case_bam_file);
    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(case_bam_file);

    // Step 2. Make sure the chromosome IDs and names are the same
    for control_bam_file in control_bam_files.iter() {
        let chromosome_names_map_: BiMap<Box<str>,u16> = create_chromosome_names_map(control_bam_file);
        let chromosome_lengths_: HashMap<Box<str>,usize> = get_chromosome_lengths(control_bam_file);
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

    // Step 3. Get read IDs maps
    log_info!("Converting read names to IDs.");
    let case_read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
        case_bam_file,
        case_bam_bai_file,
        num_threads
    );
    let mut control_read_names_map: HashMap<Box<str>,BiMap<Box<str>,usize>> = HashMap::new();
    for (index,control_bam_file) in control_bam_files.iter().enumerate() {
        let control_bam_bai_file: &str = control_bam_bai_files[index];
        let read_names_map: BiMap<Box<str>,usize> = create_read_names_map(
            control_bam_file,
            control_bam_bai_file,
            num_threads
        );
        control_read_names_map.insert(control_bam_file.to_string().into_boxed_str(), read_names_map);
    }

    // Step 4. Fetch the temp directory
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

    // Step 5. Identify case-specific variant calls
    log_info!("Identifying case-specific DNA variants.");
    let mut temp_files: Vec<NamedTempFile> = Vec::new();
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let mut max_distance: u32 = max(max_intrachromosomal_distance_tau, max_intrachromosomal_distance);
    max_distance = max(max_distance, max_interchromosomal_distance);
    let bin_size: u32 = 10_u32.pow((max_distance as f32).log10().floor() as u32 + 1);
    let pb = Arc::new(ProgressBar::new(chromosomes.len() as u64));
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("=>-")
    );
    let mut variant_call_idx: usize = 1;
    for chromosome in chromosomes.iter() {
        let start: usize = 1;
        let end: usize = *chromosome_lengths.get(&chromosome.to_string().into_boxed_str()).unwrap();
        let chromosome_id: u16 = chromosome_names_map.get_by_left(&chromosome.to_string().into_boxed_str()).unwrap().clone();

        log_info!("Identifying variant calls in {}:{}-{}.", chromosome, 1, end);

        // Fetch case BAM records
        let mut case_records_map: HashMap<usize,Vec<bam::Record>> = fetch_bam_records(
            case_bam_file,
            case_bam_bai_file,
            chromosome,
            start,
            end,
            &case_read_names_map,
            num_threads
        );

        // Identify case variant records
        log_info!("\tIdentifying case variant records.");
        let mut case_variant_records_: HashSet<VariantRecord> = thread_pool.install(|| {
            case_records_map
                .par_iter()
                .map(|(read_id, records)| {
                    let read_sequence: Box<str> =  get_fastx_read_sequence(records.iter().collect::<Vec<_>>().as_slice());
                    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records.iter().collect::<Vec<_>>().as_slice());
                    let alignment: Alignment = Alignment::new(
                        *read_id,
                        &*read_sequence,
                        &quality_scores,
                        records
                    );
                    let variant_records: Vec<VariantRecord> = alignment
                        .get_alignment_structure()
                        .identify_dna_variant_records(
                        min_mapping_quality,
                        min_base_quality
                    );
                    variant_records
                })
                .flatten()
                .collect()
        });
        case_records_map.clear();
        case_records_map.shrink_to_fit();
        drop(case_records_map);

        // Filter by chromosome (allow inter-chromosomal translocations)
        case_variant_records_.retain(|vr| {
            if vr.get_variant_type().clone() == VariantType::Translocation {
                if vr.get_chromosome_1() == chromosome_id || vr.get_chromosome_2() == chromosome_id {
                    true
                } else {
                    false
                }
            } else {
                if vr.get_chromosome_1() == chromosome_id {
                    true
                } else {
                    // intra-chromosomal variants in other chromosomes should not be called
                    false
                }
            }
        });

        // Cluster the case variant records
        let variant_calls_: Vec<VariantCall> = cluster_variant_records(
            case_variant_records_.into_iter().map(Arc::new).collect(),
            num_threads,
            min_size_proportion,
            max_ins_norm_edit_distance,
            max_intrachromosomal_distance_tau,
            max_intrachromosomal_distance,
            max_interchromosomal_distance,
            false
        );

        // Filter variant calls by the minimum read count
        let mut case_variant_records: Vec<Arc<VariantRecord>> = Vec::new();
        for variant_call in variant_calls_ {
            let (_, read_ids) = variant_call.get_consensus_record();
            if read_ids.len() >= min_reads {
                for variant_record in variant_call.variant_records {
                    case_variant_records.push(Arc::new(variant_record));
                }
            }
        }

        // Filter out variant records near control variant records
        for (i, control_bam_file) in control_bam_files.iter().enumerate() {
            if case_variant_records.is_empty() == false {
                // Fetch case BAM records
                let control_bam_bai_file: &str = control_bam_bai_files[i];
                let mut control_records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
                    control_bam_file,
                    control_bam_bai_file,
                    chromosome,
                    start,
                    end,
                    &control_read_names_map.get(&control_bam_file.to_string().into_boxed_str()).unwrap(),
                    num_threads
                );

                // Identify control variant records
                log_info!("\tIdentifying control variant records.");
                let mut control_variant_records: HashSet<VariantRecord> = thread_pool.install(|| {
                    control_records_map
                        .par_iter()
                        .map(|(read_id, records)| {
                            let read_sequence: Box<str> = get_fastx_read_sequence(records.iter().collect::<Vec<_>>().as_slice());
                            let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records.iter().collect::<Vec<_>>().as_slice());
                            let alignment: Alignment = Alignment::new(
                                *read_id,
                                &*read_sequence,
                                &quality_scores,
                                records
                            );
                            let variant_records: Vec<VariantRecord> = alignment
                                .get_alignment_structure()
                                .identify_dna_variant_records(
                                min_mapping_quality,
                                min_base_quality
                            );
                            variant_records
                        })
                        .flatten()
                        .collect()
                });
                control_records_map.clear();
                control_records_map.shrink_to_fit();
                drop(control_records_map);

                // Filter by chromosome (allow inter-chromosomal translocations)
                control_variant_records.retain(|vr| {
                    if vr.get_variant_type().clone() == VariantType::Translocation {
                        if vr.get_chromosome_1() == chromosome_id || vr.get_chromosome_2() == chromosome_id {
                            true
                        } else {
                            false
                        }
                    } else {
                        if vr.get_chromosome_1() == chromosome_id {
                            true
                        } else {
                            // intra-chromosomal variants in other chromosomes should not be called
                            false
                        }
                    }
                });

                // Filter out control variants
                log_info!("\tDiffing control variant records.");
                log_info!("\t{} case variant records", case_variant_records.len());
                log_info!("\t{} control variant records", control_variant_records.len());
                case_variant_records = diff_variant_records(
                    case_variant_records,
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
            }
        }

        log_info!("\tCluster variant records into variant calls.");
        let variant_calls: Vec<VariantCall> = cluster_variant_records(
            case_variant_records,
            num_threads,
            min_size_proportion,
            max_ins_norm_edit_distance,
            max_intrachromosomal_distance_tau,
            max_intrachromosomal_distance,
            max_interchromosomal_distance,
            false
        );

        log_info!("\tAdd case-specific variant calls to the set.");
        let mut variant_call_set: DNAVariantCallSet = DNAVariantCallSet::new();
        for mut variant_call in variant_calls {
            if variant_call.get_consensus_record().1.len() >= min_reads {
                // Rename the variant call ID
                variant_call.id = variant_call_idx;
                variant_call_set.add_variant_call(variant_call);
                variant_call_idx += 1;
            }
        }

        // Store variant_call_set in a temp file
        log_info!("\tStoring variant calls into a temp file.");
        let mut temp_file = NamedTempFile::new_in(dir_path).unwrap();
        let mut encoded: Vec<u8> = bincode::serialize(&variant_call_set).expect("Failed to serialize data");
        temp_file.write_all(&encoded).unwrap();
        temp_file.flush().unwrap();
        temp_files.push(temp_file);

        encoded.clear();
        encoded.shrink_to_fit();
        drop(encoded);
        variant_call_set.variant_calls.clear();
        variant_call_set.variant_calls.shrink_to_fit();
        drop(variant_call_set);
        pb.inc(1);
    }
    pb.finish_with_message("Completed identifying case-specific DNA variants.");

    // Step 6. Load all VariantCallSet objects and merge them
    log_info!("Loading all temp files and merging them into a variant call set.");
    let mut variant_calls: HashSet<VariantCall> = HashSet::new();
    for temp_file in temp_files.iter() {
        let mut file = File::open(temp_file.path().to_str().unwrap()).unwrap();
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer).unwrap();
        let variant_call_set_: DNAVariantCallSet = bincode::deserialize(&buffer).expect("Failed to deserialize data");
        buffer.clear();
        buffer.shrink_to_fit();
        for (_, variant_call) in variant_call_set_.variant_calls {
            variant_calls.insert(variant_call);
        }
    }
    let mut variant_call_set: DNAVariantCallSet = DNAVariantCallSet::new();
    let mut variant_call_id: usize = 1;
    for mut variant_call in variant_calls {
        variant_call.id = variant_call_id;
        variant_call_set.add_variant_call(variant_call);
        variant_call_id += 1;
    }
    variant_call_set.load_read_names(case_read_names_map);
    variant_call_set.load_chromosome_names(chromosome_names_map);

    variant_call_set
}
