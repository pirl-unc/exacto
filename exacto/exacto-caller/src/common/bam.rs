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
use bstr::{ByteSlice, ByteVec};
use exacto_util::prelude::*;
use noodles_bam as bam;
use noodles_bam::bai;
use noodles_bam::bai::Index;
use noodles_core::{Position, Region};
use noodles_fasta::indexed_reader::Builder;
use noodles_sam::alignment::Record;
use noodles_sam::alignment::record::cigar::{op::Kind};
use noodles_sam::alignment::record::{Cigar,Flags,QualityScores};
use noodles_sam::alignment::record::data::field::{value::Value, Tag};
use noodles_sam::Header;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::rc::Rc;
use sysinfo::System;
use tempfile::NamedTempFile;

use crate::prelude::*;
use crate::log_info;


pub fn capture_memory_usage(message: &str) {
    let mut sys = System::new_all();
    sys.refresh_all();
    let pid = sysinfo::get_current_pid().unwrap();
    if let Some(process) = sys.process(pid) {
        let memory_usage = process.memory();
        let memory_usage_gb = memory_usage as f64 / (1024.0 * 1024.0 * 1024.0);
        log_info!("{}: {:.2} GB", message, memory_usage_gb);
    } else {
        log_info!("Could not get process memory usage");
    }
}

pub fn calculate_average_base_quality_score(vec: &Vec<u8>) -> f32 {
    assert!(vec.is_empty() == false);
    let sum: usize = vec.iter().map(|&x| x as usize).sum();
    let count: f32 = vec.len() as f32;
    sum as f32 / count
}

/// Fetches all reference chromosomes in a BAM file.
///
/// Parameters:
/// - bam_file      :   BAM file.
///
/// Returns:
/// * A BiMap where the left is chromosome name and the right is chromosome ID.
pub fn create_chromosome_names_map(bam_file: &str) -> BiMap<Box<str>,u16> {
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header = reader.read_header().unwrap();
    let mut chromosome_names_map: BiMap<Box<str>,u16> = BiMap::new();
    let mut chromosome_id: u16 = 0;
    for chromosome in header.reference_sequences().iter() {
        chromosome_names_map.insert(chromosome.0.to_string().into_boxed_str(), chromosome_id);
        chromosome_id += 1;
    }
    if chromosome_id > u16::MAX {
        panic!("{} has more than {} chromosomes. Exacto supports up to {} chromosomes (u16).", bam_file, u16::MAX, u16::MAX);
    }
    chromosome_names_map
}

/// Build a BiMap of all read names in a BAM file.
///
/// # Returns:
/// * A BiMap where the left is read name and the right is read ID.
pub fn create_read_names_map(
    bam_file: &str,
    bam_bai_file: &str,
    num_threads: usize
) -> BiMap<Box<str>,usize> {
    // Step 1. Split the BAM into regions
    let chromosome_names_map: BiMap<Box<str>,u16> = create_chromosome_names_map(bam_file);
    let chromosome_names: Vec<&str> = chromosome_names_map
        .left_values()
        .map(|boxed_str| boxed_str.as_ref())
        .collect();
    let regions: HashMap<Box<str>,Vec<(usize,usize)>> = generate_regions(
        bam_file,
        &chromosome_names,
        10_000_000
    );
    let regions_flattened: Vec<(Box<str>,usize,usize)> = regions
        .into_iter()
        .flat_map(|(chromosome, intervals)| {
            intervals
                .into_iter()
                .map(move |(start, end)| (chromosome.clone(), start, end))
        })
        .collect();

    // Step 2. Identify all read names
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::read(bam_bai_file).unwrap();
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let read_names: HashSet<Box<str>> = thread_pool.install(|| {
        regions_flattened
            .par_iter()
            .flat_map(|(chromosome,start,end)| {
                let mut reader = bam::io::reader::Builder::default()
                    .build_from_path(bam_file)
                    .unwrap();
                reader.read_header();
                let start_pos = Position::new(*start as usize).unwrap();
                let end_pos = Position::new(*end as usize).unwrap();
                let region = Region::new(&**chromosome, start_pos..=end_pos);
                let mut local_reader = bam::io::reader::Builder::default()
                    .build_from_path(bam_file)
                    .unwrap();
                let query = local_reader.query(&header, &index, &region).unwrap();
                let mut read_names: HashSet<Box<str>> = HashSet::new();
                for result in query {
                    let record: bam::Record = result.unwrap();
                    read_names.insert(record.name().unwrap().to_string().into());
                }
                read_names
            })
            .collect::<HashSet<Box<str>>>()
    });

    // Step 3. Assign an ID for each read name
    let mut read_names_map: BiMap<Box<str>,usize> = BiMap::new();
    let mut read_id: usize = 1;
    let mut read_names_: Vec<&Box<str>> = read_names.iter().collect();
    read_names_.sort();
    for read_name in read_names_ {
        read_names_map.insert(read_name.clone(), read_id);
        read_id += 1;
    }
    read_names_map
}

pub fn fetch_all_bam_records(
    bam_file: &str,
    bam_bai_file: &str,
    read_names_map: &BiMap<Box<str>, usize>,
    num_threads: usize
) -> HashMap<usize, Vec<bam::Record>> {
    // Step 1. Read BAM header and index
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let _header: Header = reader.read_header().unwrap();
    let _index: Index = bai::read(bam_bai_file).unwrap();

    // Step 2. Read all records into memory
    let records: Vec<bam::Record> = reader.records()
        .filter_map(Result::ok)
        .collect();

    // Step 3. Process in parallel
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let pairs: Vec<(usize, bam::Record)> = thread_pool.install(|| {
        records.into_par_iter()
            .filter_map(|record| {
                record.clone().name().map(|name_bytes| {
                    let name = String::from_utf8(name_bytes.to_vec()).ok()?.into_boxed_str();
                    let read_id = *read_names_map.get_by_left(&name)?;
                    Some((read_id, record))
                }).flatten()
            })
            .collect()
    });

    // Step 4. Group into HashMap
    let mut records_map: HashMap<usize, Vec<bam::Record>> = HashMap::new();
    for (read_id, record) in pairs {
        records_map.entry(read_id).or_insert_with(Vec::new).push(record);
    }

    records_map
}

/// Fetch BAM records.
///
/// # Parameters:
/// * `chromosomes` is a vector of chromosome names.
///
/// # Returns:
/// * A HashMap where the key is read ID and the value is a vector of noodles_bam::Record.
pub fn fetch_bam_records(
    bam_file: &str,
    bam_bai_file: &str,
    chromosome: &str,
    start: usize,
    end: usize,
    read_names_map: &BiMap<Box<str>,usize>,
    num_threads: usize
) -> HashMap<usize,Vec<bam::Record>> {
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);

    // Step 1. Read the BAM header and index
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::read(bam_bai_file).unwrap();

    // Step 2. Collect primary records
    capture_memory_usage("\t[Memory] Before fetching primary records");
    let start_pos = Position::new(start).unwrap();
    let end_pos = Position::new(end).unwrap();
    let region = Region::new(chromosome, start_pos..=end_pos);
    let primary_records: Vec<bam::Record> = reader
        .query(&header, &index, &region)
        .unwrap()
        .filter_map(Result::ok)
        .collect();
    capture_memory_usage("\t[Memory] After fetching primary records");

    // Step 3. Identify supplementary regions
    let supplementary_regions: Vec<(Box<str>,usize,usize)> = thread_pool.install(|| {
        primary_records
            .par_iter()
            .filter_map(|record| {
                if has_tag(record, "SA") {
                    let sa_tags: Box<str> = get_tag_value(record, "SA").unwrap();
                    let sa_tags_split: Vec<&str> = sa_tags.split(';').collect();
                    let mut regions = Vec::new();
                    for sa_tag in sa_tags_split {
                        if sa_tag.is_empty() {
                            continue;
                        }
                        let sa_tag_elements: Vec<&str> = sa_tag.split(',').collect();
                        let sa_chromosome: &str = sa_tag_elements[0];
                        let sa_chromosome_length: isize = *chromosome_lengths.get(sa_chromosome).unwrap() as isize;
                        let sa_position: isize = sa_tag_elements[1].parse().unwrap();
                        let sa_start: isize = if sa_position - 100 > 0 {
                            sa_position - 100
                        } else {
                            1
                        };
                        let sa_end: isize = if sa_position + 100 < sa_chromosome_length {
                            sa_position + 100
                        } else {
                            sa_chromosome_length
                        };
                        regions.push((
                            sa_chromosome.to_string().into_boxed_str(),
                            sa_start as usize,
                            sa_end as usize,
                        ));
                    }
                    Some(regions)
                } else {
                    None
                }
            })
            .flat_map_iter(|regions| regions.into_iter())
            .collect()
    });
    capture_memory_usage("\t[Memory] After identifying supplementary regions");

    // Step 4. Merge supplementary regions
    let mut supplementary_regions_map: HashMap<Box<str>,Vec<(usize,usize)>> = HashMap::new();
    for (sa_chromosome, sa_start, sa_end) in supplementary_regions.iter() {
        supplementary_regions_map
            .entry(sa_chromosome.clone())
            .or_insert_with(Vec::new)
            .push((sa_start.clone(), sa_end.clone()));
    }
    let mut supplementary_regions_merged: Vec<(Box<str>,usize,usize)> = Vec::new();
    for (sa_chromosome, sa_regions) in supplementary_regions_map.iter() {
        let regions: Vec<(isize,isize)> = sa_regions
            .iter()
            .map(|region| {
                let sa_start: isize = region.0 as isize;
                let sa_end: isize = region.1 as isize;
                (sa_start, sa_end)
            })
            .collect();
        let merged_regions: Vec<(isize,isize)> = merge_regions(regions);
        for (sa_start, sa_end) in merged_regions.iter() {
            supplementary_regions_merged.push((sa_chromosome.clone(), *sa_start as usize, *sa_end as usize));
        }
    }

    // Step 4. Identify relevant read IDs
    let mut primary_record_read_ids: HashSet<usize> = HashSet::new();
    for record in primary_records.iter() {
        let read_name: Box<str> = record.name().unwrap().to_string().into_boxed_str();
        let read_index: usize = *read_names_map.get_by_left(&read_name).unwrap();
        primary_record_read_ids.insert(read_index);
    }
    capture_memory_usage("\t[Memory] After getting primary record IDs");

    // Step 5. Fetch supplementary records
    let supplementary_records: Vec<bam::Record> = thread_pool.install(|| {
        supplementary_regions_merged
            .par_iter()
            .flat_map(|(sa_chromosome, sa_start, sa_end)| {
                let mut local_reader = bam::io::reader::Builder::default()
                    .build_from_path(bam_file)
                    .unwrap();
                let start_pos = Position::new(*sa_start as usize).unwrap();
                let end_pos = Position::new(*sa_end as usize).unwrap();
                let region = Region::new(&**sa_chromosome, start_pos..=end_pos);
                let query = local_reader.query(&header, &index, &region).unwrap();
                query
                    .filter_map(|result| {
                        let record: bam::Record = result.unwrap();
                        if record.flags().is_unmapped() || record.flags().is_secondary() {
                            return None;
                        }
                        let read_name: Box<str> = record.name().unwrap().to_string().into_boxed_str();
                        let read_index: usize = *read_names_map.get_by_left(&read_name).unwrap();
                        if primary_record_read_ids.contains(&read_index) {
                            Some(record)
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect()
    });
    capture_memory_usage("\t[Memory] After getting supplementary record");

    // Step 6. Combine primary and supplementary records
    let mut combined_records: HashMap<usize,HashMap<(u16,usize,usize,Box<str>),&bam::Record>> = HashMap::new();
    for record in primary_records.iter().chain(supplementary_records.iter()) {
        let read_name: Box<str> = record.name().unwrap().to_string().into_boxed_str();
        let read_index: usize = *read_names_map.get_by_left(&read_name).unwrap();
        let key: (u16,usize,usize,Box<str>) = (
            record.flags().bits(),
            get_alignment_start_position(&record),
            get_alignment_end_position(&record),
            get_cigar_string(&record),
        );
        combined_records
            .entry(read_index)
            .or_insert_with(HashMap::new)
            .insert(key,record);
    }
    capture_memory_usage("\t[Memory] After combining records");

    // Step 7. Prepare output
    thread_pool.install(|| {
        combined_records
            .par_iter()
            .map(|(read_index, record_map)| {
                let mut records: Vec<bam::Record> = record_map
                    .values()
                    .map(|v| (*v).clone())
                    .collect();
                records.sort_by(|a, b| {
                    get_alignment_start_position(a).cmp(&get_alignment_start_position(b))
                });
                (*read_index, records)
            })
            .collect::<HashMap<usize,Vec<bam::Record>>>()
    })
}

/// Generate a list of regions.
///
/// Parameters:
/// - bam_file              :   BAM file.
/// - chromosomes           :   Chromosomes.
/// - chunk_size            :   Chunk size (e.g. 10_000_000).
///
/// Returns:
/// - Vec<(chromosome,start,end)>
pub fn generate_regions(
    bam_file: &str,
    chromosomes: &Vec<&str>,
    chunk_size: usize,
) -> HashMap<Box<str>,Vec<(usize,usize)>> {
    let mut regions: HashMap<Box<str>,Vec<(usize,usize)>> = HashMap::new();
    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    for chromosome in chromosomes.iter() {
        let chromosome_length: usize = *chromosome_lengths.get(&chromosome.to_string().into_boxed_str()).unwrap();
        // Divide the chromosome into chunks with overlaps
        let mut start: usize = 0;
        while start < chromosome_length {
            start += 1;

            // Compute the end of the current region
            let mut end = start + chunk_size - 1;
            if end > chromosome_length {
                end = chromosome_length;
            }

            // Add the region to the list
            regions
                .entry(chromosome.to_string().into_boxed_str())
                .or_insert_with(Vec::new)
                .push((start,end));

            // Move to the next chunk
            start = end;
        }
    }
    regions
}

/// Generate a list of regions with buffer.
///
/// Parameters:
/// - bam_file              :   BAM file.
/// - chromosomes           :   Chromosomes.
/// - chunk_size            :   Chunk size (e.g. 10_000_000).
/// - chunk_size_buffer     :   Chunk size buffer (e.g. 10_000; should be smaller than `chunk_size`).
///
/// Returns:
/// - Vec<(chromosome,start,end)>
pub fn generate_buffered_regions(
    bam_file: &str,
    chromosomes: &Vec<&str>,
    chunk_size: usize,
    chunk_size_buffer: usize
) -> HashMap<Box<str>,Vec<(usize,usize)>> {
    let mut buffered_regions: HashMap<Box<str>,Vec<(usize,usize)>> = HashMap::new();
    let chromosome_lengths: HashMap<Box<str>,usize> = get_chromosome_lengths(bam_file);
    for chromosome in chromosomes.iter() {
        let chromosome_length: usize = *chromosome_lengths.get(&chromosome.to_string().into_boxed_str()).unwrap();
        // Divide the chromosome into chunks with overlaps
        let mut start: usize = 0;
        while start < chromosome_length {
            // Compute the end of the current region
            let mut end = start + chunk_size;
            if end > chromosome_length {
                end = chromosome_length;
            }

            // Add buffer to the start and end of the region
            let buffered_start: usize = if start < chunk_size_buffer { 0 } else { start - chunk_size_buffer };
            let buffered_end: usize = if end + chunk_size_buffer > chromosome_length { chromosome_length } else { end + chunk_size_buffer };

            // Add the region to the list
            buffered_regions
                .entry(chromosome.to_string().into_boxed_str())
                .or_insert_with(Vec::new)
                .push((buffered_start,buffered_end));

            // Move to the next chunk
            start = end;
        }
    }
    buffered_regions
}

/// Get alignment end position.
pub fn get_alignment_end_position(record: &bam::Record) -> usize {
    let alignment_start_position: usize = record.alignment_start().unwrap().unwrap().get() as usize;
    let alignment_span: usize = record
        .cigar()
        .iter()
        .filter_map(|op| op.ok()) // Unwrap the Result<Op, Error> safely
        .filter(|op| matches!(op.kind(), Kind::Match | Kind::Deletion | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Skip))
        .map(|op| op.len())
        .sum();
    let alignment_last_position: usize = alignment_start_position + alignment_span as usize - 1;
    alignment_last_position
}

/// Get aligned sequence from CIGAR string.
///
/// Returns:
/// - Aligned sequence (the sequence from the original read sequence).
pub fn get_aligned_sequence_from_cigar(record: &bam::Record) -> Box<str> {
    let mut read_pos = 0;
    let mut aligned_sequence = String::new();
    for cigar in record.cigar().iter() {
        let cigar_ = cigar.unwrap();
        match cigar_.kind() {
            Kind::SequenceMatch | Kind::SequenceMismatch=> {
                for _ in 0..cigar_.len() {
                    if let Some(base) = record.sequence().get(read_pos) {
                        aligned_sequence.push(char::from(base));
                    }
                    read_pos += 1;
                }
            },
            Kind::SoftClip => {
                read_pos += cigar_.len();
            },
            Kind::Deletion |Kind::Skip => {
            },
            Kind::Insertion => {
                for _ in 0..cigar_.len() {
                    if let Some(base) = record.sequence().get(read_pos) {
                        aligned_sequence.push(char::from(base));
                    }
                    read_pos += 1;
                }
            },
            _ => {
                panic!("Unknown cigar type: {:?}", cigar_.kind());
            }
        }
    }

    if record.flags().is_reverse_complemented() {
        reverse_complement(&aligned_sequence).into()
    } else {
        aligned_sequence.into()
    }
}

/// Get alignment start position.
pub fn get_alignment_start_position(record: &bam::Record) -> usize {
    record.alignment_start().unwrap().unwrap().get() as usize
}

/// Get base quality scores for the aligned part.
pub fn get_alignment_base_quality_scores(record: &bam::Record) -> Vec<u8> {
    let quality_scores: Vec<u8> = record.quality_scores().as_ref().to_vec();
    let start: usize;
    let end: usize;
    let left_softclipping: (bool,usize) = get_left_softclipping(record);
    let right_softclipping: (bool,usize) = get_right_softclipping(record);
    if left_softclipping.0 {
        start = left_softclipping.1 as usize;
    } else {
        start = 0;
    }
    if right_softclipping.0 {
        end = quality_scores.len() - 1 - right_softclipping.1 as usize;
    } else {
        end = quality_scores.len() - 1;
    }
    quality_scores[start..=end].to_vec()
}

/// Get base quality scores
pub fn get_base_quality_scores(record: &bam::Record) -> Vec<u8> {
    record.quality_scores().as_ref().to_vec()
}

/// Fetches all reference chromosome names in a BAM file.
///
/// Parameters:
/// - bam_file      :   BAM file.
///
/// Returns:
/// - Vec<String>
pub fn get_chromosome_names(bam_file: &str) -> Vec<Box<str>> {
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header = reader.read_header().unwrap();
    let mut chromosome_names: Vec<Box<str>> = Vec::new();
    for chromosome in header.reference_sequences().iter() {
        chromosome_names.push(chromosome.0.to_string().into_boxed_str());
    }
    chromosome_names
}

/// Fetches all reference chromosomes in a BAM file.
///
/// Parameters:
/// - bam_file      :   BAM file.
///
/// Returns:
/// - HashMap<chromosome name,chromosome length>
pub fn get_chromosome_lengths(bam_file: &str) -> HashMap<Box<str>,usize> {
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header = reader.read_header().unwrap();
    let mut chromosome_lengths: HashMap<Box<str>,usize> = HashMap::new();
    for chromosome in header.reference_sequences().iter() {
        chromosome_lengths.insert(chromosome.0.to_string().into_boxed_str(), chromosome.1.length().get() as usize);
    }
    chromosome_lengths
}

pub fn get_cigar_operations(record: &bam::Record) -> Vec<(Kind,usize)> {
    record
        .cigar()
        .iter()
        .map(|cigar| {
            let cigar_ = cigar.unwrap();
            (cigar_.kind(), cigar_.len() as usize)
        })
        .collect()
}

pub fn get_cigar_string(record: &bam::Record) -> Box<str> {
    let cigar_vec: Vec<(Kind,usize)> = get_cigar_operations(record);
    cigar_vec
        .iter()
        .map(|(kind, len)| format!("{}{}", len,  kind_to_char(*kind)))
        .collect::<Vec<_>>()
        .join("")
        .into()
}

/// Get left soft-clipping.
///
/// Returns:
/// - (true if left is soft-clipped, number of soft-clipped bases)
pub fn get_left_softclipping(record: &bam::Record) -> (bool,usize) {
    let left_soft_clipped: (bool, usize) = record
        .cigar()
        .iter()
        .next()
        .and_then(|op| op.ok()) // Unwrap the Result<Op>
        .filter(|op| op.kind() == Kind::SoftClip)
        .map_or((false, 0), |op| (true, op.len() as usize));
    left_soft_clipped
}

/// Get right soft-clipping.
///
/// Returns:
/// - (true if right is soft-clipped, number of soft-clipped bases)
pub fn get_right_softclipping(record: &bam::Record) -> (bool, usize) {
    let right_soft_clipped: (bool, usize) = record
        .cigar()
        .iter()
        .last()
        .and_then(|op| op.ok()) // Unwrap the Result<Op>
        .filter(|op| op.kind() == Kind::SoftClip)
        .map_or((false, 0), |op| (true, op.len() as usize));
    right_soft_clipped
}

/// Get original FASTQ base quality scores.
pub fn get_original_base_quality_scores(records: &[&bam::Record]) -> Vec<u8> {
    for &record in records.iter() {
        if record.flags().is_supplementary() == false {
            let mut quality_scores: Vec<u8> = record.quality_scores().as_ref().to_vec();
            // let mut quality_scores: Vec<u8> = record.quality_scores().iter().collect();
            if record.flags().is_reverse_complemented() {
                quality_scores.reverse();
            }
            return quality_scores;
        }
    }
    panic!("Could not find the base quality scores.");
}

/// Get original FASTQ read sequence.
pub fn get_original_read_sequence(records: &[&bam::Record]) -> Box<str> {
    for &record in records.iter() {
        if record.flags().is_supplementary() == false {
            let s: Vec<u8> = record.sequence().iter().collect();
            let mut sequence: Box<str> = String::from_utf8(s).unwrap().into_boxed_str();
            if record.flags().is_reverse_complemented() {
                sequence = reverse_complement(&*sequence);
            }
            return sequence;
        }
    }
    panic!("Could not find the read sequence for {}.", records[0].name().unwrap());
}

/// Get primary alignment base quality scores.
pub fn get_primary_alignment_base_quality_scores(records: &[&bam::Record]) -> Vec<u8> {
    for &record in records.iter() {
        if record.flags().is_supplementary() == false {
            let quality_scores: Vec<u8> = record.quality_scores().as_ref().to_vec();
            // let quality_scores: Vec<u8> = record.quality_scores().iter().collect();
            return quality_scores;
        }
    }
    panic!("Could not find the base quality scores.");
}


/// Get primary alignment read sequence.
pub fn get_primary_alignment_read_sequence(records: &[&bam::Record]) -> Box<str> {
    for &record in records.iter() {
        if record.flags().is_supplementary() == false {
            let s: Vec<u8> = record.sequence().iter().collect();
            let sequence: String = String::from_utf8(s).unwrap();
            return sequence.into();
        }
    }
    panic!("Could not find the read sequence.");
}

/// Get read sequence.
///
/// Returns:
/// - (read_sequence, read_reverse_complemented). The `read_sequence` is the original read sequence.
/// If `read_reverse_complemented` is true, then the original read was reverse-complemented to
/// align to the reference.
pub fn get_read_sequence(record: &bam::Record) -> Box<str> {
    let s: Vec<u8> = record.sequence().iter().collect();
    let sequence: String = String::from_utf8(s).unwrap();
    sequence.into()
}

/// Get tag value.
pub fn get_tag_value(record: &bam::Record, tag: &str) -> Option<Box<str>> {
    let tag_bytes = tag.as_bytes();
    if tag_bytes.len() != 2 {
        panic!("Tag must be exactly 2 characters.");
    }
    let tag_array: [u8; 2] = [tag_bytes[0], tag_bytes[1]];
    let tag = Tag::from(tag_array);
    match record.data().get(&tag) {
        Some(Ok(value)) => {
            match value {
                Value::String(s) => Some(s.to_string().into()),
                _ => {
                    panic!("Tag is not a string.");
                }
            }
        },
        Some(Err(_)) => {
            panic!("Could not fetch the tag value.");
        },
        None => None,
    }
}

pub fn has_soft_clipping(record: &bam::Record) -> bool {
    record.cigar().iter().any(|op| matches!(op.unwrap().kind(), Kind::SoftClip))
}

pub fn has_tag(record: &bam::Record, tag: &str) -> bool {
    let tag_bytes = tag.as_bytes();
    if tag_bytes.len() != 2 {
        panic!("Tag must be exactly 2 characters.");
    }
    let tag_array: [u8; 2] = [tag_bytes[0], tag_bytes[1]];
    let tag = Tag::from(tag_array);
    match record.data().get(&tag) {
        Some(Ok(value)) => {
            true
        },
        Some(Err(_)) => {
            panic!("Could not fetch the tag value.");
        },
        None => {
            false
        }
    }
}

pub fn kind_to_char(kind: Kind) -> char {
    match kind {
        Kind::Match => 'M',
        Kind::Insertion => 'I',
        Kind::Deletion => 'D',
        Kind::Skip => 'N',
        Kind::SoftClip => 'S',
        Kind::HardClip => 'H',
        Kind::Pad => 'P',
        Kind::SequenceMatch => '=',
        Kind::SequenceMismatch => 'X',
    }
}

pub fn is_aligned_to_reverse_strand(record: &bam::Record) -> bool {
    for flag in record.flags() {
        if flag == Flags::REVERSE_COMPLEMENTED {
            return true;
        }
    }
    false
}
