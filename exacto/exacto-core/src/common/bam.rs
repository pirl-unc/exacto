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
use noodles_bam as bam;
use noodles_bam::bai as bai;
use noodles_bam::bai::Index;
use noodles_core::{Position, Region};
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::{Tag, Value};
use noodles_sam::alignment::record::Flags;
use noodles_sam::Header;
use noodles_sam::alignment::record::Record as SamRecord;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufWriter;
use rayon::ThreadPool;

use crate::prelude::*;


/// Calculates the average base quality score.
///
/// # Arguments
/// * `base_quality_scores`: Base quality scores.
///
/// # Returns
/// * Average base quality score.
pub fn calculate_average_base_quality_score(base_quality_scores: &Vec<u8>) -> f32 {
    assert!(base_quality_scores.is_empty() == false);
    let sum: usize = base_quality_scores.iter().map(|&x| x as usize).sum();
    let count: f32 = base_quality_scores.len() as f32;
    sum as f32 / count
}

/// Creates a BiMap of chromosome names and IDs in a BAM file.
///
/// # Arguments
/// * `bam_file`: BAM file.
///
/// # Returns
/// * BiMap where the left is chromosome name and the right is chromosome ID.
pub fn create_chromosome_names_map(bam_file: &str) -> BiMap<Box<str>, u16> {
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header = reader.read_header().unwrap();
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
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

/// Creates a BiMap of read names and IDs in a BAM file.
///
/// # Arguments
/// * `bam_file`: BAM file.
/// * `bam_bai_file`: BAM.BAI file.
/// * `num_threads`: Number of threads to use.
///
/// # Returns
/// * BiMap where the left is read name and the right is read ID.
pub fn create_read_names_map(
    bam_file: &str,
    bam_bai_file: &str,
    num_threads: usize
) -> BiMap<Box<str>, usize> {
    // Step 1. Split the BAM into regions
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let chromosome_names: Vec<&str> = chromosome_names_map
        .left_values()
        .map(|boxed_str| boxed_str.as_ref())
        .collect();
    let regions: HashMap<Box<str>,Vec<(usize, usize)>> = generate_regions(
        bam_file,
        &chromosome_names,
        10_000_000
    );
    let regions_flattened: Vec<(Box<str>, usize, usize)> = regions
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
    let index: Index = bai::fs::read(bam_bai_file).unwrap();
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
    let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    let mut read_id: usize = 1;
    let mut read_names_: Vec<&Box<str>> = read_names.iter().collect();
    read_names_.sort();
    for read_name in read_names_ {
        read_names_map.insert(read_name.clone(), read_id);
        read_id += 1;
    }
    read_names_map
}

/// Fetches all BAM records in a BAM file.
///
/// # Arguments
/// * `bam_file`: BAM file.
/// * `bam_bai_file`: BAM.BAI file.
/// * `read_names_map`: BiMap where the left is read name and the right is read ID.
/// * `num_threads`: Number of threads.
///
/// # Returns
/// * HashMap where the key is read ID and the value is a vector of noodles_bam::Record objects.
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
    let _index = bai::fs::read(bam_bai_file).unwrap();

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
        records_map
            .entry(read_id)
            .or_insert_with(Vec::new)
            .push(record);
    }

    records_map
}

/// Fetches BAM records in a genomic region.
///
/// # Arguments
/// * `bam_file`: BAM file.
/// * `bam_bai_file`: BAM.BAI file.
/// * `chromosome`: Chromosome name.
/// * `start`: Start position.
/// * `end`: End position.
/// * `read_names_map`: BiMap where the left is read name and the right is read ID.
/// * `num_threads`: Number of threads.
///
/// # Returns
/// * HashMap where the key is read ID and the value is a vector of noodles_bam::Record objects.
pub fn fetch_bam_records(
    bam_file: &str,
    bam_bai_file: &str,
    chromosome: &str,
    start: usize,
    end: usize,
    read_names_map: &BiMap<Box<str>, usize>,
    num_threads: usize
) -> HashMap<usize, Vec<bam::Record>> {
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let chromosome_lengths: HashMap<Box<str>, usize> = get_chromosome_lengths(bam_file);

    // Step 1. Read the BAM header and index
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index = bai::fs::read(bam_bai_file).unwrap();

    // Step 2. Collect primary records
    let start_pos = Position::new(start).unwrap();
    let end_pos = Position::new(end).unwrap();
    let region = Region::new(chromosome, start_pos..=end_pos);
    let primary_records: Vec<bam::Record> = reader
        .query(&header, &index, &region)
        .unwrap()
        .filter_map(Result::ok)
        .collect();

    // Step 3. Identify supplementary regions
    let supplementary_regions: Vec<(Box<str>, usize, usize)> = thread_pool.install(|| {
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

    // Step 6. Combine primary and supplementary records
    let mut combined_records: HashMap<usize,HashMap<(u16, usize, usize, Box<str>),&bam::Record>> = HashMap::new();
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

/// Generates a list of regions with buffer.
///
/// # Arguments
/// * `bam_file`: BAM file.
/// * `chromosomes`: Chromosomes.
/// * `chunk_size`: Chunk size (e.g. 10_000_000).
/// * `chunk_size_buffer`: Chunk size buffer (e.g. 10_000; should be smaller than `chunk_size`).
///
/// # Returns
/// * HashMap where the key is a chromosome name and the value is a vector of tuples (start, end).
pub fn generate_buffered_regions(
    bam_file: &str,
    chromosomes: &Vec<&str>,
    chunk_size: usize,
    chunk_size_buffer: usize
) -> HashMap<Box<str>, Vec<(usize, usize)>> {
    let mut buffered_regions: HashMap<Box<str>,Vec<(usize, usize)>> = HashMap::new();
    let chromosome_lengths: HashMap<Box<str>, usize> = get_chromosome_lengths(bam_file);
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

/// Generates a list of regions.
///
/// # Arguments
/// * `bam_file`: BAM file.
/// * `chromosomes`: Chromosomes.
/// * `chunk_size`: Chunk size (e.g. 10_000_000).
///
/// # Returns
/// * HashMap where the key is a chromosome name and the value is a vector of tuples (start, end).
pub fn generate_regions(
    bam_file: &str,
    chromosomes: &Vec<&str>,
    chunk_size: usize,
) -> HashMap<Box<str>, Vec<(usize, usize)>> {
    let mut regions: HashMap<Box<str>,Vec<(usize, usize)>> = HashMap::new();
    let chromosome_lengths: HashMap<Box<str>, usize> = get_chromosome_lengths(bam_file);
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

/// Gets the alignment end position.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * Alignment end position.
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

/// Gets the alignment mapping quality.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * Mapping quality.
pub fn get_alignment_mapping_quality(record: &bam::Record) -> usize {
    record.mapping_quality().unwrap().get() as usize
}

/// Gets the aligned sequence from the CIGAR string.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * Aligned sequence (the sequence from the original read sequence).
pub fn get_aligned_sequence_from_cigar(record: &bam::Record) -> Box<str> {
    let mut read_pos = 0;
    let mut aligned_sequence = String::new();
    for cigar in record.cigar().iter() {
        let cigar_ = cigar.unwrap();
        match cigar_.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch=> {
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

/// Gets the alignment start position.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * Alignment start position.
pub fn get_alignment_start_position(record: &bam::Record) -> usize {
    record.alignment_start().unwrap().unwrap().get() as usize
}

/// Gets the alignment strand.
///
/// # Arguments
/// - `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * `Strand` enum value indicating the strand of the alignment:
///     * `Strand::Reverse` if the `REVERSE_COMPLEMENTED` flag is set in the BAM record.
///     * `Strand::Forward` otherwise.
pub fn get_alignment_strand(record: &bam::Record) -> Strand {
    let strand: Strand = if record.flags().contains(Flags::REVERSE_COMPLEMENTED) {
        Strand::Reverse
    } else {
        Strand::Forward
    };
    strand
}

/// Gets the BAM header.
///
/// # Arguments
/// * `bam_file`: BAM file.
///
/// # Returns
/// * noodles_sam::Header object.
pub fn get_bam_header(bam_file: &str) -> Header {
    let mut reader = File::open(bam_file).map(bam::io::Reader::new).unwrap();
    let header: Header = reader.read_header().unwrap();
    header
}

/// Gets chromosome lengths in a BAM file.
///
/// # Arguments
/// * `bam_file`: BAM file.
///
/// # Returns
/// * HashMap where the key is a chromosome name and the value is chromosome length.
pub fn get_chromosome_lengths(bam_file: &str) -> HashMap<Box<str>, usize> {
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header = reader.read_header().unwrap();
    let mut chromosome_lengths: HashMap<Box<str>, usize> = HashMap::new();
    for chromosome in header.reference_sequences().iter() {
        chromosome_lengths.insert(chromosome.0.to_string().into_boxed_str(), chromosome.1.length().get() as usize);
    }
    chromosome_lengths
}

/// Gets chromosome names in a BAM file.
///
/// # Arguments
/// * `bam_file`: BAM file.
///
/// # Returns
/// * Vector of chromosome names.
pub fn get_chromosome_names(bam_file: &str) -> Vec<Box<str>> {
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header = reader.read_header().unwrap();
    let mut chromosome_names: Vec<Box<str>> = Vec::new();
    for chromosome in header.reference_sequences().iter() {
        chromosome_names.push(chromosome.0.to_string().into_boxed_str());
    }
    chromosome_names
}

/// Gets the CIGAR operations.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * Vector of tuples (`noodles_sam::alignment::record::cigar::op::Kind`, size).
pub fn get_cigar_operations(record: &bam::Record) -> Vec<(Kind, usize)> {
    record
        .cigar()
        .iter()
        .map(|cigar| {
            let cigar_ = cigar.unwrap();
            (cigar_.kind(), cigar_.len() as usize)
        })
        .collect()
}

/// Gets the CIGAR string.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * CIGAR string.
pub fn get_cigar_string(record: &bam::Record) -> Box<str> {
    let cigar_vec: Vec<(Kind, usize)> = get_cigar_operations(record);
    cigar_vec
        .iter()
        .map(|(kind, len)| format!("{}{}", len,  kind_to_char(*kind)))
        .collect::<Vec<_>>()
        .join("")
        .into()
}

/// Fetches the base quality scores from a list of BAM records.
///
/// # Arguments
/// * `records` - A slice of references to `bam::Record` objects of the same read name.
///
/// # Returns
/// * `Vec<u8>` representing the base quality scores of the first primary
/// (non-supplementary) read in the provided `records`. If the read is marked as
/// reverse complemented, the quality scores are reversed before being returned.
///
/// # Panics
/// This function will panic if no primary (non-supplementary) read is present in
/// the `records` slice.
///
/// # Notes
/// * The `flags` field of the `bam::Record` is used to determine whether a
///   read is supplementary or reverse complemented.
/// * The function processes records sequentially and returns the first valid
///   primary read's quality scores.
/// * Ensure that the provided `records` slice contains reads with valid quality
///   scores to avoid potential runtime errors.
pub fn get_fastx_base_quality_scores(records: &[&bam::Record]) -> Vec<u8> {
    assert!(records.len() > 0, "records must contain at least one record.");
    let read_name: Box<str> = records[0].name().unwrap().to_string().into_boxed_str();
    for &record in records.iter() {
        if read_name != record.name().unwrap().to_string().into_boxed_str() {
            panic!("All records must be from the same read name.");
        }
        if record.flags().is_supplementary() == false {
            let mut quality_scores: Vec<u8> = record.quality_scores().as_ref().to_vec();
            if record.flags().is_reverse_complemented() {
                quality_scores.reverse();
            }
            return quality_scores;
        }
    }
    panic!("Could not find the base quality scores from the primary record.");
}

/// Fetches the original read sequence.
///
/// # Arguments
/// * `records` - A slice of references to `bam::Record` objects of the same read name.
///
/// # Returns
/// * Original read sequence.
///
/// # Panics
/// * The function panics if no non-supplementary record is found in the `records` slice, or if
///   there is an issue retrieving the sequence or record name from the BAM record.
///
/// # Notes
/// * The function assumes that the input records are valid and conform to the expected BAM format.
/// * The method uses `unwrap()` on the sequence conversion and on the record's name retrieval, which
///   suggests potential panics if the string is not valid UTF-8 or if a name is missing.
pub fn get_fastx_read_sequence(records: &[&bam::Record]) -> Box<str> {
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

/// Fetches left soft-clipping information.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * Tuple (is_left_soft_clipped, soft_clip_len).
pub fn get_left_softclipping(record: &bam::Record) -> (bool, usize) {
    let left_soft_clipped: (bool, usize) = record
        .cigar()
        .iter()
        .next()
        .and_then(|op| op.ok()) // Unwrap the Result<Op>
        .filter(|op| op.kind() == Kind::SoftClip)
        .map_or((false, 0), |op| (true, op.len() as usize));
    left_soft_clipped
}

/// Fetches right soft-clipping information.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * Tuple (is_right_soft_clipped, soft_clip_len).
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

/// Fetches the primary alignment base quality scores from the first primary alignment record in the provided list.
///
/// This function iterates through a list of BAM records and returns the base quality scores associated
/// with the first record that is marked as a primary alignment (i.e., not supplementary).
/// If no primary alignment is found, the function will panic.
///
/// # Arguments
/// * `records` - A slice of references to `bam::Record` objects of the same read name.
///
/// # Returns
/// * Base quality scores of the first primary alignment found.
///
/// # Panics
/// This function panics if no primary alignment is found in the list of records. Ensure that the input
/// contains at least one primary alignment record.
pub fn get_primary_alignment_base_quality_scores(records: &[&bam::Record]) -> Vec<u8> {
    for &record in records.iter() {
        if record.flags().is_supplementary() == false {
            let quality_scores: Vec<u8> = record.quality_scores().as_ref().to_vec();
            return quality_scores;
        }
    }
    panic!("Could not find the base quality scores.");
}

/// Fetches the read sequence of the first primary alignment from a list of BAM records.
///
/// This function iterates through a list of BAM records and returns the sequence associated
/// with the first record that is marked as a primary alignment (i.e., not supplementary).
/// If no primary alignment is found, the function will panic.
///
/// # Parameters
/// * `records` - A slice of references to `bam::Record` objects of the same read name.
///
/// # Returns
/// * Read sequence of the first primary alignment found.
///
/// # Panics
/// This function panics if no primary alignment is found in the given list of records.
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

/// Fetches all read names from a BAM file.
///
/// # Arguments
/// * `bam_file`: BAM file.
/// * `bam_bai_file`: BAM.BAI file.
/// * `num_threads`: Number of threads.
///
/// # Returns
/// * HashSet of all unique read names extracted from the BAM file.
pub fn get_read_names(bam_file: &str, bam_bai_file: &str, num_threads: usize) -> HashSet<Box<str>> {
    let mut reader = bam::io::Reader::new(File::open(bam_file).unwrap());
    let _header: Header = reader.read_header().unwrap();
    let _index = bai::fs::read(bam_bai_file).unwrap();
    let records: Vec<bam::Record> = reader.records().map(|r| r.unwrap()).collect();
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let read_names: HashSet<Box<str>> = thread_pool.install(|| {
        records
            .par_iter()
            .filter_map(|record| {
                record
                    .name()
                    .and_then(|n| std::str::from_utf8(n).ok())
                    .map(|s| s.to_owned().into_boxed_str())
            })
            .collect::<HashSet<Box<str>>>()
    });
    read_names
}

/// Fetches read names in a BAM file with at least 1 record that passes the minimum mapping quality (inclusive).
///
/// # Arguments
/// * `bam_file`: BAM file.
/// * `bam_bai_file`: BAM.BAI file.
/// * `num_threads`: Number of threads.
/// * `min_mapping_quality`: Minimum mapping quality (inclusive).
///
/// # Returns
/// * HashSet of read names.
pub fn get_read_names_passing_mapping_quality(
    bam_file: &str,
    bam_bai_file: &str,
    num_threads: usize,
    min_mapping_quality: usize
) -> HashSet<Box<str>> {
    // Step 1. Read BAM header and index
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let _header: Header = reader.read_header().unwrap();
    let _index = bai::fs::read(bam_bai_file).unwrap();

    // Step 2. Read all records into memory
    let records: Vec<bam::Record> = reader.records()
        .filter_map(Result::ok)
        .collect();

    // Step 3. Create a thread pool
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    // Step 4. Identify read names that pass the minimum mapping quality
    // If one record passes the minimum mapping quality,
    // then all records of the same read name pass the minimum mapping quality
    let read_names_passing_mapq: HashSet<Box<str>> = thread_pool.install(|| {
        records.par_iter()
            .filter_map(|record| {
                if record.mapping_quality().unwrap().get() as usize >= min_mapping_quality {
                    record.name().and_then(|name_bytes| {
                        String::from_utf8(name_bytes.to_vec()).ok().map(|s| s.into_boxed_str())
                    })
                } else {
                    None
                }
            })
            .collect()
    });

    read_names_passing_mapq
}

/// Fetches read names that has at least 1 BAM record with a splicing signal
/// based on its CIGAR string or overlaps a single-exon transcript.
///
/// # Arguments
/// * `bam_file`: BAM file.
/// * `bam_bai_file`: BAM.BAI file.
/// * `gene_annotator`: Gene annotator.
/// * `num_threads`: Number of threads.
///
/// # Returns
/// * HashSet of read names that have at least 1 splicing signal or overlaps a single-exon transcript.
pub fn get_read_names_with_splicing(
    bam_file: &str,
    bam_bai_file: &str,
    gene_annotator: &(impl GeneAnnotator + Sync),
    num_threads: usize
) -> HashSet<Box<str>> {
    // Step 1. Get single-exon transcripts
    let mut single_exon_transcripts_map: HashMap<Box<str>, Vec<&Transcript>> = HashMap::new();
    for transcript in gene_annotator.get_transcripts() {
        if transcript.get_exon_ids().len() == 1 {
            single_exon_transcripts_map
                .entry(transcript.chromosome.clone())
                .or_insert_with(Vec::new)
                .push(transcript);
        }
    }

    // Step 2. Create a map of chromosome IDs and names
    let chromosomes_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);

    // Step 3. Read BAM header and index
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let _header: Header = reader.read_header().unwrap();
    let _index = bai::fs::read(bam_bai_file).unwrap();

    // Step 4. Read all records into memory
    let records: Vec<bam::Record> = reader.records()
        .filter_map(Result::ok)
        .collect();

    // Step 5. Identify read IDs that either have splicing or overlap a single-exon transcript
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let read_names_spliced: HashSet<Box<str>> = thread_pool.install(|| {
        records.par_iter()
            .filter_map(|record| {
                if has_splicing(record) {
                    record.name().and_then(|name_bytes| {
                        String::from_utf8(name_bytes.to_vec()).ok().map(|s| s.into_boxed_str())
                    })
                } else {
                    // Check if the read overlaps a single-exon transcript
                    let chromosome_id: usize = record.reference_sequence_id().unwrap().unwrap();
                    let chromosome_name: Box<str> = chromosomes_map.get_by_right(&(chromosome_id as u16)).unwrap().clone();
                    let start: usize = record.alignment_start().unwrap().unwrap().get() - 1;
                    let end: usize = record.alignment_end().unwrap().unwrap().get() - 1;

                    if let Some(transcripts) = single_exon_transcripts_map.get(&chromosome_name) {
                        for &transcript in transcripts {
                            if overlaps(transcript.start as isize, transcript.end as isize, start as isize, end as isize) {
                                return record.name().and_then(|name_bytes| {
                                    String::from_utf8(name_bytes.to_vec()).ok().map(|s| s.into_boxed_str())
                                });
                            }
                        }
                    }

                    None
                }
            })
            .collect()
    });

    read_names_spliced
}

/// Get read sequence.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * Read sequence.
pub fn get_read_sequence(record: &bam::Record) -> Box<str> {
    let s: Vec<u8> = record.sequence().iter().collect();
    let sequence: String = String::from_utf8(s).unwrap();
    sequence.into()
}

/// Retrieves the value of a specified BAM tag as a `String`.
///
/// This function takes a reference to a `bam::Record` and a tag represented as a two-character
/// string slice, and attempts to retrieve the corresponding value within the BAM record's data.
/// If the tag is present and its value is a string, the value is returned encapsulated in
/// a `Box<str>`. Otherwise, the function may return `None` or panic depending on the scenario.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
/// * `tag`: A two-character string that corresponds to the desired BAM tag.
///
/// # Returns
/// * `Some(Box<str>)` if the tag exists in the `bam::Record` and its value is a string.
/// * `None` if the tag is not found in the `bam::Record`.
///
/// # Panics
/// * Panics if the `tag` provided is not exactly two characters long.
/// * Panics if the tag exists, but its value is not a string.
/// * Panics if there is an error retrieving the tag's value.
///
/// # Notes
/// * The function expects the BAM tag to be exactly 2 characters long, as per the BAM format specification.
/// * The function is specifically designed to handle tags with string values. If the tag's value is
///   of a different type, it will panic.
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

/// Determines if a given BAM record has soft clipping in its CIGAR string.
///
/// This function inspects the CIGAR string of a BAM record and checks if any of the operations
/// include soft clipping. Soft clipping is represented by the `Kind::SoftClip` operation in the CIGAR string.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * `true` if the record contains at least one soft-clipped operation in its CIGAR string.
/// * `false` otherwise.
pub fn has_soft_clipping(record: &bam::Record) -> bool {
    record.cigar().iter().any(|op| matches!(op.unwrap().kind(), Kind::SoftClip))
}

/// Determines if a BAM record contains splicing.
///
/// This function checks the CIGAR string of a BAM record to see if it contains
/// any "SKIP" (N) operations, which indicate the presence of splicing in the
/// alignment.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * `true` if the CIGAR string contains a "SKIP" (N) operation, indicating splicing.
/// * `false` otherwise.
pub fn has_splicing(record: &bam::Record) -> bool {
    record.cigar().iter().any(|op| matches!(op.unwrap().kind(), Kind::Skip))
}

/// Checks if a BAM record contains a specific tag.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
/// * `tag`: A string slice representing the tag to search for. It must be exactly 2 characters long.
///
/// # Returns
/// * `true` - If the specified tag is found in the BAM record.
/// * `false` - If the specified tag is not found in the BAM record.
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

/// Converts a `Kind` enum variant to its corresponding single character representation.
///
/// The `Kind` enum represents different types of operations or alignment
/// events for a sequence. This function maps each `Kind` to a specific
/// character that is typically used for representing the event.
///
/// # Arguments
/// * `kind` - A `Kind` enum variant, representing the type of operation.
///
/// # Returns
/// * A `char` corresponding to the given `Kind`, as follows:
///     * `Kind::Match` -> `'M'`: Represents a match between sequences.
///     * `Kind::Insertion` -> `'I'`: Represents an insertion event.
///     * `Kind::Deletion` -> `'D'`: Represents a deletion event.
///     * `Kind::Skip` -> `'N'`: Represents a skipped region in the sequence.
///     * `Kind::SoftClip` -> `'S'`: Represents a soft-clipping event.
///     * `Kind::HardClip` -> `'H'`: Represents a hard-clipping event.
///     * `Kind::Pad` -> `'P'`: Represents a padding operation.
///     * `Kind::SequenceMatch` -> `'='`: Represents a sequence match event.
///     * `Kind::SequenceMismatch` -> `'X'`: Represents a sequence mismatch event.
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
        Kind::SequenceMismatch => 'X'
    }
}

/// Determines if a BAM record is aligned to the reverse complement strand.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * `true` if the record is aligned to the reverse complement strand (indicated by the `REVERSE_COMPLEMENTED` flag).
/// * `false` otherwise.
pub fn is_aligned_to_reverse_strand(record: &bam::Record) -> bool {
    for flag in record.flags() {
        if flag == Flags::REVERSE_COMPLEMENTED {
            return true;
        }
    }
    false
}

/// Writes a BAM file and a BAI file.
///
/// # Arguments
/// * `bam_file`: Output BAM file.
/// * `bai_file`: Output BAM.BAI file.
/// * `header`: Reference to noodles_sam::Header object.
/// * `records`: Reference to a vector of references to noodles_bam::Record objects.
pub fn write_bam_file(
    bam_file: &str,
    bai_file: &str,
    header: &Header,
    records: &Vec<&bam::Record>
) {
    // Step 1. Write BAM file
    let file = File::create(bam_file).unwrap();
    let mut writer = bam::io::Writer::new(BufWriter::new(file));
    writer.write_header(header).unwrap();
    for record in records {
        writer.write_record(header, record).unwrap();
    }
    writer.try_finish().unwrap();

    // Step 2. Explicitly drop the writer to ensure the BAM file is closed
    drop(writer);

    // Step 3: Build index
    let index = bam::fs::index(bam_file).unwrap();

    // Step 4: Save BAI
    bai::fs::write(bai_file, &index).unwrap();
}
