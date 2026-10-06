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
use bstr::ByteSlice;
use noodles_bam as bam;
use noodles_bam::bai as bai;
use noodles_bam::bai::Index;
use noodles_bgzf as bgzf;
use noodles_bgzf::io::{BufRead, Seek};
use noodles_bgzf::VirtualPosition;
use noodles_core::{Position, Region};
use noodles_util::alignment::{io, iter::Depth};
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::{Tag, Value};
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record::Record;
use noodles_sam::alignment::record::cigar::Op;
use noodles_sam::Header;
use rayon::prelude::*;
use rayon::ThreadPool;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, SeekFrom};
use std::io::{Read as _, Seek as _};
use std::num::NonZeroUsize;
use std::sync::Arc;

use crate::log_info;
use crate::prelude::*;


/// The empty block that ends every complete BGZF file (SAM/BAM specification, section 4.1.2).
const BGZF_EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
    0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
];


/// Calculates the average base quality score.
///
/// # Arguments
/// * `base_quality_scores`: Base quality scores.
///
/// # Returns
/// * Average base quality score.
pub fn calculate_average_base_quality_score(base_quality_scores: &Vec<BaseQuality>) -> f32 {
    assert!(base_quality_scores.is_empty() == false);
    let sum: usize = base_quality_scores.iter().map(|&x| x as usize).sum();
    let count: f32 = base_quality_scores.len() as f32;
    sum as f32 / count
}


/// Checks that a BAM file ends with the BGZF end-of-file block.
///
/// A BAM cut short (an interrupted copy, a full disk) otherwise reads as a complete file that
/// stops at the cut: noodles returns 0 bytes there, as it does at the real end. Every function
/// here that opens a BAM file by path calls this first.
///
/// # Arguments
/// * `bam_file`: BAM file.
///
/// # Panics
/// * If the file cannot be read or does not end with the end-of-file block.
pub fn check_bam_end_of_file(bam_file: &str) {
    let mut file: File = File::open(bam_file).unwrap_or_else(|e| panic!("Could not open {bam_file}: {e}"));
    let mut last_block: [u8; 28] = [0; 28];
    let has_eof_block: bool = file.seek(SeekFrom::End(-(BGZF_EOF_BLOCK.len() as i64))).is_ok()
        && file.read_exact(&mut last_block).is_ok()
        && last_block == BGZF_EOF_BLOCK;
    assert!(
        has_eof_block,
        "{bam_file} does not end with the BGZF end-of-file block; the file is probably truncated."
    );
}


/// Creates a BiMap of chromosome names and IDs in a BAM file.
///
/// # Arguments
/// * `bam_file`: BAM file.
///
/// # Returns
/// * BiMap where the left is chromosome name and the right is chromosome ID.
pub fn create_chromosome_names_map(bam_file: &str) -> BiMap<ReferenceChromosomeName, ReferenceChromosomeID> {
    check_bam_end_of_file(bam_file);
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header = reader.read_header().unwrap();
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    let mut chromosome_id: u32 = 0;
    for chromosome in header.reference_sequences().iter() {
        // chromosome_names_map.insert(chromosome.0.to_string().into_boxed_str(), chromosome_id);
        if chromosome_id > u16::MAX as u32 {
            panic!("{} has more than {} chromosomes. Exacto supports up to {} chromosomes (u16).", bam_file, u16::MAX, u16::MAX);
        }
        chromosome_names_map.insert(chromosome.0.to_string().into_boxed_str(), chromosome_id as u16);
        chromosome_id += 1;
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
) -> BiMap<ReadName, ReadID> {
    check_bam_end_of_file(bam_file);

    // Step 1. Split the BAM into regions
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let chromosome_names: Vec<ReferenceChromosomeName> = chromosome_names_map
        .left_values()
        .cloned()
        .collect();
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let regions: HashMap<Box<str>,Vec<(u32, u32)>> = generate_regions(
        &chromosome_names,
        &chromosome_lengths,
        10_000_000
    );
    let regions_flattened: Vec<(Box<str>, u32, u32)> = regions
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
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
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
                let start_pos: Position = Position::new(*start as usize).unwrap();
                let end_pos: Position = Position::new(*end as usize).unwrap();
                let region: Region = Region::new(&**chromosome, start_pos..=end_pos);
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
    read_names_map: &BiMap<ReadName, ReadID>,
    num_threads: usize
) -> HashMap<ReadID, Vec<bam::Record>> {
    check_bam_end_of_file(bam_file);

    // Step 1. Read BAM header and index
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header_: Header = reader.read_header().unwrap();
    let index_: Index = bai::fs::read(bam_bai_file).unwrap();

    // Step 2. Read all records into memory. Secondary records are alternative mappings of
    // the same read bases and corrupt the merged alignment model — see `index_bam_records`.
    let records: Vec<bam::Record> = reader.records()
        .map(|result| result.unwrap_or_else(|e| panic!("Could not read a record of {bam_file}: {e}")))
        .filter(|record| !record.flags().is_secondary())
        .collect();

    // Step 3. Process in parallel
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
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
pub fn fetch_bam_records<R>(
    reader: &mut bam::io::Reader<R>,
    header: &Header,
    index: &Index,
    chromosome: ReferenceChromosomeName,
    start: ReferencePosition,
    end: ReferencePosition,
    record_positions_map: &HashMap<ReadID, Vec<VirtualPosition>>,
    read_names_map: &BiMap<ReadName, ReadID>,
    max_records: ReadSupport,
    num_threads: usize
) -> HashMap<ReadID, Vec<bam::Record>>
where
    R: BufRead + Seek
{
    // Step 1. Collect the records in the region. Secondary records are left out, as in
    // `index_bam_records`: they belong to no read's records, so they bring no read into the region.
    let start_pos: Position = Position::new(start as usize).unwrap();
    let end_pos: Position = Position::new(end as usize).unwrap();
    let region: Region = Region::new(&*chromosome, start_pos..=end_pos);
    let primary_records: Vec<bam::Record> = reader
        .query(header, index, &region)
        .unwrap()
        .map(|result| result.unwrap_or_else(|e| panic!("Could not read a record in {chromosome}:{start}-{end}: {e}")))
        .filter(|record| !record.flags().is_secondary())
        .collect();

    // Step 2. Identify relevant read IDs. A read that `index_bam_records` left out (its records
    // in the file are all unmapped, or all supplementary) is skipped.
    let mut read_ids: HashSet<usize> = HashSet::new();
    for record in primary_records.iter() {
        let read_name: Box<str> = record.name().unwrap().to_string().into_boxed_str();
        if let Some(read_id) = read_names_map.get_by_left(&read_name) {
            read_ids.insert(*read_id);
        }
    }

    // Step 3. Fetch all BAM records
    let mut read_records: HashMap<usize, HashMap<(u16, u32, u32, Box<str>), bam::Record>> = HashMap::new();
    let mut record: bam::Record = bam::Record::default();
    for read_id in read_ids.iter() {
        let vps: &Vec<VirtualPosition> = record_positions_map.get(read_id).unwrap();
        if vps.len() <= max_records as usize {
            for &vp in vps {
                reader.get_mut().seek_to_virtual_position(vp).unwrap();
                let bytes_read = reader.read_record(&mut record).unwrap();
                if bytes_read == 0 {
                    continue;
                }
                let curr_read_name: Box<str> = record.name().unwrap().to_string().into_boxed_str();
                let curr_read_id: usize = *read_names_map.get_by_left(&curr_read_name).unwrap();
                if curr_read_id == *read_id {
                    let key: (u16, u32, u32, Box<str>) = (
                        record.flags().bits(),
                        get_alignment_start_position(&record),
                        get_alignment_end_position(&record),
                        get_cigar_string(&record),
                    );
                    read_records
                        .entry(curr_read_id)
                        .or_insert_with(HashMap::new)
                        .insert(key, record.clone());
                }
            }
        }
    }

    // Step 4. Sort BAM records, then output
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    thread_pool.install(|| {
        read_records
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
            .collect::<HashMap<usize, Vec<bam::Record>>>()
    })
}


/// Fetches all BAM records for one read ID by seeking to its indexed virtual positions.
///
/// `record_positions_map` comes from `index_bam_records`: it maps a read ID to the
/// virtual positions of all that read's records (primary, secondary, and supplementary),
/// wherever they align — so this gathers alignments a coordinate `query` would miss.
pub fn fetch_bam_records_for_read_id<R>(
    reader: &mut bam::io::Reader<R>,
    read_id: ReadID,
    record_positions_map: &HashMap<ReadID, Vec<VirtualPosition>>
) -> Vec<bam::Record>
where
    R: BufRead + Seek
{
    let mut records: Vec<bam::Record> = Vec::new();
    if let Some(vps) = record_positions_map.get(&read_id) {
        let mut record: bam::Record = bam::Record::default();
        for &vp in vps {
            reader.get_mut().seek_to_virtual_position(vp).unwrap();
            let bytes_read: usize = reader.read_record(&mut record).unwrap();
            if bytes_read == 0 {
                continue;
            }
            records.push(record.clone());
        }
    }
    records
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
    chromosomes: &Vec<ReferenceChromosomeName>,
    chunk_size: u32,
    chunk_size_buffer: u32
) -> HashMap<Box<str>, Vec<(u32, u32)>> {
    assert!(chunk_size > 0, "chunk_size must be > 0");
    let mut buffered_regions: HashMap<Box<str>,Vec<(u32, u32)>> = HashMap::new();
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    for chromosome in chromosomes.iter() {
        let chromosome_length: u32 = *chromosome_lengths.get(&chromosome.to_string().into_boxed_str()).unwrap();
        // Divide the chromosome into chunks with overlaps
        let mut start: u32 = 0;
        while start < chromosome_length {
            // Compute the end of the current region
            let mut end: u32 = start + chunk_size;
            if end > chromosome_length {
                end = chromosome_length;
            }

            // Add buffer to the start and end of the region
            let buffered_start: u32 = if start < chunk_size_buffer { 0 } else { start - chunk_size_buffer };
            let buffered_end: u32 = if end + chunk_size_buffer > chromosome_length {
                chromosome_length
            } else {
                end + chunk_size_buffer
            };

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
    chromosomes: &Vec<ReferenceChromosomeName>,
    chromosome_lengths: &HashMap<ReferenceChromosomeName, ReferenceChromosomeLength>,
    chunk_size: u32
) -> HashMap<ReferenceChromosomeName, Vec<(ReferencePosition, ReferencePosition)>> {
    assert!(chunk_size > 0, "chunk_size must be > 0");
    let mut regions: HashMap<Box<str>,Vec<(u32, u32)>> = HashMap::new();
    for chromosome in chromosomes.iter() {
        let chromosome_length: u32 = *chromosome_lengths.get(&chromosome.to_string().into_boxed_str()).unwrap();
        // Divide the chromosome into chunks with overlaps
        let mut start: u32 = 0;
        while start < chromosome_length {
            start += 1;

            // Compute the end of the current region
            let mut end: u32 = start + chunk_size - 1;
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
pub fn get_alignment_end_position(record: &bam::Record) -> u32 {
    let alignment_start_position: usize = record.alignment_start().unwrap().unwrap().get();
    let alignment_span: usize = record
        .cigar()
        .iter()
        .filter_map(|op| op.ok()) // Unwrap the Result<Op, Error> safely
        .filter(|op| matches!(op.kind(), Kind::Match | Kind::Deletion | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Skip))
        .map(|op| op.len())
        .sum();
    let alignment_last_position: usize = alignment_start_position + alignment_span - 1;
    alignment_last_position as u32
}

/// Gets the alignment mapping quality.
///
/// The SAM specification reserves 255 for "not available" (STAR writes it for every unique
/// mapper), and noodles reads it as `None`. It is returned as 255, the value in the file, as
/// samtools and pysam do, so such a record passes any minimum mapping quality.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * Mapping quality.
pub fn get_alignment_mapping_quality(record: &bam::Record) -> u16 {
    record.mapping_quality().map_or(255, |mapping_quality| mapping_quality.get() as u16)
}

/// Gets the aligned sequence from the CIGAR string.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * Aligned sequence (the sequence from the original read sequence).
pub fn get_aligned_sequence_from_cigar(record: &bam::Record) -> Box<str> {
    let mut read_pos: usize = 0;
    let mut aligned_sequence: String = String::new();
    for cigar in record.cigar().iter() {
        let cigar_: Op = cigar.unwrap();
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
            // A hard clip or a pad holds no base of the stored sequence.
            Kind::Deletion | Kind::Skip | Kind::HardClip | Kind::Pad => {
            },
            Kind::Insertion => {
                for _ in 0..cigar_.len() {
                    if let Some(base) = record.sequence().get(read_pos) {
                        aligned_sequence.push(char::from(base));
                    }
                    read_pos += 1;
                }
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
pub fn get_alignment_start_position(record: &bam::Record) -> u32 {
    record.alignment_start().unwrap().unwrap().get() as u32
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
    check_bam_end_of_file(bam_file);
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
pub fn get_chromosome_lengths(bam_file: &str) -> HashMap<Box<str>, u32> {
    check_bam_end_of_file(bam_file);
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header: Header = reader.read_header().unwrap();
    let mut chromosome_lengths: HashMap<Box<str>, u32> = HashMap::new();
    for chromosome in header.reference_sequences().iter() {
        chromosome_lengths.insert(chromosome.0.to_string().into_boxed_str(), chromosome.1.length().get() as u32);
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
    check_bam_end_of_file(bam_file);
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header: Header = reader.read_header().unwrap();
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
pub fn get_cigar_operations(record: &dyn Record) -> Vec<(Kind, u32)> {
    record
        .cigar()
        .iter()
        .map(|cigar| {
            let cigar_ = cigar.unwrap();
            (cigar_.kind(), cigar_.len() as u32)
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
    let cigar_vec: Vec<(Kind, u32)> = get_cigar_operations(record);
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
/// the `records` slice. `index_bam_records` leaves out reads that have no such record.
///
/// # Notes
/// * The `flags` field of the `bam::Record` is used to determine whether a
///   read is supplementary or reverse complemented.
/// * The function processes records sequentially and returns the first valid
///   primary read's quality scores.
/// * Ensure that the provided `records` slice contains reads with valid quality
///   scores to avoid potential runtime errors.
pub fn get_bam_fastx_base_quality_scores(records: &Vec<bam::Record>) -> Vec<u8> {
    assert!(records.len() > 0, "records must contain at least one record.");
    let read_name: Box<str> = records[0].name().unwrap().to_string().into_boxed_str();
    for record in records.iter() {
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
///   `index_bam_records` leaves out reads that have no non-supplementary record.
///
/// # Notes
/// * The function assumes that the input records are valid and conform to the expected BAM format.
/// * The method uses `unwrap()` on the sequence conversion and on the record's name retrieval, which
///   suggests potential panics if the string is not valid UTF-8 or if a name is missing.
pub fn get_bam_fastx_read_sequence(records: &Vec<bam::Record>) -> Box<str> {
    for record in records.iter() {
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
pub fn get_left_softclipping(record: &bam::Record) -> (bool, u32) {
    let left_soft_clipped: (bool, u32) = record
        .cigar()
        .iter()
        .next()
        .and_then(|op| op.ok()) // Unwrap the Result<Op>
        .filter(|op| op.kind() == Kind::SoftClip)
        .map_or((false, 0), |op| (true, op.len() as u32));
    left_soft_clipped
}

/// Fetches right soft-clipping information.
///
/// # Arguments
/// * `record`: Reference to a `noodles_bam::Record` object.
///
/// # Returns
/// * Tuple (is_right_soft_clipped, soft_clip_len).
pub fn get_right_softclipping(record: &bam::Record) -> (bool, u32) {
    let right_soft_clipped: (bool, u32) = record
        .cigar()
        .iter()
        .last()
        .and_then(|op| op.ok()) // Unwrap the Result<Op>
        .filter(|op| op.kind() == Kind::SoftClip)
        .map_or((false, 0), |op| (true, op.len() as u32));
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
    check_bam_end_of_file(bam_file);
    let mut reader = bam::io::Reader::new(File::open(bam_file).unwrap());
    let header_: Header = reader.read_header().unwrap();
    let index_: Index = bai::fs::read(bam_bai_file).unwrap();
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
/// * `bam_bai_file`: BAM.BAI file. Not read: the file is read from start to end.
/// * `num_threads`: Number of threads that decompress the file.
/// * `min_mapping_quality`: Minimum mapping quality (inclusive).
///
/// # Returns
/// * HashSet of read names.
///
/// # Notes
/// * The records are read one by one and only the names are held.
pub fn get_read_names_passing_mapping_quality(
    bam_file: &str,
    bam_bai_file: &str,
    num_threads: usize,
    min_mapping_quality: u16
) -> HashSet<Box<str>> {
    let _ = bam_bai_file;
    check_bam_end_of_file(bam_file);

    // Step 1. Open the BAM file
    let file: File = File::open(bam_file).unwrap_or_else(|e| panic!("Could not open {bam_file}: {e}"));
    let workers: NonZeroUsize = NonZeroUsize::new(num_threads.max(1)).unwrap();
    let mut reader = bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(workers, file));
    reader.read_header().unwrap();

    // Step 2. Identify read names that pass the minimum mapping quality
    // If one record passes the minimum mapping quality,
    // then all records of the same read name pass the minimum mapping quality
    let mut read_names_passing_mapq: HashSet<Box<str>> = HashSet::new();
    for result in reader.records() {
        let record: bam::Record = result.unwrap_or_else(|e| panic!("Could not read a record of {bam_file}: {e}"));
        if get_alignment_mapping_quality(&record) >= min_mapping_quality {
            if let Some(read_name) = record.name().and_then(|name_bytes| std::str::from_utf8(name_bytes).ok()) {
                read_names_passing_mapq.insert(read_name.into());
            }
        }
    }

    read_names_passing_mapq
}

/// Fetches read names that has at least 1 BAM record with a splicing signal
/// based on its CIGAR string or overlaps a single-exon transcript.
///
/// # Arguments
/// * `bam_file`: BAM file.
/// * `bam_bai_file`: BAM.BAI file. Not read: the file is read from start to end.
/// * `gene_annotator`: Gene annotator.
/// * `num_threads`: Number of threads that decompress the file.
///
/// # Returns
/// * HashSet of read names that have at least 1 splicing signal or overlaps a single-exon transcript.
///
/// # Notes
/// * The records are read one by one and only the names are held.
pub fn get_read_names_with_splicing(
    bam_file: &str,
    bam_bai_file: &str,
    gene_annotator: &(impl GeneAnnotator + Sync),
    num_threads: usize
) -> HashSet<Box<str>> {
    let _ = bam_bai_file;
    check_bam_end_of_file(bam_file);

    // Step 1. Get single-exon transcripts, one interval tree per chromosome
    let mut single_exon_transcripts_map: HashMap<Box<str>, IntervalTree<()>> = HashMap::new();
    for transcript in gene_annotator.get_transcripts() {
        if transcript.get_exon_ids().len() == 1 {
            single_exon_transcripts_map
                .entry(transcript.chromosome.clone())
                .or_insert_with(IntervalTree::new)
                .insert(Interval::new(transcript.start as isize, transcript.end as isize, ()));
        }
    }

    // Step 2. Create a map of chromosome IDs and names
    let chromosomes_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);

    // Step 3. Open the BAM file
    let file: File = File::open(bam_file).unwrap_or_else(|e| panic!("Could not open {bam_file}: {e}"));
    let workers: NonZeroUsize = NonZeroUsize::new(num_threads.max(1)).unwrap();
    let mut reader = bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(workers, file));
    reader.read_header().unwrap();

    // Step 4. Identify read names that either have splicing or overlap a single-exon transcript
    let mut read_names_spliced: HashSet<Box<str>> = HashSet::new();
    for result in reader.records() {
        let record: bam::Record = result.unwrap_or_else(|e| panic!("Could not read a record of {bam_file}: {e}"));
        // Skip unmapped reads: they have no reference_sequence_id or
        // alignment_start/end and should not appear in the output.
        if record.flags().is_unmapped() {
            continue;
        }
        let Some(read_name) = record.name().and_then(|name_bytes| std::str::from_utf8(name_bytes).ok()) else {
            continue;
        };
        let is_spliced: bool = has_splicing(&record) || {
            // Check if the read overlaps a single-exon transcript
            let chromosome_id: usize = record.reference_sequence_id()
                .unwrap_or_else(|| panic!("reference_sequence_id() returned None for read: {}", read_name))
                .unwrap_or_else(|e| panic!("reference_sequence_id() returned Err for read {}: {}", read_name, e));
            let chromosome_name: &Box<str> = chromosomes_map.get_by_right(&(chromosome_id as u16))
                .unwrap_or_else(|| panic!("chromosome id {} not found in chromosomes_map for read: {}", chromosome_id, read_name));
            let start: usize = record.alignment_start()
                .unwrap_or_else(|| panic!("alignment_start() returned None for read: {}", read_name))
                .unwrap_or_else(|e| panic!("alignment_start() returned Err for read {}: {}", read_name, e))
                .get() - 1;
            let end: usize = record.alignment_end()
                .unwrap_or_else(|| panic!("alignment_end() returned None for read: {}", read_name))
                .unwrap_or_else(|e| panic!("alignment_end() returned Err for read {}: {}", read_name, e))
                .get() - 1;
            single_exon_transcripts_map
                .get(chromosome_name)
                .is_some_and(|transcripts| !transcripts.overlaps(start as isize, end as isize).is_empty())
        };
        if is_spliced {
            read_names_spliced.insert(read_name.into());
        }
    }

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


/// Get BAM depths map.
///
/// `bam_bai_file` is the index of `bam_file`, wherever it is stored.
///
/// # Returns
/// * `HashMap<chromosome name, Vec<u32>>` where `depths[i]` is the read depth at
///   1-based genomic position `i + 1`. To look up 1-based position `p`, index with
///   `p - 1`. Holds 4 bytes for every base of every contig of the header; `ReadDepths` counts at
///   chosen positions only.
pub fn get_bam_depths_map(bam_file: &str, bam_bai_file: &str, num_threads: usize) -> HashMap<Box<str>, Vec<u32>> {
    check_bam_end_of_file(bam_file);
    let index: Index = bai::fs::read(bam_bai_file)
        .unwrap_or_else(|e| panic!("Could not read the index {bam_bai_file} of {bam_file}: {e}"));
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);

    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("Failed to build Rayon thread pool");

    thread_pool.install(|| {
        chromosome_lengths
            .into_par_iter()
            .map_init(|| {
                    let mut reader = io::indexed_reader::Builder::default()
                        .set_index(index.clone())
                        .build_from_path(bam_file)
                        .unwrap();
                    let header = reader.read_header().unwrap();
                    (reader, header)
            }, |(reader, header), (chromosome, length)| {
                let start = 1usize;
                let end = length as usize;

                let region: Region = format!("{chromosome}:{start}-{end}").parse().unwrap();

                let query = reader.query(header, &region).unwrap();
                let mut depth_iter = Depth::new(header, query);

                let len = end - start + 1;
                let mut depths = vec![0u32; len];

                while let Some(result) = depth_iter.next() {
                    let (pos, d) = result.unwrap();
                    let pos = pos.get();
                    let idx = pos - start;
                    if idx < len {
                        depths[idx] = d as u32;
                    }
                }

                (chromosome, depths)
            }).collect()
    })
}


pub fn get_bam_depths(bam_file: &str, bam_bai_file: &str, num_threads: usize) -> HashSet<u32> {
    let depths_map: Arc<HashMap<Box<str>, Vec<u32>>> = Arc::new(get_bam_depths_map(bam_file, bam_bai_file, num_threads));
    depths_map
        .values()
        .flat_map(|depths| depths.iter().copied())
        .collect()
}


/// Get the greatest read depth of a BAM file, over every base of every contig, as
/// `get_bam_depths_map` counts it, without holding the depths.
pub fn get_bam_max_depth(bam_file: &str, bam_bai_file: &str, num_threads: usize) -> u32 {
    check_bam_end_of_file(bam_file);
    let index: Index = bai::fs::read(bam_bai_file)
        .unwrap_or_else(|e| panic!("Could not read the index {bam_bai_file} of {bam_file}: {e}"));
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);

    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .expect("Failed to build Rayon thread pool");

    thread_pool.install(|| {
        chromosome_lengths
            .into_par_iter()
            .map_init(|| {
                    let mut reader = io::indexed_reader::Builder::default()
                        .set_index(index.clone())
                        .build_from_path(bam_file)
                        .unwrap();
                    let header = reader.read_header().unwrap();
                    (reader, header)
            }, |(reader, header), (chromosome, length)| {
                let region: Region = format!("{chromosome}:1-{length}").parse().unwrap();
                let query = reader.query(header, &region).unwrap();
                Depth::new(header, query)
                    .map(|result| result.unwrap().1 as u32)
                    .max()
                    .unwrap_or(0)
            })
            .max()
            .unwrap_or(0)
    })
}


pub struct BAMReadDepths {
    contig_lengths: HashMap<ReferenceChromosomeName, ReferenceChromosomeLength>,

    /// HashMap<contig, HashMap<position, (depth, forward reads, reverse reads)>>
    counts: HashMap<ReferenceChromosomeName, HashMap<ReferencePosition, (ReadDepth, ReadSupport, ReadSupport)>>
}

impl BAMReadDepths {
    pub fn new(
        bam_file: &str,
        bai_file: &str,
        positions: &HashMap<ReferenceChromosomeName, Vec<ReferencePosition>>,
        max_merge_distance: u32
    ) -> Self {
        // Step 1. Open the BAM file and get contig lengths.
        let mut reader = bam::io::reader::Builder::default()
            .build_from_path(bam_file)
            .unwrap_or_else(|e| panic!("Could not open {bam_file}: {e}"));
        let header: Header = reader.read_header().unwrap_or_else(|e| panic!("Could not read the header of {bam_file}: {e}"));
        let index: Index = bai::fs::read(bai_file).unwrap_or_else(|e| panic!("Could not read {bai_file}: {e}"));
        let contig_lengths: HashMap<ReferenceChromosomeName, ReferenceChromosomeLength> = header
            .reference_sequences()
            .iter()
            .map(|(name, reference_sequence)| (name.to_string().into_boxed_str(), reference_sequence.length().get() as u32))
            .collect();

        // Step 2. Initialize BAMReadDepths.
        let mut read_depths: BAMReadDepths = Self {
            contig_lengths,
            counts: HashMap::new()
        };

        // Step 3. Get depth for each reference position.
        for (contig, positions) in positions {
            let mut positions: Vec<ReferencePosition> = positions
                .iter()
                .map(|&position| read_depths.clamp(contig, position))
                .collect();
            positions.sort_unstable();
            positions.dedup();
            for group in positions.chunk_by(|a, b| b - a < max_merge_distance) {
                let (start, end): (u32, u32) = (group[0], *group.last().unwrap());
                let region: Region = Region::new(
                    &**contig,
                    Position::new(start as usize).unwrap()..=Position::new(end as usize).unwrap()
                );
                let length: usize = (end - start + 1) as usize;
                let mut depths: Vec<u32> = vec![0; length];
                let mut strand_counts: Vec<(u32, u32)> = vec![(0, 0); length];
                for result in reader.query(&header, &index, &region).unwrap() {
                    let record: bam::Record = result.unwrap_or_else(|e| panic!("Could not read a record in {contig}:{start}-{end}: {e}"));
                    BAMReadDepths::add_depths(&record, start as usize, &mut depths);
                    BAMReadDepths::add_strand_counts(&record, start as usize, &mut strand_counts);
                }
                for &position in group {
                    let offset: usize = (position - start) as usize;
                    read_depths.insert(
                        contig,
                        position,
                        depths[offset],
                        strand_counts[offset].0,
                        strand_counts[offset].1
                    );
                }
            }
        }
        read_depths
    }

    pub fn get_depth(&self, contig: &str, position: ReferencePosition) -> ReadDepth {
        self.get(&*contig, position).0
    }

    /// (forward reads, reverse reads)
    pub fn get_strands(&self, contig: &str, position: ReferencePosition) -> (ReadSupport, ReadSupport) {
        let (_, forward, reverse) = self.get(contig, position);
        (forward, reverse)
    }

    pub fn insert(&mut self, contig: &str, position: u32, depth: u32, forward: u32, reverse: u32) {
        self.counts.entry(contig.into()).or_default().insert(position, (depth, forward, reverse));
    }

    fn clamp(&self, contig: &str, position: u32) -> u32 {
        let length: u32 = *self.contig_lengths
            .get(contig)
            .unwrap_or_else(|| panic!("ReadDepths has no contig {contig}."));
        position.clamp(1, length)
    }

    fn get(&self, contig: &str, position: u32) -> (u32, u32, u32) {
        let position: u32 = self.clamp(contig, position);
        *self.counts
            .get(contig)
            .and_then(|counts| counts.get(&position))
            .unwrap_or_else(|| panic!("ReadDepths holds no count at {contig}:{position}."))
    }

    fn add_depths<R: Record + ?Sized>(record: &R, start: usize, depths: &mut [u32]) {
        let flags = record.flags().unwrap();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_qc_fail() || flags.is_duplicate() {
            return;
        }

        let end: usize = start + depths.len() - 1;
        let mut ref_pos: usize = record.alignment_start().unwrap().unwrap().get();

        for op in record.cigar().iter() {
            let op = op.unwrap();
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    for _ in 0..op.len() {
                        if ref_pos >= start && ref_pos <= end {
                            depths[ref_pos - start] += 1;
                        }
                        ref_pos += 1;
                    }
                }
                Kind::Deletion | Kind::Skip => ref_pos += op.len(),
                Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
            }
        }
    }

    fn add_strand_counts<R: Record + ?Sized>(record: &R, start: usize, strand_counts: &mut [(u32, u32)]) {
        let flags = record.flags().unwrap();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_qc_fail() || flags.is_duplicate() {
            return;
        }

        let is_reverse: bool = flags.is_reverse_complemented();
        let end: usize = start + strand_counts.len() - 1;

        let mut ref_pos: usize = record.alignment_start().unwrap().unwrap().get();

        for op in record.cigar().iter() {
            let op = op.unwrap();
            let kind = op.kind();
            let oplen = op.len();

            match kind {
                Kind::Match |
                Kind::SequenceMatch |
                Kind::SequenceMismatch => {
                    for _ in 0..oplen {
                        if ref_pos >= start && ref_pos <= end {
                            let idx = ref_pos - start;
                            if is_reverse {
                                strand_counts[idx].1 += 1;
                            } else {
                                strand_counts[idx].0 += 1;
                            }
                        }
                        ref_pos += 1;
                    }
                }
                Kind::Deletion | Kind::Skip => {
                    ref_pos += oplen;
                }
                Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {
                    // Do nothing
                }
            }
        }
    }
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
    let tag_bytes: &[u8] = tag.as_bytes();
    if tag_bytes.len() != 2 {
        panic!("Tag must be exactly 2 characters.");
    }
    let tag_array: [u8; 2] = [tag_bytes[0], tag_bytes[1]];
    let tag: Tag = Tag::from(tag_array);
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
        None => None
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
    let tag_bytes: &[u8] = tag.as_bytes();
    if tag_bytes.len() != 2 {
        panic!("Tag must be exactly 2 characters.");
    }
    let tag_array: [u8; 2] = [tag_bytes[0], tag_bytes[1]];
    let tag: Tag = Tag::from(tag_array);
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


pub fn index_bam_records(
    bam_file: &str,
    skip_unmapped: bool,
    num_threads: usize
) -> (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) {
    check_bam_end_of_file(bam_file);
    let file: File = File::open(bam_file).unwrap();
    let workers = NonZeroUsize::new(num_threads).unwrap();
    let bgzf_reader = bgzf::io::MultithreadedReader::with_worker_count(workers, file);
    let mut reader = bam::io::Reader::from(bgzf_reader);
    let header_: Header = reader.read_header().unwrap();

    let mut record_positions_map: HashMap<usize, Vec<VirtualPosition>> = HashMap::new();
    let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    let mut read_id: usize = 1;
    // Indexed by read ID (IDs start at 1): whether the read has a record that is not supplementary.
    let mut has_primary_record: Vec<bool> = vec![false];

    let mut record: bam::Record = bam::Record::default();

    loop {
        let vp: VirtualPosition = reader.get_ref().virtual_position();
        let bytes_read: usize = reader.read_record(&mut record).unwrap();
        if bytes_read == 0 {
            // End of file
            break;
        }
        let flags: Flags = record.flags();
        if (skip_unmapped && flags.is_unmapped()) || flags.is_secondary() {
            continue;
        }
        if let Some(name) = record.name() {
            if let Ok(name_str) = std::str::from_utf8(name.as_bytes()) {
                let curr_read_name: Box<str> = name_str.into();
                let curr_read_id: usize = match read_names_map.get_by_left(&curr_read_name) {
                    Some(curr_read_id) => *curr_read_id,
                    None => {
                        let curr_read_id: usize = read_id;
                        read_names_map.insert(curr_read_name.clone(), curr_read_id);
                        has_primary_record.push(false);
                        read_id += 1;
                        curr_read_id
                    }
                };
                record_positions_map.entry(curr_read_id).or_default().push(vp);
                has_primary_record[curr_read_id] |= !flags.is_supplementary();
            } else {
                // If not valid UTF-8, skip
                continue;
            }
        }
    }

    // A read whose records in the file are all supplementary (its primary record lies outside a
    // BAM cut to a region) is left out: the primary record holds the read as sequenced, which
    // every model of a read starts from. The other reads keep their IDs.
    let num_reads: usize = read_names_map.len();
    record_positions_map.retain(|read_id, _| has_primary_record[*read_id]);
    read_names_map.retain(|_, read_id| has_primary_record[*read_id]);
    if read_names_map.len() < num_reads {
        log_info!(
            "Left out {} read(s) of {} with no primary record (supplementary records only).",
            num_reads - read_names_map.len(),
            bam_file
        );
    }

    (record_positions_map, read_names_map)
}


pub fn split_regions(
    regions: &Vec<(&str, u32, u32)>,
    chunk_size: u32
) -> Vec<(Box<str>, u32, u32)> {
    assert!(chunk_size > 0, "chunk_size must be > 0");
    let mut out: Vec<(Box<str>, u32, u32)> = Vec::new();
    for (contig, start, end) in regions.iter().copied() {
        if start == 0 || end == 0 || start > end {
            continue;
        }
        let mut curr: u32 = start;
        while curr <= end {
            let mut chunk_end: u32 = curr.saturating_add(chunk_size - 1);
            if chunk_end > end {
                chunk_end = end;
            }
            out.push((
                contig.to_string().into_boxed_str(),
                curr,
                chunk_end,
            ));
            curr = chunk_end + 1;
        }
    }

    out.sort_unstable_by(|a, b| {
        a.0.cmp(&b.0)
            .then(a.1.cmp(&b.1))
    });

    out
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
