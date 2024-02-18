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


extern crate bstr;
extern crate edit_distance;
extern crate log;
extern crate noodles;
extern crate noodles_core;
extern crate noodles_sam;
extern crate rayon;
extern crate regex;
extern crate serde;
use bstr::ByteSlice;
use edit_distance::edit_distance;
use log::info;
use noodles::bam as bam;
use noodles_core::{Region, Position};
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::value::Value;
use rayon::prelude::*;
use regex::Regex;
use serde::{Serialize, Deserialize};
use std::collections::{HashMap, HashSet};
use std::process;
use crate::clustering::cluster_variant_records;
use crate::constants::*;
use crate::defaults::*;
use crate::utilities::get_chromosomes;
use crate::variant_call::VariantCall;
use crate::variant_record::VariantRecord;


pub fn identify_rna_variants_in_read(
    record: &bam::Record,
    chromosome: &str
) -> Vec<VariantRecord> {
    let mut variant_records = Vec::new();
    let read_id: String = match std::str::from_utf8(record.name().unwrap().as_bytes()) {
        Ok(s) => s.to_string(),
        Err(e) => {
            eprintln!("Could not get the read ID: {}", e);
            std::process::exit(exitcode::DATAERR);
        },
    };
    let chromosome_1: &str = chromosome;

    // Step 1. Call SNVs, insertions, and deletions using CS tag
    const CS: [u8; 2] = [b'c', b's'];
    let cs_tag: String = match record.data().get(&CS).unwrap() {
        Ok(s) => {
            match s {
                Value::String(s) => { s.to_string() },
                _ => {
                    eprintln!("cs tag is not a string");
                    std::process::exit(exitcode::DATAERR);
                }
            }
        },
        Err(e) => {
            eprintln!("Error fetching cs tag: {}", e);
            std::process::exit(exitcode::DATAERR);
        }
    };
    let re = Regex::new(r"([:\-+*~=][0-9A-Za-z]+)").unwrap(); // or ([:][0-9]+|[-+*=][A-Za-z]+)
    let mut curr_pos: u32 = match record.alignment_start().unwrap() { // curr_pos always points
        Ok(s) => { s.get() as u32 },                                  // at the end of the last variant
        Err(e) => {
            eprintln!("Error fetching start position: {}", e);
            std::process::exit(exitcode::DATAERR);
        },
    };
    curr_pos -= 1;
    for value in re.captures_iter(&cs_tag) {
        if value[0].contains(":") {
            // Increment current position
            curr_pos += 1;

            // Remove the first character
            let mut chars = value[0].chars();
            chars.next();

            // Increment current position by the number of matched nucleotides
            let num_matched_nucleotides: u32 = chars.as_str().parse::<u32>().unwrap();
            curr_pos += num_matched_nucleotides - 1; // minus 1 is necessary here
        } else if value[0].contains("*") {
            // Increment current position
            curr_pos += 1;

            // Remove the first character
            let mut value_chars = value[0].chars();
            value_chars.next();

            // Reference and alternate alleles
            let alleles: String = value_chars.as_str().to_string();
            let reference_allele: String = alleles.chars().nth(0).unwrap().to_string().to_uppercase();
            let alternate_allele: String = alleles.chars().nth(1).unwrap().to_string().to_uppercase();

            // Record a single-nucleotide variant
            let variant_record = VariantRecord::new(
                format!("{}:{}:{}", read_id, SINGLE_NUCLEOTIDE_VARIANT.to_string(), curr_pos).to_string(),
                read_id.to_string(),
                chromosome_1.to_string(),
                curr_pos,
                chromosome_1.to_string(),
                curr_pos,
                SINGLE_NUCLEOTIDE_VARIANT.to_string(),
                reference_allele.to_string(),
                alternate_allele.to_string(),
                1
            );
            variant_records.push(variant_record);
        } else if value[0].contains("+") {
            // Remove the first character
            let mut value_chars = value[0].chars();
            value_chars.next();

            // Reference and alternate alleles
            let reference_allele: String = String::from("");
            let alternate_allele: String = value_chars.as_str().to_string().to_uppercase();
            let variant_size: u32 = alternate_allele.chars().count() as u32;

            // Record an insertion
            let variant_record = VariantRecord::new(
                format!("{}:{}:{}", read_id, INSERTION.to_string(), curr_pos).to_string(),
                read_id.to_string(),
                chromosome_1.to_string(),
                curr_pos,
                chromosome_1.to_string(),
                curr_pos,
                INSERTION.to_string(),
                reference_allele.to_string(),
                alternate_allele.to_string(),
                variant_size
            );
            variant_records.push(variant_record);
        } else if value[0].contains("-") {
            // Increment current position
            curr_pos += 1;

            // Remove the first character
            let mut value_chars = value[0].chars();
            value_chars.next();

            // Reference and alternate alleles
            let reference_allele: String = value_chars.as_str().to_string().to_uppercase();
            let alternate_allele: String = String::from("");
            let sequence: String = value_chars.as_str().to_string().to_uppercase();
            let variant_size: u32 = sequence.chars().count() as u32;

            // Record a deletion
            let variant_record = VariantRecord::new(
                format!("{}:{}:{}", read_id, DELETION.to_string(), curr_pos).to_string(),
                read_id.to_string(),
                chromosome_1.to_string(),
                curr_pos,
                chromosome_1.to_string(),
                curr_pos + variant_size - 1,
                DELETION.to_string(),
                reference_allele.to_string(),
                alternate_allele.to_string(),
                variant_size
            );
            variant_records.push(variant_record);

            // Increment current position
            curr_pos += variant_size - 1;
        } else if value[0].contains("~") {
            // Increment current position
            curr_pos += 1;

            // Remove the first character
            let mut value_chars = value[0].chars();
            value_chars.next();

            // Get splicing size
            let splicing_str: String = value_chars.as_str().to_string();
            let re_splicing = Regex::new(r"\d+").unwrap(); // r"\d+" matches one or more digits
            let splicing_size: u32 = match re_splicing.find(&splicing_str) {
                Some(matched) => {
                    let num_str = &splicing_str[matched.start()..matched.end()];
                    let size: u32 = match num_str.parse::<u32>() {
                        Ok(num) => num,
                        Err(e) => {
                            eprintln!("Failed to convert found number to u32: {}", e);
                            std::process::exit(exitcode::DATAERR);
                        },
                    };
                    size
                },
                None => {
                    eprintln!("No number found in the splicing string: {}", splicing_str);
                    std::process::exit(exitcode::DATAERR);
                }
            };

            // Record a splicing
            let variant_record = VariantRecord::new(
                format!("{}:SPL:{}", read_id, curr_pos).to_string(),
                read_id.to_string(),
                chromosome_1.to_string(),
                curr_pos,
                chromosome_1.to_string(),
                curr_pos + splicing_size - 1,
                SPLICING.to_string(),
                "".to_string(),
                "".to_string(),
                splicing_size
            );
            variant_records.push(variant_record);

            // Increment current position
            curr_pos += splicing_size - 1;
        } else {
            eprintln!("Unknown delimiter for CS tag: {}", value[0].to_string());
            std::process::exit(exitcode::DATAERR);
        }
    }

    return variant_records;
}

pub fn identify_rna_variant_records(
    bam_file: &str,
    min_mapping_quality: usize,
    num_threads: usize,
    chromosomes: Vec<&str>
) -> Result<Vec<VariantRecord>, Box<dyn std::error::Error>> {
    // Step 1. Get chromosome names and sizes
    let chromosomes_map: HashMap<usize, (String, usize)> = get_chromosomes(bam_file);

    // Step 2. Prepare output data structure
    let mut variant_records_all: Vec<VariantRecord> = Vec::new();

    // Step 3. Create a thread pool
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    // Step 4. Identify variant records
    let mut reader = bam::io::indexed_reader::Builder::default().build_from_path(bam_file)?;
    let header = reader.read_header()?;
    for (chromosome_idx, chromosome) in chromosomes_map.iter() {
        let chromosome_name: &str = &chromosome.0;
        let chromosome_length: usize = chromosome.1;
        if chromosomes.contains(&chromosome_name) == false {
            continue;
        }
        let start = Position::try_from(1)?;
        let end = Position::try_from(chromosome_length)?;
        let region = Region::new(chromosome_name, start..=end);
        let query = reader.query(&header, &region)?;
        let mut records: Vec<bam::Record> = Vec::new();
        for result in query {
            let record = result?;
            if record.flags().is_unmapped() {
                continue;
            }
            if min_mapping_quality > record.mapping_quality().unwrap().get() as usize {
                continue;
            }
            records.push(record.clone());
        }
        let variant_records_list: Vec<Vec<VariantRecord>> = thread_pool.install(|| {
            records.par_iter().map(|record| {
                let chromosome_idx: usize = match record.reference_sequence_id().unwrap() {
                    Ok(s) => s as usize,
                    Err(e) => {
                        eprintln!("Failed to fetch reference sequence ID: {}", e);
                        std::process::exit(exitcode::DATAERR);
                    },
                };
                let mut variant_records_: Vec<VariantRecord> = Vec::new();
                if let Some(&(ref chromosome_name, _)) = chromosomes_map.get(&chromosome_idx) {
                    variant_records_ = identify_rna_variants_in_read(
                        record, &chromosome_name
                    );
                }
                variant_records_
            }).collect()
        });
        for variant_records in variant_records_list {
            if variant_records.is_empty() == false {
                variant_records_all.extend(variant_records);
            }
        }
    }

    Ok(variant_records_all)
}

pub fn identify_rna_variants(
    bam_file: &str,
    sample_id: &str,
    min_reads: usize,
    min_mapping_quality: usize,
    num_threads: usize,
    min_ins_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    min_del_size_proportion: f32,
    max_bnd_distance: isize,
    clustering_grid_size: isize,
    chromosomes: Vec<&str>
) -> Vec<VariantCall> {
    let mut variant_calls_all: Vec<VariantCall> = Vec::new();
    for chromosome in chromosomes.iter() {
        let mut chromosomes_: Vec<&str> = Vec::new();
        chromosomes_.push(chromosome);

        // Step 1. Call variant records
        info!("Started calling variant records in {} reads", chromosome);
        let mut variant_records: Vec<VariantRecord> = match identify_rna_variant_records(
            &bam_file, min_mapping_quality, num_threads, chromosomes_) {
            Ok(results) => results,
            Err(e) => {
                eprintln!("Error calling variants: {}", e);
                std::process::exit(exitcode::DATAERR);
            }
        };
        info!("Finished calling variant records in {} reads", chromosome);

        // Step 2. Cluster variant records into variant calls
        info!("Started clustering variant records in {}", chromosome);
        let variant_calls: Vec<VariantCall> = cluster_variant_records(
            &mut variant_records,
            sample_id,
            RNA,
            min_reads,
            num_threads,
            min_ins_size_proportion,
            max_ins_norm_edit_distance,
            min_del_size_proportion,
            max_bnd_distance,
            clustering_grid_size
        );
        info!("Finished clustering variant records in {}", chromosome);
        variant_calls_all.extend(variant_calls);
    }
    return variant_calls_all
}
