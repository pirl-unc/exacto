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
use crate::constants::*;
use crate::defaults::*;
use crate::utilities::get_chromosomes;
use crate::variant_call::VariantCall;
use crate::variant_record::VariantRecord;


pub fn call_rna_variants_in_read(
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
    let mut curr_pos: isize = match record.alignment_start().unwrap() { // curr_pos always points
        Ok(s) => { s.get() as isize },                                  // at the end of the last variant
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
            let num_matched_nucleotides: isize = chars.as_str().parse::<isize>().unwrap();
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
                format!("{}:SNV:{}", read_id, curr_pos).to_string(),
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
            let variant_size: isize = alternate_allele.chars().count() as isize;

            // Record an insertion
            let variant_record = VariantRecord::new(
                format!("{}:INS:{}", read_id, curr_pos).to_string(),
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
            let variant_size: isize = sequence.chars().count() as isize;

            // Record a deletion
            let variant_record = VariantRecord::new(
                format!("{}:DEL:{}", read_id, curr_pos).to_string(),
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
            let splicing_size: isize = match re_splicing.find(&splicing_str) {
                Some(matched) => {
                    let num_str = &splicing_str[matched.start()..matched.end()];
                    let size: isize = match num_str.parse::<isize>() {
                        Ok(num) => num,
                        Err(e) => {
                            eprintln!("Failed to convert found number to isize: {}", e);
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

    // Step 2. Call translocations (BND) using SA tag and CIGAR string
//     let mut sa_tag: String = "".to_string();
//     const SA: [u8; 2] = [b'S', b'A'];
//     let sa_tag: String = match record.data().get(&SA) {
//         Some(value) => {
//             let sa_tag: String = match record.data().get(&SA).unwrap() {
//                 Ok(s) => {
//                     match s {
//                         Value::String(s) => { s.to_string() },
//                         _ => { eprintln!("Value is not a string"); "".to_string() }
//                     }
//                 },
//                 Err(e) => { eprintln!("Error fetching SA tag: {}", e); "".to_string() }
//             };
//             sa_tag
//         },
//         None => { "".to_string() }
//     };
//
//     if sa_tag.is_empty() == false {
//         let mut position_1: isize = match record.alignment_start().unwrap() {
//             Ok(s) => { s.get() as isize },
//             Err(e) => { panic!("Failed to convert: {}", e); },
//         };
//
//         // Determine softclip direction
//         let mut cigar_idx: isize = 0;
//         let mut softclip_left: bool = false;
//         let mut softclip_right: bool = false;
//         for cigar in record.cigar().iter() {
//             let cigar = cigar.unwrap();
//             if cigar_idx == 0 && cigar.kind() == Kind::SoftClip {
//                 softclip_left = true;
//             }
//             if cigar_idx > 0 && cigar.kind() == Kind::SoftClip {
//                 softclip_right = true;
//             }
//             cigar_idx += 1;
//         }
//
//         // Determine softclip direction
//         if softclip_right {
//             for cigar in record.cigar().iter() {
//                 let cigar = cigar.unwrap();
//                 match cigar.kind() {
//                     Kind::Match
//                     | Kind::SequenceMatch
//                     | Kind::SequenceMismatch
//                     | Kind::Deletion
//                     | Kind::Skip => {
//                         position_1 += cigar.len() as isize;
//                     },
//                     _ => {} // Ignore other operations
//                 }
//             }
//         }
//
//         let supplementary_alignments: Vec<&str> = sa_tag.split(';').collect();
//         for supplementary_alignment in supplementary_alignments.iter() {
//             if supplementary_alignment.is_empty() {
//                 continue;
//             }
//             let parts: Vec<&str> = supplementary_alignment.split(',').collect();
//             let chromosome_2: String = parts[0].to_string();
//             let position_2: isize = parts[1].parse().unwrap();
//             let strand: String = parts[2].to_string();
//             let cigar: String = parts[3].to_string();
//             let mapq: f64 = parts[4].parse().unwrap();
//             let nm: isize = parts[5].parse().unwrap();
//
//             // Record a translocation (BND)
//             let variant_record = VariantRecord::new(
//                 read_id.to_string(),
//                 chromosome_1.to_string(),
//                 position_1,
//                 chromosome_2.to_string(),
//                 position_2,
//                 "BND".to_string(),
//                 "".to_string(),
//                 "".to_string(),
//                 -1
//             );
//             variant_records.push(variant_record);
//         }
//     }

    return variant_records;
}

pub fn call_rna_variant_records(
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
                    variant_records_ = call_rna_variants_in_read(
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

pub fn merge_variant_records(
    variant_records: &Vec<VariantRecord>,
    min_reads: usize,
    num_threads: usize,
    min_ins_size_proportion: f64,
    max_ins_norm_edit_distance: f64,
    min_del_size_proportion: f64
) -> Vec<VariantCall> {
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    // Step 1. Split variant records by chromosome
    // key      =   (chromosome_1,chromosome_2,variant_type)
    // value    =   Vec<VariantRecord>
    let mut variant_records_map: HashMap<(String,String,String), Vec<VariantRecord>> = HashMap::new();
    for variant_record in variant_records.iter() {
        variant_records_map
            .entry((variant_record.chromosome_1.to_string(),
                    variant_record.chromosome_2.to_string(),
                    variant_record.variant_type.to_string()))
            .or_insert(Vec::new())
            .push(variant_record.clone());
    }

    // Step 2. Sort Vec<VariantRecord> by VariantRecord.position_1
    thread_pool.install(|| {
        variant_records_map.par_iter_mut().for_each(|(_key, variant_records_)| {
            variant_records_.sort_by(|a, b| a.position_1.cmp(&b.position_1));
        });
    });

    // Step 3. Identify variant records that can be merged
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    let variant_calls_list: Vec<Vec<VariantCall>> = thread_pool.install(|| {
        variant_records_map.par_iter().map(|(key, variant_records_)| {
            let chromosome_1: &str = &key.0;
            let chromosome_2: &str = &key.1;
            let variant_type: &str = &key.2;

            // Identify variant calls
            let mut added_variant_record_ids: HashSet<&str> = HashSet::new();
            let mut variant_calls: Vec<VariantCall> = Vec::new();
            let mut variant_call_idx: isize = 1;
            for i in 0..variant_records_.len() {
                // Fetch variant record 1 data
                let id_1: &str = variant_records_[i].id.as_str();
                let pos_11: isize = variant_records_[i].position_1;
                let pos_12: isize = variant_records_[i].position_2;
                let variant_size_1: isize = variant_records_[i].variant_size;
                let reference_allele_1: &str = &variant_records_[i].reference_allele;
                let alternate_allele_1: &str = &variant_records_[i].alternate_allele;
                let alternate_allele_read_id_1: &str = &variant_records_[i].read_id;

                // Maximum clustering distance is a function of the variant size
                let mut max_distance: isize = (variant_size_1 as f64).log2().floor() as isize;

                if added_variant_record_ids.contains(id_1) {
                    // This read ID has been added previously
                    continue;
                }
                added_variant_record_ids.insert(id_1);

                // Data structures for VariantCall member variables
                let mut pos_1_values: Vec<f64> = Vec::new();
                let mut pos_2_values: Vec<f64> = Vec::new();
                let mut alternate_allele_read_ids: Vec<&str> = Vec::new();
                let mut variant_sizes: Vec<f64> = Vec::new();
                let mut reference_allele: &str = reference_allele_1;
                let mut alternate_allele: &str = alternate_allele_1;

                // Populate data structures with the current variant record
                pos_1_values.push(pos_11 as f64);
                pos_2_values.push(pos_12 as f64);
                alternate_allele_read_ids.push(alternate_allele_read_id_1);
                variant_sizes.push(variant_size_1 as f64);

                // Iterate through subsequent variant records (sorted)
                for j in (i + 1)..variant_records_.len() {
                    // Fetch variant record 2 data
                    let id_2: &str = variant_records_[j].id.as_str();
                    let pos_21: isize = variant_records_[j].position_1;
                    let pos_22: isize = variant_records_[j].position_2;
                    let variant_size_2: isize = variant_records_[j].variant_size;
                    let reference_allele_2: &str = &variant_records_[j].reference_allele;
                    let alternate_allele_2: &str = &variant_records_[j].alternate_allele;
                    let alternate_allele_read_id_2: &str = &variant_records_[j].read_id;

                    if added_variant_record_ids.contains(id_2) {
                        // This read ID has been added previously
                        continue;
                    }

                    // Calculate distances
                    let distance_1: isize = (pos_11 - pos_21).abs();
                    let distance_2: isize = (pos_12 - pos_22).abs();

                    // Break loop if distance_1 is greater than the
                    // maximum distance allowed
                    if distance_1 > max_distance {
                        break;
                    }

                    // Compute size proportion between the two variant records
                    let mut size_proportion: f64 = 1.0;
                    if variant_size_1 < variant_size_2 {
                        size_proportion = (variant_size_1 as f64) / (variant_size_2 as f64);
                    } else {
                        size_proportion = (variant_size_2 as f64) / (variant_size_1 as f64);
                    }

                    let mut cluster: bool = false;
                    if variant_type == SINGLE_NUCLEOTIDE_VARIANT {
                        if (distance_1 == 0) &&
                            (distance_2 == 0) &&
                            (alternate_allele_1 == alternate_allele_2) {
                            cluster = true;
                        }
                    } else if variant_type == INSERTION {
                        let edit_distance: usize = edit_distance(alternate_allele_1, alternate_allele_2);
                        let mut norm_edit_distance: f64 = 0.0;
                        if variant_size_1 < variant_size_2 {
                            norm_edit_distance = (edit_distance as f64) / (variant_size_2 as f64);
                        } else {
                            norm_edit_distance = (edit_distance as f64) / (variant_size_1 as f64);
                        }
                        if (distance_1 <= max_distance) &&
                            (distance_2 <= max_distance) &&
                            (size_proportion >= min_ins_size_proportion) &&
                            (norm_edit_distance <= max_ins_norm_edit_distance) {
                           cluster = true;
                        }
                    } else if variant_type == DELETION {
                        if (distance_1 <= max_distance) &&
                            (distance_2 <= max_distance) &&
                            (size_proportion >= min_del_size_proportion) {
                            cluster = true;
                        }
                    } else if variant_type == SPLICING {
                        if (distance_1 == 0) && (distance_2 == 0) {
                            cluster = true;
                        }
                    } else {
                        eprintln!("Unknown variant type: {}", variant_type);
                        std::process::exit(exitcode::DATAERR);
                    }

                    if cluster {
                        pos_1_values.push(pos_21 as f64);
                        pos_2_values.push(pos_22 as f64);
                        alternate_allele_read_ids.push(alternate_allele_read_id_2);
                        variant_sizes.push(variant_size_2 as f64);
                        added_variant_record_ids.insert(id_2);

                        if reference_allele_2.chars().count() > reference_allele.chars().count() {
                            reference_allele = reference_allele_2;
                        }
                        if alternate_allele_2.chars().count() > alternate_allele.chars().count() {
                            alternate_allele = alternate_allele_2;
                        }
                    }
                }

                if alternate_allele_read_ids.len() < min_reads {
                    // Not enough read support
                    continue;
                }

                let pos_1_sum: f64 = pos_1_values.iter().sum();
                let pos_1_count = pos_1_values.len() as f64;
                let pos_1_average: isize = (pos_1_sum / pos_1_count).round() as isize;
                let pos_2_sum: f64 = pos_2_values.iter().sum();
                let pos_2_count = pos_2_values.len() as f64;
                let pos_2_average: isize = (pos_2_sum / pos_2_count).round() as isize;
                let variant_size_sum: f64 = variant_sizes.iter().sum();
                let variant_size_count = variant_sizes.len() as f64;
                let variant_size_average: isize = (variant_size_sum / variant_size_count).round() as isize;

                let mut variant_id: String = "".to_string();
                if (variant_type == SINGLE_NUCLEOTIDE_VARIANT) ||
                    (variant_type == DELETION) ||
                    (variant_type == INSERTION) ||
                    (variant_type == SPLICING) {
                    variant_id = format!("{}_{}_{}", chromosome_1, variant_type, variant_call_idx).to_string();
                } else {
                    variant_id = format!("{}_{}_{}_{}", chromosome_1, chromosome_2, variant_type, variant_call_idx).to_string();
                }

                let mut variant_call = VariantCall::new(
                    variant_id.to_string(),
                    chromosome_1.to_string(),
                    pos_1_average,
                    chromosome_2.to_string(),
                    pos_2_average,
                    variant_type.to_string(),
                    reference_allele.to_string(),
                    alternate_allele.to_string(),
                    variant_size_average
                );

                // Add all read IDs
                for alternate_allele_read_id in alternate_allele_read_ids.iter() {
                    variant_call.add_alternate_allele_read_id(alternate_allele_read_id.to_string());
                }
                variant_calls.push(variant_call);
                variant_call_idx += 1;
            }
            variant_calls
        })
        .collect()
    });

    // Step 3. Collect all variant calls
    let mut variant_calls_all: Vec<VariantCall> = Vec::new();
    for variant_calls in variant_calls_list.iter() {
        variant_calls_all.extend(variant_calls.clone());
    }

    return variant_calls_all;
}

pub fn call_rna_variants(
    bam_file: &str,
    min_reads: usize,
    min_mapping_quality: usize,
    num_threads: usize,
    min_ins_size_proportion: f64,
    max_ins_norm_edit_distance: f64,
    min_del_size_proportion: f64,
    chromosomes: Vec<&str>
) -> Vec<VariantCall> {
    // Step 1. Call variant records
    info!("Started calling variant records in reads");
    let variant_records: Vec<VariantRecord> = match call_rna_variant_records(
        &bam_file, min_mapping_quality, num_threads, chromosomes) {
        Ok(results) => results,
        Err(e) => {
            eprintln!("Error calling variants: {}", e);
            std::process::exit(exitcode::DATAERR);
        }
    };
    info!("Finished calling variant records in reads");

    // Step 2. Merge variant records into variant calls
    info!("Started merging variant records in into variant calls");
    let variant_calls: Vec<VariantCall> = merge_variant_records(
        &variant_records,
        min_reads,
        num_threads,
        min_ins_size_proportion,
        max_ins_norm_edit_distance,
        min_del_size_proportion
    );
    info!("Finished merging variant records in into variant calls");

    return variant_calls
}
