extern crate bam;

use std::str;
use lazy_static::lazy_static;
use regex::Regex;
use polars::prelude::*;
use polars::df;


/// Identifies DNA variants from a read.
///
/// # Arguments
/// * `read_id`     - Read ID.
/// * `start_pos`   - Start position of read.
/// * `cs_tag`      - CS tag.
///
/// # Returns
/// * `PolarsResult<DataFrame>`
pub fn identify_dna_variants_in_cs_tag(read_id: &str,
                                       start_pos: i32,
                                       cs_tag: &str) -> PolarsResult<DataFrame> {
    // Step 1. Prepare vectors
    let mut vec_read_ids: Vec<String> = Vec::new();
    let mut vec_positions: Vec<i32> = Vec::new();
    let mut vec_variant_types: Vec<String> = Vec::new();
    let mut vec_reference_alleles: Vec<String> = Vec::new();
    let mut vec_alternate_alleles: Vec<String> = Vec::new();
    let mut vec_sequences: Vec<String> = Vec::new();
    let mut vec_variant_sizes: Vec<i32> = Vec::new();

    // Step 2. Parse the cs_tag
    let re = Regex::new(r"([:\-+*=][0-9A-Za-z]+)").unwrap(); // or ([:][0-9]+|[-+*=][A-Za-z]+)
    let mut curr_pos: i32 = start_pos;
    for value in re.captures_iter(cs_tag) {
        if value[0].contains(":") {
            // Remove the first character
            let mut chars = value[0].chars();
            chars.next();

            // Increment current position by the number of matched nucleotides
            let num_matched_nucleotides: i32 = chars.as_str().parse::<i32>().unwrap();
            curr_pos += num_matched_nucleotides;
        } else if value[0].contains("*") {
            // Remove the first character
            let mut value_chars = value[0].chars();
            value_chars.next();

            // Reference and alternate alleles
            let alleles: String = value_chars.as_str().to_string();
            let reference_allele: String = alleles.chars().nth(0).unwrap().to_string();
            let alternate_allele: String = alleles.chars().nth(1).unwrap().to_string();
            let sequence: String = alleles.chars().nth(1).unwrap().to_string();

            // Record a single-nucleotide variant
            vec_read_ids.push(read_id.to_string());
            vec_positions.push(curr_pos);
            vec_variant_types.push(String::from("snv"));
            vec_reference_alleles.push(reference_allele);
            vec_alternate_alleles.push(alternate_allele);
            vec_sequences.push(sequence);
            vec_variant_sizes.push(1);

            // Increment current position
            curr_pos += 1;
        } else if value[0].contains("+") {
            // Remove the first character
            let mut value_chars = value[0].chars();
            value_chars.next();

            // Reference and alternate alleles
            let reference_allele: String = String::from("");
            let alternate_allele: String = value_chars.as_str().to_string();
            let sequence: String = value_chars.as_str().to_string();
            let variant_size: i32 = sequence.chars().count() as i32;

            // Record a single-nucleotide variant
            vec_read_ids.push(read_id.to_string());
            vec_positions.push(curr_pos);
            vec_variant_types.push(String::from("insertion"));
            vec_reference_alleles.push(reference_allele);
            vec_alternate_alleles.push(alternate_allele);
            vec_sequences.push(sequence);
            vec_variant_sizes.push(variant_size);
        } else if value[0].contains("-") {
            // Remove the first character
            let mut value_chars = value[0].chars();
            value_chars.next();

            // Reference and alternate alleles
            let reference_allele: String = value_chars.as_str().to_string();
            let alternate_allele: String = String::from("");
            let sequence: String = value_chars.as_str().to_string();
            let variant_size: i32 = sequence.chars().count() as i32;

            // Record a single-nucleotide variant
            vec_read_ids.push(read_id.to_string());
            vec_positions.push(curr_pos + 1);
            vec_variant_types.push(String::from("deletion"));
            vec_reference_alleles.push(reference_allele);
            vec_alternate_alleles.push(alternate_allele);
            vec_sequences.push(sequence);
            vec_variant_sizes.push(variant_size);

            // Increment current position
            curr_pos += variant_size;
        } else {
            println!("Unknown delimiter for CS tag: {}", value[0].to_string());
            panic!("Closing program...");
        }
    }

    // Step 3. Construct a DataFrame
    let mut df: PolarsResult<DataFrame> = df!(
        "position" => &vec_positions,
        "read_id" => &vec_read_ids,
        "variant_type" => &vec_variant_types,
        "reference_allele" => &vec_reference_alleles,
        "alternate_allele" => &vec_alternate_alleles,
        "sequence" => &vec_sequences,
        "variant_size" => &vec_variant_sizes
    );
    return df;
}