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


use noodles_bam as bam;
use noodles_bam::bai;
use std::collections::HashMap;
use std::fs;

use crate::prelude::*;


pub fn check_alignment_inputs(
    bam_file: &str,
    bai_file: Option<&str>,
    reference_genome_fasta_file: &str
) {
    // Step 1. Check if the BAM file exists.
    assert_eq!(file_exists(bam_file), true, "{bam_file} does not exist.");

    // Step 2. Check if the BAI file exists.
    if let Some(bai_file) = bai_file {
        assert_eq!(file_exists(bai_file), true, "{bai_file} does not exist.");
    }

    // Step 3. Check the BAM file header.
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap_or_else(|e| panic!("{bam_file} could not be read: {e}."));
    let header = reader
        .read_header()
        .unwrap_or_else(|e| panic!("{bam_file} could not be read: {e}."));

    // Step 4. Check the BAI file.
    if let Some(bai_file) = bai_file {
        let index: bai::Index = bai::fs::read(bai_file)
            .unwrap_or_else(|e| panic!("The index {bai_file} of {bam_file} could not be read: {e}."));
    }

    // Step 5. Make sure the BAM file has cs tags.
    for result in reader.records() {
        let record: bam::Record = result.unwrap_or_else(|e| panic!("{bam_file} could not be read: {e}."));
        if record.flags().is_unmapped() {
            continue;
        }
        if !has_tag(&record, "cs") {
            panic!(
                "Read {} of {bam_file} has no cs tag; align with minimap2 --cs.",
                record.name().map(|name| name.to_string()).unwrap_or_default()
            );
        }
    }
    
    // Step 6. Check the contigs between the BAM file and the reference genome FASTA file.
    let fai_file: String = format!("{reference_genome_fasta_file}.fai");
    let fai: String = fs::read_to_string(&fai_file).unwrap_or_else(|e| panic!("{fai_file} could not be read: {e}."));

    // HashMap<chromosome, length>
    let reference_lengths: HashMap<&str, usize> = fai
        .lines()
        .filter_map(|line| {
            let mut fields = line.split('\t');
            Some((fields.next()?, fields.next()?.parse().ok()?))
        })
        .collect();

    for (i, (name, reference_sequence)) in header.reference_sequences().iter().enumerate() {
        let contig: String = name.to_string();
        let length: usize = reference_sequence.length().get();

        let reference_length: usize = *reference_lengths
            .get(contig.as_str())
            .expect("Contig {contig} not in the reference genome FASTA file.");

        if length != reference_length {
            panic!(
                "Contig {} is {} bases long in {} and {} bases long in {}.",
                contig, length, bam_file, reference_length, reference_genome_fasta_file
            );
        }
    }
}
