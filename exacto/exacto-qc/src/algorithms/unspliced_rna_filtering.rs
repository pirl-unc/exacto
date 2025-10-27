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
use exacto_core::prelude::*;
use noodles_bam as bam;
use noodles_sam::Header;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};


/// Remove unspliced RNAs.
///
/// # Arguments
/// * `bam_file`: BAM file of assembled transcripts aligned to a reference genome.
/// * `bam_bai_file`: BAM.BAI file of assembled transcripts aligned to a reference genome.
/// * `fasta_file`: FASTA file of assembled transcripts.
/// * `gene_annotator`: Gene annotator.
/// * `output_bam_file`: Output BAM file.
/// * `output_bam_bai_file`: Output BAM.BAI file.
/// * `output_fasta_file`: Output FASTA file.
/// * `num_threads`: Number of threads.
/// * `min_mapping_quality`: Minimum mapping quality (inclusive).
///
/// # Notes
/// * Unspliced RNAs are removed based on the BAM file.
/// * Transcripts that do not meet the minimum mapping quality are removed.
pub fn remove_unspliced_rnas(
    bam_file: &str,
    bam_bai_file: &str,
    fasta_file: &str,
    gene_annotator: &(impl GeneAnnotator + Sync),
    output_bam_file: &str,
    output_bam_bai_file: &str,
    output_fasta_file: &str,
    num_threads: usize,
    min_mapping_quality: usize
) {
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

    // Step 2. Get read names passing the minimum mapping quality
    let read_names_passing_mapq: HashSet<Box<str>> = get_read_names_passing_mapping_quality(
        bam_file,
        bam_bai_file,
        num_threads,
        min_mapping_quality
    );

    // Step 3. Get read names with splicing
    let read_names_spliced: HashSet<Box<str>> = get_read_names_with_splicing(
        bam_file,
        bam_bai_file,
        gene_annotator,
        num_threads
    );

    // Step 4. Identify read names to keeps
    let read_names_to_keep: HashSet<Box<str>> = read_names_passing_mapq.intersection(&read_names_spliced).cloned().collect();

    // Step 5. Read all records into memory
    let read_names_map: BiMap<Box<str>, usize> = create_read_names_map(bam_file, bam_bai_file, 1);
    let records: HashMap<usize, Vec<bam::Record>> = fetch_all_bam_records(
        bam_file,
        bam_bai_file,
        &read_names_map,
        num_threads
    );

    // Step 6. Identify BAM records to keep
    let mut records_to_keep: Vec<&bam::Record> = Vec::new();
    for (read_id, records) in records.iter() {
        let read_name: Box<str> = read_names_map.get_by_right(&read_id).unwrap().to_owned();
        if read_names_to_keep.contains(&read_name) {
            records_to_keep.extend(records);
        }
    }

    // Step 7. Write records to output BAM file
    let header: Header = get_bam_header(bam_file);
    write_bam_file(
        output_bam_file,
        output_bam_bai_file,
        &header,
        &records_to_keep
    );

    // Step 8. Write to FASTA file
    let sequence_ids: Vec<(Box<str>, usize)> = get_fasta_sequence_ids(fasta_file);
    let mut sequences: Vec<(Box<str>, Box<str>)> = Vec::new();
    for (sequence_id, sequence_length) in sequence_ids.iter() {
        if read_names_to_keep.contains(&**sequence_id) {
            let sequence: Box<str> = get_fasta_sequence(&**sequence_id, 1, *sequence_length, fasta_file);
            sequences.push((sequence_id.clone(), sequence));
        }
    }
    write_fasta_file(&sequences, output_fasta_file);
}
