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
use noodles_bam::bai as bai;
use noodles_bgzf as bgzf;
use noodles_sam::alignment::Record as _;
use noodles_sam::Header;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::BufWriter;
use std::num::NonZeroUsize;


pub fn remove_unspliced_rnas(
    bam_file: &str,
    bam_bai_file: &str,
    gene_annotator: &(impl GeneAnnotator + Sync),
    output_bam_file: &str,
    output_bam_bai_file: &str,
    num_threads: usize,
    min_mapping_quality: u16,
    write_output_bam_file: bool
) -> HashSet<Box<str>> {
    let _ = bam_bai_file;
    if write_output_bam_file {
        for output_file in [output_bam_file, output_bam_bai_file] {
            for input_file in [bam_file, bam_bai_file] {
                let is_same_file: bool = matches!(
                    (fs::canonicalize(output_file), fs::canonicalize(input_file)),
                    (Ok(output_path), Ok(input_path)) if output_path == input_path
                );
                assert!(!is_same_file, "Output file {output_file} is the input file {input_file}; write the output to a new file.");
            }
        }
    }
    check_bam_end_of_file(bam_file);
    let open_reader = || -> bam::io::Reader<bgzf::io::MultithreadedReader<File>> {
        let file: File = File::open(bam_file).unwrap_or_else(|error| panic!("Could not open {bam_file}: {error}"));
        let workers: NonZeroUsize = NonZeroUsize::new(num_threads.max(1)).unwrap();
        bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(workers, file))
    };

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
    let chromosomes_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);

    // Step 2. Decide every read in one pass: a read is kept when one of its primary or
    // supplementary records reaches the minimum mapping quality and splices or overlaps a
    // single-exon transcript.
    let mut read_names_to_keep: HashSet<Box<str>> = HashSet::new();
    let mut reader = open_reader();
    reader.read_header().unwrap();
    for result in reader.records() {
        let record: bam::Record = result.unwrap_or_else(|error| panic!("Could not read a record of {bam_file}: {error}"));
        if record.flags().is_unmapped() || record.flags().is_secondary()
            || get_alignment_mapping_quality(&record) < min_mapping_quality {
            continue;
        }
        let Some(read_name) = record.name().and_then(|name_bytes| std::str::from_utf8(name_bytes).ok()) else {
            continue;
        };
        if read_names_to_keep.contains(read_name) {
            continue;
        }
        let is_spliced: bool = has_splicing(&record) || {
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
            single_exon_transcripts_map.get(chromosome_name).is_some_and(|transcripts| {
                transcripts.iter().any(|transcript| {
                    overlaps(transcript.start as isize, transcript.end as isize, start as isize, end as isize)
                })
            })
        };
        if is_spliced {
            read_names_to_keep.insert(read_name.into());
        }
    }

    // Step 3. Write the primary and supplementary records of the kept reads, in the order of the
    // input, which is coordinate-sorted, as a BAI index requires. Unmapped records are skipped
    // even for a kept read: `bam::fs::index` rejects a record with no reference sequence.
    if write_output_bam_file {
        let mut reader = open_reader();
        let header: Header = reader.read_header().unwrap();
        let file: File = File::create(output_bam_file).unwrap_or_else(|error| panic!("Could not create {output_bam_file}: {error}"));
        let mut writer = bam::io::Writer::new(BufWriter::new(file));
        writer.write_header(&header).unwrap();
        for result in reader.records() {
            let record: bam::Record = result.unwrap_or_else(|error| panic!("Could not read a record of {bam_file}: {error}"));
            let is_kept: bool = !record.flags().is_unmapped()
                && !record.flags().is_secondary()
                && record
                    .name()
                    .and_then(|name_bytes| std::str::from_utf8(name_bytes).ok())
                    .is_some_and(|read_name| read_names_to_keep.contains(read_name));
            if is_kept {
                writer.write_record(&header, &record).unwrap();
            }
        }
        writer.try_finish().unwrap();
        drop(writer);
        let index = bam::fs::index(output_bam_file).unwrap();
        bai::fs::write(output_bam_bai_file, &index).unwrap();
    }

    read_names_to_keep
}


#[cfg(test)]
#[path = "../tests/pipeline/unspliced_rna_filtering.rs"]
mod tests;