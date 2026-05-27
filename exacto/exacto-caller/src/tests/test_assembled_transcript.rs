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
use noodles_bam::bai;
use noodles_bam::bai::Index;
use noodles_sam::Header;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::sync::Arc;

use crate::prelude::*;


fn build_assembled_transcript(
    bam_filename: &str,
    read_name: &str,
) -> (AssembledTranscript, BiMap<Box<str>, u16>) {
    let bam_path = Path::new("src/tests/data/bam").join(bam_filename);
    let bam_full_path = fs::canonicalize(&bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = bam_path.with_extension("bam.bai");
    let bam_bai_full_path = fs::canonicalize(&bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new("src/tests/data/fasta/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let (record_positions_map, read_names_map) = index_bam_records(bam_file, 2);

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17",
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1,
    );

    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: Alignment = Alignment::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map
            .get(&read_id)
            .unwrap()
            .iter()
            .map(|record| Arc::new(record.clone()))
            .collect(),
    );

    let alignment_structure: AlignmentStructure = alignment.get_alignment_structure().clone();

    let assembled_transcript: AssembledTranscript = AssembledTranscript::new(
        read_name,
        &alignment_structure,
        &chromosome_names_map,
        reference_genome_fasta_file,
    );

    (assembled_transcript, chromosome_names_map)
}

#[test]
fn test_assembled_transcript_basic_accessors() {
    let read_name = "m64012_507476_774164/1/ccs";
    let (assembled_transcript, _chromosome_names_map) = build_assembled_transcript(
        "rna-100-tumor_minimap2_mdtagged_sorted.bam",
        read_name,
    );

    assert_eq!(assembled_transcript.get_name(), read_name);
    assert_eq!(assembled_transcript.get_exons().len(), 11);
    assert_eq!(assembled_transcript.get_introns().len(), 10);
}

#[test]
fn test_assembled_transcript_basic_accessors_multi_alignment() {
    let read_name = "m64012_535544_475898/1/ccs";
    let (assembled_transcript, _chromosome_names_map) = build_assembled_transcript(
        "rna-103-tumor_minimap2_mdtagged_sorted.bam",
        read_name,
    );

    assert_eq!(assembled_transcript.get_name(), read_name);
    assert_eq!(assembled_transcript.get_exons().len(), 18);
    assert_eq!(assembled_transcript.get_introns().len(), 17);
}

#[test]
fn test_assembled_transcript_reference_positions_single_chromosome() {
    let (assembled_transcript, chromosome_names_map) = build_assembled_transcript(
        "rna-100-tumor_minimap2_mdtagged_sorted.bam",
        "m64012_507476_774164/1/ccs",
    );

    let chr17_id: u16 = *chromosome_names_map.get_by_left("chr17").unwrap();

    let (start_chr, start_pos) = assembled_transcript.get_reference_start_position();
    let (end_chr, end_pos) = assembled_transcript.get_reference_end_position();

    assert_eq!(start_chr, chr17_id);
    assert_eq!(end_chr, chr17_id);
    // The aggregate spans a positive-width window on the forward
    // strand; start must precede end.
    assert!(start_pos < end_pos);
}

#[test]
fn test_assembled_transcript_clone_is_equal() {
    let (assembled_transcript, _chromosome_names_map) = build_assembled_transcript(
        "rna-100-tumor_minimap2_mdtagged_sorted.bam",
        "m64012_507476_774164/1/ccs",
    );

    let cloned: AssembledTranscript = assembled_transcript.clone();

    assert_eq!(assembled_transcript, cloned);
    assert_eq!(assembled_transcript.get_name(), cloned.get_name());
    assert_eq!(
        assembled_transcript.get_exons().len(),
        cloned.get_exons().len()
    );
    assert_eq!(
        assembled_transcript.get_introns().len(),
        cloned.get_introns().len()
    );
}
