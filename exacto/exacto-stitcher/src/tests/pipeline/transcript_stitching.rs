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


use super::*;

use noodles_bam as bam;
use noodles_sam as sam;
use sam::alignment::io::Write;
use std::path::Path;

use crate::io::builders::build_stitched_transcript_records;
use crate::io::records::StitchedTranscriptRecord;


#[test]
fn test_stitch_reference_transcripts_one_transcript_per_mapped_read() {
    // Convert the fixture SAM into a BAM.
    let directory = tempfile::tempdir().unwrap();
    let mut reader = sam::io::reader::Builder::default()
        .build_from_path(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.sam"))
        .unwrap();
    let header = reader.read_header().unwrap();
    let bam_path = directory.path().join("reference_stitching.bam");
    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(&bam_path)
        .unwrap();
    writer.write_header(&header).unwrap();
    for record in reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    let bam_file: String = bam_path.to_str().unwrap().to_string();

    // Stitch.
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let options: StitchReferenceTranscriptsOptions = StitchReferenceTranscriptsOptions::default();
    let stitched_transcript_set: StitchedTranscriptSet = stitch_reference_transcripts(
        &bam_file,
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap(),
        &gene_annotator,
        &options,
        2
    );

    // The pipeline identifies a read by its ID, so recover the read names the same way it does.
    let (_, read_names_map): (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);

    // 26 read names, of which read 10 is unmapped and never modeled.
    assert_eq!(stitched_transcript_set.transcripts.len(), 25, "number of stitched transcripts");
    assert!(read_names_map.get_by_left("10").is_none(), "the unmapped read has no read ID");

    // The output keeps the BAM's order of first appearance.
    let read_ids: Vec<usize> = stitched_transcript_set.transcripts.iter().map(|transcript| transcript.read_id).collect();
    let mut sorted_read_ids: Vec<usize> = read_ids.clone();
    sorted_read_ids.sort();
    assert_eq!(read_ids, sorted_read_ids, "read IDs are in ascending order");

    // The fusion read 8 has two records but is one transcript.
    let num_read_8: usize = stitched_transcript_set
        .transcripts
        .iter()
        .filter(|transcript| transcript.read_id == *read_names_map.get_by_left("8").unwrap())
        .count();
    assert_eq!(num_read_8, 1, "a split read is stitched once");
}


#[test]
fn test_stitch_reference_transcripts_degraded_ends() {
    // Convert the fixture SAM into a BAM.
    let directory = tempfile::tempdir().unwrap();
    let mut reader = sam::io::reader::Builder::default()
        .build_from_path(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.sam"))
        .unwrap();
    let header = reader.read_header().unwrap();
    let bam_path = directory.path().join("reference_stitching.bam");
    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(&bam_path)
        .unwrap();
    writer.write_header(&header).unwrap();
    for record in reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    let bam_file: String = bam_path.to_str().unwrap().to_string();

    // Stitch.
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let options: StitchReferenceTranscriptsOptions = StitchReferenceTranscriptsOptions::default();
    let stitched_transcript_set: StitchedTranscriptSet = stitch_reference_transcripts(
        &bam_file,
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap(),
        &gene_annotator,
        &options,
        2
    );

    // The pipeline identifies a read by its ID, so recover the read names the same way it does.
    let (_, read_names_map): (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);

    // Reads 1 (5'-degraded, plus), 2 (5'-degraded, minus), 3 (3'-degraded), 4 (both ends)
    // and 5 (already complete) all come out as the full 200-base reference transcript.
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("1").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(100), "read 1 stitched sequence");
    assert_eq!(stitched_transcript.original_length, 150, "read 1 original length");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(50), "read 1 5' stitch end index");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(200), "read 1 3' end is already complete");

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("2").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(100), "read 2 stitched sequence");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(50), "read 2 5' stitch end index");
    assert_eq!(
        stitched_transcript.five_prime_end_resolution.as_ref().unwrap().reference_transcript_match.get_reference_transcript_id(),
        "tm",
        "read 2 5' transcript"
    );

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("3").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(100), "read 3 stitched sequence");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(0), "read 3 5' end is already complete");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(150), "read 3 3' stitch start index");

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("4").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(100), "read 4 stitched sequence");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(25), "read 4 5' stitch end index");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(175), "read 4 3' stitch start index");

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("5").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(100), "read 5 stitched sequence");
    assert_eq!(stitched_transcript.original_length, 200, "read 5 original length");
    assert_eq!(stitched_transcript.stitched_length, 200, "read 5 stitched length");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(0), "read 5 5' end is already complete");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(200), "read 5 3' end is already complete");
}


#[test]
fn test_stitch_reference_transcripts_intronic_three_prime_ends() {
    // Convert the fixture SAM into a BAM.
    let directory = tempfile::tempdir().unwrap();
    let mut reader = sam::io::reader::Builder::default()
        .build_from_path(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.sam"))
        .unwrap();
    let header = reader.read_header().unwrap();
    let bam_path = directory.path().join("reference_stitching.bam");
    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(&bam_path)
        .unwrap();
    writer.write_header(&header).unwrap();
    for record in reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    let bam_file: String = bam_path.to_str().unwrap().to_string();

    // Stitch.
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let options: StitchReferenceTranscriptsOptions = StitchReferenceTranscriptsOptions::default();
    let stitched_transcript_set: StitchedTranscriptSet = stitch_reference_transcripts(
        &bam_file,
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap(),
        &gene_annotator,
        &options,
        2
    );

    // The pipeline identifies a read by its ID, so recover the read names the same way it does.
    let (_, read_names_map): (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);

    // Internal priming rolls back and stitches: read 6 (tract downstream, no signal),
    // read 7 (minus strand, signal and tract) and read 28 (signal, tract upstream).
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("6").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(100), "read 6 stitched sequence");
    assert_eq!(stitched_transcript.original_length, 170, "read 6 original length");

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("7").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(100), "read 7 stitched sequence");
    assert_eq!(stitched_transcript.original_length, 170, "read 7 original length");

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("28").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(100), "read 28 stitched sequence");
    assert_eq!(stitched_transcript.original_length, 164, "read 28 original length");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(100), "read 28 3' stitch start index");

    // Read 15 (signal, no tract) is a genuine intronic polyadenylation site, and read 18
    // would roll back across a cassette exon. Both keep their 3' end.
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("15").unwrap())
        .unwrap();
    assert_eq!(stitched_transcript.stitched_length, 120, "read 15 stitched length");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, None, "read 15 3' end is unresolved");
    assert!(stitched_transcript.stitched_sequence.ends_with("AATAAAGCGCGCGCGCGCGC"), "read 15 keeps its intronic bases");

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("18").unwrap())
        .unwrap();
    assert_eq!(stitched_transcript.stitched_length, 140, "read 18 stitched length");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, None, "read 18 3' end is unresolved");
}


#[test]
fn test_stitch_reference_transcripts_terminal_clips() {
    // Convert the fixture SAM into a BAM.
    let directory = tempfile::tempdir().unwrap();
    let mut reader = sam::io::reader::Builder::default()
        .build_from_path(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.sam"))
        .unwrap();
    let header = reader.read_header().unwrap();
    let bam_path = directory.path().join("reference_stitching.bam");
    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(&bam_path)
        .unwrap();
    writer.write_header(&header).unwrap();
    for record in reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    let bam_file: String = bam_path.to_str().unwrap().to_string();

    // Stitch.
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let options: StitchReferenceTranscriptsOptions = StitchReferenceTranscriptsOptions::default();
    let stitched_transcript_set: StitchedTranscriptSet = stitch_reference_transcripts(
        &bam_file,
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap(),
        &gene_annotator,
        &options,
        2
    );

    // The pipeline identifies a read by its ID, so recover the read names the same way it does.
    let (_, read_names_map): (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);

    // Read 12: a 12-base polyA tail is dropped and exon 4 appended.
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("12").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(100), "read 12 stitched sequence");
    assert_eq!(stitched_transcript.original_length, 162, "read 12 original length");

    // Reads 13 and 14: a 3' clip that is not a polyA tail blocks that end only.
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("13").unwrap())
        .unwrap();
    assert_eq!(stitched_transcript.stitched_length, 159, "read 13 stitched length");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(0), "read 13 5' end is already complete");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, None, "read 13 3' end is unresolved");

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("14").unwrap())
        .unwrap();
    assert_eq!(stitched_transcript.stitched_length, 164, "read 14 stitched length");
    assert!(stitched_transcript.stitched_sequence.ends_with("TTTTGTTTTGTTTT"), "read 14 keeps its 3' clip");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, None, "read 14 3' end is unresolved");

    // Read 17: a 5' clip that spells the last 10 bases of exon 1 is kept, and the rest of
    // exon 1 is stitched before it; the stitched transcript is ts.
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("17").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(100), "read 17 stitched sequence");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(40), "read 17 5' stitch is exon 1 up to the clip");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(200), "read 17 3' end is already complete");
}


#[test]
fn test_stitch_reference_transcripts_unresolved_reads_pass_through() {
    // Convert the fixture SAM into a BAM.
    let directory = tempfile::tempdir().unwrap();
    let mut reader = sam::io::reader::Builder::default()
        .build_from_path(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.sam"))
        .unwrap();
    let header = reader.read_header().unwrap();
    let bam_path = directory.path().join("reference_stitching.bam");
    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(&bam_path)
        .unwrap();
    writer.write_header(&header).unwrap();
    for record in reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    let bam_file: String = bam_path.to_str().unwrap().to_string();

    // Stitch.
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let options: StitchReferenceTranscriptsOptions = StitchReferenceTranscriptsOptions::default();
    let stitched_transcript_set: StitchedTranscriptSet = stitch_reference_transcripts(
        &bam_file,
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap(),
        &gene_annotator,
        &options,
        2
    );

    // The pipeline identifies a read by its ID, so recover the read names the same way it does.
    let (_, read_names_map): (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);

    // Read 9 (no reference transcript) and read 23 (backsplicing, which matches no junction)
    // come out exactly as sequenced.
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("9").unwrap())
        .unwrap();
    assert_eq!(stitched_transcript.original_length, 60, "read 9 original length");
    assert_eq!(stitched_transcript.stitched_length, 60, "read 9 stitched length");
    assert!(!stitched_transcript.is_five_prime_stitched(), "read 9 5' end");
    assert!(!stitched_transcript.is_three_prime_stitched(), "read 9 3' end");

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("23").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(50), "read 23 stitched sequence");
    assert!(!stitched_transcript.is_five_prime_stitched(), "read 23 5' end");
    assert!(!stitched_transcript.is_three_prime_stitched(), "read 23 3' end");

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("27").unwrap())
        .unwrap();
    // Read 27 lies in ga and gb, which tie; their templates agree on both ends, so both
    // are stitched (20 bases each) and reported under the first in rank order.
    assert_eq!(stitched_transcript.stitched_length, 100, "read 27 stitched length");
    assert_eq!(
        stitched_transcript.five_prime_end_resolution.as_ref().unwrap().reference_transcript_match.get_reference_transcript_id(),
        "ta",
        "read 27 5' transcript"
    );
    assert!(stitched_transcript.is_three_prime_stitched(), "read 27 3' end");

    // Read 16 overruns the annotated start and read 24 starts inside intron 1. Only their
    // 5' ends are left alone.
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("16").unwrap())
        .unwrap();
    assert_eq!(stitched_transcript.stitched_length, 220, "read 16 stitched length");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, None, "read 16 5' end is unresolved");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(220), "read 16 3' end is already complete");

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("24").unwrap())
        .unwrap();
    assert_eq!(stitched_transcript.stitched_length, 180, "read 24 stitched length");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, None, "read 24 5' end is unresolved");
}


#[test]
fn test_stitch_reference_transcripts_single_exon_gene() {
    // Convert the fixture SAM into a BAM.
    let directory = tempfile::tempdir().unwrap();
    let mut reader = sam::io::reader::Builder::default()
        .build_from_path(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.sam"))
        .unwrap();
    let header = reader.read_header().unwrap();
    let bam_path = directory.path().join("reference_stitching.bam");
    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(&bam_path)
        .unwrap();
    writer.write_header(&header).unwrap();
    for record in reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    let bam_file: String = bam_path.to_str().unwrap().to_string();

    // Stitch.
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let options: StitchReferenceTranscriptsOptions = StitchReferenceTranscriptsOptions::default();
    let stitched_transcript_set: StitchedTranscriptSet = stitch_reference_transcripts(
        &bam_file,
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap(),
        &gene_annotator,
        &options,
        2
    );

    // The pipeline identifies a read by its ID, so recover the read names the same way it does.
    let (_, read_names_map): (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);

    // Read 26: 60 unspliced bases inside the 100-base single-exon gene g1. The default
    // minimum of one junction match does not apply to a single-exon reference transcript.
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("26").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(50), "stitched sequence is all of t1");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(20), "5' stitch end index");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(80), "3' stitch start index");
}


#[test]
fn test_stitch_reference_transcripts_min_num_splice_junction_matches() {
    // Convert the fixture SAM into a BAM.
    let directory = tempfile::tempdir().unwrap();
    let mut reader = sam::io::reader::Builder::default()
        .build_from_path(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.sam"))
        .unwrap();
    let header = reader.read_header().unwrap();
    let bam_path = directory.path().join("reference_stitching.bam");
    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(&bam_path)
        .unwrap();
    writer.write_header(&header).unwrap();
    for record in reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    let bam_file: String = bam_path.to_str().unwrap().to_string();

    // Stitch.
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let options: StitchReferenceTranscriptsOptions = StitchReferenceTranscriptsOptions::default();
    let stitched_transcript_set: StitchedTranscriptSet = stitch_reference_transcripts(
        &bam_file,
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap(),
        &gene_annotator,
        &options,
        2
    );

    // The pipeline identifies a read by its ID, so recover the read names the same way it does.
    let (_, read_names_map): (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);

    // Under the default minimum of 1, an end whose gene matches no junction is left alone:
    // the 3' arm of fusion read 8 (tf), of read-through read 20 (tr) and of run-on read 25 (tr).
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("8").unwrap())
        .unwrap();
    assert_eq!(stitched_transcript.stitched_length, 130, "read 8 stitched length");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(25), "read 8 5' stitch end index");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, None, "read 8 3' end is unresolved");

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("20").unwrap())
        .unwrap();
    assert_eq!(stitched_transcript.stitched_length, 170, "read 20 stitched length");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, None, "read 20 3' end is unresolved");

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("25").unwrap())
        .unwrap();
    assert_eq!(stitched_transcript.stitched_length, 300, "read 25 stitched length");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(25), "read 25 5' stitch end index");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, None, "read 25 3' end is unresolved");
}


#[test]
fn test_stitch_reference_transcripts_without_splice_junction_minimum() {
    // Convert the fixture SAM into a BAM.
    let directory = tempfile::tempdir().unwrap();
    let mut reader = sam::io::reader::Builder::default()
        .build_from_path(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.sam"))
        .unwrap();
    let header = reader.read_header().unwrap();
    let bam_path = directory.path().join("reference_stitching.bam");
    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(&bam_path)
        .unwrap();
    writer.write_header(&header).unwrap();
    for record in reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    let bam_file: String = bam_path.to_str().unwrap().to_string();

    // Stitch.
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let mut options: StitchReferenceTranscriptsOptions = StitchReferenceTranscriptsOptions::default();
    options.min_num_splice_junction_matches = 0;
    let stitched_transcript_set: StitchedTranscriptSet = stitch_reference_transcripts(
        &bam_file,
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap(),
        &gene_annotator,
        &options,
        2
    );

    // The pipeline identifies a read by its ID, so recover the read names the same way it does.
    let (_, read_names_map): (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);

    // With no minimum, each of those 3' arms takes its own gene, where it can be resolved.
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("8").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(75), "read 8 stitched sequence");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(130), "read 8 3' stitch start index");
    assert_eq!(
        stitched_transcript.three_prime_end_resolution.as_ref().unwrap().reference_transcript_match.get_reference_transcript_id(),
        "tf",
        "read 8 3' transcript"
    );

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("20").unwrap())
        .unwrap();
    // Read 20's arm ends inside gr intron 1, which holds no A run: nothing shows the end is
    // an artifact, so it is not rolled back and the read keeps its intronic bases.
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(85), "read 20 stitched sequence");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, None, "read 20 3' end is unresolved");

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("25").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(175), "read 25 stitched sequence");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(25), "read 25 5' stitch end index");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(300), "read 25 3' stitch start index");
    assert_eq!(
        stitched_transcript.five_prime_end_resolution.as_ref().unwrap().reference_transcript_match.get_reference_transcript_id(),
        "ts",
        "read 25 5' transcript"
    );
    assert_eq!(
        stitched_transcript.three_prime_end_resolution.as_ref().unwrap().reference_transcript_match.get_reference_transcript_id(),
        "tr",
        "read 25 3' transcript"
    );
}


#[test]
fn test_stitch_reference_transcripts_min_mapping_quality() {
    // Convert the fixture SAM into a BAM.
    let directory = tempfile::tempdir().unwrap();
    let mut reader = sam::io::reader::Builder::default()
        .build_from_path(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.sam"))
        .unwrap();
    let header = reader.read_header().unwrap();
    let bam_path = directory.path().join("reference_stitching.bam");
    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(&bam_path)
        .unwrap();
    writer.write_header(&header).unwrap();
    for record in reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    let bam_file: String = bam_path.to_str().unwrap().to_string();

    // Stitch.
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let mut options: StitchReferenceTranscriptsOptions = StitchReferenceTranscriptsOptions::default();
    options.min_mapping_quality = 61;
    let stitched_transcript_set: StitchedTranscriptSet = stitch_reference_transcripts(
        &bam_file,
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap(),
        &gene_annotator,
        &options,
        2
    );

    // The pipeline identifies a read by its ID, so recover the read names the same way it does.
    let (_, read_names_map): (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);

    // Every fixture record has a mapping quality of 60. A read below the minimum is not
    // modeled; it passes through as sequenced, with the reason, so that no read is lost.
    assert_eq!(read_names_map.len(), 25, "number of mapped reads");
    assert_eq!(stitched_transcript_set.transcripts.len(), 25, "one stitched transcript per mapped read");
    for stitched_transcript in stitched_transcript_set.transcripts.iter() {
        assert_eq!(stitched_transcript.unmodeled_reason, Some("mapping_quality"), "no read meets a minimum of 61");
        assert_eq!(stitched_transcript.stitched_length, stitched_transcript.original_length, "passed through");
        assert!(!stitched_transcript.is_five_prime_stitched(), "5' end");
        assert!(!stitched_transcript.is_three_prime_stitched(), "3' end");
    }
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("1").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(75), "read 1 as sequenced");
}


#[test]
fn test_stitch_reference_transcripts_isoforms_agree_or_decline() {
    // A GC genome and one gene with three isoforms:
    //   t_full      101-150, 201-250, 301-350, 401-450
    //   t_internal  271-290, 301-350, 401-450 (an internal promoter: exon 1 in intron 2)
    //   t_alt3      101-150, 201-250, 301-350, 401-430 (an earlier cleavage site)
    let directory = tempfile::tempdir().unwrap();
    let genome: String = "GC".repeat(300);
    let fasta_file = directory.path().join("genome.fa");
    std::fs::write(&fasta_file, format!(">chrT\n{genome}\n")).unwrap();
    std::fs::write(directory.path().join("genome.fa.fai"), "chrT\t600\t6\t600\t601\n").unwrap();
    let annotation_file = directory.path().join("genes.tsv");
    std::fs::write(&annotation_file, [
        "row_type\tgene_id\ttranscript_id\texon_id\texon_number\tstrand\tchromosome\tstart\tend",
        "gene\tgt\t\t\t\t+\tchrT\t101\t450",
        "transcript\tgt\tt_full\t\t\t+\tchrT\t101\t450",
        "exon\tgt\tt_full\tt_full_e1\t1\t+\tchrT\t101\t150",
        "exon\tgt\tt_full\tt_full_e2\t2\t+\tchrT\t201\t250",
        "exon\tgt\tt_full\tt_full_e3\t3\t+\tchrT\t301\t350",
        "exon\tgt\tt_full\tt_full_e4\t4\t+\tchrT\t401\t450",
        "transcript\tgt\tt_internal\t\t\t+\tchrT\t271\t450",
        "exon\tgt\tt_internal\tt_internal_e1\t1\t+\tchrT\t271\t290",
        "exon\tgt\tt_internal\tt_internal_e2\t2\t+\tchrT\t301\t350",
        "exon\tgt\tt_internal\tt_internal_e3\t3\t+\tchrT\t401\t450",
        "transcript\tgt\tt_alt3\t\t\t+\tchrT\t101\t430",
        "exon\tgt\tt_alt3\tt_alt3_e1\t1\t+\tchrT\t101\t150",
        "exon\tgt\tt_alt3\tt_alt3_e2\t2\t+\tchrT\t201\t250",
        "exon\tgt\tt_alt3\tt_alt3_e3\t3\t+\tchrT\t301\t350",
        "exon\tgt\tt_alt3\tt_alt3_e4\t4\t+\tchrT\t401\t430"
    ].join("\n") + "\n").unwrap();

    // Read "truncated": 40M50N50M from 311, exon 3 (from 311) and exon 4.
    // Read "spliced": 20M50N50M50N20M from 231, exon 2 (from 231), exon 3, exon 4 (to 420).
    let truncated: String = format!("{}{}", &genome[310..350], &genome[400..450]);
    let spliced: String = format!("{}{}{}", &genome[230..250], &genome[300..350], &genome[400..420]);
    let sam_text: String = [
        "@HD\tVN:1.6\tSO:unsorted".to_string(),
        "@SQ\tSN:chrT\tLN:600".to_string(),
        format!("truncated\t0\tchrT\t311\t60\t40M50N50M\t*\t0\t0\t{truncated}\t{}\tcs:Z::40~gc50gc:50", "I".repeat(90)),
        format!("spliced\t0\tchrT\t231\t60\t20M50N50M50N20M\t*\t0\t0\t{spliced}\t{}\tcs:Z::20~gc50gc:50~gc50gc:20", "I".repeat(90))
    ].join("\n") + "\n";
    let sam_file = directory.path().join("reads.sam");
    std::fs::write(&sam_file, sam_text).unwrap();

    // Convert the SAM into a BAM.
    let mut reader = sam::io::reader::Builder::default()
        .build_from_path(&sam_file)
        .unwrap();
    let header = reader.read_header().unwrap();
    let bam_path = directory.path().join("reads.bam");
    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(&bam_path)
        .unwrap();
    writer.write_header(&header).unwrap();
    for record in reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    let bam_file: String = bam_path.to_str().unwrap().to_string();

    // Stitch.
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(annotation_file.to_str().unwrap(), "synthetic", "v1");
    let stitched_transcript_set: StitchedTranscriptSet = stitch_reference_transcripts(
        &bam_file,
        fasta_file.to_str().unwrap(),
        &gene_annotator,
        &StitchReferenceTranscriptsOptions::default(),
        2
    );

    // The pipeline identifies a read by its ID, so recover the read names the same way it does.
    let (_, read_names_map): (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);

    // "truncated" fits t_full and t_internal equally. Upstream of it, t_full would stitch
    // exons 1 and 2 and t_internal its own exon 1, so the 5' end is left alone. Both end at
    // 450, where the read ends, so the 3' end is complete (an empty stitch).
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("truncated").unwrap())
        .unwrap();
    assert_eq!(&*stitched_transcript.stitched_sequence, truncated.as_str(), "truncated: as sequenced");
    assert!(!stitched_transcript.is_five_prime_stitched(), "truncated: the isoforms disagree upstream");
    assert!(stitched_transcript.is_three_prime_stitched(), "truncated: the isoforms agree downstream");

    // "spliced" holds the exon 2 to exon 3 junction, which t_internal lacks. t_full and t_alt3
    // agree upstream (exon 1 and exon 2 up to 230, 80 bases) and disagree downstream (to 450
    // or to 430).
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("spliced").unwrap())
        .unwrap();
    assert_eq!(stitched_transcript.stitched_length, 170, "spliced: stitched length");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(80), "spliced: 5' stitch");
    assert_eq!(
        stitched_transcript.five_prime_end_resolution.as_ref().unwrap().reference_transcript_match.get_reference_transcript_id(),
        "t_alt3",
        "spliced: the first agreeing isoform in rank order"
    );
    assert!(!stitched_transcript.is_three_prime_stitched(), "spliced: the isoforms disagree downstream");
}


#[test]
fn test_stitch_reference_transcripts_incomplete_annotation_ends() {
    // A GC genome and one transcript, tn (101-150, 201-250, 301-350, 401-450), whose 5' end
    // GENCODE could not find (mRNA_start_NF). Read "spliced": 20M50N50M50N20M from 231.
    let directory = tempfile::tempdir().unwrap();
    let genome: String = "GC".repeat(300);
    let fasta_file = directory.path().join("genome.fa");
    std::fs::write(&fasta_file, format!(">chrT\n{genome}\n")).unwrap();
    std::fs::write(directory.path().join("genome.fa.fai"), "chrT\t600\t6\t600\t601\n").unwrap();
    let gtf_file = directory.path().join("genes.gtf");
    let gene_attributes: &str = "gene_id \"gn\"; gene_type \"protein_coding\"; gene_name \"GN\"; level 2;";
    let transcript_attributes: &str = "gene_id \"gn\"; transcript_id \"tn\"; gene_type \"protein_coding\"; gene_name \"GN\"; transcript_type \"protein_coding\"; transcript_name \"GN-201\"; level 2; tag \"basic\"; tag \"mRNA_start_NF\";";
    std::fs::write(&gtf_file, [
        format!("chrT\tTEST\tgene\t101\t450\t.\t+\t.\t{gene_attributes}"),
        format!("chrT\tTEST\ttranscript\t101\t450\t.\t+\t.\t{transcript_attributes}"),
        format!("chrT\tTEST\texon\t101\t150\t.\t+\t.\t{transcript_attributes} exon_number 1; exon_id \"tn_e1\";"),
        format!("chrT\tTEST\texon\t201\t250\t.\t+\t.\t{transcript_attributes} exon_number 2; exon_id \"tn_e2\";"),
        format!("chrT\tTEST\texon\t301\t350\t.\t+\t.\t{transcript_attributes} exon_number 3; exon_id \"tn_e3\";"),
        format!("chrT\tTEST\texon\t401\t450\t.\t+\t.\t{transcript_attributes} exon_number 4; exon_id \"tn_e4\";")
    ].join("\n") + "\n").unwrap();
    let spliced: String = format!("{}{}{}", &genome[230..250], &genome[300..350], &genome[400..420]);
    let sam_text: String = [
        "@HD\tVN:1.6\tSO:unsorted".to_string(),
        "@SQ\tSN:chrT\tLN:600".to_string(),
        format!("spliced\t0\tchrT\t231\t60\t20M50N50M50N20M\t*\t0\t0\t{spliced}\t{}\tcs:Z::20~gc50gc:50~gc50gc:20", "I".repeat(90))
    ].join("\n") + "\n";
    let sam_file = directory.path().join("reads.sam");
    std::fs::write(&sam_file, sam_text).unwrap();

    // Convert the SAM into a BAM.
    let mut reader = sam::io::reader::Builder::default()
        .build_from_path(&sam_file)
        .unwrap();
    let header = reader.read_header().unwrap();
    let bam_path = directory.path().join("reads.bam");
    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(&bam_path)
        .unwrap();
    writer.write_header(&header).unwrap();
    for record in reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    let bam_file: String = bam_path.to_str().unwrap().to_string();

    // Stitch.
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_file.to_str().unwrap(), "synthetic", "v1");
    let stitched_transcript_set: StitchedTranscriptSet = stitch_reference_transcripts(
        &bam_file,
        fasta_file.to_str().unwrap(),
        &gene_annotator,
        &StitchReferenceTranscriptsOptions::default(),
        2
    );

    // tn's annotated start is where the evidence stopped, so the 5' end is left alone; its
    // 3' end is complete, so the 3' end takes exon 4 from 421 (30 bases).
    assert_eq!(stitched_transcript_set.transcripts.len(), 1, "one read");
    let stitched_transcript: &StitchedTranscript = &stitched_transcript_set.transcripts[0];
    assert!(!stitched_transcript.is_five_prime_stitched(), "mRNA_start_NF: no 5' stitch");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(90), "3' stitch start index");
    assert_eq!(stitched_transcript.stitched_length, 120, "stitched length");
}


#[test]
fn test_stitch_reference_transcripts_every_mapped_read_gives_a_row() {
    // A GC genome and one transcript, tt (101-150, 201-250, 301-350, 401-450). Three copies
    // of one read (20M50N50M50N20M from 231): "mapq_255" has no mapping quality, "no_cs_tag"
    // cannot be modeled, "mapq_5" is under a minimum of 20.
    let directory = tempfile::tempdir().unwrap();
    let genome: String = "GC".repeat(300);
    let fasta_file = directory.path().join("genome.fa");
    std::fs::write(&fasta_file, format!(">chrT\n{genome}\n")).unwrap();
    std::fs::write(directory.path().join("genome.fa.fai"), "chrT\t600\t6\t600\t601\n").unwrap();
    let annotation_file = directory.path().join("genes.tsv");
    std::fs::write(&annotation_file, [
        "row_type\tgene_id\ttranscript_id\texon_id\texon_number\tstrand\tchromosome\tstart\tend",
        "gene\tgt\t\t\t\t+\tchrT\t101\t450",
        "transcript\tgt\ttt\t\t\t+\tchrT\t101\t450",
        "exon\tgt\ttt\ttt_e1\t1\t+\tchrT\t101\t150",
        "exon\tgt\ttt\ttt_e2\t2\t+\tchrT\t201\t250",
        "exon\tgt\ttt\ttt_e3\t3\t+\tchrT\t301\t350",
        "exon\tgt\ttt\ttt_e4\t4\t+\tchrT\t401\t450"
    ].join("\n") + "\n").unwrap();
    let spliced: String = format!("{}{}{}", &genome[230..250], &genome[300..350], &genome[400..420]);
    let qualities: String = "I".repeat(90);
    let sam_text: String = [
        "@HD\tVN:1.6\tSO:unsorted".to_string(),
        "@SQ\tSN:chrT\tLN:600".to_string(),
        format!("mapq_255\t0\tchrT\t231\t255\t20M50N50M50N20M\t*\t0\t0\t{spliced}\t{qualities}\tcs:Z::20~gc50gc:50~gc50gc:20"),
        format!("no_cs_tag\t0\tchrT\t231\t60\t20M50N50M50N20M\t*\t0\t0\t{spliced}\t{qualities}"),
        format!("mapq_5\t0\tchrT\t231\t5\t20M50N50M50N20M\t*\t0\t0\t{spliced}\t{qualities}\tcs:Z::20~gc50gc:50~gc50gc:20")
    ].join("\n") + "\n";
    let sam_file = directory.path().join("reads.sam");
    std::fs::write(&sam_file, sam_text).unwrap();

    // Convert the SAM into a BAM.
    let mut reader = sam::io::reader::Builder::default()
        .build_from_path(&sam_file)
        .unwrap();
    let header = reader.read_header().unwrap();
    let bam_path = directory.path().join("reads.bam");
    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(&bam_path)
        .unwrap();
    writer.write_header(&header).unwrap();
    for record in reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    let bam_file: String = bam_path.to_str().unwrap().to_string();

    // Stitch.
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(annotation_file.to_str().unwrap(), "synthetic", "v1");
    let mut options: StitchReferenceTranscriptsOptions = StitchReferenceTranscriptsOptions::default();
    options.min_mapping_quality = 20;
    let stitched_transcript_set: StitchedTranscriptSet = stitch_reference_transcripts(
        &bam_file,
        fasta_file.to_str().unwrap(),
        &gene_annotator,
        &options,
        2
    );

    // The pipeline identifies a read by its ID, so recover the read names the same way it does.
    let (_, read_names_map): (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);
    assert_eq!(stitched_transcript_set.transcripts.len(), 3, "one stitched transcript per mapped read");

    // A missing mapping quality says nothing against the read: it is modeled and stitched
    // (exon 1 and exon 2 up to 230; exon 4 from 421).
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("mapq_255").unwrap())
        .unwrap();
    assert_eq!(stitched_transcript.unmodeled_reason, None, "mapq_255 is modeled");
    assert_eq!(stitched_transcript.stitched_length, 200, "mapq_255 is stitched");

    // The other two pass through as sequenced, with the reason.
    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("no_cs_tag").unwrap())
        .unwrap();
    assert_eq!(stitched_transcript.unmodeled_reason, Some("modeling_error"), "no_cs_tag reason");
    assert_eq!(&*stitched_transcript.stitched_sequence, spliced.as_str(), "no_cs_tag as sequenced");
    assert!(!stitched_transcript.is_five_prime_stitched() && !stitched_transcript.is_three_prime_stitched(), "no_cs_tag ends");

    let stitched_transcript: &StitchedTranscript = stitched_transcript_set
        .transcripts
        .iter()
        .find(|transcript| transcript.read_id == *read_names_map.get_by_left("mapq_5").unwrap())
        .unwrap();
    assert_eq!(stitched_transcript.unmodeled_reason, Some("mapping_quality"), "mapq_5 reason");
    assert_eq!(&*stitched_transcript.stitched_sequence, spliced.as_str(), "mapq_5 as sequenced");
}


#[test]
#[should_panic(expected = "has no cs tag on its first mapped record")]
fn test_stitch_reference_transcripts_bam_without_cs_tags_panics() {
    // A GC genome, one transcript, and one read aligned without a cs tag.
    let directory = tempfile::tempdir().unwrap();
    let genome: String = "GC".repeat(300);
    let fasta_file = directory.path().join("genome.fa");
    std::fs::write(&fasta_file, format!(">chrT\n{genome}\n")).unwrap();
    std::fs::write(directory.path().join("genome.fa.fai"), "chrT\t600\t6\t600\t601\n").unwrap();
    let annotation_file = directory.path().join("genes.tsv");
    std::fs::write(&annotation_file, [
        "row_type\tgene_id\ttranscript_id\texon_id\texon_number\tstrand\tchromosome\tstart\tend",
        "gene\tgt\t\t\t\t+\tchrT\t101\t450",
        "transcript\tgt\ttt\t\t\t+\tchrT\t101\t450",
        "exon\tgt\ttt\ttt_e1\t1\t+\tchrT\t101\t450"
    ].join("\n") + "\n").unwrap();
    let sam_text: String = [
        "@HD\tVN:1.6\tSO:unsorted".to_string(),
        "@SQ\tSN:chrT\tLN:600".to_string(),
        format!("no_cs_tag\t0\tchrT\t201\t60\t50M\t*\t0\t0\t{}\t{}", &genome[200..250], "I".repeat(50))
    ].join("\n") + "\n";
    let sam_file = directory.path().join("reads.sam");
    std::fs::write(&sam_file, sam_text).unwrap();

    // Convert the SAM into a BAM.
    let mut reader = sam::io::reader::Builder::default()
        .build_from_path(&sam_file)
        .unwrap();
    let header = reader.read_header().unwrap();
    let bam_path = directory.path().join("reads.bam");
    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(&bam_path)
        .unwrap();
    writer.write_header(&header).unwrap();
    for record in reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    let bam_file: String = bam_path.to_str().unwrap().to_string();

    // Stitch: a BAM aligned without cs tags fails the call instead of passing every read
    // through unstitched.
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(annotation_file.to_str().unwrap(), "synthetic", "v1");
    stitch_reference_transcripts(
        &bam_file,
        fasta_file.to_str().unwrap(),
        &gene_annotator,
        &StitchReferenceTranscriptsOptions::default(),
        2
    );
}


#[test]
fn test_build_stitched_transcript_records_retained_read_bases() {
    // Convert the fixture SAM into a BAM.
    let directory = tempfile::tempdir().unwrap();
    let mut reader = sam::io::reader::Builder::default()
        .build_from_path(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.sam"))
        .unwrap();
    let header = reader.read_header().unwrap();
    let bam_path = directory.path().join("reference_stitching.bam");
    let mut writer = bam::io::writer::Builder::default()
        .build_from_path(&bam_path)
        .unwrap();
    writer.write_header(&header).unwrap();
    for record in reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &record.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    let bam_file: String = bam_path.to_str().unwrap().to_string();

    // Stitch and build the records.
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let stitched_transcript_set: StitchedTranscriptSet = stitch_reference_transcripts(
        &bam_file,
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap(),
        &gene_annotator,
        &StitchReferenceTranscriptsOptions::default(),
        2
    );
    let (_, read_names_map): (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);
    let records: Vec<StitchedTranscriptRecord> = build_stitched_transcript_records(&stitched_transcript_set, &read_names_map).collect();
    assert_eq!(records.len(), 25, "one record per mapped read");

    // Read 6 ends 20 bases into intron 3, next to an A tract: it rolls back to exon 3's last
    // base, so its last 20 read bases are not in the stitched sequence, and the record says so.
    let record: &StitchedTranscriptRecord = records.iter().find(|record| &*record.read_name == "6").unwrap();
    assert_eq!(record.read_length, 170, "read 6 read length");
    assert_eq!((record.retained_read_start, record.retained_read_end), (0, 150), "read 6 keeps read bases [0, 150)");
    assert_eq!(record.three_prime_stitch_start - record.five_prime_stitch_end, 150, "read 6 retained bases in the stitched sequence");
    assert_eq!(&*record.unmodeled_reason, "", "read 6 is modeled");

    // Read 1 keeps every read base.
    let record: &StitchedTranscriptRecord = records.iter().find(|record| &*record.read_name == "1").unwrap();
    assert_eq!((record.retained_read_start, record.retained_read_end), (0, 150), "read 1 keeps read bases [0, 150)");
    assert_eq!(record.read_length, 150, "read 1 read length");
}
