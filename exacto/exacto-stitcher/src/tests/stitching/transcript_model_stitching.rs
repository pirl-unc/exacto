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
use noodles_bgzf::VirtualPosition;
use noodles_sam as sam;
use sam::alignment::io::Write;
use std::path::Path;
use std::collections::HashMap;

use crate::modeling::rna_read_modeling::model_rna_reads;


#[test]
fn test_build_stitched_transcript_both_ends_resolved() {
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

    // Model the reads.
    let (record_positions_map, read_names_map):
        (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(&bam_file);
    let fasta_map: FastaMap = FastaMap::new(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap());
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let transcript_models: Vec<TranscriptModel> = model_rna_reads(
        &bam_file,
        &read_names_map,
        &record_positions_map,
        &chromosome_names_map,
        &fasta_map,
        &gene_annotator,
        2,
        0
    );

    // Read 4: 150 bases covering the last 25 bases of exon 1 through the first 25 of exon 4.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("4").unwrap())
        .unwrap();
    let reference_transcript_match: ReferenceTranscriptMatch = transcript_model.get_reference_transcript_matches()[0].clone();
    let five_prime_end_resolution: Option<EndResolution> = Some(EndResolution {
        terminus: TranscriptTerminus::FivePrime,
        reference_transcript_match: reference_transcript_match.clone(),
        stitch_sequence: "TTTTT".into(),
        read_join_position: 0
    });
    let three_prime_end_resolution: Option<EndResolution> = Some(EndResolution {
        terminus: TranscriptTerminus::ThreePrime,
        reference_transcript_match: reference_transcript_match.clone(),
        stitch_sequence: "AAA".into(),
        read_join_position: 149
    });

    let stitched_transcript: StitchedTranscript = build_stitched_transcript(
        transcript_model,
        &five_prime_end_resolution,
        &three_prime_end_resolution
    );

    let read_sequence: String = transcript_model.get_alignment_model().get_read_sequence();
    assert_eq!(stitched_transcript.read_id, transcript_model.get_read_id(), "read ID");
    assert_eq!(stitched_transcript.original_length, 150, "original length");
    assert_eq!(stitched_transcript.stitched_length, 158, "stitched length");
    assert_eq!(&*stitched_transcript.stitched_sequence, format!("TTTTT{}AAA", read_sequence), "stitched sequence");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(5), "5' stitch end index");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(155), "3' stitch start index");
    assert_eq!(&stitched_transcript.stitched_sequence[..5], "TTTTT", "the 5' index bounds the 5' stitch");
    assert_eq!(&stitched_transcript.stitched_sequence[155..], "AAA", "the 3' index bounds the 3' stitch");
    assert!(stitched_transcript.is_five_prime_stitched(), "5' end is stitched");
    assert!(stitched_transcript.is_three_prime_stitched(), "3' end is stitched");
}


#[test]
fn test_build_stitched_transcript_no_end_resolved() {
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

    // Model the reads.
    let (record_positions_map, read_names_map):
        (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(&bam_file);
    let fasta_map: FastaMap = FastaMap::new(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap());
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let transcript_models: Vec<TranscriptModel> = model_rna_reads(
        &bam_file,
        &read_names_map,
        &record_positions_map,
        &chromosome_names_map,
        &fasta_map,
        &gene_annotator,
        2,
        0
    );

    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("4").unwrap())
        .unwrap();
    let stitched_transcript: StitchedTranscript = build_stitched_transcript(
        transcript_model,
        &None,
        &None
    );

    let read_sequence: String = transcript_model.get_alignment_model().get_read_sequence();
    assert_eq!(&*stitched_transcript.stitched_sequence, read_sequence, "the read is kept as sequenced");
    assert_eq!(stitched_transcript.original_length, 150, "original length");
    assert_eq!(stitched_transcript.stitched_length, 150, "stitched length");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, None, "5' stitch end index");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, None, "3' stitch start index");
    assert!(!stitched_transcript.is_five_prime_stitched(), "5' end is not stitched");
    assert!(!stitched_transcript.is_three_prime_stitched(), "3' end is not stitched");
}


#[test]
fn test_build_stitched_transcript_trims_outside_the_join_positions() {
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

    // Model the reads.
    let (record_positions_map, read_names_map):
        (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(&bam_file);
    let fasta_map: FastaMap = FastaMap::new(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap());
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let transcript_models: Vec<TranscriptModel> = model_rna_reads(
        &bam_file,
        &read_names_map,
        &record_positions_map,
        &chromosome_names_map,
        &fasta_map,
        &gene_annotator,
        2,
        0
    );

    // Read 4 again. Join positions inside the read drop the read bases outside them, as
    // a 3' roll-back or a dropped polyA tail does.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("4").unwrap())
        .unwrap();
    let reference_transcript_match: ReferenceTranscriptMatch = transcript_model.get_reference_transcript_matches()[0].clone();
    let five_prime_end_resolution: Option<EndResolution> = Some(EndResolution {
        terminus: TranscriptTerminus::FivePrime,
        reference_transcript_match: reference_transcript_match.clone(),
        stitch_sequence: "".into(),
        read_join_position: 10
    });
    let three_prime_end_resolution: Option<EndResolution> = Some(EndResolution {
        terminus: TranscriptTerminus::ThreePrime,
        reference_transcript_match: reference_transcript_match.clone(),
        stitch_sequence: "".into(),
        read_join_position: 129
    });

    let stitched_transcript: StitchedTranscript = build_stitched_transcript(
        transcript_model,
        &five_prime_end_resolution,
        &three_prime_end_resolution
    );

    let read_sequence: String = transcript_model.get_alignment_model().get_read_sequence();
    assert_eq!(&*stitched_transcript.stitched_sequence, &read_sequence[10..=129], "read bases 10 through 129, both inclusive");
    assert_eq!(stitched_transcript.original_length, 150, "original length");
    assert_eq!(stitched_transcript.stitched_length, 120, "stitched length");

    // A resolved end with an empty stitch has an empty range, not `None`.
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(0), "5' stitch end index");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(120), "3' stitch start index");
}


#[test]
#[should_panic(expected = "passes the 3' join position")]
fn test_build_stitched_transcript_crossed_join_positions_panic() {
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

    // Model the reads.
    let (record_positions_map, read_names_map):
        (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(&bam_file);
    let fasta_map: FastaMap = FastaMap::new(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap());
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let transcript_models: Vec<TranscriptModel> = model_rna_reads(
        &bam_file,
        &read_names_map,
        &record_positions_map,
        &chromosome_names_map,
        &fasta_map,
        &gene_annotator,
        2,
        0
    );

    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("4").unwrap())
        .unwrap();
    let reference_transcript_match: ReferenceTranscriptMatch = transcript_model.get_reference_transcript_matches()[0].clone();
    let five_prime_end_resolution: Option<EndResolution> = Some(EndResolution {
        terminus: TranscriptTerminus::FivePrime,
        reference_transcript_match: reference_transcript_match.clone(),
        stitch_sequence: "".into(),
        read_join_position: 100
    });
    let three_prime_end_resolution: Option<EndResolution> = Some(EndResolution {
        terminus: TranscriptTerminus::ThreePrime,
        reference_transcript_match: reference_transcript_match.clone(),
        stitch_sequence: "".into(),
        read_join_position: 50
    });

    build_stitched_transcript(
        transcript_model,
        &five_prime_end_resolution,
        &three_prime_end_resolution
    );
}


#[test]
fn test_stitch_transcript_model_both_ends_degraded() {
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

    // Model the reads.
    let (record_positions_map, read_names_map):
        (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(&bam_file);
    let fasta_map: FastaMap = FastaMap::new(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap());
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let transcript_models: Vec<TranscriptModel> = model_rna_reads(
        &bam_file,
        &read_names_map,
        &record_positions_map,
        &chromosome_names_map,
        &fasta_map,
        &gene_annotator,
        2,
        0
    );

    // Read 4: both ends stop 25 bases short of ts.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("4").unwrap())
        .unwrap();
    let stitched_transcript: StitchedTranscript = stitch_transcript_model(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &fasta_map,
        1,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(100), "stitched sequence is all of ts");
    assert_eq!(stitched_transcript.original_length, 150, "original length");
    assert_eq!(stitched_transcript.stitched_length, 200, "stitched length");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(25), "5' stitch end index");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(175), "3' stitch start index");
    assert_eq!(
        stitched_transcript.five_prime_end_resolution.unwrap().reference_transcript_match.get_reference_transcript_id(),
        "ts",
        "5' transcript"
    );
    assert_eq!(
        stitched_transcript.three_prime_end_resolution.unwrap().reference_transcript_match.get_reference_transcript_id(),
        "ts",
        "3' transcript"
    );
}


#[test]
fn test_stitch_transcript_model_drops_polya_tail() {
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

    // Model the reads.
    let (record_positions_map, read_names_map):
        (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(&bam_file);
    let fasta_map: FastaMap = FastaMap::new(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap());
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let transcript_models: Vec<TranscriptModel> = model_rna_reads(
        &bam_file,
        &read_names_map,
        &record_positions_map,
        &chromosome_names_map,
        &fasta_map,
        &gene_annotator,
        2,
        0
    );

    // Read 12: 150 aligned bases and a 12-base polyA tail. The tail is dropped and exon 4
    // takes its place.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("12").unwrap())
        .unwrap();
    let stitched_transcript: StitchedTranscript = stitch_transcript_model(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &fasta_map,
        1,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    assert_eq!(stitched_transcript.original_length, 162, "original length includes the tail");
    assert_eq!(stitched_transcript.stitched_length, 200, "stitched length");
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(100), "no tail base is left between the read and exon 4");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(150), "3' stitch start index");
}


#[test]
fn test_stitch_transcript_model_keeps_clip_of_unresolved_end() {
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

    // Model the reads.
    let (record_positions_map, read_names_map):
        (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(&bam_file);
    let fasta_map: FastaMap = FastaMap::new(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap());
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let transcript_models: Vec<TranscriptModel> = model_rna_reads(
        &bam_file,
        &read_names_map,
        &record_positions_map,
        &chromosome_names_map,
        &fasta_map,
        &gene_annotator,
        2,
        0
    );

    // Read 14: 150 aligned bases and a 14-base junk clip at the 3' end. That end is left
    // alone, so the clip stays in the output as evidence.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("14").unwrap())
        .unwrap();
    let stitched_transcript: StitchedTranscript = stitch_transcript_model(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &fasta_map,
        1,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    assert_eq!(stitched_transcript.original_length, 164, "original length");
    assert_eq!(stitched_transcript.stitched_length, 164, "stitched length");
    assert_eq!(
        &*stitched_transcript.stitched_sequence,
        "GC".repeat(75) + "TTTTGTTTTGTTTT",
        "stitched sequence"
    );
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(0), "5' end is already complete");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, None, "3' end is unresolved");
}


#[test]
fn test_stitch_transcript_model_intronic_rollback() {
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

    // Model the reads.
    let (record_positions_map, read_names_map):
        (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(&bam_file);
    let fasta_map: FastaMap = FastaMap::new(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap());
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let transcript_models: Vec<TranscriptModel> = model_rna_reads(
        &bam_file,
        &read_names_map,
        &record_positions_map,
        &chromosome_names_map,
        &fasta_map,
        &gene_annotator,
        2,
        0
    );

    // Read 6: 150 exonic bases and 20 bases of intron 3. The intronic bases are dropped
    // and exon 4 is appended.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("6").unwrap())
        .unwrap();
    let stitched_transcript: StitchedTranscript = stitch_transcript_model(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &fasta_map,
        1,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    assert_eq!(stitched_transcript.original_length, 170, "original length");
    assert_eq!(stitched_transcript.stitched_length, 200, "stitched length");
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(100), "stitched sequence is all of ts");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(150), "3' stitch start index");
}


#[test]
fn test_stitch_transcript_model_fusion_read() {
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

    // Model the reads.
    let (record_positions_map, read_names_map):
        (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(&bam_file);
    let fasta_map: FastaMap = FastaMap::new(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap());
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let transcript_models: Vec<TranscriptModel> = model_rna_reads(
        &bam_file,
        &read_names_map,
        &record_positions_map,
        &chromosome_names_map,
        &fasta_map,
        &gene_annotator,
        2,
        0
    );

    // Read 8: gs::gf fusion of 105 bases. Its 5' end stops 25 bases short of ts and its
    // 3' end 20 bases short of tf. The junction between the arms is never touched.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("8").unwrap())
        .unwrap();
    let gated_stitched_transcript: StitchedTranscript = stitch_transcript_model(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &fasta_map,
        1,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    // tf matches no junction, so a minimum of 1 leaves the 3' end alone.
    assert_eq!(gated_stitched_transcript.stitched_length, 130, "stitched length with the 3' end left alone");
    assert_eq!(gated_stitched_transcript.five_prime_stitch_end_index, Some(25), "5' stitch end index");
    assert_eq!(gated_stitched_transcript.three_prime_stitch_start_index, None, "3' end is unresolved");

    let stitched_transcript: StitchedTranscript = stitch_transcript_model(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &fasta_map,
        0,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    assert_eq!(stitched_transcript.stitched_length, 150, "stitched length");
    assert_eq!(&*stitched_transcript.stitched_sequence, "GC".repeat(75), "stitched sequence");
    assert_eq!(stitched_transcript.five_prime_stitch_end_index, Some(25), "5' stitch end index");
    assert_eq!(stitched_transcript.three_prime_stitch_start_index, Some(130), "3' stitch start index");
    assert_eq!(
        stitched_transcript.five_prime_end_resolution.unwrap().reference_transcript_match.get_reference_transcript_id(),
        "ts",
        "5' transcript"
    );
    assert_eq!(
        stitched_transcript.three_prime_end_resolution.unwrap().reference_transcript_match.get_reference_transcript_id(),
        "tf",
        "3' transcript"
    );
}


#[test]
fn test_stitch_rna_reads_keeps_bam_order() {
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

    // Stitch the reads.
    let (record_positions_map, read_names_map):
        (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(&bam_file);
    let fasta_map: FastaMap = FastaMap::new(Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/reference_stitching.fa").to_str().unwrap());
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-stitcher/gene_annotations.tsv").to_str().unwrap(),
        "synthetic",
        "v1"
    );
    let stitched_transcripts: Vec<StitchedTranscript> = stitch_rna_reads(
        &bam_file,
        &read_names_map,
        &record_positions_map,
        &chromosome_names_map,
        &fasta_map,
        &gene_annotator,
        0,
        1,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6,
        2
    );

    // One stitched transcript per mapped read, in read ID order, which is BAM order.
    let mut read_ids: Vec<usize> = read_names_map.right_values().copied().collect();
    read_ids.sort();
    assert_eq!(stitched_transcripts.len(), read_ids.len(), "one stitched transcript per mapped read");
    for (stitched_transcript, read_id) in stitched_transcripts.iter().zip(read_ids.iter()) {
        assert_eq!(stitched_transcript.read_id, *read_id, "read ID");
    }
}
