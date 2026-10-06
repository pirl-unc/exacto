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
fn test_identify_transcript_terminus_matches_five_prime_plus_strand() {
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

    // Read 1 (chrS, plus): 50M100N50M100N50M from 251. It lacks exon 1 of ts.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("1").unwrap())
        .unwrap();
    let end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        1
    );

    assert_eq!(end_match.len(), 1, "one best match");
    let end_match: EndMatch = end_match.into_iter().next().unwrap();
    assert_eq!(end_match.terminus, TranscriptTerminus::FivePrime, "terminus");
    assert_eq!(end_match.reference_transcript_match.get_reference_gene_id(), "gs", "gene");
    assert_eq!(end_match.reference_transcript_match.get_reference_transcript_id(), "ts", "transcript");
    assert_eq!(end_match.reference_transcript_match.num_splice_junction_matches(), 2, "junction matches");
    assert_eq!(end_match.chromosome_id, *chromosome_names_map.get_by_left("chrS").unwrap(), "chromosome");
    assert_eq!(end_match.reference_position, 251, "the first aligned base is exon 2's first base");
    assert_eq!(end_match.strand, Strand::Forward, "strand");
    assert_eq!(end_match.read_position, 0, "read position");
}


#[test]
fn test_identify_transcript_terminus_matches_three_prime_plus_strand() {
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

    // Read 3 (chrS, plus): 50M100N50M100N50M from 101. It lacks exon 4 of ts.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("3").unwrap())
        .unwrap();
    let end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );

    assert_eq!(end_match.len(), 1, "one best match");
    let end_match: EndMatch = end_match.into_iter().next().unwrap();
    assert_eq!(end_match.terminus, TranscriptTerminus::ThreePrime, "terminus");
    assert_eq!(end_match.reference_transcript_match.get_reference_transcript_id(), "ts", "transcript");
    assert_eq!(end_match.reference_position, 450, "the last aligned base is exon 3's last base");
    assert_eq!(end_match.strand, Strand::Forward, "strand");
    assert_eq!(end_match.read_position, 149, "read position");
}


#[test]
fn test_identify_transcript_terminus_matches_minus_strand() {
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

    // Read 2 (chrM, minus): 50M100N50M100N50M from 301. On the minus strand the read's
    // 5' end is the HIGHEST coordinate it covers and its 3' end the lowest.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("2").unwrap())
        .unwrap();
    let five_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        1
    );
    let three_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );

    assert_eq!(five_prime_end_match.len(), 1, "one best match");
    let five_prime_end_match: EndMatch = five_prime_end_match.into_iter().next().unwrap();
    assert_eq!(five_prime_end_match.reference_transcript_match.get_reference_transcript_id(), "tm", "5' transcript");
    assert_eq!(five_prime_end_match.chromosome_id, *chromosome_names_map.get_by_left("chrM").unwrap(), "5' chromosome");
    assert_eq!(five_prime_end_match.reference_position, 650, "5' reference position");
    assert_eq!(five_prime_end_match.strand, Strand::Reverse, "5' strand");
    assert_eq!(five_prime_end_match.read_position, 0, "5' read position");

    assert_eq!(three_prime_end_match.len(), 1, "one best match");
    let three_prime_end_match: EndMatch = three_prime_end_match.into_iter().next().unwrap();
    assert_eq!(three_prime_end_match.reference_transcript_match.get_reference_transcript_id(), "tm", "3' transcript");
    assert_eq!(three_prime_end_match.reference_position, 301, "3' reference position");
    assert_eq!(three_prime_end_match.strand, Strand::Reverse, "3' strand");
    assert_eq!(three_prime_end_match.read_position, 149, "3' read position");
}


#[test]
fn test_identify_transcript_terminus_matches_skips_soft_clips() {
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

    // Read 17: 10S50M100N50M100N50M from 251. The terminal ALIGNED base is the first
    // base after the 5' soft clip.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("17").unwrap())
        .unwrap();
    let five_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        1
    );

    assert_eq!(five_prime_end_match.len(), 1, "one best match");
    let five_prime_end_match: EndMatch = five_prime_end_match.into_iter().next().unwrap();
    assert_eq!(five_prime_end_match.reference_position, 251, "5' reference position");
    assert_eq!(five_prime_end_match.read_position, 10, "5' read position is past the 10-base clip");

    // Read 12: 50M100N50M100N50M12S from 101. The terminal aligned base is the last
    // base before the 3' soft clip.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("12").unwrap())
        .unwrap();
    let three_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );

    assert_eq!(three_prime_end_match.len(), 1, "one best match");
    let three_prime_end_match: EndMatch = three_prime_end_match.into_iter().next().unwrap();
    assert_eq!(three_prime_end_match.reference_position, 450, "3' reference position");
    assert_eq!(three_prime_end_match.read_position, 149, "3' read position is before the 12-base clip");
}


#[test]
fn test_identify_transcript_terminus_matches_no_reference_transcript() {
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

    // Read 9: 60M at chrS 810, past every annotated gene.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("9").unwrap())
        .unwrap();
    let five_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        0
    );
    let three_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        0
    );

    assert!(five_prime_end_match.is_empty(), "5' end has no reference transcript");
    assert!(three_prime_end_match.is_empty(), "3' end has no reference transcript");
}


#[test]
fn test_identify_transcript_terminus_matches_terminal_base_outside_transcript() {
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

    // Read 16: 70M100N50M100N50M100N50M from 81. It matches every junction of ts, but
    // its 5' end runs 20 bases past the annotated start (101).
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("16").unwrap())
        .unwrap();
    let five_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        1
    );
    let three_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );

    assert!(five_prime_end_match.is_empty(), "ts does not span the 5' terminal base (81)");

    assert_eq!(three_prime_end_match.len(), 1, "one best match");
    let three_prime_end_match: EndMatch = three_prime_end_match.into_iter().next().unwrap();
    assert_eq!(three_prime_end_match.reference_transcript_match.get_reference_transcript_id(), "ts", "3' transcript");
    assert_eq!(three_prime_end_match.reference_position, 600, "3' reference position");
}


#[test]
fn test_identify_transcript_terminus_matches_each_end_selects_its_own_gene() {
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

    // Read 25: 25M100N50M100N50M100N150M from 126. It runs past gs's last exon, across
    // the intergenic gap and into gr exon 1, so its ends lie in different genes.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("25").unwrap())
        .unwrap();
    let five_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        0
    );
    let three_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        0
    );

    assert_eq!(five_prime_end_match.len(), 1, "one best match");
    let five_prime_end_match: EndMatch = five_prime_end_match.into_iter().next().unwrap();
    assert_eq!(five_prime_end_match.reference_transcript_match.get_reference_transcript_id(), "ts", "5' transcript");
    assert_eq!(five_prime_end_match.reference_position, 126, "5' reference position");

    assert_eq!(three_prime_end_match.len(), 1, "one best match");
    let three_prime_end_match: EndMatch = three_prime_end_match.into_iter().next().unwrap();
    assert_eq!(three_prime_end_match.reference_transcript_match.get_reference_transcript_id(), "tr", "3' transcript");
    assert_eq!(three_prime_end_match.reference_position, 700, "3' reference position");
    assert_eq!(three_prime_end_match.read_position, 274, "3' read position");
}


#[test]
fn test_identify_transcript_terminus_matches_min_num_splice_junction_matches() {
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

    // Read 25 again. It matches three junctions of ts and none of tr, which has two exons.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("25").unwrap())
        .unwrap();
    let five_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        3
    );
    let three_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );

    assert_eq!(
        five_prime_end_match[0].reference_transcript_match.num_splice_junction_matches(),
        3,
        "ts meets a minimum of 3"
    );
    assert!(three_prime_end_match.is_empty(), "tr matches no junction, so a minimum of 1 removes it");

    let five_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        4
    );
    assert!(five_prime_end_match.is_empty(), "ts cannot meet a minimum of 4");
}


#[test]
fn test_identify_transcript_terminus_matches_single_exon_reference_transcript_bypasses_minimum() {
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

    // Read 26: 60M at chrF 521, inside the single-exon gene g1 (501-600). A single-exon
    // reference transcript has no junction to match, so the minimum does not apply.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("26").unwrap())
        .unwrap();
    let five_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        5
    );
    let three_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        5
    );

    assert_eq!(five_prime_end_match.len(), 1, "one best match");
    let five_prime_end_match: EndMatch = five_prime_end_match.into_iter().next().unwrap();
    assert_eq!(five_prime_end_match.reference_transcript_match.get_reference_transcript_id(), "t1", "5' transcript");
    assert_eq!(five_prime_end_match.reference_transcript_match.num_splice_junction_matches(), 0, "5' junction matches");
    assert_eq!(five_prime_end_match.chromosome_id, *chromosome_names_map.get_by_left("chrF").unwrap(), "5' chromosome");
    assert_eq!(five_prime_end_match.reference_position, 521, "5' reference position");

    assert_eq!(three_prime_end_match.len(), 1, "one best match");
    let three_prime_end_match: EndMatch = three_prime_end_match.into_iter().next().unwrap();
    assert_eq!(three_prime_end_match.reference_transcript_match.get_reference_transcript_id(), "t1", "3' transcript");
    assert_eq!(three_prime_end_match.reference_position, 580, "3' reference position");
    assert_eq!(three_prime_end_match.read_position, 59, "3' read position");
}


#[test]
fn test_identify_transcript_terminus_matches_returns_tied_genes() {
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

    // Read 27: 60M at chrF 671, inside ga and gb, two single-exon genes on the same
    // interval. Their matches score identically, so both are returned, in rank order (here
    // transcript ID order); the tie is never broken, the end resolution stitches only what
    // they agree on.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("27").unwrap())
        .unwrap();
    assert_eq!(transcript_model.get_reference_transcript_matches().len(), 2, "one match per gene");

    let five_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        0
    );
    let three_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        0
    );

    let five_prime_transcript_ids: Vec<&str> = five_prime_end_match
        .iter()
        .map(|m| m.reference_transcript_match.get_reference_transcript_id())
        .collect();
    let three_prime_transcript_ids: Vec<&str> = three_prime_end_match
        .iter()
        .map(|m| m.reference_transcript_match.get_reference_transcript_id())
        .collect();
    assert_eq!(five_prime_transcript_ids, vec!["ta", "tb"], "5' end: both tied genes");
    assert_eq!(three_prime_transcript_ids, vec!["ta", "tb"], "3' end: both tied genes");
}


#[test]
fn test_identify_transcript_terminus_matches_fusion_read() {
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

    // Read 8: gs::gf fusion as two records, 25M100N50M30S on chrS and 75S30M on chrF.
    // Its 5' end lies in gs and its 3' end in gf, on another chromosome.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("8").unwrap())
        .unwrap();
    let five_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        0
    );
    let three_prime_end_match: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        0
    );

    assert_eq!(five_prime_end_match.len(), 1, "one best match");
    let five_prime_end_match: EndMatch = five_prime_end_match.into_iter().next().unwrap();
    assert_eq!(five_prime_end_match.reference_transcript_match.get_reference_transcript_id(), "ts", "5' transcript");
    assert_eq!(five_prime_end_match.chromosome_id, *chromosome_names_map.get_by_left("chrS").unwrap(), "5' chromosome");
    assert_eq!(five_prime_end_match.reference_position, 126, "5' reference position");
    assert_eq!(five_prime_end_match.read_position, 0, "5' read position");

    assert_eq!(three_prime_end_match.len(), 1, "one best match");
    let three_prime_end_match: EndMatch = three_prime_end_match.into_iter().next().unwrap();
    assert_eq!(three_prime_end_match.reference_transcript_match.get_reference_transcript_id(), "tf", "3' transcript");
    assert_eq!(three_prime_end_match.chromosome_id, *chromosome_names_map.get_by_left("chrF").unwrap(), "3' chromosome");
    assert_eq!(three_prime_end_match.reference_position, 280, "3' reference position");
    assert_eq!(three_prime_end_match.read_position, 104, "3' read position");
}
