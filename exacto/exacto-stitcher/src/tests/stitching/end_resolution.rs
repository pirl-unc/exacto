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
use crate::stitching::end_matching::identify_transcript_terminus_matches;


#[test]
fn test_resolve_five_prime_end_exonic_plus_strand() {
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

    // Read 1 (chrS, plus): 50M100N50M100N50M from 251. Its first aligned base is exon 2's
    // first base, the 51st base of ts, so all of exon 1 is prepended.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("1").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_five_prime_end(transcript_model, &end_matches, &gene_annotator, &chromosome_names_map, &fasta_map);

    let end_resolution: EndResolution = end_resolution.unwrap();
    assert_eq!(end_resolution.terminus, TranscriptTerminus::FivePrime, "terminus");
    assert_eq!(end_resolution.reference_transcript_match.get_reference_transcript_id(), "ts", "transcript");
    assert_eq!(&*end_resolution.stitch_sequence, "GC".repeat(25), "stitch sequence is exon 1");
    assert_eq!(end_resolution.read_join_position, 0, "read join position");
}


#[test]
fn test_resolve_five_prime_end_exonic_minus_strand() {
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

    // Read 2 (chrM, minus): mirror of read 1. Exon 1 of tm (751-800) is prepended, reverse
    // complemented into transcript orientation.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("2").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_five_prime_end(transcript_model, &end_matches, &gene_annotator, &chromosome_names_map, &fasta_map);

    let end_resolution: EndResolution = end_resolution.unwrap();
    assert_eq!(end_resolution.reference_transcript_match.get_reference_transcript_id(), "tm", "transcript");
    assert_eq!(&*end_resolution.stitch_sequence, "GC".repeat(25), "stitch sequence is exon 1");
    assert_eq!(end_resolution.read_join_position, 0, "read join position");
}


#[test]
fn test_resolve_five_prime_end_partial_exon() {
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

    // Read 4: 25M100N50M100N50M100N25M from 126. Its first aligned base is the 26th base
    // of exon 1, so only the first 25 bases are prepended.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("4").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_five_prime_end(transcript_model, &end_matches, &gene_annotator, &chromosome_names_map, &fasta_map);

    let end_resolution: EndResolution = end_resolution.unwrap();
    assert_eq!(&*end_resolution.stitch_sequence, "GC".repeat(12) + "G", "stitch sequence is the first 25 bases of exon 1");
    assert_eq!(end_resolution.read_join_position, 0, "read join position");
}


#[test]
fn test_resolve_five_prime_end_already_complete() {
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

    // Read 5: all four exons of ts. The end resolves with nothing to add.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("5").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_five_prime_end(transcript_model, &end_matches, &gene_annotator, &chromosome_names_map, &fasta_map);

    let end_resolution: EndResolution = end_resolution.unwrap();
    assert_eq!(&*end_resolution.stitch_sequence, "", "stitch sequence");
    assert_eq!(end_resolution.read_join_position, 0, "read join position");
}


#[test]
fn test_resolve_five_prime_end_soft_clip_spelling_the_reference() {
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

    // Read 17: 10S50M100N50M100N50M from 251. The 10 clipped bases are the last 10 bases of
    // exon 1, which the aligner did not place across the intron. They spell the reference
    // transcript right before the terminal aligned base, so the read keeps them and the
    // stitch is the rest of exon 1.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("17").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_five_prime_end(transcript_model, &end_matches, &gene_annotator, &chromosome_names_map, &fasta_map);

    let end_resolution: EndResolution = end_resolution.unwrap();
    assert_eq!(end_resolution.reference_transcript_match.get_reference_transcript_id(), "ts", "transcript");
    assert_eq!(&*end_resolution.stitch_sequence, "GC".repeat(20), "stitch sequence is exon 1 up to the clip");
    assert_eq!(end_resolution.read_join_position, 0, "the read keeps its clipped bases");
}


#[test]
fn test_resolve_five_prime_end_soft_clip_within_one_edit() {
    // A GC genome and one transcript, tt (101-150, 201-300, 401-500). Three reads align
    // 100M100N100M from 201 after a 10-base 5' soft clip:
    //   "one_edit": the last 10 bases of exon 1 with one substitution (GCGCGAGCGC)
    //   "junk":     TTTTGTTTTG, which is not in tt
    //   "short":    A, one base that is not exon 1's last base (C)
    let directory = tempfile::tempdir().unwrap();
    let genome: String = "GC".repeat(300);
    let fasta_file = directory.path().join("genome.fa");
    std::fs::write(&fasta_file, format!(">chrT\n{genome}\n")).unwrap();
    std::fs::write(directory.path().join("genome.fa.fai"), "chrT\t600\t6\t600\t601\n").unwrap();
    let annotation_file = directory.path().join("genes.tsv");
    std::fs::write(&annotation_file, [
        "row_type\tgene_id\ttranscript_id\texon_id\texon_number\tstrand\tchromosome\tstart\tend",
        "gene\tgt\t\t\t\t+\tchrT\t101\t500",
        "transcript\tgt\ttt\t\t\t+\tchrT\t101\t500",
        "exon\tgt\ttt\ttt_e1\t1\t+\tchrT\t101\t150",
        "exon\tgt\ttt\ttt_e2\t2\t+\tchrT\t201\t300",
        "exon\tgt\ttt\ttt_e3\t3\t+\tchrT\t401\t500"
    ].join("\n") + "\n").unwrap();
    let aligned: String = format!("{}{}", &genome[200..300], &genome[400..500]);
    let sam_text: String = [
        "@HD\tVN:1.6\tSO:unsorted".to_string(),
        "@SQ\tSN:chrT\tLN:600".to_string(),
        format!("one_edit\t0\tchrT\t201\t60\t10S100M100N100M\t*\t0\t0\tGCGCGAGCGC{aligned}\t{}\tcs:Z::100~gc100gc:100", "I".repeat(210)),
        format!("junk\t0\tchrT\t201\t60\t10S100M100N100M\t*\t0\t0\tTTTTGTTTTG{aligned}\t{}\tcs:Z::100~gc100gc:100", "I".repeat(210)),
        format!("short\t0\tchrT\t201\t60\t1S100M100N100M\t*\t0\t0\tA{aligned}\t{}\tcs:Z::100~gc100gc:100", "I".repeat(201))
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

    // Model the reads.
    let (record_positions_map, read_names_map):
        (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(&bam_file);
    let fasta_map: FastaMap = FastaMap::new(fasta_file.to_str().unwrap());
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(annotation_file.to_str().unwrap(), "synthetic", "v1");
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

    // "one_edit" spells exon 1's last 10 bases within one edit, and no shorter or longer
    // stretch comes as close: the read keeps its 10 bases (with the error), and the first
    // 40 bases of exon 1 are stitched before them.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("one_edit").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(transcript_model, &chromosome_names_map, &gene_annotator, &TranscriptTerminus::FivePrime, 1);
    assert_eq!(end_matches.len(), 1, "end match");
    let end_resolution: Option<EndResolution> = resolve_five_prime_end(transcript_model, &end_matches, &gene_annotator, &chromosome_names_map, &fasta_map);
    let end_resolution: EndResolution = end_resolution.unwrap();
    assert_eq!(&*end_resolution.stitch_sequence, "GC".repeat(20), "one_edit: exon 1 up to the clip");
    assert_eq!(end_resolution.read_join_position, 0, "one_edit: the read keeps its clipped bases");

    // "junk" could be the remnant of a fusion arm or of aberrant splicing: left alone.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("junk").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(transcript_model, &chromosome_names_map, &gene_annotator, &TranscriptTerminus::FivePrime, 1);
    assert_eq!(end_matches.len(), 1, "end match");
    let end_resolution: Option<EndResolution> = resolve_five_prime_end(transcript_model, &end_matches, &gene_annotator, &chromosome_names_map, &fasta_map);
    assert!(end_resolution.is_none(), "junk: a clip that is not the reference is left alone");

    // "short" is one edit from both zero and one reference base: the boundary is ambiguous,
    // so the end is left alone.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("short").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(transcript_model, &chromosome_names_map, &gene_annotator, &TranscriptTerminus::FivePrime, 1);
    assert_eq!(end_matches.len(), 1, "end match");
    let end_resolution: Option<EndResolution> = resolve_five_prime_end(transcript_model, &end_matches, &gene_annotator, &chromosome_names_map, &fasta_map);
    assert!(end_resolution.is_none(), "short: an ambiguous boundary is left alone");
}


#[test]
fn test_resolve_five_prime_end_intronic() {
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

    // Read 24: 80M100N50M100N50M from 221. Its first aligned base is 30 bases inside
    // intron 1 of ts, which could be a novel exon, so the end is left alone.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("24").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::FivePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_five_prime_end(transcript_model, &end_matches, &gene_annotator, &chromosome_names_map, &fasta_map);

    assert!(end_resolution.is_none(), "an intronic 5' end is not stitched");
}


#[test]
fn test_resolve_five_prime_end_no_end_match() {
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
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("1").unwrap())
        .unwrap();
    let end_resolution: Option<EndResolution> = resolve_five_prime_end(transcript_model, &[], &gene_annotator, &chromosome_names_map, &fasta_map);

    assert!(end_resolution.is_none(), "no end match, no resolution");
}


#[test]
fn test_resolve_three_prime_end_exonic() {
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

    // Read 3: 50M100N50M100N50M from 101. Its last aligned base is exon 3's last base,
    // so all of exon 4 is appended.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("3").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    let end_resolution: EndResolution = end_resolution.unwrap();
    assert_eq!(end_resolution.terminus, TranscriptTerminus::ThreePrime, "terminus");
    assert_eq!(end_resolution.reference_transcript_match.get_reference_transcript_id(), "ts", "transcript");
    assert_eq!(&*end_resolution.stitch_sequence, "GC".repeat(25), "stitch sequence is exon 4");
    assert_eq!(end_resolution.read_join_position, 149, "read join position");
}


#[test]
fn test_resolve_three_prime_end_partial_exon() {
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

    // Read 4: 25M100N50M100N50M100N25M from 126. Its last aligned base is the 25th base
    // of exon 4, so only the last 25 bases are appended.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("4").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    let end_resolution: EndResolution = end_resolution.unwrap();
    assert_eq!(&*end_resolution.stitch_sequence, "C".to_string() + &"GC".repeat(12), "stitch sequence is the last 25 bases of exon 4");
    assert_eq!(end_resolution.read_join_position, 149, "read join position");
}


#[test]
fn test_resolve_three_prime_end_already_complete_minus_strand() {
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

    // Read 2 (chrM, minus): its last aligned base (301) is the last base of tm.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("2").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    let end_resolution: EndResolution = end_resolution.unwrap();
    assert_eq!(end_resolution.reference_transcript_match.get_reference_transcript_id(), "tm", "transcript");
    assert_eq!(&*end_resolution.stitch_sequence, "", "stitch sequence");
    assert_eq!(end_resolution.read_join_position, 149, "read join position");
}


#[test]
fn test_resolve_three_prime_end_polya_tail() {
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

    // Read 12: 50M100N50M100N50M12S from 101, where the clip is twelve A bases. A polyA
    // tail does not block stitching, and the join sits before it, so it is dropped.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("12").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    let end_resolution: EndResolution = end_resolution.unwrap();
    assert_eq!(&*end_resolution.stitch_sequence, "GC".repeat(25), "stitch sequence is exon 4");
    assert_eq!(end_resolution.read_join_position, 149, "read join position is the last aligned base");

    // Twelve A bases fail a minimum A fraction above 1.0, so the tail now blocks stitching.
    let strict_end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        1.1,
        20,
        6
    );
    assert!(strict_end_resolution.is_none(), "the A fraction threshold is applied to the clip");
}


#[test]
fn test_resolve_three_prime_end_non_adenosine_clip() {
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

    // Read 13: the 3' clip is "GCGCGCGCG", the first nine bases of exon 4.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("13").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    assert!(end_resolution.is_none(), "a 3' clip that is not a polyA tail blocks stitching");

    // Read 14: the 3' clip is "TTTTGTTTTGTTTT".
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("14").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    assert!(end_resolution.is_none(), "a 3' junk clip blocks stitching");
}


#[test]
fn test_resolve_three_prime_end_intronic_tract_downstream() {
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

    // Read 6: 50M100N50M100N70M from 101. It ends at 470, 20 bases inside intron 3 and
    // right before the 471-490 polyA tract. The read carries no polyadenylation signal,
    // so the end is an artifact: roll back to exon 3's last base and append exon 4.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("6").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    let end_resolution: EndResolution = end_resolution.unwrap();
    assert_eq!(end_resolution.reference_transcript_match.get_reference_transcript_id(), "ts", "transcript");
    assert_eq!(&*end_resolution.stitch_sequence, "GC".repeat(25), "stitch sequence is exon 4");
    assert_eq!(end_resolution.read_join_position, 149, "read join position is exon 3's last base, not the terminal base (169)");
}


#[test]
fn test_resolve_three_prime_end_intronic_minus_strand() {
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

    // Read 7 (chrM, minus): 70M100N50M100N50M from 431. It ends at 431, 20 bases inside
    // intron 3 of tm. The read carries AATAAA 16 bases before its end (genomic TTTATT at
    // 441-446), and a polyT tract lies at 411-430, which is polyA on the transcript.
    // Signal plus tract is internal priming: roll back and append exon 4.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("7").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    let end_resolution: EndResolution = end_resolution.unwrap();
    assert_eq!(end_resolution.reference_transcript_match.get_reference_transcript_id(), "tm", "transcript");
    assert_eq!(&*end_resolution.stitch_sequence, "GC".repeat(25), "stitch sequence is exon 4");
    assert_eq!(end_resolution.read_join_position, 149, "read join position");

    // The tract is 20 bases, so a minimum run of 21 no longer sees internal priming. The
    // signal then marks a genuine intronic polyadenylation site, which is left alone.
    // This outcome also results if the genomic window is not reverse complemented.
    let unprimed_end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        21
    );
    assert!(unprimed_end_resolution.is_none(), "a signal without a tract is not stitched");
}


#[test]
fn test_resolve_three_prime_end_intronic_polyadenylation() {
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

    // Read 15: 50M100N70M from 101. It ends at 320, inside intron 2, 20 bases after the
    // start of the AATAAA at 301-306, with no A run of 6 within 20 bases. A signal with
    // no tract is a genuine intronic polyadenylation site, so the end is left alone.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("15").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    assert!(end_resolution.is_none(), "intronic polyadenylation is not stitched");

    // Without a signal to search for, the end is still left alone: only internal priming
    // (an A run in the genome) shows that an intronic end is an artifact, and rolling back
    // would delete the read's 20 intronic bases.
    let unsignalled_end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &[],
        &(10..=40),
        0.8,
        20,
        6
    );
    assert!(unsignalled_end_resolution.is_none(), "an intronic end without a tract is not rolled back");
}


#[test]
fn test_resolve_three_prime_end_exonic_polyadenylation_site() {
    // A GC genome with AATAAA at 251-256, inside exon 2 of tt (101-150, 201-300, 401-500).
    // Read "tailed": 50M50N80M15S from 101, ending at 280 inside exon 2 with a soft-clipped
    // 15-base A tail. The genome holds no A run near 280, so the tail is untemplated.
    // Read "tailless": the same alignment without the tail.
    let directory = tempfile::tempdir().unwrap();
    let mut genome: Vec<u8> = "GC".repeat(300).into_bytes();
    genome[250..256].copy_from_slice(b"AATAAA");
    let genome: String = String::from_utf8(genome).unwrap();
    let fasta_file = directory.path().join("genome.fa");
    std::fs::write(&fasta_file, format!(">chrT\n{genome}\n")).unwrap();
    std::fs::write(directory.path().join("genome.fa.fai"), "chrT\t600\t6\t600\t601\n").unwrap();
    let annotation_file = directory.path().join("genes.tsv");
    std::fs::write(&annotation_file, [
        "row_type\tgene_id\ttranscript_id\texon_id\texon_number\tstrand\tchromosome\tstart\tend",
        "gene\tgt\t\t\t\t+\tchrT\t101\t500",
        "transcript\tgt\ttt\t\t\t+\tchrT\t101\t500",
        "exon\tgt\ttt\ttt_e1\t1\t+\tchrT\t101\t150",
        "exon\tgt\ttt\ttt_e2\t2\t+\tchrT\t201\t300",
        "exon\tgt\ttt\ttt_e3\t3\t+\tchrT\t401\t500"
    ].join("\n") + "\n").unwrap();
    let aligned: String = format!("{}{}", &genome[100..150], &genome[200..280]);
    let sam_text: String = [
        "@HD\tVN:1.6\tSO:unsorted".to_string(),
        "@SQ\tSN:chrT\tLN:600".to_string(),
        format!("tailed\t0\tchrT\t101\t60\t50M50N80M15S\t*\t0\t0\t{aligned}{}\t{}\tcs:Z::50~gc50gc:80", "A".repeat(15), "I".repeat(145)),
        format!("tailless\t0\tchrT\t101\t60\t50M50N80M\t*\t0\t0\t{aligned}\t{}\tcs:Z::50~gc50gc:80", "I".repeat(130))
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

    // Model the reads.
    let (record_positions_map, read_names_map):
        (HashMap<usize, Vec<VirtualPosition>>, BiMap<Box<str>, usize>) = index_bam_records(&bam_file, true, 2);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(&bam_file);
    let fasta_map: FastaMap = FastaMap::new(fasta_file.to_str().unwrap());
    let gene_annotator: TsvGeneAnnotator = TsvGeneAnnotator::new(annotation_file.to_str().unwrap(), "synthetic", "v1");
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
    let tailed: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("tailed").unwrap())
        .unwrap();
    let tailless: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("tailless").unwrap())
        .unwrap();
    let tailed_end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(tailed, &chromosome_names_map, &gene_annotator, &TranscriptTerminus::ThreePrime, 1);
    let tailless_end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(tailless, &chromosome_names_map, &gene_annotator, &TranscriptTerminus::ThreePrime, 1);
    assert!(!tailed_end_matches.is_empty(), "end match");
    assert!(!tailless_end_matches.is_empty(), "end match");

    // The signal starts at 251 and the read ends at 280: the distance from the hexamer's
    // first base to the boundary after the read's last base is 30. A signal, then an
    // untemplated tail, is a polyadenylation site: the end is left alone although exon 2
    // goes on.
    let on_offset: Option<EndResolution> = resolve_three_prime_end(tailed, &tailed_end_matches, &gene_annotator, &chromosome_names_map, &fasta_map, 40, &["AATAAA", "ATTAAA"], &(30..=30), 0.8, 20, 6);
    assert!(on_offset.is_none(), "an offset range of exactly 30 finds the signal");

    // Without the signal in range, the tail alone does not say the end is real: stitch the
    // rest of exon 2 and exon 3.
    let off_offset: Option<EndResolution> = resolve_three_prime_end(tailed, &tailed_end_matches, &gene_annotator, &chromosome_names_map, &fasta_map, 40, &["AATAAA", "ATTAAA"], &(31..=40), 0.8, 20, 6);
    let off_offset: EndResolution = off_offset.unwrap();
    assert_eq!(off_offset.stitch_sequence.len(), 120, "an offset range from 31 misses the signal, so the end is stitched");
    assert_eq!(off_offset.read_join_position, 129, "read join position is the terminal aligned base");

    // An A run in the genome would template the tail (internal priming); here the AAA of
    // the signal (254-256) counts once the window reaches it (30 bases) and the minimum run is 3.
    let primed: Option<EndResolution> = resolve_three_prime_end(tailed, &tailed_end_matches, &gene_annotator, &chromosome_names_map, &fasta_map, 40, &["AATAAA", "ATTAAA"], &(10..=40), 0.8, 30, 3);
    assert!(primed.is_some(), "an internally primed end is stitched");

    // Without a tail, the end is a truncation that happens to follow a signal.
    let tailless_end_resolution: Option<EndResolution> = resolve_three_prime_end(tailless, &tailless_end_matches, &gene_annotator, &chromosome_names_map, &fasta_map, 40, &["AATAAA", "ATTAAA"], &(10..=40), 0.8, 20, 6);
    assert!(tailless_end_resolution.is_some(), "an end without a tail is stitched");
}


#[test]
fn test_resolve_three_prime_end_intronic_tract_upstream() {
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

    // Read 28: 50M100N114M from 101. It ends at 364, the LAST base of the 345-364 polyA
    // tract, as if its tail had aligned onto the tract. AATAAA (331-336) starts 34 bases
    // before the read's end. The tract lies entirely upstream of the terminal aligned
    // base, so only a window that spans the base sees it.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("28").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    let end_resolution: EndResolution = end_resolution.unwrap();
    assert_eq!(&*end_resolution.stitch_sequence, "GC".repeat(50), "stitch sequence is exons 3 and 4");
    assert_eq!(end_resolution.read_join_position, 99, "read join position is exon 2's last base");

    // Without the tract the signal marks a genuine site. A window that only looks
    // downstream of the terminal base gives this outcome.
    let unprimed_end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        21
    );
    assert!(unprimed_end_resolution.is_none(), "a signal without a tract is not stitched");
}


#[test]
fn test_resolve_three_prime_end_rollback_onto_mismatch_base() {
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

    // Read 22: the alignment of read 15 with a mismatch (T) on exon 2's last base (300),
    // which is the roll-back target. A mismatch base is placed, so the roll-back resolves.
    // A roll-back needs internal priming: the window is widened to 30 bases so that it
    // reaches the A tract at 345.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("22").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &[],
        &(10..=40),
        0.8,
        30,
        6
    );

    let end_resolution: EndResolution = end_resolution.unwrap();
    assert_eq!(end_resolution.read_join_position, 99, "read join position");
    assert_eq!(
        transcript_model.get_alignment_model().get_base(99).get_nucleotide().as_str(),
        "T",
        "the read keeps its own base at the join"
    );
    assert_eq!(&*end_resolution.stitch_sequence, "GC".repeat(50), "stitch sequence is exons 3 and 4");
}


#[test]
fn test_resolve_three_prime_end_rollback_across_event() {
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

    // Read 18: 50M100N50M20N40M from 101. It splices from exon 2 into a cassette exon
    // (321-360) inside intron 2 of ts. Rolling back to exon 2 would trim across that
    // splice and delete the cassette exon, so the end is left alone.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("18").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        1
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &[],
        &(10..=40),
        0.8,
        20,
        6
    );

    assert!(end_resolution.is_none(), "a roll-back never crosses an alignment event");
}


#[test]
fn test_resolve_three_prime_end_read_through_arm() {
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

    // Read 20: 50M100N50M350N70M from 101. It splices from gs exon 2 into gr exon 1 and
    // ends at 720, 20 bases inside gr intron 1. The 3' end takes tr, which matches no
    // junction, so the minimum is 0 here. gr intron 1 holds no A run, so nothing shows the
    // end is an artifact and it is left alone.
    let transcript_model: &TranscriptModel = transcript_models
        .iter()
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("20").unwrap())
        .unwrap();
    let end_matches: Vec<EndMatch> = identify_transcript_terminus_matches(
        transcript_model,
        &chromosome_names_map,
        &gene_annotator,
        &TranscriptTerminus::ThreePrime,
        0
    );
    assert!(!end_matches.is_empty(), "end match");
    let end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    assert!(end_resolution.is_none(), "an intronic end without a tract is not rolled back");

    // With a minimum run of 0 every end counts as internally primed. The arm then rolls
    // back to gr exon 1's last base.
    let end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &end_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        0
    );
    let end_resolution: EndResolution = end_resolution.unwrap();
    assert_eq!(end_resolution.reference_transcript_match.get_reference_transcript_id(), "tr", "transcript");
    assert_eq!(&*end_resolution.stitch_sequence, "GC".repeat(25), "stitch sequence is gr exon 2");
    assert_eq!(end_resolution.read_join_position, 149, "read join position");
}


#[test]
fn test_resolve_three_prime_end_no_end_match() {
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
        .find(|model| model.get_read_id() == *read_names_map.get_by_left("3").unwrap())
        .unwrap();
    let end_resolution: Option<EndResolution> = resolve_three_prime_end(
        transcript_model,
        &[],
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map,
        40,
        &["AATAAA", "ATTAAA"],
        &(10..=40),
        0.8,
        20,
        6
    );

    assert!(end_resolution.is_none(), "no end match, no resolution");
}
