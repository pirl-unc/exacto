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

use super::*;


#[test]
fn scga_mini_rna_001_transcript_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    // The simulated TP53 (ENST00000269305.9, minus strand) transcript carrying the SNV at
    // chr17:7674225, error-free: one reverse record, `...343N44=1X65=568N...`, 11 blocks, 10 introns.
    let read_name: &str = "scga-mini-rna-001-tumor-1";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    assert!(alignment.get_records().len() == 1);
    assert_eq!(alignment.get_records().get(0).unwrap().reference_strand, Strand::Reverse);

    let mut num_mismatch: usize = 0;
    for base in alignment.get_bases() {
        if *base.get_kind() == AlignmentModelBaseKind::Mismatch {
            num_mismatch += 1;
        }
    }

    let alignment_model: AlignmentModel = alignment.clone();

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &identify_transcript_model_exons(&alignment_model),
        &identify_transcript_model_splice_junctions(&alignment_model),
        &gene_annotator,
        &chromosome_names_map
    );

    let mut transcript_model: TranscriptModel = TranscriptModel::new(
        alignment_model,
        reference_transcript_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map
    );

    let exons: &Vec<TranscriptModelExon> = transcript_model.get_exons();
    let splice_junctions: &Vec<TranscriptModelSpliceJunction> = transcript_model.get_splice_junctions();

    assert_eq!(num_mismatch, 1);
    assert_eq!(alignment.is_spliced(), true);
    assert_eq!(exons.len(), 11);
    assert_eq!(splice_junctions.len(), 10);
}

#[test]
fn scga_mini_rna_002_transcript_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-002-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-002-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    // The simulated TP53 transcript carrying the 12 bp insertion at chr17:7674224-7674225,
    // error-free: one reverse record, `...343N44=12I66=568N...`, 11 blocks, 10 introns.
    let read_name: &str = "scga-mini-rna-002-tumor-1";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    assert!(alignment.get_records().len() == 1);
    assert_eq!(alignment.get_records().get(0).unwrap().reference_strand, Strand::Reverse);

    let mut num_insertion: usize = 0;
    for base in alignment.get_bases() {
        if *base.get_kind() == AlignmentModelBaseKind::Insertion {
            num_insertion += 1;
        }
    }

    let alignment_model: AlignmentModel = alignment.clone();

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &identify_transcript_model_exons(&alignment_model),
        &identify_transcript_model_splice_junctions(&alignment_model),
        &gene_annotator,
        &chromosome_names_map
    );

    let mut transcript_model: TranscriptModel = TranscriptModel::new(
        alignment_model,
        reference_transcript_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map
    );

    let exons: &Vec<TranscriptModelExon> = transcript_model.get_exons();
    let splice_junctions: &Vec<TranscriptModelSpliceJunction> = transcript_model.get_splice_junctions();

    assert_eq!(num_insertion, 12);
    assert_eq!(alignment.is_spliced(), true);
    assert_eq!(exons.len(), 11);
    assert_eq!(splice_junctions.len(), 10);
}

#[test]
fn scga_mini_rna_003_transcript_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-003-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-003-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    // The simulated TP53 transcript carrying the 30 bp deletion at chr17:7674201-7674230,
    // error-free: one reverse record, `...343N20=30D60=568N...`, 11 blocks, 10 introns.
    let read_name: &str = "scga-mini-rna-003-tumor-1";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment_model: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    let mut num_deletion: usize = 0;
    for event in alignment_model.get_events().values() {
        if event.get_kind() == &AlignmentModelEventKind::Deletion {
            num_deletion += 1;
        }
    }

    assert_eq!(alignment_model.get_records().len(), 1);
    assert_eq!(alignment_model.get_records().get(0).unwrap().reference_strand, Strand::Reverse);
    assert_eq!(alignment_model.is_spliced(), true);
    assert_eq!(num_deletion, 1);

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &identify_transcript_model_exons(&alignment_model),
        &identify_transcript_model_splice_junctions(&alignment_model),
        &gene_annotator,
        &chromosome_names_map
    );

    let transcript_model: TranscriptModel = TranscriptModel::new(
        alignment_model,
        reference_transcript_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map
    );

    let exons: &Vec<TranscriptModelExon> = transcript_model.get_exons();
    let splice_junctions: &Vec<TranscriptModelSpliceJunction> = transcript_model.get_splice_junctions();

    assert_eq!(exons.len(), 11);
    assert_eq!(splice_junctions.len(), 10);
}

#[test]
fn scga_mini_rna_007_transcript_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-007-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-007-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    // The simulated ASPA-WSCD1 fusion transcript, error-free. minimap2 splits it in two: the
    // ASPA part (read 0-797, `399=5207N196=1700N94=5642N109=5066S`) and the WSCD1 part (read
    // 795-5863, `795S117=2216N...2120N4119=`). The three bases both place, read 795-797, are the
    // ones the two genes share at the junction.
    let read_name: &str = "scga-mini-rna-007-tumor-1";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment_model: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    assert_eq!(alignment_model.get_records().len(), 2);
    assert_eq!(alignment_model.get_records()[0].read_start, 0);
    assert_eq!(alignment_model.get_records()[0].read_end, 797);
    assert_eq!(alignment_model.get_records()[0].reference_strand, Strand::Forward);
    assert_eq!(alignment_model.get_records()[1].read_start, 795);
    assert_eq!(alignment_model.get_records()[1].read_end, 5863);
    assert_eq!(alignment_model.get_records()[1].reference_strand, Strand::Forward);
    assert_eq!(alignment_model.is_spliced(), true);

    let mut num_mismatch: usize = 0;
    for base in alignment_model.get_bases() {
        if *base.get_kind() == AlignmentModelBaseKind::Mismatch {
            num_mismatch += 1;
        }
    }

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &identify_transcript_model_exons(&alignment_model),
        &identify_transcript_model_splice_junctions(&alignment_model),
        &gene_annotator,
        &chromosome_names_map
    );

    let transcript_model: TranscriptModel = TranscriptModel::new(
        alignment_model,
        reference_transcript_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map
    );

    let annotation: TranscriptModelAnnotation = transcript_model.get_annotation().clone();

    let mut num_fusion_gene: usize = 0;
    for event in transcript_model.get_alignment_model().get_events().values() {
        if annotation
            .get_event(event.get_prev_read_position(), event.get_next_read_position())
            .and_then(|e| e.get_context().as_ref()) == Some(&AlignmentModelEventContext::FusionGene) {
            num_fusion_gene += 1;
        }
    }

    let options: IdentifyRNATranscriptVariantsOptions = IdentifyRNATranscriptVariantsOptions::default();
    let breakpoint_rescue: Option<RNABreakpointRescue> = options.calling.bkpt_rescue.then(|| RNABreakpointRescue {
        gene_annotator: &gene_annotator,
        chromosome_names_map: &chromosome_names_map,
        fasta_map: &fasta_map,
        min_insertion_length: options.calling.bkpt_rescue_min_ins_len,
        gap_open: options.calling.bkpt_rescue_realignment_gap_open_score,
        gap_extend: options.calling.bkpt_rescue_realignment_gap_extend_score,
        k: options.calling.bkpt_rescue_realignment_k,
        band_width: options.calling.bkpt_rescue_realignment_band_width,
        min_score_fraction: options.calling.bkpt_rescue_realignment_min_score_fraction,
        min_query_coverage: options.calling.bkpt_rescue_realignment_min_query_coverage,
        min_placed_fraction: options.calling.bkpt_rescue_realignment_min_placed_fraction,
        max_pieces: options.calling.bkpt_rescue_realignment_max_pieces
    });
    let fusion_gene_exists: bool = RNAVariantRecordCaller::new(
        options.calling.min_mapping_quality,
        0,
        options.calling.min_terminal_soft_clip_ins_len,
        breakpoint_rescue
    )
        .call(&transcript_model)
        .iter()
        .any(|variant_record| *variant_record.get_variant_type() == VariantType::FusionGene);

    let exons: &Vec<TranscriptModelExon> = transcript_model.get_exons();
    let splice_junctions: &Vec<TranscriptModelSpliceJunction> = transcript_model.get_splice_junctions();

    // Error-free, so no mismatch. The one fusion is the breakend between the two records; the
    // exons are ASPA's 4 and WSCD1's 7, the junctions the 3 + 6 introns inside them.
    assert_eq!(num_mismatch, 0);
    assert_eq!(num_fusion_gene, 1);
    assert_eq!(fusion_gene_exists, true);
    assert_eq!(exons.len(), 11);
    assert_eq!(splice_junctions.len(), 9);
}

#[test]
fn scga_mini_rna_005_transcript_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-005-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-005-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    // The simulated TP53 transcript with the 1,100 bp non-canonical splice: one reverse record,
    // `...92N100=1100N71=81N...`. The junction chr17:7673801-7674900 takes out the last 37 bases
    // of the 137 bp exon, the whole 110 bp exon and the first 42 bases of the 113 bp exon.
    let read_name: &str = "scga-mini-rna-005-tumor-1";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    assert!(alignment.get_records().len() == 1);
    assert_eq!(alignment.get_records().get(0).unwrap().reference_strand, Strand::Reverse);

    let exons: Vec<TranscriptModelExon> = identify_transcript_model_exons(&alignment);
    let splice_junctions: Vec<TranscriptModelSpliceJunction> = identify_transcript_model_splice_junctions(&alignment);

    let reference_transcript: &Transcript = gene_annotator.get_transcript("ENST00000269305.9").unwrap();

    // `identify_transcript_model_annotation` runs inside `TranscriptModel::new`,
    // so the overlay is read back off the model. The model builds the
    // ReferenceTranscriptSequence itself; the test supplies only the match.
    let transcript_model: TranscriptModel = TranscriptModel::new(
        alignment.clone(),
        vec![ReferenceTranscriptMatch::new(
            &*reference_transcript.gene_id,
            &*reference_transcript.gene_id,
            "ENST00000269305.9",
            &Vec::new(),
            0,
            0,
            0,
            0
        )],
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map
    );
    let annotation: TranscriptModelAnnotation = transcript_model.get_annotation().clone();

    let mut num_ref_bases_skipped: usize = 0;
    for event in alignment.get_events().values() {
        if let Some(event_annotation) = annotation.get_event(event.get_prev_read_position(), event.get_next_read_position()) {
            for run in event_annotation.get_skipped_reference_bases().iter() {
                num_ref_bases_skipped += run.len();
            }
        }
    }

    let mut num_noncanonical_splicing: usize = 0;
    for event in alignment.get_events().values() {
        if annotation
            .get_event(event.get_prev_read_position(), event.get_next_read_position())
            .and_then(|e| e.get_context().as_ref()) == Some(&AlignmentModelEventContext::NonCanonicalSplicing) {
            num_noncanonical_splicing += 1;
        }
    }

    // 37 + 110 + 42 exonic bases of ENST00000269305.9 skipped, one unannotated junction, and
    // the 11 exons less the skipped one: 10 exons over 9 junctions.
    assert_eq!(num_ref_bases_skipped, 189);
    assert_eq!(num_noncanonical_splicing, 1);
    assert_eq!(alignment.is_spliced(), true);
    assert_eq!(exons.len(), 10);
    assert_eq!(splice_junctions.len(), 9);
}

#[test]
fn scga_mini_rna_008_transcript_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-008-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-008-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    // The simulated TP53 transcript with the 280 bp cryptic exon chr17:7672205-7672484 in the
    // intron between the 74 bp and 107 bp exons: one reverse record, `...107=1489N280=1050N74=...`.
    let read_name: &str = "scga-mini-rna-008-tumor-1";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    assert!(alignment.get_records().len() == 1);
    assert_eq!(alignment.get_records().get(0).unwrap().reference_strand, Strand::Reverse);

    let exons: Vec<TranscriptModelExon> = identify_transcript_model_exons(&alignment);
    let splice_junctions: Vec<TranscriptModelSpliceJunction> = identify_transcript_model_splice_junctions(&alignment);

    let reference_transcript: &Transcript = gene_annotator.get_transcript("ENST00000269305.9").unwrap();

    let transcript_model: TranscriptModel = TranscriptModel::new(
        alignment.clone(),
        vec![ReferenceTranscriptMatch::new(
            &*reference_transcript.gene_id,
            &*reference_transcript.gene_id,
            "ENST00000269305.9",
            &Vec::new(),
            0,
            0,
            0,
            0
        )],
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map
    );
    let annotation: TranscriptModelAnnotation = transcript_model.get_annotation().clone();

    let mut num_ref_bases_cryptic: usize = 0;
    for base in alignment.get_bases() {
        if annotation.get_base_context(base.get_read_position()) == Some(&AlignmentModelBaseContext::Intronic) {
            num_ref_bases_cryptic += 1;
        }
    }

    let mut num_noncanonical_splicing: usize = 0;
    for event in alignment.get_events().values() {
        if annotation
            .get_event(event.get_prev_read_position(), event.get_next_read_position())
            .and_then(|e| e.get_context().as_ref()) == Some(&AlignmentModelEventContext::NonCanonicalSplicing) {
            num_noncanonical_splicing += 1;
        }
    }

    // Every base of the cryptic exon is intronic to ENST00000269305.9, and neither junction
    // into or out of it is annotated. 11 exons plus the cryptic one.
    assert_eq!(num_ref_bases_cryptic, 280);
    assert_eq!(num_noncanonical_splicing, 2);
    assert_eq!(alignment.is_spliced(), true);
    assert_eq!(exons.len(), 12);
    assert_eq!(splice_junctions.len(), 11);
}

#[test]
fn scga_mini_rna_009_transcript_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-009-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-009-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    // The simulated TP53 transcript retaining the first 8 bases of the intron after the 74 bp
    // exon, chr17:7673527-7673534: one reverse record, `...107=2811N82=92N...`, the 82 bp block
    // being the 74 bp exon and the 8 retained bases.
    let read_name: &str = "scga-mini-rna-009-tumor-1";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    assert!(alignment.get_records().len() == 1);
    assert_eq!(alignment.get_records().get(0).unwrap().reference_strand, Strand::Reverse);

    let exons: Vec<TranscriptModelExon> = identify_transcript_model_exons(&alignment);
    let splice_junctions: Vec<TranscriptModelSpliceJunction> = identify_transcript_model_splice_junctions(&alignment);

    let reference_transcript: &Transcript = gene_annotator.get_transcript("ENST00000269305.9").unwrap();

    let transcript_model: TranscriptModel = TranscriptModel::new(
        alignment.clone(),
        vec![ReferenceTranscriptMatch::new(
            &*reference_transcript.gene_id,
            &*reference_transcript.gene_id,
            "ENST00000269305.9",
            &Vec::new(),
            0,
            0,
            0,
            0
        )],
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map
    );
    let annotation: TranscriptModelAnnotation = transcript_model.get_annotation().clone();

    let mut num_ref_bases_intron: usize = 0;
    for base in alignment.get_bases() {
        if annotation.get_base_context(base.get_read_position()) == Some(&AlignmentModelBaseContext::Intronic) {
            num_ref_bases_intron += 1;
        }
    }

    let mut num_noncanonical_splicing: usize = 0;
    for event in alignment.get_events().values() {
        if annotation
            .get_event(event.get_prev_read_position(), event.get_next_read_position())
            .and_then(|e| e.get_context().as_ref()) == Some(&AlignmentModelEventContext::NonCanonicalSplicing) {
            num_noncanonical_splicing += 1;
        }
    }

    // The 8 retained bases are intronic, and the junction that now starts after them is not
    // annotated.
    assert_eq!(num_ref_bases_intron, 8);
    assert_eq!(num_noncanonical_splicing, 1);
    assert_eq!(alignment.is_spliced(), true);
    assert_eq!(exons.len(), 11);
    assert_eq!(splice_junctions.len(), 10);
}

#[test]
fn scga_mini_rna_015_transcript_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-015-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-015-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    // The simulated ASPA-WSCD1-ACAP1 transcript (two fusions), error-free. minimap2 splits it in
    // three forward records: ASPA (read 0-797), WSCD1 (read 795-1543) and ACAP1 (read 1542-1959).
    // Neighbouring records both place the bases the genes share at each junction: read 795-797
    // and read 1542-1543.
    let read_name: &str = "scga-mini-rna-015-tumor-1";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    assert!(alignment.get_records().len() == 3);
    assert_eq!(alignment.get_records().get(0).unwrap().reference_strand, Strand::Forward);
    assert_eq!(alignment.get_records().get(1).unwrap().reference_strand, Strand::Forward);
    assert_eq!(alignment.get_records().get(2).unwrap().reference_strand, Strand::Forward);
    assert_eq!((alignment.get_records()[0].read_start, alignment.get_records()[0].read_end), (0, 797));
    assert_eq!((alignment.get_records()[1].read_start, alignment.get_records()[1].read_end), (795, 1543));
    assert_eq!((alignment.get_records()[2].read_start, alignment.get_records()[2].read_end), (1542, 1959));

    let exons: Vec<TranscriptModelExon> = identify_transcript_model_exons(&alignment);
    let splice_junctions: Vec<TranscriptModelSpliceJunction> = identify_transcript_model_splice_junctions(&alignment);

    // The transcripts the simulation built the three parts from.
    let reference_transcript_1: &Transcript = gene_annotator.get_transcript("ENST00000263080.3").unwrap();
    let reference_transcript_2: &Transcript = gene_annotator.get_transcript("ENST00000317744.10").unwrap();
    let reference_transcript_3: &Transcript = gene_annotator.get_transcript("ENST00000575425.1").unwrap();

    let transcript_model: TranscriptModel = TranscriptModel::new(
        alignment.clone(),
        vec![
            ReferenceTranscriptMatch::new(
                &*reference_transcript_1.gene_id,
                &*reference_transcript_1.gene_id,
                "ENST00000263080.3",
                &Vec::new(),
                0,
                0,
                0,
                0
            ),
            ReferenceTranscriptMatch::new(
                &*reference_transcript_2.gene_id,
                &*reference_transcript_2.gene_id,
                "ENST00000317744.10",
                &Vec::new(),
                0,
                0,
                0,
                0
            ),
            ReferenceTranscriptMatch::new(
                &*reference_transcript_3.gene_id,
                &*reference_transcript_3.gene_id,
                "ENST00000575425.1",
                &Vec::new(),
                0,
                0,
                0,
                0
            )
        ],
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map
    );
    let annotation: TranscriptModelAnnotation = transcript_model.get_annotation().clone();

    let mut num_fusion_gene: usize = 0;
    for event in alignment.get_events().values() {
        if annotation
            .get_event(event.get_prev_read_position(), event.get_next_read_position())
            .and_then(|e| e.get_context().as_ref()) == Some(&AlignmentModelEventContext::FusionGene) {
            num_fusion_gene += 1;
        }
    }

    let variant_records: Vec<VariantRecord> = RNAVariantRecordCaller::new(30, 30, 0, None).call(&transcript_model);

    // 4 + 5 + 4 exons; the 3 + 4 + 3 introns inside the records are the junctions, the two
    // breakends between records are not.
    assert_eq!(num_fusion_gene, 2);
    assert_eq!(alignment.is_spliced(), true);
    assert_eq!(exons.len(), 13);
    assert_eq!(splice_junctions.len(), 10);

    // A breakend between overlapping records is stated over the bases either side of the
    // overlap, and the overlap is its sequence: read 794 (ASPA chr17:3489340) to read 798
    // (WSCD1 chr17:6087991) over `AGG`, and read 1541 (WSCD1 chr17:6110933) to read 1544
    // (ACAP1 chr17:7341948) over `AG`.
    let fusion_genes: Vec<&VariantRecord> = variant_records
        .iter()
        .filter(|variant_record| *variant_record.get_variant_type() == VariantType::FusionGene)
        .collect();
    assert_eq!(fusion_genes.len(), 2);
    assert_eq!(fusion_genes[0].get_chromosome_1(), 0);
    assert_eq!(fusion_genes[0].get_chromosome_2(), 0);
    assert_eq!(fusion_genes[0].get_position_1(), 3489340);
    assert_eq!(fusion_genes[0].get_position_2(), 6087991);
    assert_eq!(fusion_genes[0].get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(fusion_genes[0].get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(fusion_genes[0].get_sequence(), "AGG"); // overlapping alignment
    assert_eq!(fusion_genes[1].get_chromosome_1(), 0);
    assert_eq!(fusion_genes[1].get_chromosome_2(), 0);
    assert_eq!(fusion_genes[1].get_position_1(), 6110933);
    assert_eq!(fusion_genes[1].get_position_2(), 7341948);
    assert_eq!(fusion_genes[1].get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(fusion_genes[1].get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(fusion_genes[1].get_sequence(), "AG"); // overlapping alignment

    // Every other record is an exon truncation: the exons of the three transcripts the read
    // does not place, each run going to a breakend. ASPA: its last two exons (110 + 4,515) and
    // chr17:3489341-3489342 under the overlap. WSCD1: its first two exons (274 + 715), its last
    // two (201 + 4,119), chr17:6087990 and chr17:6110934-6110935 under the overlaps. ACAP1: its
    // first exon (198). Every junction is annotated, so none is non-canonical.
    assert!(variant_records.iter().all(|variant_record| matches!(
        variant_record.get_variant_type(),
        VariantType::FusionGene | VariantType::ExonTruncation
    )));
    let num_ref_bases_skipped: usize = variant_records
        .iter()
        .filter(|variant_record| *variant_record.get_variant_type() == VariantType::ExonTruncation)
        .map(|variant_record| variant_record.get_sequence().len())
        .sum();
    assert_eq!(num_ref_bases_skipped, (110 + 4515 + 2) + (274 + 715 + 201 + 4119 + 1 + 2) + 198);
}

#[test]
fn scga_mini_rna_016_transcript_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-016-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-016-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    // The simulated ASPA-WSCD1 fusion transcript with 12 untemplated bases at the junction,
    // error-free. minimap2 splits it in two forward records, ASPA (read 0-796,
    // `399=5207N196=1700N94=5642N108=5079S`) and WSCD1 (read 808-5875, `808S116=2216N...`), and
    // places neither on read 797-807. The 12th inserted base, a G, is placed by the WSCD1 record
    // on chr17:6087989, the intron base just before the exon at 6087990.
    let read_name: &str = "scga-mini-rna-016-tumor-1";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    assert!(alignment.get_records().len() == 2);
    assert_eq!(alignment.get_records().get(0).unwrap().reference_strand, Strand::Forward);
    assert_eq!(alignment.get_records().get(1).unwrap().reference_strand, Strand::Forward);
    assert_eq!((alignment.get_records()[0].read_start, alignment.get_records()[0].read_end), (0, 796));
    assert_eq!((alignment.get_records()[1].read_start, alignment.get_records()[1].read_end), (808, 5875));

    let exons: Vec<TranscriptModelExon> = identify_transcript_model_exons(&alignment);
    let splice_junctions: Vec<TranscriptModelSpliceJunction> = identify_transcript_model_splice_junctions(&alignment);

    // The transcripts the simulation built the two parts from.
    let reference_transcript_1: &Transcript = gene_annotator.get_transcript("ENST00000263080.3").unwrap();
    let reference_transcript_2: &Transcript = gene_annotator.get_transcript("ENST00000317744.10").unwrap();

    let transcript_model: TranscriptModel = TranscriptModel::new(
        alignment.clone(),
        vec![
            ReferenceTranscriptMatch::new(
                &*reference_transcript_1.gene_id,
                &*reference_transcript_1.gene_id,
                "ENST00000263080.3",
                &Vec::new(),
                0,
                0,
                0,
                0
            ),
            ReferenceTranscriptMatch::new(
                &*reference_transcript_2.gene_id,
                &*reference_transcript_2.gene_id,
                "ENST00000317744.10",
                &Vec::new(),
                0,
                0,
                0,
                0
            )
        ],
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map
    );
    let annotation: TranscriptModelAnnotation = transcript_model.get_annotation().clone();

    let mut num_fusion_gene: usize = 0;
    for event in alignment.get_events().values() {
        if annotation
            .get_event(event.get_prev_read_position(), event.get_next_read_position())
            .and_then(|e| e.get_context().as_ref()) == Some(&AlignmentModelEventContext::FusionGene) {
            num_fusion_gene += 1;
        }
    }

    let variant_records: Vec<VariantRecord> = RNAVariantRecordCaller::new(30, 30, 0, None).call(&transcript_model);

    // ASPA's 4 exons and WSCD1's 7; the 3 + 6 introns inside the records are the junctions.
    assert_eq!(num_fusion_gene, 1);
    assert_eq!(alignment.is_spliced(), true);
    assert_eq!(exons.len(), 11);
    assert_eq!(splice_junctions.len(), 9);

    // The breakend runs from ASPA's last placed base (read 796, chr17:3489342) to WSCD1's first
    // (read 808, chr17:6087989) and carries the 11 unplaced bases between them.
    let fusion_genes: Vec<&VariantRecord> = variant_records
        .iter()
        .filter(|variant_record| *variant_record.get_variant_type() == VariantType::FusionGene)
        .collect();
    assert_eq!(fusion_genes.len(), 1);
    assert_eq!(fusion_genes[0].get_chromosome_1(), 0);
    assert_eq!(fusion_genes[0].get_chromosome_2(), 0);
    assert_eq!(fusion_genes[0].get_position_1(), 3489342);
    assert_eq!(fusion_genes[0].get_position_2(), 6087989);
    assert_eq!(fusion_genes[0].get_operation_1(), &GraphOperationType::Downstream);
    assert_eq!(fusion_genes[0].get_operation_2(), &GraphOperationType::Upstream);
    assert_eq!(fusion_genes[0].get_sequence(), "CCCATCCGCCT");

    // Every junction lies in an annotated intron, so none is non-canonical. The G at chr17:6087989
    // abuts the WSCD1 exon, so it reads as a 1 bp intron retention.
    assert_eq!(
        variant_records.iter().filter(|variant_record| *variant_record.get_variant_type() == VariantType::NonCanonicalSplicing).count(),
        0
    );
    let intron_retentions: Vec<&VariantRecord> = variant_records
        .iter()
        .filter(|variant_record| *variant_record.get_variant_type() == VariantType::IntronRetention)
        .collect();
    assert_eq!(intron_retentions.len(), 1);
    assert_eq!(intron_retentions[0].get_read_position_1(), 808);
    assert_eq!(intron_retentions[0].get_read_position_2(), 808);

    // Every other record is an exon truncation: ASPA's last two exons (110 + 4,515) and WSCD1's
    // first two (274 + 715), the exons the simulation left out of the fusion.
    assert!(variant_records.iter().all(|variant_record| matches!(
        variant_record.get_variant_type(),
        VariantType::FusionGene | VariantType::IntronRetention | VariantType::ExonTruncation
    )));
    let num_ref_bases_skipped: usize = variant_records
        .iter()
        .filter(|variant_record| *variant_record.get_variant_type() == VariantType::ExonTruncation)
        .map(|variant_record| variant_record.get_sequence().len())
        .sum();
    assert_eq!(num_ref_bases_skipped, (110 + 4515) + (274 + 715));
}

#[test]
fn scga_mini_rna_011_transcript_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-011-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-011-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    // The simulated TP53 circular RNA, two laps of the 110, 137, 74 and 107 bp exons joined by
    // the back-splice chr17:7670609 -> 7674290, error-free. minimap2 splits it in three reverse
    // records: the first lap's first three exons (read 0-320), its 107 bp exon (read 318-429) and
    // the second lap (read 426-855). The records overlap on read 318-320 and read 426-429.
    let read_name: &str = "scga-mini-rna-011-tumor-1";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    assert!(alignment.get_records().len() == 3);
    assert_eq!(alignment.get_records().get(0).unwrap().reference_strand, Strand::Reverse);
    assert_eq!(alignment.get_records().get(1).unwrap().reference_strand, Strand::Reverse);
    assert_eq!(alignment.get_records().get(2).unwrap().reference_strand, Strand::Reverse);
    assert_eq!((alignment.get_records()[0].read_start, alignment.get_records()[0].read_end), (0, 320));
    assert_eq!((alignment.get_records()[1].read_start, alignment.get_records()[1].read_end), (318, 429));
    assert_eq!((alignment.get_records()[2].read_start, alignment.get_records()[2].read_end), (426, 855));

    let reference_transcript: &Transcript = gene_annotator.get_transcript("ENST00000269305.9").unwrap();

    let transcript_model: TranscriptModel = TranscriptModel::new(
        alignment.clone(),
        vec![ReferenceTranscriptMatch::new(
            &*reference_transcript.gene_id,
            &*reference_transcript.gene_id,
            "ENST00000269305.9",
            &Vec::new(),
            0,
            0,
            0,
            0
        )],
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map
    );
    let annotation: TranscriptModelAnnotation = transcript_model.get_annotation().clone();

    // The back-splice is the breakend between the second and the third record, stated over the
    // bases either side of their overlap: read 425 (first lap, 107 bp exon) and read 430 (second
    // lap, 110 bp exon). The breakend between the first and the second record, read 317 to 321,
    // is the canonical junction 7673535 -> 7670715 of the first lap; it is not pinned here: the
    // model types it a back-splice too (open bug, its two sides share the second lap's bases).
    assert_eq!(
        annotation.get_event(425, 430).and_then(|e| e.get_context().as_ref()),
        Some(&AlignmentModelEventContext::BackSplicing)
    );
    assert_eq!(alignment.is_spliced(), true);
}

#[test]
fn scga_mini_rna_014_transcript_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-014-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-014-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();

    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let (record_positions_map, read_names_map) = index_bam_records(bam_file, true, 2);

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    // A TP53 read (reverse record) whose cs reads `:834~ct918ac+g:107`: an inserted base against
    // the intron chr17:7669691-7670608, between it and the 107 bp exon chr17:7670609-7670715.
    // Every other junction of the read is an annotated TP53 intron.
    let read_name: &str = "scga-mini-rna-014-tumor_chunk_0000/306/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment_model: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    let splice_junctions: Vec<TranscriptModelSpliceJunction> = identify_transcript_model_splice_junctions(&alignment_model);

    // `~ct918ac` + `g`: on the minus strand the intron runs from its donor at chr17:7670608 down
    // to its acceptor at 7669691, and the inserted base sits at the junction rather than inside
    // it. Read off the inserted base, which is placed on the intron's last reference base, the
    // donor would come back as 7670607, one short, with a signal that is not GT.
    let junction: &TranscriptModelSpliceJunction = splice_junctions
        .iter()
        .find(|junction| junction.reference_position_1 == 7670608)
        .expect("intron at chr17:7670608 not reported");

    assert_eq!(junction.reference_position_2, 7669691);

    // The model no longer carries splice site signals, so read them off the reference: the
    // donor is the intron's first two bases and the acceptor its last two, each on the
    // junction's strand.
    let splice_site_signals = |junction: &TranscriptModelSpliceJunction| -> (Box<str>, Box<str>) {
        let chromosome_1: &str = chromosome_names_map.get_by_right(&junction.reference_chromosome_id_1).unwrap();
        let chromosome_2: &str = chromosome_names_map.get_by_right(&junction.reference_chromosome_id_2).unwrap();
        let position_1: usize = junction.reference_position_1 as usize;
        let position_2: usize = junction.reference_position_2 as usize;
        let donor: Box<str> = if junction.reference_strand_1 == Strand::Forward {
            fasta_map.get_sequence(chromosome_1, position_1, position_1 + 1).into()
        } else {
            reverse_complement(fasta_map.get_sequence(chromosome_1, position_1 - 1, position_1))
        };
        let acceptor: Box<str> = if junction.reference_strand_2 == Strand::Forward {
            fasta_map.get_sequence(chromosome_2, position_2 - 1, position_2).into()
        } else {
            reverse_complement(fasta_map.get_sequence(chromosome_2, position_2, position_2 + 1))
        };
        (donor, acceptor)
    };

    let (donor, acceptor): (Box<str>, Box<str>) = splice_site_signals(junction);
    assert_eq!(&*donor, "GT");
    assert_eq!(&*acceptor, "AG");

    // Every junction of this read reads back a canonical signal once the flanks are aligned
    // bases; before the fix this one did not.
    for junction in splice_junctions.iter() {
        let (donor, acceptor): (Box<str>, Box<str>) = splice_site_signals(junction);
        assert_eq!(&*donor, "GT", "junction {} donor", junction.splice_junction_number);
        assert_eq!(&*acceptor, "AG", "junction {} acceptor", junction.splice_junction_number);
    }
}

#[test]
fn scga_mini_rna_013_transcript_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-013-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-013-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());

    let chromosome_names: Vec<Box<str>> = get_chromosome_names(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let (record_positions_map, read_names_map) = index_bam_records(bam_file, true, 2);

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    // A TP53 read of the rna-013 sample split across the two strands: a supplementary forward
    // record (read 0-1563, `52=343N108=568N115=81N7=1180I6728N...140177N1I4=1D23=336S`) and the
    // primary reverse record (read 459-1893, `6S43=1D194=...315=918N107=2819N74=459S`).
    let read_name: &str = "scga-mini-rna-013-tumor_chunk_0000/182/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment_model: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    // Precondition: the two records overlap over read positions 459-1563. The supplementary
    // is sorted first, so the primary re-places every base in the overlap.
    assert_eq!(alignment_model.get_records().len(), 2);
    assert_eq!(alignment_model.get_records()[0].read_start, 0);
    assert_eq!(alignment_model.get_records()[0].read_end, 1563);
    assert_eq!(alignment_model.get_records()[0].reference_strand, Strand::Forward);
    assert_eq!(alignment_model.get_records()[1].read_start, 459);
    assert_eq!(alignment_model.get_records()[1].read_end, 1893);
    assert_eq!(alignment_model.get_records()[1].reference_strand, Strand::Reverse);

    // The invariant the purge exists to hold: an event's two flanking bases are placed by one
    // alignment, so they agree on chromosome and strand and lie at least two apart - one base
    // of reference is skipped, at minimum. The supplementary leaves four events whose flanks
    // the primary re-places: the intron after its 1,180 bp insertion (281,1462), its deletions
    // at (1514,1515) and (1540,1541), and the intron at (1535,1537). Kept, each would join a
    // supplementary placement to a primary one, or two adjacent primary bases.
    for ((read_position_1, read_position_2), event) in alignment_model.get_events().iter() {
        if *event.get_kind() == AlignmentModelEventKind::Breakpoint {
            continue;
        }

        let coordinate_1 = alignment_model.get_bases()[*read_position_1 as usize]
            .get_placement()
            .get_coordinate()
            .expect("event flank is placed");
        let coordinate_2 = alignment_model.get_bases()[*read_position_2 as usize]
            .get_placement()
            .get_coordinate()
            .expect("event flank is placed");

        assert_eq!(coordinate_1.0, coordinate_2.0, "event ({read_position_1},{read_position_2}) spans two chromosomes");
        assert_eq!(coordinate_1.2, coordinate_2.2, "event ({read_position_1},{read_position_2}) spans two strands");
        assert!(
            coordinate_1.1.abs_diff(coordinate_2.1) >= 2,
            "event ({read_position_1},{read_position_2}) skips no reference: {} to {}",
            coordinate_1.1, coordinate_2.1
        );
    }

    // Three splices from the supplementary (its events below the overlap), two splices and seven
    // deletions from the primary, one breakend between them over the overlap.
    assert_eq!(alignment_model.get_events().len(), 13);
    assert_eq!(alignment_model.get_event(281, 1462).is_none(), true);
    assert_eq!(alignment_model.get_event(1514, 1515).is_none(), true);
    assert_eq!(alignment_model.get_event(1535, 1537).is_none(), true);
    assert_eq!(alignment_model.get_event(1540, 1541).is_none(), true);
    assert_eq!(*alignment_model.get_event(458, 1564).unwrap().get_kind(), AlignmentModelEventKind::Breakpoint);

    // The primary's deletions are untouched, among them 3 bp at chr17:7668853-7668855 and
    // 1 bp at chr17:7668472.
    assert_eq!(*alignment_model.get_event(1472, 1473).unwrap().get_kind(), AlignmentModelEventKind::Deletion);
    assert_eq!(alignment_model.get_bases()[1472].get_placement().get_coordinate().unwrap().1, 7668856);
    assert_eq!(alignment_model.get_bases()[1473].get_placement().get_coordinate().unwrap().1, 7668852);
    assert_eq!(*alignment_model.get_event(1850, 1851).unwrap().get_kind(), AlignmentModelEventKind::Deletion);
    assert_eq!(alignment_model.get_bases()[1850].get_placement().get_coordinate().unwrap().1, 7668473);
    assert_eq!(alignment_model.get_bases()[1851].get_placement().get_coordinate().unwrap().1, 7668471);

    let mut num_splices: usize = 0;
    for ((read_position_1, read_position_2), event) in alignment_model.get_events().iter() {
        if *event.get_kind() != AlignmentModelEventKind::Splicing {
            continue;
        }
        num_splices += 1;
        for read_position in [read_position_1, read_position_2] {
            assert!(
                alignment_model.get_bases()[*read_position as usize]
                    .get_deleted_reference_bases()
                    .is_empty(),
                "splice flank at read position {read_position} was given intronic bases to restore"
            );
        }
    }
    assert!(num_splices > 0, "fixture carries no splices; the assertion above proves nothing");

    // Every base the primary owns reports the primary's strand, and none keeps an allele the
    // supplementary wrote (its `*` tokens at read 1504-1506, 1511, 1518 and 1534). The primary's
    // one `*` token is its own mismatch at read 1472.
    for read_position in 459..=1563u32 {
        let base = &alignment_model.get_bases()[read_position as usize];
        let (_, _, strand) = base.get_placement().get_coordinate().expect("overlap base is placed");
        assert_eq!(*strand, Strand::Reverse, "read position {read_position} was not re-placed by the primary");
        if read_position == 1472 {
            assert!(base.get_reference_nucleotide().is_some(), "the primary's mismatch at read position 1472 lost its reference allele");
            continue;
        }
        assert_eq!(
            base.get_reference_nucleotide(),
            None,
            "read position {read_position} kept a stale reference allele from the supplementary"
        );
    }

    let read_name: &str = "scga-mini-rna-013-tumor_chunk_0000/182/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());
    let alignment_model: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|r| Arc::new(r.clone())).collect::<Vec<_>>()
    );

    let mut lengths: Vec<usize> = Vec::new();
    for base in alignment_model.get_bases().iter() {
        if base.get_deleted_reference_bases().is_empty() {
            continue;
        }
        let (chromosome_id, position, strand) = base
            .get_placement()
            .get_coordinate()
            .expect("a deletion flank is placed");
        let contig: &str = &chromosome_names[chromosome_id as usize];
        let length: usize = base.get_deleted_reference_bases().len();
        lengths.push(length);

        // The cs walk advances the reference regardless of strand, so the deleted bases always
        // follow the flank on the reference; only the read-orientation spelling differs.
        let reference: String = fasta_map
            .get_sequence(contig, position as usize + 1, position as usize + length)
            .to_uppercase();
        let expected: String = if *strand == Strand::Reverse {
            reference.chars().rev().map(|c| match c {
                'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', other => other
            }).collect()
        } else {
            reference
        };

        let stored: String = base
            .get_deleted_reference_bases()
            .iter()
            .map(|deleted| deleted.as_str().to_uppercase())
            .collect();
        assert_eq!(
            stored, expected,
            "wrong deleted bases at read position {}", base.get_read_position()
        );

        let insertion_point: u32 = base.get_deletion_read_position().unwrap();
        if *strand == Strand::Reverse {
            assert_eq!(insertion_point, base.get_read_position());
        } else {
            assert_eq!(insertion_point, base.get_read_position() + 1);
        }
    }

    lengths.sort();
    assert_eq!(lengths, vec![1, 1, 1, 1, 1, 1, 3], "expected the primary's six 1 bp deletions and its `-aag`");


    let read_name: &str = "scga-mini-rna-013-tumor_chunk_0000/182/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let model: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|r| Arc::new(r.clone())).collect::<Vec<_>>()
    );

    let mut checked_forward: usize = 0;
    let mut checked_reverse: usize = 0;
    let mut disagree_forward: usize = 0;
    let mut disagree_reverse: usize = 0;
    let mut first_disagreement: Option<String> = None;

    for base in model.get_bases().iter() {
        // Only bases the alignment asserts sit on one reference base: Match and Mismatch.
        // Insertion and Softclip carry an anchor, not a claim.
        let kind_ok: bool = matches!(
            base.get_kind(),
            AlignmentModelBaseKind::Match | AlignmentModelBaseKind::Mismatch
        );
        if !kind_ok {
            continue;
        }
        let (chromosome_id, position, strand) = match base.get_placement().get_coordinate() {
            Some(c) => c,
            None => continue
        };
        let chromosome: &Box<str> = chromosome_names_map.get_by_right(&chromosome_id).unwrap();
        let reference: Box<str> = get_fasta_sequence(&**chromosome, position, position, fasta_file);
        let reference_base: char = reference.to_uppercase().chars().next().unwrap();

        // The read base as the aligner saw it: forward record -> the read base itself,
        // reverse record -> its complement (the read runs antiparallel to the reference).
        let read_base: char = base.get_nucleotide().as_str().to_uppercase().chars().next().unwrap();
        let expected: char = if *strand == Strand::Forward {
            read_base
        } else {
            reverse_complement(&read_base.to_string()).to_uppercase().chars().next().unwrap()
        };

        let agrees: bool = expected == reference_base;
        let is_match: bool = *base.get_kind() == AlignmentModelBaseKind::Match;

        if *strand == Strand::Forward {
            checked_forward += 1;
            if is_match && !agrees { disagree_forward += 1; }
        } else {
            checked_reverse += 1;
            if is_match && !agrees { disagree_reverse += 1; }
        }

        if is_match && !agrees && first_disagreement.is_none() {
            first_disagreement = Some(format!(
                "read position {} nucleotide {} strand {:?} -> {}:{} reference {}",
                base.get_read_position(), read_base, strand, chromosome, position, reference_base
            ));
        }
    }

    // The overlap region 459..1563: kind was flipped to Softclip by apply_breakpoints, but the
    // placement left behind is the primary's per-base coordinate. Does it still name the right
    // reference base?
    let mut overlap_checked: usize = 0;
    let mut overlap_disagree: usize = 0;
    for base in model.get_bases().iter() {
        let p: u32 = base.get_read_position();
        if p < 459 || p > 1563 {
            continue;
        }
        if *base.get_kind() != AlignmentModelBaseKind::Softclip {
            continue;
        }
        let (chromosome_id, position, strand) = match base.get_placement().get_coordinate() {
            Some(c) => c,
            None => continue
        };
        let chromosome: &Box<str> = chromosome_names_map.get_by_right(&chromosome_id).unwrap();
        let reference: Box<str> = get_fasta_sequence(&**chromosome, position, position, fasta_file);
        let reference_base: char = reference.to_uppercase().chars().next().unwrap();
        let read_base: char = base.get_nucleotide().as_str().to_uppercase().chars().next().unwrap();
        let expected: char = if *strand == Strand::Forward {
            read_base
        } else {
            reverse_complement(&read_base.to_string()).to_uppercase().chars().next().unwrap()
        };
        overlap_checked += 1;
        if expected != reference_base {
            overlap_disagree += 1;
        }
    }

    assert!(checked_forward > 0, "no forward-placed bases - the mixed-strand precondition is gone");
    assert!(checked_reverse > 0, "no reverse-placed bases - the mixed-strand precondition is gone");
    assert_eq!(disagree_forward, 0, "forward (supplementary) placements disagree with the reference");
    assert_eq!(disagree_reverse, 0, "reverse (primary) placements disagree with the reference");
}

/// Deletion immediately followed by an insertion.
///
/// All three reads carry the rna-012 duplication, a ~540 bp block after chr17:7674180 that
/// minimap2 spelled as a short deletion butted straight against an `I` op. A deletion's
/// reference cursor advances past its own gap before the next cs token is read, so the inserted
/// bases are anchored on the far side of it. Taking the literally adjacent read position as the
/// event's right flank therefore lands on an inserted base one reference position short of the
/// true flank: read 188's 1 bp deletion collapses to a zero-length span - which is what
/// `is_repeat_variant` asserts against - and the 2 bp deletions of reads 39 and 106 report 1 bp.
#[test]
fn scga_mini_rna_012_transcript_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-012-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-012-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let (record_positions_map, read_names_map) = index_bam_records(bam_file, true, 2);

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    // Every read is on the reverse strand, so this also covers the flank ordering there.
    let read_names: [&str; 3] = [
        "scga-mini-rna-012-tumor_chunk_0000/39/ccs",
        "scga-mini-rna-012-tumor_chunk_0000/106/ccs",
        "scga-mini-rna-012-tumor_chunk_0000/188/ccs"
    ];

    let mut deletions: HashMap<&str, Vec<(u32, u32)>> = HashMap::new();
    for read_name in read_names.iter() {
        let read_id: usize = *read_names_map.get_by_left(*read_name).unwrap();
        let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
        let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

        let alignment_model: AlignmentModel = AlignmentModel::new(
            read_id,
            &*read_sequence,
            &quality_scores,
            &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
        );

        assert_eq!(alignment_model.get_records()[0].reference_strand, Strand::Reverse);

        for ((read_position_1, read_position_2), event) in alignment_model.get_events().iter() {
            if *event.get_kind() == AlignmentModelEventKind::Breakpoint {
                continue;
            }

            let coordinate_1 = alignment_model.get_bases()[*read_position_1 as usize]
                .get_placement()
                .get_coordinate()
                .expect("event flank is placed");
            let coordinate_2 = alignment_model.get_bases()[*read_position_2 as usize]
                .get_placement()
                .get_coordinate()
                .expect("event flank is placed");

            // Both flanks are bases the alignment placed on the reference, so they skip at
            // least one reference base between them. An inserted base as a flank breaks this.
            assert!(
                coordinate_1.1.abs_diff(coordinate_2.1) >= 2,
                "{read_name} event ({read_position_1},{read_position_2}) skips no reference: {} to {}",
                coordinate_1.1, coordinate_2.1
            );

            if *event.get_kind() == AlignmentModelEventKind::Deletion {
                deletions
                    .entry(*read_name)
                    .or_insert_with(Vec::new)
                    .push((coordinate_1.1.min(coordinate_2.1), coordinate_1.1.max(coordinate_2.1)));
            }
        }
    }

    // `~ct343ac:13-gt` + 545 bp and + 541 bp `I`: two deleted bases, chr17:7674194-7674195.
    assert!(
        deletions.get("scga-mini-rna-012-tumor_chunk_0000/39/ccs").unwrap().contains(&(7674193, 7674196)),
        "2 bp deletion at chr17:7674194-7674195 not reported over its true flanks"
    );
    assert!(
        deletions.get("scga-mini-rna-012-tumor_chunk_0000/106/ccs").unwrap().contains(&(7674193, 7674196)),
        "2 bp deletion at chr17:7674194-7674195 not reported over its true flanks"
    );

    // `~ct343ac:8-t` + 542 bp `I`: the deleted base is chr17:7674189, flanked by 7674188 and 7674190.
    assert!(
        deletions.get("scga-mini-rna-012-tumor_chunk_0000/188/ccs").unwrap().contains(&(7674188, 7674190)),
        "1 bp deletion at chr17:7674189 not reported over its true flanks"
    );
}

/// `min_terminal_soft_clip_ins_len` floors terminal soft-clip records by raw run length:
/// 1-2 bp terminal clips are aligner anchor jitter, and as variant sites they only fragment
/// clusters downstream. Two forward reads of rna-015 as minimap2 clipped them:
/// `scga-mini-rna-015-tumor_chunk_0000/367/ccs` (`225=1D138=1D86=2S`, cs `:225-t:138-t:86`, a
/// WSCD1 last-exon fragment) ends in a 2 bp clip, and `scga-mini-rna-015-tumor_chunk_0000/737/ccs`
/// (`1S198=5160N120=207N54=87N59=904N183=`, the whole of ACAP1 ENST00000575425.1 with no
/// difference) starts with a 1 bp clip. Neither read's cs has a `+` token, so the only
/// Insertion-typed records are the clip records. The floor must remove exactly the clip records
/// below it and leave every other record untouched; 0 must behave as 1.
#[test]
fn scga_mini_rna_015_clipfloor_transcript_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-015-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-015-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    // Each read's variant records, as sorted operation strings, at floors 0 through 3.
    // Breakpoint rescue is off, so its parameters are inert.
    let mut ops_by_floor: HashMap<(&str, u32), Vec<Box<str>>> = HashMap::new();
    for read_name in ["scga-mini-rna-015-tumor_chunk_0000/367/ccs", "scga-mini-rna-015-tumor_chunk_0000/737/ccs"] {
        let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
        let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
        let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());
        let alignment_model: AlignmentModel = AlignmentModel::new(
            read_id,
            &*read_sequence,
            &quality_scores,
            &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
        );
        let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
            &identify_transcript_model_exons(&alignment_model),
            &identify_transcript_model_splice_junctions(&alignment_model),
            &gene_annotator,
            &chromosome_names_map
        );
        let transcript_model: TranscriptModel = TranscriptModel::new(
            alignment_model,
            reference_transcript_matches,
            &gene_annotator,
            &chromosome_names_map,
            &fasta_map
        );
        for min_terminal_soft_clip_ins_len in 0..=3u32 {
            let mut ops: Vec<Box<str>> = RNAVariantRecordCaller::new(0, 0, min_terminal_soft_clip_ins_len, None)
                .call(&transcript_model)
                .iter()
                .map(|variant_record| variant_record.get_graph_operation_boxed_str())
                .collect();
            ops.sort_unstable();
            ops_by_floor.insert((read_name, min_terminal_soft_clip_ins_len), ops);
        }
    }

    // Trailing 2 bp clip (`AA`, after the last aligned base chr17:6124424) alongside the read's
    // two 1 bp deletions, chr17:6124199 and 6124338: present at floors 0-2, gone at 3, and the
    // deletions untouched throughout.
    let trailing_clip: Box<str> = "0:6124424:+:D:0:6124425:+:U:AA:2:INS".into();
    let deletion_1: Box<str> = "0:6124198:+:D:0:6124200:+:U::0:DEL".into();
    let deletion_2: Box<str> = "0:6124337:+:D:0:6124339:+:U::0:DEL".into();
    for floor in 0..=2u32 {
        assert_eq!(
            ops_by_floor[&("scga-mini-rna-015-tumor_chunk_0000/367/ccs", floor)],
            vec![deletion_1.clone(), deletion_2.clone(), trailing_clip.clone()],
            "floor {floor} must keep the 2 bp trailing clip"
        );
    }
    assert_eq!(
        ops_by_floor[&("scga-mini-rna-015-tumor_chunk_0000/367/ccs", 3)],
        vec![deletion_1.clone(), deletion_2.clone()],
        "floor 3 must remove exactly the 2 bp trailing clip"
    );

    // Leading 1 bp clip (`G`, before the first aligned base chr17:7336590): 0 and 1 are the same
    // floor; 2 removes it, leaving the read with no variant record at all.
    let leading_clip: Box<str> = "0:7336589:+:D:0:7336590:+:U:G:1:INS".into();
    assert_eq!(ops_by_floor[&("scga-mini-rna-015-tumor_chunk_0000/737/ccs", 0)], vec![leading_clip.clone()]);
    assert_eq!(
        ops_by_floor[&("scga-mini-rna-015-tumor_chunk_0000/737/ccs", 1)],
        ops_by_floor[&("scga-mini-rna-015-tumor_chunk_0000/737/ccs", 0)],
        "a floor of 1 must behave exactly like 0"
    );
    assert!(ops_by_floor[&("scga-mini-rna-015-tumor_chunk_0000/737/ccs", 2)].is_empty(), "a floor of 2 must remove the 1 bp clip");
    assert!(ops_by_floor[&("scga-mini-rna-015-tumor_chunk_0000/737/ccs", 3)].is_empty());
}

#[test]
fn identify_exons_keeps_an_insertion_beside_an_intron_or_a_deletion_in_the_exon_it_is_placed_in() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;

    // Exon 1 = chr1:101-200, intron = 201-300, exon 2 = 301-400. Inserted bases are placed at the
    // reference base to their left.
    let mut state: u32 = 2_463_534_242;
    let bases: String = (0..200)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            b"ACGT"[(state % 4) as usize] as char
        })
        .collect();
    let (exon_1, exon_2): (&str, &str) = (&bases[..100], &bases[100..]);
    let cases: Vec<(&str, &str, &str, Vec<(u32, u32, u32, u32)>)> = vec![
        ("no insertion", "100=100N100=", ":100~gt100ag:100", vec![(101, 200, 0, 99), (301, 400, 100, 199)]),
        ("insertion before the intron", "100=3I100N100=", ":100+ccc~gt100ag:100", vec![(101, 200, 0, 102), (301, 400, 103, 202)]),
        ("insertion after the intron", "100=100N3I100=", ":100~gt100ag+ccc:100", vec![(101, 200, 0, 99), (301, 400, 100, 202)]),
        ("insertion before a deletion", "100=3I10D100=", ":100+ccc-acgtacgtac:100", vec![(101, 310, 0, 202)]),
        ("insertion after a deletion", "100=10D3I100=", ":100-acgtacgtac+ccc:100", vec![(101, 310, 0, 202)])
    ];
    for (name, cigar, cs, expected) in cases {
        let read: String = if cigar.contains('I') { format!("{exon_1}CCC{exon_2}") } else { format!("{exon_1}{exon_2}") };
        let qualities: String = "I".repeat(read.len());
        let sam_text: String = format!(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nread-1\t0\tchr1\t101\t60\t{cigar}\t*\t0\t0\t{read}\t{qualities}\tcs:Z:{cs}\n"
        );
        let mut sam_reader = sam::io::Reader::new(sam_text.as_bytes());
        let header = sam_reader.read_header().unwrap();
        let mut writer = bam::io::Writer::new(Vec::new());
        for record_buf in sam_reader.record_bufs(&header) {
            writer.write_alignment_record(&header, &record_buf.unwrap()).unwrap();
        }
        writer.try_finish().unwrap();
        let mut reader = bam::io::Reader::new(writer.get_ref().get_ref().as_slice());
        let mut records: Vec<bam::Record> = Vec::new();
        let mut record = bam::Record::default();
        while reader.read_record(&mut record).unwrap() > 0 {
            records.push(record.clone());
        }
        let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
        let base_quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
        let bam_records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
        let alignment_model: AlignmentModel = AlignmentModel::new(1, &*read_sequence, &base_quality_scores, &bam_records);

        let exons: Vec<(u32, u32, u32, u32)> = identify_transcript_model_exons(&alignment_model)
            .iter()
            .map(|exon| (exon.reference_start, exon.reference_end, exon.read_start_position, exon.read_end_position))
            .collect();
        assert_eq!(exons, expected, "{name}");
    }
}

#[test]
fn contextualize_types_an_event_by_shared_gene_and_a_base_by_chromosome_and_strand() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;

    // chr1 and chr2 of 3,000 bases. Gene GA (+) has exons 101-200 and 901-1000; gene GX (-)
    // overlaps its first exon at 150-250 and gene GY (-) lies in its intron at 300-400.
    let mut state: u32 = 2_463_534_242;
    let mut random_bases = |n: usize| -> String {
        (0..n)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 17;
                state ^= state << 5;
                b"ACGT"[(state % 4) as usize] as char
            })
            .collect()
    };
    let temp_dir = tempfile::tempdir().unwrap();
    let fasta_file = temp_dir.path().join("genome.fa");
    fs::write(&fasta_file, format!(">chr1\n{}\n>chr2\n{}\n", random_bases(3_000), random_bases(3_000))).unwrap();
    let gtf_file = temp_dir.path().join("genes.gtf");
    let mut gtf: String = String::new();
    for (gene, strand, exons) in [("GA", "+", vec![(101, 200), (901, 1_000)]), ("GX", "-", vec![(150, 250)]), ("GY", "-", vec![(300, 400)])] {
        let (start, end): (u32, u32) = (exons.first().unwrap().0, exons.last().unwrap().1);
        let gene_attributes: String = format!("gene_id \"{gene}\"; gene_type \"protein_coding\"; gene_name \"{gene}\"; level 2;");
        let transcript_attributes: String = format!(
            "gene_id \"{gene}\"; transcript_id \"T{gene}\"; gene_type \"protein_coding\"; gene_name \"{gene}\"; transcript_type \"protein_coding\"; transcript_name \"T{gene}\"; level 2;"
        );
        gtf.push_str(&format!("chr1\tHAVANA\tgene\t{start}\t{end}\t.\t{strand}\t.\t{gene_attributes}\n"));
        gtf.push_str(&format!("chr1\tHAVANA\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t{transcript_attributes}\n"));
        for (number, (exon_start, exon_end)) in exons.iter().enumerate() {
            gtf.push_str(&format!(
                "chr1\tHAVANA\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t{transcript_attributes} exon_number {}; exon_id \"E{gene}{}\";\n",
                number + 1,
                number + 1
            ));
        }
    }
    fs::write(&gtf_file, gtf).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_file.to_str().unwrap(), "hg38", "v41");
    let fasta_map: FastaMap = FastaMap::new(fasta_file.to_str().unwrap());
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr1".into(), 0);
    chromosome_names_map.insert("chr2".into(), 1);
    let matches: Vec<ReferenceTranscriptMatch> = vec![ReferenceTranscriptMatch::new("GA", "GA", "TGA", &Vec::new(), 0, 0, 0, 0)];

    // Each read: its sequence and its records as SAM fields (flag, chromosome, position, CIGAR, cs).
    let first_part: String = random_bases(100);
    let fold_back: String = first_part.chars().rev().map(|base| match base { 'A' => 'T', 'C' => 'G', 'G' => 'C', _ => 'A' }).collect();
    let second_part: String = random_bases(100);
    let reads: Vec<(&str, String, Vec<(u16, &str, u32, &str, &str)>)> = vec![
        // A junction inside GA from 200 to 351: GX overlaps the donor and GY the acceptor.
        ("splice inside one gene", format!("{first_part}{second_part}"), vec![(0, "chr1", 101, "100=150N100=", ":100~gt150ag:100")]),
        // The second half folds back onto the other strand of 101-200.
        ("fold-back", format!("{first_part}{fold_back}"), vec![(0, "chr1", 101, "100=100S", ":100"), (2064, "chr1", 101, "100=100S", ":100")]),
        // The second half lies in the span of GA on the other strand, then on another chromosome.
        ("other strand", format!("{first_part}{second_part}"), vec![(0, "chr1", 101, "100=100S", ":100"), (2064, "chr1", 501, "100=100S", ":100")]),
        ("other chromosome", format!("{first_part}{second_part}"), vec![(0, "chr1", 101, "100=100S", ":100"), (2048, "chr2", 501, "100S100=", ":100")])
    ];
    let mut event_contexts: Vec<Option<AlignmentModelEventContext>> = Vec::new();
    let mut second_half_contexts: Vec<AlignmentModelBaseContext> = Vec::new();
    for (name, read, fields) in reads {
        let qualities: String = "I".repeat(read.len());
        let reverse: String = read.chars().rev().map(|base| match base { 'A' => 'T', 'C' => 'G', 'G' => 'C', _ => 'A' }).collect();
        let mut sam_text: String = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:3000\n@SQ\tSN:chr2\tLN:3000\n".to_string();
        for (flag, chromosome, position, cigar, cs) in fields {
            let sequence: &str = if flag & 16 == 16 { &reverse } else { &read };
            sam_text.push_str(&format!("read-1\t{flag}\t{chromosome}\t{position}\t60\t{cigar}\t*\t0\t0\t{sequence}\t{qualities}\tcs:Z:{cs}\n"));
        }
        let mut sam_reader = sam::io::Reader::new(sam_text.as_bytes());
        let header = sam_reader.read_header().unwrap();
        let mut writer = bam::io::Writer::new(Vec::new());
        for record_buf in sam_reader.record_bufs(&header) {
            writer.write_alignment_record(&header, &record_buf.unwrap()).unwrap();
        }
        writer.try_finish().unwrap();
        let mut reader = bam::io::Reader::new(writer.get_ref().get_ref().as_slice());
        let mut records: Vec<bam::Record> = Vec::new();
        let mut record = bam::Record::default();
        while reader.read_record(&mut record).unwrap() > 0 {
            records.push(record.clone());
        }
        let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
        let base_quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
        let bam_records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
        let alignment_model: AlignmentModel = AlignmentModel::new(1, &*read_sequence, &base_quality_scores, &bam_records);
        let annotation: TranscriptModelAnnotation = identify_transcript_model_annotation(
            &alignment_model,
            &matches,
            &gene_annotator,
            &chromosome_names_map,
            &fasta_map
        );
        event_contexts.push(annotation.get_event(99, 100).and_then(|event| event.get_context().clone()));
        second_half_contexts.push(annotation.get_base_context(150).unwrap_or_else(|| panic!("{name}: read position 150 is not annotated")).clone());
    }

    // Not a fusion: both sides lie in GA.
    assert_eq!(event_contexts[0], Some(AlignmentModelEventContext::NonCanonicalSplicing));
    // Not a back-splice: the two halves share positions but not the strand. Both lie in GA.
    assert_eq!(event_contexts[1], Some(AlignmentModelEventContext::NonCanonicalSplicing));
    // Not intronic: GA is on the other strand, or on the other chromosome.
    assert_eq!(second_half_contexts[2], AlignmentModelBaseContext::Intergenic);
    assert_eq!(second_half_contexts[3], AlignmentModelBaseContext::Intergenic);
}

/// A minus-strand gene with exons 1501-1700, 1001-1200 and 501-700, and a reverse read that
/// splices from 1501 to 700, skipping exon 2, with a 2-base deletion (1502-1503) one base after
/// the junction. Read base 60 (1501) ends both the deletion and the intron. Each skipped base goes
/// to the event on its side of its anchor: all of exon 2 to the intron, so it is one exon
/// truncation, and the deleted bases to the deletion, so none of them is.
#[test]
fn identify_transcript_model_annotation_gives_a_skipped_base_to_the_event_on_its_side_of_the_anchor() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;

    let mut state: u32 = 12_345;
    let mut random_bases = |n: usize| -> String {
        (0..n)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 17;
                state ^= state << 5;
                b"ACGT"[(state % 4) as usize] as char
            })
            .collect()
    };
    let temp_dir = tempfile::tempdir().unwrap();
    let fasta_file = temp_dir.path().join("genome.fa");
    fs::write(&fasta_file, format!(">chr1\n{}\n", random_bases(2_000))).unwrap();
    let gtf_file = temp_dir.path().join("genes.gtf");
    let gene_attributes: &str = "gene_id \"GN\"; gene_type \"protein_coding\"; gene_name \"GN\"; level 2;";
    let transcript_attributes: &str = "gene_id \"GN\"; transcript_id \"TGN\"; gene_type \"protein_coding\"; gene_name \"GN\"; transcript_type \"protein_coding\"; transcript_name \"TGN\"; level 2;";
    let mut gtf: String = format!(
        "chr1\tHAVANA\tgene\t501\t1700\t.\t-\t.\t{gene_attributes}\nchr1\tHAVANA\ttranscript\t501\t1700\t.\t-\t.\t{transcript_attributes}\n"
    );
    for (number, (exon_start, exon_end)) in [(1_501, 1_700), (1_001, 1_200), (501, 700)].iter().enumerate() {
        gtf.push_str(&format!(
            "chr1\tHAVANA\texon\t{exon_start}\t{exon_end}\t.\t-\t.\t{transcript_attributes} exon_number {}; exon_id \"EGN{}\";\n",
            number + 1,
            number + 1
        ));
    }
    fs::write(&gtf_file, gtf).unwrap();
    let gene_annotator = Gencode::new_with_defaults(gtf_file.to_str().unwrap(), "hg38", "v41");
    let fasta_map: FastaMap = FastaMap::new(fasta_file.to_str().unwrap());
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr1".into(), 0);

    // 651-700, intron 701-1500, 1501, deletion 1502-1503, 1504-1563.
    let read: String = random_bases(111);
    let qualities: String = "I".repeat(111);
    let sam_text: String = format!(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:2000\nread-1\t16\tchr1\t651\t60\t50=800N1=2D60=\t*\t0\t0\t{read}\t{qualities}\tcs:Z::50~ct800ac:1-ac:60\n"
    );
    let mut sam_reader = sam::io::Reader::new(sam_text.as_bytes());
    let header = sam_reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(Vec::new());
    for record_buf in sam_reader.record_bufs(&header) {
        writer.write_alignment_record(&header, &record_buf.unwrap()).unwrap();
    }
    writer.try_finish().unwrap();
    let mut reader = bam::io::Reader::new(writer.get_ref().get_ref().as_slice());
    let mut records: Vec<bam::Record> = Vec::new();
    let mut record = bam::Record::default();
    while reader.read_record(&mut record).unwrap() > 0 {
        records.push(record.clone());
    }
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    let base_quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
    let bam_records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
    let alignment_model: AlignmentModel = AlignmentModel::new(1, &*read_sequence, &base_quality_scores, &bam_records);

    let matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &identify_transcript_model_exons(&alignment_model),
        &identify_transcript_model_splice_junctions(&alignment_model),
        &gene_annotator,
        &chromosome_names_map
    );
    let transcript_model: TranscriptModel = TranscriptModel::new(alignment_model, matches, &gene_annotator, &chromosome_names_map, &fasta_map);

    // The read runs down the reference: read base 59 is 1504, 60 is 1501 and 61 is 700.
    let skipped_runs = |read_position_1: u32, read_position_2: u32| -> Vec<(u32, u32)> {
        transcript_model
            .get_annotation()
            .get_event(read_position_1, read_position_2)
            .map(|event| event.get_skipped_runs().iter().map(|run| (run.reference_start, run.reference_end)).collect())
            .unwrap_or_default()
    };
    assert_eq!(*transcript_model.get_alignment_model().get_event(60, 61).unwrap().get_kind(), AlignmentModelEventKind::Splicing);
    assert_eq!(*transcript_model.get_alignment_model().get_event(59, 60).unwrap().get_kind(), AlignmentModelEventKind::Deletion);
    assert_eq!(skipped_runs(60, 61), vec![(1_001, 1_200)]);
    assert_eq!(skipped_runs(59, 60), vec![(1_502, 1_503)]);

    let exon_truncations: Vec<(u32, u32)> = RNAVariantRecordCaller::new(0, 0, 1_000, None)
        .call(&transcript_model)
        .iter()
        .filter(|variant_record| *variant_record.get_variant_type() == VariantType::ExonTruncation)
        .map(|variant_record| (variant_record.get_position_1(), variant_record.get_position_2()))
        .collect();
    assert_eq!(exon_truncations, vec![(1_001, 1_200)]);
}
