use exacto_core::prelude::*;
use noodles_bam as bam;
use noodles_bam::bai;
use noodles_bam::bai::Index;
use noodles_sam::Header;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::str::FromStr;
use std::sync::Arc;

use super::*;


#[test]
fn scga_mini_dna_001_alignment_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());
    let chromosome_names: Vec<Box<str>> = get_chromosome_names(bam_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

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

    // The first read in the file that carries the SNV (an `X` at chr17:7674225) in its single record.
    let read_name: &str = "scga-mini-dna-001-tumor_scga-mini-dna-001-tumor-1_1/17/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());
    let alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );
    assert_eq!(alignment.get_records().len(), 1);
    assert_eq!(alignment.iter_base_quality_scores().len(), 19057);
    assert_eq!(alignment.get_read_id(), read_id);
    assert_eq!(alignment.get_read_sequence().len(), alignment.iter_base_quality_scores().len());
    assert_eq!(alignment.get_read_sequence().len(), alignment.iter_base_quality_scores().len());

    let mut found: bool = false;
    for base in alignment.get_bases() {
        if *base.get_kind() == AlignmentModelBaseKind::Mismatch {
            found = true;
        }
    }

    assert_eq!(found, true);
    assert_eq!(alignment.is_spliced(), false);

    let mut num_forward_mismatches: usize = 0;
    let mut num_reverse_mismatches: usize = 0;

    for (read_id, records) in records_map.iter() {
        let read_sequence: Box<str> = get_bam_fastx_read_sequence(records);
        let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records);
        let alignment_model: AlignmentModel = AlignmentModel::new(
            *read_id,
            &*read_sequence,
            &quality_scores,
            &records.iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
        );

        for base in alignment_model.get_bases().iter() {
            if *base.get_kind() != AlignmentModelBaseKind::Mismatch {
                // Only a mismatch names a reference base; a match already holds it in `nucleotide`.
                assert_eq!(
                    base.get_reference_nucleotide(),
                    None,
                    "read {} position {} is {:?} yet carries a reference allele",
                    read_id, base.get_read_position(), base.get_kind()
                );
                continue;
            }

            let (chromosome_id, position, strand) = base
                .get_placement()
                .get_coordinate()
                .expect("a mismatch is placed");
            let contig: &str = &chromosome_names[chromosome_id as usize];

            let mut expected: Nucleotide = Nucleotide::from_str(
                fasta_map
                    .get_sequence(contig, position as usize, position as usize)
                    .to_uppercase()
                    .as_str()
            ).unwrap();
            if *strand == Strand::Reverse {
                expected = expected.complement();
                num_reverse_mismatches += 1;
            } else {
                num_forward_mismatches += 1;
            }

            assert_eq!(
                base.get_reference_nucleotide(),
                Some(&expected),
                "read {} position {} carries the wrong reference allele",
                read_id, base.get_read_position()
            );
            // A mismatch that agrees with the read is a complement applied on the wrong side.
            assert_ne!(
                base.get_reference_nucleotide().unwrap().as_str().to_uppercase(),
                base.get_nucleotide().as_str().to_uppercase(),
                "read {} position {} has reference == read but is kind Mismatch",
                read_id, base.get_read_position()
            );
        }
    }

    // Both arms of the strand branch must actually execute or the complement assertion is vacuous.
    assert!(num_forward_mismatches > 0, "no forward-strand mismatch exercised");
    assert!(num_reverse_mismatches > 0, "no reverse-strand mismatch exercised");
}

#[test]
fn scga_mini_dna_002_alignment_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-002-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-002-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

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

    // The first read in the file that carries the 12-base insertion after chr17:7674224.
    let read_name: &str = "scga-mini-dna-002-tumor_scga-mini-dna-002-tumor-1_1/9/ccs";
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

    let mut found: bool = false;
    for base in alignment.get_bases() {
        if *base.get_kind() == AlignmentModelBaseKind::Insertion {
            found = true;
        }
    }

    assert_eq!(found, true);
    assert_eq!(alignment.is_spliced(), false);
}

#[test]
fn scga_mini_dna_003_alignment_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-003-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-003-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

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

    // The first read in the file that carries the 30-base deletion of chr17:7674201-7674230.
    let read_name: &str = "scga-mini-dna-003-tumor_scga-mini-dna-003-tumor-1_1/1/ccs";
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

    let mut found: bool = false;
    for event in alignment.get_events().values() {
        if event.get_kind() == &AlignmentModelEventKind::Deletion {
            found = true;
        }
    }

    assert_eq!(found, true);
    assert_eq!(alignment.is_spliced(), false);
}

#[test]
fn scga_mini_dna_013_alignment_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-013-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-013-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

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

    // The one read across the whole inverted duplication, from the other strand: chr17:7687109 down
    // to 7679901, the inverted copy forward from 7673301 to 7679895, then 7679900 down to 7678934.
    let read_name: &str = "scga-mini-dna-013-tumor_scga-mini-dna-013-tumor-1_1/14/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    assert_eq!(alignment.get_records().len(), 3);
    assert_eq!(alignment.get_records().get(0).unwrap().reference_strand, Strand::Reverse);
    assert_eq!(alignment.get_records().get(1).unwrap().reference_strand, Strand::Forward);
    assert_eq!(alignment.get_records().get(2).unwrap().reference_strand, Strand::Reverse);

    let mut num_breakpoint: usize = 0;
    for event in alignment.get_events().values() {
        if event.get_kind() == &AlignmentModelEventKind::Breakpoint {
            num_breakpoint += 1;
        }
    }

    assert_eq!(num_breakpoint, 2);
    assert_eq!(alignment.is_spliced(), false);
}

#[test]
fn scga_mini_dna_013_multi_record_base_placements() {
    use noodles_sam::alignment::Record as _;
    use noodles_sam::alignment::record::cigar::op::Kind;

    let bam_path =
        Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-013-tumor_minimap2_sorted.bam");
    let read_name = "scga-mini-dna-013-tumor_scga-mini-dna-013-tumor-1_1/14/ccs";

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_path)
        .unwrap();
    reader.read_header().unwrap();

    let mut records: Vec<bam::Record> = Vec::new();
    for result in reader.records() {
        let record = result.unwrap();
        if record.flags().is_unmapped() || record.flags().is_secondary() {
            continue;
        }
        let name = std::str::from_utf8(record.name().unwrap().as_ref()).unwrap();
        if name == read_name {
            records.push(record);
        }
    }

    assert_eq!(records.len(), 3, "fixture must contain all three alignments");
    records.sort_by_key(|record| get_alignment_start_position(record));

    let chromosome_id =
        records[0].reference_sequence_id().unwrap().unwrap() as u16;
    let sequence = get_bam_fastx_read_sequence(&records);
    let qualities = get_bam_fastx_base_quality_scores(&records);
    let records: Vec<Arc<bam::Record>> =
        records.into_iter().map(Arc::new).collect();

    let model = AlignmentModel::new(0, &sequence, &qualities, &records);

    assert_eq!(model.get_bases().len(), 14_569);

    let spans: Vec<_> = model
        .get_records()
        .iter()
        .map(|record| (
            record.read_start,
            record.read_end,
            record.reference_strand.clone(),
        ))
        .collect();

    assert_eq!(
        spans,
        vec![
            (0, 7_100, Strand::Reverse),
            (7_101, 13_606, Strand::Forward),
            (13_611, 14_568, Strand::Reverse),
        ]
    );

    // Read positions are zero-based; reference positions are one-based.
    // No read base lies in two record spans, so each record's CIGAR, walked in reference
    // order, names every base it covers: `=` a match, `X` a mismatch, `I` an insertion on
    // the reference base before it. A reverse record walks the read down from its read end.
    for record in model.get_records() {
        let strand = &record.reference_strand;
        let mut reference_position = get_alignment_start_position(&record.record) - 1;
        let mut offset: u32 = 0;
        for op in record.record.cigar().iter().map(|op| op.unwrap()) {
            let kind = match op.kind() {
                Kind::SequenceMatch => AlignmentModelBaseKind::Match,
                Kind::SequenceMismatch => AlignmentModelBaseKind::Mismatch,
                Kind::Insertion => AlignmentModelBaseKind::Insertion,
                Kind::Deletion => {
                    reference_position += op.len() as u32;
                    continue;
                },
                Kind::SoftClip => continue,
                other => panic!("unexpected CIGAR operation {other:?}")
            };
            for _ in 0..op.len() {
                if kind != AlignmentModelBaseKind::Insertion {
                    reference_position += 1;
                }
                let read_position = if *strand == Strand::Forward {
                    record.read_start + offset
                } else {
                    record.read_end - offset
                };
                offset += 1;
                let base = model.get_base(read_position);

                assert_eq!(
                    base.get_kind(),
                    &kind,
                    "unexpected kind at read position {read_position}"
                );
                assert_eq!(
                    base.get_placement().get_coordinate(),
                    Some((chromosome_id, reference_position, strand)),
                    "unexpected placement at read position {read_position}"
                );
                assert!(base.is_aligned());
            }
        }
        assert_eq!(offset, record.read_end - record.read_start + 1);
    }

    // The record ends land where the records say: the reverse arm from 7,687,109 down to
    // 7,679,901, the inverted copy from 7,673,301 to 7,679,895, the reverse arm on from
    // 7,679,900 down to 7,678,934.
    for (read_position, position, strand) in [
        (0u32, 7_687_109u32, Strand::Reverse),
        (7_100, 7_679_901, Strand::Reverse),
        (7_101, 7_673_301, Strand::Forward),
        (13_606, 7_679_895, Strand::Forward),
        (13_611, 7_679_900, Strand::Reverse),
        (14_568, 7_678_934, Strand::Reverse),
    ] {
        assert_eq!(
            model.get_base(read_position).get_placement().get_coordinate(),
            Some((chromosome_id, position, &strand)),
            "unexpected placement at read position {read_position}"
        );
    }

    // The gap between the second and third records is soft-clipped
    // and has no reference placement.
    for read_position in 13_607..13_611 {
        let base = model.get_base(read_position);
        assert_eq!(base.get_kind(), &AlignmentModelBaseKind::Softclip);
        assert!(!base.get_placement().is_placed());
        assert!(!base.is_aligned());
    }

    // Regression version of the audit's base-kind census.
    assert!(model.get_bases().iter().all(|base| {
        base.get_kind() != &AlignmentModelBaseKind::Unaligned
    }));

    // Supplementary records must place bases outside the primary's span.
    let primary = model
        .get_records()
        .iter()
        .find(|record| !record.record.flags().is_supplementary())
        .unwrap();

    assert_eq!((primary.read_start, primary.read_end), (0, 7_100));
    assert!(model.get_base(7_101).is_aligned());
    assert!(model.get_base(14_568).is_aligned());
}

#[test]
fn scga_mini_dna_012_multi_record_placements_are_independent_of_input_order() {
    use noodles_sam::alignment::Record as _;

    // This read's three records share read positions twice: the primary and the first
    // supplementary at 9244, the two supplementaries at 10462-10676. Which record places
    // them depends on the order the records are applied in.
    let bam_path =
        Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-012-tumor_minimap2_sorted.bam");
    let read_name = "scga-mini-dna-012-tumor_scga-mini-dna-012-tumor-1_1/7/ccs";

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_path)
        .unwrap();
    reader.read_header().unwrap();

    let mut records: Vec<bam::Record> = Vec::new();
    for result in reader.records() {
        let record = result.unwrap();
        if record.flags().is_unmapped() || record.flags().is_secondary() {
            continue;
        }
        let name = std::str::from_utf8(record.name().unwrap().as_ref()).unwrap();
        if name == read_name {
            records.push(record);
        }
    }

    assert_eq!(records.len(), 3, "fixture must contain all three alignments");
    records.sort_by_key(|record| get_alignment_start_position(record));

    let sequence = get_bam_fastx_read_sequence(&records);
    let qualities = get_bam_fastx_base_quality_scores(&records);
    let records: Vec<Arc<bam::Record>> =
        records.into_iter().map(Arc::new).collect();

    let original = AlignmentModel::new(0, &sequence, &qualities, &records);

    let reversed_records: Vec<Arc<bam::Record>> =
        records.iter().rev().cloned().collect();
    let reversed =
        AlignmentModel::new(0, &sequence, &qualities, &reversed_records);

    assert_eq!(original.get_bases().len(), reversed.get_bases().len());

    for (position, (a, b)) in original
        .get_bases()
        .iter()
        .zip(reversed.get_bases())
        .enumerate()
    {
        assert_eq!(
            a.get_kind(),
            b.get_kind(),
            "record order changed base kind at read position {position}"
        );
        assert_eq!(
            a.get_placement(),
            b.get_placement(),
            "record order changed placement at read position {position}"
        );
    }
}

#[test]
fn scga_mini_dna_013_alignment_model_is_input_order_insensitive() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-013-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-013-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

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

    let mut num_multi_record_reads: usize = 0;
    for (read_id, records) in records_map.iter() {
        if records.len() < 2 {
            continue;
        }
        num_multi_record_reads += 1;

        let read_sequence: Box<str> = get_bam_fastx_read_sequence(records);
        let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records);
        let forward: AlignmentModel = AlignmentModel::new(
            *read_id,
            &*read_sequence,
            &quality_scores,
            &records.iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
        );
        let reversed: AlignmentModel = AlignmentModel::new(
            *read_id,
            &*read_sequence,
            &quality_scores,
            &records.iter().rev().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
        );

        // Records come out sorted by read start whatever order they went in.
        assert_eq!(forward.get_records().len(), reversed.get_records().len(), "read {read_id}");
        for (a, b) in forward.get_records().iter().zip(reversed.get_records().iter()) {
            assert_eq!(
                (a.read_start, a.read_end, &a.reference_strand),
                (b.read_start, b.read_end, &b.reference_strand),
                "read {read_id}"
            );
        }

        // Every base carries the same kind and the same placement.
        assert_eq!(forward.get_bases().len(), reversed.get_bases().len(), "read {read_id}");
        for (a, b) in forward.get_bases().iter().zip(reversed.get_bases().iter()) {
            assert_eq!(a.get_kind(), b.get_kind(), "read {} position {}", read_id, a.get_read_position());
            assert_eq!(a.get_placement(), b.get_placement(), "read {} position {}", read_id, a.get_read_position());
        }

        // Same events, keyed by the same read positions.
        let mut forward_events: Vec<((u32, u32), AlignmentModelEventKind)> = forward
            .get_events()
            .iter()
            .map(|(key, event)| (*key, event.get_kind().clone()))
            .collect();
        let mut reversed_events: Vec<((u32, u32), AlignmentModelEventKind)> = reversed
            .get_events()
            .iter()
            .map(|(key, event)| (*key, event.get_kind().clone()))
            .collect();
        forward_events.sort_by_key(|(key, _kind)| *key);
        reversed_events.sort_by_key(|(key, _kind)| *key);
        assert_eq!(forward_events, reversed_events, "read {read_id}");
    }

    // scga-mini-dna-013 carries nine split reads, eight of two records and one of three;
    // without them the comparison above is vacuous.
    assert_eq!(num_multi_record_reads, 9);
}

#[test]
fn scga_mini_dna_012_alignment_model_places_contested_bases_by_last_record() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-012-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-012-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_names: Vec<Box<str>> = get_chromosome_names(bam_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

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

    // A read across the tandem duplication in three forward records: the copy ending at chr17:7679900,
    // then the copy starting at 7673301 split in two. The primary and the first supplementary share one
    // base, and the two supplementaries share 215.
    let read_name: &str = "scga-mini-dna-012-tumor_scga-mini-dna-012-tumor-1_1/7/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    // Sorted by read start: the primary, then the two supplementaries.
    let spans: Vec<(u32, u32, Strand)> = alignment
        .get_records()
        .iter()
        .map(|record| (record.read_start, record.read_end, record.reference_strand.clone()))
        .collect();
    assert_eq!(
        spans,
        vec![
            (0, 9244, Strand::Forward),
            (9244, 10676, Strand::Forward),
            (10462, 15494, Strand::Forward)
        ]
    );
    assert!(!alignment.get_records()[0].record.flags().is_supplementary());
    assert!(alignment.get_records()[1].record.flags().is_supplementary());
    assert!(alignment.get_records()[2].record.flags().is_supplementary());
    assert_eq!(get_alignment_start_position(&alignment.get_records()[0].record), 7670657);
    assert_eq!(get_alignment_end_position(&alignment.get_records()[0].record), 7679901);
    assert_eq!(get_alignment_start_position(&alignment.get_records()[1].record), 7673300);
    assert_eq!(get_alignment_end_position(&alignment.get_records()[1].record), 7674732);
    assert_eq!(get_alignment_start_position(&alignment.get_records()[2].record), 7674519);
    assert_eq!(get_alignment_end_position(&alignment.get_records()[2].record), 7679553);

    // 216 read positions lie in two record spans: the seam at 9244, and 10462-10676.
    let mut contested: Vec<u32> = Vec::new();
    for base in alignment.get_bases().iter() {
        let p: u32 = base.get_read_position();
        let num_claimants: usize = alignment
            .get_records()
            .iter()
            .filter(|record| p >= record.read_start && p <= record.read_end)
            .count();
        if num_claimants > 1 {
            contested.push(p);
        }
    }
    assert_eq!(contested, [9244].into_iter().chain(10462..=10676).collect::<Vec<u32>>());

    // The later record in read-start order owns each shared base, placed from its own alignment start.
    // All of them are flagged soft-clips so no record can read them as aligned bases.
    for (p, expected_position) in [(9244u32, 7673300u32), (10462, 7674519), (10676, 7674732)] {
        let (chromosome_id, position, strand) = alignment.get_base(p).get_placement().get_coordinate().unwrap();
        assert_eq!(&*chromosome_names[chromosome_id as usize], "chr17");
        assert_eq!(position, expected_position);
        assert_eq!(*strand, Strand::Forward);
    }
    assert!(contested.iter().all(|p| *alignment.get_base(*p).get_kind() == AlignmentModelBaseKind::Softclip));

    // The earlier record keeps the base before each shared run, and the later one resumes after it.
    for (p, expected_position) in [(9243u32, 7679900u32), (9245, 7673301), (10461, 7674518), (10677, 7674733)] {
        let (_chromosome_id, position, strand) = alignment.get_base(p).get_placement().get_coordinate().unwrap();
        assert_eq!(position, expected_position);
        assert_eq!(*strand, Strand::Forward);
        assert_eq!(*alignment.get_base(p).get_kind(), AlignmentModelBaseKind::Match);
    }

    // Kind census, and the reach beyond the primary: both supplementary records are aligned but for
    // the shared bases.
    let mut num_match: usize = 0;
    let mut num_mismatch: usize = 0;
    let mut num_insertion: usize = 0;
    let mut num_softclip: usize = 0;
    let mut num_softclip_placed: usize = 0;
    let mut num_aligned_outside_primary: usize = 0;
    for base in alignment.get_bases().iter() {
        match base.get_kind() {
            AlignmentModelBaseKind::Match => num_match += 1,
            AlignmentModelBaseKind::Mismatch => num_mismatch += 1,
            AlignmentModelBaseKind::Insertion => num_insertion += 1,
            AlignmentModelBaseKind::Softclip => {
                num_softclip += 1;
                if base.get_placement().is_placed() {
                    num_softclip_placed += 1;
                }
            },
            other => panic!("unexpected kind {:?} at position {}", other, base.get_read_position())
        }
        let p: u32 = base.get_read_position();
        if p > 9244 && base.is_aligned() {
            num_aligned_outside_primary += 1;
        }
    }
    assert_eq!(num_match, 15270);
    assert_eq!(num_mismatch, 0);
    assert_eq!(num_insertion, 9);
    assert_eq!(num_softclip, 216);
    assert_eq!(num_softclip_placed, 216);
    assert_eq!(num_aligned_outside_primary, (15494 - 9244) - 215);

    // One breakpoint flanks each shared run. Every other event is one of the 13 deletions the three
    // CIGARs hold.
    let mut num_deletion: usize = 0;
    for event in alignment.get_events().values() {
        if *event.get_kind() == AlignmentModelEventKind::Deletion {
            num_deletion += 1;
        }
    }
    assert_eq!(alignment.get_events().len(), 15);
    assert_eq!(num_deletion, 13);
    assert_eq!(*alignment.get_event(9243, 9245).unwrap().get_kind(), AlignmentModelEventKind::Breakpoint);
    assert_eq!(*alignment.get_event(10461, 10677).unwrap().get_kind(), AlignmentModelEventKind::Breakpoint);
}

#[test]
fn scga_mini_dna_013_alignment_model_walks_reference_in_strand_direction() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-013-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-013-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

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

    let mut num_records_checked: usize = 0;
    let mut num_reverse_records_checked: usize = 0;
    for (read_id, records) in records_map.iter() {
        let read_sequence: Box<str> = get_bam_fastx_read_sequence(records);
        let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records);
        let alignment: AlignmentModel = AlignmentModel::new(
            *read_id,
            &*read_sequence,
            &quality_scores,
            &records.iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
        );

        for (index, alignment_record) in alignment.get_records().iter().enumerate() {
            let reference_start: u32 = get_alignment_start_position(&alignment_record.record);
            let reference_end: u32 = get_alignment_end_position(&alignment_record.record);
            let chromosome_id: u16 = alignment_record.record.reference_sequence_id().unwrap().unwrap() as u16;
            let strand: &Strand = &alignment_record.reference_strand;

            let mut previous_position: Option<u32> = None;
            for p in alignment_record.read_start..=alignment_record.read_end {
                // A read position inside several record spans belongs to the last of them in read-start order.
                let owner: usize = alignment
                    .get_records()
                    .iter()
                    .enumerate()
                    .filter(|(_i, record)| p >= record.read_start && p <= record.read_end)
                    .map(|(i, _record)| i)
                    .max()
                    .unwrap();
                if owner != index {
                    continue;
                }

                let (placed_chromosome_id, position, placed_strand) = alignment
                    .get_base(p)
                    .get_placement()
                    .get_coordinate()
                    .unwrap_or_else(|| panic!("read {read_id} position {p} inside record {index} is unplaced"));
                assert_eq!(placed_chromosome_id, chromosome_id, "read {read_id} position {p}");
                assert_eq!(placed_strand, strand, "read {read_id} position {p}");
                assert!(
                    position >= reference_start && position <= reference_end,
                    "read {read_id} position {p} placed at {position}, outside [{reference_start}, {reference_end}]"
                );

                // Walking the read walks the reference monotonically, in the strand's direction.
                if let Some(previous) = previous_position {
                    match strand {
                        Strand::Forward => assert!(
                            position >= previous,
                            "read {read_id} position {p}: {position} < {previous} on the forward strand"
                        ),
                        _ => assert!(
                            position <= previous,
                            "read {read_id} position {p}: {position} > {previous} on the reverse strand"
                        )
                    }
                }
                previous_position = Some(position);

                // A record's own read ends, when it keeps them, land on its reference ends.
                if p == alignment_record.read_start {
                    let expected: u32 = match strand {
                        Strand::Forward => reference_start,
                        _ => reference_end
                    };
                    assert_eq!(position, expected, "read {read_id} record {index} at its read start");
                }
                if p == alignment_record.read_end {
                    let expected: u32 = match strand {
                        Strand::Forward => reference_end,
                        _ => reference_start
                    };
                    assert_eq!(position, expected, "read {read_id} record {index} at its read end");
                }
            }

            num_records_checked += 1;
            if *strand == Strand::Reverse {
                num_reverse_records_checked += 1;
            }
        }
    }

    // 64 single-record reads, 8 split reads of two records and 1 of three.
    assert_eq!(num_records_checked, 83);
    assert!(num_reverse_records_checked > 0, "no reverse-strand record exercised");
}

#[test]
fn scga_mini_dna_013_alignment_model_leaves_no_unaligned_base() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-013-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-013-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

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

    let mut num_reads: usize = 0;
    let mut num_multi_record_reads: usize = 0;
    for (read_id, records) in records_map.iter() {
        let read_sequence: Box<str> = get_bam_fastx_read_sequence(records);
        let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records);
        let alignment: AlignmentModel = AlignmentModel::new(
            *read_id,
            &*read_sequence,
            &quality_scores,
            &records.iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
        );
        num_reads += 1;
        if records.len() > 1 {
            num_multi_record_reads += 1;
        }

        for base in alignment.get_bases().iter() {
            let p: u32 = base.get_read_position();

            // Every base starts out Unaligned; construction must resolve all of them.
            assert_ne!(*base.get_kind(), AlignmentModelBaseKind::Unaligned, "read {read_id} position {p}");

            // Aligned kinds are always placed; a soft-clip never counts as aligned.
            match base.get_kind() {
                AlignmentModelBaseKind::Match
                | AlignmentModelBaseKind::Mismatch
                | AlignmentModelBaseKind::Insertion => {
                    assert!(base.get_placement().is_placed(), "read {read_id} position {p} is {:?} but unplaced", base.get_kind());
                    assert!(base.is_aligned(), "read {read_id} position {p}");
                },
                AlignmentModelBaseKind::Softclip => {
                    assert!(!base.is_aligned(), "read {read_id} position {p}");
                },
                AlignmentModelBaseKind::Unaligned => unreachable!()
            }

            // Only a soft-clip may go unplaced.
            if !base.get_placement().is_placed() {
                assert_eq!(*base.get_kind(), AlignmentModelBaseKind::Softclip, "read {read_id} position {p}");
            }

            // A base outside every record span is a soft-clip.
            let covered: bool = alignment
                .get_records()
                .iter()
                .any(|record| p >= record.read_start && p <= record.read_end);
            if !covered {
                assert_eq!(*base.get_kind(), AlignmentModelBaseKind::Softclip, "read {read_id} position {p}");
            }
        }
    }

    assert_eq!(num_reads, 73);
    assert_eq!(num_multi_record_reads, 9);
}

#[test]
fn scga_mini_dna_016_alignment_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-016-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-016-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

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

    // The first read in the file with both records forward and its junction on the true breakends:
    // chr17:3491600 joined to chr17:6085001 across the 12 inserted bases.
    let read_name: &str = "scga-mini-dna-016-tumor_scga-mini-dna-016-tumor-3_1/29/ccs";
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

    let mut num_breakpoint: usize = 0;
    for event in alignment.get_events().values() {
        if event.get_kind() == &AlignmentModelEventKind::Breakpoint {
            num_breakpoint += 1;
        }
    }

    for ((r1,r2), event) in alignment.get_events().iter() {
        if *event.get_kind() == AlignmentModelEventKind::Breakpoint {
            let b1 = alignment.get_base(*r1);
            let b2 = alignment.get_base(*r2);
            assert!(b1.get_placement().get_coordinate().unwrap().1 == 3491600);
            assert!(b2.get_placement().get_coordinate().unwrap().1 == 6085001);
        }
    }

    assert_eq!(num_breakpoint, 1);
    assert_eq!(alignment.is_spliced(), false);
}

#[test]
fn scga_mini_dna_007_alignment_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-007-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

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

    // The first read in the file with both records forward: chr17:3491600 joined to chr17:6085001.
    let read_name: &str = "scga-mini-dna-007-tumor_scga-mini-dna-007-tumor-3_1/97/ccs";
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

    let mut num_breakpoints: usize = 0;
    for event in alignment.get_events().values() {
        if event.get_kind() == &AlignmentModelEventKind::Breakpoint {
            num_breakpoints += 1;
        }
    }

    assert_eq!(num_breakpoints, 1);
    assert_eq!(alignment.is_spliced(), false);
}

#[test]
fn scga_mini_dna_012_alignment_model_places_contested_base_on_last_record_locus() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-012-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-012-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_names: Vec<Box<str>> = get_chromosome_names(bam_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

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

    // The first read in the file across the tandem duplication with both records forward and no
    // terminal clip: the copy ending at chr17:7679900, then the copy starting at 7673301.
    let read_name: &str = "scga-mini-dna-012-tumor_scga-mini-dna-012-tumor-1_1/5/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    // The primary followed by a supplementary arm 6.6 kb upstream, sharing read position 9164: the
    // read base there matches both chr17:7679901 and chr17:7673300.
    let spans: Vec<(u32, u32, Strand)> = alignment
        .get_records()
        .iter()
        .map(|record| (record.read_start, record.read_end, record.reference_strand.clone()))
        .collect();
    assert_eq!(spans, vec![(0, 9164, Strand::Forward), (9164, 14520, Strand::Forward)]);
    assert!(!alignment.get_records()[0].record.flags().is_supplementary());
    assert!(alignment.get_records()[1].record.flags().is_supplementary());
    assert_eq!(get_alignment_start_position(&alignment.get_records()[0].record), 7670737);
    assert_eq!(get_alignment_end_position(&alignment.get_records()[0].record), 7679901);
    assert_eq!(get_alignment_start_position(&alignment.get_records()[1].record), 7673300);
    assert_eq!(get_alignment_end_position(&alignment.get_records()[1].record), 7678658);

    let mut contested: Vec<u32> = Vec::new();
    for base in alignment.get_bases().iter() {
        let p: u32 = base.get_read_position();
        let num_claimants: usize = alignment
            .get_records()
            .iter()
            .filter(|record| p >= record.read_start && p <= record.read_end)
            .count();
        if num_claimants > 1 {
            contested.push(p);
        }
    }
    assert_eq!(contested, vec![9164]);

    // The seam goes to the last record: its alignment start, not the first record's end, flagged a soft-clip.
    let (chromosome_id, position, strand) = alignment.get_base(9164).get_placement().get_coordinate().unwrap();
    assert_eq!(&*chromosome_names[chromosome_id as usize], "chr17");
    assert_eq!(position, 7673300);
    assert_eq!(*strand, Strand::Forward);
    assert_eq!(*alignment.get_base(9164).get_kind(), AlignmentModelBaseKind::Softclip);

    // The base before the seam is the first record's, one short of its alignment end.
    let (chromosome_id, position, strand) = alignment.get_base(9163).get_placement().get_coordinate().unwrap();
    assert_eq!(&*chromosome_names[chromosome_id as usize], "chr17");
    assert_eq!(position, 7679900);
    assert_eq!(*strand, Strand::Forward);
    assert_eq!(*alignment.get_base(9163).get_kind(), AlignmentModelBaseKind::Match);

    // The base after the seam continues the supplementary.
    let (chromosome_id, position, _strand) = alignment.get_base(9165).get_placement().get_coordinate().unwrap();
    assert_eq!(&*chromosome_names[chromosome_id as usize], "chr17");
    assert_eq!(position, 7673301);
    assert_eq!(*alignment.get_base(9165).get_kind(), AlignmentModelBaseKind::Match);

    // Every other base is aligned (14,504 matches, 1 mismatch, 15 inserted bases in the two CIGARs);
    // the seam is the only soft-clip, and it stays placed.
    let num_softclip: usize = alignment
        .get_bases()
        .iter()
        .filter(|base| *base.get_kind() == AlignmentModelBaseKind::Softclip)
        .count();
    let num_match: usize = alignment
        .get_bases()
        .iter()
        .filter(|base| *base.get_kind() == AlignmentModelBaseKind::Match)
        .count();
    let num_aligned: usize = alignment
        .get_bases()
        .iter()
        .filter(|base| base.is_aligned())
        .count();
    assert_eq!(num_softclip, 1);
    assert_eq!(num_match, 14504);
    assert_eq!(num_aligned, 14520);
    assert!(alignment.get_bases().iter().all(|base| base.get_placement().is_placed()));

    // The single breakpoint flanks the seam; the other events are the 15 deletions in the two CIGARs.
    let num_breakpoints: usize = alignment
        .get_events()
        .values()
        .filter(|event| *event.get_kind() == AlignmentModelEventKind::Breakpoint)
        .count();
    assert_eq!(alignment.get_events().len(), 16);
    assert_eq!(num_breakpoints, 1);
    assert_eq!(*alignment.get_event(9163, 9165).unwrap().get_kind(), AlignmentModelEventKind::Breakpoint);
}

#[test]
fn scga_mini_dna_006_alignment_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-006-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-006-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

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

    // The first read in the file across the inversion forward then reverse: chr17:7687771-7701200
    // forward, then chr17:7717000 down to 7713906.
    let read_name: &str = "scga-mini-dna-006-tumor_scga-mini-dna-006-tumor-1_1/89/ccs";
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
    assert_eq!(alignment.get_records().get(1).unwrap().reference_strand, Strand::Reverse);

    let read_name: &str = "scga-mini-dna-006-tumor_scga-mini-dna-006-tumor-1_1/89/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    let mut num_breakpoints: usize = 0;
    for event in alignment.get_events().values() {
        if event.get_kind() == &AlignmentModelEventKind::Breakpoint {
            num_breakpoints += 1;
        }
    }

    assert_eq!(num_breakpoints, 1);
    assert_eq!(alignment.is_spliced(), false);

    let mut num_breakpoints: usize = 0;
    for event in alignment.get_events().values() {
        if event.get_kind() == &AlignmentModelEventKind::Breakpoint {
            num_breakpoints += 1;
        }
    }

    assert_eq!(num_breakpoints, 1);
    assert_eq!(alignment.is_spliced(), false);
}

#[test]
fn scga_mini_dna_002_alignment_model_returns_terminal_softclips() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-002-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-002-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_names: Vec<Box<str>> = get_chromosome_names(bam_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

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

    // Every read here is a single record; nine carry a terminal soft-clip. From their CIGARs: the read
    // length, the record's read span, its reference span, the clip runs in read orientation, and the
    // number of `=` bases. They cover both ends of the read on both strands.
    let expected: HashMap<&str, (u32, (u32, u32), (u32, u32), Vec<(u32, u32)>, usize)> = HashMap::from([
        ("scga-mini-dna-002-tumor_scga-mini-dna-002-tumor-1_1/7/ccs", (18989, (2, 18988), (7668422, 7687490), vec![(0, 1)], 18917)),
        ("scga-mini-dna-002-tumor_scga-mini-dna-002-tumor-1_1/16/ccs", (16290, (0, 16286), (7668460, 7684927), vec![(16287, 16289)], 16195)),
        ("scga-mini-dna-002-tumor_scga-mini-dna-002-tumor-2_1/8/ccs", (15784, (4, 15783), (7668545, 7684305), vec![(0, 3)], 15743)),
        ("scga-mini-dna-002-tumor_scga-mini-dna-002-tumor-2_1/10/ccs", (17495, (8, 17494), (7669588, 7687044), vec![(0, 7)], 17432)),
        ("scga-mini-dna-002-tumor_scga-mini-dna-002-tumor-1_1/15/ccs", (16798, (2, 16797), (7669709, 7686506), vec![(0, 1)], 16769)),
        ("scga-mini-dna-002-tumor_scga-mini-dna-002-tumor-2_1/2/ccs", (15715, (7, 15714), (7669734, 7685429), vec![(0, 6)], 15672)),
        ("scga-mini-dna-002-tumor_scga-mini-dna-002-tumor-1_1/19/ccs", (16268, (0, 16260), (7670194, 7686462), vec![(16261, 16267)], 16223)),
        ("scga-mini-dna-002-tumor_scga-mini-dna-002-tumor-1_1/24/ccs", (14628, (0, 14624), (7671489, 7686116), vec![(14625, 14627)], 14584)),
        ("scga-mini-dna-002-tumor_scga-mini-dna-002-tumor-2_1/22/ccs", (14414, (4, 14413), (7672217, 7686608), vec![(0, 3)], 14376))
    ]);
    let mut num_reads: usize = 0;
    let mut num_clipped_forward: usize = 0;
    let mut num_clipped_reverse: usize = 0;
    for (read_id, records) in records_map.iter() {
        assert_eq!(records.len(), 1, "read {read_id}");
        num_reads += 1;

        let read_sequence: Box<str> = get_bam_fastx_read_sequence(records);
        let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records);
        let alignment: AlignmentModel = AlignmentModel::new(
            *read_id,
            &*read_sequence,
            &quality_scores,
            &records.iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
        );
        assert_eq!(alignment.get_records().len(), 1, "read {read_id}");
        assert_eq!(alignment.is_spliced(), false, "read {read_id}");

        // Maximal runs of soft-clip bases; no base may be left Unaligned.
        let mut clip_runs: Vec<(u32, u32)> = Vec::new();
        for base in alignment.get_bases().iter() {
            assert_ne!(*base.get_kind(), AlignmentModelBaseKind::Unaligned, "read {} position {}", read_id, base.get_read_position());
            if *base.get_kind() == AlignmentModelBaseKind::Softclip {
                let p: u32 = base.get_read_position();
                match clip_runs.last_mut() {
                    Some(last) if last.1 + 1 == p => last.1 = p,
                    _ => clip_runs.push((p, p))
                }
            }
        }
        if clip_runs.is_empty() {
            continue;
        }

        let read_name: &str = read_names_map.get_by_right(read_id).unwrap();
        let (num_bases, read_span, reference_span, expected_clip_runs, expected_num_match) = expected
            .get(read_name)
            .unwrap_or_else(|| panic!("read {read_name} is clipped but not expected to be"));
        let alignment_record: &AlignmentRecord = &alignment.get_records()[0];
        let strand: &Strand = &alignment_record.reference_strand;
        let reference_start: u32 = get_alignment_start_position(&alignment_record.record);
        let reference_end: u32 = get_alignment_end_position(&alignment_record.record);
        assert_eq!((alignment_record.read_start, alignment_record.read_end), *read_span, "read {read_name}");
        assert_eq!((reference_start, reference_end), *reference_span, "read {read_name}");
        assert_eq!(alignment.num_bases(), *num_bases, "read {read_name}");
        assert_eq!(clip_runs, *expected_clip_runs, "read {read_name}");

        // A terminal clip is pinned to one anchor: the base just before the alignment on the read's 5'
        // side, the alignment's last base on its 3' side. On the reverse strand the two anchors swap.
        let (expected_left_anchor, expected_right_anchor): (u32, u32) = match strand {
            Strand::Forward => (reference_start - 1, reference_end),
            _ => (reference_end, reference_start - 1)
        };
        for (start, stop) in clip_runs.iter().copied() {
            let expected_anchor: u32 = if start == 0 { expected_left_anchor } else { expected_right_anchor };
            for p in start..=stop {
                let base: &AlignmentModelBase = alignment.get_base(p);
                assert_eq!(*base.get_kind(), AlignmentModelBaseKind::Softclip, "read {read_name} position {p}");
                assert_eq!(base.is_aligned(), false, "read {read_name} position {p}");
                let (chromosome_id, position, placed_strand) = base
                    .get_placement()
                    .get_coordinate()
                    .unwrap_or_else(|| panic!("read {read_name} position {p}: terminal clip base is unplaced"));
                assert_eq!(&*chromosome_names[chromosome_id as usize], "chr17", "read {read_name} position {p}");
                assert_eq!(placed_strand, strand, "read {read_name} position {p}");
                assert_eq!(position, expected_anchor, "read {read_name} position {p}");
            }
        }

        // The aligned body sits between the clips and runs the full reference span.
        let (first_aligned_position, last_aligned_position): (u32, u32) = match strand {
            Strand::Forward => (reference_start, reference_end),
            _ => (reference_end, reference_start)
        };
        let (read_start, read_end): (u32, u32) = *read_span;
        assert_eq!(*alignment.get_base(read_start).get_kind(), AlignmentModelBaseKind::Match, "read {read_name}");
        assert_eq!(alignment.get_base(read_start).get_placement().get_coordinate().unwrap().1, first_aligned_position, "read {read_name}");
        assert_eq!(*alignment.get_base(read_end).get_kind(), AlignmentModelBaseKind::Match, "read {read_name}");
        assert_eq!(alignment.get_base(read_end).get_placement().get_coordinate().unwrap().1, last_aligned_position, "read {read_name}");
        assert!((read_start..=read_end).all(|p| alignment.get_base(p).is_aligned()), "read {read_name}");
        let num_match: usize = alignment
            .get_bases()
            .iter()
            .filter(|base| *base.get_kind() == AlignmentModelBaseKind::Match)
            .count();
        assert_eq!(num_match, *expected_num_match, "read {read_name}");
        assert!(alignment.get_bases().iter().all(|base| base.get_placement().is_placed()), "read {read_name}");

        // A terminal clip raises no event: the read's only events are the deletions in its CIGAR.
        assert!(
            alignment.get_events().values().all(|event| *event.get_kind() == AlignmentModelEventKind::Deletion),
            "read {read_name}"
        );

        if *strand == Strand::Reverse {
            num_clipped_reverse += 1;
        } else {
            num_clipped_forward += 1;
        }
    }

    assert_eq!(num_reads, 66);
    assert_eq!(num_clipped_forward + num_clipped_reverse, expected.len());
    // Both strand branches of the anchor rule must actually run.
    assert!(num_clipped_forward > 0, "no forward-strand clipped read exercised");
    assert!(num_clipped_reverse > 0, "no reverse-strand clipped read exercised");
}

#[test]
fn scga_mini_rna_013_alignment_model_returns_matches() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-013-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-013-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();

    let chromosome_names: Vec<Box<str>> = get_chromosome_names(bam_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

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

    // The one read in the file with two records.
    let read_name: &str = "scga-mini-rna-013-tumor_chunk_0000/182/ccs";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    // A spliced forward supplementary followed by a reverse primary; they share read positions 459..=1563.
    assert_eq!(alignment.num_bases(), 1900);
    assert_eq!(alignment.is_spliced(), true);
    let spans: Vec<(u32, u32, Strand)> = alignment
        .get_records()
        .iter()
        .map(|record| (record.read_start, record.read_end, record.reference_strand.clone()))
        .collect();
    assert_eq!(spans, vec![(0, 1563, Strand::Forward), (459, 1893, Strand::Reverse)]);
    assert!(alignment.get_records()[0].record.flags().is_supplementary());
    assert!(!alignment.get_records()[1].record.flags().is_supplementary());
    assert_eq!(get_alignment_start_position(&alignment.get_records()[0].record), 7673786);
    assert_eq!(get_alignment_end_position(&alignment.get_records()[0].record), 7822054);
    assert_eq!(get_alignment_start_position(&alignment.get_records()[1].record), 7668429);
    assert_eq!(get_alignment_end_position(&alignment.get_records()[1].record), 7673608);

    // The whole 1,105-base overlap belongs to the last record: reverse strand, walking down from its
    // alignment end, and flagged soft-clip so neither record can read it as aligned sequence. The
    // primary's one inserted base (read position 1037) sits on the reference position of 1038.
    let mut contested: Vec<u32> = Vec::new();
    for base in alignment.get_bases().iter() {
        let p: u32 = base.get_read_position();
        let num_claimants: usize = alignment
            .get_records()
            .iter()
            .filter(|record| p >= record.read_start && p <= record.read_end)
            .count();
        if num_claimants > 1 {
            contested.push(p);
        }
    }
    assert_eq!(contested, (459..=1563).collect::<Vec<u32>>());
    let mut previous_position: Option<u32> = None;
    let mut num_repeated_positions: usize = 0;
    for p in 459..=1563u32 {
        let base: &AlignmentModelBase = alignment.get_base(p);
        assert_eq!(*base.get_kind(), AlignmentModelBaseKind::Softclip, "position {p}");
        assert_eq!(base.is_aligned(), false, "position {p}");
        let (chromosome_id, position, strand) = base.get_placement().get_coordinate().unwrap();
        assert_eq!(&*chromosome_names[chromosome_id as usize], "chr17", "position {p}");
        assert_eq!(*strand, Strand::Reverse, "position {p}");
        if let Some(previous) = previous_position {
            assert!(position <= previous, "position {p}: {position} > {previous} on the reverse strand");
            if position == previous {
                num_repeated_positions += 1;
            }
        }
        previous_position = Some(position);
    }
    assert_eq!(num_repeated_positions, 1);
    assert_eq!(alignment.get_base(1037).get_placement(), alignment.get_base(1038).get_placement());
    assert_eq!(alignment.get_base(459).get_placement().get_coordinate().unwrap().1, 7673608);
    assert_eq!(alignment.get_base(1563).get_placement().get_coordinate().unwrap().1, 7668761);

    // Either side of the seam: the forward record's last kept base and the reverse record's first aligned one.
    let (_chromosome_id, position, strand) = alignment.get_base(458).get_placement().get_coordinate().unwrap();
    assert_eq!(*alignment.get_base(458).get_kind(), AlignmentModelBaseKind::Insertion);
    assert_eq!(position, 7675059);
    assert_eq!(*strand, Strand::Forward);
    let (_chromosome_id, position, strand) = alignment.get_base(1564).get_placement().get_coordinate().unwrap();
    assert_eq!(*alignment.get_base(1564).get_kind(), AlignmentModelBaseKind::Match);
    assert_eq!(position, 7668760);
    assert_eq!(*strand, Strand::Reverse);

    // The primary's 6-base terminal clip hangs off the base before its alignment start.
    for p in 1894..=1899u32 {
        let base: &AlignmentModelBase = alignment.get_base(p);
        assert_eq!(*base.get_kind(), AlignmentModelBaseKind::Softclip, "position {p}");
        assert_eq!(base.get_placement().get_coordinate().unwrap().1, 7668428, "position {p}");
        assert_eq!(*base.get_placement().get_coordinate().unwrap().2, Strand::Reverse, "position {p}");
    }

    // Kind census: no base is left unplaced, none is Unaligned. Every mismatch either record holds
    // lies inside the overlap, so none survives as a mismatch.
    let mut num_match: usize = 0;
    let mut num_insertion: usize = 0;
    let mut num_softclip: usize = 0;
    for base in alignment.get_bases().iter() {
        assert!(base.get_placement().is_placed(), "position {} is unplaced", base.get_read_position());
        match base.get_kind() {
            AlignmentModelBaseKind::Match => num_match += 1,
            AlignmentModelBaseKind::Insertion => num_insertion += 1,
            AlignmentModelBaseKind::Softclip => num_softclip += 1,
            other => panic!("unexpected kind {:?} at position {}", other, base.get_read_position())
        }
    }
    assert_eq!(num_match, 612);
    assert_eq!(num_insertion, 177);
    assert_eq!(num_softclip, 1105 + 6);

    // Every insertion base sits on the reference position of the base before it, so a run of inserted
    // bases shares one anchor. One run survives: the 1,180-base insertion is cut to 177 by the overlap.
    let mut insertion_runs: Vec<(u32, u32, u32)> = Vec::new();
    for base in alignment.get_bases().iter() {
        if *base.get_kind() != AlignmentModelBaseKind::Insertion {
            continue;
        }
        let p: u32 = base.get_read_position();
        let (_chromosome_id, position, _strand) = base.get_placement().get_coordinate().unwrap();
        let (_chromosome_id, anchor, _strand) = alignment.get_base(p - 1).get_placement().get_coordinate().unwrap();
        assert_eq!(position, anchor, "insertion at position {p} is not anchored to the base before it");
        match insertion_runs.last_mut() {
            Some(last) if last.1 + 1 == p && last.2 == position => last.1 = p,
            _ => insertion_runs.push((p, p, position))
        }
    }
    assert_eq!(insertion_runs, vec![(282, 458, 7675059)]);

    // Events: three splices from the forward record, one breakpoint at the seam, two splices and seven
    // deletions from the reverse record. Nothing from the forward record survives inside the overlap it lost.
    let mut events: Vec<((u32, u32), AlignmentModelEventKind)> = alignment
        .get_events()
        .iter()
        .map(|(key, event)| (*key, event.get_kind().clone()))
        .collect();
    events.sort_by_key(|(key, _kind)| *key);
    assert_eq!(
        events,
        vec![
            ((51, 52), AlignmentModelEventKind::Splicing),
            ((159, 160), AlignmentModelEventKind::Splicing),
            ((274, 275), AlignmentModelEventKind::Splicing),
            ((458, 1564), AlignmentModelEventKind::Breakpoint),
            ((532, 533), AlignmentModelEventKind::Splicing),
            ((639, 640), AlignmentModelEventKind::Splicing),
            ((954, 955), AlignmentModelEventKind::Deletion),
            ((1085, 1086), AlignmentModelEventKind::Deletion),
            ((1252, 1253), AlignmentModelEventKind::Deletion),
            ((1472, 1473), AlignmentModelEventKind::Deletion),
            ((1515, 1516), AlignmentModelEventKind::Deletion),
            ((1656, 1657), AlignmentModelEventKind::Deletion),
            ((1850, 1851), AlignmentModelEventKind::Deletion)
        ]
    );

    // On the reverse strand the deleted reference bases hang off the higher read position of each pair.
    for (lower, higher, length) in [
        (954u32, 955u32, 1usize),
        (1085, 1086, 1),
        (1252, 1253, 1),
        (1472, 1473, 3),
        (1515, 1516, 1),
        (1656, 1657, 1),
        (1850, 1851, 1)
    ] {
        assert_eq!(alignment.get_base(higher).get_deleted_reference_bases().len(), length, "position {higher}");
        assert_eq!(alignment.get_base(lower).get_deleted_reference_bases().len(), 0, "position {lower}");
    }
}

#[test]
fn scga_mini_rna_001_alignment_model_returns_terminal_softclips() {
    use noodles_sam::alignment::Record as _;

    let bam_path = Path::new(env!("EXACTO_TEST_DATA"))
        .join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(&bam_full_path)
        .unwrap();
    reader.read_header().unwrap();

    let mut records_by_name: HashMap<String, Vec<bam::Record>> = HashMap::new();

    for result in reader.records() {
        let record = result.unwrap();
        if record.flags().is_unmapped() || record.flags().is_secondary() {
            continue;
        }

        let name = std::str::from_utf8(record.name().unwrap().as_ref())
            .unwrap()
            .to_string();

        records_by_name.entry(name).or_default().push(record);
    }

    // Read lengths, expected inclusive clip ranges in original read coordinates
    // (zero-based), and their reference anchors (one-based). The first reverse
    // reads in the file with a 1-base clip at the read start (`...28=1S`,
    // anchored at the alignment end 7,676,548) and a 2-base clip at the read end
    // (`2S37=...`, anchored before the alignment start 7,668,422).
    let cases: [(&str, u32, u32, u32, u32); 2] = [
        ("scga-mini-rna-001-tumor_chunk_0000/72/ccs", 2_320, 0, 0, 7_676_548),
        ("scga-mini-rna-001-tumor_chunk_0001/5/ccs", 13_543, 13_541, 13_542, 7_668_421),
    ];

    for (read_id, (read_name, num_bases, clip_start, clip_end, anchor)) in
        cases.into_iter().enumerate()
    {
        let records = records_by_name
            .get(read_name)
            .expect("fixture read must exist");

        assert_eq!(records.len(), 1, "{read_name}: expected one record");
        assert!(!records[0].flags().is_supplementary());
        assert!(records[0].flags().is_reverse_complemented());

        let chromosome_id =
            records[0].reference_sequence_id().unwrap().unwrap() as u16;
        let read_sequence = get_bam_fastx_read_sequence(records);
        let quality_scores = get_bam_fastx_base_quality_scores(records);
        let records: Vec<Arc<bam::Record>> = records
            .iter()
            .map(|record| Arc::new(record.clone()))
            .collect();

        let alignment = AlignmentModel::new(
            read_id,
            &read_sequence,
            &quality_scores,
            &records,
        );

        assert_eq!(alignment.get_records().len(), 1);
        assert_eq!(alignment.num_bases(), num_bases);

        let mut num_softclips: usize = 0;

        for base in alignment.get_bases() {
            let position = base.get_read_position();

            assert_ne!(
                base.get_kind(),
                &AlignmentModelBaseKind::Unaligned,
                "{read_name}: unclassified base at read position {position}"
            );

            if (clip_start..=clip_end).contains(&position) {
                num_softclips += 1;

                assert_eq!(
                    base.get_kind(),
                    &AlignmentModelBaseKind::Softclip,
                    "{read_name}: expected softclip at read position {position}"
                );

                // Also fails if the base has no reference placement.
                assert_eq!(
                    base.get_placement().get_coordinate(),
                    Some((chromosome_id, anchor, &Strand::Reverse)),
                    "{read_name}: wrong clip anchor at read position {position}"
                );

                assert!(!base.is_aligned());
            } else {
                assert_ne!(
                    base.get_kind(),
                    &AlignmentModelBaseKind::Softclip,
                    "{read_name}: unexpected softclip at read position {position}"
                );
            }
        }

        assert_eq!(
            num_softclips,
            (clip_end - clip_start + 1) as usize,
            "{read_name}: incorrect terminal clip count"
        );
    }
}


/// A base below the minimum base quality is removed from a run of mismatches, and the base
/// that remains keeps its own position.
///
///   Reference position   111   112
///   Reference base       G     G
///   Read base            T     A
///   Base quality         10    40
///
/// With a minimum base quality of 20 the read supports G>A at 112 and nothing at 111.
#[test]
fn identify_variant_records_keeps_position_when_first_mismatch_has_low_base_quality() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;

    let sam_text: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nread-1\t0\tchr1\t101\t60\t10=2X8=\t*\t0\t0\tACGTACGTACTAACGTACGT\tIIIIIIIIII+IIIIIIIII\tcs:Z::10*gt*ga:8\n";
    let mut sam_reader = sam::io::Reader::new(sam_text.as_bytes());
    let header = sam_reader.read_header().unwrap();
    let record_buf = sam_reader.record_bufs(&header).next().unwrap().unwrap();
    let mut writer = bam::io::Writer::new(Vec::new());
    writer.write_alignment_record(&header, &record_buf).unwrap();
    writer.try_finish().unwrap();
    let mut reader = bam::io::Reader::new(writer.get_ref().get_ref().as_slice());
    let mut record = bam::Record::default();
    reader.read_record(&mut record).unwrap();
    let records: Vec<bam::Record> = vec![record];
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    let base_quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
    let bam_records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
    let alignment_model: AlignmentModel = AlignmentModel::new(1, &*read_sequence, &base_quality_scores, &bam_records);

    let variant_records: Vec<VariantRecord> = DNAVariantRecordCaller::new(4, 20, 0).call(&alignment_model);

    assert_eq!(variant_records.len(), 1);
    assert_eq!(*variant_records[0].get_variant_type(), VariantType::SingleNucleotideVariant);
    assert_eq!(variant_records[0].get_position_1(), 111);
    assert_eq!(variant_records[0].get_position_2(), 113);
    assert_eq!(&*variant_records[0].get_standardized_sequence(), "A");
    assert_eq!(variant_records[0].get_read_position_1(), 11);
    assert_eq!(variant_records[0].get_read_position_2(), 11);
}


/// The same read as above with the low base quality on the second mismatch.
///
///   Reference position   111   112
///   Read base            T     A
///   Base quality         40    10
///
/// The read supports G>T at 111 and nothing at 112.
#[test]
fn identify_variant_records_keeps_position_when_second_mismatch_has_low_base_quality() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;

    let sam_text: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nread-1\t0\tchr1\t101\t60\t10=2X8=\t*\t0\t0\tACGTACGTACTAACGTACGT\tIIIIIIIIIII+IIIIIIII\tcs:Z::10*gt*ga:8\n";
    let mut sam_reader = sam::io::Reader::new(sam_text.as_bytes());
    let header = sam_reader.read_header().unwrap();
    let record_buf = sam_reader.record_bufs(&header).next().unwrap().unwrap();
    let mut writer = bam::io::Writer::new(Vec::new());
    writer.write_alignment_record(&header, &record_buf).unwrap();
    writer.try_finish().unwrap();
    let mut reader = bam::io::Reader::new(writer.get_ref().get_ref().as_slice());
    let mut record = bam::Record::default();
    reader.read_record(&mut record).unwrap();
    let records: Vec<bam::Record> = vec![record];
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    let base_quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
    let bam_records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
    let alignment_model: AlignmentModel = AlignmentModel::new(1, &*read_sequence, &base_quality_scores, &bam_records);

    let variant_records: Vec<VariantRecord> = DNAVariantRecordCaller::new(4, 20, 0).call(&alignment_model);

    assert_eq!(variant_records.len(), 1);
    assert_eq!(*variant_records[0].get_variant_type(), VariantType::SingleNucleotideVariant);
    assert_eq!(variant_records[0].get_position_1(), 110);
    assert_eq!(variant_records[0].get_position_2(), 112);
    assert_eq!(&*variant_records[0].get_standardized_sequence(), "T");
    assert_eq!(variant_records[0].get_read_position_1(), 10);
    assert_eq!(variant_records[0].get_read_position_2(), 10);
}


/// A low base quality in the middle of a run of three mismatches cuts the run in two.
///
///   Reference position   111   112   113
///   Reference base       G     G     T
///   Read base            T     A     C
///   Base quality         40    10    40
///
/// The read supports G>T at 111 and T>C at 113, and nothing at 112.
#[test]
fn identify_variant_records_cuts_run_of_mismatches_at_low_base_quality() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;

    let sam_text: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nread-1\t0\tchr1\t101\t60\t10=3X7=\t*\t0\t0\tACGTACGTACTACCGTACGT\tIIIIIIIIIII+IIIIIIII\tcs:Z::10*gt*ga*tc:7\n";
    let mut sam_reader = sam::io::Reader::new(sam_text.as_bytes());
    let header = sam_reader.read_header().unwrap();
    let record_buf = sam_reader.record_bufs(&header).next().unwrap().unwrap();
    let mut writer = bam::io::Writer::new(Vec::new());
    writer.write_alignment_record(&header, &record_buf).unwrap();
    writer.try_finish().unwrap();
    let mut reader = bam::io::Reader::new(writer.get_ref().get_ref().as_slice());
    let mut record = bam::Record::default();
    reader.read_record(&mut record).unwrap();
    let records: Vec<bam::Record> = vec![record];
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    let base_quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
    let bam_records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
    let alignment_model: AlignmentModel = AlignmentModel::new(1, &*read_sequence, &base_quality_scores, &bam_records);

    let variant_records: Vec<VariantRecord> = DNAVariantRecordCaller::new(4, 20, 0).call(&alignment_model);

    assert_eq!(variant_records.len(), 2);
    assert_eq!(*variant_records[0].get_variant_type(), VariantType::SingleNucleotideVariant);
    assert_eq!(variant_records[0].get_position_1(), 110);
    assert_eq!(variant_records[0].get_position_2(), 112);
    assert_eq!(&*variant_records[0].get_standardized_sequence(), "T");
    assert_eq!(*variant_records[1].get_variant_type(), VariantType::SingleNucleotideVariant);
    assert_eq!(variant_records[1].get_position_1(), 112);
    assert_eq!(variant_records[1].get_position_2(), 114);
    assert_eq!(&*variant_records[1].get_standardized_sequence(), "C");
}


/// The first test on the reverse strand. The read runs down the reference, so the base at 112
/// comes before the base at 111 in the read.
///
///   Reference position   111   112
///   Reference base       G     G
///   Read base (forward)  T     A
///   Base quality         10    40
///
/// The read supports G>A at 112 and nothing at 111.
#[test]
fn identify_variant_records_keeps_position_when_first_mismatch_has_low_base_quality_on_reverse_strand() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;

    let sam_text: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nread-1\t16\tchr1\t101\t60\t10=2X8=\t*\t0\t0\tACGTACGTACTAACGTACGT\tIIIIIIIIII+IIIIIIIII\tcs:Z::10*gt*ga:8\n";
    let mut sam_reader = sam::io::Reader::new(sam_text.as_bytes());
    let header = sam_reader.read_header().unwrap();
    let record_buf = sam_reader.record_bufs(&header).next().unwrap().unwrap();
    let mut writer = bam::io::Writer::new(Vec::new());
    writer.write_alignment_record(&header, &record_buf).unwrap();
    writer.try_finish().unwrap();
    let mut reader = bam::io::Reader::new(writer.get_ref().get_ref().as_slice());
    let mut record = bam::Record::default();
    reader.read_record(&mut record).unwrap();
    let records: Vec<bam::Record> = vec![record];
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    let base_quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
    let bam_records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
    let alignment_model: AlignmentModel = AlignmentModel::new(1, &*read_sequence, &base_quality_scores, &bam_records);

    let variant_records: Vec<VariantRecord> = DNAVariantRecordCaller::new(4, 20, 0).call(&alignment_model);

    assert_eq!(variant_records.len(), 1);
    assert_eq!(*variant_records[0].get_variant_type(), VariantType::SingleNucleotideVariant);
    assert_eq!(*variant_records[0].get_strand_1(), Strand::Reverse);
    assert_eq!(variant_records[0].get_position_1(), 111);
    assert_eq!(variant_records[0].get_position_2(), 113);
    assert_eq!(&*variant_records[0].get_standardized_sequence(), "A");
    assert_eq!(variant_records[0].get_read_position_1(), 8);
    assert_eq!(variant_records[0].get_read_position_2(), 8);
}


/// An insertion of 12 bases after chr1:110, four of them below the minimum base quality (Q10
/// against 20). Most of its bases pass, so the record keeps all 12: the base quality filter keeps
/// or drops a sequence whole and never spells one the read did not hold.
#[test]
fn identify_variant_records_keeps_an_insertion_whole_when_most_of_its_bases_pass() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;

    let sam_text: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nread-1\t0\tchr1\t101\t60\t10=12I10=\t*\t0\t0\tACGTACGTACACGTTGCAGTCATTGACCATGG\tIIIIIIIIIII+II+II+II+IIIIIIIIIII\tcs:Z::10+acgttgcagtca:10\n";
    let mut sam_reader = sam::io::Reader::new(sam_text.as_bytes());
    let header = sam_reader.read_header().unwrap();
    let record_buf = sam_reader.record_bufs(&header).next().unwrap().unwrap();
    let mut writer = bam::io::Writer::new(Vec::new());
    writer.write_alignment_record(&header, &record_buf).unwrap();
    writer.try_finish().unwrap();
    let mut reader = bam::io::Reader::new(writer.get_ref().get_ref().as_slice());
    let mut record = bam::Record::default();
    reader.read_record(&mut record).unwrap();
    let records: Vec<bam::Record> = vec![record];
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    let base_quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
    let bam_records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
    let alignment_model: AlignmentModel = AlignmentModel::new(1, &*read_sequence, &base_quality_scores, &bam_records);

    let variant_records: Vec<VariantRecord> = DNAVariantRecordCaller::new(4, 20, 0).call(&alignment_model);

    assert_eq!(variant_records.len(), 1);
    assert_eq!(*variant_records[0].get_variant_type(), VariantType::Insertion);
    assert_eq!((variant_records[0].get_position_1(), variant_records[0].get_position_2()), (110, 111));
    assert_eq!(&*variant_records[0].get_standardized_sequence(), "ACGTTGCAGTCA");
}


/// The same insertion with seven of its 12 bases below the minimum base quality: most of its
/// bases fail, so the read gives no record rather than an insertion of the five that pass.
#[test]
fn identify_variant_records_drops_an_insertion_when_most_of_its_bases_fail() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;

    let sam_text: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nread-1\t0\tchr1\t101\t60\t10=12I10=\t*\t0\t0\tACGTACGTACACGTTGCAGTCATTGACCATGG\tIIIIIIIIII++I+I++I+II+IIIIIIIIII\tcs:Z::10+acgttgcagtca:10\n";
    let mut sam_reader = sam::io::Reader::new(sam_text.as_bytes());
    let header = sam_reader.read_header().unwrap();
    let record_buf = sam_reader.record_bufs(&header).next().unwrap().unwrap();
    let mut writer = bam::io::Writer::new(Vec::new());
    writer.write_alignment_record(&header, &record_buf).unwrap();
    writer.try_finish().unwrap();
    let mut reader = bam::io::Reader::new(writer.get_ref().get_ref().as_slice());
    let mut record = bam::Record::default();
    reader.read_record(&mut record).unwrap();
    let records: Vec<bam::Record> = vec![record];
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    let base_quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
    let bam_records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
    let alignment_model: AlignmentModel = AlignmentModel::new(1, &*read_sequence, &base_quality_scores, &bam_records);

    let variant_records: Vec<VariantRecord> = DNAVariantRecordCaller::new(4, 20, 0).call(&alignment_model);

    assert!(variant_records.is_empty());
}


/// A leading clip of 10 bases, three of them below the minimum base quality. The clip is kept
/// whole as the breakend's sequence.
#[test]
fn identify_variant_records_keeps_a_terminal_clip_whole_when_most_of_its_bases_pass() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;

    let sam_text: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nread-1\t0\tchr1\t101\t60\t10S20=\t*\t0\t0\tGGGTTTAAACACGTACGTACGTACGTACGT\tI+II+II+IIIIIIIIIIIIIIIIIIIIII\tcs:Z::20\n";
    let mut sam_reader = sam::io::Reader::new(sam_text.as_bytes());
    let header = sam_reader.read_header().unwrap();
    let record_buf = sam_reader.record_bufs(&header).next().unwrap().unwrap();
    let mut writer = bam::io::Writer::new(Vec::new());
    writer.write_alignment_record(&header, &record_buf).unwrap();
    writer.try_finish().unwrap();
    let mut reader = bam::io::Reader::new(writer.get_ref().get_ref().as_slice());
    let mut record = bam::Record::default();
    reader.read_record(&mut record).unwrap();
    let records: Vec<bam::Record> = vec![record];
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    let base_quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
    let bam_records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
    let alignment_model: AlignmentModel = AlignmentModel::new(1, &*read_sequence, &base_quality_scores, &bam_records);

    let variant_records: Vec<VariantRecord> = DNAVariantRecordCaller::new(4, 20, 4).call(&alignment_model);

    assert_eq!(variant_records.len(), 1);
    assert_eq!(*variant_records[0].get_variant_type(), VariantType::Breakpoint);
    assert_eq!(variant_records[0].get_position_1(), 101);
    assert_eq!(&*variant_records[0].get_sequence(), "GGGTTTAAAC");
}


/// minimap2 writes a match as `=ACGT` with --cs=long and as `:4` with --cs. The two forms of one
/// alignment give the same variant records.
#[test]
fn identify_variant_records_reads_the_long_form_of_the_cs_tag() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;

    let mut variant_records_by_form: Vec<Vec<Box<str>>> = Vec::new();
    for cs_tag in [":10*gt*ga:8", "=ACGTACGTAC*gt*ga=ACGTACGT"] {
        let sam_text: String = format!("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nread-1\t0\tchr1\t101\t60\t10=2X8=\t*\t0\t0\tACGTACGTACTAACGTACGT\tIIIIIIIIIIIIIIIIIIII\tcs:Z:{cs_tag}\n");
        let mut sam_reader = sam::io::Reader::new(sam_text.as_bytes());
        let header = sam_reader.read_header().unwrap();
        let record_buf = sam_reader.record_bufs(&header).next().unwrap().unwrap();
        let mut writer = bam::io::Writer::new(Vec::new());
        writer.write_alignment_record(&header, &record_buf).unwrap();
        writer.try_finish().unwrap();
        let mut reader = bam::io::Reader::new(writer.get_ref().get_ref().as_slice());
        let mut record = bam::Record::default();
        reader.read_record(&mut record).unwrap();
        let records: Vec<bam::Record> = vec![record];
        let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
        let base_quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
        let bam_records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
        let alignment_model: AlignmentModel = AlignmentModel::new(1, &*read_sequence, &base_quality_scores, &bam_records);
        variant_records_by_form.push(
            DNAVariantRecordCaller::new(4, 20, 0).call(&alignment_model)
                .iter()
                .map(|variant_record| variant_record.get_graph_operation_boxed_str())
                .collect()
        );
    }

    assert_eq!(variant_records_by_form[0], vec![Box::<str>::from("0:110:+:D:0:113:+:U:TA:2:MNV")]);
    assert_eq!(variant_records_by_form[1], variant_records_by_form[0]);
}


/// A read stored without base qualities (SAM `*`) is given Q60 on every base, so its variants
/// pass the minimum base quality.
#[test]
fn identify_variant_records_gives_a_read_without_base_qualities_q60() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;

    let sam_text: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nread-1\t0\tchr1\t101\t60\t10=2X8=\t*\t0\t0\tACGTACGTACTAACGTACGT\t*\tcs:Z::10*gt*ga:8\n";
    let mut sam_reader = sam::io::Reader::new(sam_text.as_bytes());
    let header = sam_reader.read_header().unwrap();
    let record_buf = sam_reader.record_bufs(&header).next().unwrap().unwrap();
    let mut writer = bam::io::Writer::new(Vec::new());
    writer.write_alignment_record(&header, &record_buf).unwrap();
    writer.try_finish().unwrap();
    let mut reader = bam::io::Reader::new(writer.get_ref().get_ref().as_slice());
    let mut record = bam::Record::default();
    reader.read_record(&mut record).unwrap();
    let records: Vec<bam::Record> = vec![record];
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
    let base_quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
    assert!(base_quality_scores.is_empty());
    let bam_records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
    let alignment_model: AlignmentModel = AlignmentModel::new(1, &*read_sequence, &base_quality_scores, &bam_records);

    let variant_records: Vec<VariantRecord> = DNAVariantRecordCaller::new(4, 20, 0).call(&alignment_model);

    assert!((0..20).all(|i| alignment_model.get_base(i).get_base_quality() == 60));
    assert_eq!(variant_records.len(), 1);
    assert_eq!(&*variant_records[0].get_standardized_sequence(), "TA");
}


/// A split read of 160 bases: 10 bases that align nowhere, 50 on chr1 from 501 (supplementary),
/// 100 on chr1 from 101 (primary). The supplementary record is the first in the read, and its
/// leading clip is the 10 unaligned bases. Written with a hard clip (minimap2 without -Y) and
/// with a soft clip (-Y), the record gives the same model and the same variant records: a hard
/// clip stands for read bases as a soft clip does, since the read sequence is the primary's.
#[test]
fn alignment_model_reads_a_hard_clipped_supplementary_record_as_a_soft_clipped_one() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;

    let unaligned: String = "GGGGGCCCCC".to_string();
    let supplementary_part: String = "TTGACCATGG".repeat(5);
    let primary_part: String = "ACGTTGCAAC".repeat(10);
    let read: String = format!("{unaligned}{supplementary_part}{primary_part}");
    let qualities: String = "I".repeat(160);
    let mut variant_records_by_form: Vec<Vec<Box<str>>> = Vec::new();
    for supplementary in [
        format!("read-1\t2048\tchr1\t501\t60\t10H50=100H\t*\t0\t0\t{supplementary_part}\t{}\tcs:Z::50", "I".repeat(50)),
        format!("read-1\t2048\tchr1\t501\t60\t10S50=100S\t*\t0\t0\t{read}\t{qualities}\tcs:Z::50")
    ] {
        let sam_text: String = format!(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nread-1\t0\tchr1\t101\t60\t60S100=\t*\t0\t0\t{read}\t{qualities}\tcs:Z::100\n{supplementary}\n"
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
        let mut operations: Vec<Box<str>> = DNAVariantRecordCaller::new(4, 20, 4).call(&alignment_model)
            .iter()
            .map(|variant_record| variant_record.get_graph_operation_boxed_str())
            .collect();
        operations.sort();
        variant_records_by_form.push(operations);
    }

    // The breakend between the two records, and the leading clip of 10 bases.
    assert_eq!(variant_records_by_form[1].len(), 2, "{:?}", variant_records_by_form[1]);
    assert_eq!(variant_records_by_form[0], variant_records_by_form[1]);
}

#[test]
fn alignment_model_joins_the_records_either_side_of_a_contained_record() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;

    // A read of 300 bases: A aligns read 0-199 at chr1:101, B aligns read 50-99 (inside A) at
    // chr1:151, and C aligns read 200-299 at chr1:601. The breakend joins A to C.
    let mut state: u32 = 2_463_534_242;
    let read: String = (0..300)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            b"ACGT"[(state % 4) as usize] as char
        })
        .collect();
    let qualities: String = "I".repeat(300);
    let sam_text: String = format!(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n\
         read-1\t0\tchr1\t101\t60\t200=100S\t*\t0\t0\t{read}\t{qualities}\tcs:Z::200\n\
         read-1\t2048\tchr1\t151\t60\t50S50=200S\t*\t0\t0\t{read}\t{qualities}\tcs:Z::50\n\
         read-1\t2048\tchr1\t601\t60\t200S100=\t*\t0\t0\t{read}\t{qualities}\tcs:Z::100\n"
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

    // Read bases 100-199 stay aligned by A.
    assert!((100..200).all(|i| *alignment_model.get_base(i).get_kind() == AlignmentModelBaseKind::Match));

    let breakends: Vec<(u32, u32, usize)> = DNAVariantRecordCaller::new(4, 20, 4).call(&alignment_model)
        .iter()
        .filter(|variant_record| *variant_record.get_variant_type() == VariantType::Breakpoint)
        .map(|variant_record| (variant_record.get_position_1(), variant_record.get_position_2(), variant_record.get_sequence().len()))
        .collect();
    assert_eq!(breakends, vec![(300, 601, 0)]);
}
