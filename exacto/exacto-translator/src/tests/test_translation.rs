use exacto_caller::prelude::*;
use exacto_core::prelude::*;
use flate2::read::GzDecoder;
use noodles_fastq as fastq;
use polars::prelude::*;
use std::fs;
use std::fs::File;
use std::io::{BufReader,Read};
use std::path::Path;
use tempfile::NamedTempFile;

use crate::prelude::*;


// ---------------------------------------------------------------------------
// Helper utilities for building synthetic test DataFrames
// ---------------------------------------------------------------------------

/// Build a transcript_structures DataFrame from row tuples.
/// Each row: (transcript_model_id, reference_transcript_ids, index, read_start, read_end, sequence, type, kind, context, skipped)
fn make_transcript_structures_df(
    rows: Vec<(i64, &str, i64, i64, i64, &str, &str, &str, &str, &str)>
) -> DataFrame {
    let mut col_tmi: Vec<i64> = Vec::new();
    let mut col_rti: Vec<String> = Vec::new();
    let mut col_idx: Vec<i64> = Vec::new();
    let mut col_rs: Vec<i64> = Vec::new();
    let mut col_re: Vec<i64> = Vec::new();
    let mut col_seq: Vec<String> = Vec::new();
    let mut col_type: Vec<String> = Vec::new();
    let mut col_kind: Vec<String> = Vec::new();
    let mut col_ctx: Vec<String> = Vec::new();
    let mut col_skipped: Vec<String> = Vec::new();
    for (tmi, rti, idx, rs, re, seq, rtype, kind, ctx, skipped) in rows {
        col_tmi.push(tmi);
        col_rti.push(rti.to_string());
        col_idx.push(idx);
        col_rs.push(rs);
        col_re.push(re);
        col_seq.push(seq.to_string());
        col_type.push(rtype.to_string());
        col_kind.push(kind.to_string());
        col_ctx.push(ctx.to_string());
        col_skipped.push(skipped.to_string());
    }
    DataFrame::new(vec![
        Column::from(Series::new("transcript_model_id".into(), col_tmi)),
        Column::from(Series::new("reference_transcript_ids".into(), col_rti)),
        Column::from(Series::new("index".into(), col_idx)),
        Column::from(Series::new("read_start".into(), col_rs)),
        Column::from(Series::new("read_end".into(), col_re)),
        Column::from(Series::new("sequence".into(), col_seq)),
        Column::from(Series::new("type".into(), col_type)),
        Column::from(Series::new("kind".into(), col_kind)),
        Column::from(Series::new("context".into(), col_ctx)),
        Column::from(Series::new("skipped".into(), col_skipped)),
    ]).unwrap()
}

/// Build an empty integrated_variants DataFrame with the correct schema.
fn make_empty_integrated_variants_df() -> DataFrame {
    DataFrame::new(vec![
        Column::from(Series::new("transcript_model_id".into(), Vec::<i64>::new())),
        Column::from(Series::new("reference_transcript_ids".into(), Vec::<String>::new())),
        Column::from(Series::new("rna_variant_call_id".into(), Vec::<i64>::new())),
        Column::from(Series::new("dna_variant_call_id".into(), Vec::<i64>::new())),
    ]).unwrap()
}



// ===========================================================================
// translate_rnas edge-case tests
// ===========================================================================

#[test]
fn test_translate_rnas_sample_200_normal() {
    // Step 1. Read the RNA FASTQ file
    let fastq_path = Path::new("src/tests/data/fastq/sample200normal_long_read_rna.fastq.gz");
    let fastq_full_path = fs::canonicalize(fastq_path).unwrap();
    let gzipped = is_gzipped(fastq_full_path.to_str().unwrap());
    let file = File::open(fastq_full_path.to_str().unwrap()).expect("Unable to open FASTQ file");
    let reader: Box<dyn Read> = if gzipped {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let buffered_reader = BufReader::new(reader);
    let mut fastq_reader = fastq::Reader::new(buffered_reader);
    let mut rnas: Vec<RNA> = Vec::new();
    for result in fastq_reader.records() {
        match result {
            Ok(record) => {
                let sequence_result = String::from_utf8(record.sequence().to_vec());
                let sequence: String = match sequence_result {
                    Ok(seq) => seq,
                    Err(e) => {
                        panic!("Error converting sequence to UTF-8: {}", e);
                    }
                };
                let rna: RNA = RNA::new(
                    record.name().to_string().into_boxed_str(),
                    sequence.into_boxed_str(),
                );
                rnas.push(rna);
            }
            Err(e) => {
                panic!("Error reading record: {}", e);
            }
        }
    }

    // Step 2. Translate the RNA sequences
    let translation_set: TranslationSet = translate_rnas(
        rnas,
        2
    );

    // Step 3. Check for validity
    let mut found: bool = false;
    for translation in translation_set.translations.iter() {
        if translation.rna.id == "m64012_817037_637278/74/ccs".into() {
            if translation.get_longest_orf_peptide().sequence == "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD*".into() {
                found = true;
            }
            assert!(translation.get_peptides_count() == 38);
        }
    }
    assert!(found);
    assert!(translation_set.translations.len() == 200);
}

#[test]
fn test_translate_rnas_empty_vec() {
    let result = translate_rnas(vec![], 1);
    assert_eq!(result.translations.len(), 0);
}

#[test]
fn test_translate_rnas_no_start_codon() {
    // Sequence with no AUG — should produce no peptides, so filter_map drops it
    let rna = RNA::new("no_aug".into(), "GCUGCUGCUGCUGCU".into());
    let result = translate_rnas(vec![rna], 1);
    assert_eq!(result.translations.len(), 0);
}

#[test]
fn test_translate_rnas_stop_only() {
    // Only stop codons (UAA, UAG, UGA), no start codon
    let rna = RNA::new("stops".into(), "UAAUAGUGA".into());
    let result = translate_rnas(vec![rna], 1);
    assert_eq!(result.translations.len(), 0);
}

#[test]
fn test_translate_rnas_single_codon_orf() {
    // AUG immediately followed by UAA stop → peptide = "M*"
    let rna = RNA::new("minimal_orf".into(), "AUGUAA".into());
    let result = translate_rnas(vec![rna], 1);
    assert_eq!(result.translations.len(), 1);
    let translation = &result.translations[0];
    assert_eq!(translation.get_peptides_count(), 1);
    let peptide = translation.get_longest_orf_peptide();
    assert_eq!(&*peptide.sequence, "M*");
    assert_eq!(peptide.orf_start, 0);
    assert_eq!(peptide.orf_end, 5);
}

#[test]
fn test_translate_rnas_no_stop_codon() {
    // AUG GCU GCU — no stop codon, translation runs to end
    let rna = RNA::new("no_stop".into(), "AUGGCUGCU".into());
    let result = translate_rnas(vec![rna], 1);
    assert_eq!(result.translations.len(), 1);
    let peptide = result.translations[0].get_longest_orf_peptide();
    assert_eq!(&*peptide.sequence, "MAA");
    assert_eq!(peptide.orf_start, 0);
    assert_eq!(peptide.orf_end, 8); // last codon ends at position 8
}

#[test]
fn test_translate_rnas_multiple_orfs() {
    // Two AUGs at different positions → multiple peptides
    // Frame: AUG GCU AUG GCU UAA
    //   ORF1 starts at 0: AUGGCUAUGGCUUAA → "MAMА*" (but inner AUG also starts ORF2)
    //   ORF2 starts at 6: AUGGCUUAA → "MA*"
    let rna = RNA::new("multi_orf".into(), "AUGGCUAUGGCUUAA".into());
    let result = translate_rnas(vec![rna], 1);
    assert_eq!(result.translations.len(), 1);
    let translation = &result.translations[0];
    assert!(translation.get_peptides_count() >= 2);
    // Longest ORF should be the one starting at position 0
    let longest = translation.get_longest_orf_peptide();
    assert_eq!(longest.orf_start, 0);
}

#[test]
fn test_translate_rnas_short_sequence() {
    // Sequence shorter than 3 nt — can't form any codon
    let rna = RNA::new("short".into(), "AU".into());
    let result = translate_rnas(vec![rna], 1);
    assert_eq!(result.translations.len(), 0);
}

#[test]
fn test_translate_rnas_lowercase_input() {
    // Lowercase RNA sequence — translate() converts t→u internally
    // "augcgauag" → AUG CGA UAG → M R * (stop)
    let rna = RNA::new("lowercase".into(), "augcgauag".into());
    let result = translate_rnas(vec![rna], 1);
    assert_eq!(result.translations.len(), 1);
    let peptide = result.translations[0].get_longest_orf_peptide();
    assert_eq!(&*peptide.sequence, "MR*");
}

#[test]
fn test_translate_rnas_dna_input() {
    // DNA with T instead of U — translate() does T→U conversion
    // ATGCGATAG → AUG CGA UAG → M R *
    let rna = RNA::new("dna".into(), "ATGCGATAG".into());
    let result = translate_rnas(vec![rna], 1);
    assert_eq!(result.translations.len(), 1);
    let peptide = result.translations[0].get_longest_orf_peptide();
    assert_eq!(&*peptide.sequence, "MR*");
}

#[test]
fn test_translate_rnas_multithreaded() {
    // Multiple RNAs processed with multiple threads — all should be translated
    let rnas: Vec<RNA> = (0..20).map(|i| {
        RNA::new(format!("rna_{}", i).into(), "AUGGCUGCUUAA".into())
    }).collect();
    let result = translate_rnas(rnas, 4);
    assert_eq!(result.translations.len(), 20);
    for t in result.translations.iter() {
        assert_eq!(&*t.get_longest_orf_peptide().sequence, "MAA*");
    }
}

#[test]
fn test_translate_rnas_mixed_translatable_and_not() {
    // Mix of RNAs: some with ORFs, some without
    let rnas = vec![
        RNA::new("has_orf".into(), "AUGGCUUAA".into()),   // MA*
        RNA::new("no_orf".into(), "GCUGCUGCU".into()),     // no AUG
        RNA::new("has_orf2".into(), "AUGCCCUAG".into()),   // MP*
    ];
    let result = translate_rnas(rnas, 1);
    // Only the two with ORFs should be in the result
    assert_eq!(result.translations.len(), 2);
    let ids: Vec<&str> = result.translations.iter().map(|t| &*t.rna.id).collect();
    assert!(ids.contains(&"has_orf"));
    assert!(ids.contains(&"has_orf2"));
}


// ===========================================================================
// translate_transcript_structures edge-case tests
// ===========================================================================

#[test]
fn test_tts_rna_100_with_variants() {
    let tsv_path_1 = Path::new("src/tests/data/tsv/rna-100-tumor_minimap2_mdtagged_sorted_exacto_transcript_structures.tsv");
    let tsv_full_path_1 = fs::canonicalize(tsv_path_1).unwrap();
    let tsv_file_1: &str = tsv_full_path_1.to_str().unwrap();
    let tsv_path_2 = Path::new("src/tests/data/tsv/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_calls.tsv");
    let tsv_full_path_2 = fs::canonicalize(tsv_path_2).unwrap();
    let tsv_file_2: &str = tsv_full_path_2.to_str().unwrap();
    let tsv_path_3 = Path::new("src/tests/data/tsv/rna-100_dna-001_integration.tsv");
    let tsv_full_path_3 = fs::canonicalize(tsv_path_3).unwrap();
    let tsv_file_3: &str = tsv_full_path_3.to_str().unwrap();

    let df_transcript_structures = read_tsv_file(tsv_file_1);
    let rna_variant_call_set = RNAVariantCallSet::read_tsv_file(tsv_file_2);
    let df_integrated_variants = read_tsv_file(tsv_file_3);

    let primary_structure_set = translate_transcript_structures(
        &df_transcript_structures,
        &rna_variant_call_set,
        &df_integrated_variants,
        TranslationStrategy::LongestORF,
        1
    );

    assert!(primary_structure_set.primary_structures.len() == 6);

    let output_fasta_file: NamedTempFile = NamedTempFile::new().unwrap();
    let output_tsv_file: NamedTempFile = NamedTempFile::new().unwrap();

    primary_structure_set.to_fasta_file(output_fasta_file.path().to_str().unwrap());
    primary_structure_set.to_tsv_file(output_tsv_file.path().to_str().unwrap());
}

#[test]
fn test_tts_empty_dataframe() {
    let df = make_transcript_structures_df(vec![]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);
    assert_eq!(result.primary_structures.len(), 0);
}

#[test]
fn test_tts_no_orf_found() {
    // Sequence "GCUGCUGCU" has no AUG → no ORF → empty result
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 9, "GCUGCUGCU", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);
    assert_eq!(result.primary_structures.len(), 0);
}

#[test]
fn test_tts_single_exon_simple() {
    // Simple exon: AUGGCUGCUUAA → M A A * (4 amino acids, 12 nucleotides, 4 codons)
    // read positions 1..12 (1-based)
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 12, "AUGGCUGCUUAA", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    assert_eq!(result.primary_structures.len(), 1);
    let ps = &result.primary_structures[0];

    // Should have 12 base records
    let base_records: Vec<&PrimaryStructureRecord> = ps.records.iter()
        .filter(|r| *r.get_record_type() == PrimaryStructureRecordType::Base)
        .collect();
    assert_eq!(base_records.len(), 12);

    // Codon indices cycle 0,1,2,0,1,2,...
    for (i, rec) in base_records.iter().enumerate() {
        assert_eq!(rec.get_codon_index().unwrap(), (i % 3) as u8);
    }

    // All records should be InFrame (no variants)
    for rec in base_records.iter() {
        if let Some(fs) = rec.get_frameshift_state() {
            assert_eq!(*fs, FrameshiftState::InFrame);
        }
    }

    // First 3 bases (AUG) should have amino acid "M"
    assert_eq!(base_records[0].get_amino_acid().as_ref().unwrap().as_ref(), "M");
    assert_eq!(base_records[1].get_amino_acid().as_ref().unwrap().as_ref(), "M");
    assert_eq!(base_records[2].get_amino_acid().as_ref().unwrap().as_ref(), "M");

    // Bases 3-5 (GCU) → A
    assert_eq!(base_records[3].get_amino_acid().as_ref().unwrap().as_ref(), "A");

    // Bases 9-11 (UAA) → * (stop)
    assert_eq!(base_records[9].get_amino_acid().as_ref().unwrap().as_ref(), "*");
}

#[test]
fn test_tts_orf_position_alignment() {
    // Regression test for the 0-based→1-based ORF offset bug.
    // Sequence: GGG AUG GCU UAA = 12 nucleotides
    // The ORF starts at 0-based index 3 (the A of AUG), which is read position 4 (1-based).
    // The translation should NOT include the leading GGG.
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 12, "GGGAUGGCUUAA", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    assert_eq!(result.primary_structures.len(), 1);
    let ps = &result.primary_structures[0];

    let base_records: Vec<&PrimaryStructureRecord> = ps.records.iter()
        .filter(|r| *r.get_record_type() == PrimaryStructureRecordType::Base)
        .collect();

    // ORF is AUG GCU UAA = 9 nucleotides, not 12
    assert_eq!(base_records.len(), 9);

    // First base should be at read_start=4 (the 'A' of AUG), not 1 (the 'G' of GGG)
    assert_eq!(base_records[0].get_read_start(), 4);

    // First amino acid should be M (methionine), confirming ORF alignment
    assert_eq!(base_records[0].get_amino_acid().as_ref().unwrap().as_ref(), "M");
}

#[test]
fn test_tts_insertion_frameshift() {
    // A base row with kind=insertion causes net_variant_nucleotides_count to increment.
    // Sequence layout:
    //   Row 0 (exonic match): AUG → positions 1-3
    //   Row 1 (insertion):    C   → position 4-4 (inserted base, causes frameshift)
    //   Row 2 (exonic match): GCUUAA → positions 5-10
    // Combined RNA: AUGCGCUUAA
    // ORF from translate(): AUG CGC UUA A... → ORF covers positions 1..10 after offset
    // The insertion row increments net_variant_nucleotides_count by 1 → 1 % 3 != 0 → FrameShifted
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 3, "AUG", "base", "match", "exonic", ""),
        (1, "ENST001", 1, 4, 4, "C", "base", "insertion", "exonic", ""),
        (1, "ENST001", 2, 5, 10, "GCUUAA", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    assert_eq!(result.primary_structures.len(), 1);
    let ps = &result.primary_structures[0];

    // Find the insertion record
    let base_records: Vec<&PrimaryStructureRecord> = ps.records.iter()
        .filter(|r| *r.get_record_type() == PrimaryStructureRecordType::Base)
        .collect();

    // The inserted base should have incremented net_variant_nucleotides_count
    let insertion_rec = base_records.iter().find(|r| r.get_read_start() == 4).unwrap();
    assert_eq!(insertion_rec.get_net_variant_nucleotides_count(), 1);
}

#[test]
fn test_tts_event_with_skipped() {
    // An event record with skipped coordinates decrements net_variant_nucleotides_count.
    // skipped format: "label:pos1|label:pos2"
    // skipped_length = abs_diff(pos1, pos2) + 1
    // With "exon:10|exon:12": length = |10-12| + 1 = 3, so count decremented by 3
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 9, "AUGGCUGCU", "base", "match", "exonic", ""),
        (1, "ENST001", 1, 9, 12, "", "event", "splicing", "canonical", "exon:10|exon:12"),
        (1, "ENST001", 2, 13, 15, "UAA", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    assert_eq!(result.primary_structures.len(), 1);
    let ps = &result.primary_structures[0];

    // Find the event record
    let event_records: Vec<&PrimaryStructureRecord> = ps.records.iter()
        .filter(|r| *r.get_record_type() == PrimaryStructureRecordType::Event)
        .collect();
    assert_eq!(event_records.len(), 1);
    // net_variant_nucleotides_count should be -3 (decremented by skipped_length 3)
    assert_eq!(event_records[0].get_net_variant_nucleotides_count(), -3);
}

#[test]
fn test_tts_multiple_partitions() {
    // Two different transcript_model_ids → separate primary structures
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 12, "AUGGCUGCUUAA", "base", "match", "exonic", ""),
        (2, "ENST002", 0, 1, 12, "AUGCCCCCAUAA", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    assert_eq!(result.primary_structures.len(), 2);
}

#[test]
fn test_tts_all_orfs_strategy() {
    // Sequence with two AUGs → AllORFs strategy should produce multiple primary structures
    // AUG GCU AUG GCU UAA → two AUGs
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 15, "AUGGCUAUGGCUUAA", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::AllORFs, 1);

    // Should have at least 2 primary structures (one per ORF)
    assert!(result.primary_structures.len() >= 2,
        "Expected >=2 primary structures for AllORFs, got {}", result.primary_structures.len());
}

#[test]
fn test_tts_empty_integrated_variants() {
    // Verify that an empty integrated variants DF doesn't cause issues
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 9, "AUGGCUUAA", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    assert_eq!(result.primary_structures.len(), 1);
    // All dna_variant_call_ids should be empty
    for rec in result.primary_structures[0].records.iter() {
        assert!(rec.get_dna_variant_call_ids().is_empty());
    }
}

#[test]
fn test_tts_intronic_base_increments_count() {
    // Intronic match bases should increment net_variant_nucleotides_count
    // Layout: exonic AUG (positions 1-3), intronic GCU (positions 4-6), exonic UAA (positions 7-9)
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 3, "AUG", "base", "match", "exonic", ""),
        (1, "ENST001", 1, 4, 6, "GCU", "base", "match", "intronic", ""),
        (1, "ENST001", 2, 7, 9, "UAA", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    assert_eq!(result.primary_structures.len(), 1);
    let ps = &result.primary_structures[0];
    let base_records: Vec<&PrimaryStructureRecord> = ps.records.iter()
        .filter(|r| *r.get_record_type() == PrimaryStructureRecordType::Base)
        .collect();

    // Find intronic bases (positions 4-6) and verify count incremented
    let intronic_recs: Vec<&&PrimaryStructureRecord> = base_records.iter()
        .filter(|r| r.get_read_start() >= 4 && r.get_read_start() <= 6)
        .collect();
    for rec in intronic_recs.iter() {
        assert!(rec.get_net_variant_nucleotides_count() > 0,
            "Intronic match base should increment net_variant_nucleotides_count");
    }
}

#[test]
fn test_tts_incomplete_codon_at_end() {
    // 11 nucleotides: AUG GCU GCU GC
    // translate() only produces ORF up to last complete codon (position 8, 0-based).
    // So orf_end = 8 → after offset = 9. Only 9 nucleotides (3 codons) enter the
    // primary structure; the trailing 2 bases are outside the ORF entirely.
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 11, "AUGGCUGCUGC", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    assert_eq!(result.primary_structures.len(), 1);
    let ps = &result.primary_structures[0];
    let base_records: Vec<&PrimaryStructureRecord> = ps.records.iter()
        .filter(|r| *r.get_record_type() == PrimaryStructureRecordType::Base)
        .collect();

    // Only 9 bases (3 complete codons) should be in the primary structure
    assert_eq!(base_records.len(), 9);

    // All 9 bases should have amino acids assigned (complete codons only)
    for rec in base_records.iter() {
        assert!(rec.get_amino_acid().is_some(),
            "All bases within ORF should have amino acid assigned");
    }

    // Last base should be at read position 9, not 11
    assert_eq!(base_records.last().unwrap().get_read_start(), 9);
}

#[test]
fn test_tts_multithreaded_consistency() {
    // Same data processed with 1 thread vs 4 threads should produce same results
    let rows = vec![
        (1, "ENST001", 0, 1, 12, "AUGGCUGCUUAA", "base", "match", "exonic", ""),
        (2, "ENST002", 0, 1, 9, "AUGCCCUAG", "base", "match", "exonic", ""),
        (3, "ENST003", 0, 1, 12, "AUGAAAAAGUAA", "base", "match", "exonic", ""),
    ];

    let df1 = make_transcript_structures_df(rows.clone());
    let df2 = make_transcript_structures_df(rows);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv1 = make_empty_integrated_variants_df();
    let df_iv2 = make_empty_integrated_variants_df();

    let result_1t = translate_transcript_structures(&df1, &rna_vcs, &df_iv1, TranslationStrategy::LongestORF, 1);
    let result_4t = translate_transcript_structures(&df2, &rna_vcs, &df_iv2, TranslationStrategy::LongestORF, 4);

    assert_eq!(result_1t.primary_structures.len(), result_4t.primary_structures.len());
}

#[test]
fn test_tts_same_transcript_model_different_ref_transcripts() {
    // Same transcript_model_id but different reference_transcript_ids → separate partitions
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 12, "AUGGCUGCUUAA", "base", "match", "exonic", ""),
        (1, "ENST002", 0, 1, 9, "AUGCCCUAG", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    assert_eq!(result.primary_structures.len(), 2);
}

#[test]
fn test_tts_multi_exon_spliced() {
    // Multi-exon transcript: two exon rows joined by a splicing event.
    // Read positions are contiguous in the read coordinate space:
    //   Exon 1: AUG GCU (read positions 1-6)
    //   Splicing event (read positions 6-7, no sequence)
    //   Exon 2: GCU UAA (read positions 7-12)
    // Combined RNA: AUGGCUGCUUAA → M A A *
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 6, "AUGGCU", "base", "match", "exonic", ""),
        (1, "ENST001", 1, 6, 7, "", "event", "splicing", "canonical", ""),
        (1, "ENST001", 2, 7, 12, "GCUUAA", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    assert_eq!(result.primary_structures.len(), 1);
    let ps = &result.primary_structures[0];

    let base_records: Vec<&PrimaryStructureRecord> = ps.records.iter()
        .filter(|r| *r.get_record_type() == PrimaryStructureRecordType::Base)
        .collect();
    // 12 nucleotides total across both exons
    assert_eq!(base_records.len(), 12);

    // First amino acid should be M
    assert_eq!(base_records[0].get_amino_acid().as_ref().unwrap().as_ref(), "M");
    // Last codon (bases 9,10,11) should be * (stop)
    assert_eq!(base_records[9].get_amino_acid().as_ref().unwrap().as_ref(), "*");

    // Should also have an event record between the exons
    let event_count = ps.records.iter()
        .filter(|r| *r.get_record_type() == PrimaryStructureRecordType::Event)
        .count();
    assert_eq!(event_count, 1);
}

#[test]
fn test_tts_one_leading_base_with_stop() {
    // Sequence: G AUG GCU UAA — 1 leading nucleotide before ATG, with stop codon
    // 10 nucleotides at read positions 1-10
    // translate() finds AUG at 0-based index 1 → orf_start=1, orf_end=9
    // After offset (+1): orf_start=2, orf_end=10
    // Primary structure should contain only the 9 ORF nucleotides (positions 2-10),
    // excluding the leading G at position 1.
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 10, "GAUGGCUUAA", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    assert_eq!(result.primary_structures.len(), 1);
    let ps = &result.primary_structures[0];
    let base_records: Vec<&PrimaryStructureRecord> = ps.records.iter()
        .filter(|r| *r.get_record_type() == PrimaryStructureRecordType::Base)
        .collect();

    // 9 nucleotides in ORF (AUG GCU UAA), not 10
    assert_eq!(base_records.len(), 9);

    // First base should be at read position 2 (skipping the leading G at position 1)
    assert_eq!(base_records[0].get_read_start(), 2);

    // First codon → M
    assert_eq!(base_records[0].get_amino_acid().as_ref().unwrap().as_ref(), "M");
    assert_eq!(base_records[1].get_amino_acid().as_ref().unwrap().as_ref(), "M");
    assert_eq!(base_records[2].get_amino_acid().as_ref().unwrap().as_ref(), "M");

    // Second codon → A (GCU)
    assert_eq!(base_records[3].get_amino_acid().as_ref().unwrap().as_ref(), "A");

    // Third codon → * (UAA stop)
    assert_eq!(base_records[6].get_amino_acid().as_ref().unwrap().as_ref(), "*");

    // Last base at read position 10
    assert_eq!(base_records.last().unwrap().get_read_start(), 10);

    // All should be InFrame (no variants)
    for rec in base_records.iter() {
        assert_eq!(*rec.get_frameshift_state().as_ref().unwrap(), FrameshiftState::InFrame);
    }
}

#[test]
fn test_tts_one_leading_base_no_stop_incomplete_codon() {
    // Sequence: C AUG GCU GCU GC — 1 leading nucleotide, no stop codon, ends with incomplete codon
    // 12 nucleotides at read positions 1-12
    // translate() finds AUG at 0-based index 1:
    //   codon_start=1: AUG → M, orf_end=3
    //   codon_start=4: GCU → A, orf_end=6
    //   codon_start=7: GCU → A, orf_end=9
    //   codon_start=10: GC — only 2 chars (10+2=12 >= 12), loop breaks
    // So orf_start=1, orf_end=9 (0-based). After offset (+1): orf_start=2, orf_end=10.
    // Primary structure: 9 nucleotides at positions 2-10, 3 complete codons (MAA).
    // Leading C (pos 1) and trailing GC (pos 11-12) are outside the ORF.
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 12, "CAUGGCUGCUGC", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    assert_eq!(result.primary_structures.len(), 1);
    let ps = &result.primary_structures[0];
    let base_records: Vec<&PrimaryStructureRecord> = ps.records.iter()
        .filter(|r| *r.get_record_type() == PrimaryStructureRecordType::Base)
        .collect();

    // 9 nucleotides (3 complete codons), not 12
    assert_eq!(base_records.len(), 9);

    // First base at position 2 (leading C excluded)
    assert_eq!(base_records[0].get_read_start(), 2);

    // Last base at position 10 (trailing GC at 11-12 excluded)
    assert_eq!(base_records.last().unwrap().get_read_start(), 10);

    // All 3 codons should have amino acids: M, A, A
    assert_eq!(base_records[0].get_amino_acid().as_ref().unwrap().as_ref(), "M");
    assert_eq!(base_records[3].get_amino_acid().as_ref().unwrap().as_ref(), "A");
    assert_eq!(base_records[6].get_amino_acid().as_ref().unwrap().as_ref(), "A");

    // All should be InFrame
    for rec in base_records.iter() {
        assert_eq!(*rec.get_frameshift_state().as_ref().unwrap(), FrameshiftState::InFrame);
    }
}


// ===========================================================================
// Codon-level and frameshift variant ID propagation tests
// ===========================================================================

/// Helper: build a minimal GraphOperation for test VariantRecords.
fn make_test_graph_operation() -> GraphOperation {
    GraphOperation::new(
        0, 100, Strand::Forward, GraphOperationType::Downstream,
        0, 200, Strand::Forward, GraphOperationType::Upstream,
        "A".into(), VariantType::SingleNucleotideVariant
    )
}

#[test]
fn test_tts_codon_variant_propagation() {
    // A mismatch at one nucleotide in a codon should propagate to all 3 nucleotides
    // via codon_rna_variant_call_ids.
    //
    // Layout:
    //   Row 0 (match, exonic): AUG — positions 1-3
    //   Row 1 (mismatch, exonic): G — position 4 (this is the variant)
    //   Row 2 (match, exonic): CU — positions 5-6
    //   Row 3 (match, exonic): UAA — positions 7-9
    // Combined RNA: AUGGCUUAA → M A *
    // The mismatch at position 4 changes the middle codon (GCU → amino acid A).
    // Variant 42 has read_position_1=4, read_position_2=4, matching row 1.

    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 3, "AUG", "base", "match", "exonic", ""),
        (1, "ENST001", 1, 4, 4, "G", "base", "mismatch", "exonic", ""),
        (1, "ENST001", 2, 5, 6, "CU", "base", "match", "exonic", ""),
        (1, "ENST001", 3, 7, 9, "UAA", "base", "match", "exonic", ""),
    ]);

    // Build RNAVariantCallSet with variant 42 at read positions 4-4
    let mut rna_vcs = RNAVariantCallSet::new();
    let mut vc = VariantCall::new(42);
    vc.add_variant_record(VariantRecord::new(1, 4, 4, make_test_graph_operation()));
    rna_vcs.add_variant_call(
        1,
        vec!["ENST001".into()],
        vc
    );

    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    assert_eq!(result.primary_structures.len(), 1);
    let ps = &result.primary_structures[0];
    let base_records: Vec<&PrimaryStructureRecord> = ps.records.iter()
        .filter(|r| *r.get_record_type() == PrimaryStructureRecordType::Base)
        .collect();
    assert_eq!(base_records.len(), 9);

    // Direct rna_variant_call_ids: only the nucleotide at position 4 has it
    assert!(base_records[0].get_rna_variant_call_ids().is_empty()); // pos 1
    assert!(base_records[1].get_rna_variant_call_ids().is_empty()); // pos 2
    assert!(base_records[2].get_rna_variant_call_ids().is_empty()); // pos 3
    assert_eq!(base_records[3].get_rna_variant_call_ids(), &vec![42usize]); // pos 4 — variant
    assert!(base_records[4].get_rna_variant_call_ids().is_empty()); // pos 5
    assert!(base_records[5].get_rna_variant_call_ids().is_empty()); // pos 6

    // Codon-level: all 3 nucleotides in the affected codon (positions 4,5,6) should have variant 42
    assert!(base_records[0].get_codon_rna_variant_call_ids().is_empty()); // codon 1: AUG — no variant
    assert!(base_records[1].get_codon_rna_variant_call_ids().is_empty());
    assert!(base_records[2].get_codon_rna_variant_call_ids().is_empty());
    assert_eq!(base_records[3].get_codon_rna_variant_call_ids(), &vec![42usize]); // codon 2: GCU
    assert_eq!(base_records[4].get_codon_rna_variant_call_ids(), &vec![42usize]);
    assert_eq!(base_records[5].get_codon_rna_variant_call_ids(), &vec![42usize]);
    assert!(base_records[6].get_codon_rna_variant_call_ids().is_empty()); // codon 3: UAA — no variant
    assert!(base_records[7].get_codon_rna_variant_call_ids().is_empty());
    assert!(base_records[8].get_codon_rna_variant_call_ids().is_empty());

    // Frameshift: all InFrame (mismatch doesn't change nucleotide count), so no frameshift IDs
    for rec in base_records.iter() {
        assert!(rec.get_frameshift_rna_variant_call_ids().is_empty());
    }
}

#[test]
fn test_tts_frameshift_variant_propagation() {
    // An insertion causes a frameshift; all downstream FrameShifted records should
    // carry the causal variant in frameshift_rna_variant_call_ids.
    //
    // Layout:
    //   Row 0 (match, exonic): AUG — positions 1-3
    //   Row 1 (insertion, exonic): C — position 4 (inserted, causes +1 frameshift)
    //   Row 2 (match, exonic): GCUGCUUAA — positions 5-13
    // Combined RNA: AUGCGCUGCUUAA
    // Variant 99 has read_position_1=4, read_position_2=4, matching row 1.
    //
    // Codon 1 (AUG, pos 1-3): InFrame, net_var=0
    // Codon 2 (CGC, pos 4-6): net_var at pos 6 = 1 → 1%3!=0 → FrameShifted
    // Codon 3 (UGC, pos 7-9): net_var = 1 → FrameShifted
    // Codon 4 (UUA, pos 10-12): net_var = 1 → FrameShifted

    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 3, "AUG", "base", "match", "exonic", ""),
        (1, "ENST001", 1, 4, 4, "C", "base", "insertion", "exonic", ""),
        (1, "ENST001", 2, 5, 13, "GCUGCUUAA", "base", "match", "exonic", ""),
    ]);

    let mut rna_vcs = RNAVariantCallSet::new();
    let mut vc = VariantCall::new(99);
    vc.add_variant_record(VariantRecord::new(1, 4, 4, make_test_graph_operation()));
    rna_vcs.add_variant_call(1, vec!["ENST001".into()], vc);

    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    assert_eq!(result.primary_structures.len(), 1);
    let ps = &result.primary_structures[0];
    let base_records: Vec<&PrimaryStructureRecord> = ps.records.iter()
        .filter(|r| *r.get_record_type() == PrimaryStructureRecordType::Base)
        .collect();

    // Codon 1 (positions 1-3): InFrame, no frameshift IDs
    assert_eq!(*base_records[0].get_frameshift_state().as_ref().unwrap(), FrameshiftState::InFrame);
    assert!(base_records[0].get_frameshift_rna_variant_call_ids().is_empty());
    assert!(base_records[1].get_frameshift_rna_variant_call_ids().is_empty());
    assert!(base_records[2].get_frameshift_rna_variant_call_ids().is_empty());

    // Codon 2 (positions 4-6): FrameShifted, caused by variant 99
    assert_eq!(*base_records[3].get_frameshift_state().as_ref().unwrap(), FrameshiftState::FrameShifted);
    assert_eq!(base_records[3].get_frameshift_rna_variant_call_ids(), &vec![99usize]);
    assert_eq!(base_records[4].get_frameshift_rna_variant_call_ids(), &vec![99usize]);
    assert_eq!(base_records[5].get_frameshift_rna_variant_call_ids(), &vec![99usize]);

    // Codon 3 (positions 7-9): Still FrameShifted, still carries variant 99
    assert_eq!(*base_records[6].get_frameshift_state().as_ref().unwrap(), FrameshiftState::FrameShifted);
    assert_eq!(base_records[6].get_frameshift_rna_variant_call_ids(), &vec![99usize]);
    assert_eq!(base_records[7].get_frameshift_rna_variant_call_ids(), &vec![99usize]);
    assert_eq!(base_records[8].get_frameshift_rna_variant_call_ids(), &vec![99usize]);

    // Codon-level: variant 99 should be on codon 2 (where the insertion physically is)
    assert!(base_records[0].get_codon_rna_variant_call_ids().is_empty()); // codon 1
    assert_eq!(base_records[3].get_codon_rna_variant_call_ids(), &vec![99usize]); // codon 2
    assert!(base_records[6].get_codon_rna_variant_call_ids().is_empty()); // codon 3
}

#[test]
fn test_tts_no_variants_empty_propagation_columns() {
    // When no variants exist, all 4 new columns should be empty on every record.
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 12, "AUGGCUGCUUAA", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    let ps = &result.primary_structures[0];
    for rec in ps.records.iter() {
        assert!(rec.get_codon_rna_variant_call_ids().is_empty());
        assert!(rec.get_codon_dna_variant_call_ids().is_empty());
        assert!(rec.get_frameshift_rna_variant_call_ids().is_empty());
        assert!(rec.get_frameshift_dna_variant_call_ids().is_empty());
    }
}

#[test]
fn test_tts_event_record_empty_propagation_columns() {
    // Event records should have empty codon and frameshift variant ID fields.
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 6, "AUGGCU", "base", "match", "exonic", ""),
        (1, "ENST001", 1, 6, 7, "", "event", "splicing", "canonical", ""),
        (1, "ENST001", 2, 7, 12, "GCUUAA", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);

    let ps = &result.primary_structures[0];
    let event_records: Vec<&PrimaryStructureRecord> = ps.records.iter()
        .filter(|r| *r.get_record_type() == PrimaryStructureRecordType::Event)
        .collect();
    assert_eq!(event_records.len(), 1);
    assert!(event_records[0].get_codon_rna_variant_call_ids().is_empty());
    assert!(event_records[0].get_codon_dna_variant_call_ids().is_empty());
    assert!(event_records[0].get_frameshift_rna_variant_call_ids().is_empty());
    assert!(event_records[0].get_frameshift_dna_variant_call_ids().is_empty());
}

#[test]
fn test_tts_dataframe_has_new_columns() {
    // Verify the to_dataframe() output includes the 4 new columns (20 total).
    let df = make_transcript_structures_df(vec![
        (1, "ENST001", 0, 1, 9, "AUGGCUUAA", "base", "match", "exonic", ""),
    ]);
    let rna_vcs = RNAVariantCallSet::new();
    let df_iv = make_empty_integrated_variants_df();
    let result = translate_transcript_structures(&df, &rna_vcs, &df_iv, TranslationStrategy::LongestORF, 1);
    let out_df = result.to_dataframe();

    assert_eq!(out_df.width(), 21);
    assert!(out_df.column("codon_rna_variant_call_ids").is_ok());
    assert!(out_df.column("codon_dna_variant_call_ids").is_ok());
    assert!(out_df.column("frameshift_rna_variant_call_ids").is_ok());
    assert!(out_df.column("frameshift_dna_variant_call_ids").is_ok());
}
