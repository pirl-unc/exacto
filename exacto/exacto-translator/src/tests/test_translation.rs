use exacto_caller::prelude::*;
use exacto_integrator::prelude::*;
use polars::prelude::*;
use std::collections::HashSet;
use std::fs;
use std::path::Path;

use crate::prelude::*;


/// Wrap a `(id, sequence)` pair into the boxed-str shape that
/// `translate_sequences` consumes.
fn make_input(id: &str, sequence: &str) -> Vec<(Box<str>, Box<str>)> {
    vec![(id.to_string().into_boxed_str(), sequence.to_string().into_boxed_str())]
}

/// Convenience wrapper: translate one synthetic sequence with the default
/// AUG start-codon set on a single thread. Used by every edge-case test
/// below so the boilerplate stays in one place.
fn translate_one(sequence: &str, strategy: TranslationStrategy) -> TranscriptSet {
    translate_sequences(
        make_input("test", sequence),
        strategy,
        HashSet::from_iter(vec!["AUG"]),
        1,
    )
}

/// Collect the amino-acid string from a `PrimaryStructure`. Mirrors what
/// the FASTA writer would emit, so tests can assert against literal peptides.
fn amino_acid_sequence(ps: &PrimaryStructure) -> String {
    ps.amino_acids.iter().map(|aa| aa.get_amino_acid()).collect()
}


#[test]
fn test_translate_transcripts_1() {
    let ts: TranscriptSet = translate_sequences(
        vec![("1".to_string().into_boxed_str(), "GCGAUGGCUGAAAAACUGACUGGCCAUUAA".to_string().into_boxed_str())],
        TranslationStrategy::AllORFs,
        HashSet::from_iter(vec!["AUG"]),
        1
    );

    let ps: &PrimaryStructure = ts.transcripts
        .get(0)
        .unwrap()
        .primary_structures
        .get(0)
        .unwrap();

    assert_eq!(ts.len(), 1);
    assert_eq!(ps.orf_start, 3);
    assert_eq!(ps.orf_end, 29);
}

#[test]
fn test_translate_transcripts_2() {
    let assembly_tsv_path = Path::new("src/tests/data/tsv/rna-100_transcriptome_assembly_read_support.tsv");
    let assembly_tsv_full_path = fs::canonicalize(assembly_tsv_path).unwrap();
    let assembly_tsv_file: &str = assembly_tsv_full_path.to_str().unwrap();

    let transcript_structure_tsv_path = Path::new("src/tests/data/tsv/rna-100-tumor_minimap2_mdtagged_sorted_exacto_transcript_structures.tsv");
    let transcript_structure_tsv_full_path = fs::canonicalize(transcript_structure_tsv_path).unwrap();
    let transcript_structure_tsv_file: &str = transcript_structure_tsv_full_path.to_str().unwrap();

    let rna_variants_tsv_path = Path::new("src/tests/data/tsv/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_calls.tsv");
    let rna_variants_tsv_full_path = fs::canonicalize(rna_variants_tsv_path).unwrap();
    let rna_variants_tsv_file: &str = rna_variants_tsv_full_path.to_str().unwrap();

    let dna_variants_tsv_path = Path::new("src/tests/data/tsv/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv");
    let dna_variants_tsv_full_path = fs::canonicalize(dna_variants_tsv_path).unwrap();
    let dna_variants_tsv_file: &str = dna_variants_tsv_full_path.to_str().unwrap();

    let integrated_variants_tsv_path = Path::new("src/tests/data/tsv/rna-100_dna-001_integration.tsv");
    let integrated_variants_tsv_full_path = fs::canonicalize(integrated_variants_tsv_path).unwrap();
    let integrated_variants_tsv_file: &str = integrated_variants_tsv_full_path.to_str().unwrap();

    let assembly_support_records: Vec<AssembledTranscriptSupportRecord> = load_assembled_transcript_support_records(&assembly_tsv_file);
    let transcript_structure_records: Vec<TranscriptModelStructureRecord> = load_transcript_model_structure_records(&transcript_structure_tsv_file);
    let rna_variant_records: Vec<RNAVariantRecord> = load_rna_variant_records(&rna_variants_tsv_file);
    let dna_variant_records: Vec<DNAVariantRecord> = load_dna_variant_records(&dna_variants_tsv_file);
    let integrated_variant_records: Vec<IntegratedVariantRecord> = load_integrated_variant_records(&integrated_variants_tsv_file);

    let ts: TranscriptSet = translate_structures(
        &assembly_support_records,
        &transcript_structure_records,
        &rna_variant_records,
        &dna_variant_records,
        &integrated_variant_records,
        TranslationStrategy::LongestORF,
        HashSet::from_iter(vec!["AUG"]),
        1
    );

    // ------------------------------------------------------------------
    // Linkage assertion: the integration fixture links exactly one RNA
    // variant ↔ DNA variant pair (rna_variant_id=1 ↔ dna_variant_id=1, one
    // SNV at chr17:7674224, RNA read position 879). The same linkage must
    // surface in all four downstream representations:
    //   1. NucleotideRecord     (one row per ORF nucleotide)
    //   2. PrimaryStructureRecord (one row per ORF, ids aggregated)
    //   3. NucleotideRecord DataFrame
    //   4. PrimaryStructureRecord DataFrame
    // ------------------------------------------------------------------

    // --- (1) NucleotideRecord: exactly one nucleotide carries the link ---
    let nucleotide_records: Vec<NucleotideRecord> = build_nucleotide_records(&ts).collect();
    let snv_records: Vec<&NucleotideRecord> = nucleotide_records
        .iter()
        .filter(|r| r.rna_variant_id == Some(1))
        .collect();
    assert_eq!(
        snv_records.len(), 1,
        "expected exactly one nucleotide tagged with rna_variant_id=1, found {}",
        snv_records.len()
    );

    let snv: &NucleotideRecord = snv_records[0];
    let expected_dna_ids: HashSet<u32> = HashSet::from([1]);
    assert_eq!(
        snv.dna_variant_ids.as_ref(),
        Some(&expected_dna_ids),
        "SNV nucleotide should carry the linked DNA variant id from the integration fixture"
    );
    assert!(snv.is_nucleotide_variant, "SNV nucleotide should be flagged as a variant");

    // --- (2) PrimaryStructureRecord: the parent PS carries both ids ---
    let ps_records: Vec<PrimaryStructureRecord> = build_primary_structure_records(&ts).collect();
    let ps_with_rna_1: Vec<&PrimaryStructureRecord> = ps_records
        .iter()
        .filter(|r| r.rna_variant_ids.split(',').any(|id| id == "1"))
        .collect();
    assert_eq!(
        ps_with_rna_1.len(), 1,
        "expected exactly one PrimaryStructureRecord carrying rna_variant_id=1, found {}",
        ps_with_rna_1.len()
    );
    let ps_record: &PrimaryStructureRecord = ps_with_rna_1[0];
    assert!(
        ps_record.dna_variant_ids.split(',').any(|id| id == "1"),
        "PrimaryStructureRecord carrying rna_variant_id=1 should also carry dna_variant_id=1 \
         (got dna_variant_ids={:?})",
        ps_record.dna_variant_ids
    );
    assert_eq!(
        ps_record.primary_structure_id, snv.primary_structure_id as usize,
        "PS-record linkage row should belong to the same primary structure as the \
         nucleotide-record SNV row"
    );

    // --- (3) NucleotideRecord DataFrame: same linkage survives the conversion ---
    let nucleotide_df: DataFrame = nucleotide_records_to_dataframe(build_nucleotide_records(&ts));
    let nuc_rna_col = nucleotide_df.column("rna_variant_id").unwrap().u32().unwrap();
    let nuc_dna_col = nucleotide_df.column("dna_variant_ids").unwrap().str().unwrap();
    let nuc_linked_rows: Vec<usize> = nuc_rna_col
        .iter()
        .enumerate()
        .filter_map(|(i, v)| (v == Some(1)).then_some(i))
        .collect();
    assert_eq!(
        nuc_linked_rows.len(), 1,
        "nucleotide DataFrame should have exactly one rna_variant_id=1 row, found {}",
        nuc_linked_rows.len()
    );
    assert_eq!(
        nuc_dna_col.get(nuc_linked_rows[0]),
        Some("1"),
        "nucleotide DataFrame's linked row should have dna_variant_ids = '1'"
    );

    // --- (4) PrimaryStructureRecord DataFrame: same linkage survives the conversion ---
    let ps_df: DataFrame = primary_structure_records_to_dataframe(build_primary_structure_records(&ts));
    let ps_rna_col = ps_df.column("rna_variant_ids").unwrap().str().unwrap();
    let ps_dna_col = ps_df.column("dna_variant_ids").unwrap().str().unwrap();
    let ps_linked_rows: Vec<usize> = ps_rna_col
        .iter()
        .enumerate()
        .filter_map(|(i, v)| {
            v.and_then(|s| s.split(',').any(|id| id == "1").then_some(i))
        })
        .collect();
    assert_eq!(
        ps_linked_rows.len(), 1,
        "PrimaryStructure DataFrame should have exactly one row carrying rna_variant_id=1, found {}",
        ps_linked_rows.len()
    );
    let ps_dna_str: &str = ps_dna_col.get(ps_linked_rows[0]).unwrap();
    assert!(
        ps_dna_str.split(',').any(|id| id == "1"),
        "PS DataFrame's linked row should carry dna_variant_id=1 (got dna_variant_ids={:?})",
        ps_dna_str
    );
}


#[test]
fn test_translate_transcripts_2_duplicate_gov_canonicalization() {
    // Regression test: RNA variant calling can emit one row per matching
    // reference_transcript_id, so the same biological variant (same
    // GraphOperationView) appears with several distinct variant_ids in the
    // RNA-variants TSV. Before the canonicalization fix in
    // `build_transcript_set`, the `BiMap<variant_id, GOV>` would evict every
    // prior pair and retain only the last id — orphaning any integration
    // record keyed on an earlier id and leaving `dna_variant_ids` empty.
    //
    // This test reuses the rna-100 fixtures (which contain a single RNA
    // variant id=1) and appends a *duplicate* row with the same GOV but a
    // different variant_id (=999). Integration still references id=1.
    // Translation must surface rna_variant_id=1 ↔ dna_variant_id=1 on the
    // SNV nucleotide, proving the canonicalization remap holds.

    let assembly_tsv_full_path = fs::canonicalize(
        Path::new("src/tests/data/tsv/rna-100_transcriptome_assembly_read_support.tsv")
    ).unwrap();
    let transcript_structure_tsv_full_path = fs::canonicalize(
        Path::new("src/tests/data/tsv/rna-100-tumor_minimap2_mdtagged_sorted_exacto_transcript_structures.tsv")
    ).unwrap();
    let rna_variants_tsv_full_path = fs::canonicalize(
        Path::new("src/tests/data/tsv/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_calls.tsv")
    ).unwrap();
    let dna_variants_tsv_full_path = fs::canonicalize(
        Path::new("src/tests/data/tsv/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv")
    ).unwrap();
    let integrated_variants_tsv_full_path = fs::canonicalize(
        Path::new("src/tests/data/tsv/rna-100_dna-001_integration.tsv")
    ).unwrap();

    let assembly_support_records: Vec<AssembledTranscriptSupportRecord> =
        load_assembled_transcript_support_records(assembly_tsv_full_path.to_str().unwrap());
    let transcript_structure_records: Vec<TranscriptModelStructureRecord> =
        load_transcript_model_structure_records(transcript_structure_tsv_full_path.to_str().unwrap());
    let mut rna_variant_records: Vec<RNAVariantRecord> =
        load_rna_variant_records(rna_variants_tsv_full_path.to_str().unwrap());
    let dna_variant_records: Vec<DNAVariantRecord> =
        load_dna_variant_records(dna_variants_tsv_full_path.to_str().unwrap());
    let integrated_variant_records: Vec<IntegratedVariantRecord> =
        load_integrated_variant_records(integrated_variants_tsv_full_path.to_str().unwrap());

    // Inject a duplicate RNA variant row: same GOV as id=1, new id=999.
    // Ordering matters — id=1 must precede id=999 so that "first id wins"
    // canonicalization keeps id=1 (matching the integration TSV's reference).
    assert_eq!(rna_variant_records.len(), 1, "fixture should start with exactly one RNA variant");
    let mut duplicate: RNAVariantRecord = rna_variant_records[0].clone();
    duplicate.variant_id = 999;
    duplicate.reference_transcript_id = "ENST_DUPLICATE.1".into();
    rna_variant_records.push(duplicate);

    let ts: TranscriptSet = translate_structures(
        &assembly_support_records,
        &transcript_structure_records,
        &rna_variant_records,
        &dna_variant_records,
        &integrated_variant_records,
        TranslationStrategy::LongestORF,
        HashSet::from_iter(vec!["AUG"]),
        1,
    );

    let nucleotide_records: Vec<NucleotideRecord> = build_nucleotide_records(&ts).collect();

    // The SNV nucleotide must carry the *canonical* rna_variant_id=1
    // (not the duplicate id=999), and the DNA linkage must still resolve.
    let snv_records: Vec<&NucleotideRecord> = nucleotide_records
        .iter()
        .filter(|r| r.is_nucleotide_variant)
        .collect();
    assert_eq!(
        snv_records.len(), 1,
        "expected exactly one variant nucleotide, found {}",
        snv_records.len()
    );
    let snv = snv_records[0];
    assert_eq!(
        snv.rna_variant_id, Some(1),
        "SNV nucleotide should carry the canonical (first-seen) rna_variant_id=1, \
         not the duplicate id=999 — got {:?}",
        snv.rna_variant_id
    );
    let expected_dna_ids: HashSet<u32> = HashSet::from([1]);
    assert_eq!(
        snv.dna_variant_ids.as_ref(),
        Some(&expected_dna_ids),
        "SNV nucleotide should still carry dna_variant_id=1 after canonicalization"
    );

    // The aggregated PrimaryStructureRecord row should also carry both ids.
    let ps_records: Vec<PrimaryStructureRecord> = build_primary_structure_records(&ts).collect();
    let ps = ps_records.iter()
        .find(|r| r.rna_variant_ids.split(',').any(|id| id == "1"))
        .expect("expected a PrimaryStructureRecord carrying rna_variant_id=1");
    assert!(
        ps.dna_variant_ids.split(',').any(|id| id == "1"),
        "PS row should carry dna_variant_id=1 (got dna_variant_ids={:?})",
        ps.dna_variant_ids
    );
    // The canonical-id rule means rna_variant_ids should be "1" alone —
    // the duplicate id=999 shares the same GOV and is never registered.
    assert_eq!(
        ps.rna_variant_ids, "1",
        "PS row should expose only the canonical rna_variant_id=1, not 999"
    );
}


// ===========================================================================
// `translate_sequences` edge-case tests
//
// These mirror the pre-refactor `translate_rnas_*` tests that lived below
// (now removed). They exercise the synthetic-sequence path — the cheapest
// way to verify the core ORF-detection / codon-table / strategy logic
// without spinning up the full record-driven pipeline.
//
// Conventions:
//   * One input → one Transcript in `ts.transcripts`, regardless of whether
//     an ORF was found. If no ORF, `transcript.primary_structures` is empty.
//   * `orf_start` / `orf_end` are 0-indexed; `orf_end` is the position of
//     the last nucleotide of the last emitted codon (inclusive).
//   * Stop codons (UAA / UAG / UGA) become "*" and terminate the peptide.
// ===========================================================================

#[test]
fn test_translate_sequences_empty_vec() {
    // No inputs → no transcripts.
    let ts: TranscriptSet = translate_sequences(
        Vec::new(),
        TranslationStrategy::LongestORF,
        HashSet::from_iter(vec!["AUG"]),
        1,
    );
    assert_eq!(ts.transcripts.len(), 0);
}

#[test]
fn test_translate_sequences_no_start_codon() {
    // Sequence without AUG → Transcript is still created, but with no PSs.
    let ts: TranscriptSet = translate_one("GCUGCUGCUGCUGCU", TranslationStrategy::LongestORF);
    assert_eq!(ts.transcripts.len(), 1);
    assert_eq!(ts.transcripts[0].primary_structures.len(), 0);
}

#[test]
fn test_translate_sequences_stop_only() {
    // Only stop codons, no AUG → no ORFs.
    let ts: TranscriptSet = translate_one("UAAUAGUGA", TranslationStrategy::LongestORF);
    assert_eq!(ts.transcripts.len(), 1);
    assert_eq!(ts.transcripts[0].primary_structures.len(), 0);
}

#[test]
fn test_translate_sequences_single_codon_orf() {
    // Minimal ORF: AUG immediately followed by UAA → peptide "M*".
    let ts: TranscriptSet = translate_one("AUGUAA", TranslationStrategy::LongestORF);
    assert_eq!(ts.transcripts.len(), 1);
    assert_eq!(ts.transcripts[0].primary_structures.len(), 1);

    let ps: &PrimaryStructure = &ts.transcripts[0].primary_structures[0];
    assert_eq!(ps.orf_start, 0);
    assert_eq!(ps.orf_end, 5); // last nt of UAA at position 5
    assert_eq!(amino_acid_sequence(ps), "M*");
}

#[test]
fn test_translate_sequences_no_stop_codon() {
    // AUG GCU GCU — no stop codon. Translation runs to the end of the last
    // complete codon.
    let ts: TranscriptSet = translate_one("AUGGCUGCU", TranslationStrategy::LongestORF);
    assert_eq!(ts.transcripts.len(), 1);
    assert_eq!(ts.transcripts[0].primary_structures.len(), 1);

    let ps: &PrimaryStructure = &ts.transcripts[0].primary_structures[0];
    assert_eq!(ps.orf_start, 0);
    assert_eq!(ps.orf_end, 8); // last nt of the third codon
    assert_eq!(amino_acid_sequence(ps), "MAA");
}

#[test]
fn test_translate_sequences_multiple_orfs() {
    // Two AUGs at different positions → AllORFs strategy yields >=2 PSs;
    // LongestORF picks the one starting at position 0.
    let sequence: &str = "AUGGCUAUGGCUUAA";

    let ts_all: TranscriptSet = translate_one(sequence, TranslationStrategy::AllORFs);
    assert!(
        ts_all.transcripts[0].primary_structures.len() >= 2,
        "AllORFs should emit at least one PS per start codon, got {}",
        ts_all.transcripts[0].primary_structures.len()
    );

    let ts_longest: TranscriptSet = translate_one(sequence, TranslationStrategy::LongestORF);
    assert_eq!(ts_longest.transcripts[0].primary_structures.len(), 1);
    assert_eq!(ts_longest.transcripts[0].primary_structures[0].orf_start, 0);
}

#[test]
fn test_translate_sequences_short_sequence() {
    // Sequence shorter than a codon → no ORF possible.
    let ts: TranscriptSet = translate_one("AU", TranslationStrategy::LongestORF);
    assert_eq!(ts.transcripts.len(), 1);
    assert_eq!(ts.transcripts[0].primary_structures.len(), 0);
}

#[test]
fn test_translate_sequences_dna_input() {
    // DNA input (T instead of U). `translate()` normalizes T→U internally;
    // ATGCGATAG → AUG CGA UAG → "MR*".
    let ts: TranscriptSet = translate_one("ATGCGATAG", TranslationStrategy::LongestORF);
    assert_eq!(ts.transcripts[0].primary_structures.len(), 1);
    assert_eq!(amino_acid_sequence(&ts.transcripts[0].primary_structures[0]), "MR*");
}

#[test]
fn test_translate_sequences_lowercase_input() {
    // Lowercase input. `translate()` uppercases each codon before lookup;
    // "augcgauag" → AUG CGA UAG → "MR*".
    let ts: TranscriptSet = translate_one("augcgauag", TranslationStrategy::LongestORF);
    assert_eq!(ts.transcripts[0].primary_structures.len(), 1);
    assert_eq!(amino_acid_sequence(&ts.transcripts[0].primary_structures[0]), "MR*");
}

#[test]
fn test_translate_sequences_multithreaded() {
    // Multi-thread translation should be order-stable per input and produce
    // identical peptides for identical sequences.
    let inputs: Vec<(Box<str>, Box<str>)> = (0..20)
        .map(|i| (
            format!("rna_{}", i).into_boxed_str(),
            "AUGGCUGCUUAA".to_string().into_boxed_str(),
        ))
        .collect();

    let ts: TranscriptSet = translate_sequences(
        inputs,
        TranslationStrategy::LongestORF,
        HashSet::from_iter(vec!["AUG"]),
        4,
    );

    assert_eq!(ts.transcripts.len(), 20);
    for t in ts.transcripts.iter() {
        assert_eq!(t.primary_structures.len(), 1);
        // AUG GCU GCU UAA → M A A *
        assert_eq!(amino_acid_sequence(&t.primary_structures[0]), "MAA*");
    }
}

#[test]
fn test_translate_sequences_mixed_translatable_and_not() {
    // Mixed batch: each input maps to exactly one Transcript, but only the
    // sequences that contain an AUG produce a primary structure.
    let inputs: Vec<(Box<str>, Box<str>)> = vec![
        ("has_orf".to_string().into_boxed_str(),  "AUGGCUUAA".to_string().into_boxed_str()),   // MA*
        ("no_orf".to_string().into_boxed_str(),   "GCUGCUGCU".to_string().into_boxed_str()),    // no AUG
        ("has_orf2".to_string().into_boxed_str(), "AUGCCCUAG".to_string().into_boxed_str()),   // MP*
    ];

    let ts: TranscriptSet = translate_sequences(
        inputs,
        TranslationStrategy::LongestORF,
        HashSet::from_iter(vec!["AUG"]),
        1,
    );

    assert_eq!(ts.transcripts.len(), 3);

    // Per-input PS counts: 1, 0, 1.
    let ps_counts: Vec<usize> = ts.transcripts.iter()
        .map(|t| t.primary_structures.len())
        .collect();
    assert_eq!(ps_counts, vec![1, 0, 1]);

    // Verify each translatable transcript's peptide.
    assert_eq!(amino_acid_sequence(&ts.transcripts[0].primary_structures[0]), "MA*");
    assert_eq!(amino_acid_sequence(&ts.transcripts[2].primary_structures[0]), "MP*");
}
