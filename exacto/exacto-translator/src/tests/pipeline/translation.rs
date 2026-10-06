use bimap::BiMap;
use exacto_caller::prelude::*;
use exacto_core::prelude::*;
use exacto_integrator::prelude::*;
use polars::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

use super::*;


/// Wrap a `(id, sequence)` pair into the boxed-str shape that
/// `translate_sequences` consumes.
fn make_input(id: &str, sequence: &str) -> Vec<(Box<str>, Box<str>)> {
    vec![(id.to_string().into_boxed_str(), sequence.to_string().into_boxed_str())]
}

/// Convenience wrapper: translate one synthetic sequence with the default
/// AUG start-codon set on a single thread. Used by every edge-case test
/// below so the boilerplate stays in one place.
fn translate_one(sequence: &str, strategy: TranslationStrategy) -> AssembledTranscriptSet {
    translate_sequences(
        make_input("test", sequence),
        strategy,
        &HashSet::from_iter(vec!["AUG"]),
        1,
    ).unwrap()
}

/// Collect the amino-acid string from a `Peptide`. Mirrors what
/// the FASTA writer would emit, so tests can assert against literal peptides.
fn amino_acid_sequence(ps: &Proteoform) -> String {
    ps.get_sequence().to_string()
}


#[test]
fn test_translate_transcripts_1() {
    let ts: AssembledTranscriptSet = translate_sequences(
        vec![("1".to_string().into_boxed_str(), "GCGAUGGCUGAAAAACUGACUGGCCAUUAA".to_string().into_boxed_str())],
        TranslationStrategy::AllORFs,
        &HashSet::from_iter(vec!["AUG"]),
        1
    ).unwrap();

    let ps: &Proteoform = ts.transcripts
        .get(0)
        .unwrap()
        .proteoforms
        .get(0)
        .unwrap();

    assert_eq!(ts.len(), 1);
    assert_eq!(ps.orf_start, 3);
    assert_eq!(ps.orf_end, 29);
}

#[test]
fn test_translate_transcripts_2() {
    let assembly_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-translator/rna-100_transcriptome_assembly_read_support.tsv");
    let assembly_tsv_full_path = fs::canonicalize(assembly_tsv_path).unwrap();
    let assembly_tsv_file: &str = assembly_tsv_full_path.to_str().unwrap();

    let transcript_structure_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-translator/rna-100-tumor_minimap2_mdtagged_sorted_exacto_transcript_structures.tsv");
    let transcript_structure_tsv_full_path = fs::canonicalize(transcript_structure_tsv_path).unwrap();
    let transcript_structure_tsv_file: &str = transcript_structure_tsv_full_path.to_str().unwrap();

    let rna_variants_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-translator/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_calls.tsv");
    let rna_variants_tsv_full_path = fs::canonicalize(rna_variants_tsv_path).unwrap();
    let rna_variants_tsv_file: &str = rna_variants_tsv_full_path.to_str().unwrap();

    let dna_variants_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-translator/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv");
    let dna_variants_tsv_full_path = fs::canonicalize(dna_variants_tsv_path).unwrap();
    let dna_variants_tsv_file: &str = dna_variants_tsv_full_path.to_str().unwrap();

    let integrated_variants_tsv_path = Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-translator/rna-100_dna-001_integration.tsv");
    let integrated_variants_tsv_full_path = fs::canonicalize(integrated_variants_tsv_path).unwrap();
    let integrated_variants_tsv_file: &str = integrated_variants_tsv_full_path.to_str().unwrap();

    let assembly_support_records: Vec<AssembledTranscriptSupportRecord> = load_assembled_transcript_support_records(&assembly_tsv_file).unwrap();
    let transcript_structure_records: Vec<AssembledTranscriptModelAlignmentRecord> = load_assembled_transcript_model_alignment_records(&transcript_structure_tsv_file);
    let rna_variant_records: Vec<AssembledTranscriptVariantRecord> = load_assembled_transcript_variant_records(&rna_variants_tsv_file);
    let dna_variant_records: Vec<DNAVariantRecord> = load_dna_variant_records(&dna_variants_tsv_file);
    let integrated_variant_records: Vec<IntegratedVariantRecord> = load_integrated_variant_records(&integrated_variants_tsv_file);

    let ts: AssembledTranscriptSet = translate_transcripts(
        &assembly_support_records,
        &transcript_structure_records,
        &rna_variant_records,
        &dna_variant_records,
        &integrated_variant_records,
        &Vec::new(),
        TranslationStrategy::LongestORF,
        &HashSet::from_iter(vec!["AUG"]),
        1
    ).unwrap();

    // ------------------------------------------------------------------
    // Linkage assertion: the integration fixture links exactly one RNA
    // variant ↔ DNA variant pair (assembled_transcript_variant_id=1 ↔ dna_variant_id=1, one
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
        .filter(|r| r.assembled_transcript_variant_id == Some(1))
        .collect();
    assert_eq!(
        snv_records.len(), 1,
        "expected exactly one nucleotide tagged with assembled_transcript_variant_id=1, found {}",
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
    let ps_records: Vec<ProteoformRecord> = build_proteoform_records(&ts).collect();
    let ps_with_rna_1: Vec<&ProteoformRecord> = ps_records
        .iter()
        .filter(|r| r.assembled_transcript_variant_ids.split(LIST_SEPARATOR).any(|id| id == "1"))
        .collect();
    assert_eq!(
        ps_with_rna_1.len(), 1,
        "expected exactly one PrimaryStructureRecord carrying assembled_transcript_variant_id=1, found {}",
        ps_with_rna_1.len()
    );
    let ps_record: &ProteoformRecord = ps_with_rna_1[0];
    assert!(
        ps_record.dna_variant_ids.split(LIST_SEPARATOR).any(|id| id == "1"),
        "PrimaryStructureRecord carrying assembled_transcript_variant_id=1 should also carry dna_variant_id=1 \
         (got dna_variant_ids={:?})",
        ps_record.dna_variant_ids
    );
    assert_eq!(
        ps_record.proteoform_id, snv.proteoform_id as usize,
        "PS-record linkage row should belong to the same primary structure as the \
         nucleotide-record SNV row"
    );

    // --- (3) NucleotideRecord DataFrame: same linkage survives the conversion ---
    let nucleotide_df: DataFrame = nucleotide_records_to_dataframe(build_nucleotide_records(&ts));
    let nuc_rna_col = nucleotide_df.column("assembled_transcript_variant_id").unwrap().u32().unwrap();
    let nuc_dna_col = nucleotide_df.column("dna_variant_ids").unwrap().str().unwrap();
    let nuc_linked_rows: Vec<usize> = nuc_rna_col
        .iter()
        .enumerate()
        .filter_map(|(i, v)| (v == Some(1)).then_some(i))
        .collect();
    assert_eq!(
        nuc_linked_rows.len(), 1,
        "nucleotide DataFrame should have exactly one assembled_transcript_variant_id=1 row, found {}",
        nuc_linked_rows.len()
    );
    assert_eq!(
        nuc_dna_col.get(nuc_linked_rows[0]),
        Some("1"),
        "nucleotide DataFrame's linked row should have dna_variant_ids = '1'"
    );

    // --- (4) PrimaryStructureRecord DataFrame: same linkage survives the conversion ---
    let ps_df: DataFrame = proteoform_records_to_dataframe(build_proteoform_records(&ts));
    let ps_rna_col = ps_df.column("assembled_transcript_variant_ids").unwrap().str().unwrap();
    let ps_dna_col = ps_df.column("dna_variant_ids").unwrap().str().unwrap();
    let ps_linked_rows: Vec<usize> = ps_rna_col
        .iter()
        .enumerate()
        .filter_map(|(i, v)| {
            v.and_then(|s| s.split(LIST_SEPARATOR).any(|id| id == "1").then_some(i))
        })
        .collect();
    assert_eq!(
        ps_linked_rows.len(), 1,
        "PrimaryStructure DataFrame should have exactly one row carrying assembled_transcript_variant_id=1, found {}",
        ps_linked_rows.len()
    );
    let ps_dna_str: &str = ps_dna_col.get(ps_linked_rows[0]).unwrap();
    assert!(
        ps_dna_str.split(LIST_SEPARATOR).any(|id| id == "1"),
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
    // Translation must surface assembled_transcript_variant_id=1 ↔ dna_variant_id=1 on the
    // SNV nucleotide, proving the canonicalization remap holds.

    let assembly_tsv_full_path = fs::canonicalize(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-translator/rna-100_transcriptome_assembly_read_support.tsv")
    ).unwrap();
    let transcript_structure_tsv_full_path = fs::canonicalize(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-translator/rna-100-tumor_minimap2_mdtagged_sorted_exacto_transcript_structures.tsv")
    ).unwrap();
    let rna_variants_tsv_full_path = fs::canonicalize(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-translator/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_calls.tsv")
    ).unwrap();
    let dna_variants_tsv_full_path = fs::canonicalize(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-translator/dna-001-tumor_minimap2_mdtagged_sorted_exacto_somatic_variants.tsv")
    ).unwrap();
    let integrated_variants_tsv_full_path = fs::canonicalize(
        Path::new(env!("EXACTO_TEST_DATA")).join("exacto/exacto-translator/rna-100_dna-001_integration.tsv")
    ).unwrap();

    let assembly_support_records: Vec<AssembledTranscriptSupportRecord> =
        load_assembled_transcript_support_records(assembly_tsv_full_path.to_str().unwrap()).unwrap();
    let transcript_structure_records: Vec<AssembledTranscriptModelAlignmentRecord> =
        load_assembled_transcript_model_alignment_records(transcript_structure_tsv_full_path.to_str().unwrap());
    let mut rna_variant_records: Vec<AssembledTranscriptVariantRecord> =
        load_assembled_transcript_variant_records(rna_variants_tsv_full_path.to_str().unwrap());
    let dna_variant_records: Vec<DNAVariantRecord> =
        load_dna_variant_records(dna_variants_tsv_full_path.to_str().unwrap());
    let integrated_variant_records: Vec<IntegratedVariantRecord> =
        load_integrated_variant_records(integrated_variants_tsv_full_path.to_str().unwrap());

    // Inject a duplicate RNA variant row: same GOV as id=1, new id=999.
    // Ordering matters — id=1 must precede id=999 so that "first id wins"
    // canonicalization keeps id=1 (matching the integration TSV's reference).
    assert_eq!(rna_variant_records.len(), 1, "fixture should start with exactly one RNA variant");
    let mut duplicate: AssembledTranscriptVariantRecord = rna_variant_records[0].clone();
    duplicate.variant_id = 999;
    duplicate.reference_transcript_id = "ENST_DUPLICATE.1".into();
    rna_variant_records.push(duplicate);

    let ts: AssembledTranscriptSet = translate_transcripts(
        &assembly_support_records,
        &transcript_structure_records,
        &rna_variant_records,
        &dna_variant_records,
        &integrated_variant_records,
        &Vec::new(),
        TranslationStrategy::LongestORF,
        &HashSet::from_iter(vec!["AUG"]),
        1,
    ).unwrap();

    let nucleotide_records: Vec<NucleotideRecord> = build_nucleotide_records(&ts).collect();

    // The SNV nucleotide must carry the *canonical* assembled_transcript_variant_id=1
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
        snv.assembled_transcript_variant_id, Some(1),
        "SNV nucleotide should carry the canonical (first-seen) assembled_transcript_variant_id=1, \
         not the duplicate id=999 — got {:?}",
        snv.assembled_transcript_variant_id
    );
    let expected_dna_ids: HashSet<u32> = HashSet::from([1]);
    assert_eq!(
        snv.dna_variant_ids.as_ref(),
        Some(&expected_dna_ids),
        "SNV nucleotide should still carry dna_variant_id=1 after canonicalization"
    );

    // The aggregated PrimaryStructureRecord row should also carry both ids.
    let ps_records: Vec<ProteoformRecord> = build_proteoform_records(&ts).collect();
    let ps = ps_records.iter()
        .find(|r| r.assembled_transcript_variant_ids.split(LIST_SEPARATOR).any(|id| id == "1"))
        .expect("expected a PrimaryStructureRecord carrying assembled_transcript_variant_id=1");
    assert!(
        ps.dna_variant_ids.split(LIST_SEPARATOR).any(|id| id == "1"),
        "PS row should carry dna_variant_id=1 (got dna_variant_ids={:?})",
        ps.dna_variant_ids
    );
    // The canonical-id rule means assembled_transcript_variant_ids should be "1" alone —
    // the duplicate id=999 shares the same GOV and is never registered.
    assert_eq!(
        ps.assembled_transcript_variant_ids, "1",
        "PS row should expose only the canonical assembled_transcript_variant_id=1, not 999"
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
    let ts: AssembledTranscriptSet = translate_sequences(
        Vec::new(),
        TranslationStrategy::LongestORF,
        &HashSet::from_iter(vec!["AUG"]),
        1,
    ).unwrap();
    assert_eq!(ts.transcripts.len(), 0);
}

#[test]
fn test_translate_sequences_no_start_codon() {
    // Sequence without AUG → Transcript is still created, but with no PSs.
    let ts: AssembledTranscriptSet = translate_one("GCUGCUGCUGCUGCU", TranslationStrategy::LongestORF);
    assert_eq!(ts.transcripts.len(), 1);
    assert_eq!(ts.transcripts[0].proteoforms.len(), 0);
}

#[test]
fn test_translate_sequences_stop_only() {
    // Only stop codons, no AUG → no ORFs.
    let ts: AssembledTranscriptSet = translate_one("UAAUAGUGA", TranslationStrategy::LongestORF);
    assert_eq!(ts.transcripts.len(), 1);
    assert_eq!(ts.transcripts[0].proteoforms.len(), 0);
}

#[test]
fn test_translate_sequences_single_codon_orf() {
    // Minimal ORF: AUG immediately followed by UAA → peptide "M*".
    let ts: AssembledTranscriptSet = translate_one("AUGUAA", TranslationStrategy::LongestORF);
    assert_eq!(ts.transcripts.len(), 1);
    assert_eq!(ts.transcripts[0].proteoforms.len(), 1);

    let ps: &Proteoform = &ts.transcripts[0].proteoforms[0];
    assert_eq!(ps.orf_start, 0);
    assert_eq!(ps.orf_end, 5); // last nt of UAA at position 5
    assert_eq!(amino_acid_sequence(ps), "M*");
}

#[test]
fn test_translate_sequences_no_stop_codon() {
    // AUG GCU GCU — no stop codon. An ORF that never reaches a stop codon is an
    // incomplete CDS and must be excluded entirely (no primary structure).
    let ts: AssembledTranscriptSet = translate_one("AUGGCUGCU", TranslationStrategy::LongestORF);
    assert_eq!(ts.transcripts.len(), 1);
    assert_eq!(
        ts.transcripts[0].proteoforms.len(), 0,
        "an ORF lacking a proper stop codon must not produce a primary structure"
    );
}

#[test]
fn test_translate_sequences_multiple_orfs() {
    // Two AUGs at different positions → AllORFs strategy yields >=2 PSs;
    // LongestORF picks the one starting at position 0.
    let sequence: &str = "AUGGCUAUGGCUUAA";

    let ts_all: AssembledTranscriptSet = translate_one(sequence, TranslationStrategy::AllORFs);
    assert!(
        ts_all.transcripts[0].proteoforms.len() >= 2,
        "AllORFs should emit at least one PS per start codon, got {}",
        ts_all.transcripts[0].proteoforms.len()
    );

    let ts_longest: AssembledTranscriptSet = translate_one(sequence, TranslationStrategy::LongestORF);
    assert_eq!(ts_longest.transcripts[0].proteoforms.len(), 1);
    assert_eq!(ts_longest.transcripts[0].proteoforms[0].orf_start, 0);
}

#[test]
fn test_translate_sequences_short_sequence() {
    // Sequence shorter than a codon → no ORF possible.
    let ts: AssembledTranscriptSet = translate_one("AU", TranslationStrategy::LongestORF);
    assert_eq!(ts.transcripts.len(), 1);
    assert_eq!(ts.transcripts[0].proteoforms.len(), 0);
}

#[test]
fn test_translate_sequences_dna_input() {
    // DNA input (T instead of U). `translate()` normalizes T→U internally;
    // ATGCGATAG → AUG CGA UAG → "MR*".
    let ts: AssembledTranscriptSet = translate_one("ATGCGATAG", TranslationStrategy::LongestORF);
    assert_eq!(ts.transcripts[0].proteoforms.len(), 1);
    assert_eq!(amino_acid_sequence(&ts.transcripts[0].proteoforms[0]), "MR*");
}

#[test]
fn test_translate_sequences_lowercase_input() {
    // Lowercase input. `translate()` uppercases each codon before lookup;
    // "augcgauag" → AUG CGA UAG → "MR*".
    let ts: AssembledTranscriptSet = translate_one("augcgauag", TranslationStrategy::LongestORF);
    assert_eq!(ts.transcripts[0].proteoforms.len(), 1);
    assert_eq!(amino_acid_sequence(&ts.transcripts[0].proteoforms[0]), "MR*");
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

    let ts: AssembledTranscriptSet = translate_sequences(
        inputs,
        TranslationStrategy::LongestORF,
        &HashSet::from_iter(vec!["AUG"]),
        4,
    ).unwrap();

    assert_eq!(ts.transcripts.len(), 20);
    for t in ts.transcripts.iter() {
        assert_eq!(t.proteoforms.len(), 1);
        // AUG GCU GCU UAA → M A A *
        assert_eq!(amino_acid_sequence(&t.proteoforms[0]), "MAA*");
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

    let ts: AssembledTranscriptSet = translate_sequences(
        inputs,
        TranslationStrategy::LongestORF,
        &HashSet::from_iter(vec!["AUG"]),
        1,
    ).unwrap();

    assert_eq!(ts.transcripts.len(), 3);

    // Per-input PS counts: 1, 0, 1.
    let ps_counts: Vec<usize> = ts.transcripts.iter()
        .map(|t| t.proteoforms.len())
        .collect();
    assert_eq!(ps_counts, vec![1, 0, 1]);

    // Verify each translatable transcript's peptide.
    assert_eq!(amino_acid_sequence(&ts.transcripts[0].proteoforms[0]), "MA*");
    assert_eq!(amino_acid_sequence(&ts.transcripts[2].proteoforms[0]), "MP*");
}


/// Builds a minimal `TranscriptModelStructureRecord` for an `event` row whose
/// `kind`/`context` can be overridden. Coordinate/descriptor fields use values
/// copied from a real splicing-event fixture row so the GraphOperationView parses.
fn event_structure_record(kind: &str, context: &str) -> AssembledTranscriptModelAlignmentRecord {
    AssembledTranscriptModelAlignmentRecord {
        assembled_transcript_name: "m64012_507476_774164/1/ccs".into(),
        reference_gene_name: "TP53".into(),
        reference_transcript_id: "ENST00000269305.9".into(),
        index: 1,
        read_start: 113,
        read_end: 114,
        sequence: "".into(),
        record_type: "event".into(),
        kind: kind.into(),
        context: context.into(),
        chromosome_1: "chr17".into(),
        position_1: 7676622,
        operation_1: "D".into(),
        strand_1: "-".into(),
        chromosome_2: "chr17".into(),
        position_2: 7687377,
        operation_2: "U".into(),
        strand_2: "-".into(),
        reference_gene_id_1: "ENSG00000141510.18".into(),
        reference_transcript_id_1: "ENST00000269305.9".into(),
        reference_exon_id_1: "ENSE00002667911.1".into(),
        reference_gene_id_2: "ENSG00000141510.18".into(),
        reference_transcript_id_2: "ENST00000269305.9".into(),
        reference_exon_id_2: "ENSE00003753508.2".into(),
        skipped: "".into()
    }
}

/// An unclassified event is serialized with an empty `context`; building a
/// structure item from it must yield `Event { context: None }` rather than
/// panicking (regression for the `from_str("")` panic in translate-structs).
#[test]
fn test_build_transcript_structure_item_empty_event_context() {
    let record = event_structure_record("breakpoint", "");
    let item = build_transcript_alignment_item(&record).unwrap();
    match item.item_type {
        TranscriptAlignmentItemType::Event { kind, context } => {
            assert_eq!(kind, AlignmentModelEventKind::Breakpoint);
            assert_eq!(context, None, "empty event context must parse to None");
        }
        _ => panic!("expected an Event item type")
    }
}

/// A populated event context still round-trips to `Some(..)`.
#[test]
fn test_build_transcript_structure_item_populated_event_context() {
    let record = event_structure_record("splicing", "canonical");
    let item = build_transcript_alignment_item(&record).unwrap();
    match item.item_type {
        TranscriptAlignmentItemType::Event { kind, context } => {
            assert_eq!(kind, AlignmentModelEventKind::Splicing);
            assert_eq!(context, Some(AlignmentModelEventContext::CanonicalSplicing));
        }
        _ => panic!("expected an Event item type")
    }
}

/// Builds a minimal `DNAVariantRecord` from a variant grammar descriptor.
fn dna_variant_record(
    chromosome_1: &str,
    position_1: u32,
    operation_1: &str,
    chromosome_2: &str,
    position_2: u32,
    operation_2: &str
) -> DNAVariantRecord {
    DNAVariantRecord {
        origin: DNAVariantOrigin::Somatic.as_str().into(),
        variant_id: 1,
        chromosome_1: chromosome_1.into(),
        position_1,
        strand_1: "+".into(),
        operation_1: operation_1.into(),
        chromosome_2: chromosome_2.into(),
        position_2,
        strand_2: "+".into(),
        operation_2: operation_2.into(),
        sequence: "".into(),
        variant_size: Some(0),
        variant_type: "".into(),
        consensus_read_names: "".into(),
        num_consensus_read_names: 0,
        read_names: "".into(),
        num_read_names: 0
    }
}

/// A cycle-creating DNA variant (Upstream -> Downstream, same chromosome,
/// position_1 < position_2 — e.g. a tandem duplication) must build a
/// GraphOperationView with a defaulted cycle count instead of panicking.
/// Regression for the `GraphOperationView::new` cycle-count assertion firing
/// during translate-structs DNA variant loading.
#[test]
fn test_graph_operation_view_from_dna_variant_record_cycle_defaults_num_cycles() {
    // The exact shape from the reported panic: chr7:66126215:U -> chr7:66206600:D.
    let record = dna_variant_record("chr7", 66126215, "U", "chr7", 66206600, "D");
    let gov = graph_operation_view_from_dna_variant_record(&record, "t").unwrap();
    assert_eq!(
        gov.get_num_cycles(), &Some(1),
        "a cycle-creating DNA operation must default num_cycles to Some(1)"
    );
}

/// A non-cycle DNA variant (e.g. a deletion: Downstream -> Upstream) leaves
/// num_cycles as None — the default must apply only to the cycle-creating shape.
#[test]
fn test_graph_operation_view_from_dna_variant_record_non_cycle_keeps_none() {
    let record = dna_variant_record("chr17", 7674224, "D", "chr17", 7674226, "U");
    let gov = graph_operation_view_from_dna_variant_record(&record, "t").unwrap();
    assert_eq!(
        gov.get_num_cycles(), &None,
        "a non-cycle DNA operation must keep num_cycles as None"
    );
}

/// Builds a minimal `NucleotideRecord` with an overridable `dna_variant_ids` set.
fn nucleotide_record(primary_structure_id: u32, dna_variant_ids: Option<HashSet<u32>>) -> NucleotideRecord {
    NucleotideRecord {
        proteoform_id: primary_structure_id,
        assembled_transcript_name: "t".into(),
        amino_acid_index: 0,
        amino_acid: "M".into(),
        codon_index: 0,
        nucleotide: "A".into(),
        is_amino_acid_variant: dna_variant_ids.is_some(),
        is_nucleotide_variant: dna_variant_ids.is_some(),
        assembled_transcript_read_position: 0,
        assembled_transcript_alignment_index: None,
        assembled_transcript_variant_id: None,
        assembled_transcript_variant: None,
        dna_variant_ids,
        dna_variant: None,
        preceding_event_assembled_transcript_variant_id: None,
        preceding_event_assembled_transcript_variant: None,
        preceding_event_dna_variant_ids: None,
        preceding_event_dna_variant: None,
        is_reference_stitched: false
    }
}

/// Regression: serializing a `NucleotideRecord` whose `dna_variant_ids` is a
/// MULTI-element set must produce a single semicolon-joined field, not one field
/// per element. Previously the raw `HashSet` expanded into multiple CSV columns,
/// changing the record width and making the writer fail on the first such row —
/// which silently truncated the nucleotides TSV (only records up to the first
/// multi-element set were written).
#[test]
fn test_nucleotide_records_write_with_multielement_dna_variant_ids() {
    let records: Vec<NucleotideRecord> = vec![
        nucleotide_record(1, None),
        nucleotide_record(2, Some(HashSet::from([6026]))),       // single element
        nucleotide_record(3, Some(HashSet::from([12756, 5271]))) // multi element: the failing case
    ];

    let tmp = tempfile::NamedTempFile::new().unwrap();
    let result = exacto_core::prelude::write_tsv_file(records.into_iter(), tmp.path());
    assert!(
        result.is_ok(),
        "write_tsv_file must not fail on a multi-element dna_variant_ids set: {:?}",
        result.err()
    );

    // Header + all three data rows must be present (no truncation).
    let content = std::fs::read_to_string(tmp.path()).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines.len(), 4, "expected header + 3 data rows, got {}", lines.len());

    // The multi-element set is rendered sorted and semicolon-joined in one field.
    assert!(
        lines[3].contains("5271;12756"),
        "multi-element dna_variant_ids must serialize as one joined field; row was: {}",
        lines[3]
    );
}

/// Build a single-ORF Transcript whose codon at amino-acid index 2 carries the
/// given RNA-variant GraphOperationView, translate it, and return its one
/// PrimaryStructureRecord.
///
/// Sequence `AUG AAA AAA AAA AAA UAA` translates to `MKKKK*` — ORF `[0, 17]`,
/// six residues at indices 0..=5. The variant (id 1) spans transcript position 7
/// (the middle base of the codon at positions 6..=8, i.e. amino-acid index 2):
/// a single base whose read span `[7, 7]` is registered both in `rna_variants`
/// (for descriptor rendering / frameshift lookup) and in the read-span index
/// (so `get_nucleotide` resolves the variant id by read coordinate there).
fn primary_structure_record_with_variant_at_index_2(gov: GraphOperationView) -> ProteoformRecord {
    let annotation = TranscriptAlignmentAnnotation {
        position_1_annotation: Annotation { gene_id: None, transcript_id: None, exon_id: None },
        position_2_annotation: Annotation { gene_id: None, transcript_id: None, exon_id: None }
    };
    let item = AssembledTranscriptAlignmentRecord::new(
        0,
        7, // read_start
        7, // read_end (single base => amino-acid index 2, codon positions 6..=8)
        TranscriptAlignmentItemType::Base {
            kind: AlignmentModelBaseKind::Insertion,
            context: AlignmentModelBaseContext::Exonic
        },
        gov.clone(),
        annotation
    );
    let mut transcript_structure = AssembledTranscriptAlignment::new(1, Vec::new(), Vec::new());
    transcript_structure.add_item(item);

    let mut rna_variants: BiMap<u32, GraphOperationView> = BiMap::new();
    rna_variants.insert(1, gov);

    let mut transcript = AssembledTranscript::new(
        "t".into(),
        "AUGAAAAAAAAAAAAUAA".into(),
        Vec::new(),
        "t".into(),
        transcript_structure,
        rna_variants,
        BiMap::new(),
        HashMap::new(),
        HashMap::new(),
        vec![(7, 7, 1)], // read span [7,7] -> assembled_transcript_variant_id 1
        HashMap::from([(1, VariantType::Insertion)]),
        Vec::new()
    );
    transcript.translate(&TranslationStrategy::LongestORF, &HashSet::from_iter(vec!["AUG"]));

    let ts = AssembledTranscriptSet::new(vec![transcript]);
    let records: Vec<ProteoformRecord> = build_proteoform_records(&ts).collect();
    assert_eq!(records.len(), 1, "expected exactly one primary structure");
    records.into_iter().next().unwrap()
}

/// A frame-shifting RNA variant (net +1 base insertion) must flag the whole
/// downstream region as mutant, from its codon (amino-acid index 2) through the
/// final residue of the ORF (index 5) — not just the codon it lands in.
#[test]
fn test_mutant_intervals_frameshift_flags_downstream_region() {
    // Net delta +1: insertion of a single base between adjacent anchors.
    let frameshift_gov = GraphOperationView::new(
        "chr1",
        100,
        GraphOperationType::Downstream,
        Strand::Forward,
        "chr1",
        101,
        GraphOperationType::Upstream,
        Strand::Forward,
        "A",
        None
    );
    assert!(frameshift_gov.is_frameshift(), "fixture must be a frameshift");

    let record = primary_structure_record_with_variant_at_index_2(frameshift_gov);
    assert_eq!(
        record.mutant_amino_acid_intervals, "2:5",
        "a frameshift must flag from its codon to the ORF end"
    );
    assert_eq!(record.num_mutant_amino_acids, 4);
}

/// An in-frame indel (net +3 bases) is NOT a frameshift: only the codon it lands
/// in (amino-acid index 2) is flagged; the downstream region is untouched.
#[test]
fn test_mutant_intervals_in_frame_indel_flags_only_site() {
    // Net delta +3: insertion of three bases, divisible by 3 -> no frameshift.
    let in_frame_gov = GraphOperationView::new(
        "chr1",
        100,
        GraphOperationType::Downstream,
        Strand::Forward,
        "chr1",
        101,
        GraphOperationType::Upstream,
        Strand::Forward,
        "ACG",
        None
    );
    assert!(!in_frame_gov.is_frameshift(), "fixture must not be a frameshift");

    let record = primary_structure_record_with_variant_at_index_2(in_frame_gov);
    assert_eq!(
        record.mutant_amino_acid_intervals, "2:2",
        "an in-frame indel must flag only its own codon"
    );
    assert_eq!(record.num_mutant_amino_acids, 1);
}

/// A cryptic-exon / read-through `Match` base is linked to its CRX/UTR RNA
/// variant by READ-COORDINATE containment (the GOVs never match: the structure
/// base carries the matched genomic sequence, the variant call an empty one).
/// The base therefore carries the variant id, is flagged mutant, and the id
/// surfaces in the primary-structure record. Regression for SPRED1
/// (hunknown_294149), whose downstream cryptic exons were dropped from the
/// mutant set because Match bases were never linked to their variant.
#[test]
fn test_mutant_intervals_cryptic_exon_match_base_linked_by_read_span() {
    let annotation = TranscriptAlignmentAnnotation {
        position_1_annotation: Annotation { gene_id: None, transcript_id: None, exon_id: None },
        position_2_annotation: Annotation { gene_id: None, transcript_id: None, exon_id: None }
    };
    // Structure base carries the matched genomic sequence (non-empty)...
    let structure_gov = GraphOperationView::new(
        "chr1",
        100,
        GraphOperationType::Include,
        Strand::Forward,
        "chr1",
        200,
        GraphOperationType::Include,
        Strand::Forward,
        "ACGTACGT",
        None
    );
    // ...while the CRX/UTR variant call (id 5) carries an EMPTY sequence, so the
    // two GOVs are NOT equal — only the read span links them.
    let variant_gov = GraphOperationView::new(
        "chr1",
        100,
        GraphOperationType::Include,
        Strand::Forward,
        "chr1",
        200,
        GraphOperationType::Include,
        Strand::Forward,
        "",
        None
    );
    assert_ne!(structure_gov, variant_gov, "GOVs must differ (sequence field)");

    let item = AssembledTranscriptAlignmentRecord::new(
        0,
        7, // read_start
        7, // read_end -> amino-acid index 2
        TranscriptAlignmentItemType::Base {
            kind: AlignmentModelBaseKind::Match,
            context: AlignmentModelBaseContext::Intergenic
        },
        structure_gov,
        annotation
    );
    let mut transcript_structure = AssembledTranscriptAlignment::new(1, Vec::new(), Vec::new());
    transcript_structure.add_item(item);

    let mut rna_variants: BiMap<u32, GraphOperationView> = BiMap::new();
    rna_variants.insert(5, variant_gov);

    let mut transcript = AssembledTranscript::new(
        "t".into(),
        "AUGAAAAAAAAAAAAUAA".into(),
        Vec::new(),
        "t".into(),
        transcript_structure,
        rna_variants,
        BiMap::new(),
        HashMap::new(),
        HashMap::new(),
        // The variant's read span [5, 9] contains the cryptic-exon base at 7.
        vec![(5, 9, 5)],
        HashMap::from([(5, VariantType::CrypticExon)]),
        Vec::new()
    );
    transcript.translate(&TranslationStrategy::LongestORF, &HashSet::from_iter(vec!["AUG"]));
    let ts = AssembledTranscriptSet::new(vec![transcript]);

    // Primary-structure level: the cryptic-exon codons are flagged and the CRX
    // variant id surfaces. Read 5 (codon 1) has no Base row, so it is matched to
    // the variant span on its own position. A cryptic exon is not known to keep
    // the reading frame, so the mutant region runs from its first codon to the
    // end of the ORF.
    let ps_records: Vec<ProteoformRecord> = build_proteoform_records(&ts).collect();
    assert_eq!(ps_records.len(), 1);
    assert_eq!(
        ps_records[0].mutant_amino_acid_intervals, "1:5",
        "a cryptic-exon base must flag its amino acid, and every one after it, as mutant"
    );
    assert_eq!(ps_records[0].num_mutant_amino_acids, 5);
    assert_eq!(
        ps_records[0].assembled_transcript_variant_ids, "5",
        "the cryptic-exon base must carry its CRX/UTR variant id"
    );

    // Nucleotide level: the base at read 7 carries the variant id and is flagged.
    let nuc_records: Vec<NucleotideRecord> = build_nucleotide_records(&ts).collect();
    let cryptic: Vec<&NucleotideRecord> = nuc_records.iter()
        .filter(|r| r.assembled_transcript_read_position == 7)
        .collect();
    assert_eq!(cryptic.len(), 1);
    assert_eq!(cryptic[0].assembled_transcript_variant_id, Some(5));
    assert!(cryptic[0].is_nucleotide_variant, "cryptic-exon nucleotide must be a variant");
    assert!(cryptic[0].is_amino_acid_variant);

    // A codon with no covering structure base (and thus no variant span) stays
    // un-flagged — read-span matching does not over-report.
    let unflagged: Vec<&NucleotideRecord> = nuc_records.iter()
        .filter(|r| r.assembled_transcript_read_position == 0)
        .collect();
    assert_eq!(unflagged.len(), 1);
    assert!(!unflagged[0].is_nucleotide_variant);
    assert_eq!(unflagged[0].assembled_transcript_variant_id, None);
}


/// A minimal aligned block for `assembled_transcript_name`: one exon covering the whole of
/// `sequence`, the sequence its support row carries. Only the name, sequence and read span
/// matter to the join under test, so everything else is left blank.
fn alignment_record_named(assembled_transcript_name: &str, sequence: &str) -> AssembledTranscriptModelAlignmentRecord {
    AssembledTranscriptModelAlignmentRecord {
        assembled_transcript_name: assembled_transcript_name.into(),
        reference_gene_name: "".into(),
        reference_transcript_id: "".into(),
        index: 0,
        read_start: 0,
        read_end: sequence.len() as u32 - 1,
        sequence: sequence.into(),
        record_type: "base".into(),
        kind: "match".into(),
        context: "exonic".into(),
        chromosome_1: "chr17".into(),
        position_1: 7_687_377,
        operation_1: "I".into(),
        strand_1: "+".into(),
        chromosome_2: "chr17".into(),
        position_2: 7_687_377 + sequence.len() as u32 - 1,
        operation_2: "I".into(),
        strand_2: "+".into(),
        reference_gene_id_1: "".into(),
        reference_transcript_id_1: "".into(),
        reference_exon_id_1: "".into(),
        reference_gene_id_2: "".into(),
        reference_transcript_id_2: "".into(),
        reference_exon_id_2: "".into(),
        skipped: "".into()
    }
}


/// Both `read_names` conventions parse: the external assembly pipeline joins with `,`, every
/// table exacto writes itself uses `LIST_SEPARATOR`.
#[test]
fn test_split_read_names_accepts_both_producer_conventions() {
    let expected: Vec<Box<str>> = vec!["read-a/1/ccs".into(), "read-b/1/ccs".into()];
    assert_eq!(split_read_names("read-a/1/ccs,read-b/1/ccs"), expected);
    assert_eq!(
        split_read_names(&["read-a/1/ccs", "read-b/1/ccs"].join(LIST_SEPARATOR)),
        expected
    );

    // Empty field means "no reads recorded", not one nameless read.
    assert!(split_read_names("").is_empty());
    // Trailing/duplicated separators and padding must not mint empty read names.
    assert_eq!(split_read_names("read-a/1/ccs, read-b/1/ccs,"), expected);
}


/// A consensus row becomes a support row whose name is the cluster id as a string — which is
/// exactly how `determine-rna-consensus` names sequences in the FASTA that goes on to be
/// aligned and called, so the join against the alignments is exact.
#[test]
fn test_support_record_from_consensus_names_by_cluster_id() {
    let record: AssembledTranscriptSupportRecord = AssembledTranscriptSupportRecord::from_consensus(
        12,
        "ACGT".into(),
        ["read-a/1/ccs", "read-b/1/ccs"].join(LIST_SEPARATOR).into_boxed_str()
    );
    assert_eq!(record.assembled_transcript_name.as_ref(), "12");
    assert_eq!(record.sequence.as_ref(), "ACGT");
    assert_eq!(
        split_read_names(&record.read_names),
        vec![Box::<str>::from("read-a/1/ccs"), Box::<str>::from("read-b/1/ccs")]
    );
}


/// One alignment row and a support row that names a different transcript: the classic
/// consensus-vs-external-assembler naming mismatch (`"0"` against `cid_0_0`). It used to
/// produce a transcript with a stitched-together sequence and no reads; it must now abort.
#[test]
fn test_build_transcript_set_rejects_unmatched_transcript_names() {
    let support_records: Vec<AssembledTranscriptSupportRecord> = vec![
        AssembledTranscriptSupportRecord::from_consensus(0, "ACGT".into(), "read-a/1/ccs".into())
    ];
    let alignment_records: Vec<AssembledTranscriptModelAlignmentRecord> =
        vec![alignment_record_named("cid_0_0", "ACGT")];

    let result = build_transcript_set(
        &support_records,
        &alignment_records,
        &Vec::new(),
        &Vec::new(),
        &Vec::new(),
        &Vec::new()
    );
    assert!(matches!(
        result,
        Err(TranslatorError::UnmatchedTranscripts { num_unmatched: 1, num_transcripts: 1, .. })
    ));
}


/// The matching case, so the test above is not passing merely because the builder rejects
/// everything: the same alignment row named by cluster id joins and carries its reads through.
#[test]
fn test_build_transcript_set_accepts_consensus_named_transcripts() {
    let support_records: Vec<AssembledTranscriptSupportRecord> = vec![
        AssembledTranscriptSupportRecord::from_consensus(
            0,
            "ACGT".into(),
            ["read-a/1/ccs", "read-b/1/ccs"].join(LIST_SEPARATOR).into_boxed_str()
        )
    ];
    let alignment_records: Vec<AssembledTranscriptModelAlignmentRecord> =
        vec![alignment_record_named("0", "ACGT")];

    let transcript_set: AssembledTranscriptSet = build_transcript_set(
        &support_records,
        &alignment_records,
        &Vec::new(),
        &Vec::new(),
        &Vec::new(),
        &Vec::new()
    ).unwrap();

    let transcripts: Vec<&AssembledTranscript> = transcript_set.iter().collect();
    assert_eq!(transcripts.len(), 1);
    assert_eq!(transcripts[0].get_assembled_transcript_name(), "0");
    assert_eq!(transcripts[0].sequence.as_ref(), "ACGT");
    assert_eq!(transcripts[0].get_read_ids().len(), 2);
}


/// DNA read names are deduplicated per read, not per variant.
///
/// `build_dna_variant_records` joins this field with `LIST_SEPARATOR`, so splitting it on ','
/// matched nothing and inserted each variant's whole list into the dedup set as one blob. A
/// proteoform integrating two DNA variants that share a read then reported that read twice.
#[test]
fn test_split_read_names_deduplicates_shared_dna_reads() {
    let first: Box<str> = ["readA", "readB"].join(LIST_SEPARATOR).into_boxed_str();
    let second: Box<str> = ["readB", "readC"].join(LIST_SEPARATOR).into_boxed_str();

    let mut names: HashSet<&str> = HashSet::new();
    for blob in [&first, &second] {
        for name in split_read_names_borrowed(blob) {
            names.insert(name);
        }
    }
    let mut sorted: Vec<&str> = names.into_iter().collect();
    sorted.sort_unstable();

    assert_eq!(sorted, vec!["readA", "readB", "readC"]);
    assert_eq!(sorted.join(LIST_SEPARATOR), "readA;readB;readC");
}


// ===========================================================================
// Reference-stitched provenance
//
// `stitch-reference-rnas` pastes reference-transcript bases onto consensus
// termini and records the boundaries as stitch points. These tests pin the
// whole provenance path: stitch points -> `ReferenceStitchedSpans` intervals
// -> per-nucleotide / per-amino-acid flags -> the interval column and the
// per-nucleotide boolean in the output records.
// ===========================================================================

/// Stitch-point arithmetic: `[0, fp)` and `[tp, len)` are reference-stitched,
/// empty spans are dropped, and a pass-through (fp=0, tp=len) yields none.
#[test]
fn test_reference_stitched_spans_from_stitch_points() {
    // 3'-only stitch (the 14411 shape: fp=0, tp=1434, len=5262).
    let spans = ReferenceStitchedSpans::from_stitch_points(14411, 0, 1434, 5262).unwrap();
    assert_eq!(spans.assembled_transcript_name.as_ref(), "14411");
    assert_eq!(spans.intervals, vec![(1434, 5262)]);

    // Both ends stitched.
    let spans = ReferenceStitchedSpans::from_stitch_points(7, 5, 10, 12).unwrap();
    assert_eq!(spans.intervals, vec![(0, 5), (10, 12)]);

    // Pass-through: nothing stitched.
    let spans = ReferenceStitchedSpans::from_stitch_points(3, 0, 9, 9).unwrap();
    assert!(spans.intervals.is_empty());
}

/// Inverted stitch points must be an error rather than mint a negative span.
#[test]
fn test_reference_stitched_spans_rejects_inverted_points() {
    assert!(matches!(
        ReferenceStitchedSpans::from_stitch_points(1, 10, 5, 12),
        Err(TranslatorError::InvalidStitchPoints { five_prime: 10, three_prime: 5, stitched_length: 12, .. })
    ));
}

/// End-to-end through `build_transcript_set` + translation: a 5'-stitched
/// prefix flags exactly the nucleotides and amino acids it covers, in both
/// output record streams.
///
/// Sequence `AUGAAAAAAAAAAAAUAA` (18 nt, ORF [0,17], peptide MKKKK*): a 5'
/// stitch point of 7 marks read positions 0..=6 as reference-stitched, i.e.
/// codons 0 and 1 entirely and codon 2's first base — so amino acids 0..=2
/// are stitched (any-nucleotide reduction) and 3..=5 are not.
#[test]
fn test_reference_stitched_provenance_flows_to_records() {
    let support_records: Vec<AssembledTranscriptSupportRecord> = vec![
        AssembledTranscriptSupportRecord::from_consensus(
            0,
            "AUGAAAAAAAAAAAAUAA".into(),
            "read-a/1/ccs".into()
        )
    ];
    let alignment_records: Vec<AssembledTranscriptModelAlignmentRecord> =
        vec![alignment_record_named("0", "AUGAAAAAAAAAAAAUAA")];
    let stitched_spans: Vec<ReferenceStitchedSpans> =
        vec![ReferenceStitchedSpans::from_stitch_points(0, 7, 18, 18).unwrap()];

    let ts: AssembledTranscriptSet = translate_transcripts(
        &support_records,
        &alignment_records,
        &Vec::new(),
        &Vec::new(),
        &Vec::new(),
        &stitched_spans,
        TranslationStrategy::LongestORF,
        &HashSet::from_iter(vec!["AUG"]),
        1
    ).unwrap();

    // Per-nucleotide: stitched iff read position < 7.
    let nucleotide_records: Vec<NucleotideRecord> = build_nucleotide_records(&ts).collect();
    assert_eq!(nucleotide_records.len(), 18);
    for record in nucleotide_records.iter() {
        assert_eq!(
            record.is_reference_stitched,
            record.assembled_transcript_read_position < 7,
            "read position {} has the wrong is_reference_stitched flag",
            record.assembled_transcript_read_position
        );
    }

    // Per-proteoform: amino acids 0..=2 stitched, rendered in the uniform
    // start:end interval spelling with a matching count.
    let ps_records: Vec<ProteoformRecord> = build_proteoform_records(&ts).collect();
    assert_eq!(ps_records.len(), 1);
    assert_eq!(ps_records[0].reference_stitched_amino_acid_intervals, "0:2");
    assert_eq!(ps_records[0].num_reference_stitched_amino_acids, 3);
}

/// No stitch annotations -> the columns are inert: empty intervals, zero
/// count, every nucleotide unflagged.
#[test]
fn test_reference_stitched_provenance_absent_without_annotations() {
    let support_records: Vec<AssembledTranscriptSupportRecord> = vec![
        AssembledTranscriptSupportRecord::from_consensus(0, "AUGGCUUAA".into(), "read-a/1/ccs".into())
    ];
    let alignment_records: Vec<AssembledTranscriptModelAlignmentRecord> =
        vec![alignment_record_named("0", "AUGGCUUAA")];

    let ts: AssembledTranscriptSet = translate_transcripts(
        &support_records,
        &alignment_records,
        &Vec::new(),
        &Vec::new(),
        &Vec::new(),
        &Vec::new(),
        TranslationStrategy::LongestORF,
        &HashSet::from_iter(vec!["AUG"]),
        1
    ).unwrap();

    let ps_records: Vec<ProteoformRecord> = build_proteoform_records(&ts).collect();
    assert_eq!(ps_records[0].reference_stitched_amino_acid_intervals, "");
    assert_eq!(ps_records[0].num_reference_stitched_amino_acids, 0);
    assert!(build_nucleotide_records(&ts).all(|r| !r.is_reference_stitched));
}

/// A stitch annotation whose stitched_length disagrees with the transcript's
/// actual sequence length is a stale pairing (annotations from a different
/// stitch run) and must be an error, not silently mislabel provenance.
#[test]
fn test_build_transcript_set_rejects_mismatched_stitched_length() {
    let support_records: Vec<AssembledTranscriptSupportRecord> = vec![
        AssembledTranscriptSupportRecord::from_consensus(0, "AUGGCUUAA".into(), "read-a/1/ccs".into())
    ];
    let alignment_records: Vec<AssembledTranscriptModelAlignmentRecord> =
        vec![alignment_record_named("0", "AUGGCUUAA")];
    // Sequence is 9 nt; the annotation claims 12.
    let stitched_spans: Vec<ReferenceStitchedSpans> =
        vec![ReferenceStitchedSpans::from_stitch_points(0, 3, 12, 12).unwrap()];

    let result = build_transcript_set(
        &support_records,
        &alignment_records,
        &Vec::new(),
        &Vec::new(),
        &Vec::new(),
        &stitched_spans
    );
    assert!(matches!(
        result,
        Err(TranslatorError::StitchedLengthMismatch { sequence_length: 9, stitched_length: 12, .. })
    ));
}


// ===========================================================================
// Review fixes of 2026-09-30 (tasks/exacto-translator-code-review.pdf)
// ===========================================================================

/// F-01. The mutant region after an unannotated junction opens whatever the intron length.
/// The sequence is the same for both introns: exon 1 is ATG + 19 GCA (read 0-59), exon 2 is
/// 30 GCA + TAA + 10 bases (read 60-162), so the ORF is read 0-152, 51 residues, and the
/// junction precedes residue 20. The old rule took the intron length modulo 3 as the frame
/// change: 1,001 nt marked 20:50, 1,002 nt only 20:20.
#[test]
fn test_mutant_region_after_junction_does_not_depend_on_intron_length() {
    let sequence: String = format!("ATG{}{}TAAGCAGCAGCAG", "GCA".repeat(19), "GCA".repeat(30));
    for intron_length in [1001u32, 1002] {
        let acceptor: u32 = 1060 + intron_length;
        let support_records: Vec<AssembledTranscriptSupportRecord> = vec![AssembledTranscriptSupportRecord {
            assembled_transcript_name: "t".into(),
            sequence: sequence.clone().into_boxed_str(),
            read_names: "r1".into()
        }];
        let alignment_records: Vec<AssembledTranscriptModelAlignmentRecord> = vec![
            AssembledTranscriptModelAlignmentRecord {
                assembled_transcript_name: "t".into(), index: 0, read_start: 0, read_end: 59,
                sequence: sequence[0..60].into(), record_type: "base".into(), kind: "match".into(), context: "exonic".into(),
                chromosome_1: "chr1".into(), position_1: 1000, operation_1: "I".into(), strand_1: "+".into(),
                chromosome_2: "chr1".into(), position_2: 1059, operation_2: "I".into(), strand_2: "+".into(),
                ..Default::default()
            },
            AssembledTranscriptModelAlignmentRecord {
                assembled_transcript_name: "t".into(), index: 1, read_start: 59, read_end: 60,
                sequence: "".into(), record_type: "event".into(), kind: "splicing".into(), context: "noncanonical".into(),
                chromosome_1: "chr1".into(), position_1: 1059, operation_1: "D".into(), strand_1: "+".into(),
                chromosome_2: "chr1".into(), position_2: acceptor, operation_2: "U".into(), strand_2: "+".into(),
                ..Default::default()
            },
            AssembledTranscriptModelAlignmentRecord {
                assembled_transcript_name: "t".into(), index: 2, read_start: 60, read_end: 162,
                sequence: sequence[60..163].into(), record_type: "base".into(), kind: "match".into(), context: "exonic".into(),
                chromosome_1: "chr1".into(), position_1: acceptor, operation_1: "I".into(), strand_1: "+".into(),
                chromosome_2: "chr1".into(), position_2: acceptor + 102, operation_2: "I".into(), strand_2: "+".into(),
                ..Default::default()
            }
        ];
        let variant_records: Vec<AssembledTranscriptVariantRecord> = vec![AssembledTranscriptVariantRecord {
            variant_id: 1, assembled_transcript_name: "t".into(), variant_type: "NCS".into(),
            chromosome_1: "chr1".into(), position_1: 1059, operation_1: "D".into(), strand_1: "+".into(),
            chromosome_2: "chr1".into(), position_2: acceptor, operation_2: "U".into(), strand_2: "+".into(),
            read_start: 59, read_end: 60,
            ..Default::default()
        }];

        let ts: AssembledTranscriptSet = translate_transcripts(
            &support_records, &alignment_records, &variant_records, &Vec::new(), &Vec::new(), &Vec::new(),
            TranslationStrategy::LongestORF, &HashSet::from_iter(vec!["AUG"]), 1
        ).unwrap();

        let records: Vec<ProteoformRecord> = build_proteoform_records(&ts).collect();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].amino_acid_sequence_length, 51);
        assert_eq!(records[0].mutant_amino_acid_intervals, "20:50", "intron of {} nt", intron_length);
        assert_eq!(records[0].num_mutant_amino_acids, 31);
    }
}

/// F-01. A retained intron (an Include-to-Include region call) opens the mutant region: the
/// 100 retained bases put everything after them in another frame than the spliced transcript's.
/// The old rule never opened a region for a region call, and marked only residues 20:53.
#[test]
fn test_mutant_region_opens_at_retained_intron() {
    // Exon 1 = ATG + 19 GCA (read 0-59), retained intron = 33 GCA + G (read 60-159), exon 2 =
    // GC + 20 GCA + TAA (read 160-224). The retained G puts exon 2 one base out of the frame it
    // has in the spliced transcript; the ORF is read 0-224, 75 residues.
    let sequence: String = format!("ATG{}{}GGC{}TAA", "GCA".repeat(19), "GCA".repeat(33), "GCA".repeat(20));
    let support_records: Vec<AssembledTranscriptSupportRecord> = vec![AssembledTranscriptSupportRecord {
        assembled_transcript_name: "t".into(),
        sequence: sequence.clone().into_boxed_str(),
        read_names: "r1".into()
    }];
    let alignment_records: Vec<AssembledTranscriptModelAlignmentRecord> = [(0usize, 59usize, "exonic", 1000u32), (60, 159, "intronic", 1060), (160, 224, "exonic", 1160)]
        .into_iter()
        .enumerate()
        .map(|(index, (read_start, read_end, context, position_1))| AssembledTranscriptModelAlignmentRecord {
            assembled_transcript_name: "t".into(), index: index as u32, read_start: read_start as u32, read_end: read_end as u32,
            sequence: sequence[read_start..=read_end].into(), record_type: "base".into(), kind: "match".into(), context: context.into(),
            chromosome_1: "chr1".into(), position_1, operation_1: "I".into(), strand_1: "+".into(),
            chromosome_2: "chr1".into(), position_2: position_1 + (read_end - read_start) as u32, operation_2: "I".into(), strand_2: "+".into(),
            ..Default::default()
        })
        .collect();
    let variant_records: Vec<AssembledTranscriptVariantRecord> = vec![AssembledTranscriptVariantRecord {
        variant_id: 1, assembled_transcript_name: "t".into(), variant_type: "IRT".into(),
        chromosome_1: "chr1".into(), position_1: 1060, operation_1: "I".into(), strand_1: "+".into(),
        chromosome_2: "chr1".into(), position_2: 1159, operation_2: "I".into(), strand_2: "+".into(),
        read_start: 60, read_end: 159,
        ..Default::default()
    }];

    let ts: AssembledTranscriptSet = translate_transcripts(
        &support_records, &alignment_records, &variant_records, &Vec::new(), &Vec::new(), &Vec::new(),
        TranslationStrategy::LongestORF, &HashSet::from_iter(vec!["AUG"]), 1
    ).unwrap();

    let records: Vec<ProteoformRecord> = build_proteoform_records(&ts).collect();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].amino_acid_sequence_length, 75);
    assert_eq!(records[0].mutant_amino_acid_intervals, "20:74");
}

/// F-01. An SNV keeps the reading frame: only its own residue is marked.
#[test]
fn test_mutant_region_not_opened_by_snv() {
    // ATG + 9 GCA + TAA; the SNV turns read 16 (codon 5) from C to G.
    let sequence: &str = "ATGGCAGCAGCAGCAGGAGCAGCAGCAGCATAA";
    let support_records: Vec<AssembledTranscriptSupportRecord> = vec![AssembledTranscriptSupportRecord {
        assembled_transcript_name: "t".into(),
        sequence: sequence.into(),
        read_names: "r1".into()
    }];
    let alignment_records: Vec<AssembledTranscriptModelAlignmentRecord> = [(0usize, 15usize, "match"), (16, 16, "mismatch"), (17, 32, "match")]
        .into_iter()
        .enumerate()
        .map(|(index, (read_start, read_end, kind))| AssembledTranscriptModelAlignmentRecord {
            assembled_transcript_name: "t".into(), index: index as u32, read_start: read_start as u32, read_end: read_end as u32,
            sequence: sequence[read_start..=read_end].into(), record_type: "base".into(), kind: kind.into(), context: "exonic".into(),
            chromosome_1: "chr1".into(), position_1: 1000 + read_start as u32, operation_1: "I".into(), strand_1: "+".into(),
            chromosome_2: "chr1".into(), position_2: 1000 + read_end as u32, operation_2: "I".into(), strand_2: "+".into(),
            ..Default::default()
        })
        .collect();
    let variant_records: Vec<AssembledTranscriptVariantRecord> = vec![AssembledTranscriptVariantRecord {
        variant_id: 1, assembled_transcript_name: "t".into(), variant_type: "SNV".into(),
        chromosome_1: "chr1".into(), position_1: 1016, operation_1: "I".into(), strand_1: "+".into(),
        chromosome_2: "chr1".into(), position_2: 1016, operation_2: "I".into(), strand_2: "+".into(),
        sequence: "G".into(), read_start: 16, read_end: 16,
        ..Default::default()
    }];

    let ts: AssembledTranscriptSet = translate_transcripts(
        &support_records, &alignment_records, &variant_records, &Vec::new(), &Vec::new(), &Vec::new(),
        TranslationStrategy::LongestORF, &HashSet::from_iter(vec!["AUG"]), 1
    ).unwrap();

    let records: Vec<ProteoformRecord> = build_proteoform_records(&ts).collect();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].amino_acid_sequence, "MAAAAGAAAA*");
    assert_eq!(records[0].mutant_amino_acid_intervals, "5:5");
}

/// F-03. A terminal soft-clip insertion has no Base row: the aligner clipped its bases, and
/// call-rna-vars called them. Its bases are matched to it by their own position, so the
/// residue they are read in is marked.
#[test]
fn test_soft_clipped_insertion_marks_its_residue() {
    // The ATG at read 0-2 is the clipped insertion; the first Base row starts at read 3.
    let sequence: &str = "ATGGCAGCAGCATAA";
    let support_records: Vec<AssembledTranscriptSupportRecord> = vec![AssembledTranscriptSupportRecord {
        assembled_transcript_name: "t".into(),
        sequence: sequence.into(),
        read_names: "r1".into()
    }];
    let alignment_records: Vec<AssembledTranscriptModelAlignmentRecord> = vec![AssembledTranscriptModelAlignmentRecord {
        assembled_transcript_name: "t".into(), index: 0, read_start: 3, read_end: 14,
        sequence: sequence[3..15].into(), record_type: "base".into(), kind: "match".into(), context: "exonic".into(),
        chromosome_1: "chr1".into(), position_1: 1001, operation_1: "I".into(), strand_1: "+".into(),
        chromosome_2: "chr1".into(), position_2: 1012, operation_2: "I".into(), strand_2: "+".into(),
        ..Default::default()
    }];
    let variant_records: Vec<AssembledTranscriptVariantRecord> = vec![AssembledTranscriptVariantRecord {
        variant_id: 33, assembled_transcript_name: "t".into(), variant_type: "INS".into(),
        chromosome_1: "chr1".into(), position_1: 1000, operation_1: "D".into(), strand_1: "+".into(),
        chromosome_2: "chr1".into(), position_2: 1001, operation_2: "U".into(), strand_2: "+".into(),
        sequence: "ATG".into(), read_start: 0, read_end: 2,
        ..Default::default()
    }];

    let ts: AssembledTranscriptSet = translate_transcripts(
        &support_records, &alignment_records, &variant_records, &Vec::new(), &Vec::new(), &Vec::new(),
        TranslationStrategy::LongestORF, &HashSet::from_iter(vec!["AUG"]), 1
    ).unwrap();

    let records: Vec<ProteoformRecord> = build_proteoform_records(&ts).collect();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].mutant_amino_acid_intervals, "0:0");
    assert_eq!(records[0].assembled_transcript_variant_ids, "33");
    let nucleotide_records: Vec<NucleotideRecord> = build_nucleotide_records(&ts).collect();
    assert_eq!(
        nucleotide_records.iter().map(|r| r.assembled_transcript_variant_id).collect::<Vec<Option<u32>>>()[..4],
        [Some(33), Some(33), Some(33), None]
    );
    assert_eq!(nucleotide_records[0].assembled_transcript_alignment_index, None);
}

/// F-04. An integrated DNA variant that the DNA variants file does not hold is an error that
/// names it, raised while the set is built; it used to panic only when the tables were written.
#[test]
fn test_build_transcript_set_rejects_unknown_dna_variant() {
    let support_records: Vec<AssembledTranscriptSupportRecord> = vec![AssembledTranscriptSupportRecord {
        assembled_transcript_name: "t".into(),
        sequence: "ATGGCATAA".into(),
        read_names: "r1".into()
    }];
    let alignment_records: Vec<AssembledTranscriptModelAlignmentRecord> = vec![AssembledTranscriptModelAlignmentRecord {
        assembled_transcript_name: "t".into(), index: 0, read_start: 0, read_end: 8,
        sequence: "ATGGCATAA".into(), record_type: "base".into(), kind: "match".into(), context: "exonic".into(),
        chromosome_1: "chr1".into(), position_1: 1000, operation_1: "I".into(), strand_1: "+".into(),
        chromosome_2: "chr1".into(), position_2: 1008, operation_2: "I".into(), strand_2: "+".into(),
        ..Default::default()
    }];
    let dna_variant_records: Vec<DNAVariantRecord> = vec![DNAVariantRecord {
        variant_id: 8, origin: "somatic".into(),
        chromosome_1: "chr1".into(), position_1: 1004, operation_1: "I".into(), strand_1: "+".into(),
        chromosome_2: "chr1".into(), position_2: 1004, operation_2: "I".into(), strand_2: "+".into(),
        sequence: "G".into(),
        ..Default::default()
    }];
    let integrated_variant_records: Vec<IntegratedVariantRecord> = vec![IntegratedVariantRecord {
        assembled_transcript_name: "t".into(),
        reference_gene_name: "".into(),
        reference_transcript_id: "".into(),
        rna_variant_id: 1,
        dna_variant_id: 7,
        distance: 0,
        rna_variant_position: "position_1".into(),
        dna_variant_position: "position_1".into()
    }];

    let result = build_transcript_set(
        &support_records, &alignment_records, &Vec::new(), &dna_variant_records, &integrated_variant_records, &Vec::new()
    );
    match result {
        Err(TranslatorError::UnknownDnaVariant { transcript, dna_variant_id: 7 }) => assert_eq!(transcript.as_ref(), "t"),
        other => panic!("expected UnknownDnaVariant for id 7, got {:?}", other.map(|set| set.len()))
    }
}

/// F-04. A Base row that does not spell the support sequence at its span means the sequence
/// comes from another run, and is an error. A minus-strand row carries its bases in forward
/// orientation and matches the reverse complement; U reads as T.
#[test]
fn test_build_transcript_set_checks_rows_against_the_sequence() {
    let support_records: Vec<AssembledTranscriptSupportRecord> = vec![AssembledTranscriptSupportRecord {
        assembled_transcript_name: "t".into(),
        sequence: "AUGGCAUAA".into(),
        read_names: "r1".into()
    }];
    let row = |sequence: &str, strand: &str| AssembledTranscriptModelAlignmentRecord {
        assembled_transcript_name: "t".into(), index: 4, read_start: 0, read_end: 8,
        sequence: sequence.into(), record_type: "base".into(), kind: "match".into(), context: "exonic".into(),
        chromosome_1: "chr1".into(), position_1: 1000, operation_1: "I".into(), strand_1: strand.into(),
        chromosome_2: "chr1".into(), position_2: 1008, operation_2: "I".into(), strand_2: strand.into(),
        ..Default::default()
    };

    for (sequence, strand) in [("ATGGCATAA", "+"), ("TTATGCCAT", "-")] {
        assert!(
            build_transcript_set(&support_records, &vec![row(sequence, strand)], &Vec::new(), &Vec::new(), &Vec::new(), &Vec::new()).is_ok(),
            "row {} on {} spells the sequence", sequence, strand
        );
    }
    for (sequence, strand) in [("ATGGCATAG", "+"), ("ATGGCATAA", "-"), ("ATGGCATA", "+")] {
        assert!(
            matches!(
                build_transcript_set(&support_records, &vec![row(sequence, strand)], &Vec::new(), &Vec::new(), &Vec::new(), &Vec::new()),
                Err(TranslatorError::SequenceMismatch { index: 4, read_start: 0, read_end: 8, .. })
            ),
            "row {} on {} does not spell the sequence", sequence, strand
        );
    }
}

/// F-05. A base outside ACGTUN (IUPAC R) is read as N: its codon translates to X, as in
/// exacto_core::translate. It used to panic in the nucleotide lookup.
#[test]
fn test_translate_transcripts_reads_iupac_base_as_n() {
    let support_records: Vec<AssembledTranscriptSupportRecord> = vec![AssembledTranscriptSupportRecord {
        assembled_transcript_name: "t".into(),
        sequence: "ATGRCATAA".into(),
        read_names: "r1".into()
    }];
    let alignment_records: Vec<AssembledTranscriptModelAlignmentRecord> = vec![AssembledTranscriptModelAlignmentRecord {
        assembled_transcript_name: "t".into(), index: 0, read_start: 0, read_end: 8,
        sequence: "ATGRCATAA".into(), record_type: "base".into(), kind: "match".into(), context: "exonic".into(),
        chromosome_1: "chr1".into(), position_1: 1000, operation_1: "I".into(), strand_1: "+".into(),
        chromosome_2: "chr1".into(), position_2: 1008, operation_2: "I".into(), strand_2: "+".into(),
        ..Default::default()
    }];

    let ts: AssembledTranscriptSet = translate_transcripts(
        &support_records, &alignment_records, &Vec::new(), &Vec::new(), &Vec::new(), &Vec::new(),
        TranslationStrategy::LongestORF, &HashSet::from_iter(vec!["AUG"]), 1
    ).unwrap();

    let records: Vec<ProteoformRecord> = build_proteoform_records(&ts).collect();
    assert_eq!(records[0].amino_acid_sequence, "MX*");
    let nucleotide_records: Vec<NucleotideRecord> = build_nucleotide_records(&ts).collect();
    assert_eq!(nucleotide_records[3].nucleotide.as_ref(), "N");
    assert_eq!(nucleotide_records[3].amino_acid.as_ref(), "X");
}

/// F-05. translate_fastx_file returns an error naming the file or the record instead of
/// panicking, and translates an IUPAC base to X.
#[test]
fn test_translate_fastx_file_errors_and_iupac() {
    let directory = tempfile::tempdir().unwrap();
    let path = |name: &str| directory.path().join(name).to_str().unwrap().to_string();

    // A missing input file.
    assert!(matches!(
        translate_fastx_file(&path("missing.fasta"), &path("out.fasta.gz"), &path("out.tsv"),
                             TranslationStrategy::LongestORF, HashSet::from_iter(vec!["AUG"]), 1),
        Err(TranslatorError::File { .. })
    ));

    // A FASTQ record whose name is not UTF-8.
    fs::write(path("bad_name.fastq"), b"@r\xff1\nATGGCATAA\n+\nIIIIIIIII\n").unwrap();
    assert!(matches!(
        translate_fastx_file(&path("bad_name.fastq"), &path("out.fasta.gz"), &path("out.tsv"),
                             TranslationStrategy::LongestORF, HashSet::from_iter(vec!["AUG"]), 1),
        Err(TranslatorError::Record { record: 1, .. })
    ));

    // An IUPAC base in the second record.
    fs::write(path("iupac.fasta"), ">r1\nATGGCATAA\n>r2\nATGRCATAA\n").unwrap();
    translate_fastx_file(&path("iupac.fasta"), &path("out.fasta.gz"), &path("out.tsv"),
                         TranslationStrategy::LongestORF, HashSet::from_iter(vec!["AUG"]), 1).unwrap();
    let table: String = fs::read_to_string(path("out.tsv")).unwrap();
    let peptides: Vec<&str> = table.lines().skip(1).map(|line| line.split('\t').nth(3).unwrap()).collect();
    assert_eq!(peptides, vec!["MA*", "MX*"]);
}
