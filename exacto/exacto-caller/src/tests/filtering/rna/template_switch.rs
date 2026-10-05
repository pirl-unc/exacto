use bimap::BiMap;
use exacto_core::prelude::*;
use std::collections::HashSet;
use tempfile::TempDir;

use super::*;


/// Forward/forward junction with a repeat planted on both flanks.
///
/// chrA (leaves at 30, '+'): positions 21-30 = "GGGGGGGCAT", 31-40 = "ACGTATTTTT".
/// chrB (enters at 51, '+'): positions 41-50 = "TTTTTTTCAT", 51-60 = "ACGTACCCCC".
/// h_3prime = LCP("ACGTATTTTT", "ACGTACCCCC") = 5; h_5prime = common suffix of
/// "GGGGGGGCAT" and "TTTTTTTCAT" = 3 ("CAT"). Positions are 1-based inclusive, the
/// `FastaMap::get_sequence` convention: an off-by-one in a flank window shifts a length.
#[test]
fn compute_junction_homology_forward_forward_junction_measures_planted_repeat() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}", "T".repeat(20), "GGGGGGGCAT", "ACGTATTTTT", "G".repeat(40)).into()),
        ("chrB".into(), format!("{}{}{}{}", "C".repeat(40), "TTTTTTTCAT", "ACGTACCCCC", "A".repeat(40)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    let graph_operation: GraphOperation = GraphOperation::new(
        0, 30, Strand::Forward, GraphOperationType::Downstream,
        1, 51, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    let homology: JunctionHomology = compute_junction_homology(
        &graph_operation,
        &chromosome_names_map,
        &fasta_map,
        10
    ).unwrap();

    assert_eq!(homology.h_3prime, 5, "prefix of the abandoned continuation vs the entered head");
    assert_eq!(homology.h_5prime, 3, "suffix of the aligned tail vs the entering side's preceding bases");
    assert_eq!(homology.total(), 8);
    assert_eq!(homology.canonical_splice, false);
}

/// The same junction geometry with the ENTERING side on the minus strand: chrB is the
/// reverse complement of the contig above, so the read enters at its genomic END. The
/// homology must be identical to the forward/forward case; the biology does not care how
/// the aligner spelled the second side.
#[test]
fn compute_junction_homology_reverse_entering_side_matches_forward_measurement() {
    // chrB is 100 bases, so the entry point maps from position 51 to 100 - 51 + 1 = 50.
    let chr_b_forward: String = format!("{}{}{}{}", "C".repeat(40), "TTTTTTTCAT", "ACGTACCCCC", "A".repeat(40));
    let chr_b_reverse: String = chr_b_forward.chars().rev().map(|base| match base {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        other => other
    }).collect();

    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}", "T".repeat(20), "GGGGGGGCAT", "ACGTATTTTT", "G".repeat(40)).into()),
        ("chrB".into(), chr_b_reverse.into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    // Reverse-strand read entering at genomic position 50 of the reverse-complemented
    // contig, junction on its genomic 3' side: Downstream + Reverse = Enters.
    let graph_operation: GraphOperation = GraphOperation::new(
        0, 30, Strand::Forward, GraphOperationType::Downstream,
        1, 50, Strand::Reverse, GraphOperationType::Downstream,
        "".into(), VariantType::Breakpoint
    );
    let homology: JunctionHomology = compute_junction_homology(
        &graph_operation,
        &chromosome_names_map,
        &fasta_map,
        10
    ).unwrap();

    assert_eq!((homology.h_5prime, homology.h_3prime), (3, 5));
    assert_eq!(homology.canonical_splice, false);
}

/// GT-AG: the leaving side continues "GT", the entering side is preceded by "AG".
#[test]
fn compute_junction_homology_detects_canonical_splice_motif() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        // Positions 31-32 = "GT".
        ("chrA".into(), format!("{}{}", "A".repeat(30), "GTCCCCCCCC").into()),
        // Positions 49-50 = "AG"; the read enters at 51.
        ("chrB".into(), format!("{}{}{}", "C".repeat(48), "AG", "T".repeat(10)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    let graph_operation: GraphOperation = GraphOperation::new(
        0, 30, Strand::Forward, GraphOperationType::Downstream,
        1, 51, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    let homology: JunctionHomology = compute_junction_homology(
        &graph_operation,
        &chromosome_names_map,
        &fasta_map,
        10
    ).unwrap();

    assert_eq!(homology.canonical_splice, true);
    // The motif is not homology: nothing is shared across the junction.
    assert_eq!(homology.total(), 0);
}

/// N never extends homology: an unknown base is not evidence of a shared motif.
#[test]
fn compute_junction_homology_stops_at_n_bases() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}", "T".repeat(30), "ACNTAGGGGG").into()),
        ("chrB".into(), format!("{}{}", "C".repeat(50), "ACNTACCCCC").into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    let graph_operation: GraphOperation = GraphOperation::new(
        0, 30, Strand::Forward, GraphOperationType::Downstream,
        1, 51, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    let homology: JunctionHomology = compute_junction_homology(
        &graph_operation,
        &chromosome_names_map,
        &fasta_map,
        10
    ).unwrap();

    assert_eq!(homology.h_3prime, 2, "match must stop at the N, not skip it");
    assert_eq!(homology.h_5prime, 0);
}

/// Contig-edge clamping shortens the window instead of failing.
#[test]
fn compute_junction_homology_shortens_window_at_contig_edge() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), "ACGTACGT".into()),
        ("chrB".into(), format!("{}{}", "T".repeat(50), "ACGTAAAAAA").into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    // Leaves chrA at its last base: the continuation window is entirely off-contig, and
    // the aligned-tail window is clipped to the 8 bases the contig has.
    let graph_operation: GraphOperation = GraphOperation::new(
        0, 8, Strand::Forward, GraphOperationType::Downstream,
        1, 51, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    let homology: JunctionHomology = compute_junction_homology(
        &graph_operation,
        &chromosome_names_map,
        &fasta_map,
        10
    ).unwrap();

    assert_eq!(homology.h_3prime, 0, "no continuation exists past the contig end");
    // The clipped tail "ACGTACGT" still shares its final T with the entering side's "TTTTTTTTTT".
    assert_eq!(homology.h_5prime, 1);
}

/// Two sides with the same traversal role cannot be oriented; the fold-back shapes land
/// here by construction and are the shape predicate's business instead.
#[test]
fn compute_junction_homology_returns_none_for_same_role_sides() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), "ACGT".repeat(30).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);

    // Forward + Downstream on both sides: both leave, neither enters.
    let graph_operation: GraphOperation = GraphOperation::new(
        0, 30, Strand::Forward, GraphOperationType::Downstream,
        0, 80, Strand::Forward, GraphOperationType::Downstream,
        "".into(), VariantType::Breakpoint
    );

    assert!(compute_junction_homology(&graph_operation, &chromosome_names_map, &fasta_map, 10).is_none());
}

#[test]
fn is_foldback_shape_returns_matches() {
    // Same chromosome, same breakpoint operation on both sides, opposite strands.
    let foldback: GraphOperation = GraphOperation::new(
        0, 1000, Strand::Forward, GraphOperationType::Downstream,
        0, 1400, Strand::Reverse, GraphOperationType::Downstream,
        "".into(), VariantType::Breakpoint
    );
    assert_eq!(is_foldback_shape(&foldback), true);

    // Same strand: a plain deletion-like junction.
    let linear: GraphOperation = GraphOperation::new(
        0, 1000, Strand::Forward, GraphOperationType::Downstream,
        0, 1400, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    assert_eq!(is_foldback_shape(&linear), false);

    // Opposite strands but different chromosomes: a translocation, not a fold.
    let translocation: GraphOperation = GraphOperation::new(
        0, 1000, Strand::Forward, GraphOperationType::Downstream,
        1, 1400, Strand::Reverse, GraphOperationType::Downstream,
        "".into(), VariantType::Translocation
    );
    assert_eq!(is_foldback_shape(&translocation), false);
}

#[test]
fn is_junction_variant_type_returns_matches() {
    assert_eq!(is_junction_variant_type(&VariantType::Breakpoint), true);
    assert_eq!(is_junction_variant_type(&VariantType::Translocation), true);
    assert_eq!(is_junction_variant_type(&VariantType::FusionGene), true);
    assert_eq!(is_junction_variant_type(&VariantType::CircularRNA), true);
    assert_eq!(is_junction_variant_type(&VariantType::NonCanonicalSplicing), true);
    assert_eq!(is_junction_variant_type(&VariantType::Insertion), false);
    assert_eq!(is_junction_variant_type(&VariantType::Deletion), false);
    assert_eq!(is_junction_variant_type(&VariantType::SingleNucleotideVariant), false);
    assert_eq!(is_junction_variant_type(&VariantType::MultiNucleotideVariant), false);
}

/// The LEAVING side on the minus strand: chrA is the reverse complement of the planted
/// contig, so the read leaves at genomic position 80 - 30 + 1 = 51 and Upstream + Reverse
/// = Leaves. The other reverse branch of the flank fetch; same 3 and 5 expected.
#[test]
fn compute_junction_homology_reverse_leaving_side_matches_forward_measurement() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let chr_a_forward: String = format!("{}{}{}{}", "T".repeat(20), "GGGGGGGCAT", "ACGTATTTTT", "G".repeat(40));
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), reverse_complement(&chr_a_forward)),
        ("chrB".into(), format!("{}{}{}{}", "C".repeat(40), "TTTTTTTCAT", "ACGTACCCCC", "A".repeat(40)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    let graph_operation: GraphOperation = GraphOperation::new(
        0, 51, Strand::Reverse, GraphOperationType::Upstream,
        1, 51, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    let homology: JunctionHomology = compute_junction_homology(
        &graph_operation,
        &chromosome_names_map,
        &fasta_map,
        10
    ).unwrap();

    assert_eq!((homology.h_5prime, homology.h_3prime), (3, 5));
    assert_eq!(homology.canonical_splice, false);
}

/// Both sides on the minus strand: every flank is fetched reverse-complemented.
#[test]
fn compute_junction_homology_both_sides_reverse_matches_forward_measurement() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let chr_a_forward: String = format!("{}{}{}{}", "T".repeat(20), "GGGGGGGCAT", "ACGTATTTTT", "G".repeat(40));
    let chr_b_forward: String = format!("{}{}{}{}", "C".repeat(40), "TTTTTTTCAT", "ACGTACCCCC", "A".repeat(40));
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), reverse_complement(&chr_a_forward)),
        ("chrB".into(), reverse_complement(&chr_b_forward))
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    // chrA leaves at 80 - 30 + 1 = 51; chrB enters at 100 - 51 + 1 = 50.
    let graph_operation: GraphOperation = GraphOperation::new(
        0, 51, Strand::Reverse, GraphOperationType::Upstream,
        1, 50, Strand::Reverse, GraphOperationType::Downstream,
        "".into(), VariantType::Breakpoint
    );
    let homology: JunctionHomology = compute_junction_homology(
        &graph_operation,
        &chromosome_names_map,
        &fasta_map,
        10
    ).unwrap();

    assert_eq!((homology.h_5prime, homology.h_3prime), (3, 5));
    assert_eq!(homology.canonical_splice, false);
}

/// `GraphOperation::new` orders sides by position on one chromosome, so an operation
/// spelled entering-side-first comes back leaving-side-first. Roles are assigned per side,
/// so the measurement does not depend on which slot each side was given in.
#[test]
fn compute_junction_homology_ignores_side_order() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    // 21-30 GGGGGGGCAT, 31-40 ACGTATTTTT, 41-50 TTTTTTTCAT, 51-60 ACGTACCCCC: a deletion-shaped
    // junction leaving at 30 and entering at 51 on one contig.
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}{}{}", "T".repeat(20), "GGGGGGGCAT", "ACGTATTTTT", "TTTTTTTCAT", "ACGTACCCCC", "A".repeat(40)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);

    let graph_operation: GraphOperation = GraphOperation::new(
        0, 51, Strand::Forward, GraphOperationType::Upstream,
        0, 30, Strand::Forward, GraphOperationType::Downstream,
        "".into(), VariantType::Breakpoint
    );
    assert_eq!(graph_operation.get_position_1(), 30, "standardization put the leaving side first");

    let homology: JunctionHomology = compute_junction_homology(
        &graph_operation,
        &chromosome_names_map,
        &fasta_map,
        10
    ).unwrap();

    assert_eq!((homology.h_5prime, homology.h_3prime), (3, 5));
}

/// A repeat longer than the flank reports exactly the flank; a flank longer than the
/// repeat reports the repeat. The flank is a ceiling on the measurement, never a floor.
#[test]
fn compute_junction_homology_caps_each_side_at_flank() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    // chrA: G at 16-30, C at 31-45; leaves at 30.  chrB: G at 36-50, C at 51-65; enters at 51.
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}", "A".repeat(15), "G".repeat(15), "C".repeat(15), "A".repeat(15)).into()),
        ("chrB".into(), format!("{}{}{}{}", "T".repeat(35), "G".repeat(15), "C".repeat(15), "T".repeat(15)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    let graph_operation: GraphOperation = GraphOperation::new(
        0, 30, Strand::Forward, GraphOperationType::Downstream,
        1, 51, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );

    let homology: JunctionHomology = compute_junction_homology(&graph_operation, &chromosome_names_map, &fasta_map, 5).unwrap();
    assert_eq!((homology.h_5prime, homology.h_3prime), (5, 5));

    let homology: JunctionHomology = compute_junction_homology(&graph_operation, &chromosome_names_map, &fasta_map, 10).unwrap();
    assert_eq!((homology.h_5prime, homology.h_3prime), (10, 10));
    assert_eq!(homology.total(), 20);

    // The planted repeat is 15 on each side; a 20-base flank runs into the mismatching filler.
    let homology: JunctionHomology = compute_junction_homology(&graph_operation, &chromosome_names_map, &fasta_map, 20).unwrap();
    assert_eq!((homology.h_5prime, homology.h_3prime), (15, 15));
}

#[test]
fn compute_junction_homology_returns_none_for_non_junction_variant_types() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}", "T".repeat(20), "GGGGGGGCAT", "ACGTATTTTT", "G".repeat(40)).into()),
        ("chrB".into(), format!("{}{}{}{}", "C".repeat(40), "TTTTTTTCAT", "ACGTACCCCC", "A".repeat(40)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    // Same geometry as the planted-repeat case, so the only reason for None is the type gate.
    for variant_type in [VariantType::Insertion, VariantType::Deletion, VariantType::SingleNucleotideVariant] {
        let graph_operation: GraphOperation = GraphOperation::new(
            0, 30, Strand::Forward, GraphOperationType::Downstream,
            1, 51, Strand::Forward, GraphOperationType::Upstream,
            "".into(), variant_type.clone()
        );
        assert!(
            compute_junction_homology(&graph_operation, &chromosome_names_map, &fasta_map, 10).is_none(),
            "{variant_type:?} is not a junction"
        );
    }
}

#[test]
fn compute_junction_homology_returns_none_for_zero_flank() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}", "T".repeat(20), "GGGGGGGCAT", "ACGTATTTTT", "G".repeat(40)).into()),
        ("chrB".into(), format!("{}{}{}{}", "C".repeat(40), "TTTTTTTCAT", "ACGTACCCCC", "A".repeat(40)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    let graph_operation: GraphOperation = GraphOperation::new(
        0, 30, Strand::Forward, GraphOperationType::Downstream,
        1, 51, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );

    assert!(compute_junction_homology(&graph_operation, &chromosome_names_map, &fasta_map, 0).is_none());
    assert!(compute_junction_homology(&graph_operation, &chromosome_names_map, &fasta_map, 1).is_some());
}

/// Position 0 is outside 1-based coordinates: there is no base to leave from or enter at.
#[test]
fn compute_junction_homology_returns_none_at_position_zero() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}", "T".repeat(20), "GGGGGGGCAT", "ACGTATTTTT", "G".repeat(40)).into()),
        ("chrB".into(), format!("{}{}{}{}", "C".repeat(40), "TTTTTTTCAT", "ACGTACCCCC", "A".repeat(40)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    let leaving_at_zero: GraphOperation = GraphOperation::new(
        0, 0, Strand::Forward, GraphOperationType::Downstream,
        1, 51, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    assert!(compute_junction_homology(&leaving_at_zero, &chromosome_names_map, &fasta_map, 10).is_none());

    let entering_at_zero: GraphOperation = GraphOperation::new(
        0, 30, Strand::Forward, GraphOperationType::Downstream,
        1, 0, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    assert!(compute_junction_homology(&entering_at_zero, &chromosome_names_map, &fasta_map, 10).is_none());
}

/// hg38 is soft-masked: repeats are lowercase. The comparison must be case-insensitive.
#[test]
fn compute_junction_homology_measures_soft_masked_repeat() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}", "T".repeat(20), "gggggggcat", "acgtattttt", "G".repeat(40)).into()),
        ("chrB".into(), format!("{}{}{}{}", "C".repeat(40), "tttttttcat", "acgtaccccc", "A".repeat(40)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    let graph_operation: GraphOperation = GraphOperation::new(
        0, 30, Strand::Forward, GraphOperationType::Downstream,
        1, 51, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    let homology: JunctionHomology = compute_junction_homology(
        &graph_operation,
        &chromosome_names_map,
        &fasta_map,
        10
    ).unwrap();

    assert_eq!((homology.h_5prime, homology.h_3prime), (3, 5));
}

/// The 5' suffix walk stops at an N just as the 3' prefix walk does. Without the N the
/// two 5' flanks would be identical over all 10 bases.
#[test]
fn compute_junction_homology_stops_at_n_bases_on_five_prime_side() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}", "T".repeat(20), "GGGGGGNCAT", "ACGTATTTTT", "G".repeat(40)).into()),
        ("chrB".into(), format!("{}{}{}{}", "C".repeat(40), "GGGGGGNCAT", "ACGTACCCCC", "A".repeat(40)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    let graph_operation: GraphOperation = GraphOperation::new(
        0, 30, Strand::Forward, GraphOperationType::Downstream,
        1, 51, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    let homology: JunctionHomology = compute_junction_homology(
        &graph_operation,
        &chromosome_names_map,
        &fasta_map,
        10
    ).unwrap();

    assert_eq!(homology.h_5prime, 3, "suffix match must stop at the N, not skip it");
    assert_eq!(homology.h_3prime, 5);
}

/// A minus-strand splice shows CT...AC on the genome. The motif check runs on
/// read-oriented flanks, so it must still read as GT-AG.
#[test]
fn compute_junction_homology_detects_canonical_splice_motif_on_reverse_strand() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let chr_a_forward: String = format!("{}{}", "A".repeat(30), "GTCCCCCCCC");
    let chr_b_forward: String = format!("{}{}{}", "C".repeat(48), "AG", "T".repeat(10));
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), reverse_complement(&chr_a_forward)),
        ("chrB".into(), reverse_complement(&chr_b_forward))
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    // chrA (40 bases) leaves at 40 - 30 + 1 = 11; chrB (60 bases) enters at 60 - 51 + 1 = 10.
    let graph_operation: GraphOperation = GraphOperation::new(
        0, 11, Strand::Reverse, GraphOperationType::Upstream,
        1, 10, Strand::Reverse, GraphOperationType::Downstream,
        "".into(), VariantType::Breakpoint
    );
    let homology: JunctionHomology = compute_junction_homology(
        &graph_operation,
        &chromosome_names_map,
        &fasta_map,
        10
    ).unwrap();

    assert_eq!(homology.canonical_splice, true);
    assert_eq!(homology.total(), 0);
}

/// GT without AG, and AG without GT, are not a splice.
#[test]
fn compute_junction_homology_requires_both_halves_of_the_splice_motif() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}", "A".repeat(30), "GTCCCCCCCC").into()),           // GT at 31-32
        ("chrB".into(), format!("{}{}", "C".repeat(50), "T".repeat(10)).into()),          // no AG before 51
        ("chrC".into(), "A".repeat(40).into()),                                           // no GT after 30
        ("chrD".into(), format!("{}{}{}", "C".repeat(48), "AG", "T".repeat(10)).into())   // AG at 49-50
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);
    chromosome_names_map.insert("chrC".into(), 2);
    chromosome_names_map.insert("chrD".into(), 3);

    let donor_only: GraphOperation = GraphOperation::new(
        0, 30, Strand::Forward, GraphOperationType::Downstream,
        1, 51, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    let homology: JunctionHomology = compute_junction_homology(&donor_only, &chromosome_names_map, &fasta_map, 10).unwrap();
    assert_eq!(homology.canonical_splice, false);

    let acceptor_only: GraphOperation = GraphOperation::new(
        2, 30, Strand::Forward, GraphOperationType::Downstream,
        3, 51, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    let homology: JunctionHomology = compute_junction_homology(&acceptor_only, &chromosome_names_map, &fasta_map, 10).unwrap();
    assert_eq!(homology.canonical_splice, false);
}

/// Leaving 3 bases into a contig: the aligned-tail window is clipped to those 3 bases and
/// still measures the repeat they carry.
#[test]
fn compute_junction_homology_clamps_leaving_window_at_contig_start() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}", "CAT", "ACGTATTTTT", "G".repeat(40)).into()),
        ("chrB".into(), format!("{}{}{}{}", "C".repeat(40), "TTTTTTTCAT", "ACGTACCCCC", "A".repeat(40)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    let graph_operation: GraphOperation = GraphOperation::new(
        0, 3, Strand::Forward, GraphOperationType::Downstream,
        1, 51, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    let homology: JunctionHomology = compute_junction_homology(
        &graph_operation,
        &chromosome_names_map,
        &fasta_map,
        10
    ).unwrap();

    assert_eq!((homology.h_5prime, homology.h_3prime), (3, 5));
}

/// Entering at a contig's first base: nothing precedes it on that template, so the 5'
/// window is empty and only the 3' side can measure.
#[test]
fn compute_junction_homology_clamps_entering_window_at_contig_start() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}", "T".repeat(20), "GGGGGGGCAT", "ACGTATTTTT", "G".repeat(40)).into()),
        ("chrB".into(), format!("{}{}", "ACGTACCCCC", "A".repeat(40)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    let graph_operation: GraphOperation = GraphOperation::new(
        0, 30, Strand::Forward, GraphOperationType::Downstream,
        1, 1, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    let homology: JunctionHomology = compute_junction_homology(
        &graph_operation,
        &chromosome_names_map,
        &fasta_map,
        10
    ).unwrap();

    assert_eq!((homology.h_5prime, homology.h_3prime), (0, 5));
}

/// A hairpin planted on one contig: eight A's, a three-base loop, eight T's. Positions
/// 51-58 = A, 59-61 = C, 62-69 = T, G filler either side so the stem cannot extend.
#[test]
fn compute_foldback_stem_measures_planted_hairpin() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}{}", "G".repeat(50), "A".repeat(8), "CCC", "T".repeat(8), "G".repeat(50)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);

    // Leaves the plus strand at 51, re-enters the minus strand at 69: fold-back shaped.
    let graph_operation: GraphOperation = GraphOperation::new(
        0, 51, Strand::Forward, GraphOperationType::Downstream,
        0, 69, Strand::Reverse, GraphOperationType::Downstream,
        "".into(), VariantType::Breakpoint
    );
    let stem: FoldbackStem = compute_foldback_stem(
        &graph_operation,
        &chromosome_names_map,
        &fasta_map,
        10,
        5
    ).unwrap();

    assert_eq!(stem, FoldbackStem { span: 18, stem: 8, loop_length: 3 });
}

/// The loop cap is a hard gate: a three-base loop is invisible when at most two are allowed,
/// and nothing else in the window pairs.
#[test]
fn compute_foldback_stem_caps_the_loop_at_max_loop_length() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}{}", "G".repeat(50), "A".repeat(8), "CCC", "T".repeat(8), "G".repeat(50)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);

    let graph_operation: GraphOperation = GraphOperation::new(
        0, 51, Strand::Forward, GraphOperationType::Downstream,
        0, 69, Strand::Reverse, GraphOperationType::Downstream,
        "".into(), VariantType::Breakpoint
    );

    let stem: FoldbackStem = compute_foldback_stem(&graph_operation, &chromosome_names_map, &fasta_map, 10, 2).unwrap();
    assert_eq!(stem, FoldbackStem { span: 18, stem: 0, loop_length: 0 });

    // The cap is inclusive: a loop of exactly three is allowed at three.
    let stem: FoldbackStem = compute_foldback_stem(&graph_operation, &chromosome_names_map, &fasta_map, 10, 3).unwrap();
    assert_eq!(stem, FoldbackStem { span: 18, stem: 8, loop_length: 3 });
}

/// An N in the stem ends the pairing walk: the four A's nearest the loop pair, the N does not.
#[test]
fn compute_foldback_stem_stops_at_n_bases() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}{}", "G".repeat(50), "AAANAAAA", "CCC", "T".repeat(8), "G".repeat(50)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);

    let graph_operation: GraphOperation = GraphOperation::new(
        0, 51, Strand::Forward, GraphOperationType::Downstream,
        0, 69, Strand::Reverse, GraphOperationType::Downstream,
        "".into(), VariantType::Breakpoint
    );
    let stem: FoldbackStem = compute_foldback_stem(
        &graph_operation,
        &chromosome_names_map,
        &fasta_map,
        10,
        5
    ).unwrap();

    assert_eq!(stem, FoldbackStem { span: 18, stem: 4, loop_length: 3 });
}

#[test]
fn compute_foldback_stem_returns_none_for_non_foldback_shape() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        ("chrA".into(), format!("{}{}{}{}{}", "G".repeat(50), "A".repeat(8), "CCC", "T".repeat(8), "G".repeat(50)).into()),
        ("chrB".into(), format!("{}{}{}{}{}", "G".repeat(50), "A".repeat(8), "CCC", "T".repeat(8), "G".repeat(50)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    // The same hairpin, but the read stays on one strand: a deletion-shaped junction.
    let same_strand: GraphOperation = GraphOperation::new(
        0, 51, Strand::Forward, GraphOperationType::Downstream,
        0, 69, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    assert!(compute_foldback_stem(&same_strand, &chromosome_names_map, &fasta_map, 10, 5).is_none());

    // Opposite strands on different chromosomes: a translocation, not a fold.
    let across_chromosomes: GraphOperation = GraphOperation::new(
        0, 51, Strand::Forward, GraphOperationType::Downstream,
        1, 69, Strand::Reverse, GraphOperationType::Downstream,
        "".into(), VariantType::Translocation
    );
    assert!(compute_foldback_stem(&across_chromosomes, &chromosome_names_map, &fasta_map, 10, 5).is_none());
}

#[test]
fn classify_template_switch_flags_hairpin_foldback() {
    let evidence: TemplateSwitchEvidence = TemplateSwitchEvidence {
        junction_homology: None,
        foldback_stem: Some(FoldbackStem { span: 100, stem: 8, loop_length: 3 }),
        at_exon_boundaries: false,
        has_exit_junction: false,
        dispersion_beyond_homology: None,
        num_reads: 3
    };

    assert_eq!(
        classify_template_switch(&evidence, 8, 4, 10, 1000, 6),
        TemplateSwitchVerdict::Flagged(TemplateSwitchReason::Foldback)
    );
}

/// A fold-back is judged on its stem, never on homology: the two flanks are one locus read
/// both ways. An unresolved fold-back is Clear even with a repeat panel A would flag.
#[test]
fn classify_template_switch_clears_unresolved_foldback_without_consulting_homology() {
    let evidence: TemplateSwitchEvidence = TemplateSwitchEvidence {
        junction_homology: Some(JunctionHomology { h_5prime: 10, h_3prime: 10, canonical_splice: false }),
        foldback_stem: Some(FoldbackStem { span: 100, stem: 2, loop_length: 0 }),
        at_exon_boundaries: false,
        has_exit_junction: false,
        dispersion_beyond_homology: None,
        num_reads: 3
    };

    assert_eq!(classify_template_switch(&evidence, 8, 4, 10, 1000, 6), TemplateSwitchVerdict::Clear);
}

#[test]
fn classify_template_switch_clears_foldback_on_exon_boundaries() {
    let evidence: TemplateSwitchEvidence = TemplateSwitchEvidence {
        junction_homology: None,
        foldback_stem: Some(FoldbackStem { span: 100, stem: 8, loop_length: 3 }),
        at_exon_boundaries: true,
        has_exit_junction: false,
        dispersion_beyond_homology: None,
        num_reads: 3
    };

    assert_eq!(classify_template_switch(&evidence, 8, 4, 10, 1000, 6), TemplateSwitchVerdict::Clear);
}

#[test]
fn classify_template_switch_clears_foldback_with_exit_junction() {
    let evidence: TemplateSwitchEvidence = TemplateSwitchEvidence {
        junction_homology: None,
        foldback_stem: Some(FoldbackStem { span: 100, stem: 8, loop_length: 3 }),
        at_exon_boundaries: false,
        has_exit_junction: true,
        dispersion_beyond_homology: None,
        num_reads: 3
    };

    assert_eq!(classify_template_switch(&evidence, 8, 4, 10, 1000, 6), TemplateSwitchVerdict::Clear);
}

/// Beyond `foldback_max_distance` the fold-back panel steps aside and homology decides.
/// The distance gate is inclusive.
#[test]
fn classify_template_switch_hands_distant_foldback_to_homology() {
    let mut evidence: TemplateSwitchEvidence = TemplateSwitchEvidence {
        junction_homology: Some(JunctionHomology { h_5prime: 5, h_3prime: 5, canonical_splice: false }),
        foldback_stem: Some(FoldbackStem { span: 1001, stem: 8, loop_length: 3 }),
        at_exon_boundaries: false,
        has_exit_junction: false,
        dispersion_beyond_homology: None,
        num_reads: 3
    };
    assert_eq!(
        classify_template_switch(&evidence, 8, 4, 10, 1000, 6),
        TemplateSwitchVerdict::Flagged(TemplateSwitchReason::StrongHomology)
    );

    evidence.foldback_stem = Some(FoldbackStem { span: 1000, stem: 8, loop_length: 3 });
    assert_eq!(
        classify_template_switch(&evidence, 8, 4, 10, 1000, 6),
        TemplateSwitchVerdict::Flagged(TemplateSwitchReason::Foldback)
    );
}

#[test]
fn classify_template_switch_foldback_min_stem_is_inclusive() {
    let mut evidence: TemplateSwitchEvidence = TemplateSwitchEvidence {
        junction_homology: None,
        foldback_stem: Some(FoldbackStem { span: 100, stem: 6, loop_length: 0 }),
        at_exon_boundaries: false,
        has_exit_junction: false,
        dispersion_beyond_homology: None,
        num_reads: 3
    };
    assert_eq!(
        classify_template_switch(&evidence, 8, 4, 10, 1000, 6),
        TemplateSwitchVerdict::Flagged(TemplateSwitchReason::Foldback)
    );

    evidence.foldback_stem = Some(FoldbackStem { span: 100, stem: 5, loop_length: 0 });
    assert_eq!(classify_template_switch(&evidence, 8, 4, 10, 1000, 6), TemplateSwitchVerdict::Clear);
}

#[test]
fn classify_template_switch_returns_not_assessed_without_homology() {
    let evidence: TemplateSwitchEvidence = TemplateSwitchEvidence {
        junction_homology: None,
        foldback_stem: None,
        at_exon_boundaries: false,
        has_exit_junction: false,
        dispersion_beyond_homology: Some(50),
        num_reads: 3
    };

    assert_eq!(classify_template_switch(&evidence, 8, 4, 10, 1000, 6), TemplateSwitchVerdict::NotAssessed);
}

#[test]
fn classify_template_switch_clears_canonical_splice_despite_strong_homology() {
    let evidence: TemplateSwitchEvidence = TemplateSwitchEvidence {
        junction_homology: Some(JunctionHomology { h_5prime: 10, h_3prime: 10, canonical_splice: true }),
        foldback_stem: None,
        at_exon_boundaries: false,
        has_exit_junction: false,
        dispersion_beyond_homology: Some(50),
        num_reads: 3
    };

    assert_eq!(classify_template_switch(&evidence, 8, 4, 10, 1000, 6), TemplateSwitchVerdict::Clear);
}

#[test]
fn classify_template_switch_clears_exon_boundary_junction_despite_strong_homology() {
    let evidence: TemplateSwitchEvidence = TemplateSwitchEvidence {
        junction_homology: Some(JunctionHomology { h_5prime: 10, h_3prime: 10, canonical_splice: false }),
        foldback_stem: None,
        at_exon_boundaries: true,
        has_exit_junction: false,
        dispersion_beyond_homology: Some(50),
        num_reads: 3
    };

    assert_eq!(classify_template_switch(&evidence, 8, 4, 10, 1000, 6), TemplateSwitchVerdict::Clear);
}

/// `min_homology` is inclusive; one base short falls into the soft band.
#[test]
fn classify_template_switch_min_homology_is_inclusive() {
    let mut evidence: TemplateSwitchEvidence = TemplateSwitchEvidence {
        junction_homology: Some(JunctionHomology { h_5prime: 4, h_3prime: 4, canonical_splice: false }),
        foldback_stem: None,
        at_exon_boundaries: false,
        has_exit_junction: false,
        dispersion_beyond_homology: None,
        num_reads: 3
    };
    assert_eq!(
        classify_template_switch(&evidence, 8, 4, 10, 1000, 6),
        TemplateSwitchVerdict::Flagged(TemplateSwitchReason::StrongHomology)
    );

    evidence.junction_homology = Some(JunctionHomology { h_5prime: 4, h_3prime: 3, canonical_splice: false });
    assert_eq!(classify_template_switch(&evidence, 8, 4, 10, 1000, 6), TemplateSwitchVerdict::Clear);
}

/// The soft band flags only with measured dispersion strictly beyond the maximum. No
/// measurement, dispersion at the maximum, or homology below the soft floor all clear.
#[test]
fn classify_template_switch_flags_soft_homology_only_with_dispersion_beyond_maximum() {
    let mut evidence: TemplateSwitchEvidence = TemplateSwitchEvidence {
        junction_homology: Some(JunctionHomology { h_5prime: 2, h_3prime: 2, canonical_splice: false }),
        foldback_stem: None,
        at_exon_boundaries: false,
        has_exit_junction: false,
        dispersion_beyond_homology: Some(11),
        num_reads: 3
    };
    assert_eq!(
        classify_template_switch(&evidence, 8, 4, 10, 1000, 6),
        TemplateSwitchVerdict::Flagged(TemplateSwitchReason::SoftHomologyWithDispersion)
    );

    evidence.dispersion_beyond_homology = Some(10);
    assert_eq!(classify_template_switch(&evidence, 8, 4, 10, 1000, 6), TemplateSwitchVerdict::Clear);

    evidence.dispersion_beyond_homology = None;
    assert_eq!(classify_template_switch(&evidence, 8, 4, 10, 1000, 6), TemplateSwitchVerdict::Clear);

    evidence.dispersion_beyond_homology = Some(11);
    evidence.junction_homology = Some(JunctionHomology { h_5prime: 2, h_3prime: 1, canonical_splice: false });
    assert_eq!(classify_template_switch(&evidence, 8, 4, 10, 1000, 6), TemplateSwitchVerdict::Clear);
}

#[test]
fn is_template_switch_returns_matches() {
    use TemplateSwitchReason::*;
    use TemplateSwitchVerdict::*;

    // Only the pooled breakend family is ever suppressed.
    assert_eq!(is_template_switch(Flagged(StrongHomology), &VariantType::Breakpoint), true);
    assert_eq!(is_template_switch(Flagged(Foldback), &VariantType::Translocation), true);
    assert_eq!(is_template_switch(Flagged(SoftHomologyWithDispersion), &VariantType::Breakpoint), true);
    assert_eq!(is_template_switch(Flagged(StrongHomology), &VariantType::FusionGene), false);
    assert_eq!(is_template_switch(Flagged(StrongHomology), &VariantType::CircularRNA), false);
    assert_eq!(is_template_switch(Flagged(StrongHomology), &VariantType::NonCanonicalSplicing), false);

    // Hub evidence is reported, never suppressed.
    assert_eq!(is_template_switch(Flagged(SoftHomologyAtHub), &VariantType::Breakpoint), false);

    assert_eq!(is_template_switch(Clear, &VariantType::Breakpoint), false);
    assert_eq!(is_template_switch(NotAssessed, &VariantType::Breakpoint), false);
}

/// Two identical members fix the consensus; two more spell the breakpoint a few bases
/// away on each side. The spread is the larger of the two sides.
#[test]
fn pooled_breakpoint_spread_returns_the_larger_side_spread() {
    let mut variant_records: HashSet<VariantRecord> = HashSet::new();
    for (read_id, position_1, position_2) in [(1, 1000, 2000), (2, 1000, 2000), (3, 1003, 2010), (4, 995, 2000)] {
        variant_records.insert(VariantRecord::new(
            read_id, 0, 1,
            GraphOperation::new(
                0, position_1, Strand::Forward, GraphOperationType::Downstream,
                0, position_2, Strand::Forward, GraphOperationType::Upstream,
                "".into(), VariantType::Breakpoint
            )
        ));
    }
    let variant_call: VariantCall = VariantCall::from_variant_records(1, variant_records, 0, 4, 6, 2);
    assert_eq!(variant_call.get_consensus_graph_operation().get_position_1(), 1000);

    // Side 1 spans 995..1003 = 8; side 2 spans 2000..2010 = 10.
    assert_eq!(pooled_breakpoint_spread(&variant_call), Some(10));
}

#[test]
fn pooled_breakpoint_spread_measures_translocations() {
    let mut variant_records: HashSet<VariantRecord> = HashSet::new();
    for (read_id, position_1, position_2) in [(1, 1000, 2000), (2, 1000, 2000), (3, 1004, 2001)] {
        variant_records.insert(VariantRecord::new(
            read_id, 0, 1,
            GraphOperation::new(
                0, position_1, Strand::Forward, GraphOperationType::Downstream,
                1, position_2, Strand::Forward, GraphOperationType::Upstream,
                "".into(), VariantType::Translocation
            )
        ));
    }
    let variant_call: VariantCall = VariantCall::from_variant_records(1, variant_records, 0, 4, 6, 2);

    assert_eq!(pooled_breakpoint_spread(&variant_call), Some(4));
}

/// An unresolved member (Noop second side) carries no partner position; it must not
/// stretch either side's spread.
#[test]
fn pooled_breakpoint_spread_excludes_unresolved_members() {
    let mut variant_records: HashSet<VariantRecord> = HashSet::new();
    for (read_id, position_1, position_2) in [(1, 1000, 2000), (2, 1000, 2000), (3, 1003, 2010)] {
        variant_records.insert(VariantRecord::new(
            read_id, 0, 1,
            GraphOperation::new(
                0, position_1, Strand::Forward, GraphOperationType::Downstream,
                0, position_2, Strand::Forward, GraphOperationType::Upstream,
                "".into(), VariantType::Breakpoint
            )
        ));
    }
    // Spelled the way the caller spells an unanchored clip: same chromosome, same position,
    // Noop on the second side.
    variant_records.insert(VariantRecord::new(
        4, 0, 1,
        GraphOperation::new(
            0, 5000, Strand::Forward, GraphOperationType::Downstream,
            0, 5000, Strand::Forward, GraphOperationType::Noop,
            "".into(), VariantType::Breakpoint
        )
    ));
    let variant_call: VariantCall = VariantCall::from_variant_records(1, variant_records, 0, 4, 6, 2);

    assert_eq!(pooled_breakpoint_spread(&variant_call), Some(10));
}

/// One resolved spelling has no spread to measure, which is not a measured spread of zero.
#[test]
fn pooled_breakpoint_spread_returns_none_with_fewer_than_two_resolved_members() {
    let mut variant_records: HashSet<VariantRecord> = HashSet::new();
    for read_id in [1, 2] {
        variant_records.insert(VariantRecord::new(
            read_id, 0, 1,
            GraphOperation::new(
                0, 5000, Strand::Forward, GraphOperationType::Downstream,
                0, 5000, Strand::Forward, GraphOperationType::Noop,
                "".into(), VariantType::Breakpoint
            )
        ));
    }
    let variant_call: VariantCall = VariantCall::from_variant_records(1, variant_records, 0, 4, 6, 2);
    assert_eq!(*variant_call.get_consensus_graph_operation().get_variant_type(), VariantType::Breakpoint);

    assert_eq!(pooled_breakpoint_spread(&variant_call), None);
}

#[test]
fn pooled_breakpoint_spread_returns_none_outside_the_breakend_family() {
    let mut variant_records: HashSet<VariantRecord> = HashSet::new();
    for (read_id, position_1, position_2) in [(1, 1000, 2000), (2, 1000, 2000), (3, 1004, 2001)] {
        variant_records.insert(VariantRecord::new(
            read_id, 0, 1,
            GraphOperation::new(
                0, position_1, Strand::Forward, GraphOperationType::Downstream,
                1, position_2, Strand::Forward, GraphOperationType::Upstream,
                "".into(), VariantType::FusionGene
            )
        ));
    }
    let variant_call: VariantCall = VariantCall::from_variant_records(1, variant_records, 0, 4, 6, 2);

    assert_eq!(pooled_breakpoint_spread(&variant_call), None);
}


/// One junction between two loci, with a GT..AG intron motif and 10 bases of homology across it,
/// spelled once as a read in the gene's orientation gives it and once as a read of the other
/// strand does (both strands flipped). The second reads the motif as CT..AC. It is the same
/// junction, so it has the same homology, and the splice exemption clears it either way.
#[test]
fn classify_template_switch_exempts_a_canonical_splice_junction_spelled_from_either_strand() {
    let directory: TempDir = tempfile::tempdir().unwrap();
    let fasta_file: String = directory.path().join("reference.fa").to_str().unwrap().to_string();
    let sequences: Vec<(Box<str>, Box<str>)> = vec![
        // Leaving side: the read leaves after 30; positions 31-32 = "GT".
        ("chrA".into(), format!("{}{}{}", "A".repeat(30), "GTCCCCCCCCGGGGGGGGGG", "A".repeat(30)).into()),
        // Entering side: the read enters at 51; positions 49-50 = "AG", and 51-60 repeat 31-40.
        ("chrB".into(), format!("{}{}{}{}", "C".repeat(40), "TTTTTTTTAG", "GTCCCCCCCCAAAAAAAAAA", "C".repeat(30)).into())
    ];
    write_fasta_file(&sequences, &fasta_file);
    let fasta_map: FastaMap = FastaMap::new(&fasta_file);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chrA".into(), 0);
    chromosome_names_map.insert("chrB".into(), 1);

    let sense: GraphOperation = GraphOperation::new(
        0, 30, Strand::Forward, GraphOperationType::Downstream,
        1, 51, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );
    let antisense: GraphOperation = GraphOperation::new(
        0, 30, Strand::Reverse, GraphOperationType::Downstream,
        1, 51, Strand::Reverse, GraphOperationType::Upstream,
        "".into(), VariantType::Breakpoint
    );

    for graph_operation in [&sense, &antisense] {
        let homology: JunctionHomology = compute_junction_homology(graph_operation, &chromosome_names_map, &fasta_map, 20).unwrap();
        assert_eq!(homology.canonical_splice, true);
        assert_eq!(homology.total(), 10);
        let evidence: TemplateSwitchEvidence = TemplateSwitchEvidence {
            junction_homology: Some(homology),
            foldback_stem: None,
            at_exon_boundaries: false,
            has_exit_junction: false,
            dispersion_beyond_homology: None,
            num_reads: 5
        };
        let verdict: TemplateSwitchVerdict = classify_template_switch(&evidence, 8, 5, 0, 1_000, 5);
        assert_eq!(verdict, TemplateSwitchVerdict::Clear);
        assert!(!is_template_switch(verdict, graph_operation.get_variant_type()));
    }
}
