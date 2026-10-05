use abpoa_rs::AlignmentMode;
use std::collections::{HashMap, HashSet};
use std::io::Write;

use crate::prelude::*;


fn cluster(reads: &[&str]) -> Vec<Box<str>> {
    reads.iter().map(|read| Box::<str>::from(*read)).collect()
}

/// abPOA takes the magnitude of the match and mismatch scores, so either sign gives the same
/// consensus sequence, and a gap open score of 0 selects linear gaps. The gap scores that end
/// abPOA are refused (`refuses_negative_gap_open_score`, `refuses_gap_extend_score_below_1`).
///
///   Match   Mismatch   Gap open   Gap extend   Consensus sequence
///   0       4          6          2            the majority (the defaults)
///   0       -4         6          2            the majority
///   2       4          6          2            the majority
///   -2      -4         6          2            the majority
///   0       4          0          2            the majority
#[test]
fn accepts_match_and_mismatch_scores_of_either_sign_and_linear_gaps() {
    let reads = cluster(&[
        "ACGTACGTACGTTTGACCATGCTAGCTAGCTA",
        "ACGTACGTACGTTTGACCATGCTAGCTAGCTA",
        "ACGTACGTACGTTTGACCATGCTAGCTAGCTA",
        "ACGTACGTACGTTTGACCATGCTAGCTAGCTA",
        "ACGTACGTACGTTTGAACCATGCTAGCTAGCTA" // one inserted A, minority of one
    ]);

    let options: IdentifyConsensusSequencesOptions = Default::default();
    assert_eq!(
        (options.poa_match_score, options.poa_mismatch_score, options.poa_gap_open_score, options.poa_gap_extend_score),
        (0, 4, 6, 2)
    );

    for (match_score, mismatch_score, gap_open_score, gap_extend_score) in [
        (0, 4, 6, 2),
        (0, -4, 6, 2),
        (2, 4, 6, 2),
        (-2, -4, 6, 2),
        (0, 4, 0, 2)
    ] {
        let consensus: Box<str> = perform_partial_order_alignment(
            &reads,
            AlignmentMode::Global,
            match_score,
            mismatch_score,
            gap_open_score,
            gap_extend_score
        );

        assert_eq!(
            &*consensus,
            "ACGTACGTACGTTTGACCATGCTAGCTAGCTA",
            "scores ({}, {}, {}, {})",
            match_score,
            mismatch_score,
            gap_open_score,
            gap_extend_score
        );
    }
}

#[test]
fn recovers_majority_consensus() {
    let truth = "ACGTACGTACGTACGT";
    let reads = cluster(&[
        "ACGTACGTACGTACGT",
        "ACGTACGTACGTACGT",
        "ACGTTCGTACGTACGT", // position 4 A->T
        "ACGTACGAACGTACGT", // position 7 T->A
        "ACGTACGTACGTACGT",
    ]);

    let options: IdentifyConsensusSequencesOptions = Default::default();

    let consensus = perform_partial_order_alignment(
        &reads,
        AlignmentMode::Global,
        options.poa_match_score,
        options.poa_mismatch_score,
        options.poa_gap_open_score,
        options.poa_gap_extend_score
    );

    assert_eq!(&*consensus, truth);
}

#[test]
fn is_invariant_to_read_order() {
    let forward = [
        "ACGTACGTACGTACGT",
        "ACGTTCGTACGTACGT",
        "ACGTACGAACGTACGT",
        "ACGTACGTACGTACGT",
        "ACGTACGTACGTACGT",
    ];
    let mut backward = forward;
    backward.reverse();

    let options: IdentifyConsensusSequencesOptions = Default::default();

    let out_fwd = perform_partial_order_alignment(
        &cluster(&forward),
        AlignmentMode::Global,
        options.poa_match_score,
        options.poa_mismatch_score,
        options.poa_gap_open_score,
        options.poa_gap_extend_score
    );

    let out_bwd = perform_partial_order_alignment(
        &cluster(&backward),
        AlignmentMode::Global,
        options.poa_match_score,
        options.poa_mismatch_score,
        options.poa_gap_open_score,
        options.poa_gap_extend_score
    );

    assert_eq!(out_fwd, out_bwd);
}

#[test]
fn handles_empty_and_singleton_clusters() {
    let empty = Vec::<Box<str>>::new();

    let options: IdentifyConsensusSequencesOptions = Default::default();

    let out = perform_partial_order_alignment(
        &empty,
        AlignmentMode::Global,
        options.poa_match_score,
        options.poa_mismatch_score,
        options.poa_gap_open_score,
        options.poa_gap_extend_score
    );

    assert_eq!(&*out, "");

    let singleton = cluster(&["ACGTACGT"]);

    let out = perform_partial_order_alignment(
        &singleton,
        AlignmentMode::Global,
        options.poa_match_score,
        options.poa_mismatch_score,
        options.poa_gap_open_score,
        options.poa_gap_extend_score
    );

    assert_eq!(&*out, "ACGTACGT");
}

/// Deterministic LCG, so the synthetic cluster below is reproducible.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

/// A truth sequence plus `count` noisy copies carrying substitutions and deletions,
/// standing in for a long-read cluster.
fn noisy_cluster(length: usize, count: usize, sub_pct: u64, del_pct: u64) -> (String, Vec<Box<str>>) {
    const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut rng = Lcg(42);

    let truth: Vec<u8> = (0..length).map(|_| BASES[(rng.next() % 4) as usize]).collect();
    let reads = (0..count)
        .map(|_| {
            let mut read = Vec::with_capacity(length);
            for &base in &truth {
                match rng.next() % 100 {
                    roll if roll < del_pct => continue,
                    roll if roll < del_pct + sub_pct => read.push(BASES[(rng.next() % 4) as usize]),
                    _ => read.push(base),
                }
            }
            Box::<str>::from(String::from_utf8(read).unwrap())
        })
        .collect();

    (String::from_utf8(truth).unwrap(), reads)
}

/// Levenshtein distance, so the assertion below tolerates indels rather than
/// forcing a position-by-position comparison that one deletion would derail.
fn edit_distance(a: &str, b: &str) -> usize {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut curr = vec![0; b.len() + 1];

    for i in 1..=a.len() {
        curr[0] = i;
        for j in 1..=b.len() {
            let substitution = prev[j - 1] + usize::from(a[i - 1] != b[j - 1]);
            curr[j] = substitution.min(prev[j] + 1).min(curr[j - 1] + 1);
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    prev[b.len()]
}

/// Regression guard for the abPOA `ag_backtrack` abort: reads must be added to the
/// graph base-coded, not as raw ASCII. When they are not, abPOA reads past its 5x5
/// scoring matrix and calls `exit(1)`, which takes the whole test binary down rather
/// than failing this test. The short, near-identical reads in the tests above never
/// reach that path -- it needs reads long and divergent enough to force a non-trivial
/// backtrack, as real long-read clusters are.
#[test]
fn handles_long_divergent_reads() {
    let (truth, reads) = noisy_cluster(600, 12, 8, 3);

    let options: IdentifyConsensusSequencesOptions = Default::default();

    let consensus = perform_partial_order_alignment(
        &reads,
        AlignmentMode::Global,
        options.poa_match_score,
        options.poa_mismatch_score,
        options.poa_gap_open_score,
        options.poa_gap_extend_score
    );

    // Read errors are independent, so the consensus should recover the truth up to a
    // few uncorrected deletions in homopolymer runs.
    let distance = edit_distance(&consensus, &truth);
    assert!(
        distance <= truth.len() / 100,
        "consensus diverges from truth by {distance} edits (len {} vs {})",
        consensus.len(),
        truth.len()
    );
}

/// Write a minimal FASTQ (uniform quality strings) for testing.
fn write_fastq(path: &std::path::Path, records: &[(&str, &str)]) {
    let mut file = std::fs::File::create(path).unwrap();
    for (name, sequence) in records {
        writeln!(file, "@{}\n{}\n+\n{}", name, sequence, "I".repeat(sequence.len())).unwrap();
    }
}

/// `identify_consensus_sequences` fetches each cluster's read sequences from the FASTQ by
/// name, computes the consensus, and attaches the cluster's read names. Reads not named in
/// any cluster (the decoy) are ignored.
#[test]
fn identifies_consensus_from_fastq_by_name() {
    let dir = tempfile::tempdir().unwrap();
    let fastq_path = dir.path().join("reads.fastq");
    write_fastq(&fastq_path, &[
        ("r0", "ACGTACGTACGTACGT"),
        ("r1", "ACGTACGTACGTACGT"),
        ("r2", "ACGTTCGTACGTACGT"),
        ("r3", "ACGTACGAACGTACGT"),
        ("r4", "ACGTACGTACGTACGT"),
        ("decoy", "TTTTTTTTTTTTTTTT"), // not in any cluster
    ]);

    let mut clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::new();
    clusters.insert(
        0usize,
        ["r0", "r1", "r2", "r3", "r4"].iter().map(|name| Box::<str>::from(*name)).collect()
    );

    let options: IdentifyConsensusSequencesOptions = Default::default();

    let out = identify_consensus_sequences(
        &clusters,
        fastq_path.to_str().unwrap(),
        &options,
        1
    );

    let cluster_0 = out.sequences.iter().find(|s| s.get_cluster_id() == 0).unwrap();
    assert_eq!(cluster_0.get_consensus_sequence(), "ACGTACGTACGTACGT");
    assert_eq!(cluster_0.get_read_names().len(), 5); // decoy excluded
}


/// Soft-clip sequences for one insertion event have ragged ends: each clip is as long as
/// what was left of its read, so a cluster mixes truncated and full-length spellings.
///
/// `AlignmentMode::Extend` looks like the natural fit for that and is not — its free end
/// gaps stop the heaviest-bundle traversal from behaving like a support-weighted majority.
/// `Global` handles truncation correctly and wins wherever the two disagree, so clip
/// consensus uses `Global`. These are the discriminating cases; keep them passing.
#[test]
fn ragged_clip_consensus_follows_support_under_global() {
    let options: IdentifyConsensusSequencesOptions = Default::default();

    let observed = cluster(&[
        "ATCGC",
        "ATCGCTATC",    // majority
        "ATCGCTATC",    // majority
        "ATCGCTATC",    // majority
        "ATCGCTATC",    // majority
        "ATCGCTATC",    // majority
        "ATCGCTATC",    // majority
        "ATCGCTATC",    // majority
        "ATCGCTATC",    // majority
        "ATCGCTATC",    // majority
        "ATCGCTATC",    // majority
        "CATCGCTATC"
    ]);

    let out = perform_partial_order_alignment(
        &observed,
        AlignmentMode::Global,
        options.poa_match_score,
        options.poa_mismatch_score,
        options.poa_gap_open_score,
        options.poa_gap_extend_score
    );

    assert_eq!(&*out, "ATCGCTATC");

    let short_majority = cluster(&[
        "ATCGC",
        "ATCGC",
        "ATCGC",
        "ATCGC",
        "ATCGCTATCGGATCCTG"
    ]);

    let out = perform_partial_order_alignment(
        &short_majority,
        AlignmentMode::Global,
        options.poa_match_score,
        options.poa_mismatch_score,
        options.poa_gap_open_score,
        options.poa_gap_extend_score
    );

    assert_eq!(&*out, "ATCGC");

    let noisy = cluster(&[
        "ATCGCT",
        "ATCGCTATCGGATCCTG",    // majority
        "ATCGCTATCGGATCCTG",    // majority
        "ATCGCTATCGGATCCTG",    // majority
        "ATCGCTTTCGGATCCTG"     // A>T
    ]);

    let out = perform_partial_order_alignment(
        &noisy,
        AlignmentMode::Global,
        options.poa_match_score,
        options.poa_mismatch_score,
        options.poa_gap_open_score,
        options.poa_gap_extend_score
    );

    assert_eq!(&*out, "ATCGCTATCGGATCCTG");
}

#[test]
fn test_perform_partial_order_alignment() {
    let options: IdentifyConsensusSequencesOptions = Default::default();

    let truth = "TCGAATTCCCGA";

    let reads = cluster(&[
        "TCGAATTCCCG",  // missing last A
        "TTCCGGA",      // C>G
        "AATTCCCGA",    // truncated
        "CGAATTCCCGA"   // truncated
    ]);

    let consensus = perform_partial_order_alignment(
        &reads,
        AlignmentMode::Global,
        options.poa_match_score,
        options.poa_mismatch_score,
        options.poa_gap_open_score,
        options.poa_gap_extend_score
    );

    assert_eq!(&*consensus, truth);
}

/// abPOA ends the process on an empty read, so empty reads are skipped. The reads are sorted
/// shortest first before alignment, which puts an empty read first. A cluster of empty reads
/// has no consensus sequence.
#[test]
fn skips_empty_reads() {
    let options: IdentifyConsensusSequencesOptions = Default::default();

    let reads = cluster(&["", "ACGTACGTAC", "ACGTACGTAC"]);

    let out = perform_partial_order_alignment(
        &reads,
        AlignmentMode::Global,
        options.poa_match_score,
        options.poa_mismatch_score,
        options.poa_gap_open_score,
        options.poa_gap_extend_score
    );

    assert_eq!(&*out, "ACGTACGTAC");

    let reads = cluster(&["", ""]);

    let out = perform_partial_order_alignment(
        &reads,
        AlignmentMode::Global,
        options.poa_match_score,
        options.poa_mismatch_score,
        options.poa_gap_open_score,
        options.poa_gap_extend_score
    );

    assert_eq!(&*out, "");
}

/// abPOA ends the process on a negative gap open penalty, so it is refused before abPOA sees
/// it.
#[test]
#[should_panic(expected = "The gap open score must be at least 0, got -6.")]
fn refuses_negative_gap_open_score() {
    let reads = cluster(&["ACGTACGTAC", "ACGTACGTAC"]);

    perform_partial_order_alignment(&reads, AlignmentMode::Global, 0, 4, -6, 2);
}

/// abPOA ends the process on a gap extend penalty below 1, so it is refused before abPOA sees
/// it.
#[test]
#[should_panic(expected = "The gap extend score must be at least 1, got 0.")]
fn refuses_gap_extend_score_below_1() {
    let reads = cluster(&["ACGTACGTAC", "ACGTACGTAC"]);

    perform_partial_order_alignment(&reads, AlignmentMode::Global, 0, 4, 6, 0);
}
