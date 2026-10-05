use exacto_core::prelude::*;
use std::fs;
use std::path::Path;

use super::*;


/// Two spellings of one donor converge on the same canonical position.
///
/// On `scga-mini-rna-015-tumor` the aligner placed this boundary at six positions across 289 reads.
/// The two modal ones became clusters 3 and 4, one transcript reported as two, and the remaining
/// 19 reads fell below `min_reads_per_cluster` in their minority spelling and were dropped. The
/// spellings differ by exactly the reference bases they roll into the front of the insertion:
/// chr17:6,110,797-6,110,799 is `CTG`, and cluster 4's insertion is `CTG` + cluster 3's.
///
/// The canonical position is the annotated GENCODE donor, chr17:6,110,936, whose motif is a
/// canonical `GT`; every placement the aligner chose was non-canonical.
#[test]
fn normalise_junction_boundary_converges_two_donor_spellings() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());

    // The 467 bp insertion of cluster 3, whose leading 136 bases are exactly the reference bases
    // chr17:6,110,800-6,110,935, which is what lets the boundary slide.
    let cluster_3_insertion: &str = "CCTAACAAATCCAAAGTGTTTGTGGCTTTGTCAAGCTTCCCAGGAGCCGGGAACACATGGGCACGGCACCTCATTGAGCATGCCACTGGCTTCTATACAGGGAGCTACTACTTTGATGGAACCCTCTACAACAAAGCTCCTGAAACTGGGCACTGGTCTCCTGGAAAGTGGGCGCCATTACCTTGCTGCCAGCCGCGCCTTCGTTGTCGGCATTTGTGACCTGGCCCGCCTGGGTCCACCAGAGCCCATGATGGCGGAGTGTCTGGAAAAATTCACCGTGAGCCTGAACCACAAGCTGGACAGCCATGCGGAGCTTCTAGATGCCACCCAACACACACTGCAGCAGCAGATCCAGACCCTGGTCAAGGAAGGTCTGCGGGGTTTCCGAGAGGCTCGCCGGGATTTCTGGCGGGGGGCTGAGAGCCTGGAGGCTGCCCTGACCCACAACGCAGAGGTTCCCAGGCGCC";
    let cluster_4_insertion: String = format!("CTG{}", cluster_3_insertion);

    // position_1 is the first intron base, so cluster 3's exon ends at 6,110,799 and cluster 4's at
    // 6,110,796. position_2 is the acceptor side and is the same for both.
    let from_cluster_3: u32 = normalise_junction_boundary(
        6_110_800, 6_165_488, cluster_3_insertion, JunctionBoundary::Donor, "chr17", &fasta_map
    );
    let from_cluster_4: u32 = normalise_junction_boundary(
        6_110_797, 6_165_488, cluster_4_insertion.as_str(), JunctionBoundary::Donor, "chr17", &fasta_map
    );

    assert_eq!(from_cluster_3, 6_110_936, "cluster 3's spelling did not reach the annotated donor");
    assert_eq!(from_cluster_4, 6_110_936, "cluster 4's spelling did not reach the annotated donor");
    assert_eq!(
        from_cluster_3, from_cluster_4,
        "the two spellings of one boundary did not converge, so they would still cluster apart"
    );
}

/// A boundary with nothing to slide against stays put.
///
/// Both directions of the guard matter: normalisation must be inert on the ordinary case, or it
/// would move every junction in the library rather than the ambiguous ones.
#[test]
fn normalise_junction_boundary_leaves_unambiguous_junctions_alone() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());

    // No insertion at all: the overwhelmingly common case.
    assert_eq!(
        normalise_junction_boundary(6_110_800, 6_165_488, "", JunctionBoundary::Donor, "chr17", &fasta_map),
        6_110_800
    );

    // An insertion whose leading base does not match the reference base the boundary would cross.
    // chr17:6,110,800 is `C`, so an insertion starting `A` pins the boundary where it is.
    assert_eq!(
        normalise_junction_boundary(6_110_800, 6_165_488, "AAAAA", JunctionBoundary::Donor, "chr17", &fasta_map),
        6_110_800
    );
}

/// The rule is written in transcript order, so the reverse strand slides the other way.
///
/// A junction with `position_1 > position_2` is on the reverse strand, and the bases it is compared
/// against are complemented. chr17:6,110,798-6,110,800 reads `TGC` forward, so in transcript order
/// the boundary crosses `G`, then `C`, then `A` as it slides toward lower coordinates.
#[test]
fn normalise_junction_boundary_slides_reverse_strand_toward_lower_coordinates() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());

    // Slides three bases, from 6,110,800 down to 6,110,797.
    assert_eq!(
        normalise_junction_boundary(6_110_800, 6_110_000, "GCA", JunctionBoundary::Donor, "chr17", &fasta_map),
        6_110_797
    );

    // The forward-strand spelling of the same bases must not slide it: `GCA` is not what the
    // forward strand reads there, so a strand mix-up would show up here.
    assert_eq!(
        normalise_junction_boundary(6_110_800, 6_165_488, "GCA", JunctionBoundary::Donor, "chr17", &fasta_map),
        6_110_800
    );
}

/// An insertion on the other side of the intron makes the acceptor ambiguous, and it converges too.
///
/// `scga-mini-rna-013-tumor` is this case. One boundary was spelled at chr17:7,669,665 / 7,669,670 /
/// 7,669,675 with insertions of 751 / 746 / 741 bases, acceptor plus insertion length constant at
/// 7,670,416, which is what identifies them as one event rather than three, and every other
/// junction in the three chains was identical. They became three clusters, two of which assembled
/// byte-identical transcripts, because a donor-only rule cannot see this ambiguity at all.
///
/// The reference here is `...GACTGACCCTTTTT...` on the forward strand, and the transcript is on the
/// reverse, so the acceptor slides across the complements `T G A C T G` reading up from 7,669,665.
/// Those are the insertion's trailing bases, since it is the following exon that grows backwards
/// over the boundary: the mirror of the donor, which consumes leading bases.
#[test]
fn normalise_junction_boundary_converges_two_acceptor_spellings() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());

    // Two spellings of one molecule: the earlier acceptor carries six more inserted bases, which are
    // exactly the six reference bases between the two placements.
    let spelled_early: u32 = normalise_junction_boundary(
        7_673_700, 7_669_665, "AAAGTCAGT", JunctionBoundary::Acceptor, "chr17", &fasta_map
    );
    let spelled_late: u32 = normalise_junction_boundary(
        7_673_700, 7_669_671, "AAA", JunctionBoundary::Acceptor, "chr17", &fasta_map
    );

    assert_eq!(spelled_early, 7_669_671, "the acceptor did not slide across its microhomology");
    assert_eq!(spelled_late, 7_669_671, "the already-canonical acceptor moved with nothing to slide against");
    assert_eq!(
        spelled_early, spelled_late,
        "the two spellings of one acceptor did not converge, so they would still cluster apart"
    );
}

/// The two boundaries are placed independently, and neither rule fires on the other's side.
///
/// This is the guard that keeps the fix a strict generalisation: asking for the acceptor of a
/// junction whose insertion sits at the donor must return the acceptor untouched, and vice versa.
/// Both use `scga-mini-rna-015-tumor`'s donor-side insertion, whose leading bases match at
/// chr17:6,110,800 and whose trailing bases have no reason to match at the acceptor.
#[test]
fn normalise_junction_boundary_places_each_boundary_independently() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());

    let cluster_3_insertion: &str = "CCTAACAAATCCAAAGTGTTTGTGGCTTTGTCAAGCTTCCCAGGAGCCGGGAACACATGGGCACGGCACCTCATTGAGCATGCCACTGGCTTCTATACAGGGAGCTACTACTTTGATGGAACCCTCTACAACAAAGCTCCTGAAACTGGGCACTGGTCTCCTGGAAAGTGGGCGCCATTACCTTGCTGCCAGCCGCGCCTTCGTTGTCGGCATTTGTGACCTGGCCCGCCTGGGTCCACCAGAGCCCATGATGGCGGAGTGTCTGGAAAAATTCACCGTGAGCCTGAACCACAAGCTGGACAGCCATGCGGAGCTTCTAGATGCCACCCAACACACACTGCAGCAGCAGATCCAGACCCTGGTCAAGGAAGGTCTGCGGGGTTTCCGAGAGGCTCGCCGGGATTTCTGGCGGGGGGCTGAGAGCCTGGAGGCTGCCCTGACCCACAACGCAGAGGTTCCCAGGCGCC";

    // The donor still reaches the annotated GENCODE donor, exactly as before the acceptor rule existed.
    assert_eq!(
        normalise_junction_boundary(
            6_110_800, 6_165_488, cluster_3_insertion, JunctionBoundary::Donor, "chr17", &fasta_map
        ),
        6_110_936,
        "adding the acceptor case changed the donor's answer"
    );

    // And the acceptor of that same junction does not move, because nothing abuts it.
    assert_eq!(
        normalise_junction_boundary(
            6_110_800, 6_165_488, cluster_3_insertion, JunctionBoundary::Acceptor, "chr17", &fasta_map
        ),
        6_165_488,
        "the acceptor slid on an insertion that sits at the donor"
    );

    // An empty insertion is inert on both ends.
    assert_eq!(
        normalise_junction_boundary(7_673_700, 7_669_665, "", JunctionBoundary::Acceptor, "chr17", &fasta_map),
        7_669_665
    );
}
/// The slide is exactly as long as the insertion's matching prefix: it stops at the first
/// inserted base that disagrees with the reference, and it cannot outrun a short insertion.
///
/// chr17:6,110,800-6,110,806 reads `CCTAACA` on the forward strand.
#[test]
fn normalise_junction_boundary_slides_exactly_the_matching_prefix() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());

    // Five bases match, the sixth (`G` against `C`) does not: the donor moves five.
    assert_eq!(
        normalise_junction_boundary(6_110_800, 6_165_488, "CCTAAG", JunctionBoundary::Donor, "chr17", &fasta_map),
        6_110_805
    );

    // Three inserted bases, all matching: the donor moves three and the insertion is spent.
    assert_eq!(
        normalise_junction_boundary(6_110_800, 6_165_488, "CCT", JunctionBoundary::Donor, "chr17", &fasta_map),
        6_110_803
    );
}

/// A forward-strand acceptor slides toward lower coordinates, consuming the insertion's
/// trailing bases against the reference read downward from the acceptor.
///
/// The mirror of the reverse-strand acceptor case, on the same `CCTAAC` stretch: reading down
/// from 6,110,805 the reference is `C A A T C`, so an insertion ending `TAAC` moves it four.
#[test]
fn normalise_junction_boundary_slides_forward_strand_acceptor_toward_lower_coordinates() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());

    assert_eq!(
        normalise_junction_boundary(6_110_000, 6_110_805, "GGGTAAC", JunctionBoundary::Acceptor, "chr17", &fasta_map),
        6_110_801
    );

    // The same insertion's leading bases are not what the donor side reads, so it stays.
    assert_eq!(
        normalise_junction_boundary(6_110_000, 6_110_805, "GGGTAAC", JunctionBoundary::Donor, "chr17", &fasta_map),
        6_110_000
    );
}

/// A boundary never slides past the other end of its intron: reaching it would dissolve the
/// junction. Both boundaries stop at the limit even with matching insertion bases to spare.
#[test]
fn normalise_junction_boundary_stops_at_the_far_end_of_the_intron() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());

    // The cluster-3 insertion matches 136 reference bases from 6,110,800, but the intron is
    // only three bases long: the donor halts on the acceptor.
    let cluster_3_insertion: &str = "CCTAACAAATCCAAAGTGTTTGTGGCTTTGTCAAGCTTCCCAGGAGCCGGGAACACATGGGCACGGCACCTCATTGAGCATGCCACTGGCTTCTATACAGGGAGCTACTACTTTGATGGAACCCTCTACAACAAAGCTCCTGAAACTGGGCACTGGTCTCCTGGAAAGTGGGCGCCATTACCTTGCTGCCAGCCGCGCCTTCGTTGTCGGCATTTGTGACCTGGCCCGCCTGGGTCCACCAGAGCCCATGATGGCGGAGTGTCTGGAAAAATTCACCGTGAGCCTGAACCACAAGCTGGACAGCCATGCGGAGCTTCTAGATGCCACCCAACACACACTGCAGCAGCAGATCCAGACCCTGGTCAAGGAAGGTCTGCGGGGTTTCCGAGAGGCTCGCCGGGATTTCTGGCGGGGGGCTGAGAGCCTGGAGGCTGCCCTGACCCACAACGCAGAGGTTCCCAGGCGCC";
    assert_eq!(
        normalise_junction_boundary(6_110_800, 6_110_803, cluster_3_insertion, JunctionBoundary::Donor, "chr17", &fasta_map),
        6_110_803
    );

    // The reverse-strand acceptor of the scga-mini-rna-013 case would slide six; with the donor
    // three bases up the intron, it halts on the donor.
    assert_eq!(
        normalise_junction_boundary(7_669_668, 7_669_665, "AAAGTCAGT", JunctionBoundary::Acceptor, "chr17", &fasta_map),
        7_669_668
    );
}

/// Equal positions are not an intron. The guard returns the boundary untouched on both sides,
/// even with an insertion that would otherwise slide it.
#[test]
fn normalise_junction_boundary_leaves_zero_length_introns_alone() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());

    assert_eq!(
        normalise_junction_boundary(6_110_800, 6_110_800, "CCTAAC", JunctionBoundary::Donor, "chr17", &fasta_map),
        6_110_800
    );
    assert_eq!(
        normalise_junction_boundary(6_110_800, 6_110_800, "CCTAAC", JunctionBoundary::Acceptor, "chr17", &fasta_map),
        6_110_800
    );
}

/// Insertion bases are compared case-insensitively, so a lowercase spelling of an insertion
/// slides exactly as far as its uppercase spelling.
#[test]
fn normalise_junction_boundary_matches_insertion_bases_case_insensitively() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());

    let cluster_3_insertion: &str = "CCTAACAAATCCAAAGTGTTTGTGGCTTTGTCAAGCTTCCCAGGAGCCGGGAACACATGGGCACGGCACCTCATTGAGCATGCCACTGGCTTCTATACAGGGAGCTACTACTTTGATGGAACCCTCTACAACAAAGCTCCTGAAACTGGGCACTGGTCTCCTGGAAAGTGGGCGCCATTACCTTGCTGCCAGCCGCGCCTTCGTTGTCGGCATTTGTGACCTGGCCCGCCTGGGTCCACCAGAGCCCATGATGGCGGAGTGTCTGGAAAAATTCACCGTGAGCCTGAACCACAAGCTGGACAGCCATGCGGAGCTTCTAGATGCCACCCAACACACACTGCAGCAGCAGATCCAGACCCTGGTCAAGGAAGGTCTGCGGGGTTTCCGAGAGGCTCGCCGGGATTTCTGGCGGGGGGCTGAGAGCCTGGAGGCTGCCCTGACCCACAACGCAGAGGTTCCCAGGCGCC";
    let lowercase_insertion: String = cluster_3_insertion.to_lowercase();

    assert_eq!(
        normalise_junction_boundary(6_110_800, 6_165_488, lowercase_insertion.as_str(), JunctionBoundary::Donor, "chr17", &fasta_map),
        6_110_936
    );
    assert_eq!(
        normalise_junction_boundary(6_110_800, 6_110_000, "gca", JunctionBoundary::Donor, "chr17", &fasta_map),
        6_110_797
    );
}


/// A read retains an intron when one of its exon blocks holds the exonic base on either side of
/// it, or when it holds the exonic base on one side and ends among the bases of the intron. A
/// read that ends before reaching the intron is uninformative, and so is a block that stops
/// inside the intron where the read goes on.
///
/// Two introns of chr18 overlap: A is 1,000,101-1,000,900 and B is 1,000,501-1,002,000. The
/// exonic bases beside A are 1,000,100 and 1,000,901, and the exonic base before B is 1,000,500.
///
/// The read goes on past both ends of its blocks:
///
///   Read                               Exon blocks                                Retains
///   covers A and B                     1,000,001-1,003,000                        A, B
///   covers A, stops inside B           1,000,001-1,001,000                        A
///   holds the two bases beside A       1,000,100-1,000,901                        A
///   starts at the first base of A      1,000,101-1,000,901                        none
///   stops at the last base of A        1,000,100-1,000,900                        none
///   stops inside A from below          1,000,001-1,000,500                        none
///   stops inside A from above          1,000,500-1,001,000                        none
///   stops at the last base before A    1,000,001-1,000,100                        none
///   starts inside A, stops inside B    1,000,200-1,000,800                        none
///   splices A                          1,000,001-1,000,100, 1,000,901-1,001,000   none
///   splices A, covers it elsewhere     1,000,001-1,001,000                        none
///
/// The last read splices A in one alignment record and covers it in another: it spliced it.
///
/// The read ends at the last base of its blocks, at the first, or at both:
///
///   Read                               Exon blocks            Ends at      Retains
///   ends inside A from below           1,000,001-1,000,500    last         A
///   ends at the first base of A        1,000,001-1,000,101    last         A
///   ends at the last base of A         1,000,100-1,000,900    last         A, B
///   ends at the last base before A     1,000,001-1,000,100    last         none
///   covers A, ends inside B            1,000,001-1,001,000    last         A, B
///   starts inside A                    1,000,700-1,001,000    first        A
///   starts at the first base of A      1,000,101-1,000,901    first        A
///   starts at the base after A         1,000,901-1,001,000    first        none
///   starts inside A, ends inside B     1,000,200-1,000,800    both         B
///   starts inside A, stops inside B    1,000,200-1,000,800    first        none
///
/// The read that starts inside A and ends inside B holds the exonic base before B and no
/// exonic base beside A.
#[test]
fn identify_retained_introns_returns_matches() {
    let introns: Vec<(u16, u32, u32)> = vec![(1, 1_000_101, 1_000_900), (1, 1_000_501, 1_002_000)];
    let introns_by_end: Vec<(u16, u32, u32)> = vec![(1, 1_000_900, 1_000_101), (1, 1_002_000, 1_000_501)];

    // (exon blocks, splice junctions, ends at the first base, ends at the last base, retained introns)
    let cases: Vec<(Vec<(u16, u32, u32)>, Vec<SpliceJunction>, bool, bool, Vec<usize>)> = vec![
        (vec![(1, 1_000_001, 1_003_000)], Vec::new(), false, false, vec![0, 1]),
        (vec![(1, 1_000_001, 1_001_000)], Vec::new(), false, false, vec![0]),
        (vec![(1, 1_000_100, 1_000_901)], Vec::new(), false, false, vec![0]),
        (vec![(1, 1_000_101, 1_000_901)], Vec::new(), false, false, Vec::new()),
        (vec![(1, 1_000_100, 1_000_900)], Vec::new(), false, false, Vec::new()),
        (vec![(1, 1_000_001, 1_000_500)], Vec::new(), false, false, Vec::new()),
        (vec![(1, 1_000_500, 1_001_000)], Vec::new(), false, false, Vec::new()),
        (vec![(1, 1_000_001, 1_000_100)], Vec::new(), false, false, Vec::new()),
        (vec![(1, 1_000_200, 1_000_800)], Vec::new(), false, false, Vec::new()),
        (
            vec![(1, 1_000_001, 1_000_100), (1, 1_000_901, 1_001_000)],
            vec![SpliceJunction::new(1, 1, 1_000_101, 1_000_900, Strand::Forward, Strand::Forward)],
            false,
            false,
            Vec::new()
        ),
        (
            vec![(1, 1_000_001, 1_001_000)],
            vec![SpliceJunction::new(1, 1, 1_000_101, 1_000_900, Strand::Forward, Strand::Forward)],
            false,
            false,
            Vec::new()
        ),
        (vec![(1, 1_000_001, 1_000_500)], Vec::new(), false, true, vec![0]),
        (vec![(1, 1_000_001, 1_000_101)], Vec::new(), false, true, vec![0]),
        (vec![(1, 1_000_100, 1_000_900)], Vec::new(), false, true, vec![0, 1]),
        (vec![(1, 1_000_001, 1_000_100)], Vec::new(), false, true, Vec::new()),
        (vec![(1, 1_000_001, 1_001_000)], Vec::new(), false, true, vec![0, 1]),
        (vec![(1, 1_000_700, 1_001_000)], Vec::new(), true, false, vec![0]),
        (vec![(1, 1_000_101, 1_000_901)], Vec::new(), true, false, vec![0]),
        (vec![(1, 1_000_901, 1_001_000)], Vec::new(), true, false, Vec::new()),
        (vec![(1, 1_000_200, 1_000_800)], Vec::new(), true, true, vec![1]),
        (vec![(1, 1_000_200, 1_000_800)], Vec::new(), true, false, Vec::new())
    ];
    for (exons, splice_junctions, ends_at_reference_start, ends_at_reference_end, expected) in cases {
        let summary: RNAReadCharacterizationSummary = RNAReadCharacterizationSummary {
            read_id: 0,
            chromosome: Some(1),
            reference_start: exons.first().unwrap().1,
            reference_end: exons.last().unwrap().2,
            ends_at_reference_start,
            ends_at_reference_end,
            splice_junctions,
            exons: exons.clone(),
            offset: 0,
            length: 0
        };
        assert_eq!(
            summary.identify_retained_introns(&introns, &introns_by_end),
            expected,
            "exon blocks {:?}, ends at the first base {}, at the last base {}",
            exons, ends_at_reference_start, ends_at_reference_end
        );
    }
}
