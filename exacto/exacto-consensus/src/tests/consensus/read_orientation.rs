use super::*;


/// The reads in the orientation of fewer reads are reverse-complemented.
///
/// The transcript is 60 bases. The longest read is the transcript.
///
///   Read                               Orientation      Reverse-complemented
///   transcript                         longest read     no
///   bases 1-50                         same             no
///   bases 11-60                        same             no
///   bases 6-55, reverse complement     opposite         yes
///   bases 1-40, reverse complement     opposite         yes
#[test]
fn orient_reads_reverse_complements_reads_of_the_orientation_of_fewer_reads() {
    let mut reads: Vec<Box<str>> = vec![
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC".into(),
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCG".into(),
        "TTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC".into(),
        "TTCATCGGATACCGTTAAGCTGCAACTTGGATCCTGAACTGCTAAGGTAC".into(),
        "TAAGCTGCAACTTGGATCCTGAACTGCTAAGGTACGCCAT".into()
    ];

    let num_reversed: usize = orient_reads(&mut reads, 15);

    assert_eq!(num_reversed, 2);
    assert_eq!(reads, vec![
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC".into(),
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCG".into(),
        "TTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC".into(),
        "GTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAA".into(),
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTA".into()
    ] as Vec<Box<str>>);
}


/// The orientation of more reads wins over that of the longest read, and that of the longest
/// read wins when as many reads are in each.
///
///   Reads                                                      Reverse-complemented
///   transcript reverse-complemented (longest), 3 reads same    the longest read
///   as the transcript
///   transcript (longest), bases 1-50 reverse-complemented      bases 1-50
#[test]
fn orient_reads_keeps_the_orientation_of_more_reads_then_of_the_longest_read() {
    let mut reads: Vec<Box<str>> = vec![
        "GCAAGTTCATCGGATACCGTTAAGCTGCAACTTGGATCCTGAACTGCTAAGGTACGCCAT".into(),
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCG".into(),
        "TTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC".into(),
        "GTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAA".into()
    ];

    let num_reversed: usize = orient_reads(&mut reads, 15);

    assert_eq!(num_reversed, 1);
    assert_eq!(reads, vec![
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC".into(),
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCG".into(),
        "TTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC".into(),
        "GTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAA".into()
    ] as Vec<Box<str>>);

    let mut reads: Vec<Box<str>> = vec![
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC".into(),
        "CGGATACCGTTAAGCTGCAACTTGGATCCTGAACTGCTAAGGTACGCCAT".into()
    ];

    let num_reversed: usize = orient_reads(&mut reads, 15);

    assert_eq!(num_reversed, 1);
    assert_eq!(reads, vec![
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC".into(),
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCG".into()
    ] as Vec<Box<str>>);
}


/// A read that shares as many k-mers with the longest read in either orientation is left as it
/// is: one shorter than the k-mers, and one that does not overlap the longest read. When the
/// longest read is its own reverse complement, as a fold-back chimera is, no k-mer tells the
/// orientations apart and every read is left as it is. So is every read with a k-mer size of 0.
#[test]
fn orient_reads_leaves_reads_it_cannot_orient() {
    let mut reads: Vec<Box<str>> = vec![
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC".into(),
        "GCAAGTTCATCGGATACCGTTAAGCTGCAACTTGGATCCTGAACTGCTAA".into(),     // bases 11-60, reverse complement
        "GTTCATCGGA".into(),                                              // 10 bases
        "CCCTTTGGGAAACCCTTTGGGAAACCCTTTGGG".into()                        // another sequence
    ];

    let num_reversed: usize = orient_reads(&mut reads, 15);

    assert_eq!(num_reversed, 1);
    assert_eq!(reads, vec![
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC".into(),
        "TTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC".into(),
        "GTTCATCGGA".into(),
        "CCCTTTGGGAAACCCTTTGGGAAACCCTTTGGG".into()
    ] as Vec<Box<str>>);

    let folded_back: Vec<Box<str>> = vec![
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGCGCAAGTTCATCGGATACCGTTAAGCTGCAACTTGGATCCTGAACTGCTAAGGTACGCCAT".into(),
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCG".into(),
        "GCAAGTTCATCGGATACCGTTAAGCTGCAACTTGGATCCTGAACTGCTAA".into()
    ];
    let mut reads: Vec<Box<str>> = folded_back.clone();

    let num_reversed: usize = orient_reads(&mut reads, 15);

    assert_eq!(num_reversed, 0);
    assert_eq!(reads, folded_back);

    let mixed: Vec<Box<str>> = vec![
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC".into(),
        "GCAAGTTCATCGGATACCGTTAAGCTGCAACTTGGATCCTGAACTGCTAA".into()
    ];
    let mut reads: Vec<Box<str>> = mixed.clone();

    let num_reversed: usize = orient_reads(&mut reads, 0);

    assert_eq!(num_reversed, 0);
    assert_eq!(reads, mixed);
}
