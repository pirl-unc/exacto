use std::fs;
use std::path::Path;

use super::*;


/// Retention is tested once per transcript strand, and a read is returned when an intron it
/// retains falls short of its minimum read support.
///
/// The reads sit on three exons of TP53 (chr17, reverse strand) and the two introns between
/// them, 7,674,972-7,675,052 and 7,675,237-7,675,993:
///
///    14 reads  retain the first intron and splice the second
///   286 reads  splice the first intron and retain the second
///     2 reads  hold no junction and retain both introns
///
/// The spliced reads go to the reverse strand, and the junction-less reads to both. No read of
/// the forward strand splices an intron, so nothing is tested there. On the reverse strand each
/// intron is retained or spliced by all 302 reads, a minimum of 46 at an error rate of 0.01 and a
/// max FPR of 1e-6. The first intron is retained by 16 reads, which are returned; the second by
/// 288, which are not.
#[test]
fn identify_nascent_rna_read_ids_per_strand_returns_reads_retaining_intron_below_min_read_support() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());
    let chromosome_names_map: BiMap<Box<str>, u16> = BiMap::from_iter([
        (Box::<str>::from("chr17"), 0u16),
        (Box::<str>::from("chr18"), 1u16)
    ]);

    let mut summaries: Vec<RNAReadCharacterizationSummary> = Vec::new();
    for read_id in 0..14usize {
        summaries.push(RNAReadCharacterizationSummary {
            read_id,
            chromosome: Some(0),
            reference_start: 7_676_272,
            reference_end: 7_674_859,
            ends_at_reference_start: true,
            ends_at_reference_end: true,
            splice_junctions: vec![
                SpliceJunction::new(0, 0, 7_675_993, 7_675_237, Strand::Reverse, Strand::Reverse)
            ],
            exons: vec![(0, 7_675_994, 7_676_272), (0, 7_674_859, 7_675_236)],
            offset: 0,
            length: 0
        });
    }
    for read_id in 14..300usize {
        summaries.push(RNAReadCharacterizationSummary {
            read_id,
            chromosome: Some(0),
            reference_start: 7_676_272,
            reference_end: 7_674_859,
            ends_at_reference_start: true,
            ends_at_reference_end: true,
            splice_junctions: vec![
                SpliceJunction::new(0, 0, 7_675_052, 7_674_972, Strand::Reverse, Strand::Reverse)
            ],
            exons: vec![(0, 7_675_053, 7_676_272), (0, 7_674_859, 7_674_971)],
            offset: 0,
            length: 0
        });
    }
    for read_id in 300..302usize {
        summaries.push(RNAReadCharacterizationSummary {
            read_id,
            chromosome: Some(0),
            reference_start: 7_676_272,
            reference_end: 7_674_859,
            ends_at_reference_start: true,
            ends_at_reference_end: true,
            splice_junctions: Vec::new(),
            exons: vec![(0, 7_674_859, 7_676_272)],
            offset: 0,
            length: 0
        });
    }

    let nascent_read_ids: HashSet<usize> = identify_nascent_rna_read_ids_per_strand(
        &mut summaries,
        &chromosome_names_map,
        &fasta_map,
        0.01,
        1e-6,
        2
    );

    assert_eq!(nascent_read_ids, (0..14).chain(300..302).collect::<HashSet<usize>>());
    // Junction-less reads first, then the reverse strand, each in read order.
    assert_eq!(
        summaries.iter().map(|summary| summary.read_id).collect::<Vec<usize>>(),
        (300..302).chain(0..300).collect::<Vec<usize>>()
    );
}
