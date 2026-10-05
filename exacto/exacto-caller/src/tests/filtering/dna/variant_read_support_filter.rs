use bimap::BiMap;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

use super::*;


/// A translocation of six reads from chr17:7,674,224 (depth 30) to chr18:5,170,100, which no read
/// covers. The minimum read support is 6 at depth 30 and 11 at depth 60. Without a depth for an
/// uncovered flank the call cannot clear it; with one, that flank alone is judged at that depth.
#[test]
fn dna_variant_read_support_filter_judges_an_uncovered_flank_at_depth_when_uncovered() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-dna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());
    let chromosome_names_map: BiMap<Box<str>, u16> = vec![(Box::<str>::from("chr17"), 0u16), (Box::<str>::from("chr18"), 1u16)].into_iter().collect();

    // Contig lengths from the header; the counts are set by hand.
    let mut read_depths: BAMReadDepths = BAMReadDepths::new(bam_file, &format!("{bam_file}.bai"), &HashMap::new(), 1_000);
    read_depths.insert("chr17", 7_674_224, 30, 15, 15);
    read_depths.insert("chr18", 5_170_100, 0, 0, 0);

    let variant_records: HashSet<VariantRecord> = (1..=6)
        .map(|read_id| VariantRecord::new(
            read_id,
            100,
            101,
            GraphOperation::new(
                0, 7_674_224, Strand::Forward, GraphOperationType::Downstream,
                1, 5_170_100, Strand::Forward, GraphOperationType::Upstream,
                "".into(),
                VariantType::Translocation
            )
        ))
        .collect();
    let mut variant_call: VariantCall = VariantCall::from_variant_records(1, variant_records, 0, 4, 6, 2);
    variant_call.set_total_depth(30);

    let read_support_index: DNAVariantReadSupportIndex = DNAVariantReadSupportIndex::new(60, 30, 0.5, 0.001, 0.01, 0.02, 1e-6);
    assert_eq!(read_support_index.get_min_read_support((0, 30)), 6);
    assert_eq!(read_support_index.get_min_read_support((0, 60)), 11);

    let germline: DNAVariantReadSupportFilter = DNAVariantReadSupportFilter::new(
        4, 4, 4, 4, 30, None, &read_support_index, &read_depths, &chromosome_names_map, &fasta_map
    );
    assert!(!germline.passes(&variant_call));

    let at_depth_30: DNAVariantReadSupportFilter = DNAVariantReadSupportFilter::new(
        4, 4, 4, 4, 30, Some(30), &read_support_index, &read_depths, &chromosome_names_map, &fasta_map
    );
    assert!(at_depth_30.passes(&variant_call));

    let at_depth_60: DNAVariantReadSupportFilter = DNAVariantReadSupportFilter::new(
        4, 4, 4, 4, 30, Some(60), &read_support_index, &read_depths, &chromosome_names_map, &fasta_map
    );
    assert!(!at_depth_60.passes(&variant_call));

    // A covered flank keeps its own depth: at 60 it needs 11 reads, whatever an uncovered one needs.
    read_depths.insert("chr17", 7_674_224, 60, 30, 30);
    let covered_at_depth_60: DNAVariantReadSupportFilter = DNAVariantReadSupportFilter::new(
        4, 4, 4, 4, 30, Some(30), &read_support_index, &read_depths, &chromosome_names_map, &fasta_map
    );
    assert!(!covered_at_depth_60.passes(&variant_call));
}
