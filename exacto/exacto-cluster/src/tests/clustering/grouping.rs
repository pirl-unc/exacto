use bimap::BiMap;
use exacto_core::prelude::*;
use std::fs;
use std::path::Path;

use super::*;
use crate::options::ClusterRNAReadsOptions;


/// Two loci are grouped together when the reads with a junction in both reach the minimum read
/// support at the number of reads of the deeper locus.
///
/// The reads are of WSCD1 and ACAP1 (chr17, forward strand, 1.2 Mb apart). A bridging read
/// holds the last junction of WSCD1 and three junctions of ACAP1, so it anchors on ACAP1.
///
///   WSCD1 reads   ACAP1 reads   Bridging reads   Reads of a locus   Minimum   Groups
///   40            40            7                47                 8         40 and 47
///   40            40            8                48                 8         88
///   40            40            8                48                 9         40 and 48
///
/// The minimum of the last row is `min_reads`, which a bridge is never held to less than.
#[test]
fn group_transcript_models_by_shared_junctions_joins_loci_bridged_by_min_read_support() {
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");
    let chromosome_names_map: BiMap<Box<str>, u16> = BiMap::from_iter([(Box::from("chr17"), 0u16)]);
    let mut options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;
    // The documented 7/8-read boundary is specifically for this tail cutoff.
    options.filtering.max_fpr = 1e-4;

    let wscd1_introns: Vec<(u32, u32)> = vec![
        (6_070_653, 6_080_370),
        (6_081_086, 6_087_989),
        (6_088_105, 6_090_320),
        (6_090_506, 6_095_101),
        (6_095_224, 6_109_606),
        (6_109_767, 6_110_770),
        (6_110_936, 6_117_987),
        (6_118_189, 6_120_308)
    ];
    let acap1_introns: Vec<(u32, u32)> = vec![
        (7_336_788, 7_341_947),
        (7_342_068, 7_342_274),
        (7_342_329, 7_342_415),
        (7_342_475, 7_343_378)
    ];
    let bridging_introns: Vec<(u32, u32)> = vec![
        (6_109_767, 6_110_770),
        (7_342_068, 7_342_274),
        (7_342_329, 7_342_415),
        (7_342_475, 7_343_378)
    ];

    // (number of bridging reads, min reads, read IDs of every group)
    let cases: Vec<(usize, usize, Vec<Vec<usize>>)> = vec![
        (7, 3, vec![(0..40).collect(), (40..87).collect()]),
        (8, 3, vec![(0..88).collect()]),
        (8, 9, vec![(0..40).collect(), (40..88).collect()])
    ];
    for (num_bridging_reads, min_reads, expected) in cases {
        let mut summaries: Vec<RNAReadCharacterizationSummary> = Vec::new();
        for (num_reads, introns) in [(40, &wscd1_introns), (40, &acap1_introns), (num_bridging_reads, &bridging_introns)] {
            for _ in 0..num_reads {
                summaries.push(RNAReadCharacterizationSummary {
                    read_id: summaries.len(),
                    chromosome: Some(0),
                    reference_start: introns.first().unwrap().0 - 100,
                    reference_end: introns.last().unwrap().1 + 100,
                    ends_at_reference_start: true,
                    ends_at_reference_end: true,
                    splice_junctions: introns
                        .iter()
                        .map(|&(intron_start, intron_end)| SpliceJunction::new(
                            0,
                            0,
                            intron_start,
                            intron_end,
                            Strand::Forward,
                            Strand::Forward
                        ))
                        .collect(),
                    exons: Vec::new(),
                    offset: 0,
                    length: 0
                });
            }
        }

        let groups: Vec<Vec<usize>> = group_transcript_models_by_shared_junctions(
            &summaries,
            &SpliceJunctionAnnotationIndex::new(&gene_annotator, &chromosome_names_map),
            options.max_locus_gap,
            options.unspliced_bin_size,
            min_reads,
            options.filtering.max_slippage_repeat_len,
            options.error_model.sequencing_error,
            options.error_model.slippage_prob,
            options.filtering.max_fpr,
            2
        )
            .into_iter()
            .map(|group| group.into_iter().map(|summary| summary.read_id).collect())
            .collect();
        assert_eq!(
            groups,
            expected,
            "{} bridging reads, min reads {}",
            num_bridging_reads,
            min_reads
        );
    }
}


#[test]
fn group_transcript_models_by_shared_junctions_joins_novel_junction_to_annotated_intron_it_overlaps() {
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");
    let chromosome_names_map: BiMap<Box<str>, u16> = BiMap::from_iter([(Box::from("chr17"), 0u16)]);
    let options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;

    let tp53_introns: Vec<(u32, u32)> = vec![(7_687_376, 7_676_623), (7_676_520, 7_676_404)];
    // (junction of the fragments, read IDs of every group)
    let cases: Vec<((u32, u32), Vec<Vec<usize>>)> = vec![
        ((7_687_377, 7_676_623), vec![(0..65).collect()]),
        ((7_700_000, 7_676_623), vec![(0..60).collect(), (60..65).collect()])
    ];
    for (fragment_junction, expected) in cases {
        let mut summaries: Vec<RNAReadCharacterizationSummary> = Vec::new();
        for (num_reads, introns) in [(60, tp53_introns.clone()), (5, vec![fragment_junction])] {
            for _ in 0..num_reads {
                summaries.push(RNAReadCharacterizationSummary {
                    read_id: summaries.len(),
                    chromosome: Some(0),
                    reference_start: introns.last().unwrap().1 - 100,
                    reference_end: introns.first().unwrap().0 + 100,
                    ends_at_reference_start: true,
                    ends_at_reference_end: true,
                    splice_junctions: introns
                        .iter()
                        .map(|&(position_1, position_2)| SpliceJunction::new(
                            0,
                            0,
                            position_1,
                            position_2,
                            Strand::Reverse,
                            Strand::Reverse
                        ))
                        .collect(),
                    exons: Vec::new(),
                    offset: 0,
                    length: 0
                });
            }
        }

        let groups: Vec<Vec<usize>> = group_transcript_models_by_shared_junctions(
            &summaries,
            &SpliceJunctionAnnotationIndex::new(&gene_annotator, &chromosome_names_map),
            options.max_locus_gap,
            options.unspliced_bin_size,
            options.junction.min_reads,
            options.filtering.max_slippage_repeat_len,
            options.error_model.sequencing_error,
            options.error_model.slippage_prob,
            options.filtering.max_fpr,
            2
        )
            .into_iter()
            .map(|group| group.into_iter().map(|summary| summary.read_id).collect())
            .collect();
        assert_eq!(groups, expected, "fragments splice {:?}", fragment_junction);
    }
}
