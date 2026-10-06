use super::*;
use crate::options::ClusterRNAReadsOptions;
use std::path::Path;

/// Reads 10, 100, 103 and 104 of scga-mini-rna-007 align the ASPA-WSCD1 fusion across two
/// records; read 154 holds 62 bases of the ASPA arm as a leading clip, the aligner places the
/// fusion's untemplated AGG at the end of WSCD1's intron, and the read splices only WSCD1. Every
/// cluster below is keyed by read 154's chain.
#[test]
fn breakpoint_clip_needs_aligned_support_and_moves_a_read_only_to_one_cluster() {
    use exacto_core::prelude::{create_chromosome_names_map, index_bam_records, Gencode};
    let root = Path::new(env!("EXACTO_TEST_DATA"));
    let bam = root.join("alignment/scga-mini-rna-007-tumor_minimap2_sorted.bam");
    let bam = bam.to_str().unwrap();
    let fasta = FastaMap::new(root.join("references/hg38_chr17-18.fa.gz").to_str().unwrap());
    let annotator = Gencode::new_with_defaults(
        root.join("references/gencode.v41.annotation.chr17-18.gtf.gz").to_str().unwrap(),
        "hg38",
        "v41",
    );
    let chromosomes = create_chromosome_names_map(bam);
    let (positions, names) = index_bam_records(bam, true, 1);
    let ids: Vec<usize> = [10, 100, 103, 104, 154]
        .iter()
        .map(|n| *names.get_by_left(format!("scga-mini-rna-007-tumor_chunk_0000/{n}/ccs").as_str()).unwrap())
        .collect();
    let mut options = ClusterRNAReadsOptions::PACBIO_HIFI;
    let c = &options.calling;
    let temp = tempfile::tempdir().unwrap();
    let (file, summaries, _) = characterize_rna_reads(
        bam,
        &ids,
        &names,
        &positions,
        &chromosomes,
        &fasta,
        &annotator,
        1,
        c.min_mapping_quality,
        c.min_terminal_soft_clip_ins_len,
        c.bkpt_rescue,
        c.bkpt_rescue_min_ins_len,
        c.bkpt_rescue_realignment_gap_open_score,
        c.bkpt_rescue_realignment_gap_extend_score,
        c.bkpt_rescue_realignment_k,
        c.bkpt_rescue_realignment_band_width,
        c.bkpt_rescue_realignment_min_score_fraction,
        c.bkpt_rescue_realignment_min_query_coverage,
        c.bkpt_rescue_realignment_min_placed_fraction,
        c.bkpt_rescue_realignment_max_pieces,
        None,
        c.chunk_size,
        temp.path().to_str().unwrap(),
    );
    let original = load_transcript_models(&file, &summaries);
    let clipped = |models: &[TranscriptModel]| -> Vec<(VariantType, u32, u32)> {
        models
            .iter()
            .find(|m| m.get_read_id() == ids[4])
            .unwrap()
            .get_variant_records()
            .iter()
            .filter(|r| r.get_read_position_1() < 70 && matches!(r.get_variant_type(), VariantType::Insertion | VariantType::FusionGene))
            .map(|r| (r.get_variant_type().clone(), r.get_read_position_1(), r.get_read_position_2()))
            .collect()
    };
    assert_eq!(clipped(&original), vec![(VariantType::Insertion, 0, 61)]);
    let chain = original.iter().find(|m| m.get_read_id() == ids[4]).unwrap().get_splice_junctions_key();
    let cluster = |id, reads: &[usize], junctions: Vec<SpliceJunction>| {
        SpliceJunctionCluster::new(id, junctions, reads.iter().copied().collect(), vec![])
    };
    let fusion_and_clip = || vec![cluster(0, &ids[..4], chain.clone()), cluster(1, &ids[4..], chain.clone())];

    // Four aligned reads clear a floor of 3: the clip becomes the fusion junction between ASPA's
    // last base (read position 61) and WSCD1's first (65), and the read joins the fusion.
    let mut models = original.clone();
    let mut clusters = fusion_and_clip();
    assert_eq!(resolve_breakpoint_clips(&mut models, &mut clusters, &chromosomes, &fasta, &options.calling, &options.filtering, |_| 3), 1);
    assert_eq!(clipped(&models), vec![(VariantType::FusionGene, 61, 65)]);
    assert!(clusters[0].read_ids.contains(&ids[4]));
    assert!(clusters[1].read_ids.is_empty());

    // A floor of 5: the four aligned reads support nothing, and a clip never seeds a junction.
    let mut models = original.clone();
    assert_eq!(resolve_breakpoint_clips(&mut models, &mut fusion_and_clip(), &chromosomes, &fasta, &options.calling, &options.filtering, |_| 5), 0);
    assert_eq!(clipped(&models), clipped(&original));

    // No cluster supporting the fusion holds the read's chain: the clip is still the fusion's
    // arm, and the read stays where its chain put it.
    let mut models = original.clone();
    let mut clusters = fusion_and_clip();
    clusters[0].splice_junctions.clear();
    assert_eq!(resolve_breakpoint_clips(&mut models, &mut clusters, &chromosomes, &fasta, &options.calling, &options.filtering, |_| 3), 1);
    assert_eq!(clipped(&models), vec![(VariantType::FusionGene, 61, 65)]);
    assert!(clusters[1].read_ids.contains(&ids[4]));

    // Two clusters support the fusion and hold the chain: one junction, so the clip is resolved,
    // but no cluster is chosen for the read.
    options.filtering.min_reads = 1;
    options.filtering.min_total_depth = 1;
    let mut models = original.clone();
    let mut clusters = vec![
        cluster(0, &ids[..2], chain.clone()),
        cluster(1, &ids[2..4], chain.clone()),
        cluster(2, &ids[4..], chain),
    ];
    assert_eq!(resolve_breakpoint_clips(&mut models, &mut clusters, &chromosomes, &fasta, &options.calling, &options.filtering, |_| 1), 1);
    assert_eq!(clipped(&models), vec![(VariantType::FusionGene, 61, 65)]);
    assert!(clusters[2].read_ids.contains(&ids[4]));
}
