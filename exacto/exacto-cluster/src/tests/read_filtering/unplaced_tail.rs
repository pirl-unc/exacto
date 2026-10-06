use bimap::BiMap;
use exacto_core::prelude::{create_chromosome_names_map, index_bam_records, FastaMap, Gencode, Strand};
use std::fs;
use std::path::Path;

use super::*;
use crate::prelude::ClusterRNAReadsOptions;


/// A large insertion, a splice junction, a short aligned block and a terminal soft clip are a
/// stretch the aligner left unplaced when the block is shorter than the insertion and the clip.
///
/// The first read is of a fusion transcript joining WSCD1 exon 7 to ACAP1 exon 2 (chr17, forward
/// strand). minimap2 chained 17 bases of ACAP1 exon 5 to a copy of them at 6,165,489, 54,689
/// bases past the 29th base of WSCD1 exon 7, and spelled the end of the read as
/// `29= 467I 54689N 17= 68S`. The other reads vary one part of it each.
///
///   Read length   Insertion     Junction      Block         Clip          Identified
///   1,910         1,358-1,824   1,357-1,825   1,825-1,841   1,842-1,909   1,358-1,909
///   1,910         85-551        84-552        68-84         0-67          0-551, at the start of the read
///   1,852         1,358-1,824   1,357-1,825   1,825-1,841   1,842-1,851   none, the block is longer than the clip
///   1,542         1,358-1,456   1,357-1,457   1,457-1,473   1,474-1,541   none, the insertion is under 100 bases
///   1,910         1,358-1,824   none          1,825-1,841   1,842-1,909   none, no junction at the insertion
///   1,842         1,358-1,824   1,357-1,825   1,825-1,841   none          none, the read ends aligned
///   1,910         none          1,357-1,358   1,358-1,841   1,842-1,909   none, a clip past an ordinary junction
#[test]
fn identify_unplaced_tail_returns_matches() {
    // (read length, insertion, junction, clip, unplaced tail), all in read positions
    let cases: Vec<(u32, Option<(u32, u32)>, Option<(u32, u32)>, Option<(u32, u32)>, Option<(u32, u32)>)> = vec![
        (1_910, Some((1_358, 1_824)), Some((1_357, 1_825)), Some((1_842, 1_909)), Some((1_358, 1_909))),
        (1_910, Some((85, 551)), Some((84, 552)), Some((0, 67)), Some((0, 551))),
        (1_852, Some((1_358, 1_824)), Some((1_357, 1_825)), Some((1_842, 1_851)), None),
        (1_542, Some((1_358, 1_456)), Some((1_357, 1_457)), Some((1_474, 1_541)), None),
        (1_910, Some((1_358, 1_824)), None, Some((1_842, 1_909)), None),
        (1_842, Some((1_358, 1_824)), Some((1_357, 1_825)), None, None),
        (1_910, None, Some((1_357, 1_358)), Some((1_842, 1_909)), None)
    ];
    for (read_length, insertion, junction, clip, expected) in cases {
        // The insertion sits at the end of WSCD1 exon 7 as aligned, the clip at the end of the block.
        let mut variant_records: Vec<VariantRecord> = Vec::new();
        for (read_positions, position) in [(insertion, 6_110_799u32), (clip, 6_165_505u32)] {
            let Some((read_position_1, read_position_2)) = read_positions else {
                continue;
            };
            variant_records.push(VariantRecord::new(
                0,
                read_position_1,
                read_position_2,
                GraphOperation::new(
                    0,
                    position,
                    Strand::Forward,
                    GraphOperationType::Downstream,
                    0,
                    position + 1,
                    Strand::Forward,
                    GraphOperationType::Upstream,
                    "A".repeat((read_position_2 - read_position_1 + 1) as usize).into(),
                    VariantType::Insertion
                )
            ));
        }

        // An ordinary junction of WSCD1 ahead of the stretch, and the junction into the block.
        let mut splice_junctions: Vec<TranscriptModelSpliceJunction> = vec![
            TranscriptModelSpliceJunction::new(
                0,
                0,
                6_109_767,
                6_110_770,
                Strand::Forward,
                Strand::Forward,
                1,
                700,
                701
            )
        ];
        if let Some((read_position_1, read_position_2)) = junction {
            splice_junctions.push(TranscriptModelSpliceJunction::new(
                0,
                0,
                6_110_800,
                6_165_488,
                Strand::Forward,
                Strand::Forward,
                2,
                read_position_1,
                read_position_2
            ));
        }

        assert_eq!(
            identify_unplaced_tail(&variant_records, &splice_junctions, read_length, 100),
            expected,
            "read length {}, insertion {:?}, junction {:?}, clip {:?}",
            read_length,
            insertion,
            junction,
            clip
        );
    }
}


/// The stretch is identified on a read as the transcript model spells it, and the model holds
/// its junctions numbered in read order once the junction into the block is taken out.
///
/// `chunk_0000/121/ccs` of scga-mini-rna-015 is a read of 1,618 bases of its ASPA-WSCD1-ACAP1
/// fusion, aligned in two records: ASPA on a supplementary one, and WSCD1 on the primary one,
/// which ends `29= 467I 54689N 17= 68S`. It holds 8 junctions, the last of them the junction into
/// the 17 bases, 6,110,800-6,165,488.
#[test]
fn scga_mini_rna_015_identify_unplaced_tail_returns_stretch_of_read() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-015-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_map: FastaMap = FastaMap::new(fasta_full_path.to_str().unwrap());
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let (record_positions_map, read_names_map) = index_bam_records(bam_file, true, 2);
    let read_id: usize = *read_names_map.get_by_left("scga-mini-rna-015-tumor_chunk_0000/121/ccs").unwrap();
    let options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;
    let temp_dir = tempfile::tempdir().unwrap();
    let (temp_file, summaries, _) = characterize_rna_reads(
        bam_file,
        &vec![read_id],
        &read_names_map,
        &record_positions_map,
        &chromosome_names_map,
        &fasta_map,
        &gene_annotator,
        1,
        options.calling.min_mapping_quality,
        options.calling.min_terminal_soft_clip_ins_len,
        options.calling.bkpt_rescue,
        options.calling.bkpt_rescue_min_ins_len,
        options.calling.bkpt_rescue_realignment_gap_open_score,
        options.calling.bkpt_rescue_realignment_gap_extend_score,
        options.calling.bkpt_rescue_realignment_k,
        options.calling.bkpt_rescue_realignment_band_width,
        options.calling.bkpt_rescue_realignment_min_score_fraction,
        options.calling.bkpt_rescue_realignment_min_query_coverage,
        options.calling.bkpt_rescue_realignment_min_placed_fraction,
        options.calling.bkpt_rescue_realignment_max_pieces,
        None,
        options.calling.chunk_size,
        temp_dir.path().to_str().unwrap()
    );
    let mut transcript_models: Vec<TranscriptModel> = load_transcript_models(&temp_file, summaries.iter());
    assert_eq!(transcript_models.len(), 1);
    let transcript_model: &mut TranscriptModel = &mut transcript_models[0];
    let read_length: u32 = transcript_model.get_alignment_model().num_bases();
    assert_eq!(read_length, 1_618);
    assert_eq!(transcript_model.get_splice_junctions().len(), 8);

    // The stretch runs from the first base of the insertion to the last base of the read.
    let unplaced_tail: Option<(u32, u32)> = identify_unplaced_tail(
        transcript_model.get_variant_records(),
        transcript_model.get_splice_junctions(),
        read_length,
        options.calling.bkpt_rescue_min_ins_len
    );
    assert_eq!(unplaced_tail, Some((1_066, 1_617)));

    // Without the junction into the block, the 7 junctions left keep their numbers.
    transcript_model.retain_splice_junctions(|splice_junction| {
        splice_junction.read_position_2 < 1_066 || splice_junction.read_position_1 > 1_617
    });
    let splice_junctions: Vec<(u16, u32, u32)> = transcript_model
        .get_splice_junctions()
        .iter()
        .map(|splice_junction| (
            splice_junction.splice_junction_number,
            splice_junction.reference_position_1,
            splice_junction.reference_position_2
        ))
        .collect();
    assert_eq!(
        splice_junctions,
        vec![
            (1, 3_476_396, 3_481_602),
            (2, 3_481_799, 3_483_498),
            (3, 3_483_593, 3_489_234),
            (4, 6_088_105, 6_090_320),
            (5, 6_090_506, 6_095_101),
            (6, 6_095_224, 6_109_606),
            (7, 6_109_767, 6_110_770)
        ]
    );

    // Without the first junction, the 6 junctions left are numbered from 1 again.
    transcript_model.retain_splice_junctions(|splice_junction| splice_junction.reference_position_1 != 3_476_396);
    let splice_junctions: Vec<(u16, u32, u32)> = transcript_model
        .get_splice_junctions()
        .iter()
        .map(|splice_junction| (
            splice_junction.splice_junction_number,
            splice_junction.reference_position_1,
            splice_junction.reference_position_2
        ))
        .collect();
    assert_eq!(
        splice_junctions,
        vec![
            (1, 3_481_799, 3_483_498),
            (2, 3_483_593, 3_489_234),
            (3, 6_088_105, 6_090_320),
            (4, 6_090_506, 6_095_101),
            (5, 6_095_224, 6_109_606),
            (6, 6_109_767, 6_110_770)
        ]
    );
}
