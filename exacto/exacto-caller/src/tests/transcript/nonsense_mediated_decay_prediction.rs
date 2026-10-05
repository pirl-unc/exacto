use bimap::BiMap;
use exacto_core::prelude::*;
use noodles_bam as bam;
use noodles_bam::bai;
use noodles_bam::bai::Index;
use noodles_sam::Header;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::sync::Arc;

use super::*;


#[test]
fn scga_mini_rna_001_predict_nonsense_mediated_decay_returns_matches() {
    use std::collections::HashSet;
    use NonsenseMediatedDecayVerdict::*;

    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bam_bai_path = Path::new(env!("EXACTO_TEST_DATA")).join("simulation/ground_truth/scga-mini-rna-001-tumor_minimap2_sorted.bam.bai");
    let bam_bai_full_path = fs::canonicalize(bam_bai_path).unwrap();
    let bam_bai_file: &str = bam_bai_full_path.to_str().unwrap();
    let reference_genome_fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let reference_genome_fasta_full_path = fs::canonicalize(reference_genome_fasta_path).unwrap();
    let reference_genome_fasta_file: &str = reference_genome_fasta_full_path.to_str().unwrap();
    let gencode_gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gencode_gtf_full_path = fs::canonicalize(gencode_gtf_path).unwrap();
    let gencode_gtf_file: &str = gencode_gtf_full_path.to_str().unwrap();

    let fasta_map: FastaMap = FastaMap::new(reference_genome_fasta_file);
    let chromosome_lengths: HashMap<Box<str>, u32> = get_chromosome_lengths(bam_file);
    let chromosome_names_map: BiMap<Box<str>, u16> = create_chromosome_names_map(bam_file);
    let end: u32 = *chromosome_lengths.get("chr17").unwrap();

    let gene_annotator = Gencode::new_with_defaults(
        gencode_gtf_file,
        "hg38",
        "v41"
    );

    let (record_positions_map, read_names_map) = index_bam_records(
        bam_file,
        true,
        2
    );

    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(bam_file)
        .unwrap();
    let header: Header = reader.read_header().unwrap();
    let index: Index = bai::fs::read(bam_bai_file).unwrap();

    let records_map: HashMap<usize, Vec<bam::Record>> = fetch_bam_records(
        &mut reader,
        &header,
        &index,
        "chr17".into(),
        1,
        end,
        &record_positions_map,
        &read_names_map,
        7,
        1
    );

    // The simulated TP53 (ENST00000269305.9, minus strand) transcript carrying the SNV at
    // chr17:7674225, error-free: 2,512 bases in 11 exons. Its last junction follows read position
    // 1241, the end of the 107 bp exon (114 + 102 + 22 + 279 + 184 + 113 + 110 + 137 + 74 + 107 =
    // 1242 bases).
    let read_name: &str = "scga-mini-rna-001-tumor-1";
    let read_id: usize = *read_names_map.get_by_left(read_name).unwrap();
    let read_sequence: Box<str> = get_bam_fastx_read_sequence(records_map.get(&read_id).unwrap());
    let quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(records_map.get(&read_id).unwrap());

    let alignment_model: AlignmentModel = AlignmentModel::new(
        read_id,
        &*read_sequence,
        &quality_scores,
        &records_map.get(&read_id).unwrap().iter().map(|record| Arc::new(record.clone())).collect::<Vec<_>>()
    );

    let reference_transcript_matches: Vec<ReferenceTranscriptMatch> = identify_reference_transcript_matches(
        &identify_transcript_model_exons(&alignment_model),
        &identify_transcript_model_splice_junctions(&alignment_model),
        &gene_annotator,
        &chromosome_names_map
    );

    let transcript_model: TranscriptModel = TranscriptModel::new(
        alignment_model,
        reference_transcript_matches,
        &gene_annotator,
        &chromosome_names_map,
        &fasta_map
    );

    let last_junction: Option<u32> = transcript_model
        .get_splice_junctions()
        .iter()
        .map(|junction| junction.read_position_1)
        .max();
    assert_eq!(transcript_model.get_alignment_model().num_bases(), 2512);
    assert_eq!(transcript_model.get_splice_junctions().len(), 10);
    assert_eq!(last_junction, Some(1241));
    assert!(transcript_model.get_nmd_predictions().is_empty());

    // The longest frame is the TP53 CDS. It starts 28 bases into the 102 bp exon, at read 142, and
    // its 394 codons end at read 1323, in the last exon: no junction follows its stop codon.
    let start_codons: HashSet<&str> = HashSet::from(["AUG"]);
    let calls: Vec<NonsenseMediatedDecayCall> = NonsenseMediatedDecayPredictor {
        translation_strategy: &TranslationStrategy::LongestORF,
        start_codons: &start_codons,
        distance_threshold: 50
    }.predict(&transcript_model);
    assert_eq!(
        calls,
        vec![NonsenseMediatedDecayCall { orf_start: 142, orf_end: 1323, verdict: NotPredicted { distance_to_last_junction: None } }]
    );

    // Every AUG that reaches a stop codon, in read order: 35 frames. Nine stop more than 50 bases
    // before the last junction. The other 26, the CDS among them, stop in the last exon.
    let calls: Vec<NonsenseMediatedDecayCall> = NonsenseMediatedDecayPredictor {
        translation_strategy: &TranslationStrategy::AllORFs,
        start_codons: &start_codons,
        distance_threshold: 50
    }.predict(&transcript_model);
    let mut predicted: Vec<(u32, u32, u32)> = Vec::new();
    let mut not_predicted: Vec<(u32, u32, Option<u32>)> = Vec::new();
    for call in calls.iter() {
        match call.verdict {
            Predicted { distance_to_last_junction } => predicted.push((call.orf_start, call.orf_end, distance_to_last_junction)),
            NotPredicted { distance_to_last_junction } => not_predicted.push((call.orf_start, call.orf_end, distance_to_last_junction))
        }
    }
    assert_eq!(
        predicted,
        vec![
            (207, 224, 1017),
            (263, 271, 970),
            (297, 308, 933),
            (323, 508, 733),
            (677, 1174, 67),
            (698, 1174, 67),
            (761, 1174, 67),
            (800, 1174, 67),
            (1112, 1174, 67)
        ]
    );
    assert_eq!(not_predicted.len(), 26);
    assert!(not_predicted.iter().all(|(_, orf_end, distance_to_last_junction)| *orf_end > 1241 && distance_to_last_junction.is_none()));
    assert!(not_predicted.contains(&(142, 1323, None)));
    // No frame is longer than the CDS, the one the longest-frame strategy returns.
    assert!(calls.iter().all(|call| call.orf_start < call.orf_end));
    assert_eq!(calls.iter().map(|call| call.orf_end - call.orf_start).max(), Some(1323 - 142));

    // The five frames that stop at read 1174 lie 67 bases before the last junction. Decay is
    // predicted only beyond the threshold, so at a threshold of 67 they keep their distance.
    let calls: Vec<NonsenseMediatedDecayCall> = NonsenseMediatedDecayPredictor {
        translation_strategy: &TranslationStrategy::AllORFs,
        start_codons: &start_codons,
        distance_threshold: 67
    }.predict(&transcript_model);
    let at_threshold: Vec<(u32, u32)> = calls
        .iter()
        .filter(|call| call.verdict == NotPredicted { distance_to_last_junction: Some(67) })
        .map(|call| (call.orf_start, call.orf_end))
        .collect();
    assert_eq!(calls.len(), 35);
    assert_eq!(at_threshold, vec![(677, 1174), (698, 1174), (761, 1174), (800, 1174), (1112, 1174)]);
    assert_eq!(calls.iter().filter(|call| matches!(call.verdict, Predicted { .. })).count(), 4);
}

#[test]
fn predict_nonsense_mediated_decay_needs_a_stop_codon_more_than_the_threshold_before_the_last_junction() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;
    use std::collections::HashSet;
    use NonsenseMediatedDecayVerdict::*;

    // Gene GN (+) has three 200-base exons, chr1:101-300, 401-600 and 701-900. Each read is 600
    // bases of GCC repeated (no A or T) around one frame, ATG (GCC x 8) TAA: the only AUG and the
    // only stop codon. Spliced, the read's exons are read positions 0-199, 200-399 and 400-599, so
    // the last junction follows read position 399.
    let filler: String = "GCC".repeat(200);
    let read_with_stop_codon_ending_at = |orf_end: usize| -> String {
        format!("{}ATG{}TAA{}", &filler[..orf_end - 29], "GCC".repeat(8), &filler[orf_end + 1..])
    };
    let cases: Vec<(&str, bool, String, Vec<NonsenseMediatedDecayCall>)> = vec![
        (
            "51 bases before the last junction", true, read_with_stop_codon_ending_at(348),
            vec![NonsenseMediatedDecayCall { orf_start: 319, orf_end: 348, verdict: Predicted { distance_to_last_junction: 51 } }]
        ),
        (
            "50 bases before the last junction", true, read_with_stop_codon_ending_at(349),
            vec![NonsenseMediatedDecayCall { orf_start: 320, orf_end: 349, verdict: NotPredicted { distance_to_last_junction: Some(50) } }]
        ),
        (
            "at the last junction", true, read_with_stop_codon_ending_at(399),
            vec![NonsenseMediatedDecayCall { orf_start: 370, orf_end: 399, verdict: NotPredicted { distance_to_last_junction: Some(0) } }]
        ),
        (
            "across the last junction", true, read_with_stop_codon_ending_at(400),
            vec![NonsenseMediatedDecayCall { orf_start: 371, orf_end: 400, verdict: NotPredicted { distance_to_last_junction: None } }]
        ),
        (
            "in the last exon", true, read_with_stop_codon_ending_at(402),
            vec![NonsenseMediatedDecayCall { orf_start: 373, orf_end: 402, verdict: NotPredicted { distance_to_last_junction: None } }]
        ),
        // Measured to the last junction, not to the first one, 49 bases downstream.
        (
            "before both junctions", true, read_with_stop_codon_ending_at(150),
            vec![NonsenseMediatedDecayCall { orf_start: 121, orf_end: 150, verdict: Predicted { distance_to_last_junction: 249 } }]
        ),
        (
            "unspliced", false, read_with_stop_codon_ending_at(150),
            vec![NonsenseMediatedDecayCall { orf_start: 121, orf_end: 150, verdict: NotPredicted { distance_to_last_junction: None } }]
        ),
        // A frame that runs off the 3' end is not a complete ORF.
        ("no stop codon", true, format!("{}ATG{}", &filler[..100], &filler[103..]), Vec::new())
    ];

    let temp_dir = tempfile::tempdir().unwrap();
    let gtf_file = temp_dir.path().join("genes.gtf");
    let gene_attributes: &str = "gene_id \"GN\"; gene_type \"protein_coding\"; gene_name \"GN\"; level 2;";
    let transcript_attributes: &str = "gene_id \"GN\"; transcript_id \"TGN\"; gene_type \"protein_coding\"; gene_name \"GN\"; transcript_type \"protein_coding\"; transcript_name \"TGN\"; level 2;";
    let mut gtf: String = format!(
        "chr1\tHAVANA\tgene\t101\t900\t.\t+\t.\t{gene_attributes}\nchr1\tHAVANA\ttranscript\t101\t900\t.\t+\t.\t{transcript_attributes}\n"
    );
    for (number, (exon_start, exon_end)) in [(101, 300), (401, 600), (701, 900)].iter().enumerate() {
        gtf.push_str(&format!(
            "chr1\tHAVANA\texon\t{exon_start}\t{exon_end}\t.\t+\t.\t{transcript_attributes} exon_number {}; exon_id \"EGN{}\";\n",
            number + 1,
            number + 1
        ));
    }
    fs::write(&gtf_file, gtf).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_file.to_str().unwrap(), "hg38", "v41");
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr1".into(), 0);

    for (i, (name, spliced, read, expected)) in cases.into_iter().enumerate() {
        // chr1 holds the read at 101, around two GT..AG introns of 100 bases when it is spliced.
        let intron: String = format!("GT{}AG", "C".repeat(96));
        let (cigar, cs, chromosome): (&str, &str, String) = if spliced {
            (
                "200=100N200=100N200=",
                ":200~gt100ag:200~gt100ag:200",
                format!("{}{}{intron}{}{intron}{}{}", "C".repeat(100), &read[..200], &read[200..400], &read[400..], "C".repeat(100))
            )
        } else {
            ("600=", ":600", format!("{}{read}{}", "C".repeat(100), "C".repeat(300)))
        };
        let fasta_file = temp_dir.path().join(format!("genome_{i}.fa"));
        fs::write(&fasta_file, format!(">chr1\n{chromosome}\n")).unwrap();
        let fasta_map: FastaMap = FastaMap::new(fasta_file.to_str().unwrap());

        let qualities: String = "I".repeat(read.len());
        let sam_text: String = format!(
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nread-1\t0\tchr1\t101\t60\t{cigar}\t*\t0\t0\t{read}\t{qualities}\tcs:Z:{cs}\n"
        );
        let mut sam_reader = sam::io::Reader::new(sam_text.as_bytes());
        let header = sam_reader.read_header().unwrap();
        let mut writer = bam::io::Writer::new(Vec::new());
        for record_buf in sam_reader.record_bufs(&header) {
            writer.write_alignment_record(&header, &record_buf.unwrap()).unwrap();
        }
        writer.try_finish().unwrap();
        let mut reader = bam::io::Reader::new(writer.get_ref().get_ref().as_slice());
        let mut records: Vec<bam::Record> = Vec::new();
        let mut record = bam::Record::default();
        while reader.read_record(&mut record).unwrap() > 0 {
            records.push(record.clone());
        }
        let read_sequence: Box<str> = get_bam_fastx_read_sequence(&records);
        let base_quality_scores: Vec<u8> = get_bam_fastx_base_quality_scores(&records);
        let bam_records: Vec<Arc<bam::Record>> = records.into_iter().map(Arc::new).collect();
        let alignment_model: AlignmentModel = AlignmentModel::new(1, &*read_sequence, &base_quality_scores, &bam_records);
        let transcript_model: TranscriptModel = TranscriptModel::new(
            alignment_model,
            Vec::new(),
            &gene_annotator,
            &chromosome_names_map,
            &fasta_map
        );

        let junctions: Vec<u32> = transcript_model.get_splice_junctions().iter().map(|junction| junction.read_position_1).collect();
        assert_eq!(junctions, if spliced { vec![199, 399] } else { Vec::new() }, "{name}");

        let calls: Vec<NonsenseMediatedDecayCall> = NonsenseMediatedDecayPredictor {
            translation_strategy: &TranslationStrategy::AllORFs,
            start_codons: &HashSet::from(["AUG"]),
            distance_threshold: 50
        }.predict(&transcript_model);
        assert_eq!(calls, expected, "{name}");
    }
}

#[test]
fn classify_stop_codon_returns_matches() {
    use NonsenseMediatedDecayVerdict::*;

    assert_eq!(classify_stop_codon(10, None, 50), NotPredicted { distance_to_last_junction: None });           // unspliced
    assert_eq!(classify_stop_codon(120, Some(100), 50), NotPredicted { distance_to_last_junction: None });     // stop in the last exon
    assert_eq!(classify_stop_codon(50, Some(100), 50), NotPredicted { distance_to_last_junction: Some(50) });  // at the threshold
    assert_eq!(classify_stop_codon(49, Some(100), 50), Predicted { distance_to_last_junction: 51 });
}
