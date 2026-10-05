use bimap::BiMap;
use exacto_caller::prelude::{GraphOperation, GraphOperationType, VariantCall, VariantRecord, VariantType};
use exacto_cluster::prelude::{RNAReadCluster, RNAReadClusterSet};
use exacto_core::prelude::*;
use noodles_bam as bam;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;
use tempfile::tempdir;

use super::*;


/// The corrected FASTQ has to be readable *and indexable* by the next stage.
///
/// `determine-rna-consensus` resolves sequences by seeking to a read's offset, and
/// `build_fastq_offset_index` panics outright on plain gzip because a single deflate stream
/// cannot be seeked. So writing the output with `flate2` would produce a `.fastq.gz` that opens
/// fine in every other tool and then breaks the one stage that consumes it. Indexing it here is
/// what pins the output to BGZF.
#[test]
fn test_corrected_fastq_is_indexable_by_the_consensus_stage() {
    let directory = tempdir().unwrap();
    let fastq_file = directory.path().join("corrected.fastq.gz");
    let fastq_file: &str = fastq_file.to_str().unwrap();

    let records: Vec<(Box<str>, Box<str>, Box<str>)> = vec![
        ("read-a/1/ccs".into(), "ACGTACGT".into(), "IIIIIIII".into()),
        ("read-b/1/ccs".into(), "ACGTACGT".into(), "IIIIIIII".into()),
        ("read-c/1/ccs".into(), "TTTT".into(), "IIII".into())
    ];
    write_fastq_file(&records, fastq_file);

    // Plain gzip would satisfy this but fail everything below it.
    assert!(is_bgzipped(fastq_file), "corrected FASTQ is not BGZF-compressed.");

    let index = build_fastq_offset_index(fastq_file);
    assert_eq!(index.len(), 3);

    let mut reader: FastqReader = FastqReader::open(fastq_file);
    let sequence: Box<str> = reader.get_read_sequence(
        "read-c/1/ccs",
        *index.get("read-c/1/ccs").unwrap()
    );
    assert_eq!(&*sequence, "TTTT");
}


/// The cluster set numbers chromosomes over the contigs its tables name; the BAM numbers them
/// by header order. Correction compares calls against the model, which uses BAM ids, so a call
/// keyed on the wrong id never matches and is neither honoured nor applied. Single-contig
/// fixtures give id 0 on both sides and cannot see this, so the maps here disagree on purpose.
#[test]
fn test_cluster_operations_are_translated_to_bam_chromosome_ids() {
    // Cluster set: chr17 -> 0, chr2 -> 1, chrUn -> 2 (chrUn is not in the BAM).
    let mut cluster_chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    cluster_chromosome_names_map.insert("chr17".into(), 0u16);
    cluster_chromosome_names_map.insert("chr2".into(), 1u16);
    cluster_chromosome_names_map.insert("chrUn".into(), 2u16);
    // BAM header: chr2 -> 0, chr17 -> 1, chrM -> 2.
    let mut bam_chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    bam_chromosome_names_map.insert("chr2".into(), 0u16);
    bam_chromosome_names_map.insert("chr17".into(), 1u16);
    bam_chromosome_names_map.insert("chrM".into(), 2u16);

    let snv = |chromosome: u16, position: u32| -> GraphOperation {
        GraphOperation::new(
            chromosome, position - 1, Strand::Forward, GraphOperationType::Downstream,
            chromosome, position + 1, Strand::Forward, GraphOperationType::Upstream,
            "T".into(), VariantType::SingleNucleotideVariant
        )
    };
    // One record per call, so the consensus is that record's own operation.
    let call = |id: usize, chromosome: u16, position: u32| -> VariantCall {
        VariantCall::from_variant_records(
            id, HashSet::from([VariantRecord::new(0, 0, 0, snv(chromosome, position))]), 2, -4, 4, 2
        )
    };
    let mut cluster_set = RNAReadClusterSet::new(BiMap::new(), cluster_chromosome_names_map);
    cluster_set.add_cluster(RNAReadCluster::new(
        7,
        HashSet::new(),
        Vec::new(),
        vec![call(0, 1, 100), call(1, 0, 200), call(2, 2, 300)],   // chr2, chr17, chrUn
        HashMap::new(),
        HashSet::new()
    ));
    cluster_set.add_cluster(RNAReadCluster::new(
        8, HashSet::new(), Vec::new(), Vec::new(), HashMap::new(), HashSet::new()
    ));

    let translated: HashMap<usize, HashSet<GraphOperation>> =
        translate_cluster_operations(&cluster_set, &bam_chromosome_names_map);

    // chr2 is cluster id 1 but BAM id 0; chr17 the reverse. chrUn has nowhere to go.
    assert_eq!(translated.len(), 2);
    assert_eq!(translated[&7], HashSet::from([snv(0, 100), snv(1, 200)]));
    assert!(translated[&8].is_empty(), "a cluster without calls must still get an entry");
}


/// Every clustered read is written exactly once, in `(cluster, read name)` order, with a
/// non-empty sequence and a quality string of the same length. The FASTQ is the only output;
/// it is what the next stage reads.
#[test]
fn test_correction_emits_every_clustered_read_once_in_cluster_order() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let bam_file = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_file: &str = bam_file.to_str().unwrap();

    let mut read_names: Vec<Box<str>> = index_bam_records(bam_file, true, 1).1.left_values().cloned().collect();
    read_names.sort();
    assert!(read_names.len() >= 3, "fixture must hold at least three reads to make two ordered clusters");

    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), 0u16);
    chromosome_names_map.insert("chr18".into(), 1u16);
    let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    for (read_id, read_name) in read_names.iter().enumerate() {
        read_names_map.insert(read_name.clone(), read_id);
    }

    // The first read by name goes to the higher cluster id, so cluster order is not name order.
    let mut cluster_set = RNAReadClusterSet::new(read_names_map.clone(), chromosome_names_map);
    cluster_set.add_cluster(RNAReadCluster::new(
        1, HashSet::from([0usize]), Vec::new(), Vec::new(), HashMap::new(), HashSet::new()
    ));
    cluster_set.add_cluster(RNAReadCluster::new(
        0, (1..read_names.len()).collect(), Vec::new(), Vec::new(), HashMap::new(), HashSet::new()
    ));

    let directory = tempdir().unwrap();
    let output_fastq_file = directory.path().join("corrected.fastq.gz");
    let output_fastq_file: &str = output_fastq_file.to_str().unwrap();
    correct_rna_reads(bam_file, output_fastq_file, &cluster_set, None, &options, 2, 1000).unwrap();

    // Cluster 0's reads by name, then cluster 1's single read: not name order.
    let expected_order: Vec<Box<str>> = read_names[1..]
        .iter()
        .cloned()
        .chain([read_names[0].clone()])
        .collect();
    let written: Vec<(Box<str>, Box<str>, Vec<u8>)> = open_fastq_reader(output_fastq_file)
        .records()
        .map(|record| {
            let record = record.unwrap();
            (
                std::str::from_utf8(record.name()).unwrap().into(),
                std::str::from_utf8(record.sequence()).unwrap().into(),
                record.quality_scores().iter().map(|&q| q - b'!').collect()
            )
        })
        .collect();
    let written_order: Vec<Box<str>> = written.iter().map(|(read_name, _, _)| read_name.clone()).collect();
    assert_eq!(written_order, expected_order, "FASTQ order is not (cluster, read name) order");
    for (read_name, sequence, quality) in written.iter() {
        assert!(!sequence.is_empty(), "read {read_name} was emitted empty");
        assert_eq!(sequence.len(), quality.len(), "read {read_name}: quality string fell out of step");
    }
}


/// A read split across two records is corrected across *both*, its fusion arm survives, and the
/// called-variants list is consulted on that path.
///
/// Read 4 of `scga-mini-rna-015` is 1814 bp, with a primary at chr17:3,476,144 (`cs :230-t:21…`,
/// then `1165S`) and a supplementary at chr17:6,087,988. Reading the primary alone leaves those
/// 1165 bases as clip — untouchable, and never corrected. Three things have to hold at once:
///
/// * the supplementary's small events (three one-base insertions inside the primary's clip) are
///   applied, so the clipped bases change;
/// * its `475I` run is left alone — that is a chained fusion arm, and correcting an insertion means
///   dropping those bases out of the read, so a wrong predicate here deletes 475 bp of real
///   sequence and nothing downstream can recover it;
/// * a deletion the cluster called is kept as a gap rather than filled.
#[test]
fn test_correction_reaches_a_supplementary_record_without_draining_its_fusion_arm() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let bam_file = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-015-tumor_minimap2_sorted.bam");
    let bam_file: &str = bam_file.to_str().unwrap();

    let read_name: Box<str> = "scga-mini-rna-015-tumor_chunk_0000/4/ccs".into();
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), 0u16);
    let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    read_names_map.insert(read_name.clone(), 0usize);

    // Each variant as a one-record call, so the consensus is that record's own operation.
    // Correct into a temp FASTQ and read the one record back.
    let run = |variants: Vec<GraphOperation>| -> (Box<str>, Box<str>) {
        let variant_calls: Vec<VariantCall> = variants
            .into_iter()
            .enumerate()
            .map(|(id, variant)| VariantCall::from_variant_records(
                id, HashSet::from([VariantRecord::new(0, 0, 0, variant)]), 2, -4, 4, 2
            ))
            .collect();
        let mut cluster_set = RNAReadClusterSet::new(read_names_map.clone(), chromosome_names_map.clone());
        cluster_set.add_cluster(RNAReadCluster::new(
            0, HashSet::from([0usize]), Vec::new(), variant_calls, HashMap::new(), HashSet::new()
        ));
        let directory = tempdir().unwrap();
        let output_fastq_file = directory.path().join("corrected.fastq.gz");
        let output_fastq_file: &str = output_fastq_file.to_str().unwrap();
        correct_rna_reads(bam_file, output_fastq_file, &cluster_set, None, &options, 1, 1000).unwrap();
        let mut reader = open_fastq_reader(output_fastq_file);
        let mut records = reader.records();
        let record = records
            .next()
            .expect("correction must emit every read the clusters table names")
            .unwrap();
        assert!(records.next().is_none(), "one clustered read must give one corrected read");
        assert_eq!(record.name(), read_name.as_bytes());
        (
            std::str::from_utf8(record.sequence()).unwrap().into(),
            std::str::from_utf8(record.quality_scores()).unwrap().into()
        )
    };
    let (sequence, quality): (Box<str>, Box<str>) = run(Vec::new());

    // The original read is 1814 bp. Dropping the 475 bp arm is the failure this test exists to
    // catch, and it would be plainly visible in the length.
    let length: usize = sequence.len();
    assert!(
        length >= 1814 - options.max_correctable_event_len,
        "read collapsed from 1814 to {length}: a structural-scale event was dropped"
    );
    assert!(length <= 1814 + options.max_correctable_event_len, "read grew implausibly, to {length}");
    assert_eq!(sequence.len(), quality.len(), "quality string fell out of step");

    // The read's last 1165 bases are soft-clipped in the primary, so a primary-only correction
    // could not have touched them. A change there proves the multi-record path is live.
    let (record_positions_map, bam_read_names_map) = index_bam_records(bam_file, true, 1);
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    reader.read_header().unwrap();
    let records: Vec<bam::Record> = fetch_bam_records_for_read_id(
        &mut reader, *bam_read_names_map.get_by_left(&read_name).unwrap(), &record_positions_map
    );
    let original: Box<str> = get_bam_fastx_read_sequence(&records);
    assert_ne!(
        &sequence[sequence.len() - 1165..], &original[original.len() - 1165..],
        "the primary's clipped bases are unchanged: the supplementary was never reached"
    );

    // The primary's `-t` after 230 matched bases from 3,476,144 deletes reference base 3,476,374.
    // Called as a trusted deletion, that gap is kept instead of filled: the output loses exactly
    // the one restored base and nothing else changes.
    let trusted_deletion: GraphOperation = GraphOperation::new(
        0, 3_476_373, Strand::Forward, GraphOperationType::Downstream,
        0, 3_476_375, Strand::Forward, GraphOperationType::Upstream,
        "".into(), VariantType::Deletion
    );
    let (kept, _): (Box<str>, Box<str>) = run(vec![trusted_deletion]);
    assert_eq!(kept.len() + 1, sequence.len(), "the trusted deletion was not honoured");
    let differs_at: usize = sequence
        .bytes()
        .zip(kept.bytes())
        .position(|(filled, kept)| filled != kept)
        .unwrap_or(kept.len());
    assert_eq!(
        format!("{}{}", &sequence[..differs_at], &sequence[differs_at + 1..]), &*kept,
        "more than the trusted deletion changed"
    );
}


/// Correction must be a pure function of its inputs — same BAM and cluster tables, same
/// corrected BYTES — across runs.
///
/// Two in-process calls are a real perturbation: `RandomState` reseeds for every `HashMap`
/// correction builds internally. rna-001 is the multi-read fixture (a single-read fixture cannot
/// perturb iteration order); reads are pooled into one cluster with no variant list, which
/// maximizes the corrected surface and therefore the bytes at stake.
#[test]
fn test_correction_bytes_are_deterministic_across_runs() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let bam_file = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_file: &str = bam_file.to_str().unwrap();

    let read_names: HashSet<Box<str>> = index_bam_records(bam_file, true, 1).1.left_values().cloned().collect();
    assert!(read_names.len() > 1, "fixture must be multi-read to perturb iteration order");

    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), 0u16);
    chromosome_names_map.insert("chr18".into(), 1u16);
    let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    for (read_id, read_name) in read_names.iter().enumerate() {
        read_names_map.insert(read_name.clone(), read_id);
    }

    let run = || -> HashMap<Box<str>, (Box<str>, Box<str>)> {
        let mut cluster_set = RNAReadClusterSet::new(
            read_names_map.clone(),
            chromosome_names_map.clone()
        );
        cluster_set.add_cluster(RNAReadCluster::new(
            0,
            read_names_map.right_values().copied().collect(),
            Vec::new(),
            Vec::new(),
            HashMap::new(),
            HashSet::new()
        ));
        let directory = tempdir().unwrap();
        let output_fastq_file = directory.path().join("corrected.fastq.gz");
        let output_fastq_file: &str = output_fastq_file.to_str().unwrap();
        correct_rna_reads(bam_file, output_fastq_file, &cluster_set, None, &options, 2, 1000).unwrap();
        open_fastq_reader(output_fastq_file)
            .records()
            .map(|record| {
                let record = record.unwrap();
                (
                    std::str::from_utf8(record.name()).unwrap().into(),
                    (
                        std::str::from_utf8(record.sequence()).unwrap().into(),
                        std::str::from_utf8(record.quality_scores()).unwrap().into()
                    )
                )
            })
            .collect()
    };

    let (first, second) = (run(), run());
    assert!(!first.is_empty(), "correction emitted nothing; the assertion below would be vacuous");
    assert_eq!(
        first, second,
        "two correction runs produced different corrected bytes for identical input"
    );
}


/// A read shared between two clusters is written once, in the cluster it was assigned to: a
/// second copy in the FASTQ would be a second read to every later stage. The second read by name
/// is assigned to cluster 0 and shared with cluster 1, so it is written among cluster 0's reads
/// and nowhere else.
#[test]
fn test_correction_emits_a_shared_read_once_in_its_assigned_cluster() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let bam_file = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_file: &str = bam_file.to_str().unwrap();

    let mut read_names: Vec<Box<str>> = index_bam_records(bam_file, true, 1).1.left_values().cloned().collect();
    read_names.sort();

    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), 0u16);
    chromosome_names_map.insert("chr18".into(), 1u16);
    let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    for (read_id, read_name) in read_names.iter().enumerate() {
        read_names_map.insert(read_name.clone(), read_id);
    }

    let mut cluster_set = RNAReadClusterSet::new(read_names_map, chromosome_names_map);
    cluster_set.add_cluster(RNAReadCluster::new(
        0, HashSet::from([0usize, 1]), Vec::new(), Vec::new(), HashMap::new(), HashSet::new()
    ));
    cluster_set.add_cluster(
        RNAReadCluster::new(1, (1..read_names.len()).collect(), Vec::new(), Vec::new(), HashMap::new(), HashSet::new())
            .with_shared_read_ids(HashSet::from([1usize]))
    );

    let directory = tempdir().unwrap();
    let output_fastq_file = directory.path().join("corrected.fastq.gz");
    let output_fastq_file: &str = output_fastq_file.to_str().unwrap();
    correct_rna_reads(bam_file, output_fastq_file, &cluster_set, None, &options, 2, 1000).unwrap();

    let written: Vec<Box<str>> = open_fastq_reader(output_fastq_file)
        .records()
        .map(|record| std::str::from_utf8(record.unwrap().name()).unwrap().into())
        .collect();
    assert_eq!(written, read_names, "each read once, the shared read among cluster 0's");
}


/// Correction works one chunk of whole clusters at a time. Two clusters over rna-001 with a
/// chunk size of one force one chunk per cluster; the bytes, and their order, must match a
/// single-chunk run.
#[test]
fn test_correction_bytes_are_independent_of_chunking() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let bam_file = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_file: &str = bam_file.to_str().unwrap();

    let mut read_names: Vec<Box<str>> = index_bam_records(bam_file, true, 1).1.left_values().cloned().collect();
    read_names.sort();
    assert!(read_names.len() >= 2, "fixture must hold at least two reads to make two clusters");

    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), 0u16);
    chromosome_names_map.insert("chr18".into(), 1u16);
    let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    for (read_id, read_name) in read_names.iter().enumerate() {
        read_names_map.insert(read_name.clone(), read_id);
    }

    let run = |chunk_size: usize| -> Vec<(Box<str>, Box<str>, Box<str>)> {
        let mut cluster_set = RNAReadClusterSet::new(
            read_names_map.clone(),
            chromosome_names_map.clone()
        );
        // The first read alone, the rest together.
        cluster_set.add_cluster(RNAReadCluster::new(
            0, HashSet::from([0usize]), Vec::new(), Vec::new(), HashMap::new(), HashSet::new()
        ));
        cluster_set.add_cluster(RNAReadCluster::new(
            1, (1..read_names.len()).collect(), Vec::new(), Vec::new(), HashMap::new(), HashSet::new()
        ));
        let directory = tempdir().unwrap();
        let output_fastq_file = directory.path().join("corrected.fastq.gz");
        let output_fastq_file: &str = output_fastq_file.to_str().unwrap();
        correct_rna_reads(bam_file, output_fastq_file, &cluster_set, None, &options, 2, chunk_size).unwrap();
        open_fastq_reader(output_fastq_file)
            .records()
            .map(|record| {
                let record = record.unwrap();
                (
                    std::str::from_utf8(record.name()).unwrap().into(),
                    std::str::from_utf8(record.sequence()).unwrap().into(),
                    std::str::from_utf8(record.quality_scores()).unwrap().into()
                )
            })
            .collect()
    };

    let (one_chunk, one_chunk_per_cluster) = (run(1000), run(1));
    assert_eq!(one_chunk.len(), read_names.len(), "every clustered read must be emitted");
    assert_eq!(one_chunk, one_chunk_per_cluster, "corrected bytes changed with the chunking");
}


/// A read the clusters name but the BAM lacks is an input error, found before any read is
/// corrected: the run returns it with the count and a name, and writes no output, not even a
/// partial one.
#[test]
fn test_a_read_missing_from_the_bam_stops_the_run_before_any_output() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions::default();
    let bam_file = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_file: &str = bam_file.to_str().unwrap();

    let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    for (read_id, read_name) in index_bam_records(bam_file, true, 1).1.left_values().enumerate() {
        read_names_map.insert(read_name.clone(), read_id);
    }
    let num_reads: usize = read_names_map.len() + 1;
    read_names_map.insert("not-in-the-bam".into(), num_reads - 1);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr17".into(), 0u16);
    let mut cluster_set = RNAReadClusterSet::new(read_names_map, chromosome_names_map);
    cluster_set.add_cluster(RNAReadCluster::new(
        0, (0..num_reads).collect(), Vec::new(), Vec::new(), HashMap::new(), HashSet::new()
    ));

    let directory = tempdir().unwrap();
    let output_fastq_file = directory.path().join("corrected.fastq.gz");
    let output_fastq_file: &str = output_fastq_file.to_str().unwrap();
    let result = correct_rna_reads(bam_file, output_fastq_file, &cluster_set, None, &options, 2, 1000);
    match result {
        Err(CorrectionError::ReadsNotInBam { num_missing, num_reads: total, read_name, .. }) => {
            assert_eq!((num_missing, total, &*read_name), (1, num_reads, "not-in-the-bam"));
        },
        other => panic!("expected ReadsNotInBam, got {other:?}")
    }
    assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 0, "a failed run left a file behind");
}


/// A BAM stored without base qualities (SAM `*`, what minimap2 writes for FASTA input) is
/// corrected like any other: kept bases get the Q60 the alignment model gives them, a rewritten
/// base the requested quality.
#[test]
fn test_a_bam_without_base_qualities_is_corrected() {
    use noodles_sam as sam;
    use sam::alignment::io::Write;

    let directory = tempdir().unwrap();
    let bam_file = directory.path().join("no_qualities.bam");
    let sam_text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nread-1\t0\tchr1\t101\t60\t6M\t*\t0\t0\tACGTAC\t*\tcs:Z::2*ta:3\n";
    let mut reader = sam::io::Reader::new(sam_text.as_bytes());
    let header = reader.read_header().unwrap();
    let record = reader.record_bufs(&header).next().unwrap().unwrap();
    let mut writer = bam::io::Writer::new(File::create(&bam_file).unwrap());
    writer.write_header(&header).unwrap();
    writer.write_alignment_record(&header, &record).unwrap();
    writer.try_finish().unwrap();
    drop(writer);
    let bam_file: &str = bam_file.to_str().unwrap();

    let mut read_names_map: BiMap<Box<str>, usize> = BiMap::new();
    read_names_map.insert("read-1".into(), 0usize);
    let mut chromosome_names_map: BiMap<Box<str>, u16> = BiMap::new();
    chromosome_names_map.insert("chr1".into(), 0u16);
    let mut cluster_set = RNAReadClusterSet::new(read_names_map, chromosome_names_map);
    cluster_set.add_cluster(RNAReadCluster::new(
        0, HashSet::from([0usize]), Vec::new(), Vec::new(), HashMap::new(), HashSet::new()
    ));

    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions { corrected_base_quality: 30, ..CorrectRNAReadsOptions::default() };
    let output_fastq_file = directory.path().join("corrected.fastq.gz");
    let output_fastq_file: &str = output_fastq_file.to_str().unwrap();
    correct_rna_reads(bam_file, output_fastq_file, &cluster_set, None, &options, 1, 1000).unwrap();

    let records: Vec<(Vec<u8>, Vec<u8>)> = open_fastq_reader(output_fastq_file)
        .records()
        .map(|record| {
            let record = record.unwrap();
            (record.sequence().to_vec(), record.quality_scores().iter().map(|&q| q - b'!').collect())
        })
        .collect();
    assert_eq!(records, vec![(b"ACTTAC".to_vec(), vec![60, 60, 30, 60, 60, 60])]);
}


/// Trimming transcript ends reads the gene annotation, so asking for it without one is an
/// error, returned before the BAM is read.
#[test]
fn test_trimming_transcript_ends_without_an_annotation_is_an_error() {
    let options: CorrectRNAReadsOptions = CorrectRNAReadsOptions { trim_transcript_ends: true, ..CorrectRNAReadsOptions::default() };
    let cluster_set = RNAReadClusterSet::new(BiMap::new(), BiMap::new());
    let directory = tempdir().unwrap();
    let output_fastq_file = directory.path().join("corrected.fastq.gz");
    let result = correct_rna_reads(
        "no-such.bam", output_fastq_file.to_str().unwrap(), &cluster_set, None, &options, 1, 1000
    );
    assert!(matches!(result, Err(CorrectionError::NoGeneAnnotation)), "got {result:?}");
}
