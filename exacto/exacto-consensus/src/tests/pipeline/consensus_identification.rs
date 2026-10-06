use std::io::Write;

use super::*;


/// The clusters are sorted by their number of reads, largest first, and clusters with the same
/// number of reads by their cluster ID.
///
///   Cluster ID   Reads   Place
///   4            2       3
///   9            5       1
///   1            2       2
///   7            1       4
#[test]
fn sort_clusters_by_num_reads_returns_largest_cluster_first() {
    let clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::from([
        (4, HashSet::from(["read-a".into(), "read-b".into()])),
        (9, HashSet::from(["read-c".into(), "read-d".into(), "read-e".into(), "read-f".into(), "read-g".into()])),
        (1, HashSet::from(["read-h".into(), "read-i".into()])),
        (7, HashSet::from(["read-j".into()]))
    ]);

    let sorted_clusters: Vec<(usize, &HashSet<Box<str>>)> = sort_clusters_by_num_reads(&clusters);

    let cluster_ids: Vec<usize> = sorted_clusters.iter().map(|&(cluster_id, _)| cluster_id).collect();
    assert_eq!(cluster_ids, vec![9, 1, 4, 7]);
    for (cluster_id, read_names) in sorted_clusters {
        assert_eq!(read_names, &clusters[&cluster_id]);
    }
}


/// Each cluster gets the consensus of its own reads, and the set holds the clusters by their
/// cluster ID, on one thread and on several.
///
/// The threads take the clusters from one queue, largest first, and finish them in no fixed
/// order, so the consensus sequences have to find their clusters by cluster ID.
///
///   Cluster ID   Reads   Consensus          Reads apart from the consensus
///   11           5       ACGTACGTACGTACGT   read-11-3 (A>T at 5), read-11-4 (T>A at 8)
///   3            4       TTGACCATGCTAGCTA   read-3-3 (C>G at 6)
///   20           3       GGGCCCAAATTTGGGC
///   8            2       CATGCATGGTACGTAC
///   5            1       TTTTACGTCCCCACGT
///   2            1       AGAGAGTCTCTCAGAG
///
/// The read `decoy` is in no cluster.
#[test]
fn identify_consensus_sequences_returns_consensus_of_each_cluster() {
    let temp_dir = tempfile::tempdir().unwrap();
    let fastq_path = temp_dir.path().join("reads.fastq");
    let records: Vec<(&str, &str)> = vec![
        ("read-5-0", "TTTTACGTCCCCACGT"),
        ("read-11-0", "ACGTACGTACGTACGT"),
        ("read-3-0", "TTGACCATGCTAGCTA"),
        ("read-11-3", "ACGTTCGTACGTACGT"),
        ("read-20-0", "GGGCCCAAATTTGGGC"),
        ("decoy", "TTTTTTTTTTTTTTTT"),
        ("read-8-0", "CATGCATGGTACGTAC"),
        ("read-11-1", "ACGTACGTACGTACGT"),
        ("read-3-3", "TTGACGATGCTAGCTA"),
        ("read-20-1", "GGGCCCAAATTTGGGC"),
        ("read-2-0", "AGAGAGTCTCTCAGAG"),
        ("read-3-1", "TTGACCATGCTAGCTA"),
        ("read-11-4", "ACGTACGAACGTACGT"),
        ("read-8-1", "CATGCATGGTACGTAC"),
        ("read-20-2", "GGGCCCAAATTTGGGC"),
        ("read-3-2", "TTGACCATGCTAGCTA"),
        ("read-11-2", "ACGTACGTACGTACGT")
    ];
    let mut fastq_file = std::fs::File::create(&fastq_path).unwrap();
    for (read_name, sequence) in records.iter() {
        writeln!(fastq_file, "@{}\n{}\n+\n{}", read_name, sequence, "I".repeat(sequence.len())).unwrap();
    }
    drop(fastq_file);

    let clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::from([
        (11, HashSet::from(["read-11-0".into(), "read-11-1".into(), "read-11-2".into(), "read-11-3".into(), "read-11-4".into()])),
        (3, HashSet::from(["read-3-0".into(), "read-3-1".into(), "read-3-2".into(), "read-3-3".into()])),
        (20, HashSet::from(["read-20-0".into(), "read-20-1".into(), "read-20-2".into()])),
        (8, HashSet::from(["read-8-0".into(), "read-8-1".into()])),
        (5, HashSet::from(["read-5-0".into()])),
        (2, HashSet::from(["read-2-0".into()]))
    ]);
    let options: IdentifyConsensusSequencesOptions = IdentifyConsensusSequencesOptions::DEFAULT;

    for num_threads in [1usize, 2, 4] {
        let consensus_set: ConsensusSequenceSet = identify_consensus_sequences(
            &clusters,
            fastq_path.to_str().unwrap(),
            &options,
            num_threads
        );

        // Vec<(cluster ID, consensus sequence)>
        let consensus_sequences: Vec<(usize, &str)> = consensus_set.sequences
            .iter()
            .map(|sequence| (sequence.get_cluster_id(), sequence.get_consensus_sequence()))
            .collect();
        assert_eq!(
            consensus_sequences,
            vec![
                (2, "AGAGAGTCTCTCAGAG"),
                (3, "TTGACCATGCTAGCTA"),
                (5, "TTTTACGTCCCCACGT"),
                (8, "CATGCATGGTACGTAC"),
                (11, "ACGTACGTACGTACGT"),
                (20, "GGGCCCAAATTTGGGC")
            ],
            "{} threads",
            num_threads
        );
        for sequence in consensus_set.sequences.iter() {
            assert_eq!(sequence.get_read_names(), &clusters[&sequence.get_cluster_id()], "{} threads", num_threads);
        }
    }
}


/// A cluster of more than `max_reads_per_cluster` reads gets its consensus sequence from that
/// many of its reads, drawn at random. The result is the same on one thread and on several.
///
/// Cluster 7 holds 30 reads: 24 of the transcript and 6 with one base changed each. Cluster 2
/// holds 5 reads. The read `decoy` is in no cluster.
///
///   max_reads_per_cluster   Cluster 7                         Cluster 2
///   10                      10 reads, the transcript          all reads
///   30                      all reads                         all reads
///   0                       all reads                         all reads
///   1                       1 read, its sequence              1 read, its sequence
#[test]
fn identify_consensus_sequences_subsamples_clusters_over_max_reads_per_cluster() {
    let temp_dir = tempfile::tempdir().unwrap();
    let fastq_path = temp_dir.path().join("reads.fastq");
    let mut records: Vec<(String, &str)> = Vec::new();
    for i in 0..24 {
        records.push((format!("read-7-{}", i), "ACGTACGTACGTTTGACCATGCTAGCTAGCTA"));
    }
    records.push(("read-7-24".to_string(), "ACGTTCGTACGTTTGACCATGCTAGCTAGCTA"));     // A>T at 5
    records.push(("read-7-25".to_string(), "ACGTACGAACGTTTGACCATGCTAGCTAGCTA"));     // T>A at 8
    records.push(("read-7-26".to_string(), "ACGTACGTACGATTGACCATGCTAGCTAGCTA"));     // T>A at 12
    records.push(("read-7-27".to_string(), "ACGTACGTACGTTTCACCATGCTAGCTAGCTA"));     // G>C at 15
    records.push(("read-7-28".to_string(), "ACGTACGTACGTTTGACCATGGTAGCTAGCTA"));     // C>G at 22
    records.push(("read-7-29".to_string(), "ACGTACGTACGTTTGACCATGCTAGCTAGCAA"));     // T>A at 31
    records.push(("decoy".to_string(), "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    for i in 0..4 {
        records.push((format!("read-2-{}", i), "GGGCCCAAATTTGGGCCATGCATGGTACGTAC"));
    }
    records.push(("read-2-4".to_string(), "GGGCCCAAATTTGGGCCATGCATCGTACGTAC"));      // G>C at 24
    let mut fastq_file = std::fs::File::create(&fastq_path).unwrap();
    for (read_name, sequence) in records.iter() {
        writeln!(fastq_file, "@{}\n{}\n+\n{}", read_name, sequence, "I".repeat(sequence.len())).unwrap();
    }
    drop(fastq_file);

    let clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::from([
        (7, (0..30).map(|i| format!("read-7-{}", i).into_boxed_str()).collect()),
        (2, (0..5).map(|i| format!("read-2-{}", i).into_boxed_str()).collect())
    ]);

    for max_reads_per_cluster in [10usize, 30, 0] {
        let mut options: IdentifyConsensusSequencesOptions = IdentifyConsensusSequencesOptions::DEFAULT;
        options.max_reads_per_cluster = max_reads_per_cluster;

        for num_threads in [1usize, 2, 4] {
            let consensus_set: ConsensusSequenceSet = identify_consensus_sequences(
                &clusters,
                fastq_path.to_str().unwrap(),
                &options,
                num_threads
            );

            // Vec<(cluster ID, consensus sequence, number of reads)>
            let consensus_sequences: Vec<(usize, &str, usize)> = consensus_set.sequences
                .iter()
                .map(|sequence| (
                    sequence.get_cluster_id(),
                    sequence.get_consensus_sequence(),
                    sequence.get_read_names().len()
                ))
                .collect();
            assert_eq!(
                consensus_sequences,
                vec![
                    (2, "GGGCCCAAATTTGGGCCATGCATGGTACGTAC", 5),
                    (7, "ACGTACGTACGTTTGACCATGCTAGCTAGCTA", 30)
                ],
                "max reads per cluster {}, {} threads",
                max_reads_per_cluster,
                num_threads
            );
        }
    }

    let mut options: IdentifyConsensusSequencesOptions = IdentifyConsensusSequencesOptions::DEFAULT;
    options.max_reads_per_cluster = 1;
    // Vec<consensus sequence>
    let mut first: Vec<Box<str>> = Vec::new();
    for num_threads in [1usize, 2, 4] {
        let consensus_set: ConsensusSequenceSet = identify_consensus_sequences(
            &clusters,
            fastq_path.to_str().unwrap(),
            &options,
            num_threads
        );

        let consensus_sequences: Vec<Box<str>> = consensus_set.sequences
            .iter()
            .map(|sequence| sequence.consensus_sequence.clone())
            .collect();
        for (sequence, consensus_sequence) in consensus_set.sequences.iter().zip(consensus_sequences.iter()) {
            assert!(
                records.iter().any(|(read_name, read_sequence)| {
                    sequence.get_read_names().contains(read_name.as_str()) && *read_sequence == &**consensus_sequence
                }),
                "cluster {}, {} threads",
                sequence.get_cluster_id(),
                num_threads
            );
        }
        if first.is_empty() {
            first = consensus_sequences;
        } else {
            assert_eq!(consensus_sequences, first, "{} threads", num_threads);
        }
    }
}


/// A read of a cluster that is not in the FASTQ stops the run before any alignment.
///
/// Cluster 7 holds read-7-1, which the FASTQ does not hold.
#[test]
#[should_panic(expected = "assertion `left == right` failed")]
fn identify_consensus_sequences_panics_for_reads_not_in_fastq() {
    let temp_dir = tempfile::tempdir().unwrap();
    let fastq_path = temp_dir.path().join("reads.fastq");
    let mut fastq_file = std::fs::File::create(&fastq_path).unwrap();
    writeln!(fastq_file, "@read-7-0\nACGTACGTAC\n+\nIIIIIIIIII").unwrap();
    drop(fastq_file);

    let clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::from([
        (7, HashSet::from(["read-7-0".into(), "read-7-1".into()]))
    ]);

    identify_consensus_sequences(
        &clusters,
        fastq_path.to_str().unwrap(),
        &IdentifyConsensusSequencesOptions::DEFAULT,
        1
    );
}


/// The reads of a cluster are turned to the orientation of most of them before alignment: the
/// clusters TSV and the FASTQ hold no strand, and abPOA aligns each read as it is given.
///
/// Cluster 4 holds 6 reads of a 60-base transcript, half of them reverse-complemented, as an
/// unoriented library gives them. The longest read decides the orientation of the consensus
/// sequence. With an `orientation_kmer_size` of 0 the reads are aligned as they are, and the
/// consensus sequence joins pieces of both strands.
///
///   Read       Bases    Orientation
///   read-4-0   1-60     transcript
///   read-4-1   3-58     reverse complement
///   read-4-2   1-55     transcript
///   read-4-3   1-50     reverse complement
///   read-4-4   6-60     transcript
///   read-4-5   11-60    reverse complement
#[test]
fn identify_consensus_sequences_orients_the_reads_of_a_cluster() {
    let temp_dir = tempfile::tempdir().unwrap();
    let fastq_path = temp_dir.path().join("reads.fastq");
    let records: Vec<(&str, &str)> = vec![
        ("read-4-0", "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC"),
        ("read-4-1", "AAGTTCATCGGATACCGTTAAGCTGCAACTTGGATCCTGAACTGCTAAGGTACGCC"),
        ("read-4-2", "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAA"),
        ("read-4-3", "CGGATACCGTTAAGCTGCAACTTGGATCCTGAACTGCTAAGGTACGCCAT"),
        ("read-4-4", "GTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC"),
        ("read-4-5", "GCAAGTTCATCGGATACCGTTAAGCTGCAACTTGGATCCTGAACTGCTAA")
    ];
    let mut fastq_file = std::fs::File::create(&fastq_path).unwrap();
    for (read_name, sequence) in records.iter() {
        writeln!(fastq_file, "@{}\n{}\n+\n{}", read_name, sequence, "I".repeat(sequence.len())).unwrap();
    }
    drop(fastq_file);

    let clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::from([
        (4, (0..6).map(|i| format!("read-4-{}", i).into_boxed_str()).collect())
    ]);
    let options: IdentifyConsensusSequencesOptions = IdentifyConsensusSequencesOptions::DEFAULT;

    let consensus_set: ConsensusSequenceSet = identify_consensus_sequences(
        &clusters,
        fastq_path.to_str().unwrap(),
        &options,
        1
    );

    assert_eq!(
        consensus_set.sequences[0].get_consensus_sequence(),
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC"
    );

    let mut options: IdentifyConsensusSequencesOptions = IdentifyConsensusSequencesOptions::DEFAULT;
    options.orientation_kmer_size = 0;

    let consensus_set: ConsensusSequenceSet = identify_consensus_sequences(
        &clusters,
        fastq_path.to_str().unwrap(),
        &options,
        1
    );

    assert_ne!(
        consensus_set.sequences[0].get_consensus_sequence(),
        "ATGGCGTACCTTAGCAGTTCAGGATCCAAGTTGCAGCTTAACGGTATCCGATGAACTTGC"
    );
}


/// An empty record in the FASTQ does not end the process: abPOA would, so the read is skipped.
///
/// Cluster 6 holds 3 reads, one of them empty.
#[test]
fn identify_consensus_sequences_skips_empty_reads() {
    let temp_dir = tempfile::tempdir().unwrap();
    let fastq_path = temp_dir.path().join("reads.fastq");
    let records: Vec<(&str, &str)> = vec![
        ("read-6-0", "TTGACCATGCTAGCTA"),
        ("read-6-1", ""),
        ("read-6-2", "TTGACCATGCTAGCTA")
    ];
    let mut fastq_file = std::fs::File::create(&fastq_path).unwrap();
    for (read_name, sequence) in records.iter() {
        writeln!(fastq_file, "@{}\n{}\n+\n{}", read_name, sequence, "I".repeat(sequence.len())).unwrap();
    }
    drop(fastq_file);

    let clusters: HashMap<usize, HashSet<Box<str>>> = HashMap::from([
        (6, HashSet::from(["read-6-0".into(), "read-6-1".into(), "read-6-2".into()]))
    ]);
    let options: IdentifyConsensusSequencesOptions = IdentifyConsensusSequencesOptions::DEFAULT;

    let consensus_set: ConsensusSequenceSet = identify_consensus_sequences(
        &clusters,
        fastq_path.to_str().unwrap(),
        &options,
        1
    );

    assert_eq!(consensus_set.sequences[0].get_consensus_sequence(), "TTGACCATGCTAGCTA");
    assert_eq!(consensus_set.sequences[0].get_read_names().len(), 3);
}
