use exacto_core::prelude::Gencode;
use exacto_qc::prelude::remove_unspliced_rnas;
use noodles_bam as bam;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::fs::File;
use std::path::Path;

use crate::prelude::*;

use super::*;

/// Five truncated fusion reads used to contribute an ASPA terminal clip as an insertion
/// to the ordinary WSCD1 cluster. Resolve their partner evidence before allele calling.
#[test]
fn scga_mini_rna_007_rescues_clipped_fusion_reads() {
    let root = Path::new(env!("EXACTO_TEST_DATA"));
    let bam = root.join("alignment/scga-mini-rna-007-tumor_minimap2_sorted.bam");
    let fasta = root.join("references/hg38_chr17-18.fa.gz");
    let gtf = root.join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let annotator = Gencode::new_with_defaults(gtf.to_str().unwrap(), "hg38", "v41");
    let temp = tempfile::tempdir().unwrap();
    let clusters = cluster_rna_reads(bam.to_str().unwrap(), &format!("{}.bai", bam.display()),
        fasta.to_str().unwrap(), &annotator, None, None, &ClusterRNAReadsOptions::PACBIO_HIFI,
        2, temp.path().to_str().unwrap(), temp.path().to_str().unwrap(), "rna007");
    assert_eq!(clusters.get_clusters().len(), 3);
    let fusion = clusters.get_clusters().into_iter().find(|cluster| cluster.get_variant_calls().iter()
        .any(|call| *call.get_consensus_graph_operation().get_variant_type() == VariantType::FusionGene)).unwrap();
    let names = fusion.get_read_names(clusters.get_read_names_map());
    assert_eq!(names.len(), 114, "the five resolved fragments join the original 109 fusion reads");
    for number in [85, 122, 138, 141, 154] {
        let name = format!("scga-mini-rna-007-tumor_chunk_0000/{number}/ccs");
        assert!(names.contains(name.as_str()), "missing {number}");
        let id = *clusters.get_read_names_map().get_by_left(name.as_str()).unwrap();
        assert!(fusion.get_variant_calls().iter().any(|call| {
            let op = call.get_consensus_graph_operation();
            *op.get_variant_type() == VariantType::FusionGene
                && (op.get_position_1(), op.get_position_2()) == (3_489_340, 6_087_991)
                && call.get_read_ids().contains(&id)
        }), "rescued read must support the fusion allele, not only its chain");
    }
    assert!(!clusters.get_clusters().iter().flat_map(|cluster| cluster.get_variant_calls())
        .any(|call| *call.get_consensus_graph_operation().get_variant_type() == VariantType::Insertion));
    let total: usize = clusters.get_clusters().iter().map(|cluster| cluster.get_read_ids().len()).sum();
    assert_eq!(total, 390, "rescue preserves all assigned reads");
    let unique: HashSet<_> = clusters.get_clusters().iter().flat_map(|cluster| cluster.get_read_ids()).collect();
    assert_eq!(unique.len(), total, "no read is assigned twice");
}

/// At a junction read floor of 1 every read lands in one cluster at most, and every read but two
/// lands in one. `chunk_0000/91` and `chunk_0000/188` each hold a deletion against a junction,
/// which their models spell as a junction one base off that no other read spells; at a maximum
/// FPR of 1 a novel splice junction still needs 2 reads.
#[test]
fn scga_mini_rna_001_cluster_rna_reads_returns_clusters() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bai_file: String = format!("{}.bai", bam_file);
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");
    let temp_dir = tempfile::tempdir().unwrap();

    let mut options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;
    options.calling.min_mapping_quality = 20;
    options.calling.min_terminal_soft_clip_ins_len = 0;
    options.junction.min_reads = 1;
    options.filtering.max_fpr = 1.0;

    let cluster_set: RNAReadClusterSet = cluster_rna_reads(
        bam_file,
        bai_file.as_str(),
        fasta_file,
        &gene_annotator,
        None,
        None,
        &options,
        2,
        temp_dir.path().to_str().unwrap(),
        temp_dir.path().to_str().unwrap(),
        "test"
    );

    let mut seen: HashMap<usize, usize> = HashMap::new();
    for cluster in cluster_set.get_clusters() {
        for read_id in cluster.get_read_ids().iter() {
            *seen.entry(*read_id).or_insert(0) += 1;
        }
    }
    let mut unclustered: Vec<&str> = cluster_set
        .get_read_names_map()
        .iter()
        .filter(|(_, read_id)| !seen.contains_key(read_id))
        .map(|(read_name, _)| &**read_name)
        .collect();
    unclustered.sort_unstable();

    assert!(seen.values().all(|count| *count == 1), "no read lands in two clusters");
    assert_eq!(
        unclustered,
        vec!["scga-mini-rna-001-tumor_chunk_0000/188/ccs", "scga-mini-rna-001-tumor_chunk_0000/91/ccs"]
    );
}


/// A circular RNA clusters apart from its linear host.
///
/// It cannot do so on splice junctions alone. A circular read's *linear* junctions are a
/// **contiguous** run of the host transcript's chain — entries 7, 8 and 9 of TP53's ten — so
/// Step 8's degradation fold would commit every one of them to the linear backbone as a 5'/3'
/// truncation fragment, and the locus would come back as a single cluster mixing both transcripts.
///
/// Phasing cannot recover it either, so this is not a threshold that more depth would clear. The
/// CIR site spans 3,677 bp across introns and a non-carrier's exon cannot contain all of it, so
/// non-carriers genotype `NotCovered` rather than `Reference`. That reaches the MEC matrix as
/// `None`, which the distance and the consensus both skip — no zeros, cost 0 under every
/// partition, invisible to BIC.
///
/// Chaining on the back-splice is what separates them: the circular chain then holds a key the
/// linear chain does not, so neither contains the other and both stand as maximal chains.
///
/// scga-mini-rna-011 simulates 200 reads of the circular RNA (`chunk_0000/1` to `/200`, the first
/// row of its `.transcript` file) and 200 of the linear host (`/201` to `/400`).
#[test]
fn scga_mini_rna_011_cluster_rna_reads_separates_circular_rna() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-011-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bai_file: String = format!("{}.bai", bam_file);
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");
    let temp_dir = tempfile::tempdir().unwrap();

    let cluster_set: RNAReadClusterSet = cluster_rna_reads(
        bam_file,
        bai_file.as_str(),
        fasta_file,
        &gene_annotator,
        None,
        None,
        &ClusterRNAReadsOptions::DEFAULT,
        2,
        temp_dir.path().to_str().unwrap(),
        temp_dir.path().to_str().unwrap(),
        "test"
    );
    let is_circular_read = |read_name: &str| -> bool {
        let read_number: usize = read_name
            .strip_prefix("scga-mini-rna-011-tumor_chunk_0000/")
            .and_then(|rest| rest.strip_suffix("/ccs"))
            .map_or(0, |number| number.parse().unwrap());
        (1..=200).contains(&read_number)
    };

    // The cluster carrying the back-splice holds every read carrying it, and nothing else, and
    // every one of them is a read of the circular RNA. Asserting both directions matters: the
    // failure this guards against is a *merge*, which leaves the circular reads clustered.
    let circular_cluster: &RNAReadCluster = cluster_set
        .get_clusters()
        .iter()
        .find(|cluster| cluster.get_variant_calls().iter().any(|call| {
            *call.get_consensus_graph_operation().get_variant_type() == VariantType::CircularRNA
        }))
        .expect("No cluster carries the back-splice.");
    let carriers: HashSet<usize> = circular_cluster.get_variant_calls()[0].get_read_ids().into_iter().collect();
    let read_ids: HashSet<usize> = circular_cluster.get_read_ids().iter().copied().collect();
    assert_eq!(read_ids, carriers, "The circular cluster holds reads beyond those carrying the back-splice.");
    assert_eq!(read_ids.len(), 151);
    assert!(
        circular_cluster.get_read_names(&cluster_set.get_read_names_map()).iter().all(|read_name| is_circular_read(read_name)),
        "A read of the linear host is with the circular reads."
    );

    // The host transcript keeps its own cluster, with its first intron.
    assert_eq!(cluster_set.get_clusters().len(), 2);
    let linear_cluster: &RNAReadCluster = cluster_set
        .get_clusters()
        .iter()
        .find(|cluster| cluster.get_id() != circular_cluster.get_id())
        .unwrap();
    assert!(linear_cluster.get_variant_calls().is_empty());
    assert!(linear_cluster
        .get_splice_junctions()
        .iter()
        .any(|junction| (junction.position_1, junction.position_2) == (7_687_376, 7_676_623)));
}

/// The back-splice chains but is never published as a splice junction.
///
/// `RNAReadCluster` feeds the splice-junctions TSV, where every row is read as a linear intron. A
/// cycle is not one — on a reverse-strand transcript it is spelled `pos1 < pos2`, the reverse of
/// the convention every other row follows — so a consumer would read it as a strand error. It is
/// already reported as a `CIR` variant call, which is where it belongs.
#[test]
fn scga_mini_rna_011_cluster_rna_reads_omits_cycles_from_splice_junctions() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-011-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bai_file: String = format!("{}.bai", bam_file);
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");
    let temp_dir = tempfile::tempdir().unwrap();

    let cluster_set: RNAReadClusterSet = cluster_rna_reads(
        bam_file,
        bai_file.as_str(),
        fasta_file,
        &gene_annotator,
        None,
        None,
        &ClusterRNAReadsOptions::DEFAULT,
        2,
        temp_dir.path().to_str().unwrap(),
        temp_dir.path().to_str().unwrap(),
        "test"
    );

    // The back-splice of the circular RNA.
    const CYCLE_POSITION_1: u32 = 7_670_611;
    const CYCLE_POSITION_2: u32 = 7_674_288;

    for cluster in cluster_set.get_clusters() {
        for junction in cluster.get_splice_junctions().iter() {
            assert!(
                !(junction.position_1 == CYCLE_POSITION_1 && junction.position_2 == CYCLE_POSITION_2),
                "Cluster {} published the back-splice as a splice junction.",
                cluster.get_id()
            );
        }
    }

    // The cluster that chains on it still exists and calls it, so the assertion above is not
    // passing merely because the circular reads went missing.
    assert!(
        cluster_set.get_clusters().iter().flat_map(|cluster| cluster.get_variant_calls().iter()).any(|call| {
            let graph_operation: &GraphOperation = call.get_consensus_graph_operation();
            *graph_operation.get_variant_type() == VariantType::CircularRNA
                && graph_operation.get_position_1() == CYCLE_POSITION_1
                && graph_operation.get_position_2() == CYCLE_POSITION_2
        }),
        "The back-splice is not called; the check above proved nothing."
    );
}


/// Reads that disagree about where an insertion/junction boundary sits still cluster together.
///
/// scga-mini-rna-015's fusion joins WSCD1 exon 7 to ACAP1 exon 2. On reads too short at the
/// ACAP1 end to place the join, minimap2 chains 17 bases of ACAP1 to a copy of them at 6,165,489
/// and spells the end of the read as a large insertion and a junction into the 17 bases. It places
/// that junction in several spots, because the insertion's leading 136 bases are exactly the
/// reference bases chr17:6,110,800-6,110,935 and rolling them into the insertion costs nothing.
///
/// `normalise_junction_boundary` slides each read to the rightmost sequence-preserving placement,
/// which collapses the fork. The slide is an exact match and stops at the first sequencing error
/// inside the microhomology, so a read whose slide stops early would be left in a cluster of its
/// own.
///
/// The placement they converge on, chr17:6,110,936, is the annotated GENCODE donor and reads
/// `GTAAGT`. The stretch `identify_unplaced_tail` takes out of a read together with this junction
/// is left in: the junction is the one this test pools.
#[test]
fn scga_mini_rna_015_cluster_rna_reads_pools_insertion_boundary_spellings() {
    use noodles_sam::alignment::record::cigar::op::Kind;

    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-015-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bai_file: String = format!("{}.bai", bam_file);
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");
    let temp_dir = tempfile::tempdir().unwrap();

    // The first intron base of the junction into the 17 bases, as the aligner spelled it on each
    // primary record. HashMap<read name, first intron base>
    let mut spellings: HashMap<Box<str>, usize> = HashMap::new();
    let mut reader = bam::io::reader::Builder::default().build_from_path(bam_file).unwrap();
    let header: noodles_sam::Header = reader.read_header().unwrap();
    for result in reader.record_bufs(&header) {
        let record = result.unwrap();
        if record.flags().is_unmapped() || record.flags().is_secondary() || record.flags().is_supplementary() {
            continue;
        }
        let mut position: usize = record.alignment_start().unwrap().get();
        for op in record.cigar().as_ref().iter() {
            if op.kind() == Kind::Skip && position + op.len() - 1 == 6_165_488 {
                spellings.insert(record.name().unwrap().to_string().into(), position);
            }
            if matches!(op.kind(), Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch) {
                position += op.len();
            }
        }
    }
    let mut options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;
    options.calling.unplaced_tail_removal = false;
    let cluster_set: RNAReadClusterSet = cluster_rna_reads(
        bam_file,
        bai_file.as_str(),
        fasta_file,
        &gene_annotator,
        None,
        None,
        &options,
        2,
        temp_dir.path().to_str().unwrap(),
        temp_dir.path().to_str().unwrap(),
        "test"
    );

    // One cluster holds the junction, at the annotated donor, and every read spelling it that
    // reaches a cluster, from more than one placement.
    let clusters: Vec<&RNAReadCluster> = cluster_set
        .get_clusters()
        .into_iter()
        .filter(|cluster| cluster.get_splice_junctions().iter().any(|junction| junction.position_2 == 6_165_488))
        .collect();
    assert_eq!(clusters.len(), 1, "The spellings of one boundary are still separate clusters.");
    assert!(
        clusters[0].get_splice_junctions().iter().any(|junction| (junction.position_1, junction.position_2) == (6_110_936, 6_165_488)),
        "The pooled cluster's junction is not at the annotated donor 6,110,936: {:?}",
        clusters[0].get_splice_junctions().iter().map(|junction| (junction.position_1, junction.position_2)).collect::<Vec<_>>()
    );
    let read_names: HashSet<Box<str>> = clusters[0].get_read_names(&cluster_set.get_read_names_map());
    let starts: HashSet<usize> = spellings
        .iter()
        .filter(|(read_name, _)| read_names.contains(*read_name))
        .map(|(_, start)| *start)
        .collect();
    assert!(starts.len() > 1, "The pooled reads were spelled in one place only: {:?}", starts);
    for cluster in cluster_set.get_clusters() {
        if cluster.get_id() != clusters[0].get_id() {
            let other_read_names: HashSet<Box<str>> = cluster.get_read_names(&cluster_set.get_read_names_map());
            assert!(
                spellings.keys().all(|read_name| !other_read_names.contains(read_name)),
                "A read spelling the junction is in cluster {}.",
                cluster.get_id()
            );
        }
    }
}


/// A phantom fork of a mixed-strand (inversion) cluster is removed by the junction read floor
/// while its parent survives it.
///
/// scga-mini-rna-013's inversion reads hold a forward arm and an SA-linked reverse-strand arm. 81
/// of them spell the junction into the inverted segment as 7,673,700-7,669,660; three spell it
/// 7,673,700-7,669,670, the alignment-respelling fork observed on rna-006. Same-length chains
/// differing in one junction cannot fold, so the phantom chain forms a cluster of its own whenever
/// the floor lets it.
///
/// The `min_reads = 3` run guards against a vacuous pass: the phantom chain must actually reach
/// a cluster when the floor admits it, or the assertion at `min_reads = 4` would hold with
/// split-read modelling itself broken.
///
/// The respelled junction is a novel splice junction 3 reads splice, and at the default
/// maximum FPR its reads are removed before the floor sees them. The floor is the subject
/// here, so the maximum FPR is 1, at which the minimum read support of every novel splice
/// junction is 2.
#[test]
fn scga_mini_rna_013_cluster_rna_reads_drops_phantom_fork_below_junction_read_floor() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-013-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bai_file: String = format!("{}.bai", bam_file);
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");

    let mut options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;
    options.filtering.max_fpr = 1.0;

    // The clusters holding a junction, as (reads, number of breakend calls).
    let clusters_with = |cluster_set: &RNAReadClusterSet, junction: (u32, u32)| -> Vec<(usize, usize)> {
        cluster_set
            .get_clusters()
            .iter()
            .filter(|cluster| cluster.get_splice_junctions().iter().any(|splice_junction| (splice_junction.position_1, splice_junction.position_2) == junction))
            .map(|cluster| (
                cluster.get_read_ids().len(),
                cluster
                    .get_variant_calls()
                    .iter()
                    .filter(|call| *call.get_consensus_graph_operation().get_variant_type() == VariantType::Breakpoint)
                    .count()
            ))
            .collect()
    };

    let mut found: Vec<(Vec<(usize, usize)>, Vec<(usize, usize)>)> = Vec::new();
    for min_reads in [3, 4] {
        let temp_dir = tempfile::tempdir().unwrap();
        let mut options: ClusterRNAReadsOptions = options.clone();
        options.junction.min_reads = min_reads;
        let cluster_set: RNAReadClusterSet = cluster_rna_reads(
            bam_file,
            bai_file.as_str(),
            fasta_file,
            &gene_annotator,
            None,
            None,
            &options,
            2,
            temp_dir.path().to_str().unwrap(),
            temp_dir.path().to_str().unwrap(),
            "test"
        );
        found.push((clusters_with(&cluster_set, (7_673_700, 7_669_670)), clusters_with(&cluster_set, (7_673_700, 7_669_660))));
    }

    // Guard: with the floor at 3, the phantom chain reaches a cluster of its own — the fork
    // exists and cannot fold. The point: with the floor at 4 the phantom is gone and the parent,
    // with both breakend calls, is untouched.
    assert_eq!(found[0], (vec![(3, 2)], vec![(81, 2)]), "floor 3: (phantom, parent)");
    assert_eq!(found[1], (vec![], vec![(81, 2)]), "floor 4: (phantom, parent)");
}


/// Clustering must not depend on `HashMap`/`HashSet` iteration order.
///
/// Rust seeds `RandomState` from the OS once per process and then bumps the key for
/// every `HashMap` created, so two `cluster_rna_reads` calls in the SAME process already
/// hash differently — which is exactly the perturbation this test needs, and why it does
/// not have to shell out to a second process to be meaningful.
///
/// Asserted on the full (id, reads, junctions, calls) tuple rather than on a cluster count:
/// a count is invariant under exactly the permutation this is meant to catch.
#[test]
fn cluster_rna_reads_is_deterministic_across_runs() {
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");

    // At the default options each of these samples already leaves hundreds of failed calls,
    // several cells and a phasing per cluster in play, the permutation surface this test needs.
    let options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;

    // A digest that is sensitive to cluster IDs, membership, the id->cluster pairing, the
    // junction chain, the published variant call spellings, the reference transcript picks,
    // and the failed calls. A run-dependent variant spelling would change correct-rna-reads'
    // edits while leaving ids, membership and junctions identical, so those fields are
    // load-bearing here, not decoration.
    type ClusterDigest = (
        usize,                        // cluster id
        Vec<Box<str>>,                // sorted read names
        Vec<SpliceJunction>,          // junction chain
        Vec<Box<str>>,                // sorted variant call spellings
        Vec<(Box<str>, Box<str>)>     // sorted reference (gene, transcript) ids
    );
    let digest = |cluster_set: &RNAReadClusterSet| -> (Vec<ClusterDigest>, Vec<(Box<str>, &'static str)>) {
        let mut rows: Vec<ClusterDigest> = cluster_set
            .get_clusters()
            .iter()
            .map(|cluster| {
                let mut read_names: Vec<Box<str>> = cluster.get_read_names(&cluster_set.get_read_names_map()).into_iter().collect();
                read_names.sort();
                let mut variants: Vec<Box<str>> = cluster
                    .get_variant_calls()
                    .iter()
                    .map(|variant_call| variant_call.get_consensus_graph_operation().as_boxed_str())
                    .collect();
                variants.sort();
                let mut reference_ids: Vec<(Box<str>, Box<str>)> =
                    cluster.get_reference_gene_transcript_ids().iter().cloned().collect();
                reference_ids.sort();
                (
                    cluster.get_id(),
                    read_names,
                    cluster.get_splice_junctions().clone(),
                    variants,
                    reference_ids
                )
            })
            .collect();
        rows.sort();
        let mut failed: Vec<(Box<str>, &'static str)> = cluster_set
            .get_failed_variant_calls()
            .iter()
            .map(|failed_variant_call| (
                failed_variant_call.graph_operation.as_boxed_str(),
                failed_variant_call.failure.as_str()
            ))
            .collect();
        failed.sort();
        (rows, failed)
    };

    // rna-001 holds spliced, unspliced and partially spliced reads of one gene; rna-011 a
    // circular RNA beside its host; rna-015 a fusion of three genes, whose truncated reads are
    // fragments of more than one backbone, the Step 8 ordering that has to be total.
    for fixture in ["scga-mini-rna-001", "scga-mini-rna-011", "scga-mini-rna-015"] {
        let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join(format!("alignment/{}-tumor_minimap2_sorted.bam", fixture));
        let bam_full_path = fs::canonicalize(bam_path).unwrap();
        let bam_file: &str = bam_full_path.to_str().unwrap();
        let bai_file: String = format!("{}.bai", bam_file);

        let temp_dir_first = tempfile::tempdir().unwrap();
        let first: RNAReadClusterSet = cluster_rna_reads(
            bam_file,
            bai_file.as_str(),
            fasta_file,
            &gene_annotator,
            None,
            None,
            &options,
            2,
            temp_dir_first.path().to_str().unwrap(),
            temp_dir_first.path().to_str().unwrap(),
            "test"
        );
        let temp_dir_second = tempfile::tempdir().unwrap();
        let second: RNAReadClusterSet = cluster_rna_reads(
            bam_file,
            bai_file.as_str(),
            fasta_file,
            &gene_annotator,
            None,
            None,
            &options,
            2,
            temp_dir_second.path().to_str().unwrap(),
            temp_dir_second.path().to_str().unwrap(),
            "test"
        );

        let (first_digest, second_digest) = (digest(&first), digest(&second));
        assert!(
            !first_digest.0.is_empty(),
            "fixture {} produced no clusters; the assertion below would be vacuous",
            fixture
        );
        assert_eq!(
            first_digest, second_digest,
            "two clustering runs of {} disagreed, so HashMap iteration order is still reaching the output",
            fixture
        );
    }
}


/// A variant call the DNA list vouches for is kept regardless of the read-support gates.
///
/// scga-mini-rna-013's inversion is called as two breakends, each on the 81 reads of its
/// cluster. At a read floor of 82 both fail it, land among the failed calls with `TooFewReads`,
/// and the cluster is published without them. Listing one breakend as a DNA variant is the
/// difference: that call is kept on the cluster and nothing is left in the failed calls for it,
/// while the breakend the list does not name still fails.
///
/// The DNA row is spelled from the failed call itself rather than from constants, so the test
/// pins the gate and not a coordinate that a change in breakend placement would move.
#[test]
fn scga_mini_rna_013_cluster_rna_reads_keeps_dna_listed_breakend_call() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-013-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bai_file: String = format!("{}.bai", bam_file);
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");

    let mut options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;
    options.filtering.min_reads = 82;

    let temp_dir_unlisted = tempfile::tempdir().unwrap();
    let unlisted: RNAReadClusterSet = cluster_rna_reads(
        bam_file,
        bai_file.as_str(),
        fasta_file,
        &gene_annotator,
        None,
        None,
        &options,
        2,
        temp_dir_unlisted.path().to_str().unwrap(),
        temp_dir_unlisted.path().to_str().unwrap(),
        "test"
    );
    let breakends = |cluster_set: &RNAReadClusterSet| -> (Vec<Box<str>>, Vec<Box<str>>) {
        let mut kept: Vec<Box<str>> = cluster_set
            .get_clusters()
            .iter()
            .flat_map(|cluster| cluster.get_variant_calls().iter())
            .filter(|call| *call.get_consensus_graph_operation().get_variant_type() == VariantType::Breakpoint)
            .map(|call| call.get_consensus_graph_operation().as_boxed_str())
            .collect();
        kept.sort();
        let mut failed: Vec<Box<str>> = cluster_set
            .get_failed_variant_calls()
            .iter()
            .filter(|failed| *failed.graph_operation.get_variant_type() == VariantType::Breakpoint)
            .map(|failed| failed.graph_operation.as_boxed_str())
            .collect();
        failed.sort();
        (kept, failed)
    };

    // Guard: both breakends fail the read floor rather than being absent for some other reason.
    let failed_breakends: Vec<&FailedRNAVariantCall> = unlisted
        .get_failed_variant_calls()
        .iter()
        .filter(|failed| *failed.graph_operation.get_variant_type() == VariantType::Breakpoint)
        .collect();
    assert_eq!(failed_breakends.len(), 2, "Without a DNA list both breakends must be among the failed calls.");
    assert!(failed_breakends.iter().all(|failed| failed.num_reads == 81 && matches!(failed.failure, VariantCallFailure::TooFewReads { .. })));
    assert!(breakends(&unlisted).0.is_empty());
    let breakend_op: &GraphOperation = &failed_breakends[0].graph_operation;
    let other_op: &GraphOperation = &failed_breakends[1].graph_operation;

    // The first breakend as one DNA row, spelled the way the consensus call reports it.
    let dna_variant_records: Vec<DNAVariantRecord> = vec![DNAVariantRecord {
        origin: "somatic".into(),
        variant_id: 1,
        chromosome_1: unlisted.get_chromosome_names_map().get_by_right(&breakend_op.get_chromosome_1()).unwrap().clone(),
        position_1: breakend_op.get_position_1(),
        strand_1: breakend_op.get_strand_1().as_str().into(),
        operation_1: breakend_op.get_operation_type_1().as_str().into(),
        chromosome_2: unlisted.get_chromosome_names_map().get_by_right(&breakend_op.get_chromosome_2()).unwrap().clone(),
        position_2: breakend_op.get_position_2(),
        strand_2: breakend_op.get_strand_2().as_str().into(),
        operation_2: breakend_op.get_operation_type_2().as_str().into(),
        sequence: breakend_op.get_sequence().into(),
        variant_size: Some(0),
        variant_type: "BND".into(),
        consensus_read_names: "".into(),
        num_consensus_read_names: 0,
        read_names: "".into(),
        num_read_names: 0
    }];

    let temp_dir_listed = tempfile::tempdir().unwrap();
    let listed: RNAReadClusterSet = cluster_rna_reads(
        bam_file,
        bai_file.as_str(),
        fasta_file,
        &gene_annotator,
        Some(&dna_variant_records),
        None,
        &options,
        2,
        temp_dir_listed.path().to_str().unwrap(),
        temp_dir_listed.path().to_str().unwrap(),
        "test"
    );

    // The point: the listed breakend is kept and no longer fails; the other still does.
    assert_eq!(
        breakends(&listed),
        (vec![breakend_op.as_boxed_str()], vec![other_op.as_boxed_str()])
    );
}


/// The reads of one fusion transcript are one cluster, however the aligner spelled them.
///
/// scga-mini-rna-015 holds ASPA, WSCD1 and ACAP1 (chr17, forward strand, more than 1 Mb apart),
/// their three ordinary transcripts, and one fusion transcript joining ASPA exon 4 to WSCD1 exon 3
/// and WSCD1 exon 7 to ACAP1 exon 2. pbsim3 numbers the reads of each row of its `.transcript`
/// file in turn: `chunk_0000/1` to `/200` are the fusion, `/201` to `/400` WSCD1, `/401` to `/600`
/// ASPA and `/601` to `/800` ACAP1.
///
/// On reads too short at the ACAP1 end to place the join, minimap2 spells the end of the read as
/// a large insertion, a junction into 17 bases at 6,165,489 and a clip. That junction would key a
/// chain of its own; the stretch is taken out of the read and the read joins the fusion.
///
/// A fusion read holding junctions of only one gene carries no evidence of the fusion and goes
/// with that gene's transcript, so the fusion cluster is all fusion reads but not every fusion
/// read is in it.
#[test]
fn scga_mini_rna_015_cluster_rna_reads_keeps_fusion_reads_in_one_cluster() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-015-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bai_file: String = format!("{}.bai", bam_file);
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");
    let temp_dir = tempfile::tempdir().unwrap();
    let options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;

    let cluster_set: RNAReadClusterSet = cluster_rna_reads(
        bam_file,
        bai_file.as_str(),
        fasta_file,
        &gene_annotator,
        None,
        None,
        &options,
        2,
        temp_dir.path().to_str().unwrap(),
        temp_dir.path().to_str().unwrap(),
        "test"
    );
    // The transcript a read was simulated from: 1 fusion, 2 WSCD1, 3 ASPA, 4 ACAP1, 0 any other.
    let transcript_of = |read_name: &str| -> usize {
        let read_number: usize = read_name
            .strip_prefix("scga-mini-rna-015-tumor_chunk_0000/")
            .and_then(|rest| rest.strip_suffix("/ccs"))
            .map_or(0, |number| number.parse().unwrap());
        if (1..=800).contains(&read_number) { (read_number + 199) / 200 } else { 0 }
    };
    let transcripts = |cluster: &RNAReadCluster| -> HashSet<usize> {
        cluster.get_read_names(&cluster_set.get_read_names_map()).iter().map(|read_name| transcript_of(read_name)).collect()
    };

    // The fusion cluster holds fusion reads only, with the chain of the whole transcript and
    // without the junction into the 17 bases.
    let fusion_cluster: &RNAReadCluster = cluster_set
        .get_clusters()
        .iter()
        .find(|cluster| cluster.get_variant_calls().iter().any(|call| *call.get_consensus_graph_operation().get_variant_type() == VariantType::FusionGene))
        .expect("No cluster calls the fusion.");
    assert_eq!(transcripts(fusion_cluster), HashSet::from([1]));
    assert_eq!(fusion_cluster.get_read_ids().len(), 181);
    let junctions: Vec<(u32, u32)> = fusion_cluster
        .get_splice_junctions()
        .iter()
        .map(|junction| (junction.position_1, junction.position_2))
        .collect();
    assert_eq!(
        junctions,
        vec![
            (3_476_396, 3_481_602),
            (3_481_799, 3_483_498),
            (3_483_593, 3_489_234),
            (6_088_105, 6_090_320),
            (6_090_506, 6_095_101),
            (6_095_224, 6_109_606),
            (6_109_767, 6_110_770),
            (7_342_068, 7_342_274),
            (7_342_329, 7_342_415),
            (7_342_475, 7_343_378)
        ]
    );
    assert!(!cluster_set
        .get_clusters()
        .iter()
        .any(|cluster| cluster.get_splice_junctions().iter().any(|junction| junction.position_2 == 6_165_488)));

    // Its calls are the two joins of the fusion. The first is carried by the reads holding ASPA,
    // four of them only as a clip; the second by the reads the aligner placed on ACAP1. The
    // insertion, the clip and the join to 6,165,489 of the others are no calls.
    let mut calls: Vec<(VariantType, u32, u32, usize)> = fusion_cluster
        .get_variant_calls()
        .iter()
        .map(|variant_call| {
            let graph_operation: &GraphOperation = variant_call.get_consensus_graph_operation();
            let read_ids: HashSet<usize> = variant_call.get_read_ids().into_iter().collect();
            (
                graph_operation.get_variant_type().clone(),
                graph_operation.get_position_1(),
                graph_operation.get_position_2(),
                read_ids.len()
            )
        })
        .collect();
    calls.sort_by_key(|(_, position_1, position_2, _)| (*position_1, *position_2));
    assert_eq!(
        calls,
        vec![
            (VariantType::FusionGene, 3_489_340, 6_087_991, 172),
            (VariantType::FusionGene, 6_110_933, 7_341_948, 68)
        ]
    );

    // The three ordinary transcripts are the only other clusters, one each, holding every read
    // of their own transcript.
    assert_eq!(cluster_set.get_clusters().len(), 4);
    for transcript in [2, 3, 4] {
        let holding: Vec<&RNAReadCluster> = cluster_set
            .get_clusters()
            .into_iter()
            .filter(|cluster| transcripts(cluster).contains(&transcript))
            .collect();
        assert_eq!(holding.len(), 1, "transcript {}", transcript);
    }
}


/// The stretch a read holds unplaced stays in the read when the removal is off, or when its
/// insertion is under the minimum length.
///
/// 120 reads of scga-mini-rna-015's fusion end in a large insertion, a junction into 17 bases at
/// 6,165,489 and a clip, such as `29= 467I 54689N 17= 68S`: 44 on their primary record, 76 on a
/// supplementary one. Their insertions run from 419 to 521 bases.
///
///   Removal   Minimum insertion length   Clusters   Junction into the 17 bases
///   on        400                        4          in no cluster; the reads join the fusion
///   on        500                        5          keys a cluster of its own
///   off       100                        5          keys a cluster of its own
///
/// Left in, the junction into the 17 bases keys a chain of its own, and the insertion, the clip
/// and the join to 6,165,489 are its calls.
#[test]
fn scga_mini_rna_015_cluster_rna_reads_keeps_unplaced_tails_under_min_ins_len() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-015-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bai_file: String = format!("{}.bai", bam_file);
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");

    // (unplaced tail removal, unplaced tail min ins len, number of clusters, stretch taken out)
    let cases: Vec<(bool, u32, usize, bool)> = vec![
        (true, 400, 4, true),
        (true, 500, 5, false),
        (false, 100, 5, false)
    ];
    for (unplaced_tail_removal, unplaced_tail_min_ins_len, num_clusters, expected) in cases {
        let temp_dir = tempfile::tempdir().unwrap();
        let mut options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;
        options.calling.unplaced_tail_removal = unplaced_tail_removal;
        options.calling.unplaced_tail_min_ins_len = unplaced_tail_min_ins_len;

        let cluster_set: RNAReadClusterSet = cluster_rna_reads(
            bam_file,
            bai_file.as_str(),
            fasta_file,
            &gene_annotator,
            None,
            None,
            &options,
            2,
            temp_dir.path().to_str().unwrap(),
            temp_dir.path().to_str().unwrap(),
            "test"
        );
        assert_eq!(
            cluster_set.get_clusters().len(),
            num_clusters,
            "removal {}, min ins len {}",
            unplaced_tail_removal,
            unplaced_tail_min_ins_len
        );

        // The cluster holding the junction into the 17 bases, when the stretch is left in, and
        // the calls spelling the stretch with it.
        let decoy_cluster: Option<&RNAReadCluster> = cluster_set
            .get_clusters()
            .into_iter()
            .find(|cluster| cluster.get_splice_junctions().iter().any(|junction| junction.position_2 == 6_165_488));
        assert_eq!(
            decoy_cluster.is_none(),
            expected,
            "removal {}, min ins len {}",
            unplaced_tail_removal,
            unplaced_tail_min_ins_len
        );
        if let Some(decoy_cluster) = decoy_cluster {
            let mut calls: Vec<(VariantType, u32, u32)> = decoy_cluster
                .get_variant_calls()
                .iter()
                .map(|variant_call| {
                    let graph_operation: &GraphOperation = variant_call.get_consensus_graph_operation();
                    (
                        graph_operation.get_variant_type().clone(),
                        graph_operation.get_position_1(),
                        graph_operation.get_position_2()
                    )
                })
                .collect();
            calls.sort_by_key(|(_, position_1, position_2)| (*position_1, *position_2));
            assert_eq!(
                calls,
                vec![
                    (VariantType::FusionGene, 3_489_340, 6_087_991),
                    (VariantType::Insertion, 6_110_799, 6_110_800),
                    (VariantType::FusionGene, 6_110_799, 6_165_489),
                    (VariantType::Insertion, 6_165_505, 6_165_506)
                ]
            );
        }
    }
}


/// A terminal soft clip that spells the reference across an intron makes no call, unless the
/// removal is off or the boundary of the intron is further than the maximum distance.
///
/// Nine reads of scga-mini-rna-011's linear TP53 (chr17, reverse strand) start 9 bases into the
/// exon above the intron 7,674,972-7,675,052, too few for the aligner to splice. The intron
/// starts with the base the exon does, so the alignment runs one base into it, to 7,674,972, and
/// the 9 bases are clipped. The clips spell the exon with no edit.
///
///   Removal   Boundary within   Bases per edit   Calls of the linear cluster
///   on        3                 8                none
///   off       3                 8                INS of 9 bases after 7,674,972, 9 reads
///   on        0                 8                INS of 9 bases after 7,674,972, 9 reads
///   on        3                 0                none
#[test]
fn scga_mini_rna_011_cluster_rna_reads_removes_soft_clips_across_intron() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-011-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bai_file: String = format!("{}.bai", bam_file);
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");

    // (soft clip removal, max boundary distance, bases per edit, calls of the linear cluster)
    let cases: Vec<(bool, u32, u32, Vec<(VariantType, u32, u32, usize, u32)>)> = vec![
        (true, 3, 8, Vec::new()),
        (false, 3, 8, vec![(VariantType::Insertion, 7_674_972, 7_674_973, 9, 9)]),
        (true, 0, 8, vec![(VariantType::Insertion, 7_674_972, 7_674_973, 9, 9)]),
        (true, 3, 0, Vec::new())
    ];
    for (soft_clip_removal, soft_clip_max_boundary_distance, soft_clip_bases_per_edit, expected) in cases {
        let temp_dir = tempfile::tempdir().unwrap();
        let mut options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;
        options.calling.soft_clip_removal = soft_clip_removal;
        options.calling.soft_clip_max_boundary_distance = soft_clip_max_boundary_distance;
        options.calling.soft_clip_bases_per_edit = soft_clip_bases_per_edit;

        let cluster_set: RNAReadClusterSet = cluster_rna_reads(
            bam_file,
            bai_file.as_str(),
            fasta_file,
            &gene_annotator,
            None,
            None,
            &options,
            2,
            temp_dir.path().to_str().unwrap(),
            temp_dir.path().to_str().unwrap(),
            "test"
        );

        // The linear cluster is the one holding TP53's first intron; the circular RNA's does not.
        let linear_cluster: &RNAReadCluster = cluster_set
            .get_clusters()
            .iter()
            .find(|cluster| cluster.get_splice_junctions().iter().any(|junction| (junction.position_1, junction.position_2) == (7_687_376, 7_676_623)))
            .unwrap();
        let calls: Vec<(VariantType, u32, u32, usize, u32)> = linear_cluster
            .get_variant_calls()
            .iter()
            .map(|variant_call| {
                let graph_operation: &GraphOperation = variant_call.get_consensus_graph_operation();
                (
                    graph_operation.get_variant_type().clone(),
                    graph_operation.get_position_1(),
                    graph_operation.get_position_2(),
                    graph_operation.get_sequence().len(),
                    variant_call.get_num_reads()
                )
            })
            .collect();
        assert_eq!(
            calls,
            expected,
            "removal {}, boundary within {}, {} bases per edit",
            soft_clip_removal,
            soft_clip_max_boundary_distance,
            soft_clip_bases_per_edit
        );
    }
}


/// The read support index is the same on any number of threads, and the same built from a
/// worker of another thread pool. There the builder calculates its cells on the worker and
/// builds no pool: a pool built on every call from the workers of another pool exhausts the
/// threads of the process.
///
/// Row 0 is a variant outside a repeat, at a sequencing error of 0.01 and a max FPR of 1e-6:
///
///   depth        0    1    2   30   100   1,000   10,000
///   minimum      2    2    3    9    19     138    1,319
#[test]
fn rna_variant_read_support_index_returns_matches_on_any_thread_pool() {
    let depths: HashSet<u32> = [0u32, 1, 2, 30, 100, 1_000, 10_000].into_iter().collect();

    let on_one_thread: RNAVariantReadSupportIndex = RNAVariantReadSupportIndex::new(
        &depths,
        30,
        0.01,
        0.01,
        1e-6,
        1
    );
    let on_four_threads: RNAVariantReadSupportIndex = RNAVariantReadSupportIndex::new(
        &depths,
        30,
        0.01,
        0.01,
        1e-6,
        4
    );
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap();
    let on_worker: Vec<RNAVariantReadSupportIndex> = thread_pool.install(|| {
        (0..16)
            .into_par_iter()
            .map(|_| RNAVariantReadSupportIndex::new(
                &depths,
                30,
                0.01,
                0.01,
                1e-6,
                4
            ))
            .collect()
    });

    for (depth, min_read_support) in [(0u32, 2u32), (1, 2), (2, 3), (30, 9), (100, 19), (1_000, 138), (10_000, 1_319)] {
        assert_eq!(on_one_thread.get_min_read_support((0, depth)), min_read_support, "depth {}", depth);
    }
    // Rows 0 and 2 to 30, at every depth: a missing cell panics.
    for repeat_length in std::iter::once(0).chain(2..=30) {
        for depth in depths.iter() {
            on_one_thread.get_min_read_support((repeat_length, *depth));
        }
    }
    assert_eq!(on_four_threads, on_one_thread);
    for index in on_worker.iter() {
        assert_eq!(index, &on_one_thread);
    }
}


/// The tables are written to the output directory, under the output prefix, and hold the
/// clusters the function returns.
#[test]
fn scga_mini_rna_001_cluster_rna_reads_writes_tables_to_output_dir() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bai_file: String = format!("{}.bai", bam_file);
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");
    let temp_dir = tempfile::tempdir().unwrap();
    let output_dir = temp_dir.path().join("rna_clusters");

    let mut options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;
    options.calling.min_mapping_quality = 20;
    options.calling.min_terminal_soft_clip_ins_len = 0;
    options.junction.min_reads = 4;

    let cluster_set: RNAReadClusterSet = cluster_rna_reads(
        bam_file,
        bai_file.as_str(),
        fasta_file,
        &gene_annotator,
        None,
        None,
        &options,
        2,
        temp_dir.path().to_str().unwrap(),
        output_dir.to_str().unwrap(),
        "scga-mini-rna-001"
    );

    for suffix in ["", "_summary", "_reference_transcripts", "_splice_junctions", "_variants_passed", "_variants_failed", "_template_switch"] {
        let tsv_file = output_dir.join(format!("scga-mini-rna-001_exacto_rna_clusters{}.tsv", suffix));
        assert!(tsv_file.exists(), "{} was not written.", tsv_file.display());
    }

    let loaded_cluster_set: RNAReadClusterSet = load_rna_read_cluster_set(
        output_dir.join("scga-mini-rna-001_exacto_rna_clusters.tsv").to_str().unwrap(),
        output_dir.join("scga-mini-rna-001_exacto_rna_clusters_reference_transcripts.tsv").to_str().unwrap(),
        output_dir.join("scga-mini-rna-001_exacto_rna_clusters_splice_junctions.tsv").to_str().unwrap(),
        output_dir.join("scga-mini-rna-001_exacto_rna_clusters_variants_passed.tsv").to_str().unwrap()
    );
    assert!(!cluster_set.get_clusters().is_empty());
    assert_eq!(loaded_cluster_set.get_cluster_reads_names(), cluster_set.get_cluster_reads_names());
    for cluster in cluster_set.get_clusters() {
        let loaded_splice_junctions: Vec<(u32, u32)> = loaded_cluster_set
            .get_cluster(cluster.get_id())
            .get_splice_junctions()
            .iter()
            .map(|junction| (junction.position_1, junction.position_2))
            .collect();
        let splice_junctions: Vec<(u32, u32)> = cluster
            .get_splice_junctions()
            .iter()
            .map(|junction| (junction.position_1, junction.position_2))
            .collect();
        assert_eq!(loaded_splice_junctions, splice_junctions, "cluster {}", cluster.get_id());
    }
}

/// A DNA row spells its sequence in forward orientation. An RNA call on the reverse strand holds
/// its sequence in read orientation, so the DNA list is looked up by the call's standardized
/// sequence.
///
/// scga-mini-rna-001's TP53 reads (reverse strand) carry the SNV chr17:7674225 C>A, which the
/// reads hold as T, 151 of them. At a read floor of 152 the call fails without the DNA list.
#[test]
fn scga_mini_rna_001_cluster_rna_reads_keeps_dna_listed_snv_on_the_reverse_strand() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bai_file: String = format!("{}.bai", bam_file);
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");

    let mut options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;
    options.filtering.min_reads = 152;

    let dna_variant_records: Vec<DNAVariantRecord> = vec![DNAVariantRecord {
        origin: "somatic".into(),
        variant_id: 1,
        chromosome_1: "chr17".into(),
        position_1: 7_674_224,
        strand_1: "+".into(),
        operation_1: "D".into(),
        chromosome_2: "chr17".into(),
        position_2: 7_674_226,
        strand_2: "+".into(),
        operation_2: "U".into(),
        sequence: "A".into(),
        variant_size: Some(1),
        variant_type: "SNV".into(),
        consensus_read_names: "".into(),
        num_consensus_read_names: 0,
        read_names: "".into(),
        num_read_names: 0
    }];
    let is_snv = |variant_call: &VariantCall| -> bool {
        let op: &GraphOperation = variant_call.get_consensus_graph_operation();
        *op.get_variant_type() == VariantType::SingleNucleotideVariant
            && op.get_position_1() == 7_674_224
            && op.get_position_2() == 7_674_226
            && op.get_standardized_sequence() == "A"
    };

    let mut kept: Vec<bool> = Vec::new();
    for dna_list in [None, Some(&dna_variant_records)] {
        let temp_dir = tempfile::tempdir().unwrap();
        let cluster_set: RNAReadClusterSet = cluster_rna_reads(
            bam_file,
            bai_file.as_str(),
            fasta_file,
            &gene_annotator,
            dna_list,
            None,
            &options,
            2,
            temp_dir.path().to_str().unwrap(),
            temp_dir.path().to_str().unwrap(),
            "test"
        );
        kept.push(cluster_set.get_clusters().iter().any(|cluster| cluster.get_variant_calls().iter().any(|call| is_snv(call))));
    }

    assert_eq!(kept, vec![false, true]);
}


#[test]
fn scga_mini_rna_001_cluster_rna_reads_keys_junctions_by_transcript_strand() {
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::Flags;

    let file = |path: &str| fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join(path)).unwrap().to_str().unwrap().to_owned();
    let fasta = file("references/hg38_chr17-18.fa.gz");
    let gtf = file("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let annotator = Gencode::new_with_defaults(&gtf, "hg38", "v41");
    let temp = tempfile::tempdir().unwrap();
    let bam: String = temp.path().join("mixed-orientation.bam").to_str().unwrap().to_owned();
    let mut reader = bam::io::reader::Builder::default().build_from_path(file("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam")).unwrap();
    let header: noodles_sam::Header = reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(File::create(&bam).unwrap());
    writer.write_header(&header).unwrap();
    for (i, result) in reader.record_bufs(&header).enumerate() {
        let mut record = result.unwrap();
        if i % 2 == 1 {
            record.flags_mut().toggle(Flags::REVERSE_COMPLEMENTED);
        }
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    bam::bai::fs::write(format!("{bam}.bai"), &bam::fs::index(&bam).unwrap()).unwrap();

    let set = cluster_rna_reads(&bam, &format!("{bam}.bai"), &fasta, &annotator, None, None, &ClusterRNAReadsOptions::DEFAULT,
        2, temp.path().to_str().unwrap(), temp.path().to_str().unwrap(), "mixed-orientation");
    // (reads assigned to the cluster, its calls); reads shared with another cell are left out.
    let mut clusters: Vec<(usize, Vec<(u32, u32, u32)>)> = set
        .get_clusters()
        .iter()
        .map(|cluster| (
            cluster.get_read_ids().len() - cluster.get_shared_read_ids().len(),
            cluster
                .get_variant_calls()
                .iter()
                .map(|call| (call.get_consensus_graph_operation().get_position_1(), call.get_consensus_graph_operation().get_position_2(), call.get_num_reads()))
                .collect()
        ))
        .collect();
    clusters.sort();
    assert_eq!(clusters, vec![(154, vec![]), (171, vec![(7_674_224, 7_674_226, 151)])]);
    assert!(set.get_clusters().iter().flat_map(|cluster| cluster.get_splice_junctions().iter()).all(|junction| junction.strand_1 == Strand::Reverse));
}


#[test]
fn scga_mini_rna_001_cluster_rna_reads_writes_header_only_tables_for_no_cluster() {
    let file = |path: &str| fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join(path)).unwrap().to_str().unwrap().to_owned();
    let bam = file("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let fasta = file("references/hg38_chr17-18.fa.gz");
    let gtf = file("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let annotator = Gencode::new_with_defaults(&gtf, "hg38", "v41");
    let temp = tempfile::tempdir().unwrap();
    let directory: &str = temp.path().to_str().unwrap();
    let mut options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;
    options.junction.min_reads = 1_000;
    let set = cluster_rna_reads(&bam, &format!("{bam}.bai"), &fasta, &annotator, None, None, &options,
        2, directory, directory, "empty");
    assert!(set.get_clusters().is_empty());
    let table = |suffix: &str| -> String { format!("{directory}/empty_exacto_rna_clusters{suffix}.tsv") };
    for suffix in ["", "_summary", "_reference_transcripts", "_splice_junctions", "_variants_passed", "_variants_failed", "_template_switch"] {
        let text: String = fs::read_to_string(table(suffix)).unwrap();
        assert!(text.starts_with("cluster_id\t") && text.lines().count() == 1, "{suffix}: {text:?}");
    }
    let loaded = load_rna_read_cluster_set(&table(""), &table("_reference_transcripts"), &table("_splice_junctions"), &table("_variants_passed"));
    assert!(loaded.get_clusters().is_empty());
}


#[test]
fn scga_mini_rna_001_cluster_rna_reads_leaves_out_dna_records_on_absent_contigs() {
    let file = |path: &str| fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join(path)).unwrap().to_str().unwrap().to_owned();
    let bam = file("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let fasta = file("references/hg38_chr17-18.fa.gz");
    let gtf = file("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let annotator = Gencode::new_with_defaults(&gtf, "hg38", "v41");
    let temp = tempfile::tempdir().unwrap();
    let dna_variant_records: Vec<DNAVariantRecord> = vec![DNAVariantRecord {
        variant_id: 1,
        origin: "somatic".into(),
        chromosome_1: "chrEBV".into(),
        position_1: 1_000,
        strand_1: "+".into(),
        operation_1: "D".into(),
        chromosome_2: "chrEBV".into(),
        position_2: 1_002,
        strand_2: "+".into(),
        operation_2: "U".into(),
        sequence: "A".into(),
        variant_size: Some(1),
        variant_type: "SNV".into(),
        consensus_read_names: "".into(),
        num_consensus_read_names: 0,
        read_names: "".into(),
        num_read_names: 0
    }];
    let set = cluster_rna_reads(&bam, &format!("{bam}.bai"), &fasta, &annotator, Some(&dna_variant_records), None, &ClusterRNAReadsOptions::DEFAULT,
        2, temp.path().to_str().unwrap(), temp.path().to_str().unwrap(), "absent-contig");
    // Reads assigned to each cell; reads shared with the other cell are left out.
    let mut sizes: Vec<usize> = set
        .get_clusters()
        .iter()
        .map(|cluster| cluster.get_read_ids().len() - cluster.get_shared_read_ids().len())
        .collect();
    sizes.sort_unstable();
    assert_eq!(sizes, vec![154, 171]);
}


#[test]
fn scga_mini_rna_001_cluster_rna_reads_counts_duplicates_in_depth_and_keeps_mapq_255() {
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::Flags;

    let file = |path: &str| fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join(path)).unwrap().to_str().unwrap().to_owned();
    let fasta = file("references/hg38_chr17-18.fa.gz");
    let gtf = file("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let annotator = Gencode::new_with_defaults(&gtf, "hg38", "v41");
    let temp = tempfile::tempdir().unwrap();
    let bam: String = temp.path().join("duplicate-flag.bam").to_str().unwrap().to_owned();
    let mut reader = bam::io::reader::Builder::default().build_from_path(file("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam")).unwrap();
    let header: noodles_sam::Header = reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(File::create(&bam).unwrap());
    writer.write_header(&header).unwrap();
    for result in reader.record_bufs(&header) {
        let mut record = result.unwrap();
        record.flags_mut().insert(Flags::DUPLICATE);
        *record.mapping_quality_mut() = None;
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    bam::bai::fs::write(format!("{bam}.bai"), &bam::fs::index(&bam).unwrap()).unwrap();

    let mut options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;
    options.calling.min_mapping_quality = 20;
    let set = cluster_rna_reads(&bam, &format!("{bam}.bai"), &fasta, &annotator, None, None, &options,
        2, temp.path().to_str().unwrap(), temp.path().to_str().unwrap(), "duplicate-flag");
    let calls: Vec<(u32, u32, i32)> = set
        .get_clusters()
        .iter()
        .flat_map(|cluster| cluster.get_variant_calls().iter())
        .map(|call| (call.get_consensus_graph_operation().get_position_1(), call.get_num_reads(), call.get_total_depth()))
        .collect();
    assert_eq!(calls, vec![(7_674_224, 151, 305)]);
}


#[test]
fn remove_unspliced_rnas_writes_records_in_input_order() {
    let file = |path: &str| fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join(path)).unwrap().to_str().unwrap().to_owned();
    let bam = file("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let gtf = file("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let annotator = Gencode::new_with_defaults(&gtf, "hg38", "v41");
    let read_names = |path: &str| -> Vec<String> {
        let mut reader = bam::io::reader::Builder::default().build_from_path(path).unwrap();
        reader.read_header().unwrap();
        reader
            .records()
            .map(|record| String::from_utf8(record.unwrap().name().unwrap().to_vec()).unwrap())
            .collect()
    };
    let input: Vec<String> = read_names(&bam);
    for run in 0..3 {
        let temp = tempfile::tempdir().unwrap();
        let output: String = temp.path().join("kept.bam").to_str().unwrap().to_owned();
        let kept = remove_unspliced_rnas(&bam, &format!("{bam}.bai"), &annotator, &output, &format!("{output}.bai"), 2, 0, true);
        let expected: Vec<String> = input.iter().filter(|name| kept.contains(name.as_str())).cloned().collect();
        assert_eq!(expected.len(), 358, "run {run}");
        assert_eq!(read_names(&output), expected, "run {run}");
    }
}


#[test]
fn scga_mini_rna_013_cluster_rna_reads_counts_carriers_in_breakend_depth() {
    let file = |path: &str| fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join(path)).unwrap().to_str().unwrap().to_owned();
    let bam = file("alignment/scga-mini-rna-013-tumor_minimap2_sorted.bam");
    let fasta = file("references/hg38_chr17-18.fa.gz");
    let gtf = file("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let annotator = Gencode::new_with_defaults(&gtf, "hg38", "v41");
    let temp = tempfile::tempdir().unwrap();
    let set = cluster_rna_reads(&bam, &format!("{bam}.bai"), &fasta, &annotator, None, None, &ClusterRNAReadsOptions::PACBIO_HIFI,
        2, temp.path().to_str().unwrap(), temp.path().to_str().unwrap(), "inversion-breakend");
    let mut calls: Vec<(String, u32, usize)> = set
        .get_clusters()
        .iter()
        .flat_map(|cluster| {
            cluster
                .get_variant_calls()
                .iter()
                .map(|call| (call.get_consensus_graph_operation().get_variant_type().as_str().to_string(), call.get_num_reads(), cluster.get_read_ids().len()))
        })
        .collect();
    calls.sort();
    assert_eq!(calls, vec![("BND".to_string(), 81, 81), ("BND".to_string(), 81, 81)]);
}


/// A read that does not tell its cell from the other is a read of both.
///
/// At the PacBio HiFi options scga-mini-rna-001 is one junction cluster phased into two cells by
/// the SNV at chr17:7,674,225: cluster 1 carries it (transcript row 1, `chunk_0000/1-200`) and
/// cluster 2 does not (row 2, `chunk_0000/201-400`). 20 reads are NotCovered at the SNV, so they
/// fit both cells. MEC assigns each to cluster 1, the lower-numbered cell, the 9 reads of row 2
/// that cluster 1 holds among them; each is also a read of cluster 2, which for those 9 is their
/// own transcript. The assignments stay a partition, and the calls are those of the cells alone.
#[test]
fn scga_mini_rna_001_cluster_rna_reads_shares_reads_that_fit_both_cells() {
    let bam_path = Path::new(env!("EXACTO_TEST_DATA")).join("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let bam_full_path = fs::canonicalize(bam_path).unwrap();
    let bam_file: &str = bam_full_path.to_str().unwrap();
    let bai_file: String = format!("{}.bai", bam_file);
    let fasta_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/hg38_chr17-18.fa.gz");
    let fasta_full_path = fs::canonicalize(fasta_path).unwrap();
    let fasta_file: &str = fasta_full_path.to_str().unwrap();
    let gtf_path = Path::new(env!("EXACTO_TEST_DATA")).join("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let gtf_full_path = fs::canonicalize(gtf_path).unwrap();
    let gene_annotator: Gencode = Gencode::new_with_defaults(gtf_full_path.to_str().unwrap(), "hg38", "v41");
    let temp_dir = tempfile::tempdir().unwrap();

    let cluster_set: RNAReadClusterSet = cluster_rna_reads(
        bam_file,
        bai_file.as_str(),
        fasta_file,
        &gene_annotator,
        None,
        None,
        &ClusterRNAReadsOptions::PACBIO_HIFI,
        2,
        temp_dir.path().to_str().unwrap(),
        temp_dir.path().to_str().unwrap(),
        "test"
    );
    let mut clusters: Vec<&RNAReadCluster> = cluster_set.get_clusters();
    clusters.sort_by_key(|cluster| cluster.get_id());
    assert_eq!(clusters.iter().map(|cluster| cluster.get_id()).collect::<Vec<usize>>(), vec![1, 2]);
    let (snv_cell, other_cell): (&RNAReadCluster, &RNAReadCluster) = (clusters[0], clusters[1]);
    assert_eq!(snv_cell.get_splice_junctions(), other_cell.get_splice_junctions(), "one junction cluster");
    let is_row_2 = |read_id: &usize| -> bool {
        let read_name: &str = cluster_set.get_read_names_map().get_by_right(read_id).unwrap();
        let read_number: usize = read_name.split('/').nth(1).unwrap().parse().unwrap();
        read_name.contains("_chunk_0000/") && (201..=400).contains(&read_number)
    };

    // The assignments: a partition of 171 and 154 reads.
    let assigned = |cluster: &RNAReadCluster| -> HashSet<usize> {
        cluster.get_read_ids().difference(cluster.get_shared_read_ids()).copied().collect()
    };
    assert_eq!((assigned(snv_cell).len(), assigned(other_cell).len()), (171, 154));
    assert!(assigned(snv_cell).is_disjoint(&assigned(other_cell)));

    // The 20 reads assigned to the SNV cell are shared with the other, every row-2 read of the
    // SNV cell among them; nothing is shared the other way.
    let shared: &HashSet<usize> = other_cell.get_shared_read_ids();
    assert_eq!(shared.len(), 20);
    assert!(shared.is_subset(&assigned(snv_cell)));
    assert!(snv_cell.get_shared_read_ids().is_empty());
    let row_2_in_snv_cell: HashSet<usize> = assigned(snv_cell).into_iter().filter(|read_id| is_row_2(read_id)).collect();
    assert_eq!(row_2_in_snv_cell.len(), 9);
    assert!(row_2_in_snv_cell.is_subset(shared));

    // Each shared read is NotCovered at the SNV, the site that tells the cells apart.
    let snv_calls: Vec<(u32, u32, usize)> = snv_cell
        .get_variant_calls()
        .iter()
        .map(|variant_call| {
            let op: &GraphOperation = variant_call.get_consensus_graph_operation();
            (op.get_position_1(), op.get_position_2(), variant_call.get_id())
        })
        .collect();
    assert_eq!(snv_calls.iter().map(|(p1, p2, _)| (*p1, *p2)).collect::<Vec<(u32, u32)>>(), vec![(7_674_224, 7_674_226)]);
    let alleles: &HashMap<usize, Allele> = &snv_cell.get_genotypes()[&snv_calls[0].2];
    assert!(shared.iter().all(|read_id| matches!(alleles[read_id], Allele::NotCovered)));
    assert!(other_cell.get_variant_calls().is_empty());
}


/// A run restricted to an allow list judges only the calls on it. scga-mini-rna-001 is two cells
/// of one junction cluster, told apart by the SNV at chr17:7,674,225. With a list that holds the
/// SNV one base off, the run is the one without a list. With an empty list the SNV fails as not
/// listed, so nothing is phased on and the cluster is one cell.
#[test]
fn scga_mini_rna_001_cluster_rna_reads_judges_only_allowed_variants() {
    let file = |path: &str| fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join(path)).unwrap().to_str().unwrap().to_owned();
    let bam = file("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam");
    let fasta = file("references/hg38_chr17-18.fa.gz");
    let gtf = file("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let annotator = Gencode::new_with_defaults(&gtf, "hg38", "v41");
    let snv = |position_1: u32, position_2: u32| -> RNAReadClusterVariantRecord {
        RNAReadClusterVariantRecord {
            cluster_id: 1,
            chromosome_1: "chr17".into(),
            position_1,
            strand_1: "-".into(),
            operation_type_1: "D".into(),
            chromosome_2: "chr17".into(),
            position_2,
            strand_2: "-".into(),
            operation_type_2: "U".into(),
            sequence: "T".into(),
            variant_type: "SNV".into()
        }
    };
    // (reads assigned to each cell, the cell's calls), and the failed calls at the SNV.
    let run = |allowed: Option<&Vec<RNAReadClusterVariantRecord>>| -> (Vec<(usize, Vec<(u32, u32, u32)>)>, Vec<&'static str>) {
        let temp = tempfile::tempdir().unwrap();
        let set = cluster_rna_reads(&bam, &format!("{bam}.bai"), &fasta, &annotator, None, allowed, &ClusterRNAReadsOptions::DEFAULT,
            2, temp.path().to_str().unwrap(), temp.path().to_str().unwrap(), "allowed");
        let mut cells: Vec<(usize, Vec<(u32, u32, u32)>)> = set
            .get_clusters()
            .iter()
            .map(|cluster| (
                cluster.get_read_ids().len() - cluster.get_shared_read_ids().len(),
                cluster
                    .get_variant_calls()
                    .iter()
                    .map(|call| (call.get_consensus_graph_operation().get_position_1(), call.get_consensus_graph_operation().get_position_2(), call.get_num_reads()))
                    .collect()
            ))
            .collect();
        cells.sort();
        let failures: Vec<&'static str> = set
            .get_failed_variant_calls()
            .iter()
            .filter(|failed| failed.graph_operation.get_position_1() == 7_674_224)
            .map(|failed| failed.failure.as_str())
            .collect();
        (cells, failures)
    };

    let (without_list, failures) = run(None);
    assert_eq!(without_list, vec![(154, vec![]), (171, vec![(7_674_224, 7_674_226, 151)])]);
    assert!(failures.is_empty());
    assert_eq!(run(Some(&vec![snv(7_674_225, 7_674_227)])), (without_list, vec![]));

    let (cells, failures) = run(Some(&vec![]));
    assert_eq!(cells, vec![(325, vec![])]);
    assert_eq!(failures, vec!["not_in_allowed_list"]);
}


/// The reads of the ASPA-WSCD1 fusion sequenced from either end are one transcript: every record
/// of every odd-numbered read of scga-mini-rna-007 is turned to the other strand, as a library of
/// unoriented reads holds them. The fusion's junction is then spelled with both strands flipped by
/// half its reads, and its clipped reads carry the ASPA arm at either end; it is still one call,
/// one cluster of the 109 fusion reads and the five clipped ones, and no insertion.
#[test]
fn scga_mini_rna_007_cluster_rna_reads_keeps_fusion_whole_across_read_orientations() {
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::Flags;

    let file = |path: &str| fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join(path)).unwrap().to_str().unwrap().to_owned();
    let fasta = file("references/hg38_chr17-18.fa.gz");
    let gtf = file("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let annotator = Gencode::new_with_defaults(&gtf, "hg38", "v41");
    let temp = tempfile::tempdir().unwrap();
    let bam: String = temp.path().join("mixed-orientation.bam").to_str().unwrap().to_owned();
    let mut reader = bam::io::reader::Builder::default().build_from_path(file("alignment/scga-mini-rna-007-tumor_minimap2_sorted.bam")).unwrap();
    let header: noodles_sam::Header = reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(File::create(&bam).unwrap());
    writer.write_header(&header).unwrap();
    for result in reader.record_bufs(&header) {
        let mut record = result.unwrap();
        let name: String = record.name().unwrap().to_string();
        let number: usize = name.split('/').nth(1).unwrap().parse().unwrap();
        if number % 2 == 1 {
            record.flags_mut().toggle(Flags::REVERSE_COMPLEMENTED);
        }
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    bam::bai::fs::write(format!("{bam}.bai"), &bam::fs::index(&bam).unwrap()).unwrap();

    let set = cluster_rna_reads(&bam, &format!("{bam}.bai"), &fasta, &annotator, None, None, &ClusterRNAReadsOptions::PACBIO_HIFI,
        2, temp.path().to_str().unwrap(), temp.path().to_str().unwrap(), "mixed-orientation");
    assert_eq!(set.get_clusters().len(), 3);
    let fusions: Vec<&VariantCall> = set
        .get_clusters()
        .iter()
        .flat_map(|cluster| cluster.get_variant_calls())
        .filter(|call| *call.get_consensus_graph_operation().get_variant_type() == VariantType::FusionGene)
        .collect();
    assert_eq!(fusions.len(), 1);
    assert_eq!(fusions[0].get_num_reads(), 114);
    assert!(!set.get_clusters().iter().flat_map(|cluster| cluster.get_variant_calls())
        .any(|call| *call.get_consensus_graph_operation().get_variant_type() == VariantType::Insertion));
}


/// The reads of the TP53 circular RNA sequenced from either end key one chain: every record of
/// every odd-numbered read of scga-mini-rna-011 is turned to the other strand. A read from the other
/// end names the back-splice by the same two bases with both strands flipped, so the circle is one
/// cluster of 151 reads beside the linear transcript, as with every read in transcript orientation.
#[test]
fn scga_mini_rna_011_cluster_rna_reads_keeps_circle_whole_across_read_orientations() {
    use noodles_sam::alignment::io::Write as _;
    use noodles_sam::alignment::record::Flags;

    let file = |path: &str| fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join(path)).unwrap().to_str().unwrap().to_owned();
    let fasta = file("references/hg38_chr17-18.fa.gz");
    let gtf = file("references/gencode.v41.annotation.chr17-18.gtf.gz");
    let annotator = Gencode::new_with_defaults(&gtf, "hg38", "v41");
    let temp = tempfile::tempdir().unwrap();
    let bam: String = temp.path().join("mixed-orientation.bam").to_str().unwrap().to_owned();
    let mut reader = bam::io::reader::Builder::default().build_from_path(file("alignment/scga-mini-rna-011-tumor_minimap2_sorted.bam")).unwrap();
    let header: noodles_sam::Header = reader.read_header().unwrap();
    let mut writer = bam::io::Writer::new(File::create(&bam).unwrap());
    writer.write_header(&header).unwrap();
    for result in reader.record_bufs(&header) {
        let mut record = result.unwrap();
        let name: String = record.name().unwrap().to_string();
        let number: usize = name.split('/').nth(1).unwrap().parse().unwrap();
        if number % 2 == 1 {
            record.flags_mut().toggle(Flags::REVERSE_COMPLEMENTED);
        }
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    bam::bai::fs::write(format!("{bam}.bai"), &bam::fs::index(&bam).unwrap()).unwrap();

    let set = cluster_rna_reads(&bam, &format!("{bam}.bai"), &fasta, &annotator, None, None, &ClusterRNAReadsOptions::PACBIO_HIFI,
        2, temp.path().to_str().unwrap(), temp.path().to_str().unwrap(), "mixed-orientation");
    assert_eq!(set.get_clusters().len(), 2);
    let circles: Vec<(usize, u32)> = set
        .get_clusters()
        .iter()
        .flat_map(|cluster| cluster.get_variant_calls().iter().map(move |call| (cluster.get_read_ids().len(), call)))
        .filter(|(_, call)| {
            let operation: &GraphOperation = call.get_consensus_graph_operation();
            (operation.get_position_1(), operation.get_position_2()) == (7_670_611, 7_674_288)
        })
        .map(|(cluster_reads, call)| (cluster_reads, call.get_num_reads()))
        .collect();
    assert_eq!(circles, vec![(151, 151)]);
}
