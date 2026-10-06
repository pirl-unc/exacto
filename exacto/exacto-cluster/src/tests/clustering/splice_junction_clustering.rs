use bimap::BiMap;
use exacto_caller::prelude::*;
use exacto_core::prelude::{create_chromosome_names_map, index_bam_records, load_bincode_temp_files, FastaMap, Gencode};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;

use super::*;

use crate::prelude::ClusterRNAReadsOptions;


fn data(relative: &str) -> String {
    fs::canonicalize(Path::new(env!("EXACTO_TEST_DATA")).join(relative)).unwrap().to_str().unwrap().to_string()
}

/// Every read of `bam_file` modelled as the pipeline models it, plus the maps the clustering
/// takes. The BAM header declares chr17 and chr18, so the chr17-18 FASTA and GENCODE slice are
/// the reference.
fn transcript_models(
    bam_file: &str
) -> (Vec<TranscriptModel>, BiMap<Box<str>, usize>, BiMap<Box<str>, u16>, Gencode, FastaMap) {
    let fasta_map = FastaMap::new(&data("references/hg38_chr17-18.fa.gz"));
    let gene_annotator = Gencode::new_with_defaults(&data("references/gencode.v41.annotation.chr17-18.gtf.gz"), "hg38", "v41");
    let (record_positions_map, read_names_map) = index_bam_records(bam_file, true, 2);
    let chromosome_names_map = create_chromosome_names_map(bam_file);
    let mut read_ids: Vec<usize> = read_names_map.right_values().copied().collect();
    read_ids.sort_unstable();
    let options: ClusterRNAReadsOptions = ClusterRNAReadsOptions::DEFAULT;
    let directory = tempfile::tempdir().unwrap();
    let (temp_file, summaries, _) = characterize_rna_reads(
        bam_file, 
        &read_ids, 
        &read_names_map, 
        &record_positions_map,
        &chromosome_names_map, 
        &fasta_map,
        &gene_annotator, 
        2,
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
        directory.path().to_str().unwrap()
    );
    let models: Vec<TranscriptModel> = load_transcript_models(&temp_file, &summaries);
    (models, read_names_map, chromosome_names_map, gene_annotator, fasta_map)
}

/// scga-mini-rna-001's spliced and unspliced reads of TP53: the unspliced reads form one
/// gene-keyed cluster with no junctions, and no read is duplicated. Two reads are lost:
/// `chunk_0000/91` and `chunk_0000/188` each hold a deletion against a junction, which their
/// models spell as a junction one base off that no other read spells, and at a maximum FPR of 1 a
/// novel splice junction still needs 2 reads.
#[test]
fn scga_mini_rna_001_cluster_transcript_models_by_splice_junctions_keys_unspliced_reads_by_gene() {
    let (models, read_names_map, chromosome_names_map, gene_annotator, fasta_map) =
        transcript_models(&data("alignment/scga-mini-rna-001-tumor_minimap2_sorted.bam"));
    let clusters = cluster_transcript_models_by_splice_junctions(&models, &SpliceJunctionAnnotationIndex::new(&gene_annotator, &chromosome_names_map), &chromosome_names_map, &fasta_map, 1, 0.01, 1.0, 0, None);

    let unspliced: Vec<&SpliceJunctionCluster> = clusters.iter().filter(|c| c.splice_junctions.is_empty()).collect();
    assert_eq!(unspliced.len(), 1);
    assert_eq!(unspliced[0].read_ids.len(), 44);
    assert!(unspliced[0].reference_gene_transcript_ids.iter().all(|(_, transcript_id)| transcript_id.is_empty()));

    let mut seen: HashMap<usize, usize> = HashMap::new();
    for cluster in clusters.iter() {
        for read_id in cluster.read_ids.iter() {
            *seen.entry(*read_id).or_insert(0) += 1;
        }
    }
    let mut unclustered: Vec<&str> = read_names_map
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


/// A read the DNA list vouches for exempts its whole chain from the junction read floor.
///
/// scga-mini-rna-013's 3-read phantom chain, which spells the junction into the inverted segment
/// as 7,673,700-7,669,670 (see the pipeline test
/// `scga_mini_rna_013_cluster_rna_reads_drops_phantom_fork_below_junction_read_floor`), dies at a
/// floor of 4 on its own. Flagging one of its reads as a DNA carrier brings the chain back at the
/// same floor, while a carrier on the parent changes nothing: the floor exists to suppress noise,
/// and DNA evidence is exactly what distinguishes a rare real transcript from it.
#[test]
fn scga_mini_rna_013_cluster_transcript_models_by_splice_junctions_exempts_dna_carriers_from_read_floor() {
    let (models, _, chromosome_names_map, gene_annotator, fasta_map) =
        transcript_models(&data("alignment/scga-mini-rna-013-tumor_minimap2_sorted.bam"));
    let annotation_index = SpliceJunctionAnnotationIndex::new(&gene_annotator, &chromosome_names_map);
    let cluster = |min_reads: usize, carriers: Option<&HashSet<usize>>| -> Vec<SpliceJunctionCluster> {
        cluster_transcript_models_by_splice_junctions(&models, &annotation_index, &chromosome_names_map, &fasta_map, min_reads, 0.01, 1.0, 0, carriers)
    };
    // The reads of the chains holding a junction.
    let reads_with = |clusters: &Vec<SpliceJunctionCluster>, junction: (u32, u32)| -> Vec<HashSet<usize>> {
        clusters
            .iter()
            .filter(|cluster| cluster.splice_junctions.iter().any(|splice_junction| (splice_junction.position_1, splice_junction.position_2) == junction))
            .map(|cluster| cluster.read_ids.iter().copied().collect())
            .collect()
    };
    let phantom: (u32, u32) = (7_673_700, 7_669_670);
    let parent: (u32, u32) = (7_673_700, 7_669_660);

    // Guard: at a floor of 3 the phantom chain forms its own cluster; at 4 it does not.
    let admitted: Vec<HashSet<usize>> = reads_with(&cluster(3, None), phantom);
    assert_eq!(admitted.iter().map(|read_ids| read_ids.len()).collect::<Vec<usize>>(), vec![3]);
    let floored = cluster(4, None);
    assert!(reads_with(&floored, phantom).is_empty(), "The 3-read phantom chain should fall below a floor of 4.");

    // One phantom read flagged as a DNA carrier brings its whole chain back at the same floor.
    let phantom_carrier: usize = *admitted[0].iter().min().unwrap();
    assert_eq!(reads_with(&cluster(4, Some(&HashSet::from([phantom_carrier]))), phantom), admitted);

    // A carrier on the parent chain, which clears the floor anyway, rescues nothing else.
    let parent_carrier: usize = *reads_with(&floored, parent)[0].iter().min().unwrap();
    assert!(reads_with(&cluster(4, Some(&HashSet::from([parent_carrier]))), phantom).is_empty());
}
