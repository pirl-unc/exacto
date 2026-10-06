use exacto_core::prelude::write_tsv_table;
use crate::prelude::load_dna_variant_records;
use std::fs;
use std::path::Path;

use super::*;


#[test]
fn every_dataframe_has_the_columns_of_its_table_file_also_when_empty() {
    let temp_dir = tempfile::tempdir().unwrap();
    let path = temp_dir.path().join("table.tsv");
    // An empty table is written as its header line alone.
    let header = |path: &Path| -> Vec<String> {
        let text: String = fs::read_to_string(path).unwrap();
        assert_eq!(text.lines().count(), 1, "{text}");
        text.trim_end().split('\t').map(String::from).collect()
    };
    let columns = |dataframe: DataFrame| -> Vec<String> {
        dataframe.get_column_names().iter().map(|name| name.to_string()).collect()
    };

    write_tsv_table(Vec::<AssembledTranscriptRecord>::new(), &path).unwrap();
    assert_eq!(header(&path), columns(assembled_transcript_records_to_dataframe(Vec::new())));
    write_tsv_table(Vec::<AssembledTranscriptExonRecord>::new(), &path).unwrap();
    assert_eq!(header(&path), columns(assembled_transcript_exon_records_to_dataframe(Vec::new())));
    write_tsv_table(Vec::<AssembledTranscriptSpliceJunctionRecord>::new(), &path).unwrap();
    assert_eq!(header(&path), columns(assembled_transcript_splice_junction_records_to_dataframe(Vec::new())));
    write_tsv_table(Vec::<AssembledTranscriptReferenceTranscriptMatchRecord>::new(), &path).unwrap();
    assert_eq!(header(&path), columns(assembled_transcript_reference_transcript_match_records_to_dataframe(Vec::new())));
    write_tsv_table(Vec::<AssembledTranscriptFilterStatusRecord>::new(), &path).unwrap();
    assert_eq!(header(&path), columns(assembled_transcript_filter_status_records_to_dataframe(Vec::new())));
    write_tsv_table(Vec::<AssembledTranscriptNonsenseMediatedDecayRecord>::new(), &path).unwrap();
    assert_eq!(header(&path), columns(assembled_transcript_nonsense_mediated_decay_records_to_dataframe(Vec::new())));
    write_tsv_table(Vec::<AssembledTranscriptModelAlignmentRecord>::new(), &path).unwrap();
    assert_eq!(header(&path), columns(assembled_transcript_model_alignment_records_to_dataframe(Vec::new())));
    write_tsv_table(Vec::<AssembledTranscriptVariantRecord>::new(), &path).unwrap();
    assert_eq!(header(&path), columns(assembled_transcript_variant_records_to_dataframe(Vec::new())));
    write_tsv_table(Vec::<DNAVariantRecord>::new(), &path).unwrap();
    assert_eq!(header(&path), columns(dna_variant_records_to_dataframe(Vec::new())));

    // The loaders read an empty table back as no rows.
    assert!(load_dna_variant_records(path.to_str().unwrap()).is_empty());
}
