use std::io::Write;
use std::path::Path;
use tempfile::NamedTempFile;

use super::*;


#[test]
fn load_dna_variant_records_returns_matches_1() {
    // A bare-bones file: the grammar columns only, with ids that are not 1, 2, ...
    let mut file: NamedTempFile = NamedTempFile::new().unwrap();
    write!(
        file,
        "origin\tvariant_id\tchromosome_1\tposition_1\tstrand_1\toperation_1\tchromosome_2\tposition_2\tstrand_2\toperation_2\tsequence\n\
         somatic\t5\tchr1\t528516\t+\tD\tchr1\t528518\t+\tU\tT\n\
         somatic\t42\tchr7\t100000\t+\tD\tchr7\t100002\t+\tU\tG"
    ).unwrap();
    let records: Vec<DNAVariantRecord> = load_dna_variant_records(file.path().to_str().unwrap());

    assert_eq!(records.len(), 2, "expected 2 records");

    // Row 1: variant_id is preserved from the file (5), not renumbered.
    let r1 = &records[0];
    assert_eq!(r1.variant_id, 5, "variant_id must be taken from the file, not renumbered");
    assert_eq!(&*r1.chromosome_1, "chr1");
    assert_eq!(r1.position_1, 528516);
    assert_eq!(&*r1.strand_1, "+");
    assert_eq!(&*r1.operation_1, "D");
    assert_eq!(&*r1.chromosome_2, "chr1");
    assert_eq!(r1.position_2, 528518);
    assert_eq!(&*r1.strand_2, "+");
    assert_eq!(&*r1.operation_2, "U");
    assert_eq!(&*r1.sequence, "T");
    // Columns absent from the bare-bones grammar default to empty / 0 / no size.
    assert_eq!(r1.variant_size, None);
    assert_eq!(&*r1.variant_type, "");
    assert_eq!(&*r1.consensus_read_names, "");
    assert_eq!(r1.num_consensus_read_names, 0);
    assert_eq!(&*r1.read_names, "");
    assert_eq!(r1.num_read_names, 0);

    // Row 2: a different, non-sequential id is likewise preserved.
    let r2 = &records[1];
    assert_eq!(r2.variant_id, 42);
    assert_eq!(&*r2.chromosome_1, "chr7");
    assert_eq!(r2.position_1, 100000);
    assert_eq!(&*r2.sequence, "G");
}

#[test]
fn load_dna_variant_records_returns_matches_2() {
    // A table written by call-germline-dna-vars: the simulated SNV of scga-mini-dna-001.
    let path = Path::new(env!("EXACTO_TEST_DATA")).join("variant_calling/dna/scga-mini-dna-001-tumor_exacto_germline_dna_variants.tsv");
    let records: Vec<DNAVariantRecord> = load_dna_variant_records(path.to_str().unwrap());

    assert_eq!(records.len(), 1);
    let r = &records[0];
    assert_eq!(r.variant_id, 1);
    assert_eq!(&*r.origin, "germline");
    assert_eq!(&*r.chromosome_1, "chr17");
    assert_eq!(r.position_1, 7674224);
    assert_eq!(&*r.chromosome_2, "chr17");
    assert_eq!(r.position_2, 7674226);
    assert_eq!(&*r.sequence, "A");
    assert_eq!(r.variant_size, Some(1));
    assert_eq!(&*r.variant_type, "SNV");
    // Read-name columns are populated (not defaulted) when present, and each count is the
    // number of names in its column.
    assert!(r.num_consensus_read_names > 0);
    assert_eq!(r.num_consensus_read_names as usize, r.consensus_read_names.split(';').count());
    assert!(r.num_read_names > 0);
    assert_eq!(r.num_read_names as usize, r.read_names.split(';').count());
}

#[test]
fn load_assembled_transcript_variant_records_returns_matches() {
    let mut file: NamedTempFile = NamedTempFile::new().unwrap();
    write!(
        file,
        "variant_id\tassembled_transcript_name\treference_gene_name\treference_transcript_id\tchromosome_1\tposition_1\tstrand_1\toperation_1\tchromosome_2\tposition_2\tstrand_2\toperation_2\tvariant_size\tvariant_type\tsequence\tread_start\tread_end\n\
         1\tcid_1_1\tTP53\tENST00000269305.9\tchr17\t7674220\t+\tD\tchr17\t7674222\t+\tU\t1\tSNV\tA\t100\t101\n\
         2\tcid_1_1\tTP53;BRCA1\tENST00000269305.9;ENST00000357654.9\tchr17\t7674900\t+\tD\tchr17\t43044295\t+\tU\t\tFusion\t\t200\t201\n"
    ).unwrap();
    let records: Vec<AssembledTranscriptVariantRecord> =
        load_assembled_transcript_variant_records(file.path().to_str().unwrap());

    assert_eq!(records.len(), 2);
    let r1 = &records[0];
    assert_eq!(r1.variant_id, 1);
    assert_eq!(&*r1.assembled_transcript_name, "cid_1_1");
    assert_eq!(&*r1.chromosome_1, "chr17");
    assert_eq!(r1.position_1, 7674220);
    assert_eq!(&*r1.variant_type, "SNV");
}
