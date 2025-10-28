use exacto_caller::prelude::*;
use exacto_core::prelude::*;
use flate2::read::GzDecoder;
use noodles_fastq as fastq;
use std::fs;
use std::fs::File;
use std::io::{BufReader,Read};
use std::path::Path;
use tempfile::NamedTempFile;

use crate::prelude::*;


#[test]
fn test_translation_1() {
    // Step 1. Read the RNA FASTQ file
    let fastq_path = Path::new("src/tests/data/fastq/sample200normal_long_read_rna.fastq.gz");
    let fastq_full_path = fs::canonicalize(fastq_path).unwrap();
    let gzipped = is_gzipped(fastq_full_path.to_str().unwrap());
    let file = File::open(fastq_full_path.to_str().unwrap()).expect("Unable to open FASTQ file");
    let reader: Box<dyn Read> = if gzipped {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let buffered_reader = BufReader::new(reader);
    let mut fastq_reader = fastq::Reader::new(buffered_reader);
    let mut rnas: Vec<RNA> = Vec::new();
    for result in fastq_reader.records() {
        match result {
            Ok(record) => {
                let sequence_result = String::from_utf8(record.sequence().to_vec());
                let sequence: String = match sequence_result {
                    Ok(seq) => seq,
                    Err(e) => {
                        panic!("Error converting sequence to UTF-8: {}", e);
                    }
                };
                let rna: RNA = RNA::new(
                    record.name().to_string().into_boxed_str(),
                    sequence.into_boxed_str(),
                );
                rnas.push(rna);
            }
            Err(e) => {
                panic!("Error reading record: {}", e);
            }
        }
    }

    // Step 2. Translate the RNA sequences
    let translation_set: TranslationSet = translate_rnas(
        rnas,
        2
    );

    // Step 3. Check for validity
    let mut found: bool = false;
    for translation in translation_set.translations.iter() {
        if translation.rna.id == "m64012_817037_637278/74/ccs".into() {
            if translation.get_longest_orf_peptide().sequence == "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD*".into() {
                found = true;
            }
            assert!(translation.get_peptides_count() == 38);
        }
    }
    assert!(found);
    assert!(translation_set.translations.len() == 200);
}

#[test]
fn test_translation_2() {
    let tsv_path_1 = Path::new("src/tests/data/tsv/rna-100-tumor_minimap2_mdtagged_sorted_exacto_transcript_structures.tsv");
    let tsv_full_path_1 = fs::canonicalize(tsv_path_1).unwrap();
    let tsv_file_1: &str = tsv_full_path_1.to_str().unwrap();
    let tsv_path_2 = Path::new("src/tests/data/tsv/rna-100-tumor_minimap2_mdtagged_sorted_exacto_rna_variant_calls.tsv");
    let tsv_full_path_2 = fs::canonicalize(tsv_path_2).unwrap();
    let tsv_file_2: &str = tsv_full_path_2.to_str().unwrap();
    let tsv_path_3 = Path::new("src/tests/data/tsv/rna-100_dna-001_integration.tsv");
    let tsv_full_path_3 = fs::canonicalize(tsv_path_3).unwrap();
    let tsv_file_3: &str = tsv_full_path_3.to_str().unwrap();

    let df_transcript_structures = read_tsv_file(tsv_file_1);
    let rna_variant_call_set = RNAVariantCallSet::read_tsv_file(tsv_file_2);
    let df_integrated_variants = read_tsv_file(tsv_file_3);

    let primary_structure_set = translate_transcript_structures(
        &df_transcript_structures,
        &rna_variant_call_set,
        &df_integrated_variants,
        TranslationStrategy::LongestORF,
        1
    );

    assert!(primary_structure_set.primary_structures.len() == 6);

    let output_fasta_file: NamedTempFile = NamedTempFile::new().unwrap();
    let output_tsv_file: NamedTempFile = NamedTempFile::new().unwrap();

    primary_structure_set.to_fasta_file(output_fasta_file.path().to_str().unwrap());
    primary_structure_set.to_tsv_file(output_tsv_file.path().to_str().unwrap());
}
