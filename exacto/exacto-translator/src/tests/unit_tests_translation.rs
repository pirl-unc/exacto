use exacto_util::common::files::is_gzipped;
use flate2::read::GzDecoder;
use noodles_fastq as fastq;
use std::fs;
use std::fs::File;
use std::io::{BufReader,Read};
use std::path::Path;

use crate::algorithms::translation::translate;
use crate::structs::rna::RNA;
use crate::structs::translation_set::TranslationSet;


#[test]
fn translate_rna_fastq_1() {
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
    let translation_set: TranslationSet = translate(
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
            assert!(translation.get_peptides_count() == 19);
        }
    }
    assert!(found);
    assert!(translation_set.translations.len() == 200);
}

