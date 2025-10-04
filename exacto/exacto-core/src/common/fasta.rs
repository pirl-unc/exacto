// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


use flate2::read::MultiGzDecoder;
use noodles_core::{Position, Region};
use noodles_fasta as fasta;
use noodles_fasta::{fai, Record};
use noodles_fasta::io::indexed_reader::Builder;
use noodles_fasta::io::Reader as FastaReader;
use noodles_fasta::record::{Definition, Sequence};
use std::error::Error;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::Command;


/// Create a FAI index file for a FASTA file.
///
/// # Arguments
/// * `fasta_file`: FASTA file.
///
/// # Notes
/// - For `.fa`, this function writes `<file>.fai` using noodles-fasta.
/// - For `.fa.gz`, this function shells out to `samtools faidx` to produce `<file>.fa.gz.fai` and `<file>.fa.gz.gzi`.
///
/// # Requirements
/// - `.fa` path: no extra tools required.
/// - `.fa.gz` path: `samtools` must be in PATH (needed to create both .fai and .gzi correctly).
pub fn create_fai_file<P: AsRef<Path>>(fasta_file: P) -> Result<PathBuf, Box<dyn Error>> {
    let fasta_path = fasta_file.as_ref();

    // Case 1: Uncompressed FASTA (.fa / .fasta)
    let is_gz = fasta_path
        .extension()
        .map(|e| e == "gz")
        .unwrap_or(false);

    if !is_gz {
        let index = fasta::fs::index(fasta_path)?;
        let fai_out = fasta_path.with_extension(
            format!(
                "{}.fai",
                fasta_path
                    .extension()
                    .and_then(|s| s.to_str())
                    .unwrap_or_default()
            )
        );
        let fai_out = {
            let mut s = fasta_path.as_os_str().to_owned();
            s.push(".fai");
            PathBuf::from(s)
        };
        fai::fs::write(&fai_out, &index)?;
        return Ok(fai_out);
    }

    // Case 2: Gzipped FASTA (.fa.gz)
    let status = Command::new("samtools")
        .args(["faidx", fasta_path.to_str().ok_or("Invalid UTF-8 path")?])
        .status();
    match status {
        Ok(s) if s.success() => {
            let mut fai_out = fasta_path.as_os_str().to_owned();
            fai_out.push(".fai");
            Ok(PathBuf::from(fai_out))
        }
        Ok(s) => Err(format!(
            "samtools faidx exited with status {} when indexing {:?}. \
             Ensure the file is BGZF-compressed (bgzip) or decompress to .fa and try again.",
            s.code().unwrap_or(-1),
            fasta_path
        ).into()),
        Err(e) => Err(format!(
            "Could not run `samtools faidx` for gzipped FASTA {:?}: {}. \
             Either install samtools (recommended for .fa.gz) or decompress to .fa and use noodles_fasta::fs::index.",
            fasta_path, e
        ).into()),
    }
}

/// Checks if an index files exists for a FASTA file.
///
/// # Arguments
/// * `fasta_file`: FASTA file.
///
/// # Returns
/// * true if an index file exists.
/// * false otherwise.
pub fn fasta_index_exists(fasta_file: &str) -> bool {
    let fasta_path = Path::new(fasta_file);
    let possible_suffixes = [
        ".fai",
        ".fa.fai",
        ".fasta.fai",
        ".fa.gz.fai",
        ".fasta.gz.fai",
    ];

    possible_suffixes.iter().any(|suffix| {
        let candidate = fasta_path.with_file_name(
            format!("{}{}", fasta_path.file_name().unwrap().to_string_lossy(), suffix)
        );
        candidate.exists()
    })
}

/// Fetches a sequence from a FASTA file.
///
/// # Arguments
/// * `sequence_id`: Sequence ID.
/// * `start`: Start position.
/// * `end`: End position.
/// * `fasta_file`: FASTA file.
///
/// # Returns
/// * Sequence.
pub fn get_fasta_sequence(
    sequence_id: &str,
    start: usize,
    end: usize,
    fasta_file: &str
) -> Box<str> {
    if fasta_index_exists(fasta_file) == false {
        create_fai_file(fasta_file);
    }

    let mut fasta_reader = Builder::default()
        .build_from_path(fasta_file)
        .unwrap();

    let position_start = Position::try_from(start).unwrap();
    let position_end = Position::try_from(end).unwrap();
    let region = Region::new(sequence_id, position_start..=position_end);
    let ref_record = fasta_reader.query(&region).unwrap();
    let ref_sequence_bytes: &[u8] = ref_record.sequence().as_ref();
    let sequence: &str = std::str::from_utf8(ref_sequence_bytes).expect("Failed to convert sequence to UTF-8");
    sequence.into()
}

/// Fetches all sequence IDs in a FASTA file.
///
/// # Arguments
/// * `fasta_file`: FASTA file.
///
/// # Returns
/// * Vector of tuples (sequence ID, sequence length).
pub fn get_fasta_sequence_ids(fasta_file: &str) -> Vec<(Box<str>, usize)> {
    let file = File::open(fasta_file).expect("Failed to open FASTA file");
    let reader: Box<dyn Read> = if fasta_file.ends_with(".gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let buf_reader = BufReader::new(reader);
    let mut fasta_reader = FastaReader::new(buf_reader);

    fasta_reader
        .records()
        .map(|res| res.expect("Failed to read FASTA record"))
        .map(|record| {
            let name: String = std::str::from_utf8(record.definition().name().as_ref())
                .expect("Invalid UTF-8 in FASTA header")
                .to_string();
            let len = record.sequence().len();  // total bases for this record
            (name.into_boxed_str(), len)
        })
        .collect()
}

/// Writes a FASTA file.
///
/// # Arguments
/// * `sequences`: Vector of tuples (sequence ID, sequence).
/// * `output_fasta_file`: Output FASTA file.
pub fn write_fasta_file(sequences: &Vec<(Box<str>, Box<str>)>, output_fasta_file: &str) {
    let mut writer = fasta::io::writer::Builder::default()
        .set_line_base_count(80)
        .build_from_path(output_fasta_file)
        .unwrap();

    for (name, sequence) in sequences {
        let definition: Definition = Definition::new(name.as_ref(), None);
        let record: Record = Record::new(definition, Sequence::from(sequence.as_bytes().to_owned()));
        writer.write_record(&record).unwrap();
    }

    create_fai_file(output_fasta_file);
}
