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


use bio::io::bed;
use noodles_bgzf as bgzf;
use noodles_fasta as fasta;
use polars::prelude::*;
use std::fs::File;
use std::str;
use std::sync::Arc;
use std::io::{BufReader, Read};
use std::path::Path;
use tempfile::TempPath;


pub fn file_exists(file: &str) -> bool {
    Path::new(file).exists() && Path::new(file).is_file()
}


/// True only for BGZF (block gzip, as written by `bgzip`): gzip magic + FEXTRA flag +
/// the BGZF `BC` extra subfield. Plain gzip (from `gzip`) also starts `1F 8B` but lacks
/// the `BC` subfield and is *not* randomly seekable, so this returns false for it.
pub fn is_bgzipped(file: &str) -> bool {
    let mut header = [0u8; 16];
    let mut file = match File::open(file) {
        Ok(file) => file,
        Err(_) => return false
    };
    if file.read_exact(&mut header).is_err() {
        return false;                                       // shorter than a BGZF header
    }
    header[0] == 0x1F                                       // gzip magic
        && header[1] == 0x8B
        && header[2] == 0x08                                // CM = deflate
        && (header[3] & 0x04) != 0                          // FLG.FEXTRA set
        && header[12] == b'B'                               // BGZF 'BC' subfield id
        && header[13] == b'C'
}


/// Check if file is gzipped.
pub fn is_gzipped(file_path: &str) -> bool {
    if let Ok(mut file) = File::open(file_path) {
        let mut buffer = [0u8; 2];
        if file.read_exact(&mut buffer).is_ok() {
            return buffer == [0x1F, 0x8B];
        }
    }
    false
}


pub fn load_bincode_temp_files<T: serde::de::DeserializeOwned>(temp_paths: &[TempPath]) -> Vec<T> {
    let mut items: Vec<T> = Vec::new();
    for temp_path in temp_paths.iter() {
        let file: File = File::open(temp_path).unwrap();
        let mut reader: BufReader<File> = BufReader::new(file);
        let chunk_items: Vec<T> = bincode::deserialize_from(&mut reader)
            .expect("Failed to deserialize data.");
        items.extend(chunk_items);
    }
    items
}


/// Read a BED file.
pub fn read_bed_file(bed_file: &str) -> Vec<bed::Record> {
    let file = File::open(bed_file).unwrap();
    let reader = BufReader::new(file);
    let mut bed_reader = bed::Reader::new(reader);
    let mut records: Vec<bed::Record> = Vec::new();
    for record in bed_reader.records() {
        let record = record.unwrap();
        records.push(record);
    }
    records
}


/// Read a TSV file.
///
/// tsv_file can be gzipped or not.
pub fn read_tsv_file(tsv_file: &str) -> DataFrame {
    let parse_options = CsvParseOptions {
        separator: b'\t',
        ..Default::default()
    };
    let options = CsvReadOptions {
        parse_options: Arc::new(parse_options),
        ..Default::default()
    };
    let df: DataFrame = options
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap();
    df
}


pub fn read_fasta_file(fasta_file: &str) -> Vec<(Box<str>, Box<str>)> {
    let mut sequences: Vec<(Box<str>, Box<str>)> = Vec::new();

    // Check if file is gzip compressed by reading magic bytes
    let is_gzipped: bool = {
        let mut file: File = File::open(fasta_file).unwrap();
        let mut magic = [0; 2];
        match file.read_exact(&mut magic) {
            Ok(()) => magic == [0x1F, 0x8B],
            Err(_) => false
        }
    };

    if is_gzipped {
        // Handle gzipped FASTA file
        let mut reader = File::open(fasta_file)
            .map(bgzf::io::Reader::new)
            .map(fasta::io::Reader::new).unwrap();
        for result in reader.records() {
            let record = result.unwrap();
            let sequence_id: Box<str> = str::from_utf8(record.name()).unwrap().to_string().into_boxed_str();
            let sequence: Box<str> = str::from_utf8(record.sequence().as_ref()).unwrap()
                .to_string()
                .into_boxed_str();
            sequences.push((sequence_id, sequence));
        }
    } else {
        // Handle uncompressed FASTA file
        let mut reader = File::open(fasta_file)
            .map(BufReader::new)
            .map(fasta::io::Reader::new).unwrap();
        for result in reader.records() {
            let record = result.unwrap();
            let sequence_id: Box<str> = str::from_utf8(record.name()).unwrap().to_string().into_boxed_str();
            let sequence: Box<str> = str::from_utf8(record.sequence().as_ref()).unwrap()
                .to_string()
                .into_boxed_str();
            sequences.push((sequence_id, sequence));
        }
    }

    sequences
}
