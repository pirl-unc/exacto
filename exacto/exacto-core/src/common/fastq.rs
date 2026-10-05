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
use noodles_bgzf as bgzf;
use noodles_fastq as fastq;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Seek, SeekFrom, Write};

use crate::prelude::*;


/// Where a read's record starts, in whichever coordinate space the FASTQ needs: a raw
/// byte offset for an uncompressed file, or a BGZF virtual position for a bgzip file.
/// Pairs with a [`FastqReader`] opened on the same file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FastqOffset {
    /// Byte offset from the start of an uncompressed FASTQ.
    Plain(u64),
    /// BGZF virtual position into a bgzip-compressed FASTQ.
    Bgzf(bgzf::VirtualPosition)
}


/// A seekable FASTQ reader over either an uncompressed or a bgzip-compressed file.
/// Give each worker thread its own reader so seeks don't contend on a shared cursor.
pub enum FastqReader {
    Plain(fastq::io::Reader<BufReader<File>>),
    Bgzf(fastq::io::Reader<bgzf::io::Reader<File>>)
}

impl FastqReader {
    /// Open `fastq_file`, auto-detecting BGZF by its gzip magic plus the BGZF `BC`
    /// extra subfield.
    ///
    /// # Panics
    /// Panics if the file cannot be opened, or if it is *plain* gzip — plain gzip is a
    /// single deflate stream and is not randomly seekable, so re-compress it with
    /// `bgzip` or leave it uncompressed.
    pub fn open(fastq_file: &str) -> Self {
        if is_bgzipped(fastq_file) {
            let reader = File::open(fastq_file)
                .map(bgzf::io::Reader::new)
                .map(fastq::io::Reader::new)
                .unwrap();
            FastqReader::Bgzf(reader)
        } else if is_gzipped(fastq_file) {
            panic!(
                "'{}' is plain gzip, which is not randomly seekable; \
                 re-compress it with bgzip or decompress it first",
                fastq_file
            );
        } else {
            let reader = File::open(fastq_file)
                .map(BufReader::new)
                .map(fastq::io::Reader::new)
                .unwrap();
            FastqReader::Plain(reader)
        }
    }

    /// Read one read's sequence and quality string by seeking to the `offset` that
    /// [`build_fastq_offset_index`] recorded for `read_name`.
    ///
    /// # Panics
    /// Panics on I/O error, if `offset`'s kind does not match the reader's compression,
    /// if the record at `offset` is not `read_name` (a stale index), or if the sequence
    /// or quality is not valid UTF-8.
    pub fn get_read_record(&mut self, read_name: &str, offset: FastqOffset) -> (Box<str>, Box<str>) {
        let mut record = fastq::Record::default();
        match self {
            FastqReader::Plain(reader) => {
                let FastqOffset::Plain(byte_offset) = offset else {
                    panic!("uncompressed FASTQ needs a byte offset for read '{}'", read_name);
                };
                reader.get_mut().seek(SeekFrom::Start(byte_offset)).unwrap();
                reader.read_record(&mut record).unwrap();
            }
            FastqReader::Bgzf(reader) => {
                let FastqOffset::Bgzf(virtual_position) = offset else {
                    panic!("bgzip FASTQ needs a virtual position for read '{}'", read_name);
                };
                reader.get_mut().seek(virtual_position).unwrap();
                reader.read_record(&mut record).unwrap();
            }
        }

        // The offset must land exactly on this read's record.
        assert_eq!(
            std::str::from_utf8(record.name()).unwrap(),
            read_name,
            "offset does not point at read '{}'", read_name
        );

        (
            std::str::from_utf8(record.sequence()).expect("read sequence is valid UTF-8").into(),
            std::str::from_utf8(record.quality_scores()).expect("read quality is valid UTF-8").into()
        )
    }

    /// Read one read's sequence; see [`FastqReader::get_read_record`].
    pub fn get_read_sequence(&mut self, read_name: &str, offset: FastqOffset) -> Box<str> {
        self.get_read_record(read_name, offset).0
    }
}


/// Record the start position of *every* read in `fastq_file` in one streaming pass, keyed
/// by read name. Offsets are byte offsets for an uncompressed FASTQ and BGZF virtual
/// positions for a bgzip FASTQ, so they pair with a [`FastqReader`] opened on the same
/// file.
///
/// # Panics
/// Panics if the file cannot be opened, or if it is plain (non-BGZF) gzip.
pub fn build_fastq_offset_index(fastq_file: &str) -> HashMap<Box<str>, FastqOffset> {
    if is_bgzipped(fastq_file) {
        build_bgzf_offset_index(fastq_file)
    } else if is_gzipped(fastq_file) {
        panic!(
            "'{}' is plain gzip, which is not randomly seekable; \
             re-compress it with bgzip or decompress it first",
            fastq_file
        );
    } else {
        build_plain_offset_index(fastq_file)
    }
}


/// Index an uncompressed FASTQ by the raw byte offset of each read's '@' line.
fn build_plain_offset_index(fastq_file: &str) -> HashMap<Box<str>, FastqOffset> {
    let mut index = HashMap::new();
    let mut reader = BufReader::new(File::open(fastq_file).unwrap());
    let mut offset: u64 = 0;
    let mut line: Vec<u8> = Vec::new();
    loop {
        let record_start = offset;                     // byte position of the '@' line
        line.clear();
        let n = reader.read_until(b'\n', &mut line).unwrap();
        if n == 0 { break; }                                // EOF
        offset += n as u64;
        if let Some(name) = parse_fastq_name(&line) {
            index.insert(name.into(), FastqOffset::Plain(record_start));
        }
        for _ in 0..3 {                                     // consume seq, '+', qual
            line.clear();
            offset += reader.read_until(b'\n', &mut line).unwrap() as u64;
        }
    }
    index
}


/// Index a bgzip FASTQ by the BGZF virtual position of each read's '@' line.
/// The BGZF reader tracks its uncompressed position as it consumes blocks, so the
/// virtual position captured before a record can be seeked back to exactly.
fn build_bgzf_offset_index(fastq_file: &str) -> HashMap<Box<str>, FastqOffset> {
    let mut index = HashMap::new();
    let mut reader = File::open(fastq_file).map(bgzf::io::Reader::new).unwrap();
    let mut line: Vec<u8> = Vec::new();
    loop {
        let record_start = reader.virtual_position();       // BGZF position of the '@' line
        line.clear();
        let n = reader.read_until(b'\n', &mut line).unwrap();
        if n == 0 { break; }                                // EOF
        if let Some(name) = parse_fastq_name(&line) {
            index.insert(name.into(), FastqOffset::Bgzf(record_start));
        }
        for _ in 0..3 {                                     // consume seq, '+', qual
            line.clear();
            reader.read_until(b'\n', &mut line).unwrap();
        }
    }
    index
}


/// Open `fastq_file` for a single streaming pass, transparently handling uncompressed, plain
/// gzip and BGZF.
///
/// The counterpart to [`FastqReader`], which seeks and therefore rejects plain gzip. Nothing
/// here seeks, so plain gzip is fine — reach for this whenever every record is going to be
/// visited anyway, and for [`FastqReader`] only when a few reads are wanted out of many.
///
/// # Panics
/// Panics if the file cannot be opened.
pub fn open_fastq_reader(fastq_file: &str) -> fastq::io::Reader<Box<dyn BufRead>> {
    let file: File = File::open(fastq_file).expect("Failed to open FASTQ file");
    let reader: Box<dyn BufRead> = if is_bgzipped(fastq_file) {
        Box::new(BufReader::new(bgzf::io::Reader::new(file)))
    } else if is_gzipped(fastq_file) {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    fastq::io::Reader::new(reader)
}


/// Write `records` — (read name, sequence, quality) — as a BGZF-compressed FASTQ.
///
/// BGZF rather than plain gzip is deliberate. [`build_fastq_offset_index`] and [`FastqReader`]
/// index a FASTQ by seeking into it, and plain gzip is a single deflate stream that cannot be
/// seeked — so a `.fastq.gz` written with plain gzip reads fine in most tools and then panics
/// the moment a downstream exacto stage tries to index it. Writing BGZF keeps the output usable
/// as the next stage's input.
///
/// # Panics
/// Panics if the file cannot be created, on I/O error, or if any record's sequence and quality
/// differ in length — which is what a correction step that edits bases without editing the
/// quality string in step would produce.
pub fn write_fastq_file(records: &[(Box<str>, Box<str>, Box<str>)], output_fastq_file: &str) {
    let mut writer = open_fastq_writer(output_fastq_file);
    for (read_name, sequence, quality) in records.iter() {
        write_fastq_record(&mut writer, read_name, sequence.as_bytes(), quality.as_bytes());
    }
    writer.finish().unwrap();
}


/// Open `output_fastq_file` for writing as BGZF — see [`write_fastq_file`] for why BGZF. Call
/// `finish` on the writer when done.
///
/// # Panics
/// Panics if the file cannot be created.
pub fn open_fastq_writer(output_fastq_file: &str) -> bgzf::io::Writer<File> {
    bgzf::io::Writer::new(File::create(output_fastq_file).expect("Failed to create FASTQ file"))
}


/// Write one FASTQ record.
///
/// # Panics
/// Panics on I/O error, or if `sequence` and `quality` differ in length — which is what a
/// correction step that edits bases without editing the quality string in step would produce.
pub fn write_fastq_record<W: Write>(writer: &mut W, read_name: &str, sequence: &[u8], quality: &[u8]) {
    assert_eq!(
        sequence.len(),
        quality.len(),
        "read '{}' has {} base(s) but {} quality value(s).",
        read_name,
        sequence.len(),
        quality.len()
    );
    writeln!(writer, "@{}", read_name).unwrap();
    writer.write_all(sequence).unwrap();
    writer.write_all(b"\n+\n").unwrap();
    writer.write_all(quality).unwrap();
    writer.write_all(b"\n").unwrap();
}


/// Parse the read name from a FASTQ header line: the bytes after '@', up to the first
/// ASCII whitespace. Returns `None` if the line does not start with '@' or the name is
/// not valid UTF-8.
fn parse_fastq_name(header_line: &[u8]) -> Option<&str> {
    let after_at = header_line.strip_prefix(b"@")?;
    let end = after_at
        .iter()
        .position(|&byte| byte == b' ' || byte == b'\t' || byte == b'\n' || byte == b'\r')
        .unwrap_or(after_at.len());
    std::str::from_utf8(&after_at[..end]).ok()
}

