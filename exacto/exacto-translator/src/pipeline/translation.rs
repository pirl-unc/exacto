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


use exacto_caller::prelude::{
    AssembledTranscriptModelAlignmentRecord,
    AssembledTranscriptVariantRecord,
    DNAVariantRecord
};
use exacto_integrator::prelude::IntegratedVariantRecord;
use noodles_bgzf as bgzf;
use seq_io::fasta::{Reader as FastaReader, Record as FastaRecord};
use seq_io::fastq::{Reader as FastqReader, Record as FastqRecord};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::num::NonZeroUsize;

use crate::prelude::*;


pub fn translate_transcripts(
    assembled_transcript_support_records: &Vec<AssembledTranscriptSupportRecord>,
    assembled_transcript_model_alignment_records: &Vec<AssembledTranscriptModelAlignmentRecord>,
    assembled_transcript_variant_records: &Vec<AssembledTranscriptVariantRecord>,
    dna_variant_records: &Vec<DNAVariantRecord>,
    integrated_variant_records: &Vec<IntegratedVariantRecord>,
    reference_stitched_spans: &Vec<ReferenceStitchedSpans>,
    translation_strategy: TranslationStrategy,
    start_codons: &HashSet<&str>,
    num_threads: usize
) -> Result<AssembledTranscriptSet, TranslatorError> {
    let mut ts: AssembledTranscriptSet = build_transcript_set(
        assembled_transcript_support_records,
        assembled_transcript_model_alignment_records,
        assembled_transcript_variant_records,
        dna_variant_records,
        integrated_variant_records,
        reference_stitched_spans
    )?;
    
    ts.translate(
        translation_strategy,
        start_codons,
        num_threads
    )?;
    
    Ok(ts)
}


pub fn translate_sequences(
        sequences: Vec<(Box<str>, Box<str>)>,
        translation_strategy: TranslationStrategy,
        start_codons: &HashSet<&str>,
        num_threads: usize
) -> Result<AssembledTranscriptSet, TranslatorError> {
    let mut transcripts: Vec<AssembledTranscript> = Vec::with_capacity(sequences.len());
    for (transcript_id, sequence) in sequences.into_iter() {
        // Codons are read by byte.
        if !sequence.is_ascii() {
            return Err(TranslatorError::NonAsciiSequence { transcript: transcript_id });
        }
        let transcript: AssembledTranscript = AssembledTranscript::new_with_default(transcript_id, sequence);
        transcripts.push(transcript);
    }

    let mut ts: AssembledTranscriptSet = AssembledTranscriptSet::new(transcripts);

    ts.translate(
        translation_strategy,
        start_codons,
        num_threads
    )?;

    Ok(ts)
}


fn translate_and_write_record(
    transcript_id: &str,
    transcript_sequence: &str,
    translation_strategy: &TranslationStrategy,
    start_codons: &HashSet<&str>,
    peptide_sequence_id: &mut u32,
    fasta_writer: &mut impl Write,
    tsv_writer: &mut impl Write,
) -> std::io::Result<()> {
    // Step 1. Translate
    let mut transcript: AssembledTranscript =
        AssembledTranscript::new_with_default(transcript_id.into(), transcript_sequence.into());
    transcript.translate(translation_strategy, start_codons);

    // Step 2. Change the primary sequence ID
    for primary_structure in transcript.proteoforms.iter_mut() {
        primary_structure.id = *peptide_sequence_id;
        *peptide_sequence_id += 1;
    }

    // Step 3. Write to the output FASTA and TSV files
    for primary_structure in transcript.proteoforms.iter() {
        let peptide_id: String = primary_structure.id.to_string();
        let peptide_sequence: &str = primary_structure.get_sequence();
        writeln!(
            tsv_writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            transcript_id,
            transcript_sequence,
            peptide_id,
            peptide_sequence,
            peptide_sequence.len() as u32,
            primary_structure.orf_start,
            primary_structure.orf_end
        )?;
        seq_io::fasta::write_to(&mut *fasta_writer, peptide_id.as_bytes(), peptide_sequence.as_bytes())?;
    }
    Ok(())
}


pub fn translate_fastx_file(
    fastx_file: &str,
    output_fasta_file: &str,
    output_tsv_file: &str,
    translation_strategy: TranslationStrategy,
    start_codons: HashSet<&str>,
    num_threads: usize,
) -> Result<(), TranslatorError> {
    let file_error = |file: &str, error: &dyn std::fmt::Display| TranslatorError::File {
        file: file.into(),
        reason: error.to_string().into_boxed_str()
    };
    let record_error = |record: usize, reason: &dyn std::fmt::Display| TranslatorError::Record {
        file: fastx_file.into(),
        record,
        reason: reason.to_string().into_boxed_str()
    };

    // Step 1. Open and sniff the header: BGZF vs regular gzip vs plain.
    let mut input = BufReader::new(File::open(fastx_file).map_err(|e| file_error(fastx_file, &e))?);
    let (is_bgzf, is_gzip) = {
        let head = input.fill_buf().map_err(|e| file_error(fastx_file, &e))?;
        let is_gzip = head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b;
        // BGZF is gzip carrying the "BC" extra subfield (SAM spec §4.1).
        let is_bgzf = head.len() >= 18
            && head[0] == 0x1f && head[1] == 0x8b && head[2] == 0x08
            && (head[3] & 0x04) != 0
            && head[12] == 0x42 && head[13] == 0x43;
        (is_bgzf, is_gzip)
    };

    // Step 2. Decompress: multithreaded BGZF, single-threaded gzip, or plain.
    let decoded: Box<dyn Read> = if is_bgzf {
        let worker_count = NonZeroUsize::new(num_threads.max(1)).unwrap();
        Box::new(bgzf::io::MultithreadedReader::with_worker_count(worker_count, input))
    } else if is_gzip {
        Box::new(flate2::read::MultiGzDecoder::new(input))
    } else {
        Box::new(input)
    };

    // Step 3. Peek the first DECOMPRESSED byte to pick FASTA ('>') vs FASTQ ('@').
    let mut decoded = BufReader::new(decoded);
    let is_fasta = matches!(decoded.fill_buf().map_err(|e| file_error(fastx_file, &e))?.first(), Some(b'>'));

    // Step 4. Output writers
    let mut fasta_writer = File::create(output_fasta_file)
        .map(bgzf::io::Writer::new)
        .map_err(|e| file_error(output_fasta_file, &e))?;
    let mut tsv_writer = BufWriter::new(File::create(output_tsv_file).map_err(|e| file_error(output_tsv_file, &e))?);
    writeln!(tsv_writer, "transcript_id\ttranscript_sequence\tpeptide_id\tpeptide_sequence\tpeptide_length\torf_start\torf_end")
        .map_err(|e| file_error(output_tsv_file, &e))?;

    // Step 5. Parse with the matching reader; translate + write each record. A name or a
    // sequence that is not ASCII text is an error naming the record (1-based).
    let mut peptide_sequence_id: u32 = 1;
    let mut record_number: usize = 0;
    let mut translate_and_write = |id: Result<&str, std::str::Utf8Error>, seq: &[u8], record_number: usize| -> Result<(), TranslatorError> {
        let id: &str = id.map_err(|e| record_error(record_number, &e))?;
        let seq: &str = std::str::from_utf8(seq).map_err(|e| record_error(record_number, &e))?;
        if !seq.is_ascii() {
            return Err(record_error(record_number, &"the sequence holds a character outside ASCII"));
        }
        translate_and_write_record(
            id,
            seq,
            &translation_strategy,
            &start_codons,
            &mut peptide_sequence_id,
            &mut fasta_writer,
            &mut tsv_writer
        ).map_err(|e| file_error(output_tsv_file, &e))
    };
    if is_fasta {
        let mut reader = FastaReader::new(decoded);
        while let Some(rec) = reader.next() {
            record_number += 1;
            let record = rec.map_err(|e| record_error(record_number, &e))?;
            translate_and_write(record.id(), &record.full_seq(), record_number)?; // full_seq strips line breaks
        }
    } else {
        let mut reader = FastqReader::new(decoded);
        while let Some(rec) = reader.next() {
            record_number += 1;
            let record = rec.map_err(|e| record_error(record_number, &e))?;
            translate_and_write(record.id(), record.seq(), record_number)?;
        }
    }

    tsv_writer.flush().map_err(|e| file_error(output_tsv_file, &e))?;
    fasta_writer.try_finish().map_err(|e| file_error(output_fasta_file, &e))?;
    Ok(())
}


#[cfg(test)]
#[path = "../tests/pipeline/translation.rs"]
mod tests;