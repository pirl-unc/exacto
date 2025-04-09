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


use flate2::Compression;
use flate2::write::GzEncoder;
use polars::prelude::*;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use serde::{Serialize, Deserialize};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Mutex;

use crate::structs::mutant_peptide::MutantPeptide;


#[derive(Debug,Serialize,Deserialize)]
pub struct MutantPeptidesSet {
    pub mutant_peptides: Vec<MutantPeptide>
}

impl MutantPeptidesSet {
    pub fn new() -> Self {
        Self {
            mutant_peptides: Vec::new()
        }
    }

    pub fn add_mutant_peptide(&mut self, mutant_peptide: MutantPeptide) {
        self.mutant_peptides.push(mutant_peptide);
    }

    pub fn get_size(&self) -> usize {
        self.mutant_peptides.len()
    }

    pub fn to_dataframe(self, num_threads: usize) -> DataFrame {
        let thread_pool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();
        let chunk_size = (self.mutant_peptides.len() + num_threads - 1) / num_threads;
        let rows = thread_pool.install(|| {
            self.mutant_peptides
                .par_chunks(chunk_size)
                .flat_map_iter(|chunk| {
                    chunk.iter().map(|mutant_peptide| {
                        (
                            mutant_peptide.id.as_ref(),
                            mutant_peptide.peptide_sequence.as_ref(),
                            mutant_peptide.get_peptide_sequence_read_names_string(),
                            mutant_peptide.k as u32,
                            mutant_peptide.get_kmer_mask_string(),
                            mutant_peptide.get_rna_variants_string(),
                            mutant_peptide.get_rna_variants_read_names_string(),
                            mutant_peptide.get_dna_variants_string(),
                            mutant_peptide.get_dna_variants_read_names_string()
                        )
                    })
                })
                .collect::<Vec<_>>()
        });

        let peptide_id_column = rows.iter().map(|r| r.0).collect::<Vec<_>>();
        let peptide_sequence_column = rows.iter().map(|r| r.1).collect::<Vec<_>>();
        let peptide_sequence_read_names_column = rows.iter().map(|r| r.2.clone()).collect::<Vec<_>>();
        let k_column = rows.iter().map(|r| r.3).collect::<Vec<_>>();
        let kmer_mask_column = rows.iter().map(|r| r.4.clone()).collect::<Vec<_>>();
        let rna_variants_column = rows.iter().map(|r| r.5.clone()).collect::<Vec<_>>();
        let rna_variants_read_names_column = rows.iter().map(|r| r.6.clone()).collect::<Vec<_>>();
        let dna_variants_column = rows.iter().map(|r| r.7.clone()).collect::<Vec<_>>();
        let dna_variants_read_names_column = rows.iter().map(|r| r.8.clone()).collect::<Vec<_>>();

        // Build the DataFrame
        DataFrame::new(vec![
            Column::from(Series::new("peptide_id".into(), peptide_id_column)),
            Column::from(Series::new("peptide_sequence".into(), peptide_sequence_column)),
            Column::from(Series::new("peptide_sequence_read_names".into(), peptide_sequence_read_names_column)),
            Column::from(Series::new("k".into(), k_column)),
            Column::from(Series::new("kmer_mask".into(), kmer_mask_column)),
            Column::from(Series::new("rna_variants".into(), rna_variants_column)),
            Column::from(Series::new("rna_variants_read_names".into(), rna_variants_read_names_column)),
            Column::from(Series::new("dna_variants".into(), dna_variants_column)),
            Column::from(Series::new("dna_variants_read_names".into(), dna_variants_read_names_column))
        ]).unwrap()
    }

    pub fn to_tsv(&self, file: &str, buffer_size: usize, num_threads: usize, gzip: bool) {
        let file = File::create(file).unwrap();

        let thread_pool = ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();

        let header: String = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            "peptide_id",
            "peptide_sequence",
            "peptide_read_names",
            "k",
            "kmer_mask",
            "rna_variants",
            "rna_variants_read_names",
            "dna_variants",
            "dna_variants_read_names"
        );

        let mutant_peptides: Vec<&MutantPeptide> = self.mutant_peptides.iter().collect();

        if gzip {
            let buf_writer = BufWriter::new(file);
            let writer = Arc::new(Mutex::new(GzEncoder::new(buf_writer, Compression::default())));
            writer.lock().unwrap().write_all(header.as_bytes()).unwrap();
            thread_pool.scope(|s| {
                for (chunk_idx, chunk) in mutant_peptides.chunks(buffer_size).enumerate() {
                    let writer = writer.clone();
                    s.spawn(move |_| {
                        let mut local_buffer = String::new();
                        for (i, mutant_peptide) in chunk.iter().enumerate() {
                            let row = mutant_peptide.to_tsv_string();
                            local_buffer.push_str(&row);
                        }
                        let mut writer_guard = writer.lock().unwrap();
                        writer_guard.write_all(local_buffer.as_bytes()).unwrap();
                    });
                }
            });
            writer.lock().unwrap().flush().unwrap();
        } else {
            let writer = Arc::new(Mutex::new(BufWriter::new(file)));
            writer.lock().unwrap().write_all(header.as_bytes()).unwrap();
            thread_pool.scope(|s| {
                for (chunk_idx, chunk) in mutant_peptides.chunks(buffer_size).enumerate() {
                    let writer = writer.clone();
                    s.spawn(move |_| {
                        let mut local_buffer = String::new();
                        for (i, mutant_peptide) in chunk.iter().enumerate() {
                            let row = mutant_peptide.to_tsv_string();
                            local_buffer.push_str(&row);
                        }
                        let mut writer_guard = writer.lock().unwrap();
                        writer_guard.write_all(local_buffer.as_bytes()).unwrap();
                    });
                }
            });
            writer.lock().unwrap().flush().unwrap();
        }
    }
}
