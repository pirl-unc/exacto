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


use std::collections::{HashMap, HashSet};
use std::str::FromStr;
use polars::prelude::DataFrame;

use crate::prelude::*;


#[derive(Debug)]
pub struct TsvGeneAnnotator {
    assembly: Box<str>,
    version: Box<str>,
    genes: HashMap<Box<str>, Gene>,                                 // key = gene ID
    transcript_gene_ids: HashMap<Box<str>, Box<str>>,               // key = transcript ID, value = gene ID
    exon_transcript_ids: HashMap<Box<str>, HashSet<Box<str>>>,      // key = exon ID, value = transcript IDs
    gene_itrees_map: HashMap<Box<str>, IntervalTree<Box<str>>>      // key = chromosome, interval tree has gene IDs
}

impl TsvGeneAnnotator {
    /// Constructor for TsvGeneAnnotator object.
    ///
    /// # Arguments
    /// * `tsv_file`: TSV file.
    /// * `assembly`: Assembly (e.g. 'human').
    /// * `version`: Gene annotation version.
    pub fn new(
        tsv_file: &str,
        assembly: &str,
        version: &str
    ) -> Self {
        let df: DataFrame = read_tsv_file(tsv_file);

        // Step 1. Prepare data structures
        let mut genes: HashMap<Box<str>,Gene> = HashMap::new();
        let mut transcript_gene_ids: HashMap<Box<str> ,Box<str>> = HashMap::new();
        let mut exon_transcript_ids: HashMap<Box<str>, HashSet<Box<str>>> = HashMap::new();

        // Step 2. Get DataFrame columns
        let col_row_type = df.column("row_type").unwrap().str().unwrap();
        let col_gene_id = df.column("gene_id").unwrap().str().unwrap();
        let col_transcript_id = df.column("transcript_id").unwrap().str().unwrap();
        let col_exon_id = df.column("exon_id").unwrap().str().unwrap();
        let col_exon_number = df.column("exon_number").unwrap().i64().unwrap();
        let col_strand = df.column("strand").unwrap().str().unwrap();
        let col_chromosome = df.column("chromosome").unwrap().str().unwrap();
        let col_start = df.column("start").unwrap().i64().unwrap();
        let col_end = df.column("end").unwrap().i64().unwrap();

        // Step 3. Load genes
        for i in 0..df.height() {
            let row_type: &str = col_row_type.get(i).unwrap_or("");
            let gene_id: &str = col_gene_id.get(i).unwrap_or("");
            let strand: Strand = Strand::from_str(col_strand.get(i).unwrap()).unwrap();
            let chromosome: &str = col_chromosome.get(i).unwrap();
            let start: u32 = col_start.get(i).unwrap() as u32;
            let end: u32 = col_end.get(i).unwrap() as u32;

            assert!(row_type == "gene" || row_type == "transcript" || row_type == "exon");

            if row_type == "gene" {
                let gene = Gene::new(
                    gene_id,
                    "",
                    chromosome,
                    start,
                    end,
                    strand,
                    "",
                    0,
                    ""
                );
                genes.insert(gene_id.into(), gene);
            }
        }

        // Step 4. Load transcripts
        for i in 0..df.height() {
            let row_type: &str = col_row_type.get(i).unwrap_or("");
            let gene_id: &str = col_gene_id.get(i).unwrap_or("");
            let transcript_id: &str = col_transcript_id.get(i).unwrap_or("");
            let strand: Strand = Strand::from_str(col_strand.get(i).unwrap()).unwrap();
            let chromosome: &str = col_chromosome.get(i).unwrap();
            let start: u32 = col_start.get(i).unwrap() as u32;
            let end: u32 = col_end.get(i).unwrap() as u32;

            assert!(row_type == "gene" || row_type == "transcript" || row_type == "exon");

            if row_type == "transcript" {
                let transcript = Transcript::new(
                    gene_id,
                    transcript_id,
                    "",
                    chromosome,
                    start,
                    end,
                    strand,
                    0,
                    "",
                    "",
                    "",
                    HashSet::new()
                );
                genes.get_mut(gene_id).unwrap().add_transcript(transcript);
                transcript_gene_ids.insert(transcript_id.into(), gene_id.into());
            }
        }

        // Step 5. Load exons
        for i in 0..df.height() {
            let row_type: &str = col_row_type.get(i).unwrap_or("");
            let gene_id: &str = col_gene_id.get(i).unwrap_or("");
            let transcript_id: &str = col_transcript_id.get(i).unwrap_or("");
            let exon_id: &str = col_exon_id.get(i).unwrap_or("");
            let exon_number: u16 = col_exon_number.get(i).unwrap_or(0) as u16;
            let strand: Strand = Strand::from_str(col_strand.get(i).unwrap()).unwrap();
            let chromosome: &str = col_chromosome.get(i).unwrap();
            let start: u32 = col_start.get(i).unwrap() as u32;
            let end: u32 = col_end.get(i).unwrap() as u32;

            assert!(row_type == "gene" || row_type == "transcript" || row_type == "exon");

            if row_type == "exon" {
                let exon: Exon = Exon::new(
                    gene_id,
                    transcript_id,
                    exon_id,
                    "",
                    chromosome,
                    start,
                    end,
                    strand,
                    0,
                    exon_number
                );
                genes.get_mut(gene_id).unwrap().get_mut_transcript(transcript_id).unwrap().add_exon(exon);
                exon_transcript_ids
                    .entry(exon_id.into())
                    .or_insert_with(HashSet::new)
                    .insert(transcript_id.into());
            }
        }

        // Step 6. Build interval trees
        let mut gene_itrees_map: HashMap<Box<str>,IntervalTree<Box<str>>> = HashMap::new();
        for gene in genes.values() {
            let interval = Interval::new(
                gene.start as isize,
                gene.end as isize,
                gene.gene_id.clone()
            );
            gene_itrees_map
                .entry(gene.chromosome.clone())
                .or_insert_with(IntervalTree::new)
                .insert(interval);
        }

        Self {
            assembly: assembly.to_string().into(),
            version: version.to_string().into(),
            genes,
            transcript_gene_ids,
            exon_transcript_ids,
            gene_itrees_map
        }
    }
}


impl GeneAnnotator for TsvGeneAnnotator {
    fn get_assembly(&self) -> &str {
        &*self.assembly
    }

    fn get_version(&self) -> &str {
        &*self.version
    }

    fn get_gene_ids_at_locus(&self, chromosome: &str, position: u32) -> Vec<Box<str>> {
        if let Some(itree) = self.gene_itrees_map.get(chromosome) {
            let results: Vec<&Box<str>> = itree.overlaps(position as isize, position as isize);
            if results.is_empty() {
                Vec::new()
            } else {
                results.iter().map(|b| (**b).clone()).collect()
            }
        } else {
            Vec::new()
        }
    }

    fn get_gene_ids_overlapping_region(
        &self,
        chromosome: &str,
        start: u32,
        end: u32
    ) -> Vec<Box<str>> {
        if let Some(itree) = self.gene_itrees_map.get(chromosome) {
            let results: Vec<&Box<str>> = itree.overlaps(start as isize, end as isize);
            if results.is_empty() {
                Vec::new()
            } else {
                results.iter().map(|b| (**b).clone()).collect()
            }
        } else {
            Vec::new()
        }
    }

    fn get_transcript_ids_overlapping_region(&self, chromosome: &str, start: u32, end: u32) -> Vec<Box<str>> {
        let mut transcript_ids: Vec<Box<str>> = Vec::new();
        let gene_ids: Vec<Box<str>> = self.get_gene_ids_overlapping_region(chromosome, start, end);
        let start_: isize = start as isize;
        let end_: isize = end as isize;;
        for gene_id in gene_ids.iter() {
            let gene: &Gene = self.get_gene(&**gene_id).unwrap();
            for transcript in gene.transcripts.values() {
                let transcript_start: isize = transcript.start as isize;
                let transcript_end: isize = transcript.end as isize;
                if overlaps(start_, end_, transcript_start, transcript_end) {
                    transcript_ids.push(transcript.transcript_id.clone());
                }
            }
        }
        transcript_ids
    }

    fn get_exon_ids_overlapping_region(&self, chromosome: &str, start: u32, end: u32) -> Vec<Box<str>> {
        let mut exon_ids: Vec<Box<str>> = Vec::new();
        let transcript_ids: Vec<Box<str>> = self.get_transcript_ids_overlapping_region(chromosome, start, end);
        let start_: isize = start as isize;
        let end_: isize = end as isize;;
        for transcript_id in transcript_ids.iter() {
            let transcript: &Transcript = self.get_transcript(transcript_id).unwrap();
            for exon in transcript.exons.values() {
                let exon_start: isize = exon.start as isize;
                let exon_end: isize = exon.end as isize;
                if overlaps(start_, end_, exon_start, exon_end) {
                    exon_ids.push(exon.exon_id.clone());
                }
            }
        }
        exon_ids
    }

    fn get_gene(&self, gene_id: &str) -> Option<&Gene> {
        if self.genes.contains_key(gene_id) {
            Some(self.genes.get(gene_id).unwrap())
        } else {
            None
        }
    }

    fn get_genes(&self) -> Vec<&Gene> {
        self.genes.values().collect()
    }

    fn get_transcript(&self, transcript_id: &str) -> Option<&Transcript> {
        let gene_id = self.transcript_gene_ids.get(transcript_id)?;
        self.genes.get(gene_id)?.transcripts.get(transcript_id)
    }

    fn get_transcripts(&self) -> Vec<&Transcript> {
        let mut transcripts: Vec<&Transcript> = Vec::new();
        for gene in self.genes.values() {
            for transcript in gene.transcripts.values() {
                transcripts.push(transcript);
            }
        }
        transcripts
    }

    fn get_exon(&self, transcript_id: &str, exon_id: &str) -> Option<&Exon> {
        let gene_id = self.transcript_gene_ids.get(transcript_id)?;
        self.genes.get(gene_id)?.transcripts.get(transcript_id)?.exons.get(exon_id)
    }
    
    fn get_exons(&self) -> Vec<&Exon> {
        let mut exons: Vec<&Exon> = Vec::new();
        for gene in self.genes.values() {
            for transcript in gene.transcripts.values() {
                for exon in transcript.exons.values() {
                    exons.push(exon);
                }
            }
        }
        exons
    }

    fn rank_transcripts<'a>(&self, transcripts: Vec<&'a Transcript>) -> Vec<&'a Transcript> {
        transcripts // no ranking information available, return as-is
    }
}
