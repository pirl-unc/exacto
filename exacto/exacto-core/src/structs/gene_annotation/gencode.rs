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
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufRead,BufReader};
use std::str::FromStr;

use crate::prelude::*;


#[derive(Debug)]
pub struct Gencode {
    pub assembly: Box<str>,
    pub version: Box<str>,
    pub genes: HashMap<Box<str>,Gene>,                                      // key = gene ID
    pub transcript_gene_ids: HashMap<Box<str>,Box<str>>,                    // key = transcript ID, value = gene ID
    pub exon_transcript_ids: HashMap<Box<str>,Box<str>>,                    // key = exon ID, value = transcript ID
    pub gene_itrees_map: HashMap<Box<str>,IntervalTree<Box<str>>>,          // key = chromosome, interval tree has gene IDs
}

impl Gencode {
    /// Constructor for Gencode object.
    ///
    /// # Arguments
    /// * `gtf_file`: GTF file.
    /// * `assembly`: Assembly (e.g. 'human').
    /// * `version`: GENCODE version.
    /// * `gene_types`: Gene types (optional). If None, all gene types are loaded.
    /// * `gene_levels`: Gene levels (optional). If None, all gene levels are loaded.
    /// * `transcript_types`: Transcript types (optional). If None, all transcript types are loaded.
    /// * `transcript_levels`: Transcript levels (optional). If None, all transcript levels are loaded.
    pub fn new(
        gtf_file: &str,
        assembly: &str,
        version: &str,
        gene_types: Option<HashSet<&str>>,
        gene_levels: Option<HashSet<u8>>,
        transcript_types: Option<HashSet<&str>>,
        transcript_levels: Option<HashSet<u8>>
    ) -> Self {
        fn extract_attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
            let key_prefix = format!("{key} ");
            for field in attrs.split(';') {
                let field = field.trim();
                if field.starts_with(&key_prefix) {
                    let value = field[key_prefix.len()..].trim();
                    return Some(value.trim_matches('"'));
                }
            }
            None
        }
        // Step 1. Prepare data structures
        let mut genes: HashMap<Box<str>,Gene> = HashMap::new();
        let mut transcript_gene_ids: HashMap<Box<str>,Box<str>> = HashMap::new();
        let mut exon_transcript_ids: HashMap<Box<str>,Box<str>> = HashMap::new();

        // Step 2. Load genes
        let file = File::open(gtf_file).expect("Unable to reopen GTF file");
        let reader: Box<dyn BufRead> = if is_gzipped(gtf_file) {
            let decoder = MultiGzDecoder::new(file);
            Box::new(BufReader::new(decoder))
        } else {
            Box::new(BufReader::new(file))
        };
        for line in reader.lines().flatten() {
            if line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() != 9 {
                continue;
            }
            let feature_type = fields[2];
            if feature_type != "gene" {
                continue;
            }

            let seq_name: &str = fields[0];
            let source: &str = fields[1];
            let start: u32 = fields[3].parse::<u32>().unwrap_or(0);
            let end: u32 = fields[4].parse::<u32>().unwrap_or(0);
            let strand: Strand = Strand::from_str(fields[6]).unwrap();

            let attrs: &str = fields[8];
            let gene_id: &str = extract_attr(attrs, "gene_id").unwrap_or("");
            let gene_name: &str = extract_attr(attrs, "gene_name").unwrap_or("");
            let gene_type: &str = extract_attr(attrs, "gene_type").unwrap_or("");
            let level: u8 = extract_attr(attrs, "level").and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);

            if gene_types.is_some() && !gene_types.as_ref().unwrap().contains(&gene_type) {
                continue;
            }

            if gene_levels.is_some() && !gene_levels.as_ref().unwrap().contains(&level) {
                continue;
            }

            let gene: Gene = Gene::new(
                gene_id,
                source,
                seq_name,
                start,
                end,
                strand,
                gene_name,
                level,
                gene_type
            );

            genes.insert(gene_id.into(), gene);
        }

        // Step 3. Load transcripts
        let file = File::open(gtf_file).expect("Unable to reopen GTF file");
        let reader: Box<dyn BufRead> = if is_gzipped(gtf_file) {
            let decoder = MultiGzDecoder::new(file);
            Box::new(BufReader::new(decoder))
        } else {
            Box::new(BufReader::new(file))
        };
        for line in reader.lines().flatten() {
            if line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() != 9 {
                continue;
            }
            let feature_type = fields[2];
            if feature_type != "transcript" {
                continue;
            }

            let seq_name: &str = fields[0];
            let source: &str = fields[1];
            let start: u32 = fields[3].parse::<u32>().unwrap_or(0);
            let end: u32 = fields[4].parse::<u32>().unwrap_or(0);
            let strand: Strand = Strand::from_str(fields[6]).unwrap();

            let attrs: &str = fields[8];
            let gene_id: &str = extract_attr(attrs, "gene_id").unwrap_or("");
            let transcript_id: &str = extract_attr(attrs, "transcript_id").unwrap_or("");
            let transcript_name: &str = extract_attr(attrs, "transcript_name").unwrap_or("");
            let transcript_type: &str = extract_attr(attrs, "transcript_type").unwrap_or("");
            let transcript_support_level: &str = extract_attr(attrs, "transcript_support_level").unwrap_or("");
            let level: u8 = extract_attr(attrs, "level").and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);

            if transcript_types.is_some() && !transcript_types.as_ref().unwrap().contains(&transcript_type) {
                continue;
            }

            if transcript_levels.is_some() && !transcript_levels.as_ref().unwrap().contains(&level) {
                continue;
            }

            let transcript: Transcript = Transcript::new(
                gene_id,
                transcript_id,
                source,
                seq_name,
                start,
                end,
                strand,
                level,
                transcript_name,
                transcript_type,
                transcript_support_level,
            );

            if genes.contains_key(gene_id) {
                genes
                    .get_mut(gene_id).unwrap()
                    .add_transcript(transcript);
                transcript_gene_ids.insert(transcript_id.into(), gene_id.into());
            }
        }

        // Step 4. Load exons, UTRs start codon, stop codon, CDS
        let file = File::open(gtf_file).expect("Unable to reopen GTF file");
        let reader: Box<dyn BufRead> = if is_gzipped(gtf_file) {
            let decoder = MultiGzDecoder::new(file);
            Box::new(BufReader::new(decoder))
        } else {
            Box::new(BufReader::new(file))
        };
        for line in reader.lines().flatten() {
            if line.starts_with('#') {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() != 9 {
                continue;
            }
            let feature_type = fields[2];
            if feature_type == "gene" {
                continue;
            }
            if feature_type == "transcript" {
                continue;
            }

            let seq_name: &str = fields[0];
            let source: &str = fields[1];
            let start: u32 = fields[3].parse::<u32>().unwrap_or(0);
            let end: u32 = fields[4].parse::<u32>().unwrap_or(0);
            let strand: Strand = Strand::from_str(fields[6]).unwrap();

            let attrs: &str = fields[8];
            let gene_id: &str = extract_attr(attrs, "gene_id").unwrap_or("");
            let transcript_id: &str = extract_attr(attrs, "transcript_id").unwrap_or("");
            let exon_id: &str = extract_attr(attrs, "exon_id").unwrap_or("");
            let level: u8 = extract_attr(attrs, "level").and_then(|s| s.parse::<u8>().ok()).unwrap_or(0);
            let exon_number = extract_attr(attrs, "exon_number").and_then(|s| s.parse::<u16>().ok()).unwrap_or(0);

            if genes.contains_key(gene_id) == false {
                continue;
            }

            if genes.get(gene_id).unwrap().transcripts.contains_key(transcript_id) == false {
                continue;
            }

            match feature_type {
                "exon" => {
                    let exon = Exon::new(
                        gene_id,
                        transcript_id,
                        exon_id,
                        source,
                        seq_name,
                        start,
                        end,
                        strand,
                        level,
                        exon_number
                    );
                    genes.get_mut(gene_id).unwrap()
                        .get_mut_transcript(transcript_id).unwrap()
                        .add_exon(exon);
                    exon_transcript_ids.insert(exon_id.into(), transcript_id.into());
                },
                "UTR" => {
                    let utr = UTR::new(
                        gene_id,
                        transcript_id,
                        exon_id,
                        source,
                        seq_name,
                        start,
                        end,
                        strand,
                        level,
                        exon_number
                    );
                    genes.get_mut(gene_id).unwrap()
                        .get_mut_transcript(transcript_id).unwrap()
                        .add_utr(utr);
                },
                "start_codon" => {
                    let start_codon = StartCodon::new(
                        gene_id,
                        transcript_id,
                        exon_id,
                        source,
                        seq_name,
                        start,
                        end,
                        strand,
                        level,
                        exon_number
                    );
                    genes.get_mut(gene_id).unwrap()
                        .get_mut_transcript(transcript_id).unwrap()
                        .add_start_codon(start_codon);
                },
                "stop_codon" => {
                    let stop_codon = StopCodon::new(
                        gene_id,
                        transcript_id,
                        exon_id,
                        source,
                        seq_name,
                        start,
                        end,
                        strand,
                        level,
                        exon_number
                    );
                    genes.get_mut(gene_id).unwrap()
                        .get_mut_transcript(transcript_id).unwrap()
                        .add_stop_codon(stop_codon);
                },
                "CDS" => {
                    let cds = CDS::new(
                        gene_id,
                        transcript_id,
                        exon_id,
                        source,
                        seq_name,
                        start,
                        end,
                        strand,
                        level,
                        exon_number
                    );
                    genes.get_mut(gene_id).unwrap()
                        .get_mut_transcript(transcript_id).unwrap()
                        .add_cds(cds);
                },
                _ => {}
            }
        }

        // Step 5. Build interval trees
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

    pub fn new_with_defaults(gtf_file: &str, assembly: &str, version: &str) -> Self {
        Self::new(
            gtf_file,
            assembly,
            version,
            None,
            None,
            None,
            None
        )
    }
}


impl GeneAnnotator for Gencode {
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
        let gene_id = self.transcript_gene_ids.get(transcript_id).unwrap();
        self.genes.get(gene_id).unwrap().transcripts.get(transcript_id)
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

    fn get_exon(&self, exon_id: &str) -> Option<&Exon> {
        let transcript_id = self.exon_transcript_ids.get(exon_id).unwrap();
        let gene_id = self.transcript_gene_ids.get(transcript_id).unwrap();
        self.genes.get(gene_id).unwrap().transcripts.get(transcript_id).unwrap().exons.get(exon_id)
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
}
