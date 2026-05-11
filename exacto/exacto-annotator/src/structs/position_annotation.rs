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


use exacto_core::prelude::GenicRegion;
use serde::{Serialize, Deserialize};
use std::collections::{HashMap, HashSet};
use std::hash::Hasher;


#[derive(Debug,Eq,PartialEq,Serialize,Deserialize)]
pub struct PositionAnnotation {
    pub genic_region: GenicRegion,
    pub gene_ids: HashSet<Box<str>>,

    /// Mapping between gene ID and transcript ID
    pub transcript_ids: HashMap<Box<str>, HashSet<Box<str>>>,

    /// Mapping between transcript ID and exon ID
    pub exon_ids: HashMap<Box<str>, Box<str>>
}

impl PositionAnnotation {
    pub fn new(genic_region: GenicRegion) -> Self {
        Self {
            genic_region: genic_region,
            gene_ids: HashSet::new(),
            transcript_ids: HashMap::new(),
            exon_ids: HashMap::new()
        }
    }

    pub fn add_exon_id(&mut self, transcript_id: &str, exon_id: &str) {
        self.exon_ids.insert(transcript_id.into(), exon_id.into());
    }

    pub fn add_gene_id(&mut self, gene_id: &str) {
        self.gene_ids.insert(gene_id.into());
    }
    
    pub fn add_transcript_id(&mut self, gene_id: &str, transcript_id: &str) {
        self.transcript_ids
            .entry(gene_id.into())
            .or_insert_with(HashSet::new)
            .insert(transcript_id.into());
    }

    pub fn from_string(s: &str) -> Self {
        let parts: Vec<&str> = s.split(';').collect();

        // Gene IDs
        let gene_ids: HashSet<Box<str>> = parts.get(0)
            .map(|p| p.split(',').map(|s| s.into()).collect())
            .unwrap_or_default();

        // Transcript IDs
        let mut transcript_ids: HashMap<Box<str>, HashSet<Box<str>>> = HashMap::new();
        if let Some(p) = parts.get(1) {
            for pair in p.split(',') {
                let mut split = pair.splitn(2, '-');
                if let (Some(gene_id), Some(transcript_id)) = (split.next(), split.next()) {
                    transcript_ids
                        .entry(gene_id.into())
                        .or_insert_with(HashSet::new)
                        .insert(transcript_id.into());
                }
            }
        }
        
        // Exon IDs
        let exon_ids: HashMap<Box<str>, Box<str>> = parts.get(2)
            .map(|p| {
                p.split(',')
                    .filter_map(|pair| {
                        let mut split = pair.splitn(2, '-');
                        let transcript_id = split.next()?;
                        let exon_id = split.next()?;
                        Some((transcript_id.into(), exon_id.into()))
                    })
                    .collect()
            })
            .unwrap_or_default();

        // Determine genic region
        let genic_region = if !exon_ids.is_empty() {
            GenicRegion::Exonic
        } else if !transcript_ids.is_empty() {
            GenicRegion::Intronic
        } else {
            GenicRegion::Intergenic
        };

        Self {
            genic_region,
            gene_ids,
            transcript_ids,
            exon_ids
        }
    }

    pub fn to_string(&self) -> String {
        let mut parts: Vec<String> = Vec::new();

        // Gene IDs
        let mut gene_ids: Vec<_> = self.gene_ids.iter().map(|s| s.as_ref()).collect();
        gene_ids.sort();
        parts.push(gene_ids.join(","));

        // Transcript IDs
        let mut transcripts: Vec<String> = self.transcript_ids
            .iter()
            .flat_map(|(gene_id, transcript_ids)| {
                transcript_ids.iter().map(move |tid| format!("{}-{}", gene_id, tid))
            })
            .collect::<Vec<_>>();
        transcripts.sort();
        let transcript_ids: String = transcripts.join(",");
        parts.push(transcript_ids);

        // Exon IDs
        let mut exons: Vec<String> = self.exon_ids
            .iter()
            .map(|(transcript_id, exon_id)| format!("{}-{}", transcript_id, exon_id))
            .collect::<Vec<_>>();
        exons.sort();
        let exon_ids: String = exons.join(",");
        parts.push(exon_ids);

        parts.join(";")
    }
}

impl Clone for PositionAnnotation {
    fn clone(&self) -> Self {
        PositionAnnotation {
            genic_region: self.genic_region.clone(),
            gene_ids: self.gene_ids.clone(),
            transcript_ids: self.transcript_ids.clone(),
            exon_ids: self.exon_ids.clone()
        }
    }
}
