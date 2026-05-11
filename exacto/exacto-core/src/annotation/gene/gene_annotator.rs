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


use crate::prelude::*;


pub trait GeneAnnotator {
    fn get_assembly(&self) -> &str;
    fn get_version(&self) -> &str;
    fn get_gene_ids_at_locus(&self, chromosome: &str, position: u32) -> Vec<Box<str>>;
    fn get_gene_ids_overlapping_region(&self, chromosome: &str, start: u32, end: u32) -> Vec<Box<str>>;
    fn get_transcript_ids_overlapping_region(&self, chromosome: &str, start: u32, end: u32) -> Vec<Box<str>>;
    fn get_exon_ids_overlapping_region(&self, chromosome: &str, start: u32, end: u32) -> Vec<Box<str>>;
    fn get_gene(&self, gene_id: &str) -> Option<&Gene>;
    fn get_genes(&self) -> Vec<&Gene>;
    fn get_transcript(&self, transcript_id: &str) -> Option<&Transcript>;
    fn get_transcripts(&self) -> Vec<&Transcript>;
    fn get_exon(&self, transcript_id: &str, exon_id: &str) -> Option<&Exon>;
    fn get_exons(&self) -> Vec<&Exon>;
    fn rank_transcripts<'a>(&self, transcripts: Vec<&'a Transcript>) -> Vec<&'a Transcript>;
}
