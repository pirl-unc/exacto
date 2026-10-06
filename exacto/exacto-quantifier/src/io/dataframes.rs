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


use polars::prelude::*;

use crate::prelude::*;


pub fn quantification_cluster_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = QuantificationClusterRecord>
{
    let mut cluster_id: Vec<u32> = Vec::new();
    let mut gene_annotation: Vec<String> = Vec::new();
    let mut reference_gene_id: Vec<String> = Vec::new();
    let mut transcript_annotation: Vec<String> = Vec::new();
    let mut reference_transcript_id: Vec<String> = Vec::new();
    let mut allele: Vec<String> = Vec::new();
    let mut num_reads: Vec<u32> = Vec::new();
    let mut cpm: Vec<f64> = Vec::new();

    for r in records {
        cluster_id.push(r.cluster_id as u32);
        gene_annotation.push(r.gene_annotation.to_string());
        reference_gene_id.push(r.reference_gene_id.to_string());
        transcript_annotation.push(r.transcript_annotation.to_string());
        reference_transcript_id.push(r.reference_transcript_id.to_string());
        allele.push(r.allele.to_string());
        num_reads.push(r.num_reads as u32);
        cpm.push(r.cpm);
    }

    DataFrame::new(vec![
        Column::from(Series::new("cluster_id".into(), cluster_id)),
        Column::from(Series::new("gene_annotation".into(), gene_annotation)),
        Column::from(Series::new("reference_gene_id".into(), reference_gene_id)),
        Column::from(Series::new("transcript_annotation".into(), transcript_annotation)),
        Column::from(Series::new("reference_transcript_id".into(), reference_transcript_id)),
        Column::from(Series::new("allele".into(), allele)),
        Column::from(Series::new("num_reads".into(), num_reads)),
        Column::from(Series::new("cpm".into(), cpm))
    ]).unwrap()
}


pub fn quantification_reference_transcript_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = QuantificationReferenceTranscriptRecord>
{
    let mut reference_gene_id: Vec<String> = Vec::new();
    let mut reference_transcript_id: Vec<String> = Vec::new();
    let mut num_clusters: Vec<u32> = Vec::new();
    let mut num_reads: Vec<u32> = Vec::new();
    let mut cpm: Vec<f64> = Vec::new();

    for r in records {
        reference_gene_id.push(r.reference_gene_id.to_string());
        reference_transcript_id.push(r.reference_transcript_id.to_string());
        num_clusters.push(r.num_clusters as u32);
        num_reads.push(r.num_reads as u32);
        cpm.push(r.cpm);
    }

    DataFrame::new(vec![
        Column::from(Series::new("reference_gene_id".into(), reference_gene_id)),
        Column::from(Series::new("reference_transcript_id".into(), reference_transcript_id)),
        Column::from(Series::new("num_clusters".into(), num_clusters)),
        Column::from(Series::new("num_reads".into(), num_reads)),
        Column::from(Series::new("cpm".into(), cpm))
    ]).unwrap()
}
