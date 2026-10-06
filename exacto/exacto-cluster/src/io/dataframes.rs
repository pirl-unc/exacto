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


pub fn rna_read_cluster_id_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = RNAReadClusterIDRecord>
{
    let mut cluster_id: Vec<u32> = Vec::new();
    let mut read_name: Vec<String> = Vec::new();
    let mut is_shared: Vec<bool> = Vec::new();

    for r in records {
        cluster_id.push(r.cluster_id as u32);
        read_name.push(r.read_name.to_string());
        is_shared.push(r.is_shared);
    }

    DataFrame::new(vec![
        Column::from(Series::new("cluster_id".into(), cluster_id)),
        Column::from(Series::new("read_name".into(), read_name)),
        Column::from(Series::new("is_shared".into(), is_shared))
    ]).unwrap()
}


pub fn rna_read_cluster_summary_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = RNAReadClusterSummaryRecord>
{
    let mut cluster_id: Vec<u32> = Vec::new();
    let mut num_reads: Vec<u32> = Vec::new();

    for r in records {
        cluster_id.push(r.cluster_id as u32);
        num_reads.push(r.num_reads as u32);
    }

    DataFrame::new(vec![
        Column::from(Series::new("cluster_id".into(), cluster_id)),
        Column::from(Series::new("num_reads".into(), num_reads))
    ]).unwrap()
}


pub fn rna_read_cluster_reference_gene_transcript_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item =RNAReadClusterReferenceTranscriptRecord>
{
    let mut cluster_id: Vec<u32> = Vec::new();
    let mut reference_gene_id: Vec<String> = Vec::new();
    let mut reference_transcript_id: Vec<String> = Vec::new();

    for r in records {
        cluster_id.push(r.cluster_id as u32);
        reference_gene_id.push(r.reference_gene_id.to_string());
        reference_transcript_id.push(r.reference_transcript_id.to_string());
    }

    DataFrame::new(vec![
        Column::from(Series::new("cluster_id".into(), cluster_id)),
        Column::from(Series::new("reference_gene_id".into(), reference_gene_id)),
        Column::from(Series::new("reference_transcript_id".into(), reference_transcript_id))
    ]).unwrap()
}


pub fn rna_read_cluster_splice_junction_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = RNAReadClusterSpliceJunctionRecord>
{
    let mut cluster_id: Vec<u32> = Vec::new();
    let mut chromosome_1: Vec<String> = Vec::new();
    let mut chromosome_2: Vec<String> = Vec::new();
    let mut position_1: Vec<u32> = Vec::new();
    let mut position_2: Vec<u32> = Vec::new();
    let mut strand_1: Vec<String> = Vec::new();
    let mut strand_2: Vec<String> = Vec::new();

    for r in records {
        cluster_id.push(r.cluster_id as u32);
        chromosome_1.push(r.chromosome_1.to_string());
        chromosome_2.push(r.chromosome_2.to_string());
        position_1.push(r.position_1);
        position_2.push(r.position_2);
        strand_1.push(r.strand_1.to_string());
        strand_2.push(r.strand_2.to_string());
    }

    DataFrame::new(vec![
        Column::from(Series::new("cluster_id".into(), cluster_id)),
        Column::from(Series::new("chromosome_1".into(), chromosome_1)),
        Column::from(Series::new("position_1".into(), position_1)),
        Column::from(Series::new("strand_1".into(), strand_1)),
        Column::from(Series::new("chromosome_2".into(), chromosome_2)),
        Column::from(Series::new("position_2".into(), position_2)),
        Column::from(Series::new("strand_2".into(), strand_2))
    ]).unwrap()
}


pub fn rna_read_cluster_variant_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = RNAReadClusterVariantRecord>
{
    let mut cluster_id: Vec<u32> = Vec::new();
    let mut chromosome_1: Vec<String> = Vec::new();
    let mut chromosome_2: Vec<String> = Vec::new();
    let mut position_1: Vec<u32> = Vec::new();
    let mut position_2: Vec<u32> = Vec::new();
    let mut strand_1: Vec<String> = Vec::new();
    let mut strand_2: Vec<String> = Vec::new();
    let mut operation_type_1: Vec<String> = Vec::new();
    let mut operation_type_2: Vec<String> = Vec::new();
    let mut sequence: Vec<String> = Vec::new();
    let mut variant_type: Vec<String> = Vec::new();
    
    for r in records {
        cluster_id.push(r.cluster_id as u32);
        chromosome_1.push(r.chromosome_1.to_string());
        chromosome_2.push(r.chromosome_2.to_string());
        position_1.push(r.position_1);
        position_2.push(r.position_2);
        strand_1.push(r.strand_1.to_string());
        strand_2.push(r.strand_2.to_string());
        operation_type_1.push(r.operation_type_1.to_string());
        operation_type_2.push(r.operation_type_2.to_string());
        sequence.push(r.sequence.to_string());
        variant_type.push(r.variant_type.to_string());
    }

    DataFrame::new(vec![
        Column::from(Series::new("cluster_id".into(), cluster_id)),
        Column::from(Series::new("chromosome_1".into(), chromosome_1)),
        Column::from(Series::new("position_1".into(), position_1)),
        Column::from(Series::new("strand_1".into(), strand_1)),
        Column::from(Series::new("operation_type_1".into(), operation_type_1)),
        Column::from(Series::new("chromosome_2".into(), chromosome_2)),
        Column::from(Series::new("position_2".into(), position_2)),
        Column::from(Series::new("strand_2".into(), strand_2)),
        Column::from(Series::new("operation_type_2".into(), operation_type_2)),
        Column::from(Series::new("sequence".into(), sequence)),
        Column::from(Series::new("variant_type".into(), variant_type))
    ]).unwrap()
}


pub fn rna_read_cluster_failed_variant_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = RNAReadClusterFailedVariantRecord>
{
    let mut cluster_id: Vec<Option<u32>> = Vec::new();
    let mut chromosome_1: Vec<String> = Vec::new();
    let mut position_1: Vec<u32> = Vec::new();
    let mut strand_1: Vec<String> = Vec::new();
    let mut operation_type_1: Vec<String> = Vec::new();
    let mut chromosome_2: Vec<String> = Vec::new();
    let mut position_2: Vec<u32> = Vec::new();
    let mut strand_2: Vec<String> = Vec::new();
    let mut operation_type_2: Vec<String> = Vec::new();
    let mut sequence: Vec<String> = Vec::new();
    let mut variant_type: Vec<String> = Vec::new();
    let mut failure_reason: Vec<String> = Vec::new();
    let mut alt_count: Vec<u32> = Vec::new();
    let mut reference_count: Vec<u32> = Vec::new();
    let mut coverage: Vec<Option<u32>> = Vec::new();
    let mut p_value: Vec<Option<f64>> = Vec::new();
    let mut p_value_cutoff: Vec<Option<f64>> = Vec::new();
    let mut expected_error_rate: Vec<Option<f64>> = Vec::new();
    let mut slippage_repeat_length: Vec<Option<u32>> = Vec::new();
    let mut junction_homology: Vec<Option<u32>> = Vec::new();

    for r in records {
        cluster_id.push(r.cluster_id.map(|id| id as u32));
        chromosome_1.push(r.chromosome_1.to_string());
        position_1.push(r.position_1);
        strand_1.push(r.strand_1.to_string());
        operation_type_1.push(r.operation_type_1.to_string());
        chromosome_2.push(r.chromosome_2.to_string());
        position_2.push(r.position_2);
        strand_2.push(r.strand_2.to_string());
        operation_type_2.push(r.operation_type_2.to_string());
        sequence.push(r.sequence.to_string());
        variant_type.push(r.variant_type.to_string());
        failure_reason.push(r.failure_reason.to_string());
        alt_count.push(r.alt_count);
        reference_count.push(r.reference_count);
        coverage.push(r.coverage);
        p_value.push(r.p_value);
        p_value_cutoff.push(r.p_value_cutoff);
        expected_error_rate.push(r.expected_error_rate);
        slippage_repeat_length.push(r.slippage_repeat_length);
        junction_homology.push(r.junction_homology);
    }

    DataFrame::new(vec![
        Column::from(Series::new("cluster_id".into(), cluster_id)),
        Column::from(Series::new("chromosome_1".into(), chromosome_1)),
        Column::from(Series::new("position_1".into(), position_1)),
        Column::from(Series::new("strand_1".into(), strand_1)),
        Column::from(Series::new("operation_type_1".into(), operation_type_1)),
        Column::from(Series::new("chromosome_2".into(), chromosome_2)),
        Column::from(Series::new("position_2".into(), position_2)),
        Column::from(Series::new("strand_2".into(), strand_2)),
        Column::from(Series::new("operation_type_2".into(), operation_type_2)),
        Column::from(Series::new("sequence".into(), sequence)),
        Column::from(Series::new("variant_type".into(), variant_type)),
        Column::from(Series::new("failure_reason".into(), failure_reason)),
        Column::from(Series::new("alt_count".into(), alt_count)),
        Column::from(Series::new("reference_count".into(), reference_count)),
        Column::from(Series::new("coverage".into(), coverage)),
        Column::from(Series::new("p_value".into(), p_value)),
        Column::from(Series::new("p_value_cutoff".into(), p_value_cutoff)),
        Column::from(Series::new("expected_error_rate".into(), expected_error_rate)),
        Column::from(Series::new("slippage_repeat_length".into(), slippage_repeat_length)),
        Column::from(Series::new("junction_homology".into(), junction_homology))
    ]).unwrap()
}


pub fn rna_read_cluster_template_switch_records_to_dataframe<I>(records: I) -> DataFrame
where
    I: IntoIterator<Item = RNAReadClusterTemplateSwitchRecord>
{
    let mut cluster_id: Vec<u32> = Vec::new();
    let mut chromosome_1: Vec<String> = Vec::new();
    let mut position_1: Vec<u32> = Vec::new();
    let mut strand_1: Vec<String> = Vec::new();
    let mut operation_type_1: Vec<String> = Vec::new();
    let mut chromosome_2: Vec<String> = Vec::new();
    let mut position_2: Vec<u32> = Vec::new();
    let mut strand_2: Vec<String> = Vec::new();
    let mut operation_type_2: Vec<String> = Vec::new();
    let mut sequence: Vec<String> = Vec::new();
    let mut variant_type: Vec<String> = Vec::new();
    let mut homology_left: Vec<Option<u32>> = Vec::new();
    let mut homology_right: Vec<Option<u32>> = Vec::new();
    let mut canonical_splice: Vec<Option<bool>> = Vec::new();
    let mut breakpoint_dispersion: Vec<Option<u32>> = Vec::new();
    let mut foldback: Vec<bool> = Vec::new();
    let mut partner_count: Vec<Option<u32>> = Vec::new();
    let mut num_reads: Vec<u32> = Vec::new();
    let mut flagged: Vec<bool> = Vec::new();
    let mut suppressed: Vec<bool> = Vec::new();

    for r in records {
        cluster_id.push(r.cluster_id as u32);
        chromosome_1.push(r.chromosome_1.to_string());
        position_1.push(r.position_1);
        strand_1.push(r.strand_1.to_string());
        operation_type_1.push(r.operation_type_1.to_string());
        chromosome_2.push(r.chromosome_2.to_string());
        position_2.push(r.position_2);
        strand_2.push(r.strand_2.to_string());
        operation_type_2.push(r.operation_type_2.to_string());
        sequence.push(r.sequence.to_string());
        variant_type.push(r.variant_type.to_string());
        homology_left.push(r.homology_left);
        homology_right.push(r.homology_right);
        canonical_splice.push(r.canonical_splice);
        breakpoint_dispersion.push(r.breakpoint_dispersion);
        foldback.push(r.foldback);
        partner_count.push(r.partner_count);
        num_reads.push(r.num_reads);
        flagged.push(r.flagged);
        suppressed.push(r.suppressed);
    }

    DataFrame::new(vec![
        Column::from(Series::new("cluster_id".into(), cluster_id)),
        Column::from(Series::new("chromosome_1".into(), chromosome_1)),
        Column::from(Series::new("position_1".into(), position_1)),
        Column::from(Series::new("strand_1".into(), strand_1)),
        Column::from(Series::new("operation_type_1".into(), operation_type_1)),
        Column::from(Series::new("chromosome_2".into(), chromosome_2)),
        Column::from(Series::new("position_2".into(), position_2)),
        Column::from(Series::new("strand_2".into(), strand_2)),
        Column::from(Series::new("operation_type_2".into(), operation_type_2)),
        Column::from(Series::new("sequence".into(), sequence)),
        Column::from(Series::new("variant_type".into(), variant_type)),
        Column::from(Series::new("homology_left".into(), homology_left)),
        Column::from(Series::new("homology_right".into(), homology_right)),
        Column::from(Series::new("canonical_splice".into(), canonical_splice)),
        Column::from(Series::new("breakpoint_dispersion".into(), breakpoint_dispersion)),
        Column::from(Series::new("foldback".into(), foldback)),
        Column::from(Series::new("partner_count".into(), partner_count)),
        Column::from(Series::new("num_reads".into(), num_reads)),
        Column::from(Series::new("flagged".into(), flagged)),
        Column::from(Series::new("suppressed".into(), suppressed))
    ]).unwrap()
}
