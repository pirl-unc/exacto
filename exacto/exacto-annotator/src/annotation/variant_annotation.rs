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


use exacto_core::prelude::*;
use polars::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use crate::prelude::*;


pub fn read_variant_calls_tsv_file(tsv_file: &str) -> DataFrame {
    let schema_overwrite: Schema = Schema::from_iter(vec![
        Field::new("chromosome_1".into(), DataType::String),
        Field::new("chromosome_2".into(), DataType::String)
    ]);
    let parse_options = CsvParseOptions::default()
        .with_separator(b'\t');
    CsvReadOptions::default()
        .with_parse_options(parse_options)
        .with_has_header(true)
        .with_raise_if_empty(false)
        .with_schema_overwrite(Some(Arc::new(schema_overwrite)))
        .try_into_reader_with_file_path(Some(tsv_file.into()))
        .unwrap()
        .finish()
        .unwrap()
}


pub(crate) fn annotate_position(
    reference_chromosome_name: &str,
    reference_position: ReferencePosition,
    gene_annotator: &(impl GeneAnnotator + Sync)
) -> PositionAnnotation {
    let gene_ids: Vec<Box<str>> = gene_annotator.get_gene_ids_at_locus(&*reference_chromosome_name, reference_position);
    let mut transcript_ids: HashMap<Box<str>,HashSet<Box<str>>> = HashMap::new();
    let mut exon_ids: HashMap<Box<str>,Box<str>> = HashMap::new();

    let pos: isize = reference_position as isize;
    for gene_id in &gene_ids {
        let Some(gene) = gene_annotator.get_gene(&*gene_id) else { continue };
        for transcript_id in gene.get_transcript_ids() {
            let Some(transcript) = gene_annotator.get_transcript(&*transcript_id) else { continue };

            let mut record_transcript: bool = overlaps(pos, pos, transcript.start as isize, transcript.end as isize);

            for exon_id in transcript.get_exon_ids() {
                let Some(exon) = gene_annotator.get_exon(&*transcript_id, &*exon_id) else { continue };
                if overlaps(pos, pos, exon.start as isize, exon.end as isize) {
                    record_transcript = true;
                    exon_ids.insert(transcript_id.clone(), exon_id.clone());
                }
            }

            if record_transcript {
                transcript_ids
                    .entry(gene_id.clone())
                    .or_default()
                    .insert(transcript_id.clone());
            }
        }
    }

    let mut genic_region: GenicRegion = GenicRegion::Intergenic;
    if !transcript_ids.is_empty() {
        genic_region = GenicRegion::Intronic;
    }
    if !exon_ids.is_empty() {
        genic_region = GenicRegion::Exonic;
    }

    let mut position_annotation: PositionAnnotation = PositionAnnotation::new(genic_region);

    for gene_id in gene_ids.iter() {
        position_annotation.add_reference_gene_id(gene_id.clone());
    }
    for (gene_id, transcript_ids) in transcript_ids.iter() {
        for transcript_id in transcript_ids.iter() {
            position_annotation.add_reference_transcript_id(gene_id.clone(), transcript_id.clone());
        }
    }
    for (transcript_id, exon_id) in exon_ids.iter() {
        position_annotation.add_reference_exon_id(transcript_id.clone(), exon_id.clone());
    }

    position_annotation
}
