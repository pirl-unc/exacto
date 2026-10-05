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


use csv::Error;
use exacto_core::prelude::write_tsv_table;
use std::fs;
use std::path::{Path, PathBuf};

use crate::prelude::*;


pub fn write_rna_read_cluster_set(
    cluster_set: &RNAReadClusterSet,
    output_dir: &str,
    output_prefix: &str
) -> Result<(), Error> {
    fs::create_dir_all(output_dir)?;
    let tsv_file = |suffix: &str| -> PathBuf {
        Path::new(output_dir).join(format!("{}_exacto_rna_clusters{}.tsv", output_prefix, suffix))
    };

    write_tsv_table(build_rna_read_cluster_id_records(cluster_set), &tsv_file(""))?;
    write_tsv_table(build_rna_read_cluster_summary_records(cluster_set), &tsv_file("_summary"))?;
    write_tsv_table(build_rna_read_cluster_reference_gene_transcript_records(cluster_set), &tsv_file("_reference_transcripts"))?;
    write_tsv_table(build_rna_read_cluster_splice_junction_records(cluster_set), &tsv_file("_splice_junctions"))?;
    write_tsv_table(build_rna_read_cluster_variant_records(cluster_set), &tsv_file("_variants_passed"))?;
    write_tsv_table(build_rna_read_cluster_failed_variant_records(cluster_set), &tsv_file("_variants_failed"))?;
    write_tsv_table(build_rna_read_cluster_template_switch_records(cluster_set), &tsv_file("_template_switch"))?;
    Ok(())
}
