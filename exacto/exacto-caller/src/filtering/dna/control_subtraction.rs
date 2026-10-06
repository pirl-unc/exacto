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


use exacto_core::prelude::ReadSupport;
use rayon::prelude::*;
use rayon::ThreadPool;
use std::sync::Arc;

use crate::prelude::*;
use crate::filtering::dna::control_subtraction_filter::DNAControlSubtractionFilter;
use crate::filtering::dna::control_variant_index::DNAControlVariantIndex;


pub(crate) fn diff_dna_variant_records(
    a: Vec<Arc<VariantRecord>>,
    b: Vec<Arc<VariantRecord>>,
    bin_size: u32,
    num_threads: usize,
    min_size_proportion: f64,
    max_ins_norm_edit_distance: f64,
    max_clustering_distance: u32,
    sequencing_error: f64,
    min_b_read_support: ReadSupport,
    apply_infinite_sites_assumption: bool,
    stranded: bool
) -> Vec<Arc<VariantRecord>> {
    if b.is_empty() {
        return a;
    }
    
    let control_variant_index: DNAControlVariantIndex = DNAControlVariantIndex::new(
        &b,
        bin_size,
        min_b_read_support
    );
    
    let control_subtraction_filter: DNAControlSubtractionFilter = DNAControlSubtractionFilter::new(
        &control_variant_index,
        min_size_proportion,
        max_ins_norm_edit_distance,
        max_clustering_distance,
        sequencing_error,
        apply_infinite_sites_assumption,
        stranded
    );
    
    let thread_pool: ThreadPool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();
    
    thread_pool.install(|| {
        a.into_par_iter()
            .with_max_len(1)
            .filter(|variant_record| control_subtraction_filter.passes(variant_record))
            .collect()
    })
}


#[cfg(test)]
#[path = "../../tests/filtering/dna/control_subtraction.rs"]
mod tests;
