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

use crate::pipeline::options::QuantifyRNAAbundancesOptions;
use crate::prelude::*;
use crate::quantification::expectation_maximization::run_finite_mixture_em;


pub fn quantify_rna_abundances(
    clusters: &HashMap<usize, HashSet<Box<str>>>,
    options: &QuantifyRNAAbundancesOptions
) -> ClusterQuantificationSet {
    let (converged, cluster_cpm, loglik) = run_finite_mixture_em(
        clusters,
        options.pseudo_count,
        options.max_iter,
        options.tol
    );

    let mut set: ClusterQuantificationSet = ClusterQuantificationSet::new();
    for (cluster_id, cpm) in cluster_cpm.iter()  {
        let quantification: ClusterQuantification = ClusterQuantification::new(
            *cluster_id,
            *cpm,
            clusters.get(cluster_id).unwrap().clone()
        );
        set.add(quantification);
    }

    set
}


#[cfg(test)]
#[path = "../tests/pipeline/rna_quantification.rs"]
mod tests;