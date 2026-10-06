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


use std::collections::{BTreeSet, HashMap, HashSet};


pub(crate) fn run_finite_mixture_em(
    clusters: &HashMap<usize, HashSet<Box<str>>>,
    pseudo_count: f64,
    max_iter: usize,
    tol: f64
) -> (bool, Vec<(usize, f64)>, Vec<f64>) {
    // Step 1. Get the cluster IDs.
    let cluster_ids: Vec<usize> = clusters.keys().map(|k| *k).collect();
    let k: usize = cluster_ids.len();

    // Step 2. Create an index of the cluster IDs.
    // HashMap<cluster ID, cluster index>
    let cluster_index: HashMap<usize, usize> = cluster_ids
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, c)| (c, i))
        .collect();

    // Step 3. Build each read's compatibility set S_r (as cluster indices).
    // The isoforms the read could have originated from:
    // |S_r| = 1 unique
    // |S_r| > 1 ambiguous
    // HashMap<read name, BTreeSet<cluster index>>
    let mut read_compatibility_sets: HashMap<Box<str>, BTreeSet<usize>> = HashMap::new();
    for (cluster_id, read_names) in clusters.iter() {
        for read_name in read_names.iter() {
            read_compatibility_sets
                .entry(read_name.clone())
                .or_insert(BTreeSet::new())
                .insert(*cluster_index.get(cluster_id).unwrap());
        }
    }

    // Step 4. Collapse reads that share an identical compatibility set into
    // equivalence classes, counting n_S reads per class.
    // HashMap<cluster indices, number of reads>
    let mut equivalence_classes: HashMap<Vec<usize>, f64> = HashMap::new();
    for set in read_compatibility_sets.values() {
        let key: Vec<usize> = set.iter().cloned().collect();
        *equivalence_classes.entry(key).or_insert(0.0) += 1.0;
    }

    // Sort for reproducible floating-point summation order.
    // Vec<(Vec<cluster index>, number of reads)>
    let mut classes: Vec<(Vec<usize>, f64)> = equivalence_classes.into_iter().collect();
    classes.sort_by(|a, b| a.0.cmp(&b.0));

    // Get the total number of reads.
    let num_total_reads: f64 = classes.iter().map(|(_, n)| *n).sum();

    // Step 5. Count each cluster's unique support: reads whose class has |S| = 1
    // (reads that have come from single cluster).
    let mut unique_count: Vec<f64> = vec![0.0f64; k];
    for (set, num_reads) in &classes {
        if set.len() == 1 {
            unique_count[set[0]] += *num_reads;
        }
    }

    // Step 6. Initialize theta from unique support (a warm start).
    // The pseudocount keeps clusters with zero unique reads alive, and collapses θ to uniform
    // when no cluster has any unique support.
    let init_sum: f64 = unique_count.iter().sum::<f64>() + k as f64 * pseudo_count;
    let mut theta: Vec<f64> = unique_count.iter().map(|&u| (u + pseudo_count) / init_sum).collect();

    // Step 7. Run expectation maximization (EM).
    let mut converged: bool = false;
    let mut num_iterations: usize = 0;
    let mut log_likelihood_values: Vec<f64> = Vec::new();
    for _ in 0..max_iter {
        num_iterations += 1;

        // Compute the objective at the current theta.
        log_likelihood_values.push(log_likelihood(&classes, &theta));

        // E-step: responsibilities under the current theta.
        let gamma: Vec<Vec<f64>> = e_step(&classes, &theta);

        // M-step: expected counts -> renormalized new theta.
        let theta_new: Vec<f64> = m_step(&classes, &gamma, k, num_total_reads);

        // Check for convergence: compute the largest per-cluster change in theta.
        let delta = theta_new
            .iter()
            .zip(&theta)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);

        theta = theta_new;

        if delta < tol {
            converged = true;
            break;
        }
    }

    // Abundance alpha_c = N * theta_c (expected reads per cluster; equals the final M-step's alpha).
    let abundance: Vec<f64> = theta.iter().map(|&t| t * num_total_reads).collect();

    // Counts per million: CPM_c = alpha_c / N * 1e6 = theta_c * 1e6.
    let cpm: Vec<f64> = abundance.iter().map(|&a| a / num_total_reads * 1e6).collect();

    // Step 8. Pair each cluster's CPM with its original cluster ID (index c <-> cluster_ids[c]),
    // and return it with the convergence flag and the per-iteration log-likelihoods.
    let cluster_cpm: Vec<(usize, f64)> = (0..k)
        .map(|c| (cluster_ids[c], cpm[c]))
        .collect();

    (converged, cluster_cpm, log_likelihood_values)
}


fn e_step(
    classes: &Vec<(Vec<usize>, f64)>,
    theta: &Vec<f64>
) -> Vec<Vec<f64>> {
    let mut gamma: Vec<Vec<f64>> = Vec::new();
    for (set, num_reads) in classes.iter() {
        // Read total probability P(r) = sum(theta_c).
        // This is the responsibility normalizer.
        let z: f64 = set.iter().map(|&c| theta[c]).sum();
        if z > 0.0 {
            // Responsibility γ_c = theta_c / total probability for each c in the class.
            gamma.push(set.iter().map(|&c| theta[c] / z).collect());
        } else {
            // Degenerate (should not happen given the pseudocount): uniform.
            gamma.push(vec![1.0 / set.len() as f64; set.len()]);
        }
    }
    gamma
}


fn m_step(
    classes: &Vec<(Vec<usize>, f64)>,
    gamma: &Vec<Vec<f64>>,
    k: usize,
    num_total_reads: f64
) -> Vec<f64> {
    // Expected counts alpha_c: each class hands its n_S reads to its compatible
    // clusters, split by the E-step responsibilities.
    let mut alpha: Vec<f64> = vec![0.0f64; k];
    for ((set, n), weights) in classes.iter().zip(gamma) {
        for (&c, &g) in set.iter().zip(weights) {
            alpha[c] += *n * g;
        }
    }

    // Normalize the expected counts onto the simplex (lambda = N): theta_c = alpha_c / N.
    alpha.iter().map(|&a| a / num_total_reads).collect()
}


fn log_likelihood(classes: &Vec<(Vec<usize>, f64)>, theta: &Vec<f64>) -> f64 {
    classes
        .iter()
        .map(|(set, n)| {
            let z: f64 = set.iter().map(|&c| theta[c]).sum();
            if z > 0.0 {
                *n * z.ln()
            } else {
                0.0
            }
        })
        .sum()
}


#[cfg(test)]
#[path = "../tests/quantification/expectation_maximization.rs"]
mod tests;