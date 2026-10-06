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


use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64;


pub type ReadMatrix = Vec<Vec<Option<f64>>>;


/// Binarize a present cell value (robust to any 0/1-ish encoding).
#[inline]
fn allele(x: f64) -> u8 {
    if x >= 0.5 {
        1
    } else {
        0
    }
}

/// Mismatches between one read and a consensus, counting ONLY the columns the
/// read actually observes (None cells are skipped).
fn distance(read: &[Option<f64>], consensus: &[u8]) -> u32 {
    let mut d = 0;
    for (c, cell) in read.iter().enumerate() {
        if let Some(v) = cell {
            if allele(*v) != consensus[c] {
                d += 1;
            }
        }
    }
    d
}

/// Column-wise majority vote over observed cells of the given member reads.
/// Ties and all-unknown columns default to 1 (an arbitrary but deterministic
/// tie-break; it does not change the MEC value at a tie).
fn consensus_of(matrix: &ReadMatrix, members: &[usize], n_vars: usize) -> Vec<u8> {
    let mut cons = vec![1u8; n_vars];
    for c in 0..n_vars {
        let (mut ones, mut zeros) = (0i32, 0i32);
        for &i in members {
            if let Some(v) = matrix[i][c] {
                if allele(v) == 1 {
                    ones += 1;
                } else {
                    zeros += 1;
                }
            }
        }
        cons[c] = if ones >= zeros { 1 } else { 0 };
    }
    cons
}

/// Assign every read to its closest centroid; return (labels, total_cost).
fn assign(matrix: &ReadMatrix, centroids: &[Vec<u8>]) -> (Vec<usize>, u32) {
    let n = matrix.len();
    let k = centroids.len();
    let mut labels = vec![0usize; n];
    let mut total = 0;
    for i in 0..n {
        let mut best = 0usize;
        let mut best_d = u32::MAX;
        for c in 0..k {
            let d = distance(&matrix[i], &centroids[c]);
            if d < best_d {
                best_d = d;
                best = c;
            }
        }
        labels[i] = best;
        total += best_d;
    }
    (labels, total)
}

/// One alternating-minimization run. Returns (mec_cost, labels, centroids).
///
/// `max_iter` caps the number of Lloyd (assign/update) iterations.
fn mec_one_run(
    matrix: &ReadMatrix,
    k: usize,
    n_vars: usize,
    max_iter: usize,
    rng: &mut Pcg64,
) -> (u32, Vec<usize>, Vec<Vec<u8>>) {
    let n = matrix.len();

    // ---- k-means++-style seeding on Hamming distance over observed cells ----
    // A seed's centroid is its own row binarized, computed once when the seed is chosen.
    let mut seeds: Vec<usize> = vec![rng.random_range(0..n)];
    let mut centroids: Vec<Vec<u8>> = vec![consensus_of(matrix, &[seeds[0]], n_vars)];
    while seeds.len() < k {
        // distance of each read to its nearest chosen seed
        let d: Vec<f64> = matrix
            .iter()
            .map(|read| centroids.iter().map(|cs| distance(read, cs)).min().unwrap() as f64)
            .collect();
        let sum: f64 = d.iter().sum();
        let pick = if sum == 0.0 {
            // all remaining reads identical to a seed: pick any unused index
            let mut cand = rng.random_range(0..n);
            while seeds.contains(&cand) {
                cand = rng.random_range(0..n);
            }
            cand
        } else {
            let mut r = rng.random::<f64>() * sum;
            let mut idx = n - 1;
            for i in 0..n {
                r -= d[i];
                if r <= 0.0 {
                    idx = i;
                    break;
                }
            }
            idx
        };
        seeds.push(pick);
        centroids.push(consensus_of(matrix, &[pick], n_vars));
    }

    // ---- Lloyd iterations ----
    let mut prev: Option<Vec<usize>> = None;
    let (mut labels, mut cost) = assign(matrix, &centroids);
    for _ in 0..max_iter {
        // One pass tallies every cluster's column votes; the majority rule is consensus_of's.
        let mut ones: Vec<Vec<i32>> = vec![vec![0; n_vars]; k];
        let mut zeros: Vec<Vec<i32>> = vec![vec![0; n_vars]; k];
        let mut sizes: Vec<usize> = vec![0; k];
        for (read, &label) in matrix.iter().zip(labels.iter()) {
            sizes[label] += 1;
            for (c, cell) in read.iter().enumerate() {
                if let Some(v) = cell {
                    if allele(*v) == 1 {
                        ones[label][c] += 1;
                    } else {
                        zeros[label][c] += 1;
                    }
                }
            }
        }
        for c in 0..k {
            if sizes[c] == 0 {
                // empty cluster: reseed to the current worst-fit read so we keep k clusters
                let mut worst = 0usize;
                let mut worst_d = 0u32;
                for i in 0..n {
                    let dd = distance(&matrix[i], &centroids[labels[i]]);
                    if dd >= worst_d {
                        worst_d = dd;
                        worst = i;
                    }
                }
                centroids[c] = consensus_of(matrix, &[worst], n_vars);
            } else {
                centroids[c] = (0..n_vars).map(|v| if ones[c][v] >= zeros[c][v] { 1 } else { 0 }).collect();
            }
        }
        // reassign
        let (new_labels, new_cost) = assign(matrix, &centroids);
        labels = new_labels;
        cost = new_cost;
        if prev.as_ref() == Some(&labels) {
            break;
        }
        prev = Some(labels.clone());
    }
    (cost, labels, centroids)
}


/// Best-of-`restarts` MEC solution for a fixed k.
///
/// One centroid is the column majority of every read whatever the seed, so k = 1 needs one
/// run. `restarts`, `max_iter` and `base_seed` behave as before: one Pcg64 per (k, r), the
/// lowest cost kept, the lowest r among equal costs.
fn solve_mec(
    matrix: &ReadMatrix,
    k: usize,
    n_vars: usize,
    restarts: usize,
    max_iter: usize,
    base_seed: u64,
) -> (u32, Vec<usize>, Vec<Vec<u8>>) {
    let restarts: usize = if k == 1 { restarts.min(1) } else { restarts };
    let mut best: Option<(u32, Vec<usize>, Vec<Vec<u8>>)> = None;
    for r in 0..restarts {
        let seed = base_seed ^ ((k as u64) << 32) ^ (r as u64).wrapping_mul(0x9E37_79B1);
        let mut rng = Pcg64::seed_from_u64(seed);
        let cand = mec_one_run(matrix, k, n_vars, max_iter, &mut rng);
        if best.as_ref().map_or(true, |b| cand.0 < b.0) {
            best = Some(cand);
        }
    }
    best.unwrap()
}


// --------------------------------------------------------------------------- //
// Model selection: BIC (primary), plus AIC for reference
// --------------------------------------------------------------------------- //
/// Full result for a single value of k. `select_clusters` returns one of these
/// per k, so callers can inspect every candidate — not just the winner.
pub struct KResult {
    pub k: usize,
    pub mec: u32,
    pub eps_hat: f64,
    pub bic: f64,
    pub aic: f64,
    pub labels: Vec<usize>,      // labels[i] = cluster index assigned to read i
    pub consensus: Vec<Vec<u8>>, // consensus[c] = haplotype of cluster c (1/0 per variant)
}

/// Under a per-cell Bernoulli error model, log-likelihood at a given MEC is
///     L = (N - MEC) * ln(1 - eps) + MEC * ln(eps),   eps_hat = MEC / N,
/// with the convention MEC = 0 -> L = 0. Free parameters p = k * n_vars + 1
/// (the consensus bits plus the error rate). Returns (eps_hat, bic, aic).
fn score_stats(k: usize, mec: u32, n_observed: usize, n_vars: usize) -> (f64, f64, f64) {
    let n = n_observed as f64;
    let mec_f = mec as f64;
    let eps_hat = if n > 0.0 { mec_f / n } else { 0.0 };
    let loglik = if mec == 0 || eps_hat <= 0.0 {
        0.0 // perfect fit; 0 * ln(0) taken as 0
    } else {
        let eps = eps_hat.min(0.5);
        (n - mec_f) * (1.0 - eps).ln() + mec_f * eps.ln()
    };
    let p = (k * n_vars + 1) as f64;
    let bic = -2.0 * loglik + p * n.ln();
    let aic = -2.0 * loglik + 2.0 * p;
    (eps_hat, bic, aic)
}


/// Index of the BIC-minimizing result (the selected number of clusters).
/// MEC decreases monotonically in k, so we pick by BIC, not by minimum MEC.
pub fn best_by_bic(results: &[KResult]) -> usize {
    results
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.bic.partial_cmp(&b.1.bic).unwrap())
        .map(|(i, _)| i)
        .unwrap()
}


/// Solve MEC for every k = 1..=max_k and return one `KResult` per k, each
/// including that k's MEC, BIC/AIC scores, read->cluster assignments (`labels`),
/// and per-cluster consensus haplotypes. This function does NOT print or pick a
/// winner — call `best_by_bic` / `print_report` on the returned vector for that.
///
/// `restarts`, `max_iter`, `base_seed` are forwarded to the per-k solver.
pub fn perform_minimum_error_correction(
    matrix: &ReadMatrix,
    max_k: usize,
    num_restarts: usize,
    max_iter: usize,
    seed: u64
) -> Vec<KResult> {
    assert!(max_k > 0, "max_k must be > 0");
    assert!(num_restarts > 0, "num_restarts must be > 0");
    let n_reads = matrix.len();
    let n_vars = matrix.first().map_or(0, |r| r.len());
    let n_observed: usize = matrix.iter().flatten().filter(|c| c.is_some()).count();
    let k_hi = max_k.min(n_reads); // can't have more clusters than reads

    let mut results = Vec::with_capacity(k_hi);
    for k in 1..=k_hi {
        let (mec, labels, consensus) = solve_mec(
            matrix,
            k,
            n_vars,
            num_restarts,
            max_iter,
            seed
        );
        let (eps_hat, bic, aic) = score_stats(k, mec, n_observed, n_vars);
        results.push(KResult {
            k,
            mec,
            eps_hat,
            bic,
            aic,
            labels,
            consensus,
        });
    }
    results
}