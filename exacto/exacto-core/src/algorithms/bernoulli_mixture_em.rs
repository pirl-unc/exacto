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


#[derive(Clone, Debug)]
struct Params {
    /// `mu[k][d] = P(alt allele at site d | cell k)`. Shape: `K × D`.
    mu: Vec<Vec<f64>>,
    /// `pi[k]` = mixing weight of cell `k` (prior P(a read is from cell k)).
    /// Sums to 1. Length: `K`.
    pi: Vec<f64>,
}

impl Params {
    fn n_components(&self) -> usize {
        self.pi.len()
    }
}

/// The result of one EM run.
#[derive(Clone, Debug)]
pub struct EmResult {
    pub params: Params,
    pub resp: Vec<Vec<f64>>,
    pub loglik_history: Vec<f64>
}

impl EmResult {
    /// Final (best) observed-data log-likelihood.
    fn final_loglik(&self) -> f64 {
        *self.loglik_history.last().unwrap()
    }
}


fn logsumexp(xs: &[f64]) -> f64 {
    let m = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY {
        return f64::NEG_INFINITY; // all -inf
    }
    let sum: f64 = xs.iter().map(|&x| (x - m).exp()).sum();
    m + sum.ln()
}


fn read_loglik(read: &[Option<f64>], mu_k: &[f64]) -> f64 {
    let mut ll = 0.0;
    for (d, &call) in read.iter().enumerate() {
        if let Some(x) = call {
            // log(μ) and log(1-μ); `ln_1p(-μ)` == log(1-μ) but more accurate as μ→1.
            ll += x * mu_k[d].ln() + (1.0 - x) * (-mu_k[d]).ln_1p();
        }
        // `None` => missing => contributes nothing (marginalized away).
    }
    ll
}


fn e_step(data: &ReadMatrix, params: &Params) -> (Vec<Vec<f64>>, f64) {
    let n = data.len();
    let k = params.n_components();

    let mut resp = vec![vec![0.0; k]; n];
    let mut total_loglik = 0.0;

    for i in 0..n {
        let mut log_joint = vec![0.0; k];
        for c in 0..k {
            log_joint[c] = params.pi[c].ln() + read_loglik(&data[i], &params.mu[c]);
        }

        let log_norm = logsumexp(&log_joint);

        for c in 0..k {
            resp[i][c] = (log_joint[c] - log_norm).exp();
        }

        total_loglik += log_norm;
    }

    (resp, total_loglik)
}


fn m_step(data: &ReadMatrix, resp: &[Vec<f64>]) -> Params {
    let eps = 1e-6;
    let n = data.len();
    let k = resp[0].len();
    let d_dim = data[0].len();

    let mut mu = vec![vec![0.0; d_dim]; k];
    for c in 0..k {
        for d in 0..d_dim {
            let mut num = 0.0; // Σ γ_ik x_id  over observed reads
            let mut den = 0.0; // Σ γ_ik       over observed reads
            for i in 0..n {
                if let Some(x) = data[i][d] {
                    num += resp[i][c] * x;
                    den += resp[i][c];
                }
            }
            let raw = if den > 0.0 { num / den } else { 0.5 };
            mu[c][d] = raw.clamp(eps, 1.0 - eps);
        }
    }

    let mut pi = vec![0.0; k];
    for c in 0..k {
        let s: f64 = (0..n).map(|i| resp[i][c]).sum();
        pi[c] = s / n as f64;
    }

    Params { mu, pi }
}


fn run_em(data: &ReadMatrix, init: Params, max_iter: usize, tol: f64) -> EmResult {
    let mut params = init;
    let mut prev_ll = f64::NEG_INFINITY;
    let mut history: Vec<f64> = Vec::new();
    let mut resp: Vec<Vec<f64>> = Vec::new();

    for _ in 0..max_iter {
        // E-step: soft labels + current log-likelihood.
        let (r, ll) = e_step(data, &params);
        resp = r;
        history.push(ll);

        // Converged? (log-likelihood stopped rising)
        if ll - prev_ll < tol {
            break;
        }
        prev_ll = ll;

        // M-step: re-estimate parameters from the soft labels.
        params = m_step(data, &resp);
    }

    EmResult {
        params,
        resp,
        loglik_history: history,
    }
}


fn run_em_restarts(
    data: &ReadMatrix,
    k: usize,
    n_starts: usize,
    seed: u64,
    max_iter: usize,
    tol: f64,
) -> EmResult {
    let d_dim = data[0].len();
    let mut rng = Pcg64::seed_from_u64(seed);
    let mut best: Option<EmResult> = None;

    for _ in 0..n_starts {
        let mu: Vec<Vec<f64>> = (0..k)
            .map(|_| (0..d_dim).map(|_| rng.random_range(0.2..0.8)).collect())
            .collect();
        let pi = vec![1.0 / k as f64; k];

        let result = run_em(data, Params { mu, pi }, max_iter, tol);

        let keep = match &best {
            None => true,
            Some(b) => result.final_loglik() > b.final_loglik(),
        };
        if keep {
            best = Some(result);
        }
    }

    best.expect("n_starts must be >= 1")
}


fn n_free_params(k: usize, d: usize) -> usize {
    k * d + (k - 1)
}


fn assignment_entropy(resp: &[Vec<f64>]) -> f64 {
    let mut en = 0.0;
    for row in resp {
        for &g in row {
            if g > 0.0 {
                en -= g * g.ln();
            }
        }
    }
    en
}


#[derive(Debug)]
pub struct ModelScores {
    pub k: usize,
    pub loglik: f64,
    pub n_params: usize,
    pub bic: f64,
    pub aic: f64,
    pub icl: f64
}


fn score_model(k: usize, fit: &EmResult, n: usize, d: usize) -> ModelScores {
    let loglik = fit.final_loglik();
    let p = n_free_params(k, d);
    let log_n = (n as f64).ln();
    let bic = -2.0 * loglik + p as f64 * log_n;
    let aic = -2.0 * loglik + 2.0 * p as f64;
    let icl = bic + 2.0 * assignment_entropy(&fit.resp);
    ModelScores { k, loglik, n_params: p, bic, aic, icl }
}


pub fn run_bernoulli_mixture_em(
    data: &ReadMatrix,
    kmax: usize,
    n_starts: usize,
    seed: u64,
    max_iter: usize,
    tol: f64
) -> Vec<(EmResult, ModelScores)> {
    let n = data.len();
    let d = data[0].len();
    (1..=kmax)
        .map(|k| {
            let fit = run_em_restarts(data, k, n_starts, seed, max_iter, tol);
            let scores = score_model(k, &fit, n, d);
            (fit, scores)
        })
        .collect()
}


pub fn argmin<I: IntoIterator<Item = f64>>(iter: I) -> usize {
    let mut best_i = 0;
    let mut best_v = f64::INFINITY;
    for (i, v) in iter.into_iter().enumerate() {
        if v < best_v {
            best_v = v;
            best_i = i;
        }
    }
    best_i
}


fn l1(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).sum()
}


fn align_to_truth(mu_hat: &[Vec<f64>], mu_true: &[Vec<f64>]) -> Vec<usize> {
    assert_eq!(mu_hat.len(), 2, "alignment helper assumes K = 2");
    let identity = l1(&mu_hat[0], &mu_true[0]) + l1(&mu_hat[1], &mu_true[1]);
    let swapped = l1(&mu_hat[0], &mu_true[1]) + l1(&mu_hat[1], &mu_true[0]);
    if identity <= swapped {
        vec![0, 1]
    } else {
        vec![1, 0]
    }
}


pub fn hard_label(resp_row: &[f64]) -> usize {
    let mut best = 0;
    for c in 1..resp_row.len() {
        if resp_row[c] > resp_row[best] {
            best = c;
        }
    }
    best
}


fn n_observed(read: &[Option<f64>]) -> usize {
    read.iter().filter(|c| c.is_some()).count()
}
