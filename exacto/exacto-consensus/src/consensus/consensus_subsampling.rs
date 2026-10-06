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


use rand::SeedableRng;
use rand_pcg::Pcg64;
use std::collections::HashSet;


pub fn subsample_read_names(
    read_names: &HashSet<Box<str>>,
    max_reads_per_cluster: usize,
    seed: u64
) -> Vec<&str> {
    let mut sorted_read_names: Vec<&str> = read_names.iter().map(|read_name| &**read_name).collect();
    sorted_read_names.sort_unstable();
    if max_reads_per_cluster == 0 || sorted_read_names.len() <= max_reads_per_cluster {
        return sorted_read_names;
    }
    let mut rng: Pcg64 = Pcg64::seed_from_u64(seed);
    let mut indices: Vec<usize> = rand::seq::index::sample(
        &mut rng,
        sorted_read_names.len(),
        max_reads_per_cluster
    ).into_vec();
    indices.sort_unstable();
    indices.into_iter().map(|index| sorted_read_names[index]).collect()
}


#[cfg(test)]
#[path = "../tests/consensus/consensus_subsampling.rs"]
mod tests;