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
use exacto_caller::prelude::*;
use std::collections::{HashMap, HashSet};


fn to_read_matrix(
    variants_to_reads: &HashMap<usize, Vec<(usize, Allele)>>
) -> (Vec<usize>, Vec<usize>, Vec<Vec<Option<f64>>>) {
    // Deterministic column and row order (HashMap iteration order is not stable).
    let mut col_variants: Vec<usize> = variants_to_reads.keys().copied().collect();
    col_variants.sort_unstable();
    let mut row_read_ids: Vec<usize> = variants_to_reads
        .values()
        .flatten()
        .map(|(read_id, _)| *read_id)
        .collect();
    row_read_ids.sort_unstable();
    row_read_ids.dedup();
    let row_of: HashMap<usize, usize> = row_read_ids
        .iter()
        .enumerate()
        .map(|(row, read_id)| (*read_id, row))
        .collect();

    // Fill each column straight from its genotype rows; a missing or NotCovered cell stays None.
    let mut matrix: Vec<Vec<Option<f64>>> = vec![vec![None; col_variants.len()]; row_read_ids.len()];
    for (column, variant_id) in col_variants.iter().enumerate() {
        for (read_id, allele) in variants_to_reads[variant_id].iter() {
            matrix[row_of[read_id]][column] = match allele {
                Allele::Alternate => Some(1.0),
                Allele::Reference => Some(0.0),
                Allele::NotCovered => None
            };
        }
    }

    (row_read_ids, col_variants, matrix)
}


fn tells(row: &[Option<f64>], a: usize, b: usize, consensus: &[Vec<u8>], is_observed: &[Vec<bool>]) -> bool {
    row.iter().enumerate().any(|(site, allele)| {
        allele.is_some() && is_observed[a][site] && is_observed[b][site] && consensus[a][site] != consensus[b][site]
    })
}


pub(crate) fn phase_genotypes(
    read_ids: &HashSet<usize>,
    genotypes: &HashMap<usize, Vec<(usize, Allele)>>,
    min_reads: usize,
    max_k: usize,
    num_restarts: usize,
    max_iter: usize,
    seed: u64
) -> (HashMap<usize, HashSet<usize>>, HashMap<usize, HashSet<usize>>) {
    // Per-cluster and pure: build matrix -> fit + select K -> hard-assign.
    let (row_read_ids, _col_variants, matrix) = to_read_matrix(
        genotypes
    );

    // The sites that could set two cells apart, and the number of cells the floor leaves room
    // for. `min_reads` of 0 counts as 1.
    let max_k: usize = max_k.min(read_ids.len() / min_reads.max(1));
    let phased_sites: Vec<usize> = (0..matrix.first().map_or(0, |row| row.len()))
        .filter(|&site| {
            let num_alternate: usize = matrix.iter().filter(|row| row[site] == Some(1.0)).count();
            let num_reference: usize = matrix.iter().filter(|row| row[site] == Some(0.0)).count();
            num_alternate >= min_reads && num_reference >= min_reads
        })
        .collect();
    let matrix: Vec<Vec<Option<f64>>> = matrix
        .iter()
        .map(|row| phased_sites.iter().map(|&site| row[site]).collect())
        .collect();

    // Nothing to phase (no reads, no site that could part two cells, or no room for two).
    // With no sites there is no evidence that any two reads differ so the cluster is one just cell.
    if row_read_ids.is_empty() || matrix[0].is_empty() || max_k < 2 {
        return (HashMap::from([(0, read_ids.clone())]), HashMap::new());
    }

    // Every cluster read has a row, NotCovered included, so no read can be dropped on
    // the way through the matrix. Genotyping guarantees this; check it where it matters.
    assert!(
        row_read_ids.len() == read_ids.len() && row_read_ids.iter().all(|read_id| read_ids.contains(read_id)),
        "Genotype rows must be exactly the cluster's reads."
    );

    // Perform minimum error correction.
    let results: Vec<KResult> = perform_minimum_error_correction(
        &matrix,
        max_k,
        num_restarts,
        max_iter,
        seed
    );

    // Group reads by their assigned cell, for every number of cells: cell index -> read IDs.
    // (BIC, result, per cell and site whether a read of the cell observes the site, cells)
    let mut best: Option<(f64, &KResult, Vec<Vec<bool>>, HashMap<usize, HashSet<usize>>)> = None;
    for result in results.iter() {
        let mut cells: HashMap<usize, HashSet<usize>> = HashMap::new();
        for (read_id, cell) in row_read_ids.iter().zip(result.labels.iter()) {
            cells.entry(*cell).or_default().insert(*read_id);
        }
        // Vec<per cell, per site, whether a read of the cell observes the site>
        let mut is_observed: Vec<Vec<bool>> = vec![vec![false; phased_sites.len()]; result.consensus.len()];
        for (cell, row) in result.labels.iter().zip(matrix.iter()) {
            for (site, allele) in row.iter().enumerate() {
                is_observed[*cell][site] |= allele.is_some();
            }
        }
        // HashMap<cell index, number of reads telling it from every other cell>
        let mut num_telling_reads: HashMap<usize, usize> = HashMap::new();
        for (cell, row) in result.labels.iter().zip(matrix.iter()) {
            let is_telling: bool = cells
                .keys()
                .filter(|other| *other != cell)
                .all(|other| tells(row, *cell, *other, &result.consensus, &is_observed));
            if is_telling {
                *num_telling_reads.entry(*cell).or_insert(0) += 1;
            }
        }
        let is_under_min_reads: bool = cells
            .keys()
            .any(|cell| num_telling_reads.get(cell).copied().unwrap_or(0) < min_reads);
        if result.k > 1 && is_under_min_reads {
            continue;
        }
        if best.as_ref().map_or(true, |(bic, _, _, _)| result.bic < *bic) {
            best = Some((result.bic, result, is_observed, cells));
        }
    }
    let (_, result, is_observed, cells) = best.unwrap();

    // HashMap<cell index, HashSet<read ID of another cell that this cell fits as well>>
    let mut shared: HashMap<usize, HashSet<usize>> = HashMap::new();
    for ((read_id, cell), row) in row_read_ids.iter().zip(result.labels.iter()).zip(matrix.iter()) {
        for other in cells.keys().filter(|other| *other != cell) {
            if !tells(row, *cell, *other, &result.consensus, &is_observed) {
                shared.entry(*other).or_default().insert(*read_id);
            }
        }
    }

    (cells, shared)
}


#[cfg(test)]
#[path = "../tests/phasing/genotype_phasing.rs"]
mod tests;