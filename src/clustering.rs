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


extern crate bstr;
extern crate edit_distance;
extern crate log;
extern crate noodles;
extern crate noodles_core;
extern crate noodles_sam;
extern crate rayon;
extern crate regex;
extern crate serde;
use bstr::ByteSlice;
use edit_distance::edit_distance;
use log::info;
use noodles::bam as bam;
use noodles_core::{Region, Position};
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::value::Value;
use rayon::prelude::*;
use regex::Regex;
use serde::{Serialize, Deserialize};
use std::collections::{HashMap, HashSet};
use std::process;
// use std::sync::{Arc, atomic::{AtomicUsize, Ordering}};
use crate::constants::*;
use crate::defaults::*;
use crate::utilities::get_chromosomes;
use crate::variant_call::VariantCall;
use crate::variant_record::VariantRecord;


struct UnionFind {
    // key      =   child ID
    // value    =   parent ID
    parents: HashMap<String, String>,

    // key      =   parent ID
    // value    =   number of children
    sizes: HashMap<String, u32>
}

impl UnionFind {
    fn new() -> Self {
        UnionFind {
            parents: HashMap::new(),
            sizes: HashMap::new()
        }
    }

    fn get_size(&self, x: &str) -> u32 {
        *self.sizes.get(x).unwrap_or(&0)
    }

    fn find(&mut self, x: &str) -> String {
        if self.parents.contains_key(x) == false {
            self.parents.insert(x.to_string(), x.to_string());
            self.sizes.insert(x.to_string(), 1);
            return x.to_string();
        }

        // Compress path
        let mut path = Vec::new();
        let mut current = x;
        while let Some(parent) = self.parents.get(current) {
            if parent == current {
                break;
            }
            path.push(current.to_string());
            current = parent;
        }

        let root = current.to_string();
        for ancestor in path {
            self.parents.insert(ancestor, root.clone());
        }
        return root;
    }

    fn union(&mut self, x: &str, y: &str) {
        let root_x = self.find(x);
        let root_y = self.find(y);
        if root_x != root_y {
            let size_x = self.get_size(&root_x);
            let size_y = self.get_size(&root_y);
            if size_x < size_y {
                self.parents.insert(root_x.to_string(), root_y.to_string());
                self.sizes.entry(root_y.to_string()).and_modify(|e| *e += size_x);
            } else {
                self.parents.insert(root_y.clone(), root_x.clone());
                self.sizes.entry(root_x.clone()).and_modify(|e| *e += size_y);
            }
        }
    }

    fn clusters(&mut self) -> Vec<Vec<String>> {
        let keys: Vec<String> = self.parents.keys().cloned().collect();
        let mut map: HashMap<String, Vec<String>> = HashMap::new();
        for key in keys {
            let parent = self.find(&key);
            map.entry(parent)
                .or_insert_with(Vec::new)
                .push(key);
        }
        map.into_values().collect()
    }
}

pub fn can_cluster_variant_records(
    variant_record_1: &VariantRecord,
    variant_record_2: &VariantRecord,
    max_distance: isize,
    min_ins_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    min_del_size_proportion: f32) -> bool {
    // Make sure the variant types match
    if variant_record_1.variant_type != variant_record_2.variant_type {
        return false;
    }
    // Calculate distances
    let distance_1: isize = (variant_record_1.position_1 as isize - variant_record_2.position_1 as isize).abs();
    let distance_2: isize = (variant_record_1.position_2 as isize - variant_record_2.position_2 as isize).abs();

    // Compute size proportion between the two variant records
    let mut size_proportion: f32 = 1.0;
    if variant_record_1.variant_size < variant_record_2.variant_size {
        size_proportion = (variant_record_1.variant_size as f32) / (variant_record_1.variant_size as f32);
    } else {
        size_proportion = (variant_record_2.variant_size as f32) / (variant_record_1.variant_size as f32);
    }

    let mut cluster: bool = false;
    if variant_record_1.variant_type == SINGLE_NUCLEOTIDE_VARIANT {
        if (distance_1 == 0) &&
            (distance_2 == 0) &&
            (variant_record_1.alternate_allele == variant_record_2.alternate_allele) {
            cluster = true;
        }
    } else if variant_record_1.variant_type == INSERTION {
        let edit_distance: usize = edit_distance(&variant_record_1.alternate_allele, &variant_record_2.alternate_allele);
        let mut norm_edit_distance: f32 = 0.0;
        if variant_record_1.variant_size < variant_record_2.variant_size {
            norm_edit_distance = (edit_distance as f32) / (variant_record_2.variant_size as f32);
        } else {
            norm_edit_distance = (edit_distance as f32) / (variant_record_1.variant_size as f32);
        }
        if (distance_1 <= max_distance) &&
            (distance_2 <= max_distance) &&
            (size_proportion >= min_ins_size_proportion) &&
            (norm_edit_distance <= max_ins_norm_edit_distance) {
           cluster = true;
        }
    } else if variant_record_1.variant_type == DELETION {
        if (distance_1 <= max_distance) &&
            (distance_2 <= max_distance) &&
            (size_proportion >= min_del_size_proportion) {
            cluster = true;
        }
    } else if variant_record_1.variant_type == SPLICING {
        if (distance_1 == 0) && (distance_2 == 0) {
            cluster = true;
        }
    } else if variant_record_1.variant_type == BREAKPOINT {
        if (distance_1 <= max_distance) &&
            (distance_2 <= max_distance) {
            cluster = true;
        }
    } else {
        eprintln!("Unknown variant type: {}", variant_record_1.variant_type);
        std::process::exit(exitcode::DATAERR);
    }
    return cluster;
}

pub fn cluster_variant_records(
    variant_records: &mut Vec<VariantRecord>,
    sample_id: &str,
    nucleic_acid: &str,
    min_reads: usize,
    num_threads: usize,
    min_ins_size_proportion: f32,
    max_ins_norm_edit_distance: f32,
    min_del_size_proportion: f32,
    max_bnd_distance: isize,
    grid_size: isize
) -> Vec<VariantCall> {
    let thread_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()
        .unwrap();

    // Step 1. Create grids
    // key      =   (position_1,position_2,variant_type)
    // value    =   Vec<&VariantRecord>
    info!("Started creating a grid map");
    let mut grid_map: HashMap<(isize,isize,String), Vec<&VariantRecord>> = HashMap::new();
    for variant_record in variant_records.iter() {
        let cell_key: (isize,isize,String) = (variant_record.position_1 as isize / grid_size,
                                              variant_record.position_2 as isize / grid_size,
                                              variant_record.variant_type.to_string());
        grid_map.entry(cell_key)
                .or_insert(Vec::new())
                .push(variant_record);
    }
    let num_cells: usize = grid_map.len();
    info!("Finished creating a grid map of {} cells (keys)", num_cells);

    // Step 2. Sort Vec<VariantRecord> by VariantRecord.position_1
    thread_pool.install(|| {
        grid_map.par_iter_mut().for_each(|(_key, variant_records_)| {
            variant_records_.sort_by(|a, b| a.position_1.cmp(&b.position_1));
        });
    });

    // Step 3. Create VariantCall objects
    info!("Started creating VariantCall objects");
//     let counter = Arc::new(AtomicUsize::new(0)); // Shared counter
    let variant_calls_list: Vec<Vec<VariantCall>> = thread_pool.install(|| {
        grid_map.par_iter().map(|(cell_key, variant_records_)| {
            let variant_type = &cell_key.2;
            // Create a HashMap of variant records
            // key      =   VariantRecord.id
            // value    =   VariantRecord
            let mut map: HashMap<&str, &VariantRecord> = HashMap::new();
            let mut uf = UnionFind::new();
            // Cluster within the cell
            for i in 0..variant_records_.len() {
                map.insert(variant_records_[i].id.as_str(), variant_records_[i]);
                uf.union(&variant_records_[i].id, &variant_records_[i].id);
                // Maximum clustering distance is a function of the variant size
                // if the variant type is SNV, INS or DEL. If the variant type
                // is BND, then use max_bnd_distance.
                let mut max_distance: isize = 0;
                if variant_type == BREAKPOINT {
                    max_distance = max_bnd_distance;
                } else if variant_type == SINGLE_NUCLEOTIDE_VARIANT || variant_type == SPLICING {
                    max_distance = 0;
                } else {
                    max_distance = (variant_records_[i].variant_size as f32).log2().floor() as isize;
                }
                for j in (i + 1)..variant_records_.len() {
                    if max_distance == 0 {
                        if variant_records_[i].position_1 != variant_records_[j].position_1 {
                            break;
                        }
                    }
                    // Break loop if distance_1 is greater than the
                    // maximum distance allowed
                    let distance_1: isize = (variant_records_[i].position_1 as isize - variant_records_[j].position_1 as isize).abs();
                    if distance_1 > max_distance {
                        break;
                    }
                    if can_cluster_variant_records(
                        variant_records_[i],
                        variant_records_[j],
                        max_distance,
                        min_ins_size_proportion,
                        max_ins_norm_edit_distance,
                        min_del_size_proportion) {
                        uf.union(&variant_records_[i].id, &variant_records_[j].id);
                    }
                }
            }
            // Cluster with directly adjacent cells
            // VariantRecord objects in boundaries of cells could be merged
            // [0,-1], [0,1], [1,-1], [1,0], [1,1]
            let grid_row: isize = cell_key.0;
            let grid_col: isize = cell_key.1;
            let variant_type: &str = cell_key.2.as_str();
            if variant_type != SINGLE_NUCLEOTIDE_VARIANT && variant_type != SPLICING {
                for i in 0..=1 {
                    for j in -1..=1 {
                        if i == 0 && j == 0 {
                            continue;
                        }
                        let adj_cell_key: (isize,isize,String) = (grid_row + i, grid_col + j, variant_type.to_string());
                        if grid_map.contains_key(&adj_cell_key) == false {
                            continue;
                        }
                        if let Some(adj_variant_records) = grid_map.get(&adj_cell_key) {
                            for variant_record_1 in variant_records_.iter() {
                                // Maximum clustering distance is a function of the variant size
                                // if the variant type is SNV, INS or DEL. If the variant type
                                // is BND, then use max_bnd_distance.
                                let mut max_distance: isize = 0;
                                if variant_type == BREAKPOINT {
                                    max_distance = max_bnd_distance;
                                } else {
                                    max_distance = (variant_record_1.variant_size as f32).log2().floor() as isize;
                                }
                                for variant_record_2 in adj_variant_records.iter() {
                                    // Break loop if distance_1 is greater than the
                                    // maximum distance allowed
                                    let distance_1: isize = (variant_record_1.position_1 as isize - variant_record_2.position_1 as isize).abs();
                                    if distance_1 > max_distance {
                                        break;
                                    }
                                    if can_cluster_variant_records(
                                        variant_record_1,
                                        variant_record_2,
                                        max_distance,
                                        min_ins_size_proportion,
                                        max_ins_norm_edit_distance,
                                        min_del_size_proportion) {
                                        map.insert(variant_record_2.id.as_str(), variant_record_2);
                                        uf.union(&variant_record_1.id, &variant_record_2.id);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // Create common VariantCall objects
            let mut variant_calls: Vec<VariantCall> = Vec::new();
            for variant_record_ids in uf.clusters().iter() {
                if variant_record_ids.len() < min_reads {
                    // Not enough read support
                    continue;
                }
                // Create data structures for VariantCall member variables
                let mut position_1_values: Vec<f32> = Vec::new();
                let mut position_2_values: Vec<f32> = Vec::new();
                let mut alternate_allele_read_ids: Vec<&str> = Vec::new();
                let mut variant_sizes: Vec<f32> = Vec::new();
                let mut reference_allele: &str = "";
                let mut alternate_allele: &str = "";
                // Iterate over clustered VariantRecord objects
                // and append necessary data
                for variant_record_id in variant_record_ids.iter() {
                    let variant_record: &VariantRecord = map.get(variant_record_id.as_str()).unwrap();
                    position_1_values.push(variant_record.position_1 as f32);
                    position_2_values.push(variant_record.position_2 as f32);
                    alternate_allele_read_ids.push(variant_record.read_id.as_str());
                    variant_sizes.push(variant_record.variant_size as f32);
                    // NEED TO UPDATE THIS LATER TO GET ALL ALLELES
                    // OR CALL THE CONSENSUS ALLELE
                    if variant_record.reference_allele.chars().count() > reference_allele.chars().count() {
                        reference_allele = variant_record.reference_allele.as_str();
                    }
                    if variant_record.alternate_allele.chars().count() > alternate_allele.chars().count() {
                        alternate_allele = variant_record.alternate_allele.as_str();
                    }
                }
                // Compute average values
                let position_1_sum: f32 = position_1_values.iter().sum();
                let position_1_count = position_1_values.len() as f32;
                let position_1_average: isize = (position_1_sum / position_1_count).round() as isize;
                let position_2_sum: f32 = position_2_values.iter().sum();
                let position_2_count = position_2_values.len() as f32;
                let position_2_average: isize = (position_2_sum / position_2_count).round() as isize;
                let variant_size_sum: f32 = variant_sizes.iter().sum();
                let variant_size_count = variant_sizes.len() as f32;
                let variant_size_average: isize = (variant_size_sum / variant_size_count).round() as isize;
                // Create VariantCall ID
                let chromosome_1: &str = &map.get(variant_record_ids[0].as_str()).unwrap().chromosome_1;
                let chromosome_2: &str = &map.get(variant_record_ids[0].as_str()).unwrap().chromosome_2;
                let variant_type: &str = &map.get(variant_record_ids[0].as_str()).unwrap().variant_type;
                let variant_call_id: String = format!(
                    "{}_{}_exacto_{}:{}_{}:{}_{}",
                    sample_id,
                    nucleic_acid,
                    chromosome_1.to_string(),
                    position_1_average.to_string(),
                    chromosome_2.to_string(),
                    position_2_average.to_string(),
                    variant_type.to_string(),
                ).to_string();
                // Create a VariantCall object
                let mut variant_call = VariantCall::new(
                    variant_call_id.to_string(),
                    sample_id.to_string(),
                    nucleic_acid.to_string(),
                    chromosome_1.to_string(),
                    position_1_average as u32,
                    chromosome_2.to_string(),
                    position_2_average as u32,
                    variant_type.to_string(),
                    reference_allele.to_string(),
                    alternate_allele.to_string(),
                    variant_size_average as u32
                );
                // Add all read IDs
                for alternate_allele_read_id in alternate_allele_read_ids.iter() {
                    variant_call.add_alternate_allele_read_id(alternate_allele_read_id.to_string());
                }
                variant_calls.push(variant_call);
            }
//             let current_count = counter.fetch_add(1, Ordering::SeqCst) + 1; // fetch_add returns the previous value
//             if current_count % 10000 == 0 {
//                 info!("Processed {}/{} cells", current_count, num_cells);
//             }
            variant_calls
        })
        .collect()
    });
    info!("Finished creating VariantCall objects");

    // Consolidate VariantCall objects
    info!("Started consolidating VariantCall objects");
    let mut unique_variant_call_ids: HashSet<&str> = HashSet::new();
    let mut unique_variant_calls: Vec<VariantCall> = Vec::new();
    for variant_calls in variant_calls_list.iter() {
        for variant_call in variant_calls.iter() {
            if unique_variant_call_ids.contains(variant_call.id.as_str()) == false {
                unique_variant_calls.push(variant_call.clone());
                unique_variant_call_ids.insert(variant_call.id.as_str());
            }
        }
    }
//     let unique_variant_calls: Vec<VariantCall> = variant_calls_list
//         .into_iter()
//         .flat_map(|vec| vec.into_iter())                                    // Flatten the Vec<Vec<VariantCall>> to Vec<VariantCall>
//         .filter(|variant_call| unique_ids.insert(variant_call.id.clone()))  // Filter unique VariantCalls by id
//         .collect();                                                         // Collect the unique VariantCalls into a Vec

    return unique_variant_calls;
}