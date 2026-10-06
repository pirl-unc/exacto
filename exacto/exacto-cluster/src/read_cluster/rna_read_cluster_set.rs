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


use bimap::BiMap;
use exacto_caller::prelude::*;
use std::collections::{HashMap, HashSet};

use crate::prelude::*;


#[derive(Debug, Clone, PartialEq)]
pub struct FailedRNAVariantCall {
    pub cluster_id: usize,
    pub graph_operation: GraphOperation,
    pub num_reads: u32,
    pub total_depth: i32,
    pub failure: VariantCallFailure
}


#[derive(Debug,Clone)]
pub struct RNAReadClusterSet {
    /// HashMap<cluster ID, RNAReadCluster>
    clusters: HashMap<usize, RNAReadCluster>,

    /// BiMap<read name, read ID>
    read_names_map: BiMap<Box<str>, usize>,

    /// BiMap<chromosome name, chromosome ID>
    chromosome_names_map: BiMap<Box<str>, u16>,

    failed_variant_calls: Vec<FailedRNAVariantCall>
}

impl RNAReadClusterSet {
    pub fn new(
        read_names_map: BiMap<Box<str>, usize>,
        chromosome_names_map: BiMap<Box<str>, u16>
    ) -> Self {
        Self {
            clusters: HashMap::new(),
            read_names_map: read_names_map,
            chromosome_names_map: chromosome_names_map,
            failed_variant_calls: Vec::new()
        }
    }

    pub fn add_cluster(&mut self, cluster: RNAReadCluster) {
        assert_eq!(self.clusters.contains_key(&cluster.get_id()), false);
        self.clusters.insert(cluster.get_id(), cluster);
    }

    pub fn add_failed_variant_call(&mut self, failed_variant_call: FailedRNAVariantCall) {
        self.failed_variant_calls.push(failed_variant_call);
    }

    pub fn get_chromosome_names_map(&self) -> &BiMap<Box<str>, u16> {
        &self.chromosome_names_map
    }

    pub fn get_cluster(&self, id: usize) -> &RNAReadCluster {
        self.clusters.get(&id).unwrap()
    }
    
    pub fn get_clusters(&self) -> Vec<&RNAReadCluster> {
        self.clusters.values().collect()
    }
    
    pub fn get_cluster_reads_names(&self) -> HashMap<usize, HashSet<Box<str>>> {
        let mut cluster_reads_names: HashMap<usize, HashSet<Box<str>>> = HashMap::new();
        for (cluster_id, cluster) in self.clusters.iter() {
            let read_names: HashSet<Box<str>> = cluster.get_read_names(self.get_read_names_map());
            cluster_reads_names.insert(*cluster_id, read_names);
        }
        cluster_reads_names
    }
    
    pub fn get_failed_variant_calls(&self) -> &Vec<FailedRNAVariantCall> {
        &self.failed_variant_calls
    }

    pub fn get_introns(&self) -> HashMap<usize, HashMap<Box<str>, Vec<(u32, u32)>>> {
        // HashMap<cluster ID, HashMap<chromosome, Vec<(start, end)>>>
        let mut cluster_introns: HashMap<usize, HashMap<Box<str>, Vec<(u32, u32)>>> = HashMap::new();
        for cluster in self.get_clusters().iter() {
            // Get cluster introns.
            let introns = cluster_introns.entry(cluster.get_id()).or_default();
            for splice_junction in cluster.get_splice_junctions().iter() {
                if splice_junction.chromosome_1 == splice_junction.chromosome_2 {
                    let chromosome: Box<str> = self.chromosome_names_map
                        .get_by_right(&splice_junction.chromosome_1)
                        .unwrap()
                        .clone();
                    introns.entry(chromosome).or_default().push((
                        splice_junction.position_1.min(splice_junction.position_2),
                        splice_junction.position_1.max(splice_junction.position_2)
                    ));
                }
            }
            for spans in introns.values_mut() {
                spans.sort_unstable();
                spans.dedup();
            }
        }
        cluster_introns
    }
    
    pub fn get_read_clusters_map(&self) -> HashMap<Box<str>, usize> {
        let mut map: HashMap<Box<str>, usize> = HashMap::new();
        for cluster in self.get_clusters().iter() {
            // Make sure the assignments are disjoint.
            for read_id in cluster.get_read_ids().difference(cluster.get_shared_read_ids()) {
                let read_name: &Box<str> = self.get_read_names_map().get_by_right(read_id).unwrap();
                let previous = map.insert(read_name.clone(), cluster.get_id());
                assert!(
                    previous.is_none(),
                    "Read {} is assigned to clusters {} and {}; assignments must be disjoint.",
                    read_name, previous.unwrap(), cluster.get_id()
                );
            }
        }
        map
    }
    
    pub fn set_reference_gene_transcript_ids(
        &mut self,
        cluster_id: usize,
        reference_gene_transcript_ids: HashSet<(Box<str>, Box<str>)>
    ) {
        if let Some(cluster) = self.clusters.get_mut(&cluster_id) {
            cluster.set_reference_gene_transcript_ids(reference_gene_transcript_ids);
        }
    }

    pub fn get_read_names_map(&self) -> &BiMap<Box<str>, usize> {
        &self.read_names_map
    }
}
