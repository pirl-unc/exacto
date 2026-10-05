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
use exacto_core::prelude::{ClusterID, ReadID, ReadName, VariantID};
use std::collections::{HashMap, HashSet};

use crate::prelude::*;


pub struct RNAReadVariantCallSet {
    /// BiMap<read name, read ID>
    pub read_names_map: BiMap<ReadName, ReadID>,

    /// HashMap<phased cluster ID, HashSet<read ID>
    pub cluster_read_ids: HashMap<ClusterID, HashSet<ReadID>>,

    /// HashMap<read ID, phased cluster ID>
    pub read_cluster_ids: HashMap<ReadID, ClusterID>,

    /// HashMap<phased cluster ID, Vec<SpliceJunction>>
    pub cluster_junctions: HashMap<ClusterID, Vec<SpliceJunction>>,

    /// HashMap<phased cluster ID, Vec<VariantCall>>
    pub cluster_variants: HashMap<ClusterID, Vec<VariantCall>>,

    /// HashMap<variant call ID, Vec<(read ID, Allele)>>
    pub genotypes: HashMap<VariantID, Vec<(ReadID, Allele)>>
}

impl RNAReadVariantCallSet {
    pub fn new() -> Self {
        RNAReadVariantCallSet {
            read_names_map: BiMap::new(),
            cluster_read_ids: HashMap::new(),
            read_cluster_ids: HashMap::new(),
            cluster_junctions: HashMap::new(),
            cluster_variants: HashMap::new(),
            genotypes: HashMap::new()
        }
    }

    pub fn add_genotype(
        &mut self,
        variant_call_id: VariantID,
        read_id: ReadID,
        allele: Allele
    ) {
        self.genotypes
            .entry(variant_call_id)
            .or_insert_with(Vec::new)
            .push((read_id, allele));
    }

    pub fn add_junctions(&mut self, cluster_id: ClusterID, junctions: Vec<SpliceJunction>) {
        self.cluster_junctions.insert(cluster_id, junctions);
    }

    pub fn add_read(&mut self, read_name: &str, read_id: ReadID, cluster_id: ClusterID) {
        assert_eq!(self.read_names_map.contains_left(read_name), false);
        assert_eq!(self.read_names_map.contains_right(&read_id), false);
        self.read_names_map.insert(read_name.into(), read_id);
        self.cluster_read_ids.entry(cluster_id).or_insert(HashSet::new()).insert(read_id);
        self.read_cluster_ids.insert(read_id, cluster_id);
    }

    pub fn add_variant(&mut self, variant_call: VariantCall, cluster_id: ClusterID) {
        self.cluster_variants
            .entry(cluster_id)
            .or_insert(Vec::new())
            .push(variant_call);
    }
}

impl Clone for RNAReadVariantCallSet {
    fn clone(&self) -> Self {
        RNAReadVariantCallSet {
            read_names_map: self.read_names_map.clone(),
            cluster_read_ids: self.cluster_read_ids.clone(),
            read_cluster_ids: self.read_cluster_ids.clone(),
            cluster_junctions: self.cluster_junctions.clone(),
            cluster_variants: self.cluster_variants.clone(),
            genotypes: self.genotypes.clone()
        }
    }
}
