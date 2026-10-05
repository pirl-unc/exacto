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


use serde::{Serialize, Deserialize};
use std::collections::HashSet;
use std::hash::Hasher;


#[derive(Debug,Serialize,Deserialize)]
pub struct ClusterQuantification {
    pub cluster_id: usize,

    /// Counts per million.
    pub cpm: f64,

    /// Reads in this cluster.
    pub read_names: HashSet<Box<str>>
}

impl ClusterQuantification {
    pub fn new(
        cluster_id: usize,
        cpm: f64,
        read_names: HashSet<Box<str>>
    ) -> Self {
        Self {
            cluster_id: cluster_id,
            cpm: cpm,
            read_names: read_names
        }
    }

    pub fn get_cluster_id(&self) -> usize {
        self.cluster_id
    }

    pub fn get_cpm(&self) -> f64 {
        self.cpm
    }

    pub fn get_read_names(&self) -> &HashSet<Box<str>> {
        &self.read_names
    }
    
    pub fn get_read_count(&self) -> usize {
        self.read_names.len()
    }
}

impl Clone for ClusterQuantification {
    fn clone(&self) -> Self {
        ClusterQuantification {
            cluster_id: self.cluster_id,
            cpm: self.cpm,
            read_names: self.read_names.clone()
        }
    }
}
