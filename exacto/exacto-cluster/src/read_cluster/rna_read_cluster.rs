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
use serde::{Serialize, Deserialize};
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};

use crate::prelude::*;


#[derive(Debug,Eq,Serialize,Deserialize)]
pub struct RNAReadCluster {
    id: usize,

    read_ids: HashSet<usize>,

    /// The reads of `read_ids` assigned to another cluster that this one fits as well.
    shared_read_ids: HashSet<usize>,
    
    splice_junctions: Vec<SpliceJunction>,

    variant_calls: Vec<VariantCall>,
    
    /// HashMap<variant call ID, HashMap<read ID, Allele>>
    genotypes: HashMap<usize, HashMap<usize, Allele>>,

    /// HashSet<(reference gene ID, reference transcript ID)>
    reference_gene_transcript_ids: HashSet<(Box<str>, Box<str>)>
}

impl PartialEq for RNAReadCluster {
    fn eq(&self, other: &Self) -> bool {
        if self.id == other.id {
            true
        } else {
            false
        }
    }
}

impl Hash for RNAReadCluster {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

impl RNAReadCluster {
    pub fn new(
        id: usize,
        read_ids: HashSet<usize>,
        splice_junctions: Vec<SpliceJunction>,
        variant_calls: Vec<VariantCall>,
        genotypes: HashMap<usize, HashMap<usize, Allele>>,
        reference_gene_transcript_ids: HashSet<(Box<str>, Box<str>)>
    ) -> Self {
        Self {
            id: id,
            read_ids: read_ids,
            shared_read_ids: HashSet::new(),
            splice_junctions: splice_junctions,
            variant_calls: variant_calls,
            genotypes: genotypes,   
            reference_gene_transcript_ids: reference_gene_transcript_ids
        }
    }

    pub fn with_shared_read_ids(mut self, shared_read_ids: HashSet<usize>) -> Self {
        assert!(shared_read_ids.is_subset(&self.read_ids), "Shared reads of cluster {} must be its reads.", self.id);
        self.shared_read_ids = shared_read_ids;
        self
    }

    pub fn get_id(&self) -> usize {
        self.id
    }
    
    pub fn get_genotypes(&self) -> &HashMap<usize, HashMap<usize, Allele>> {
        &self.genotypes
    }

    pub fn get_read_ids(&self) -> &HashSet<usize> {
        &self.read_ids
    }

    pub fn get_shared_read_ids(&self) -> &HashSet<usize> {
        &self.shared_read_ids
    }

    pub fn get_read_names(&self, read_names_map: &BiMap<Box<str>, usize>) -> HashSet<Box<str>> {
        let mut read_names = HashSet::new();
        for read_id in self.read_ids.iter() {
            let read_name = read_names_map.get_by_right(read_id).unwrap();
            read_names.insert(read_name.clone());
        }
        read_names
    }

    pub fn set_reference_gene_transcript_ids(&mut self, reference_gene_transcript_ids: HashSet<(Box<str>, Box<str>)>) {
        self.reference_gene_transcript_ids = reference_gene_transcript_ids;
    }

    pub fn get_reference_gene_transcript_ids(&self) -> &HashSet<(Box<str>, Box<str>)> {
        &self.reference_gene_transcript_ids
    }

    pub fn get_splice_junctions(&self) -> &Vec<SpliceJunction> {
        &self.splice_junctions
    }

    pub fn get_variant_calls(&self) -> &Vec<VariantCall> {
        &self.variant_calls
    }
}

impl Clone for RNAReadCluster {
    fn clone(&self) -> Self {
        RNAReadCluster {
            id: self.id,
            read_ids: self.read_ids.clone(),
            shared_read_ids: self.shared_read_ids.clone(),
            splice_junctions: self.splice_junctions.clone(),
            variant_calls: self.variant_calls.clone(),
            reference_gene_transcript_ids: self.reference_gene_transcript_ids.clone(),
            genotypes: self.genotypes.clone()
        }
    }
}
