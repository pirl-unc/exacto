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


use std::any::Any;
use std::collections::HashMap;
use std::sync::Arc;

use crate::prelude::*;


#[derive(Debug)]
pub struct VarGraphReferenceNode {
    chromosome: Box<str>,
    start: u32,
    end: u32,
    sequence: Box<str>,
    annotations: HashMap<Box<str>, Arc<dyn Any + Send + Sync + 'static>>
}

impl VarGraphReferenceNode {
    pub fn new(chromosome: &str, start: u32, end: u32, sequence: &str) -> Self {
        assert!(start <= end, "{} <= {} not satisfied.", start, end);
        assert!(sequence.len() as u32 == (end - start + 1), "{} != {} - {} + 1", sequence.len(), end, start);
        VarGraphReferenceNode {
            chromosome: chromosome.into(),
            start: start,
            end: end,
            sequence: sequence.into(),
            annotations: HashMap::new()
        }
    }

    pub fn add_annotation(&mut self, key: &str, value: Arc<dyn Any + Send + Sync + 'static>) {
        self.annotations.insert(key.into(), value);
    }

    pub fn annotation_key_exists(&self, key: &str) -> bool {
        self.annotations.contains_key(key)
    }

    pub fn get_annotation(&self, key: &str) -> &Arc<dyn Any + Send + Sync + 'static> {
        self.annotations.get(key).expect("Annotation not found.")
    }

    pub fn get_annotations(&self) -> &HashMap<Box<str>, Arc<dyn Any + Send + Sync + 'static>> {
        &self.annotations
    }

    pub fn get_chromosome(&self) -> &str {
        &*self.chromosome
    }

    pub fn get_start(&self) -> u32 {
        self.start
    }

    pub fn get_end(&self) -> u32 {
        self.end
    }

    pub fn get_sequence(&self) -> &str {
        &*self.sequence
    }
    
    pub fn get_sequence_length(&self) -> u32 {
        self.sequence.len() as u32
    }

    pub fn set_annotations(&mut self, annotations: &HashMap<Box<str>, Arc<dyn Any + Send + Sync + 'static>>) {
        for (key, value) in annotations {
            self.annotations.insert(key.clone(), Arc::clone(value));
        }
    }
}

impl Clone for VarGraphReferenceNode {
    fn clone(&self) -> VarGraphReferenceNode {
        VarGraphReferenceNode {
            chromosome: self.chromosome.clone(),
            start: self.start,
            end: self.end,
            sequence: self.sequence.clone(),
            annotations: self.annotations.clone()
        }
    }
}
